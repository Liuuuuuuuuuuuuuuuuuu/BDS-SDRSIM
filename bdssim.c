/* bdssim.c : 8.184 Msps 北斗 B1I 基帶產生器 */
#include <stdio.h>
#include <string.h>
#include <stdint.h>
#include <math.h>
#include <omp.h>
#include "bdssim.h"
#include "channel.h"
#include "coord.h"
#include "orbits.h"
#include "navbits.h"
#include "timeconv.h"
#include "path.h"
#define OMEGA_E   7.2921150e-5

/* 旋轉靜態使用者位置到當前 t，單位秒 ------------------------- */
static void update_static_xyz(coord_t *u, int week, double sow)
{
    double lat = u->llh[0] * (M_PI/180.0);
    double lon0= u->llh[1] * (M_PI/180.0);
    double h   = u->llh[2];
    double sinp= sin(lat); double cosp= cos(lat);
    double N = WGS_A / sqrt(1.0 - WGS_E2*sinp*sinp);
    double dt = (week - u->week)*604800.0 + (sow - u->sow);
    double lon = lon0 + OMEGA_E*dt;
    u->xyz[0] = (N + h)*cosp*cos(lon);
    u->xyz[1] = (N + h)*cosp*sin(lon);
    u->xyz[2] = ((1.0-WGS_E2)*N + h)*sinp;
}

/* 使用者在 ECEF 系下的速度（僅地球自轉） -------------------- */
static void user_ecef_velocity(const coord_t *u,double v[3])
{
    v[0] = -OMEGA_E*u->xyz[1];
    v[1] =  OMEGA_E*u->xyz[0];
    v[2] =  0.0;
}
#define FSAMP     8.184e6
#define SAMP_1MS  8184

static const int geo[] ={1,2,3,4,5,59,60,61,62,63};
static int is_geo(int p){for(int i=0;i<10;i++) if(p==geo[i]) return 1; return 0;}

/*----------------------------------------------------*/
int select_channels(channel_t *ch,int *n,const coord_t*u)
{
    struct cand{int prn;double elev,rho,rdot;} c[63]; int m=0;
    double uv[3]; user_ecef_velocity(u,uv);
    for(int prn=1;prn<=63;++prn){
        if(is_geo(prn))continue;
        double sat[3],vel[3]; calc_sat_position_velocity(prn,u->week,u->sow,sat,vel);
        double enu[3]; ecef2enu(u,sat,enu);
        double el=enu_elevation_deg(enu); if(el<10)continue;
        double dx=sat[0]-u->xyz[0], dy=sat[1]-u->xyz[1], dz=sat[2]-u->xyz[2];
        double rho=hypot(hypot(dx,dy),dz);
        double rdot=(dx*(vel[0]-uv[0]) + dy*(vel[1]-uv[1]) + dz*(vel[2]-uv[2]))/rho;
        c[m++] = (struct cand){prn,el,rho,rdot};
    }
    /* sort by elev desc */
    for(int i=0;i<m-1;++i) for(int j=i+1;j<m;++j)
        if(c[j].elev>c[i].elev){struct cand t=c[i];c[i]=c[j];c[j]=t;}
    *n = m<MAX_CH?m:MAX_CH;
    for(int i=0;i<*n;++i) channel_reset(&ch[i],c[i].prn);
    return *n;
}
/*----------------------------------------------------*/
void generate_signal(const sim_config_t *cfg)
{
    coord_t usr={0};
    path_t path={0};
    if(cfg->path_type==1)       load_path_xyz(cfg->path_file,&path);
    else if(cfg->path_type==2)  load_path_llh(cfg->path_file,&path);
    else if(cfg->path_type==3)  load_path_nmea(cfg->path_file,&path);
    if(cfg->path_type!=0 && path.n==0){fputs("path read error\n",stderr);return;}
    if(cfg->path_type==0)      llh2xyz(cfg->llh,&usr);
    else { interpolate_path(&path,0.0,&usr); xyz2llh(usr.xyz,&usr); }
    if(utc_to_bdt(cfg->time_start,&usr.week,&usr.sow)!=0){fputs("UTC format err\n",stderr);return;}

    navbits_init();

    channel_t ch[MAX_CH]; int n_ch; select_channels(ch,&n_ch,&usr);

    /* 首次幾何 – 初始化振幅/NCO */
    double uvel0[3]; user_ecef_velocity(&usr,uvel0);
    for(int i=0;i<n_ch;++i){
        double sat[3],vel[3]; calc_sat_position_velocity(ch[i].prn,usr.week,usr.sow,sat,vel);
        double dx=sat[0]-usr.xyz[0],dy=sat[1]-usr.xyz[1],dz=sat[2]-usr.xyz[2];
        double rho=hypot(hypot(dx,dy),dz);
        double rdot=(dx*(vel[0]-uvel0[0]) + dy*(vel[1]-uvel0[1]) + dz*(vel[2]-uvel0[2]))/rho;
        update_channel_dynamics(&ch[i],rho,rdot,n_ch);
    }

    FILE *fp=fopen("beidou_b1i.bin","wb"); if(!fp){perror("bin");return;}
    int16_t tmpI[MAX_CH][SAMP_1MS],tmpQ[MAX_CH][SAMP_1MS];
    int32_t accI[SAMP_1MS],accQ[SAMP_1MS];

    const int STEP_MS = cfg->step_ms;
    const uint64_t total_ms=(uint64_t)cfg->duration*1000;
    for(uint64_t ms=0; ms<total_ms; ms+=STEP_MS)
    {
        /* --- 幾何重算每 STEP_MS --- */
        double sow = usr.sow + ms*0.001;
        if(cfg->path_type!=0){
            interpolate_path(&path, ms/1000.0, &usr);
        }else{
            update_static_xyz(&usr, usr.week, sow);
        }
        double uvel[3]; user_ecef_velocity(&usr,uvel);
        for(int i=0;i<n_ch;++i){
            double sat[3],vel[3]; calc_sat_position_velocity(ch[i].prn,usr.week,sow,sat,vel);
            double dx=sat[0]-usr.xyz[0],dy=sat[1]-usr.xyz[1],dz=sat[2]-usr.xyz[2];
            double rho=hypot(hypot(dx,dy),dz);
            double rdot=(dx*(vel[0]-uvel[0]) + dy*(vel[1]-uvel[1]) + dz*(vel[2]-uvel[2]))/rho;
            update_channel_dynamics(&ch[i],rho,rdot,n_ch);
        }

        /* --- STEP_MS 次 1ms 取樣 --- */
        for(int step=0;step<STEP_MS;++step){
            memset(accI,0,sizeof(accI)); memset(accQ,0,sizeof(accQ));

            /* 併行各通道 */
            #pragma omp parallel for
            for(int c=0;c<n_ch;++c)
                gen_samples_1ms(&ch[c],usr.week,sow+step*0.001,
                               tmpI[c],tmpQ[c]);

            /* 歸併 */
            for(int c=0;c<n_ch;++c)
                for(int k=0;k<SAMP_1MS;++k){
                    accI[k]+=tmpI[c][k];
                    accQ[k]+=tmpQ[c][k];
                }

            /* 限幅並打包成 I/Q */
            int16_t iq[2*SAMP_1MS];
            for(int k=0;k<SAMP_1MS;++k){
                int32_t i=accI[k];
                int32_t q=accQ[k];
                if(i>32767)i=32767; else if(i<-32768)i=-32768;
                if(q>32767)q=32767; else if(q<-32768)q=-32768;
                iq[2*k]   = (int16_t)i;
                iq[2*k+1] = (int16_t)q;
            }
            fwrite(iq,sizeof(int16_t),2*SAMP_1MS,fp);

        }
        /* 進度顯示 */
        printf("\r進度: %.2f / %.2f 秒",(ms+STEP_MS)/1000.0,total_ms/1000.0);
        fflush(stdout);
    }
    puts(""); fclose(fp);
    puts("[bdssim] 完成多星基帶輸出 beidou_b1i.bin");
    free_path(&path);
}

