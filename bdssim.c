/* bdssim.c : BeiDou B1I baseband generator */
#include <stdio.h>
#include <string.h>
#include <stdint.h>
#include <math.h>
#include <stdlib.h>
#include <omp.h>
#include "bdssim.h"
#include "channel.h"
#include "coord.h"
#include "orbits.h"
#include "navbits.h"
#include "timeconv.h"
#include "path.h"
#include "globals.h"     /* nav_week */
#define OMEGA_E   7.2921150e-5
#define FSAMP_DEF 6.144e6

static const int geo[] ={1,2,3,4,5,59,60,61,62,63};
static int is_geo(int p){for(int i=0;i<10;i++) if(p==geo[i]) return 1; return 0;}

/* ---- Gaussian RNG (Box-Muller) ---- */
static double gauss_rand(void)
{
    static int have=0; static double next=0.0;
    if(have){ have=0; return next; }
    /* avoid log(0) by ensuring 0<u1<1 */
    double u1=((double)rand()+1.0)/((double)RAND_MAX+2.0);
    double u2=((double)rand()+1.0)/((double)RAND_MAX+2.0);
    double r=sqrt(-2.0*log(u1));
    double theta=2.0*M_PI*u2;
    next=r*sin(theta);
    have=1;
    return r*cos(theta);
}

/* ---- 檢查模擬起始時間與星曆 toe 差距是否過大 ---- */
static void check_ephemeris_age(int week,double sow)
{
    double t = week*604800.0 + sow;
    int warn = 0;
    for(int prn=1; prn<=63; ++prn){
        if(eph[prn].prn==0) continue;
        double toe = eph[prn].week*604800.0 + eph[prn].toe;
        double diff = fabs(t - toe);
        if(diff > 2*86400.0){
            fprintf(stderr,
                    "[warn] PRN%02d toe %.1f days away from start\n",
                    prn, diff/86400.0);
            warn = 1;
        }
    }
    if(warn)
        fputs("建議使用接近星曆 toe 的起始時間以避免 tk 錯誤\n", stderr);
}

/* ---- Rotate fixed LLH with Earth rotation ----- */
static void static_user_at(int week,double sow,const coord_t*ref,
                           coord_t*out,double vel[3])
{
    double lat=ref->llh[0]*(M_PI/180.0);
    double lon0=ref->llh[1]*(M_PI/180.0);
    double h=ref->llh[2];
    double sinp=sin(lat); double cosp=cos(lat);
    double N=WGS_A/sqrt(1.0-WGS_E2*sinp*sinp);
    double t=(week-nav_week)*604800.0 + sow;
    double lon=lon0 + OMEGA_E*t;
    double cosL=cos(lon), sinL=sin(lon);
    out->xyz[0]=(N+h)*cosp*cosL;
    out->xyz[1]=(N+h)*cosp*sinL;
    out->xyz[2]=(N*(1.0-WGS_E2)+h)*sinp;
    out->llh[0]=ref->llh[0]; out->llh[1]=ref->llh[1]; out->llh[2]=ref->llh[2];
    out->week=week; out->sow=sow;
    if(vel){
        vel[0]=-OMEGA_E*out->xyz[1];
        vel[1]= OMEGA_E*out->xyz[0];
        vel[2]=0.0;
    }
}

/*----------------------------------------------------*/
int select_channels(channel_t *ch,int *n,const coord_t*u)
{
    struct cand{int prn;double elev,rho,rdot;} c[63]; int m=0;
    double uv[3]={-OMEGA_E*u->xyz[1], OMEGA_E*u->xyz[0], 0.0};
    for(int prn=1;prn<=63;++prn){
        if(is_geo(prn))continue;
        double sat[3],vel[3];
        calc_sat_position_velocity(prn,u->week,u->sow,sat,vel);
        double enu[3]; ecef2enu(u,sat,enu);
        double el=enu_elevation_deg(enu); if(el<10)continue;
        double dx=sat[0]-u->xyz[0], dy=sat[1]-u->xyz[1], dz=sat[2]-u->xyz[2];
        double rho=hypot(hypot(dx,dy),dz);
        double rdot=(dx*(vel[0]-uv[0]) + dy*(vel[1]-uv[1]) + dz*(vel[2]-uv[2]))/rho;
        c[m++] = (struct cand){prn,el,rho,rdot};
    }
    /* sort by elevation (desc) */
    for(int i=0;i<m-1;++i) for(int j=i+1;j<m;++j)
        if(c[j].elev>c[i].elev){struct cand t=c[i];c[i]=c[j];c[j]=t;}
    *n = m<MAX_CH?m:MAX_CH;
    for(int i=0;i<*n;++i) channel_reset(&ch[i],c[i].prn);
    return *n;
}

/* 更新當前可見通道（可動態加入/移除） */
static void update_channels_dynamic(channel_t *ch,int *n,const coord_t *u,const double uvel[3],double gain)
{
    struct cand{int prn;double elev,rho,rdot;} cand[63];
    int m=0;
    double uv[3];
    if(uvel) { uv[0]=uvel[0]; uv[1]=uvel[1]; uv[2]=uvel[2]; }
    else { uv[0]=uv[1]=uv[2]=0.0; }
    for(int prn=1;prn<=63;++prn){
        if(is_geo(prn)) continue;
        double sat[3],vel[3];
        calc_sat_position_velocity(prn,u->week,u->sow,sat,vel);
        double enu[3]; ecef2enu(u,sat,enu);
        double el=enu_elevation_deg(enu); if(el<10.0) continue;
        double dx=sat[0]-u->xyz[0], dy=sat[1]-u->xyz[1], dz=sat[2]-u->xyz[2];
        double rho=hypot(hypot(dx,dy),dz);
        double rdot=(dx*(vel[0]-uv[0]) + dy*(vel[1]-uv[1]) + dz*(vel[2]-uv[2]))/rho;
        cand[m++] = (struct cand){prn,el,rho,rdot};
    }
    for(int i=0;i<m-1;++i) for(int j=i+1;j<m;++j)
        if(cand[j].elev>cand[i].elev){struct cand t=cand[i];cand[i]=cand[j];cand[j]=t;}

    int new_n = m<MAX_CH?m:MAX_CH;
    channel_t new_ch[MAX_CH];
    for(int i=0;i<new_n;++i){
        int idx=-1;
        for(int j=0;j<*n;++j) if(ch[j].prn==cand[i].prn){ idx=j; break; }
        if(idx>=0) new_ch[i]=ch[idx];
        else       channel_reset(&new_ch[i],cand[i].prn);
        update_channel_dynamics(&new_ch[i],cand[i].rho,cand[i].rdot,new_n,gain);
    }
    for(int i=0;i<new_n;++i) ch[i]=new_ch[i];
    *n = new_n;
}
/*----------------------------------------------------*/
void generate_signal(const sim_config_t *cfg)
{
    coord_t usr={0};
    path_t path={0};
    /* Initialise RNG for noise generation and random channel phase */
    srand(cfg->noise_seed);
    if(cfg->path_type==1)       load_path_xyz(cfg->path_file,&path);
    else if(cfg->path_type==2)  load_path_llh(cfg->path_file,&path);
    else if(cfg->path_type==3)  load_path_nmea(cfg->path_file,&path);
    if(cfg->path_type!=0 && path.n==0){fputs("path read error\n",stderr);return;}
    if(cfg->path_type==0)      llh2xyz(cfg->llh,&usr);
    else { interpolate_path(&path,0.0,&usr); xyz2llh(usr.xyz,&usr); }
    if(utc_to_bdt(cfg->time_start,&usr.week,&usr.sow)!=0){fputs("UTC format err\n",stderr);return;}

    coord_t ref_llh=usr;              /* 保存經緯度作旋轉基準 */
    static_user_at(usr.week,usr.sow,&ref_llh,&usr,NULL);

    navbits_init();

    /* 檢查 start 時間是否與星曆 toe 接近 */
    check_ephemeris_age(usr.week, usr.sow);

    channel_t ch[MAX_CH]; int n_ch; select_channels(ch,&n_ch,&usr);

    double uvel[3]={-OMEGA_E*usr.xyz[1], OMEGA_E*usr.xyz[0], 0.0};
    /* 首次幾何 – 初始化振幅/NCO */
    for(int i=0;i<n_ch;++i){
        double sat[3],vel[3];
        calc_sat_position_velocity(ch[i].prn,usr.week,usr.sow,sat,vel);
        double dx=sat[0]-usr.xyz[0],dy=sat[1]-usr.xyz[1],dz=sat[2]-usr.xyz[2];
        double rho=hypot(hypot(dx,dy),dz);
        double rdot=(dx*(vel[0]-uvel[0]) + dy*(vel[1]-uvel[1]) + dz*(vel[2]-uvel[2]))/rho;
        update_channel_dynamics(&ch[i],rho,rdot,n_ch,cfg->gain);
        printf("[ch%02d] rdot %.2f fd %.2fHz\n", ch[i].prn, rdot, ch[i].fd);
    }

    double fs = cfg->sample_rate ? cfg->sample_rate : FSAMP_DEF;
    int samp_per_ms = (int)(fs/1000.0 + 0.5);
    channel_set_fs(fs);

    FILE *fp=NULL, *fp8=NULL;
    if(cfg->byte_output){
        fp8=fopen("beidou_b1i_u8.bin","wb");
        if(!fp8){perror("u8"); return;}
    } else {
        fp=fopen("beidou_b1i.bin","wb");
        if(!fp){perror("bin"); return;}
    }
    int16_t tmpI[MAX_CH][samp_per_ms],tmpQ[MAX_CH][samp_per_ms];
    int32_t accI[samp_per_ms],accQ[samp_per_ms];
    double sumI=0.0,sumQ=0.0,sumI2=0.0,sumQ2=0.0;
    uint64_t samp_cnt=0;

    const int STEP_MS = cfg->step_ms;
    const uint64_t total_ms=(uint64_t)cfg->duration*1000;
    coord_t ref_usr=ref_llh; /* keep LLH */
    double start_bdt = usr.week*604800.0 + usr.sow;
    for(uint64_t ms=0; ms<total_ms; ms+=STEP_MS)
    {
        double t_abs = start_bdt + ms*0.001;
        int week = (int)(t_abs/604800.0);
        double sow = t_abs - week*604800.0;

        double uvel[3];
        if(cfg->path_type==0){
            static_user_at(week,sow,&ref_usr,&usr,uvel);
        } else {
            coord_t prev=usr;
            interpolate_path(&path, ms/1000.0, &usr);
            xyz2llh(usr.xyz,&usr);
            double dt=STEP_MS*0.001;
            uvel[0]=(usr.xyz[0]-prev.xyz[0])/dt;
            uvel[1]=(usr.xyz[1]-prev.xyz[1])/dt;
            uvel[2]=(usr.xyz[2]-prev.xyz[2])/dt;
        }

        update_channels_dynamic(ch,&n_ch,&usr,uvel,cfg->gain);

        /* --- STEP_MS 次 1ms 取樣 --- */
        for(int step=0;step<STEP_MS;++step){
            memset(accI,0,sizeof(accI)); memset(accQ,0,sizeof(accQ));

            /* 併行各通道 */
            #pragma omp parallel for
            for(int c=0;c<n_ch;++c)
                gen_samples_1ms(&ch[c],week,sow+step*0.001,
                               samp_per_ms,tmpI[c],tmpQ[c]);

            /* 歸併 */
            for(int c=0;c<n_ch;++c)
                for(int k=0;k<samp_per_ms;++k){
                    accI[k]+=tmpI[c][k];
                    accQ[k]+=tmpQ[c][k];
                }

            /* 限幅並打包成 I/Q (加入 AWGN) */
            int16_t iq[2*samp_per_ms];
            int8_t  i8[samp_per_ms];            /* byte output holds I only */
            for(int k=0;k<samp_per_ms;++k){
                int32_t i=accI[k];
                int32_t q=accQ[k];
                if(cfg->noise_std>0.0){
                    i += lrint(cfg->noise_std*gauss_rand());
                    q += lrint(cfg->noise_std*gauss_rand());
                }
                if(i>32767)i=32767; else if(i<-32768)i=-32768;
                if(q>32767)q=32767; else if(q<-32768)q=-32768;
                iq[2*k]   = (int16_t)i;
                iq[2*k+1] = (int16_t)q;
                i8[k]     = (int8_t)(iq[2*k]/256);
                sumI  += iq[2*k];
                sumQ  += iq[2*k+1];
                sumI2 += (double)iq[2*k]*iq[2*k];
                sumQ2 += (double)iq[2*k+1]*iq[2*k+1];
            }
            samp_cnt += samp_per_ms;
            if(fp)
                fwrite(iq,sizeof(int16_t),2*samp_per_ms,fp);
            else
                fwrite(i8,sizeof(int8_t),samp_per_ms,fp8);

        }
        /* 進度顯示 */
        printf("\r進度: %.2f / %.2f 秒",(ms+STEP_MS)/1000.0,total_ms/1000.0);
        fflush(stdout);
    }
    puts("");
    if(fp){
        fclose(fp);
        puts("[bdssim] 完成多星基帶輸出 beidou_b1i.bin");
    } else {
        fclose(fp8);
        puts("[bdssim] 完成 I-only 基帶輸出 beidou_b1i_u8.bin");
    }
    double meanI=sumI/samp_cnt, meanQ=sumQ/samp_cnt;
    double stdI = sqrt(sumI2/samp_cnt - meanI*meanI);
    double stdQ = sqrt(sumQ2/samp_cnt - meanQ*meanQ);
    printf("I mean=%.2f std=%.2f, Q mean=%.2f std=%.2f\n",
           meanI, stdI, meanQ, stdQ);
    free_path(&path);
}

