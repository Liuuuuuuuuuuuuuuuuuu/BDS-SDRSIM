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
#define FSAMP_DEF FS_OUTPUT_HZ    /* default 6.144 MHz for 16-bit I/Q output */
#define FSAMP_BYTE 25.0e6    /* 25 MHz when --byte is used */
#define FCARRIER   1561.098e6      /* B1I carrier */
#define CHIPRATE   2.046e6

/* ========= 振幅計算與 int16 飽和保護 ========= */
static inline int16_t saturate_int16(double x)
{
    if (x >  32760.0) return  32760;
    if (x < -32760.0) return -32760;
    return (int16_t)llround(x);
}

/* ---- 檢查模擬起始時間與星曆 toe 差距是否過大 ---- */
static void check_ephemeris_age(int week,double sow)
{
    double t = week*604800.0 + sow;
    int warn = 0;
    for(int prn=1; prn<=prn_max; ++prn){
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
int select_channels(channel_t *ch,int *n,const coord_t*u,
                    bool geo_first,int single_prn,bool no_geo)
{
    struct cand{int prn;double elev,rho,rdot;int pri;} c[63];
    int m=0;
    double uv[3]={-OMEGA_E*u->xyz[1], OMEGA_E*u->xyz[0], 0.0};
    for(int prn=1;prn<=prn_max;++prn){
        if(single_prn>0 && prn!=single_prn) continue;
        if(no_geo && is_d2_prn(prn)) continue;
        double sat[3],vel[3];
        calc_sat_position_velocity(prn,u->week,u->sow,sat,vel);
        double enu[3]; ecef2enu(u,sat,enu);
        double el=enu_elevation_deg(enu); if(el<5.0)continue;
        double dx=sat[0]-u->xyz[0], dy=sat[1]-u->xyz[1], dz=sat[2]-u->xyz[2];
        double rho=hypot(hypot(dx,dy),dz);
        double rdot=(dx*(vel[0]-uv[0]) + dy*(vel[1]-uv[1]) + dz*(vel[2]-uv[2]))/rho;
        int pri = (geo_first && is_d2_prn(prn)) ? 1 : 0;
        c[m++] = (struct cand){prn,el,rho,rdot,pri};
    }
    /* sort by priority and elevation (desc) */
    for(int i=0;i<m-1;++i) for(int j=i+1;j<m;++j){
        if(c[j].pri>c[i].pri ||
           (c[j].pri==c[i].pri && c[j].elev>c[i].elev)){
            struct cand t=c[i];c[i]=c[j];c[j]=t;
        }
    }
    *n = m<MAX_CH?m:MAX_CH;
    for(int i=0;i<*n;++i) channel_reset(&ch[i],c[i].prn,u->week,u->sow);
    return *n;
}

/* 更新通道但保留原始 PRN，不做重新選擇 */
static void update_channels_fixed(channel_t *ch,int n,const coord_t *u,
                                  const double uvel[3],double gain,double target_cn0,
                                  int step_ms)
{
    double uv[3];
    if(uvel){ uv[0]=uvel[0]; uv[1]=uvel[1]; uv[2]=uvel[2]; }
    else     { uv[0]=uv[1]=uv[2]=0.0; }

    double dt = step_ms * 0.001; /* seconds */

    for(int i=0;i<n;++i){
        int prn = ch[i].prn;
        double sat[3], vel[3];
        calc_sat_position_velocity(prn,u->week,u->sow,sat,vel);
        double enu[3]; ecef2enu(u,sat,enu);
        double el = enu_elevation_deg(enu);
        double dx=sat[0]-u->xyz[0], dy=sat[1]-u->xyz[1], dz=sat[2]-u->xyz[2];
        double rho=hypot(hypot(dx,dy),dz);
        double rdot=(dx*(vel[0]-uv[0]) + dy*(vel[1]-uv[1]) + dz*(vel[2]-uv[2]))/rho;
        update_channel_dynamics(&ch[i],rho,rdot,el,gain,target_cn0,n,step_ms);
        if(el < 5.0) ch[i].amp = 0.0;     /* below horizon → mute */

        /* predict next dynamics for linear interpolation */
        double u_next[3]={u->xyz[0]+uv[0]*dt,
                          u->xyz[1]+uv[1]*dt,
                          u->xyz[2]+uv[2]*dt};
        double sat2[3], vel2[3];
        double t_abs = u->week*604800.0 + u->sow + dt;
        int week2 = (int)(t_abs/604800.0);
        double sow2 = t_abs - week2*604800.0;
        calc_sat_position_velocity(prn,week2,sow2,sat2,vel2);
        double dx2=sat2[0]-u_next[0], dy2=sat2[1]-u_next[1], dz2=sat2[2]-u_next[2];
        double rho2=hypot(hypot(dx2,dy2),dz2);
        double rdot2=(dx2*(vel2[0]-uv[0]) + dy2*(vel2[1]-uv[1]) + dz2*(vel2[2]-uv[2]))/rho2;
        coord_t u2={0};
        xyz2llh(u_next,&u2);
        double enu2[3]; ecef2enu(&u2,sat2,enu2);
        double el2 = enu_elevation_deg(enu2);
        double fd2 = -FCARRIER*rdot2/299792458.0;
        double cr2 = CHIPRATE*(1.0 - rdot2/299792458.0);
        double A2 = predict_next_amp(&ch[i], el2, gain, target_cn0, n, step_ms);
        if(el2 < 5.0) A2 = 0.0;
        ch[i].fd_dot = (fd2 - ch[i].fd) / dt;
        ch[i].code_rate_dot = (cr2 - ch[i].code_rate) / dt;
        ch[i].amp_dot = (A2 - ch[i].amp) / dt;
    }
}
/*----------------------------------------------------*/
void generate_signal(const sim_config_t *cfg)
{
    coord_t usr={0};
    path_t path={0};
    /* Initialise RNG for random channel phase */
    srand(cfg->seed);
    if(cfg->path_type==1)       load_path_xyz(cfg->path_file,&path);
    else if(cfg->path_type==2)  load_path_llh(cfg->path_file,&path);
    else if(cfg->path_type==3)  load_path_nmea(cfg->path_file,&path);
    if(cfg->path_type!=0 && path.n==0){
        fputs("path read error\n",stderr);
        free_path(&path);
        return;
    }
    if(cfg->path_type==0)      llh2xyz(cfg->llh,&usr);
    else { interpolate_path(&path,0.0,&usr); xyz2llh(usr.xyz,&usr); }
    if(utc_to_bdt(cfg->time_start,&usr.week,&usr.sow)!=0){
        fputs("UTC format err\n",stderr);
        free_path(&path);
        return;
    }

    coord_t ref_llh=usr;              /* 保存經緯度作旋轉基準 */
    static_user_at(usr.week,usr.sow,&ref_llh,&usr,NULL);

    navbits_init();

    /* 檢查 start 時間是否與星曆 toe 接近 */
    check_ephemeris_age(usr.week, usr.sow);

    channel_t ch[MAX_CH];
    int n_ch;
    select_channels(ch,&n_ch,&usr,cfg->geo_first,
                    cfg->single_prn,cfg->no_geo);

    double fs = cfg->byte_output ? FSAMP_BYTE : FSAMP_DEF;
    int samp_per_ms = (int)(fs/1000.0 + 0.5);
    channel_set_fs(fs);                   /* ensure dynamics use correct Fs */

    FILE *fp=NULL, *fp8=NULL;
    if(cfg->byte_output){
        fp8=fopen("beidou_b1i_u8.bin","wb");
        if(!fp8){
            perror("u8");
            free_path(&path);
            return;
        }
    } else {
        fp=fopen("beidou_b1i.bin","wb");
        if(!fp){
            perror("bin");
            free_path(&path);
            return;
        }
    }
    int16_t tmpI[MAX_CH][samp_per_ms],tmpQ[MAX_CH][samp_per_ms];
    int32_t accI[samp_per_ms],accQ[samp_per_ms];
    double sumI=0.0,sumQ=0.0,sumI2=0.0,sumQ2=0.0;
    uint64_t samp_cnt=0;
    static uint64_t clip_cnt = 0, tot_cnt = 0;

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
        } else if(cfg->path_type==1){
            coord_t prev=usr;
            interpolate_path(&path, ms/1000.0, &usr);
            xyz2llh(usr.xyz,&usr);
            usr.week=week; usr.sow=sow;
            double prev_eci[3], cur_eci[3];
            ecef_to_eci(prev.xyz, prev.week, prev.sow, prev_eci);
            ecef_to_eci(usr.xyz, week, sow, cur_eci);
            double dt=STEP_MS*0.001;
            uvel[0]=(cur_eci[0]-prev_eci[0])/dt;
            uvel[1]=(cur_eci[1]-prev_eci[1])/dt;
            uvel[2]=(cur_eci[2]-prev_eci[2])/dt;
        } else {
            coord_t prev=usr;
            double llh[3];
            interpolate_path_llh(&path, ms/1000.0, llh);
            coord_t ref={0};
            ref.llh[0]=llh[0]; ref.llh[1]=llh[1]; ref.llh[2]=llh[2];
            static_user_at(week,sow,&ref,&usr,NULL);
            double prev_eci[3], cur_eci[3];
            ecef_to_eci(prev.xyz, prev.week, prev.sow, prev_eci);
            ecef_to_eci(usr.xyz, week, sow, cur_eci);
            double dt=STEP_MS*0.001;
            uvel[0]=(cur_eci[0]-prev_eci[0])/dt;
            uvel[1]=(cur_eci[1]-prev_eci[1])/dt;
            uvel[2]=(cur_eci[2]-prev_eci[2])/dt;
        }

        update_channels_fixed(ch,n_ch,&usr,uvel,
                              cfg->gain,g_target_cn0,STEP_MS);
        if(ms==0){
            for(int i=0;i<n_ch;++i){
                double rdot = -ch[i].fd*299792458.0/FCARRIER;
                printf("[ch%02d] rdot %.2f fd %.2fHz\n", ch[i].prn, rdot, ch[i].fd);
            }
        }

        /* --- STEP_MS 次 1ms 取樣 --- */
        for(int step=0;step<STEP_MS;++step){
            memset(accI,0,sizeof(accI)); memset(accQ,0,sizeof(accQ));

            /* 併行各通道 */
            #pragma omp parallel for
            for(int c=0;c<n_ch;++c){
                bool use_d2 = is_d2_prn(ch[c].prn);
                if(use_d2)
                    gen_samples_1ms_d2(&ch[c],week,sow+step*0.001,
                                       samp_per_ms,tmpI[c],tmpQ[c]);
                else
                    gen_samples_1ms(&ch[c],week,sow+step*0.001,
                                   samp_per_ms,tmpI[c],tmpQ[c]);
            }

            /* 歸併 */
            for(int c=0;c<n_ch;++c)
                for(int k=0;k<samp_per_ms;++k){
                    accI[k]+=tmpI[c][k];
                    accQ[k]+=tmpQ[c][k];
                }

            /* 限幅並打包成 I/Q */
            int16_t iq[2*samp_per_ms];
            int8_t  i8[samp_per_ms];            /* byte output holds I only */
            for(int k=0;k<samp_per_ms;++k){
                int32_t i=accI[k];
                int32_t q=accQ[k];
                if(fp){
                    if(i>32760 || i<-32760) clip_cnt++;
                    if(q>32760 || q<-32760) clip_cnt++;
                    tot_cnt += 2;
                    iq[2*k]   = saturate_int16((double)i);
                    iq[2*k+1] = saturate_int16((double)q);
                    sumI  += iq[2*k];
                    sumQ  += iq[2*k+1];
                    sumI2 += (double)iq[2*k]*iq[2*k];
                    sumQ2 += (double)iq[2*k+1]*iq[2*k+1];
                } else {
                    if (i>127) i=127; else if (i<-128) i=-128;
                    i8[k] = (int8_t)i;
                    sumI  += i8[k];
                    sumI2 += (double)i8[k]*i8[k];
                }
            }
            samp_cnt += samp_per_ms;
            if(fp)
                fwrite(iq,sizeof(int16_t),2*samp_per_ms,fp);
            else
                fwrite(i8,sizeof(int8_t),samp_per_ms,fp8);
            /* 每秒印一次 clipping 比例 */
            if (((int)(step+1))%1000==0 && fp){
                fprintf(stderr,"[stat] clip=%.5f%%\n", 100.0*(double)clip_cnt/(double)tot_cnt);
                clip_cnt=tot_cnt=0;
            }

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

