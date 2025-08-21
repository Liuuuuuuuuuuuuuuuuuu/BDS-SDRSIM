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
#include "globals.h"     /* nav_week, CLIGHT */
#define FSAMP_DEF FS_OUTPUT_HZ    /* default 6.144 MHz for 16-bit I/Q output */
#define FSAMP_BYTE 25.0e6    /* 25 MHz when --byte is used */
#define FCARRIER   1561.098e6      /* B1I carrier */
#define CHIPRATE   2.046e6
#define OMEGA_E    7.2921150e-5

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

/* Compute pseudorange and satellite position in ECEF */
static double compute_pseudorange(int prn,int week,double sow,
                                  const coord_t *u,double sat[3])
{
    double clk_tx[2], dummy[3];
    calc_sat_position_velocity(prn, week, sow, dummy, NULL, clk_tx);
    double t_sv = week*604800.0 + sow - clk_tx[0];
    int week_sv = (int)(t_sv/604800.0);
    double sow_sv = t_sv - week_sv*604800.0;

    double eci[3];
    calc_sat_position_velocity(prn, week_sv, sow_sv, eci, NULL, NULL);
    double ecef[3];
    eci_to_ecef(eci, week_sv, sow_sv, ecef);
    if(sat){ sat[0]=ecef[0]; sat[1]=ecef[1]; sat[2]=ecef[2]; }

    double dx = ecef[0]-u->xyz[0];
    double dy = ecef[1]-u->xyz[1];
    double dz = ecef[2]-u->xyz[2];
    double rho_geom = hypot(hypot(dx,dy),dz);
    double sag = OMEGA_E/CLIGHT * (ecef[0]*u->xyz[1] - ecef[1]*u->xyz[0]);
    return rho_geom + sag - CLIGHT * clk_tx[0];
}

/*----------------------------------------------------*/
int select_channels(channel_t *ch,int *n,const coord_t*u,
                    int single_prn,bool meo_only)
{
    struct cand{int prn;double elev,rho;} c[63];
    int m=0;
    for(int prn=1;prn<=prn_max;++prn){
        if(single_prn>0 && prn!=single_prn) continue;
        if(is_geo_prn(prn)) continue;
        if(meo_only && !is_meo_prn(prn)) continue;
        double sat[3];
        double rho0 = compute_pseudorange(prn,u->week,u->sow,u,sat);
        double enu[3]; ecef2enu(u,sat,enu);
        double el = enu_elevation_deg(enu); if(el<5.0) continue;
        c[m++] = (struct cand){prn,el,rho0};
    }
    for(int i=0;i<m-1;++i) for(int j=i+1;j<m;++j){
        if(c[j].elev>c[i].elev){
            struct cand t=c[i];c[i]=c[j];c[j]=t;
        }
    }
    *n = m<MAX_CH?m:MAX_CH;
    for(int i=0;i<*n;++i)
        channel_reset(&ch[i],c[i].prn,c[i].rho);
    return *n;
}


/* 依照使用者軌跡更新通道，使用目前與下一步的幾何資訊 */
static void update_channels_path(channel_t *ch,int n,
                                 const coord_t *u,
                                 const coord_t *u_next,
                                 const coord_t *u_next2,
                                 double gain,double target_cn0,int step_ms)
{
    double dt = step_ms * 0.001; /* seconds */

    for(int i=0;i<n;++i){
        int prn = ch[i].prn;

        double sat0[3];
        double rho0 = compute_pseudorange(prn, u->week, u->sow, u, sat0);

        double t1_abs = u->week*604800.0 + u->sow + dt;
        int week1 = (int)(t1_abs/604800.0);
        double sow1 = t1_abs - week1*604800.0;
        double sat1[3];
        double rho1 = compute_pseudorange(prn, week1, sow1, u_next, sat1);

        double rdot = (rho1 - rho0) / dt;
        channel_set_time(&ch[i], rho0, 0);
        double enu0[3]; ecef2enu(u, sat0, enu0);
        double el0 = enu_elevation_deg(enu0);
        update_channel_dynamics(&ch[i], rho0, rdot, el0, gain, target_cn0, n, step_ms);
        if(el0 < 5.0) ch[i].amp = 0.0;

        double t2_abs = t1_abs + dt;
        int week2 = (int)(t2_abs/604800.0);
        double sow2 = t2_abs - week2*604800.0;
        double rho2 = compute_pseudorange(prn, week2, sow2, u_next2, NULL);
        double rdot2 = (rho2 - rho1) / dt;
        double enu1[3]; ecef2enu(u_next, sat1, enu1);
        double el1 = enu_elevation_deg(enu1);
        double fd2 = -FCARRIER * rdot2 / CLIGHT;
        double cr2 = CHIPRATE * (1.0 - rdot2 / CLIGHT);
        double A2 = predict_next_amp(&ch[i], rho1, el1, gain, target_cn0, n, step_ms);
        if(el1 < 5.0) A2 = 0.0;
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

    /* 檢查 start 時間是否與星曆 toe 接近 */
    check_ephemeris_age(usr.week, usr.sow);

    g_t_tx = usr.week*604800.0 + usr.sow;

    channel_t ch[MAX_CH];
    int n_ch;
    select_channels(ch,&n_ch,&usr,cfg->single_prn,
                    cfg->meo_only);

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
    double start_bdt = usr.week*604800.0 + usr.sow;
    uint64_t sample_count = 0;
    for(uint64_t ms=0; ms<total_ms; ms+=STEP_MS)
    {
        g_t_tx = start_bdt + sample_count/fs;
        int week = (int)(g_t_tx/604800.0);
        double sow = g_t_tx - week*604800.0;

        double dt = STEP_MS * 0.001;
        coord_t usr_next={0}, usr2={0};
        if(cfg->path_type==0){
            usr.week = week;
            usr.sow  = sow;
            usr_next = usr;
            usr_next.sow = sow + dt;
            usr2 = usr;
            usr2.sow = sow + 2*dt;
            update_channels_path(ch,n_ch,&usr,&usr_next,&usr2,
                                 cfg->gain,g_target_cn0,STEP_MS);
        } else if(cfg->path_type==1){
            coord_t u0,u1,u2;
            interpolate_path(&path, ms/1000.0, &u0);
            interpolate_path(&path, (ms+STEP_MS)/1000.0, &u1);
            interpolate_path(&path, (ms+2*STEP_MS)/1000.0, &u2);
            xyz2llh(u0.xyz,&usr); usr.week=week; usr.sow=sow;
            xyz2llh(u1.xyz,&usr_next); usr_next.week=week; usr_next.sow=sow+dt;
            xyz2llh(u2.xyz,&usr2); usr2.week=week; usr2.sow=sow+2*dt;
            update_channels_path(ch,n_ch,&usr,&usr_next,&usr2,
                                 cfg->gain,g_target_cn0,STEP_MS);
        } else {
            double llh0[3], llh1[3], llh2[3];
            interpolate_path_llh(&path, ms/1000.0, llh0);
            interpolate_path_llh(&path, (ms+STEP_MS)/1000.0, llh1);
            interpolate_path_llh(&path, (ms+2*STEP_MS)/1000.0, llh2);
            llh2xyz(llh0,&usr); usr.week=week; usr.sow=sow;
            llh2xyz(llh1,&usr_next); usr_next.week=week; usr_next.sow=sow+dt;
            llh2xyz(llh2,&usr2); usr2.week=week; usr2.sow=sow+2*dt;
            update_channels_path(ch,n_ch,&usr,&usr_next,&usr2,
                                 cfg->gain,g_target_cn0,STEP_MS);
        }
        if(ms==0){
            for(int i=0;i<n_ch;++i){
                double rdot = -ch[i].fd*CLIGHT/FCARRIER;
                printf("[ch%02d] rdot %.2f fd %.2fHz\n", ch[i].prn, rdot, ch[i].fd);
            }
        }

        /* --- STEP_MS 次 1ms 取樣 --- */
        for(int step=0;step<STEP_MS;++step){
            g_t_tx = start_bdt + sample_count/fs;
            memset(accI,0,sizeof(accI)); memset(accQ,0,sizeof(accQ));

            /* 併行各通道 */
            #pragma omp parallel for
            for(int c=0;c<n_ch;++c){
                gen_samples_1ms(&ch[c],samp_per_ms,tmpI[c],tmpQ[c]);
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
            sample_count += samp_per_ms;
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

