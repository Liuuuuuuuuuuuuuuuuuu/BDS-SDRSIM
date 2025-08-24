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
#define FSAMP_DEF FS_OUTPUT_HZ    /* default 5.2 MHz for 16-bit I/Q output */
#define FSAMP_BYTE 25.0e6    /* 25 MHz when --byte is used */
#define IF_BYTE    1000.0    /* 1 kHz IF for --byte output */
#define F_B1I      1561.098e6      /* B1I carrier */
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

/* Compute corrected pseudorange and range rate */
static void compute_range_rate(int prn,int week,double sow,
                               const coord_t *u,const double uvel[3],
                               double sat_rx[3], double vsat_rx[3],
                               double *rho,double *rdot,double los[3])
{
    /* Receiver time in seconds since BDT epoch */
    double t_rx = week*604800.0 + sow;

    /* Satellite clock at receiver time for initial guess */
    double clk_rx[2], pos_dummy[3], vel_dummy[3];
    calc_sat_position_velocity(prn, week, sow, pos_dummy, vel_dummy, clk_rx);

    /* Initial transmit time estimate */
    double t_tx = t_rx - clk_rx[0];

    /* Iterate once with geometric range */
    int week_sv = (int)(t_tx/604800.0);
    double sow_sv = t_tx - week_sv*604800.0;
    double r_tx[3], v_tx[3], clk_sv[2];
    calc_sat_position_velocity(prn, week_sv, sow_sv, r_tx, v_tx, clk_sv);
    double dx = r_tx[0] - u->xyz[0];
    double dy = r_tx[1] - u->xyz[1];
    double dz = r_tx[2] - u->xyz[2];
    double range = hypot(hypot(dx,dy),dz);
    double tau = range / CLIGHT;

    /* Refined transmit time */
    t_tx = t_rx - tau - clk_sv[0];
    week_sv = (int)(t_tx/604800.0);
    sow_sv  = t_tx - week_sv*604800.0;
    calc_sat_position_velocity(prn, week_sv, sow_sv, r_tx, v_tx, clk_sv);

    /* Time difference between tx and rx */
    tau = t_rx - t_tx;

    /* Rotate satellite state to receiver time */
    double A = OMEGA_E * tau;
    double cA = cos(A), sA = sin(A);
    double x_rx =  cA*r_tx[0] + sA*r_tx[1];
    double y_rx = -sA*r_tx[0] + cA*r_tx[1];
    double z_rx =  r_tx[2];

    double vx_rot =  cA*v_tx[0] + sA*v_tx[1];
    double vy_rot = -sA*v_tx[0] + cA*v_tx[1];
    double vz_rot =  v_tx[2];

    double vx_rx = vx_rot + (sA*OMEGA_E)*r_tx[0] + (-cA*OMEGA_E)*r_tx[1];
    double vy_rx = vy_rot + (cA*OMEGA_E)*r_tx[0] + ( sA*OMEGA_E)*r_tx[1];
    double vz_rx = vz_rot;

    if(sat_rx){ sat_rx[0]=x_rx; sat_rx[1]=y_rx; sat_rx[2]=z_rx; }
    if(vsat_rx){ vsat_rx[0]=vx_rx; vsat_rx[1]=vy_rx; vsat_rx[2]=vz_rx; }

    dx = x_rx - u->xyz[0];
    dy = y_rx - u->xyz[1];
    dz = z_rx - u->xyz[2];
    range = hypot(hypot(dx,dy),dz);

    double ex = dx / range;
    double ey = dy / range;
    double ez = dz / range;
    if(los){ los[0]=ex; los[1]=ey; los[2]=ez; }

    if (rho) *rho = range - CLIGHT * clk_sv[0];

    if (rdot) {
        double vrx = uvel ? uvel[0] : 0.0;
        double vry = uvel ? uvel[1] : 0.0;
        double vrz = uvel ? uvel[2] : 0.0;
        double dvx = vx_rx - vrx;
        double dvy = vy_rx - vry;
        double dvz = vz_rx - vrz;
        *rdot = ex*dvx + ey*dvy + ez*dvz;  /* m/s */
    }
}

/*----------------------------------------------------*/
int select_channels(channel_t *ch,int *n,const coord_t*u,
                    int single_prn,bool meo_only)
{
    struct cand{int prn;double elev,rho,rdot;} c[63];
    int m=0;
    for(int prn=1;prn<=prn_max;++prn){
        if(single_prn>0 && prn!=single_prn) continue;
        if(is_geo_prn(prn)) continue;
        if(meo_only && !is_meo_prn(prn)) continue;
        double sat[3], rho, rdot;
        compute_range_rate(prn,u->week,u->sow,u,NULL,sat,NULL,&rho,&rdot,NULL);
        double enu[3]; ecef_to_enu(u,sat,enu);
        double el=enu_elevation_deg(enu); if(el<5.0)continue;
        c[m++] = (struct cand){prn,el,rho,rdot};
    }
    /* sort by elevation (desc) */
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
                                 const coord_t *u,const double uvel[3],
                                 const coord_t *u_next,const double uvel_next[3],
                                 double gain,double target_cn0,int step_ms)
{
    double dt = step_ms * 0.001; /* seconds */

    double t_abs = u->week*604800.0 + u->sow;
    static int last_sec = -1;
    int cur_sec = (int)t_abs;
    int print_now = 0;
    if(cur_sec != last_sec){
        last_sec = cur_sec;
        print_now = 1;
    }

    for(int i=0;i<n;++i){
        int prn = ch[i].prn;
        double sat[3], rho, rdot;
        double vsat[3], los[3];
        compute_range_rate(prn,u->week,u->sow,u,uvel,sat,vsat,&rho,&rdot,los);
        channel_set_time(&ch[i], rho, 0);
        double enu[3]; ecef_to_enu(u,sat,enu);
        double el = enu_elevation_deg(enu);
        update_channel_dynamics(&ch[i],rho,rdot,el,gain,target_cn0,n,step_ms);
        if(print_now){
            double rvx = uvel ? uvel[0] : 0.0;
            double rvy = uvel ? uvel[1] : 0.0;
            double rvz = uvel ? uvel[2] : 0.0;
            double rel_vx = vsat[0] - rvx;
            double rel_vy = vsat[1] - rvy;
            double rel_vz = vsat[2] - rvz;
            double rdot_pred = los[0]*rel_vx + los[1]*rel_vy + los[2]*rel_vz;
            double rdot_fd = CLIGHT*ch[i].fd/F_B1I;
            double delta = rdot_pred - rdot_fd;
            printf("[delta] PRN%02d %.3f %.3f %.3f\n", prn, rdot_pred, rdot_fd, delta);
        }
        if(el < 5.0) ch[i].amp = 0.0;     /* below horizon → mute */

        /* predict next dynamics using next position/velocity */
        double t_abs = u->week*604800.0 + u->sow + dt;
        int week2 = (int)(t_abs/604800.0);
        double sow2 = t_abs - week2*604800.0;
        double sat2[3], rho2, rdot2;
        compute_range_rate(prn,week2,sow2,u_next,uvel_next,sat2,NULL,&rho2,&rdot2,NULL);
        double enu2[3]; ecef_to_enu(u_next,sat2,enu2);
        double el2 = enu_elevation_deg(enu2);
        double fd2 = F_B1I*rdot2/CLIGHT;
        double cr2 = CHIPRATE*(1.0 - rdot2/CLIGHT);
        double A2 = predict_next_amp(&ch[i], rho2, el2, gain, target_cn0, n, step_ms);
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
    if(cfg->path_type==0){
        double llh_rad[3]={cfg->llh[0]*M_PI/180.0,
                           cfg->llh[1]*M_PI/180.0,
                           cfg->llh[2]};
        lla_to_ecef(llh_rad,&usr);
    } else { interpolate_path(&path,0.0,&usr); ecef_to_lla(usr.xyz,&usr); }
    if(utc_to_bdt(cfg->time_start,&usr.week,&usr.sow)!=0){
        fputs("UTC format err\n",stderr);
        free_path(&path);
        return;
    }

    coord_t ref_llh=usr;              /* 保存經緯度作旋轉基準 */
    static_user_at(usr.week,usr.sow,&ref_llh,&usr,NULL);

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

    double c_if = 1.0, s_if = 0.0;        /* mixer state for byte output */
    double cos_d = 1.0, sin_d = 0.0;      /* per-sample rotation */
    if(cfg->byte_output){
        double w_if = 2*M_PI*IF_BYTE/fs;
        cos_d = cos(w_if);
        sin_d = sin(w_if);
    }

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
        double uvel[3], uvel_next[3];
        coord_t usr_next={0};
        if(cfg->path_type==0){
            /* Static user: include Earth-rotation velocity */
            static_user_at(week, sow, &ref_llh, &usr, uvel);
            static_user_at(week, sow + dt, &ref_llh, &usr_next, uvel_next);
            update_channels_path(ch,n_ch,&usr,uvel,&usr_next,uvel_next,
                                 cfg->gain,g_target_cn0,STEP_MS);
        } else if(cfg->path_type==1){
            coord_t u0,u1,u2;
            interpolate_path(&path, ms/1000.0, &u0);
            interpolate_path(&path, (ms+STEP_MS)/1000.0, &u1);
            interpolate_path(&path, (ms+2*STEP_MS)/1000.0, &u2);
            ecef_to_lla(u0.xyz,&usr); usr.week=week; usr.sow=sow;
            ecef_to_lla(u1.xyz,&usr_next);
            ecef_to_lla(u2.xyz,&u2);
            uvel[0]=(u1.xyz[0]-u0.xyz[0])/dt;
            uvel[1]=(u1.xyz[1]-u0.xyz[1])/dt;
            uvel[2]=(u1.xyz[2]-u0.xyz[2])/dt;
            uvel_next[0]=(u2.xyz[0]-u1.xyz[0])/dt;
            uvel_next[1]=(u2.xyz[1]-u1.xyz[1])/dt;
            uvel_next[2]=(u2.xyz[2]-u1.xyz[2])/dt;
            update_channels_path(ch,n_ch,&usr,uvel,&usr_next,uvel_next,
                                 cfg->gain,g_target_cn0,STEP_MS);
        } else {
            double llh0[3], llh1[3], llh2[3];
            interpolate_path_llh(&path, ms/1000.0, llh0);
            interpolate_path_llh(&path, (ms+STEP_MS)/1000.0, llh1);
            interpolate_path_llh(&path, (ms+2*STEP_MS)/1000.0, llh2);
            coord_t ref0={0}, ref1={0}, ref2={0};
            ref0.llh[0]=llh0[0]; ref0.llh[1]=llh0[1]; ref0.llh[2]=llh0[2];
            ref1.llh[0]=llh1[0]; ref1.llh[1]=llh1[1]; ref1.llh[2]=llh1[2];
            ref2.llh[0]=llh2[0]; ref2.llh[1]=llh2[1]; ref2.llh[2]=llh2[2];
            static_user_at(week,sow,&ref0,&usr,NULL);
            static_user_at(week,sow+dt,&ref1,&usr_next,NULL);
            coord_t usr2; static_user_at(week,sow+2*dt,&ref2,&usr2,NULL);
            uvel[0]=(usr_next.xyz[0]-usr.xyz[0])/dt;
            uvel[1]=(usr_next.xyz[1]-usr.xyz[1])/dt;
            uvel[2]=(usr_next.xyz[2]-usr.xyz[2])/dt;
            uvel_next[0]=(usr2.xyz[0]-usr_next.xyz[0])/dt;
            uvel_next[1]=(usr2.xyz[1]-usr_next.xyz[1])/dt;
            uvel_next[2]=(usr2.xyz[2]-usr_next.xyz[2])/dt;
            update_channels_path(ch,n_ch,&usr,uvel,&usr_next,uvel_next,
                                 cfg->gain,g_target_cn0,STEP_MS);
        }
        if(ms==0){
            for(int i=0;i<n_ch;++i){
                double rdot = ch[i].fd*CLIGHT/F_B1I;
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
            int8_t  i8[samp_per_ms];            /* byte output holds real IF signal */
            double c = c_if, s = s_if;
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
                    double mixed = (double)i*c - (double)q*s;
                    int32_t m = (int32_t)llround(mixed);
                    if (m>127) m=127; else if (m<-128) m=-128;
                    i8[k] = (int8_t)m;
                    sumI  += i8[k];
                    sumI2 += (double)i8[k]*i8[k];
                    double tc = c*cos_d - s*sin_d;
                    double ts = s*cos_d + c*sin_d;
                    c = tc; s = ts;
                }
            }
            if(!fp){ c_if = c; s_if = s; }
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

