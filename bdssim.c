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

#define FSAMP     8.184e6
#define SAMP_1MS  8184
#define STEP_MS   10          /* 幾何更新粒度 */

static const int geo[] ={1,2,3,4,5,59,60,61,62,63};
static int is_geo(int p){for(int i=0;i<10;i++) if(p==geo[i]) return 1; return 0;}

/*----------------------------------------------------*/
int select_channels(channel_t *ch,int *n,const coord_t*u)
{
    struct cand{int prn;double elev,rho,rdot;} c[63]; int m=0;
    for(int prn=1;prn<=63;++prn){
        if(is_geo(prn))continue;
        double sat[3],vel[3]; calc_sat_position_velocity(prn,u->week,u->sow,sat,vel);
        double enu[3]; ecef2enu(u,sat,enu);
        double el=enu_elevation_deg(enu); if(el<10)continue;
        double dx=sat[0]-u->xyz[0], dy=sat[1]-u->xyz[1], dz=sat[2]-u->xyz[2];
        double rho=hypot(hypot(dx,dy),dz);
        double rdot=(dx*vel[0]+dy*vel[1]+dz*vel[2])/rho;
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
    coord_t usr={0}; llh2xyz(cfg->llh,&usr);
    if(utc_to_bdt(cfg->time_start,&usr.week,&usr.sow)!=0){fputs("UTC format err\n",stderr);return;}

    navbits_init(usr.week, usr.sow);

    channel_t ch[MAX_CH]; int n_ch; select_channels(ch,&n_ch,&usr);

    /* 首次幾何 – 初始化振幅/NCO */
    for(int i=0;i<n_ch;++i){
        double sat[3],vel[3]; calc_sat_position_velocity(ch[i].prn,usr.week,usr.sow,sat,vel);
        double dx=sat[0]-usr.xyz[0],dy=sat[1]-usr.xyz[1],dz=sat[2]-usr.xyz[2];
        double rho=hypot(hypot(dx,dy),dz);
        double rdot=(dx*vel[0]+dy*vel[1]+dz*vel[2])/rho;
        update_channel_dynamics(&ch[i],rho,rdot,n_ch);
    }

    FILE *fp=fopen("beidou_b1i.bin","wb"); if(!fp){perror("bin");return;}
    int16_t tmpI[MAX_CH][SAMP_1MS],tmpQ[MAX_CH][SAMP_1MS];
    int32_t accI[SAMP_1MS],accQ[SAMP_1MS],outI[SAMP_1MS],outQ[SAMP_1MS];

    const uint64_t total_ms=(uint64_t)cfg->duration*1000;
    for(uint64_t ms=0; ms<total_ms; ms+=STEP_MS)
    {
        /* --- 幾何重算每 STEP_MS --- */
        double sow = usr.sow + ms*0.001;
        for(int i=0;i<n_ch;++i){
            double sat[3],vel[3]; calc_sat_position_velocity(ch[i].prn,usr.week,sow,sat,vel);
            double dx=sat[0]-usr.xyz[0],dy=sat[1]-usr.xyz[1],dz=sat[2]-usr.xyz[2];
            double rho=hypot(hypot(dx,dy),dz);
            double rdot=(dx*vel[0]+dy*vel[1]+dz*vel[2])/rho;
            update_channel_dynamics(&ch[i],rho,rdot,n_ch);
        }

        /* --- STEP_MS 次 1ms 取樣 --- */
        for(int step=0;step<STEP_MS;++step){
            memset(accI,0,sizeof(accI)); memset(accQ,0,sizeof(accQ));

            /* 併行各通道 */
            #pragma omp parallel for
            for(int c=0;c<n_ch;++c)
                gen_samples_1ms(&ch[c],tmpI[c],tmpQ[c]);

            /* 歸併 */
            for(int c=0;c<n_ch;++c)
                for(int k=0;k<SAMP_1MS;++k){
                    accI[k]+=tmpI[c][k];
                    accQ[k]+=tmpQ[c][k];
                }
            for(int k=0;k<SAMP_1MS;++k){
                int32_t i=accI[k], q=accQ[k];
                if(i>32767)i=32767; else if(i<-32768)i=-32768;
                if(q>32767)q=32767; else if(q<-32768)q=-32768;
                outI[k]=i; outQ[k]=q;
            }
            fwrite(outI,2,SAMP_1MS,fp);
            fwrite(outQ,2,SAMP_1MS,fp);
        }
        /* 進度 0.01 s = STEP_MS ms */
        printf("\r進度: %.2f / %.2f 秒",(ms+STEP_MS)/1000.0,total_ms/1000.0);
        fflush(stdout);
    }
    puts(""); fclose(fp);
    puts("[bdssim] 完成多星基帶輸出 beidou_b1i.bin");
}

