/* channel.c : 單通道 1 ms B1I I/Q 產生 */
#include <string.h>
#include <stdint.h>
#include <stdlib.h>
#include <math.h>
#include "channel.h"
#include "globals.h"               /* prn_code */

#define PI2        6.2831853071795864769
#define FCARRIER   1561.098e6      /* B1I */
#define CHIPRATE   2.046e6
static double fs = 4.092e6;
#define DBM_REF   (-130.0)         /* ±1.0 → –130 dBm */

/* ---------- 32k sin LUT ---------- */
#define LUTBITS   15
#define LUTSIZE   (1u<<LUTBITS)
static float sin_lut[LUTSIZE];
__attribute__((constructor))
static void init_sin(void){
    for(unsigned i=0;i<LUTSIZE;++i) sin_lut[i]=sinf((PI2*i)/LUTSIZE);
}
static inline void fast_sincos(double ph,float*co,float*si){
    ph -= floor(ph/PI2)*PI2;
    double idx = ph*LUTSIZE/PI2;
    uint32_t i = (uint32_t)idx & (LUTSIZE-1);
    uint32_t i2=(i+1)&(LUTSIZE-1);
    float f = (float)(idx-(double)i);
    *si = sin_lut[i] + f*(sin_lut[i2]-sin_lut[i]);
    uint32_t j=(i+(LUTSIZE>>2))&(LUTSIZE-1);
    uint32_t j2=(j+1)&(LUTSIZE-1);
    *co = sin_lut[j] + f*(sin_lut[j2]-sin_lut[j]);
}

/* ---------- 振幅與功率模型 ---------- */

/*
 * 依衛星軌道類型給予不同的 EIRP。此處僅以簡化常數代表：
 *   GEO  衛星約 52 dBm
 *   IGSO 衛星約 53 dBm
 *   MEO  衛星約 55 dBm
 */

static const int geo_prn[] = {1,2,3,4,5,59,60,61,62,63};

static int is_geo_prn(int prn)
{
    for(size_t i=0;i<sizeof(geo_prn)/sizeof(geo_prn[0]);++i)
        if(prn==geo_prn[i]) return 1;
    return 0;
}

static double sat_eirp_dbm(int prn)
{
    if(is_geo_prn(prn))      return 52.0; /* GEO */
    else if(prn>=6 && prn<=10) return 53.0; /* IGSO (粗略) */
    else                      return 55.0; /* MEO  */
}

/* 大氣衰減常數 (dB) */
#define ATM_LOSS_DB    2.0

double calc_amp(int prn,double rho,double gain,double target_cn0)
{
    /* dB path-loss + 衛星 Tx-power → 線性功率，再正規化到 ±16384 */
    double lambda    = 299792458.0/FCARRIER;
    double path_loss = 20.0*log10(4.0*M_PI*rho/lambda) + ATM_LOSS_DB;
    double p_dbm     = sat_eirp_dbm(prn) - path_loss;
    double cn0_dbhz  = p_dbm - DBM_REF + 10.0*log10(fs);
    double diff_db   = target_cn0 - cn0_dbhz;
    double a         = pow(10.0,diff_db/20.0) * 16384.0;
    return gain * a;
}

/* ---------- CA cache ---------- */
static int16_t ca_wave[64][CODE_LEN];
static int     ca_ready=0;

/* BeiDou D1 Neumann-Hoffman 20-bit code (0=+1, 1=-1) */
static const uint8_t nh20_bits[20]={
    0,0,0,0,0,1,0,0,1,1,0,1,0,1,0,0,1,1,1,0
};


/* ---------- Channel helpers ---------- */
static void load_ca_once(void)
{
    if(ca_ready) return;
    for(int p=1;p<=63;++p)
        for(int i=0;i<CODE_LEN;++i)
            ca_wave[p][i] = prn_code[p][i]?+1:-1;
    ca_ready = 1;
}

void channel_set_time(channel_t *c,int week,double sow)
{
    double sf_start = floor(sow/6.0)*6.0;
    c->sf_id = ((int)(sf_start/6.0))%5 + 1;
    int ms = (int)llround((sow - sf_start)*1000.0);
    if(ms < 0)      ms = 0;
    else if(ms >= 6000) ms = 5999;
    c->bit_ptr = ms/20;
    c->ms_count = ms%20;
    get_subframe_bits(c->prn,c->sf_id,week,sf_start,6.0,c->nav_bits);

    double sf_start_d2 = floor(sow/0.6)*0.6;
    c->sf_id_d2 = ((int)(sf_start_d2/0.6))%5 + 1;
    int ms2 = (int)llround((sow - sf_start_d2)*1000.0);
    if(ms2 < 0)      ms2 = 0;
    else if(ms2 >= 600) ms2 = 599;
    c->bit_ptr_d2 = ms2/2;
    c->ms_count_d2 = ms2%2;
    get_subframe_bits(c->prn,c->sf_id_d2,week,sf_start_d2,0.6,c->nav_bits_d2);
}

void channel_reset(channel_t *c,int prn,int week,double sow){
    memset(c,0,sizeof(*c));
    c->prn   = prn;
    load_ca_once();

    /* Randomise starting carrier and code phase so I/Q averages
       are well balanced even for short captures. */
    c->carr_phase = ((double)rand()/(double)RAND_MAX)*PI2;
    c->code_phase = ((double)rand()/(double)RAND_MAX)*CODE_LEN;

    channel_set_time(c,week,sow);
}
/* 幾何→計算振幅 / 初始多普勒 */
void update_channel_dynamics(channel_t *c,double rho,double rdot,double gain,double target_cn0){
    c->amp = calc_amp(c->prn,rho,gain,target_cn0);
    c->fd  = -FCARRIER*rdot/299792458.0;               /* Doppler (Hz) */
    /*
     * Positive range rate (rdot) means the satellite is moving away
     * from the user, resulting in a lower received chipping rate.
     * The correct relationship is therefore (1 - rdot/c).
     */
    c->code_rate = CHIPRATE*(1.0 - rdot/299792458.0);  /* Code frequency (Hz) */
}

void channel_set_fs(double sample_rate)
{
    fs = sample_rate;
}

/* ---------- 產生 1 ms ---------- */
void gen_samples_1ms(channel_t *c,int week,double sow,
                     int samp_per_ms,int16_t*I,int16_t*Q)
{
    if(c->bit_ptr==0 && c->ms_count==0)
        get_subframe_bits(c->prn,c->sf_id,week,sow,6.0,c->nav_bits);

    /* Baseband output – only apply Doppler frequency */
    const double dphi = PI2*c->fd/fs;
    const double dcode = c->code_rate/fs;     /* chips per sample */
    double code_phase = c->code_phase;
    double phase = c->carr_phase;

    for(int n=0;n<samp_per_ms;++n){
        int chip = (int)code_phase;            /* 0..2045 */
        int16_t ca = ca_wave[c->prn][chip];
        uint8_t nh = nh20_bits[c->ms_count];
        int16_t nb = (c->nav_bits[c->bit_ptr]^nh)?+1:-1;
        float co,si; fast_sincos(phase,&co,&si);
        float s = c->amp*ca*nb;
        I[n]=(int16_t)lrintf(s*co);
        Q[n]=(int16_t)lrintf(s*si);

        /* NCO */
        phase += dphi;
        if(phase>=PI2)      phase-=PI2;
        else if(phase<0.0)  phase+=PI2;
        code_phase += dcode;
        if(code_phase>=CODE_LEN){
            code_phase-=CODE_LEN;
            if(++c->ms_count==20){
                c->ms_count=0;
                if(++c->bit_ptr==300){
                    c->bit_ptr=0;
                    c->sf_id=c->sf_id%5+1;
                }
            }
        }
    }
    c->carr_phase = phase;
    c->code_phase = code_phase;
}

/*
 * ======== D2 (500 bps) support ========
 *  - 資料速率 500 bps (2 ms per bit)
 *  - 不使用二次 Neumann-Hoffman 編碼
 */
void gen_samples_1ms_d2(channel_t *c, int week, double sow,
                               int samp_per_ms, int16_t *I, int16_t *Q)
{
    if(c->bit_ptr_d2==0 && c->ms_count_d2==0)
        get_subframe_bits(c->prn,c->sf_id_d2,week,sow,0.6,c->nav_bits_d2);

    const double dphi = PI2*c->fd/fs;
    const double dcode = c->code_rate/fs;
    double code_phase = c->code_phase;
    double phase = c->carr_phase;

    for(int n=0;n<samp_per_ms;++n){
        int chip = (int)code_phase;
        int16_t ca = ca_wave[c->prn][chip];
        int16_t nb = c->nav_bits_d2[c->bit_ptr_d2] ? -1:+1;
        float co,si; fast_sincos(phase,&co,&si);
        float s = c->amp*ca*nb;
        I[n]=(int16_t)lrintf(s*co);
        Q[n]=(int16_t)lrintf(s*si);

        phase += dphi;
        if(phase>=PI2)      phase-=PI2;
        else if(phase<0.0)  phase+=PI2;
        code_phase += dcode;
        if(code_phase>=CODE_LEN){
            code_phase-=CODE_LEN;
        }
    }
    c->carr_phase = phase;
    c->code_phase = code_phase;

    if(++c->ms_count_d2==2){
        c->ms_count_d2=0;
        if(++c->bit_ptr_d2==300){
            c->bit_ptr_d2=0;
            c->sf_id_d2 = c->sf_id_d2%5 + 1;
        }
    }
}

