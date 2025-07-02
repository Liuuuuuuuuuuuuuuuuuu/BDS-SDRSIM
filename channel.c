/* channel.c : 單通道 1 ms B1I I/Q 產生 */
#include <string.h>
#include <stdint.h>
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

/* ---------- 振幅 ---------- */
double calc_amp(double rho,int n_ch,double gain){
    double p = -130.0 - 20.0*log10(rho/2.0e7);
    double a = pow(10.0,(p-DBM_REF)/20.0);    /* relative field amplitude */
    /* Scale to 16‑bit output range. Each channel shares the range */
    a *= 16384.0 / n_ch;
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
void channel_reset(channel_t *c,int prn){
    memset(c,0,sizeof(*c));
    c->prn   = prn;
    c->sf_id = 1;                     /* start from subframe 1 */
    if(!ca_ready){
        for(int p=1;p<=63;++p)
            for(int i=0;i<CODE_LEN;++i)
                ca_wave[p][i] = prn_code[p][i]?+1:-1;
        ca_ready=1;
    }
}
/* 幾何→計算振幅 / 初始多普勒 */
void update_channel_dynamics(channel_t *c,double rho,double rdot,int n_ch,double gain){
    c->amp = calc_amp(rho,n_ch,gain);
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
        get_subframe_bits(c->prn,c->sf_id,week,sow,c->nav_bits);

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

