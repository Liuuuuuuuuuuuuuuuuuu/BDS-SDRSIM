/* channel.c : 單通道 1 ms B1I I/Q 產生 */
#include <string.h>
#include <stdint.h>
#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include "channel.h"
#include "bdssim.h"
#include "navbits.h"
#include "orbits.h"

#define PI2        6.2831853071795864769
static double fs = FS_OUTPUT_HZ;           /* default sample rate */

/* default target CN0, overridable via CLI */
double g_target_cn0 = 42.0;

/* (舊的接收天線圖與 multipath 模型已移除) */

/* --------------------------- 32k sin LUT ------------------------------*/
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

/* --------------------------- 載波與振幅計算 ------------------------------*/

/* PRNs 1–5 and 59–63 are GEO satellites and are skipped during
 * channel selection. */
int is_geo_prn(int prn)
{
    if(prn < 1 || prn >= MAX_SAT) return 0;
    for(unsigned i=0;i<sizeof(geo_prn)/sizeof(geo_prn[0]);++i)
        if(prn == geo_prn[i]) return 1;
    return 0;
}

int is_igso_prn(int prn)
{
    if(prn < 1 || prn >= MAX_SAT) return 0;
    const ephemeris_t *ep = &eph[prn];
    if(ep->prn == 0) return 0;
    return ep->sqrtA > 6000.0 && fabs(ep->i0) >= 0.2; /* IGSO */
}

int is_meo_prn(int prn)
{
    if(prn < 1 || prn >= MAX_SAT) return 0;
    const ephemeris_t *ep = &eph[prn];
    if(ep->prn == 0) return 0;
    return ep->sqrtA <= 6000.0; /* ~27800 km orbit */
}

/* --------------------------- 振幅模型 ------------------------------*/
static double ant_pat[37];
__attribute__((constructor))
static void init_ant_pat(void){
    for(int i=0;i<37;++i)
        ant_pat[i] = pow(10.0, -ant_pat_db[i]/20.0);
}

static inline double amp_from_geom(double rho,double elev_deg,double gain,
                                   double cn0_dBHz,int n_visible)
{
    double base = pow(10.0,(cn0_dBHz-45.0)/20.0) * 16384.0;
    double path_loss = 20200000.0 / rho; /* 模型常數 */
    int ibs = (int)((90.0 - elev_deg)/5.0);
    if(ibs < 0) ibs = 0; else if(ibs > 36) ibs = 36;
    double ant = ant_pat[ibs];
    if (n_visible < 1) n_visible = 1;
    double A = base * path_loss * ant * gain;
    return (A / sqrt((double)n_visible)) * HEADROOM_RATIO;
}

static inline double orbit_gain_amp(int prn)
{
    double dB = is_meo_prn(prn) ? GAIN_MEO_DB : GAIN_IGSO_DB;
    return pow(10.0, dB/20.0);
}

/* --------------------------- 指數平滑振幅，避免 AM 旁帶 ------------------------------*/
static inline double smooth_amp(double A_prev, double A_new, double dt_ms)
{
    double alpha = 1.0 - exp(-dt_ms / AMP_SMOOTH_TC_MS);
    if (A_prev == 0.0)
        return A_new;                  /* avoid long ramp at start */
    return A_prev + alpha * (A_new - A_prev);
}

/* 預測下一步的平滑振幅 */
double predict_next_amp(const channel_t *c,double rho_next,double elev_deg_next,
                        double gain,double target_cn0,int n_visible,double dt_ms)
{
    double A_new = amp_from_geom(rho_next,elev_deg_next,gain,
                                 target_cn0,n_visible);
    A_new *= orbit_gain_amp(c->prn);
    return smooth_amp(c->amp, A_new, dt_ms);
}

/* --------------------------- CA cache ------------------------------*/
static int16_t ca_wave[64][CODE_LEN];
static int     ca_ready=0;

/* BeiDou D1 Neumann-Hoffman 20-bit code (0=+1, 1=-1) */

/* --------------------------- Channel helpers ------------------------------*/
static void load_ca_once(void)
{
    if(ca_ready) return;
    for(int p=1;p<=prn_max;++p)
        for(int i=0;i<CODE_LEN;++i)
            ca_wave[p][i] = prn_code[p][i]?+1:-1;
    ca_ready = 1;
}

void channel_set_time(channel_t *c,double rho)
{
    /* Satellite signal left earlier by propagation delay. */
    double tx_time = g_t_tx - rho/CLIGHT; /* seconds since BDT epoch */
    int    week_tx = (int)(tx_time/604800.0);
    double sow_tx  = tx_time - week_tx*604800.0;

    /* Compute sub-ms fractional offset to align PRN/NH/data boundaries. */
    double sow_ms  = sow_tx * 1000.0;
    double frac_ms = sow_ms - floor(sow_ms);
    c->code_phase  = frac_ms * CODE_LEN;  /* 1 ms = 2046 chips */

    /* Navigation subframe alignment (6 s). */
    double sf_start = floor(sow_tx/6.0)*6.0;
    if(sf_start != c->nav_sf_start){
        c->sf_id = ((int)(sf_start/6.0))%5 + 1;  /* 1..5 */
        get_subframe_bits(c->prn,c->sf_id,week_tx,sf_start,6.0,c->nav_bits);
        c->nav_sf_start = sf_start;
    } else {
        c->sf_id = ((int)(sf_start/6.0))%5 + 1;  /* 1..5 */
    }
    int ms = (int)((sow_tx - sf_start)*1000.0);
    if(ms < 0)      ms = 0;
    else if(ms >= 6000) ms = 5999;
    c->bit_ptr = ms/20;
    c->ms_count = ms%20;

    if(ms == 0 && frac_ms < 1e-9){
        /* At subframe boundary, force perfect alignment of PRN/NH/NAV. */
        c->code_phase = 0.0;
        c->ms_count   = 0;
        c->bit_ptr    = 0;
    }

}

void channel_reset(channel_t *c,int prn,double rho){
    memset(c,0,sizeof(*c));
    c->prn   = prn;
    load_ca_once();

    /* Randomise starting carrier phase so I/Q averages are balanced. */
    c->carr_phase = ((double)rand()/(double)RAND_MAX)*PI2;

    channel_set_time(c,rho);
}
/* 幾何→計算振幅 / 初始多普勒 */
void update_channel_dynamics(channel_t *c,double rho,double rdot,double elev_deg,
                             double gain,double target_cn0,int n_visible,
                             double dt_ms)
{
    double A_new = amp_from_geom(rho,elev_deg,gain,target_cn0,n_visible);
    A_new *= orbit_gain_amp(c->prn);
    c->amp = smooth_amp(c->amp, A_new, dt_ms);
    c->amp_dot = 0.0;
    c->elev_deg = elev_deg;
    c->fd  = -F_B1I*rdot/CLIGHT;                  /* Doppler (Hz) */
    c->code_rate = CHIPRATE*(1.0 - rdot/CLIGHT);  /* Code frequency (Hz) */
}

void channel_set_fs(double sample_rate)
{
    fs = sample_rate;
}

/* 產生 1 ms */
void gen_samples_1ms(channel_t *c,int samp_per_ms,int16_t*I,int16_t*Q)
{
    double fd = c->fd;
    double code_rate = c->code_rate;
    double phase = c->carr_phase;
    double code_phase = c->code_phase;
    double dfd = c->fd_dot / fs;
    double dcode_rate = c->code_rate_dot / fs;
    double amp = c->amp;
    double damp = c->amp_dot / fs;

    for(int n=0;n<samp_per_ms;++n){
        int chip = (int)code_phase;            /* 0..2045 */
        int16_t ca = ca_wave[c->prn][chip];
        uint8_t nh = nh20_bits[c->ms_count];
        int16_t nb = (c->nav_bits[c->bit_ptr]^nh)?-1:+1;
        float co,si; fast_sincos(phase,&co,&si);
        float s = amp*ca*nb;
        I[n]=(int16_t)lrintf(s*co);
        Q[n]=(int16_t)lrintf(s*si);

        phase += PI2*(fd + 0.5*dfd)/fs;
        fd += dfd;
        if(phase>=PI2)      phase-=PI2;
        else if(phase<0.0)  phase+=PI2;
        code_phase += (code_rate + 0.5*dcode_rate)/fs;
        code_rate += dcode_rate;
        amp += damp;
        if(code_phase>=CODE_LEN){
            code_phase-=CODE_LEN;
            if(++c->ms_count==20){
                c->ms_count=0;
                if(++c->bit_ptr==300){
                    c->bit_ptr=0;
                    c->sf_id=c->sf_id%5+1;      /* 1..5 */
                }
            }
        }
    }
    c->fd = fd;
    c->code_rate = code_rate;
    c->carr_phase = phase;
    c->code_phase = code_phase;
    c->amp = amp;
}
/* ---------------------------  End  ------------------------------*/
