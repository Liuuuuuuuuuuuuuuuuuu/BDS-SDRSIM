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
static double fs = FS_OUTPUT_HZ;           /* default sample rate */

/* default target CN0, overridable via CLI */
double g_target_cn0 = 42.0;

/* (舊的接收天線圖與 multipath 模型已移除) */

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

/* ---------- 載波與振幅計算 ---------- */

static const int geo_prn[] = {1,2,3,4,5,59,60,61,62,63};

static int is_geo_prn(int prn); /* forward */

/* classify IGSO/MEO via sqrtA (semi-major axis).  BDS IGSO shares
 * the GEO semi-major axis (~42164 km, sqrtA around 6493), whereas
 * MEO satellites use a smaller semi-major axis (~27800 km, sqrtA
 * around 5282).  This allows for more robust identification even if
 * PRN assignments change. */
int is_igso_prn(int prn)
{
    if(is_geo_prn(prn))
        return 0;
    if(prn<1 || prn>=MAX_SAT) return 0;
    const ephemeris_t *ep=&eph[prn];
    if(ep->prn==0) return 0;
    return ep->sqrtA > 6000.0;  /* ~42164 km orbit */
}

int is_meo_prn(int prn)
{
    if(is_geo_prn(prn))
        return 0;
    if(prn<1 || prn>=MAX_SAT) return 0;
    const ephemeris_t *ep=&eph[prn];
    if(ep->prn==0) return 0;
    return ep->sqrtA <= 6000.0; /* ~27800 km orbit */
}

static int is_geo_prn(int prn)
{
    for(size_t i=0;i<sizeof(geo_prn)/sizeof(geo_prn[0]);++i)
        if(prn==geo_prn[i]) return 1;
    return 0;
}

int is_d2_prn(int prn)
{
    return is_geo_prn(prn);
}

/* ---- gps-sdr-sim 風格的振幅模型 ---- */
static const double ant_pat_db[37] = {
     0.00,  0.00,  0.22,  0.44,  0.67,  1.11,  1.56,  2.00,  2.44,  2.89,
     3.56,  4.22,  4.89,  5.56,  6.22,  6.89,  7.56,  8.22,  8.89,  9.78,
    10.67, 11.56, 12.44, 13.33, 14.44, 15.56, 16.67, 17.78, 18.89, 20.00,
    21.33, 22.67, 24.00, 25.56, 27.33, 29.33, 31.56
};
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
    double path_loss = 20200000.0 / rho; /* gps-sdr-sim constant */
    int ibs = (int)((90.0 - elev_deg)/5.0);
    if(ibs < 0) ibs = 0; else if(ibs > 36) ibs = 36;
    double ant = ant_pat[ibs];
    if (n_visible < 1) n_visible = 1;
    double A = base * path_loss * ant * gain;
    return (A / sqrt((double)n_visible)) * HEADROOM_RATIO;
}

static inline double orbit_gain_amp(int prn)
{
    double dB = 0.0;
    if (is_meo_prn(prn))
        dB = GAIN_MEO_DB;
    else if (is_igso_prn(prn))
        dB = GAIN_IGSO_DB;
    else /* GEO */
        dB = GAIN_GEO_DB;
    return pow(10.0, dB/20.0);
}

/* 指數平滑振幅，避免 AM 旁帶 */
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
    for(int p=1;p<=prn_max;++p)
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
    if(ms == 0){
        /* 在子幀邊界，強制對齊三個時序：PRN、NH20、NAV */
        c->code_phase = 0.0;     /* 1 ms PRN 2046 chips 從頭開始 */
        c->ms_count   = 0;       /* NH20 index 從 0 開始 */
        c->bit_ptr    = 0;       /* 50 bps NAV bit 從 0 開始 */
    }

    double sf_start_d2 = floor(sow/0.6)*0.6;      /* 每 0.6 s 一個子帧 */
    double mf_start_d2 = floor(sow/3.0)*3.0;      /* 周内秒对齐到 3 s 主帧 */
    c->sf_id_d2 = ((int)(sf_start_d2/0.6))%5 + 1;
    int ms2 = (int)llround((sow - sf_start_d2)*1000.0);
    if(ms2 < 0)      ms2 = 0;
    else if(ms2 >= 600) ms2 = 599;
    c->bit_ptr_d2 = ms2/2;
    c->ms_count_d2 = ms2%2;
    /* D2 的 SOW 以主帧(3 s) 的子帧1同步頭對齊 */
    get_subframe_bits(c->prn,c->sf_id_d2,week,mf_start_d2,3.0,c->nav_bits_d2);
    /* 只在 3 s 主帧邊界重置 D2 的碼相位與導航索引 */
    if(is_d2_prn(c->prn) && fabs(sow - mf_start_d2) < 1e-9){
        c->code_phase  = 0.0;
        c->bit_ptr_d2  = 0;
        c->ms_count_d2 = 0;
    }
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
void update_channel_dynamics(channel_t *c,double rho,double rdot,double elev_deg,
                             double gain,double target_cn0,int n_visible,
                             double dt_ms)
{
    double A_new = amp_from_geom(rho,elev_deg,gain,target_cn0,n_visible);
    A_new *= orbit_gain_amp(c->prn);
    c->amp = smooth_amp(c->amp, A_new, dt_ms);
    c->amp_dot = 0.0;
    c->elev_deg = elev_deg;
    c->fd  = -FCARRIER*rdot/299792458.0;               /* Doppler (Hz) */
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
                    c->sf_id=c->sf_id%5+1;
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

/*
 * ======== D2 (500 bps) support ========
 *  - 資料速率 500 bps (2 ms per bit)
 *  - 不使用二次 Neumann-Hoffman 編碼
 */
void gen_samples_1ms_d2(channel_t *c, int week, double sow,
                               int samp_per_ms, int16_t *I, int16_t *Q)
{
    if(c->bit_ptr_d2==0 && c->ms_count_d2==0){
        double mf_start_d2 = floor(sow/3.0)*3.0;
        get_subframe_bits(c->prn,c->sf_id_d2,week,mf_start_d2,3.0,c->nav_bits_d2);
    }

    double fd = c->fd;
    double code_rate = c->code_rate;
    double phase = c->carr_phase;
    double code_phase = c->code_phase;
    double dfd = c->fd_dot / fs;
    double dcode_rate = c->code_rate_dot / fs;
    double amp = c->amp;
    double damp = c->amp_dot / fs;

    for(int n=0;n<samp_per_ms;++n){
        int chip = (int)code_phase;
        int16_t ca = ca_wave[c->prn][chip];
        int16_t nb = c->nav_bits_d2[c->bit_ptr_d2] ? -1:+1;
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
        }
    }
    c->fd = fd;
    c->code_rate = code_rate;
    c->carr_phase = phase;
    c->code_phase = code_phase;
    c->amp = amp;

    if(++c->ms_count_d2==2){
        c->ms_count_d2=0;
        if(++c->bit_ptr_d2==300){
            c->bit_ptr_d2=0;
            c->sf_id_d2 = c->sf_id_d2%5 + 1;
        }
    }
}

