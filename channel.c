/* channel.c : 單通道 1 ms B1I I/Q 產生 */
#include <string.h>
#include <stdint.h>
#include <stdlib.h>
#include <math.h>
#include "channel.h"
#include "globals.h"               /* prn_code */
#include "coord.h"
#include "orbits.h"

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
static inline double amp_from_cn0(double cn0_dBHz, int n_visible)
{
    /* 以 gps-sdr-sim 作法：45 dB-Hz 為名目，轉線性後放大到 int16 尺度（~16384） */
    double base = pow(10.0, (cn0_dBHz - 45.0) / 20.0) * 16384.0;
    if (n_visible < 1) n_visible = 1;
    return (base / sqrt((double)n_visible)) * HEADROOM_RATIO;
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
static inline double smooth_amp(double A_prev, double A_new)
{
    /* Called every 1 ms, so dt = 1 ms; use time constant in milliseconds */
    double alpha = 1.0 - exp(-1.0 / AMP_SMOOTH_TC_MS);
    if (A_prev == 0.0)
        return A_new;                  /* avoid long ramp at start */
    return A_prev + alpha * (A_new - A_prev);
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
                             double gain,double target_cn0,int n_visible)
{
    (void)rho; /* rho currently unused in simplified amplitude model */
    double base = amp_from_cn0(target_cn0, n_visible);
    double s = sin(elev_deg * (M_PI/180.0));
    if (s < 0.0) s = 0.0;
    double A_new = base * orbit_gain_amp(c->prn) * pow(s, 1.2) * gain;
    c->amp = smooth_amp(c->amp, A_new);
    c->elev_deg = elev_deg;
    /* Instantaneous Doppler and code rate */
    c->fd        = -FCARRIER*rdot/299792458.0;         /* Doppler (Hz) */
    c->code_rate = CHIPRATE*(1.0 - rdot/299792458.0);  /* Code rate (Hz) */
    /* Initialize instantaneous state used by 10 Hz interpolation */
    c->f_inst = c->fd;
    c->fdot   = 0.0;
    c->R_inst = c->code_rate;
    c->Rdot   = 0.0;
}

/* 每 100 ms 呼叫一次：以 t 與 t+0.1s 的幾何，更新 f_inst/R_inst 與斜率 */
void update_channel_dynamics_10hz(channel_t *c, int week, double sow,
                                  const double usr_xyz[3], const double usr_vel_eci[3],
                                  double gain, double target_cn0, int n_visible)
{
    /* 1) 計算此刻與 +0.1s 的衛星位置/速度與 elev、rdot */
    double sat0[3], vel0[3], enu0[3];
    calc_sat_position_velocity(c->prn, week, sow, sat0, vel0);
    coord_t u0={0}; memcpy(u0.xyz, usr_xyz, sizeof(u0.xyz));
    /* 與通道無關的 LLH 可略過，僅為仰角計算 */
    ecef2enu(&u0, sat0, enu0);
    double el0 = enu_elevation_deg(enu0);
    double dx0=sat0[0]-usr_xyz[0], dy0=sat0[1]-usr_xyz[1], dz0=sat0[2]-usr_xyz[2];
    double rho0=hypot(hypot(dx0,dy0),dz0);
    double rdot0=(dx0*(vel0[0]-usr_vel_eci[0]) + dy0*(vel0[1]-usr_vel_eci[1]) + dz0*(vel0[2]-usr_vel_eci[2]))/rho0;

    double sow1 = sow + 0.1; int week1 = week;
    if (sow1 >= 604800.0) { sow1 -= 604800.0; week1++; }
    double sat1[3], vel1[3];
    calc_sat_position_velocity(c->prn, week1, sow1, sat1, vel1);
    double dx1=sat1[0]-usr_xyz[0], dy1=sat1[1]-usr_xyz[1], dz1=sat1[2]-usr_xyz[2];
    double rho1=hypot(hypot(dx1,dy1),dz1);
    double rdot1=(dx1*(vel1[0]-usr_vel_eci[0]) + dy1*(vel1[1]-usr_vel_eci[1]) + dz1*(vel1[2]-usr_vel_eci[2]))/rho1;

    /* 2) 幅度（平滑）：沿用原本的模型 */
    double base = amp_from_cn0(target_cn0, n_visible);
    double s = sin(el0 * (M_PI/180.0)); if (s < 0.0) s = 0.0;
    double A_new = base * orbit_gain_amp(c->prn) * pow(s, 1.2) * gain;
    c->amp = smooth_amp(c->amp, A_new);
    c->elev_deg = el0;

    /* 3) 設定 10 Hz 內插用的即時量與斜率 */
    double fd0 = -FCARRIER*rdot0/299792458.0;
    double fd1 = -FCARRIER*rdot1/299792458.0;
    c->f_inst = fd0;
    c->fdot   = (fd1 - fd0) / 0.1;          /* Hz/s */
    double R0 = CHIPRATE * (1.0 - rdot0/299792458.0);
    double R1 = CHIPRATE * (1.0 - rdot1/299792458.0);
    c->R_inst = R0;
    c->Rdot   = (R1 - R0) / 0.1;            /* chips/s^2 */
    /* 保留顯示用值 */
    c->fd        = fd0;
    c->code_rate = R0;
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

    const double dt = 1.0/fs;
    double code_phase = c->code_phase;
    double phase = c->carr_phase;
    double f_inst = c->f_inst;   /* Hz */
    double R_inst = c->R_inst;   /* chips/s */

    for(int n=0;n<samp_per_ms;++n){
        int chip = (int)code_phase;            /* 0..2045 */
        int16_t ca = ca_wave[c->prn][chip];
        uint8_t nh = nh20_bits[c->ms_count];
        int16_t nb = (c->nav_bits[c->bit_ptr]^nh)?-1:+1;
        float co,si; fast_sincos(phase,&co,&si);
        float s = c->amp*ca*nb;
        I[n]=(int16_t)lrintf(s*co);
        Q[n]=(int16_t)lrintf(s*si);
        /* NCO：只用連續內插（10 Hz 幾何推進） */
        f_inst += c->fdot * dt;
        phase  += (float)(PI2 * f_inst * dt);
        if(phase>=PI2)      phase-=PI2;
        else if(phase<0.0)  phase+=PI2;
        R_inst += c->Rdot * dt;
        code_phase += R_inst * dt;
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
    c->f_inst     = f_inst;
    c->R_inst     = R_inst;
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

    const double dt = 1.0/fs;
    double code_phase = c->code_phase;
    double phase = c->carr_phase;
    double f_inst = c->f_inst;
    double R_inst = c->R_inst;

    for(int n=0;n<samp_per_ms;++n){
        int chip = (int)code_phase;
        int16_t ca = ca_wave[c->prn][chip];
        int16_t nb = c->nav_bits_d2[c->bit_ptr_d2] ? -1:+1;
        float co,si; fast_sincos(phase,&co,&si);
        float s = c->amp*ca*nb;
        I[n]=(int16_t)lrintf(s*co);
        Q[n]=(int16_t)lrintf(s*si);

        /* 只用連續內插 */
        f_inst += c->fdot * dt;
        phase  += (float)(PI2 * f_inst * dt);
        if(phase>=PI2)      phase-=PI2;
        else if(phase<0.0)  phase+=PI2;
        R_inst += c->Rdot * dt;
        code_phase += R_inst * dt;
        if(code_phase>=CODE_LEN){
            code_phase-=CODE_LEN;
        }
    }
    c->carr_phase = phase;
    c->code_phase = code_phase;
    c->f_inst     = f_inst;
    c->R_inst     = R_inst;

    if(++c->ms_count_d2==2){
        c->ms_count_d2=0;
        if(++c->bit_ptr_d2==300){
            c->bit_ptr_d2=0;
            c->sf_id_d2 = c->sf_id_d2%5 + 1;
        }
    }
}

