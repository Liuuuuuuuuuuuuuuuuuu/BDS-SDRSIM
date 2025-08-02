/* channel.h ------------------------------------------------------------- */
#ifndef CHANNEL_H
#define CHANNEL_H
#include <stdint.h>

#define MAX_CH     8
#define CODE_LEN   2046

typedef struct {
    int     prn;
    double  amp;
    double  fd;           /* Doppler                   */
    double  code_rate;    /* chipping rate             */
    double  carr_phase;   /* rad                       */
    double  code_phase;   /* chips                     */
    double  elev_deg;     /* satellite elevation (deg) */
    /* ---- 10 Hz 幾何更新 + 樣本內線性內插 ---- */
    double  f_inst;       /* instantaneous carrier freq (Hz), incl. Doppler */
    double  fdot;         /* carrier frequency slope (Hz/s) over next 100 ms */
    double  R_inst;       /* instantaneous code rate (chips/s) */
    double  Rdot;         /* code-rate slope (chips/s^2) over next 100 ms */
    uint16_t code_ptr;
    uint16_t bit_ptr;
    uint8_t  sf_id;
    uint8_t  ms_count;      /* 0~19: ms index within data bit */
    uint8_t  nav_bits[300]; /* cached subframe bits */
    /* ---- D2 state (500 bps) ---- */
    uint16_t bit_ptr_d2;
    uint8_t  sf_id_d2;
    uint8_t  ms_count_d2;   /* 0~1 */
    uint8_t  nav_bits_d2[300];
} channel_t;

void channel_reset(channel_t *, int prn, int week, double sow);
void channel_set_time(channel_t *, int week, double sow);
void update_channel_dynamics(channel_t *, double rho, double rdot,
                             double elev_deg, double gain,
                             double target_cn0, int n_visible);
/* 10 Hz 幾何：在 (week,sow) 與 (sow+0.1) 估計 fdot/Rdot，並更新幅度 */
void update_channel_dynamics_10hz(channel_t *, int week, double sow,
                                  const double usr_xyz[3], const double usr_vel_eci[3],
                                  double gain, double target_cn0, int n_visible);
void channel_set_fs(double sample_rate);
void gen_samples_1ms(channel_t *, int week, double sow,
                     int samp_per_ms, int16_t *I, int16_t *Q);
void gen_samples_1ms_d2(channel_t *, int week, double sow,
                        int samp_per_ms, int16_t *I, int16_t *Q);
int  is_d2_prn(int prn);
int  is_igso_prn(int prn);
int  is_meo_prn(int prn);

/* Global target CN0 (dB-Hz) */
extern double g_target_cn0;

/* === 修正：用指標而非 array declarator === */
void get_subframe_bits(int prn, int sf_id, int week, double sow,
                       double frame_len, uint8_t *out);   /* out[300] */

#endif

