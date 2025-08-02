/* channel.h ------------------------------------------------------------- */
#ifndef CHANNEL_H
#define CHANNEL_H
#include <stdint.h>

#define MAX_CH     8
#define CODE_LEN   2046

typedef struct {
    int     prn;
    double  amp;
    double  amp_dot;      /* amplitude rate (per second) */
    double  fd;           /* Doppler                   */
    double  fd_dot;       /* Doppler rate (Hz/s)       */
    double  code_rate;    /* chipping rate             */
    double  code_rate_dot;/* chipping rate slope       */
    double  carr_phase;   /* rad                       */
    double  code_phase;   /* chips                     */
    double  elev_deg;     /* satellite elevation (deg) */
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
                             double target_cn0, int n_visible,
                             double dt_ms);
void channel_set_fs(double sample_rate);
void gen_samples_1ms(channel_t *, int week, double sow,
                     int samp_per_ms, int16_t *I, int16_t *Q);
void gen_samples_1ms_d2(channel_t *, int week, double sow,
                        int samp_per_ms, int16_t *I, int16_t *Q);
int  is_d2_prn(int prn);
int  is_igso_prn(int prn);
int  is_meo_prn(int prn);
double predict_next_amp(const channel_t *c,double elev_deg_next,double gain,
                        double target_cn0,int n_visible,double dt_ms);

/* Global target CN0 (dB-Hz) */
extern double g_target_cn0;

/* === 修正：用指標而非 array declarator === */
void get_subframe_bits(int prn, int sf_id, int week, double sow,
                       double frame_len, uint8_t *out);   /* out[300] */

#endif

