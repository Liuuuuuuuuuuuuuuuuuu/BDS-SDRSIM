#ifndef CHANNEL_H
#define CHANNEL_H

#include <stdint.h>

/* GEO PRNs and NH20/antenna patterns */
static const int geo_prn[] = {1,2,3,4,5,59,60,61,62,63};
static const double ant_pat_db[37] = {
     0.00,  0.00,  0.22,  0.44,  0.67,  1.11,  1.56,  2.00,  2.44,  2.89,
     3.56,  4.22,  4.89,  5.56,  6.22,  6.89,  7.56,  8.22,  8.89,  9.78,
    10.67, 11.56, 12.44, 13.33, 14.44, 15.56, 16.67, 17.78, 18.89, 20.00,
    21.33, 22.67, 24.00, 25.56, 27.33, 29.33, 31.56
};
static const uint8_t nh20_bits[20] = {
    0,0,0,0,0,1,0,0,1,1,0,1,0,1,0,0,1,1,1,0
};

typedef struct {
    int     prn;
    double  amp;
    double  amp_dot;
    double  fd;
    double  fd_dot;
    double  code_rate;
    double  code_rate_dot;
    double  carr_phase;
    double  code_phase;
    double  elev_deg;
    uint16_t bit_ptr;
    uint8_t  sf_id;
    uint8_t  ms_count;
    uint8_t  nav_bits[300];
    double  nav_sf_start;
} channel_t;

void channel_reset(channel_t *, int prn, double rho);
void channel_set_time(channel_t *, double rho);
void update_channel_dynamics(channel_t *, double rho, double rdot,
                             double elev_deg, double gain,
                             double target_cn0, int n_visible,
                             double dt_ms);
void channel_set_fs(double sample_rate);
void gen_samples_1ms(channel_t *, int samp_per_ms, int16_t *I, int16_t *Q);
int  is_geo_prn(int prn);
int  is_igso_prn(int prn);
int  is_meo_prn(int prn);
double predict_next_amp(const channel_t *c, double rho_next,
                        double elev_deg_next, double gain,
                        double target_cn0, int n_visible, double dt_ms);

#endif /* CHANNEL_H */

