/* channel.h ------------------------------------------------------------- */
#ifndef CHANNEL_H
#define CHANNEL_H
#include <stdint.h>

#define MAX_CH     12
#define CODE_LEN   2046
#define SAMP_1MS   8184

typedef struct {
    int     prn;
    double  amp;
    double  fd;           /* Doppler                   */
    double  code_rate;    /* chipping rate             */
    double  carr_phase;   /* rad                       */
    double  code_phase;   /* chips                     */
    uint16_t code_ptr;
    uint16_t bit_ptr;
    uint8_t  sf_id;
    uint8_t  ms_count;      /* 0~19: ms index within data bit */
    uint8_t  nav_bits[300]; /* cached subframe bits */
} channel_t;

void channel_reset(channel_t *, int prn);
void update_channel_dynamics(channel_t *, double rho, double rdot, int n_ch, double gain);
void gen_samples_1ms(channel_t *, int week, double sow,
                     int16_t *I, int16_t *Q);

/* === 修正：用指標而非 array declarator === */
void get_subframe_bits(int prn, int sf_id, int week, double sow,
                       uint8_t *out);   /* out[300] */

#endif

