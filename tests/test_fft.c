#include "unity.h"
#include "helpers.h"
#include "channel.h"
#include "ext/kissfft/kiss_fft.h"
#include "navbits.h"
#define CHIPRATE 2.046e6
#include <stdlib.h>
#include <math.h>
#include <string.h>


static void setup_ch(channel_t *c)
{
    channel_reset(c,1);
    c->amp = 1000.0;
    c->fd = 0.0;
    c->code_rate = CHIPRATE;
    c->carr_phase = 0.0;
    c->code_phase = 0.0;
    memset(c->nav_bits,0,sizeof(c->nav_bits));
    c->nav_bits[0]=1;
    channel_set_fs(CHIPRATE);
}

void test_fft_dc(void)
{
    channel_t c; setup_ch(&c);
    int16_t I[CODE_LEN],Q[CODE_LEN];
    gen_samples_1ms(&c,0,0,CODE_LEN,I,Q);
    kiss_fft_cfg cfg = kiss_fft_alloc(CODE_LEN,0,NULL,NULL);
    kiss_fft_cpx *in = malloc(sizeof(kiss_fft_cpx)*CODE_LEN);
    kiss_fft_cpx *out = malloc(sizeof(kiss_fft_cpx)*CODE_LEN);
    for(int i=0;i<CODE_LEN;i++){ in[i].r=I[i]; in[i].i=Q[i]; }
    kiss_fft(cfg,in,out);
    TEST_ASSERT_INT_WITHIN(40000,0,(int)out[0].r);
    free(in); free(out); free(cfg);
}

void test_fft_power(void)
{
    channel_t c; setup_ch(&c);
    int16_t I[CODE_LEN],Q[CODE_LEN];
    gen_samples_1ms(&c,0,0,CODE_LEN,I,Q);
    kiss_fft_cfg cfg = kiss_fft_alloc(CODE_LEN,0,NULL,NULL);
    kiss_fft_cpx *in = malloc(sizeof(kiss_fft_cpx)*CODE_LEN);
    kiss_fft_cpx *out = malloc(sizeof(kiss_fft_cpx)*CODE_LEN);
    double p_time=0.0;
    for(int i=0;i<CODE_LEN;i++){ in[i].r=I[i]; in[i].i=Q[i]; p_time+=I[i]*I[i]+Q[i]*Q[i]; }
    kiss_fft(cfg,in,out);
    double p_freq=0.0; for(int i=0;i<CODE_LEN;i++) p_freq+= out[i].r*out[i].r + out[i].i*out[i].i;
    free(in); free(out); free(cfg);
    TEST_ASSERT_FLOAT_WITHIN(1e-3,p_time,p_freq/CODE_LEN);
}

int test_fft_main(void)
{
    load_rinex("tests/vectors/BRDM_sample.rnx");
    navbits_init();
    RUN_TEST(test_fft_dc);
    RUN_TEST(test_fft_power);
    return 0;
}
