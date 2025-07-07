#include "unity.h"
#include "helpers.h"
#include "channel.h"
#define CHIPRATE 2.046e6
#include "navbits.h"
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

void test_amplitude_and_interleave(void)
{
    channel_t c; setup_ch(&c);
    int16_t I[CODE_LEN],Q[CODE_LEN];
    gen_samples_1ms(&c,0,0,CODE_LEN,I,Q);
    TEST_ASSERT_EQUAL_INT16(-1000,I[0]);
    TEST_ASSERT_EQUAL_INT16(0,Q[0]);
}

void test_dc_offset(void)
{
    channel_t c; setup_ch(&c);
    int16_t I[CODE_LEN],Q[CODE_LEN];
    long sumI=0,sumQ=0;
    for(int m=0;m<20;m++){
        gen_samples_1ms(&c,0,0,CODE_LEN,I,Q);
        for(int i=0;i<CODE_LEN;i++){ sumI+=I[i]; sumQ+=Q[i]; }
    }
    TEST_ASSERT_INT_WITHIN(5,0,sumI/CODE_LEN/20);
    TEST_ASSERT_INT_WITHIN(5,0,sumQ/CODE_LEN/20);
}

