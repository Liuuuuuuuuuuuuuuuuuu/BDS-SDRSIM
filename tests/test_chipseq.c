#include "unity.h"
#include "helpers.h"
#include "channel.h"
#include "navbits.h"
#define CHIPRATE 2.046e6
#include <string.h>


static void init_channel(channel_t *c)
{
    channel_reset(c,1);
    c->amp = 1.0;
    c->fd = 0.0;
    c->code_rate = CHIPRATE;
    c->carr_phase = 0.0;
    c->code_phase = 0.0;
    memset(c->nav_bits,0,sizeof(c->nav_bits));
    c->nav_bits[0]=1; /* first nav bit = 1 */
    channel_set_fs(CHIPRATE);
}

void test_chipseq_length_advances(void)
{
    channel_t c; init_channel(&c);
    int16_t I[CODE_LEN],Q[CODE_LEN];
    gen_samples_1ms(&c,0,0,CODE_LEN,I,Q);
    TEST_ASSERT_EQUAL_UINT16(1,c.ms_count);
    gen_samples_1ms(&c,0,0,CODE_LEN,I,Q);
    TEST_ASSERT_EQUAL_UINT16(2,c.ms_count);
}

void test_navbit_inversion(void)
{
    channel_t c; init_channel(&c);
    int16_t I[CODE_LEN],Q[CODE_LEN];
    gen_samples_1ms(&c,0,0,CODE_LEN,I,Q);
    int16_t first=I[0];
    for(int i=1;i<6;i++) gen_samples_1ms(&c,0,0,CODE_LEN,I,Q);
    int16_t sixth=I[0];
    TEST_ASSERT_EQUAL_INT16(-first, sixth);
}

void test_bitptr_after_20ms(void)
{
    channel_t c; init_channel(&c);
    int16_t I[CODE_LEN],Q[CODE_LEN];
    for(int i=0;i<20;i++) gen_samples_1ms(&c,0,0,CODE_LEN,I,Q);
    TEST_ASSERT_EQUAL_UINT16(1,c.bit_ptr);
}

int test_chipseq_main(void)
{
    load_rinex("tests/vectors/BRDM_sample.rnx");
    navbits_init();
    RUN_TEST(test_chipseq_length_advances);
    RUN_TEST(test_navbit_inversion);
    RUN_TEST(test_bitptr_after_20ms);
    return 0;
}
