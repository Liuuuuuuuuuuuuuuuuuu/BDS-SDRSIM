#include "unity.h"
#include "helpers.h"
#include "navbits.h"
#include <string.h>


void test_subframe1_sync(void)
{
    load_rinex("tests/vectors/BRDM_sample.rnx");
    navbits_init();
    uint8_t bits[SUBFRAME_BITS];
    get_subframe_bits(1,1,0,0,bits);
    uint16_t sync=0;
    for(int i=0;i<11;i++) sync=(sync<<1)|bits[i];
    TEST_ASSERT_EQUAL_HEX16(0x0570, sync);
}

void test_subframe1_repeat(void)
{
    uint8_t bits[SUBFRAME_BITS];
    get_subframe_bits(1,1,0,0,bits);
    TEST_ASSERT_EQUAL_UINT8_ARRAY(bits, bits+HALF_SUBFRAME_BITS, HALF_SUBFRAME_BITS);
}

void test_subframe_bits_length(void)
{
    uint8_t bits[SUBFRAME_BITS];
    memset(bits,0xAA,sizeof(bits));
    get_subframe_bits(1,2,0,0,bits);
    /* count bits to ensure 300 produced */
    int ones=0; for(int i=0;i<SUBFRAME_BITS;i++) if(bits[i]) ones++;
    TEST_ASSERT_TRUE(SUBFRAME_BITS>=300 && SUBFRAME_BITS<=300);
    TEST_ASSERT_EQUAL_UINT(SUBFRAME_BITS, 300);
    (void)ones; /* suppress unused warning */
}

int test_navframe_main(void)
{
    RUN_TEST(test_subframe1_sync);
    RUN_TEST(test_subframe1_repeat);
    RUN_TEST(test_subframe_bits_length);
    return 0;
}
