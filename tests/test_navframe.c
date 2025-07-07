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
    /* The 11-bit preamble 0x570 is defined in the BeiDou B1I spec. */
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
    int count=0;
    for(int i=0;i<SUBFRAME_BITS;i++){
        TEST_ASSERT_TRUE_MESSAGE(bits[i]==0 || bits[i]==1,
                                "bit value not 0/1");
        count++; /* every element should be filled */
    }
    TEST_ASSERT_EQUAL_INT(300, count);
}

int test_navframe_main(void)
{
    RUN_TEST(test_subframe1_sync);
    RUN_TEST(test_subframe1_repeat);
    RUN_TEST(test_subframe_bits_length);
    return 0;
}
