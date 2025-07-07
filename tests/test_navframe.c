#include "unity.h"
#include "helpers.h"
#include "navbits.h"
#include "bch.h"
#include <string.h>


void test_subframe1_sync(void)
{
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
    navbits_init();
    uint8_t bits[SUBFRAME_BITS];
    get_subframe_bits(1,1,0,0,bits);
    TEST_ASSERT_EQUAL_UINT8_ARRAY(bits, bits+HALF_SUBFRAME_BITS, HALF_SUBFRAME_BITS);
}

void test_subframe_bits_length(void)
{
    navbits_init();
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

void test_subframe1_dummy_zero(void)
{
    navbits_init();
    uint8_t bits[SUBFRAME_BITS];
    get_subframe_bits(1,1,0,0,bits);
    for(int i=SUBFRAME_BITS-6;i<SUBFRAME_BITS;i++)
        TEST_ASSERT_EQUAL_UINT8(0,bits[i]);
}

void test_subframe1_parity(void)
{
    navbits_init();
    uint8_t bits[SUBFRAME_BITS];
    get_subframe_bits(1,1,0,0,bits);
    for(int w=0; w<SUBFRAME_BITS/30; w++){
        uint32_t info=0; uint16_t enc=0;
        for(int i=0;i<15;i++){
            info = (info<<1) | bits[w*30 + 2*i];
            enc  |= bits[w*30 + 2*i + 1] << i;
        }
        uint16_t expect = bch_encode(info) & ~0x8;
        TEST_ASSERT_EQUAL_HEX16(expect, enc);
    }
}
