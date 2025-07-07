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

static void navfile_checks(const char *fname)
{
    cleanup_simulator();
    if(load_rinex(fname)!=0)
        UNITY_TEST_FAIL(__LINE__, "missing RINEX test vector");

    navbits_init();

    int prn_found=-1;
    for(int p=1;p<MAX_SAT;p++)
        if(eph[p].prn){ prn_found=p; break; }

    TEST_ASSERT_MESSAGE(prn_found>0, "no ephemeris loaded");

    ephemeris_t *e=&eph[prn_found];
    TEST_ASSERT_EQUAL_INT(prn_found, e->prn);
    TEST_ASSERT_TRUE(e->toe > 0.0);

    uint8_t bits[SUBFRAME_BITS];
    get_subframe_bits(prn_found,1,0,0,bits);
    uint32_t info=0; uint16_t enc=0;
    for(int i=0;i<15;i++){
        info = (info<<1) | bits[i*2];
        enc  |= bits[i*2 + 1] << i;
    }
    uint16_t expect = bch_encode(info) & ~0x8;
    TEST_ASSERT_EQUAL_HEX16(expect, enc);
}

void test_navframe_file_1760(void){ navfile_checks("BRDM00DLR_S_20251760000_01D_MN.rnx"); }
void test_navframe_file_1770(void){ navfile_checks("BRDM00DLR_S_20251770000_01D_MN.rnx"); }
