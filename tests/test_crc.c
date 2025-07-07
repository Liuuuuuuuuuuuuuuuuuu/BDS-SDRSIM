#include "unity.h"
#include "helpers.h"


void test_zero_word(void)
{
    uint32_t w = hamming_parity(0x0);
    TEST_ASSERT_EQUAL_HEX32(0x00000000, w);
}

void test_one_word(void)
{
    uint32_t w = hamming_parity(0x1);
    TEST_ASSERT_EQUAL_HEX32(0x80000007, w);
}

void test_full_word(void)
{
    uint32_t w = hamming_parity(0x1ffffff);
    TEST_ASSERT_EQUAL_HEX32(0xbfff7f74, w);
}

