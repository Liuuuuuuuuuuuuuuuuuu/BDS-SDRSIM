#include "unity.h"
#include "helpers.h"

void setUp(void)
{
    ensure_prn();
}

void tearDown(void){}

/* Declarations from test files */
void test_prn_sequence(void);
void test_prn_autocorr(void);
void test_prn_crosscorr(void);

void test_zero_word(void);
void test_one_word(void);
void test_full_word(void);

void test_subframe1_sync(void);
void test_subframe1_repeat(void);
void test_subframe_bits_length(void);

void test_chipseq_length_advances(void);
void test_navbit_inversion(void);
void test_bitptr_after_20ms(void);

void test_amplitude_and_interleave(void);
void test_dc_offset(void);

void test_fft_dc(void);
void test_fft_power(void);

int main(void)
{
    UNITY_BEGIN();

    RUN_TEST(test_prn_sequence);
    RUN_TEST(test_prn_autocorr);
    RUN_TEST(test_prn_crosscorr);

    RUN_TEST(test_zero_word);
    RUN_TEST(test_one_word);
    RUN_TEST(test_full_word);

    RUN_TEST(test_subframe1_sync);
    RUN_TEST(test_subframe1_repeat);
    RUN_TEST(test_subframe_bits_length);

    RUN_TEST(test_chipseq_length_advances);
    RUN_TEST(test_navbit_inversion);
    RUN_TEST(test_bitptr_after_20ms);

    RUN_TEST(test_amplitude_and_interleave);
    RUN_TEST(test_dc_offset);

    RUN_TEST(test_fft_dc);
    RUN_TEST(test_fft_power);

    return UNITY_END();
}
