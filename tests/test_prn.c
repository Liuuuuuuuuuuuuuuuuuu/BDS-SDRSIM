#include "unity.h"
#include "helpers.h"
#include "globals.h"
#include <string.h>


/* Expected first 32 chips of PRN1 */
static const uint8_t prn1_vec[32] = {
    0,1,0,1,0,1,0,1,0,0,0,0,1,1,1,0,
    1,1,0,0,0,0,0,1,1,1,1,0,0,0,0,0
};

void test_prn_sequence(void)
{
    uint8_t buf[32];
    generate_prn(1, buf, 32);
    TEST_ASSERT_EQUAL_UINT8_ARRAY(prn1_vec, buf, 32);
}

static int chip_val(uint8_t c){ return c?1:-1; }

void test_prn_autocorr(void)
{
    int sum=0;
    for(int i=0;i<CODE_LEN;i++)
        sum += chip_val(prn_code[1][i])*chip_val(prn_code[1][i]);
    TEST_ASSERT_EQUAL_INT(CODE_LEN, sum);
}

void test_prn_crosscorr(void)
{
    int max=0;
    for(int lag=0; lag<CODE_LEN; lag++){
        int sum=0;
        for(int i=0;i<CODE_LEN;i++){
            int j=(i+lag)%CODE_LEN;
            sum += chip_val(prn_code[1][i])*chip_val(prn_code[2][j]);
        }
        if(sum>max) max=sum;
    }
    TEST_ASSERT_LESS_OR_EQUAL_INT(512, max);
}

int test_prn_main(void)
{
    TEST_ASSERT_EQUAL_INT_MESSAGE(0, load_rinex("tests/vectors/BRDM_sample.rnx"), "load rinex");
    RUN_TEST(test_prn_sequence);
    RUN_TEST(test_prn_autocorr);
    RUN_TEST(test_prn_crosscorr);
    return 0;
}
