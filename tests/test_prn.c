#include "unity.h"
#include "helpers.h"
#include "globals.h"
#include <string.h>
#include <stdlib.h>


/*
 * Expected first 32 chips of PRN1.  The sequence is taken from
 * BeiDou B1I ICD Appendix A.1 and serves as a quick sanity check
 * that the Gold code generator is wired correctly.
 */
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

/*
 * Autocorrelation should equal CODE_LEN at zero lag and remain
 * below ~1/12 of that for all other shifts.  The theoretical
 * sidelobe level for 2046 length Gold codes is |R(\tau)| \le 162.
 */
void test_prn_autocorr(void)
{
    int peak=0;
    for(int i=0;i<CODE_LEN;i++)
        peak += chip_val(prn_code[1][i])*chip_val(prn_code[1][i]);
    TEST_ASSERT_EQUAL_INT(CODE_LEN, peak);

    int max_side = 0;
    for(int lag=1; lag<CODE_LEN; lag++){
        int sum=0;
        for(int i=0;i<CODE_LEN;i++){
            int j=(i+lag)%CODE_LEN;
            sum += chip_val(prn_code[1][i])*chip_val(prn_code[1][j]);
        }
        if(abs(sum) > max_side) max_side = abs(sum);
    }
    TEST_ASSERT_LESS_OR_EQUAL_INT(162, max_side);
}

/* Cross-correlation between PRN1 and PRN2 should stay below
 * about 150 chips for any lag.  This fails if the tap table or
 * feedback logic in globals.c is modified incorrectly.
 */
void test_prn_crosscorr(void)
{
    int max = 0;
    for(int lag=0; lag<CODE_LEN; lag++){
        int sum=0;
        for(int i=0;i<CODE_LEN;i++){
            int j=(i+lag)%CODE_LEN;
            sum += chip_val(prn_code[1][i])*chip_val(prn_code[2][j]);
        }
        if(abs(sum) > max) max = abs(sum);
    }
    TEST_ASSERT_LESS_OR_EQUAL_INT(162, max);
}
