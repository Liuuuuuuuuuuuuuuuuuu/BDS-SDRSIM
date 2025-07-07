#include "unity.h"

void setUp(void){}
void tearDown(void){}
extern int test_prn_main(void);
extern int test_crc_main(void);
extern int test_navframe_main(void);
extern int test_chipseq_main(void);
extern int test_scaling_main(void);
extern int test_fft_main(void);

int main(void)
{
    UNITY_BEGIN();
    test_prn_main();
    test_crc_main();
    test_navframe_main();
    test_chipseq_main();
    test_scaling_main();
    test_fft_main();
    return UNITY_END();
}
