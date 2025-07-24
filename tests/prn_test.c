#include <stdio.h>
#include "bdssim.h"

int main(void){
    sim_config_t c = {0};
    /* Use bundled RINEX file for demonstration */
    snprintf(c.rinex_file, sizeof(c.rinex_file),
             "BRDM00DLR_S_20251760000_01D_MN.rnx");
    init_simulator(&c, 0.0);

    for(int prn=1; prn<=63; ++prn){
        printf("PRN%02d:", prn);
        for(int i=0;i<16;i++)
            printf("%d", prn_code[prn][i]);
        putchar('\n');
    }
    cleanup_simulator();
    return 0;
}


