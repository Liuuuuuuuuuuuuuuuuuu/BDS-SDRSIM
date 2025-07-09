#include <stdio.h>
#include <stdint.h>
#include <assert.h>
#include "../bdssim.h"
#include "../channel.h"
#include "../globals.h"

int main(void){
    sim_config_t cfg = {0};
    snprintf(cfg.rinex_file, sizeof(cfg.rinex_file),
             "BRDM00DLR_S_20251760000_01D_MN.rnx");
    if(!init_simulator(&cfg)){
        fprintf(stderr, "init failed\n");
        return 1;
    }

    assert(is_d2_prn(3));
    assert(!is_d2_prn(8));

    channel_t c3, c8;
    channel_reset(&c3, 3, nav_week, 0.0);
    channel_reset(&c8, 8, nav_week, 0.0);
    channel_set_fs(5.12e6);
    int samp_per_ms = (int)(5.12e6/1000.0 + 0.5);
    int16_t I[samp_per_ms], Q[samp_per_ms];

    for(int i=0;i<10;i++){
        gen_samples_1ms_d2(&c3, nav_week, i*0.001, samp_per_ms, I, Q);
        gen_samples_1ms(&c8, nav_week, i*0.001, samp_per_ms, I, Q);
    }

    if(c3.bit_ptr_d2 < 1 || c3.bit_ptr_d2 != 5){
        fprintf(stderr, "PRN3 not D2 as expected (bit_ptr_d2=%u)\n", c3.bit_ptr_d2);
        return 1;
    }
    if(c8.bit_ptr != 0){
        fprintf(stderr, "PRN8 unexpected bit progress (bit_ptr=%u)\n", c8.bit_ptr);
        return 1;
    }

    cleanup_simulator();
    return 0;
}
