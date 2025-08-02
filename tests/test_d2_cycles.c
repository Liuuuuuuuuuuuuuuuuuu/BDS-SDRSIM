#include <stdio.h>
#include <stdint.h>
#include <math.h>
#include <assert.h>
#include "../bdssim.h"
#include "../channel.h"
#include "../globals.h"

int main(void)
{
    sim_config_t cfg = {0};
    snprintf(cfg.rinex_file, sizeof(cfg.rinex_file),
             "BRDM00DLR_S_20251760000_01D_MN.rnx");
    if(!init_simulator(&cfg, 0.0)){
        fprintf(stderr, "init failed\n");
        return 1;
    }

    int prn = 3; /* GEO uses D2 without NH */
    channel_t ch;
    channel_reset(&ch, prn, nav_week, 0.0);
    ch.carr_phase = 0.0;
    ch.code_phase = 0.1; /* avoid boundary ambiguity */
    update_channel_dynamics(&ch, 2.0e7, 0.0, 90.0, 1.0, cfg.target_cn0, 1);
    channel_set_fs(FS_OUTPUT_HZ);

    int samp_per_ms = (int)(FS_OUTPUT_HZ/1000.0 + 0.5);
    int16_t I[samp_per_ms], Q[samp_per_ms];

    /* generate 2 ms and verify D2 bit spans two PRN cycles */
    for(int i=0;i<2;i++){
        gen_samples_1ms_d2(&ch, nav_week, i*0.001, samp_per_ms, I, Q);
        if(i==0){
            assert(ch.ms_count_d2 == 1);
            assert(ch.bit_ptr_d2 == 0);
        } else {
            assert(ch.ms_count_d2 == 0);
            assert(ch.bit_ptr_d2 == 1);
        }
        assert(fabs(ch.code_phase - 0.1) < 1e-6);
    }

    cleanup_simulator();
    return 0;
}
