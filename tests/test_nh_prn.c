#include <stdio.h>
#include <stdint.h>
#include <math.h>
#include <assert.h>
#include "../bdssim.h"
#include "../channel.h"
#include "../globals.h"

/* Expected NH20 sequence: 0 -> +1, 1 -> -1 */
static const uint8_t nh20_bits[20]={
    0,0,0,0,0,1,0,0,1,1,0,1,0,1,0,0,1,1,1,0
};

int main(void)
{
    sim_config_t cfg = {0};
    snprintf(cfg.rinex_file, sizeof(cfg.rinex_file),
             "BRDM00DLR_S_20251760000_01D_MN.rnx");
    if(!init_simulator(&cfg, 0.0)){
        fprintf(stderr, "init failed\n");
        return 1;
    }

    int prn = 8; /* MEO/IGSO uses D1 with NH code */
    channel_t ch;
    channel_reset(&ch, prn, nav_week, 0.0, 0.0);
    ch.carr_phase = 0.0;
    ch.code_phase = 0.1; /* avoid boundary ambiguity */
    update_channel_dynamics(&ch, 2.0e7, 0.0, 90.0, 1.0, cfg.target_cn0, 1, 1.0);
    channel_set_fs(FS_OUTPUT_HZ);

    int samp_per_ms = (int)(FS_OUTPUT_HZ/1000.0 + 0.5);
    int16_t I[samp_per_ms], Q[samp_per_ms];

    /* generate 20 ms and verify NH/PRN relationship */
    int base_sign = 0;
    for(int i=0;i<20;i++){
        gen_samples_1ms(&ch, nav_week, i*0.001, samp_per_ms, I, Q);
        int sign = (I[0] >= 0) ? 1 : -1;
        if(i==0){
            base_sign = sign;               /* nav bit sign */
        } else {
            int expect = nh20_bits[i] ? -base_sign : base_sign;
            assert(sign == expect);         /* NH flips sign when bit =1 */
        }
        if(i<19){
            assert(ch.ms_count == (i+1));
            assert(ch.bit_ptr == 0);
        } else {
            assert(ch.ms_count == 0);
            assert(ch.bit_ptr == 1);
        }
        assert(fabs(ch.code_phase - 0.1) < 1e-6);
    }

    cleanup_simulator();
    return 0;
}
