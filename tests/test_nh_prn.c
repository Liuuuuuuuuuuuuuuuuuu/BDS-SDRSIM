#include <stdio.h>
#include <stdint.h>
#include <math.h>
#include <assert.h>
#include "../bdssim.h"
#include "../channel.h"
#include "../globals.h"

static const int nh20_bits[20]={
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

    int prn = 8; /* MEO/IGSO uses D1 with NH pilot */
    channel_t ch;
    channel_reset(&ch, prn, nav_week, 0.0);
    ch.carr_phase = 0.0;
    ch.code_phase = 0.0; /* start at chip 0 */
    update_channel_dynamics(&ch, 2.0e7, 0.0, 90.0, 1.0, cfg.target_cn0, 1, 1.0);
    channel_set_fs(FS_OUTPUT_HZ);

    int samp_per_ms = (int)(FS_OUTPUT_HZ/1000.0 + 0.5);
    int16_t I[samp_per_ms], Q[samp_per_ms];

    /* first ms: verify I carries D1 and Q carries NH */
    gen_samples_1ms(&ch, nav_week, 0.0, samp_per_ms, I, Q);
    int ca = prn_code[prn][0] ? +1 : -1;
    int d1 = ch.nav_bits[0] ? -1 : +1;
    int nh = nh20_bits[0] ? -1 : +1;
    assert(I[0] * (ca*d1) > 0);
    assert(Q[0] * (ca*nh) > 0);

    /* remaining 19 ms: check NH cycle and bit pointer */
    for(int i=1;i<20;i++){
        gen_samples_1ms(&ch, nav_week, i*0.001, samp_per_ms, I, Q);
        if(i<19){
            assert(ch.ms_count == (i+1));
            assert(ch.bit_ptr == 0);
        } else {
            assert(ch.ms_count == 0);
            assert(ch.bit_ptr == 1);
        }
        assert(fabs(ch.code_phase - 0.0) < 1e-6);
    }

    cleanup_simulator();
    return 0;
}
