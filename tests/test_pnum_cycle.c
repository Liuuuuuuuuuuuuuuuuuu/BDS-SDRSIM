#include <stdio.h>
#include <stdint.h>
#include <math.h>
#include <assert.h>
#include "../bdssim.h"
#include "../channel.h"
#include "../globals.h"
#include "../bch.h"
#include "../navbits.h"

static uint32_t read_word(const uint8_t *bits, int idx)
{
    uint32_t w = 0;
    for(int i=0;i<30;i++)
        w = (w<<1) | bits[idx*30 + i];
    return w;
}

int main(void)
{
    sim_config_t cfg = {0};
    snprintf(cfg.rinex_file, sizeof(cfg.rinex_file),
             "BRDM00DLR_S_20251760000_01D_MN.rnx");
    if(!init_simulator(&cfg, 0.0)){
        fprintf(stderr, "init failed\n");
        return 1;
    }

    int prn = 8; /* any D1 satellite */
    uint8_t bits[SUBFRAME_BITS];

    /* check two full cycles of page numbers */
    for(int frame=0; frame<26; ++frame){
        int pnum = frame % 24 + 1;

        double t4 = frame*30.0 + 18.0; /* subframe 4 start */
        get_subframe_bits(prn, 4, nav_week, t4, 6.0, bits);
        uint32_t w2 = read_word(bits,1);
        uint32_t sow_int = (uint32_t)(floor(t4/6.0)*6.0);
        uint32_t payload = ((sow_int & 0xFFF) << 10) | (0 << 9) | ((pnum & 0x7F) << 2) | 0x2;
        uint32_t expect = bch_interleave_22bit(payload);
        assert(w2 == expect);

        double t5 = frame*30.0 + 24.0; /* subframe 5 start */
        get_subframe_bits(prn, 5, nav_week, t5, 6.0, bits);
        w2 = read_word(bits,1);
        sow_int = (uint32_t)(floor(t5/6.0)*6.0);
        payload = ((sow_int & 0xFFF) << 10) | (0 << 9) | ((pnum & 0x7F) << 2) | 0x2;
        expect = bch_interleave_22bit(payload);
        assert(w2 == expect);
    }

    cleanup_simulator();
    return 0;
}
