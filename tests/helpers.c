#include "helpers.h"
#include "globals.h"
#include "bdssim.h"
#include <string.h>
#include <stdio.h>

int load_rinex(const char *fname)
{
    char full[512];
    snprintf(full, sizeof(full), "tests/vectors/%s", fname);

    FILE *fp = fopen(full, "r");
    if(!fp){
        return -1; /* missing local test vector */
    }
    fclose(fp);

    sim_config_t cfg = {0};
    strncpy(cfg.rinex_file, full, sizeof(cfg.rinex_file)-1);
    if(!init_simulator(&cfg))
        return -1;
    return 0;
}

static int prn_ready = 0;
void ensure_prn(void)
{
    if(!prn_ready){
        load_rinex("BRDM_sample.rnx");
        prn_ready = 1;
    }
}

void generate_prn(int prn, uint8_t *buf, size_t len)
{
    for(size_t i=0;i<len && i<CODE_LEN;i++)
        buf[i] = prn_code[prn][i];
}

uint32_t hamming_parity(uint32_t info)
{
    int parity_pos[] = {1,2,4,8,16};
    uint32_t word = 0;
    int di = 0;
    for(int bit=1; bit<=32; bit++){
        if(bit==1||bit==2||bit==4||bit==8||bit==16||bit==32) continue;
        if(info & (1u<<di)) word |= 1u<<(bit-1);
        di++;
    }
    int overall = 0;
    for(int k=0;k<5;k++){
        uint32_t mask = 0;
        for(int bit=1; bit<=32; bit++)
            if(bit & (1<<k)) mask |= 1u<<(bit-1);
        int p = __builtin_parity(word & mask);
        if(p) word |= 1u<<(parity_pos[k]-1);
        overall ^= p;
    }
    overall ^= __builtin_parity(word & ~(1u<<(32-1)));
    if(overall & 1) word |= 1u<<(32-1);
    return word;
}
