/* globals.c: global variables and PRN generation */

#include <stdio.h>
#include <string.h>
#include <stdint.h>
#include <stdbool.h>

#include "globals.h"
#include "bdssim.h"
#include "rinex.h"

/* ───────────── 全域資料區 ───────────── */
uint8_t      prn_code[MAX_SAT][CODE_LEN];   /* Gold 2046-chip */
ephemeris_t  eph[MAX_SAT];                  /* ★ 真正配置記憶體 ★ */

double nav_time_min = 0.0;
double nav_time_max = 0.0;
int    nav_week     = 0;
double iono_alpha[4] = {0};
double iono_beta[4]  = {0};
int    utc_bdt_diff  = 4;   /* default UTC->BDT offset */

int simulator_inited = 0;
int prn_max = 63;

/* Global transmit time (BDT seconds) */
double g_t_tx = 0.0;

/* ───────────── B1I PRN 產生 ───────────── */
#define ITER 2047                     /* 先跑滿 2047，再丟掉最後 1 chip */

/* 求 11-bit 內容在 mask 位置上的 XOR parity */
static inline uint8_t parity(uint16_t x)
{
    x ^= x >> 8;  x ^= x >> 4;  x ^= x >> 2;  x ^= x >> 1;
    return x & 1;
}

/* Fibonacci LFSR 單步：整串左移，回饋 bit 塞到 bit0 */
static inline uint16_t lfsr_step(uint16_t s, uint16_t mask)
{
    uint8_t fb = parity(s & mask);
    return ((s << 1) & 0x7FE) | fb;         /* 只保留 11 位 */
}

/* 生成指定 PRN (1–63) 的 2046-chip 代碼 */
static void cb1i_generate(int prn, uint8_t *dst)
{
    if (prn < 1 || prn > 63) return;

    uint16_t g1 = 0x2AA;   /* 01010101010b */
    uint16_t g2 = 0x2AA;
    const uint8_t *tap = g2_taps[prn];

    for (int i = 0; i < ITER; ++i) {
        uint8_t g1_out = (g1 >> 10) & 1;              /* stage 11 */
        uint8_t g2_out = ((g2 >> (tap[0]-1)) & 1) ^
                         ((g2 >> (tap[1]-1)) & 1) ^
                         (tap[2] ? ((g2 >> (tap[2]-1)) & 1) : 0);

        if (i < CODE_LEN) dst[i] = g1_out ^ g2_out;

        g1 = lfsr_step(g1, 0x7C1   /* 1,7,8,9,10,11 */);
        g2 = lfsr_step(g2, 0x59F   /* 1,2,3,4,5,8,9,11 */);
    }
}

static void init_prn_table(void)
{
    for (int p = 1; p <= prn_max; ++p)
        cb1i_generate(p, prn_code[p]);
    printf("[bdssim] PRN 表已產生 (1–%d)\n", prn_max);
}

/* ───────────── 對外 API ───────────── */
bool init_simulator(sim_config_t *cfg, double start_bdt)
{
    if(simulator_inited) return true;
    prn_max = cfg->prn37_only ? 37 : 63;
    init_prn_table();
    memset(eph, 0, sizeof(eph));
    if(read_rinex_nav(cfg->rinex_file, start_bdt)!=0) return false;
    nav_week = (int)(nav_time_min/604800.0);
    simulator_inited = 1;
    return true;
}
void cleanup_simulator(void){ simulator_inited = 0; }

