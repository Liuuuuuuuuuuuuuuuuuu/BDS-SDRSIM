/* globals.c : 全域變數與 PRN 產生、init/cleanup */

#include <stdio.h>
#include <string.h>
#include <stdint.h>
#include <stdbool.h>

#include "bdssim.h"
#include "channel.h"   /* 取 CODE_LEN 宏 */
#include "rinex.h"     /* read_rinex_nav() */

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

/* ───────────── B1I PRN 產生 ───────────── */
#define ITER 2047                     /* 先跑滿 2047，再丟掉最後 1 chip */

/* G2 相位抽頭表 (Table 4‑1/4‑2) */
/* 兩抽頭→tap3=0；三抽頭→全部填入，stage 編號由 1 起算 */
static const uint8_t g2_taps[64][3] = {
 /*0*/ {0,0,0},
 /*1*/ {1,3,0},{1,4,0},{1,5,0},{1,6,0},{1,8,0},
 /*6*/ {1,9,0},{1,10,0},{1,11,0},{2,7,0},{3,4,0},
 /*11*/{3,5,0},{3,6,0},{3,8,0},{3,9,0},{3,10,0},{3,11,0},
 /*17*/{4,5,0},{4,6,0},{4,8,0},{4,9,0},{4,10,0},{4,11,0},
 /*23*/{5,6,0},{5,8,0},{5,9,0},{5,10,0},{5,11,0},
 /*28*/{6,8,0},{6,9,0},{6,10,0},{6,11,0},
 /*32*/{8,9,0},{8,10,0},{8,11,0},
 /*35*/{9,10,0},{9,11,0},{10,11,0},
 /*38*/{1,2,7},{1,3,4},{1,3,6},{1,3,8},{1,3,10},
 /*43*/{1,3,11},{1,4,5},{1,4,9},{1,5,6},{1,5,8},
 /*48*/{1,5,10},{1,5,11},{1,6,9},{1,8,9},{1,9,10},{1,9,11},
 /*54*/{2,3,7},{2,5,7},{2,7,9},{3,4,5},{3,4,9},
 /*59*/{3,5,6},{3,5,8},{3,5,10},{3,5,11},{3,6,9}
};

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

        g1 = lfsr_step(g1, 0x671   /* 1,5,6,7,10,11 */);
        g2 = lfsr_step(g2, 0x59F   /* 1,2,3,4,5,8,9,11 */);
    }
}

static void init_prn_table(void)
{
    for (int p = 1; p <= 63; ++p)
        cb1i_generate(p, prn_code[p]);
    puts("[bdssim] PRN 表已產生 (1–63)");
}

/* ───────────── 對外 API ───────────── */
bool init_simulator(sim_config_t *cfg, double start_bdt)
{
    if(simulator_inited) return true;
    init_prn_table();
    memset(eph, 0, sizeof(eph));
    if(read_rinex_nav(cfg->rinex_file, start_bdt)!=0) return false;
    nav_week = (int)(nav_time_min/604800.0);
    simulator_inited = 1;
    return true;
}
void cleanup_simulator(void){ simulator_inited = 0; }

