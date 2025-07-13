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

int simulator_inited = 0;

/* ───────────── B1I PRN 產生 ───────────── */
#define LFSR_MASK 0x7FF          /* 11 bits */

static const int g1_taps[] = {1, 5, 6, 7, 10, 11};
static const int g1_tap_cnt = sizeof(g1_taps) / sizeof(g1_taps[0]);

static const int g2_taps[] = {1, 2, 3, 4, 5, 8, 9, 11};
static const int g2_tap_cnt = sizeof(g2_taps) / sizeof(g2_taps[0]);

static const int8_t g2_prn_taps[64][3] = {
    { -1, -1, -1 },
    {  1,  3, -1 }, {  1,  4, -1 }, {  1,  5, -1 }, {  1,  6, -1 },
    {  1,  8, -1 }, {  1,  9, -1 }, {  1, 10, -1 }, {  1, 11, -1 },
    {  2,  7, -1 }, {  3,  4, -1 }, {  3,  5, -1 }, {  3,  6, -1 },
    {  3,  8, -1 }, {  3,  9, -1 }, {  3, 10, -1 }, {  3, 11, -1 },
    {  4,  5, -1 }, {  4,  6, -1 }, {  4,  8, -1 }, {  4,  9, -1 },
    {  4, 10, -1 }, {  4, 11, -1 }, {  5,  6, -1 }, {  5,  8, -1 },
    {  5,  9, -1 }, {  5, 10, -1 }, {  5, 11, -1 }, {  6,  8, -1 },
    {  6,  9, -1 }, {  6, 10, -1 }, {  6, 11, -1 }, {  8,  9, -1 },
    {  8, 10, -1 }, {  8, 11, -1 }, {  9, 10, -1 }, {  9, 11, -1 },
    { 10, 11, -1 }, {  1,  2,  7 }, {  1,  3,  4 }, {  1,  3,  6 },
    {  1,  3,  8 }, {  1,  3, 10 }, {  1,  3, 11 }, {  1,  4,  5 },
    {  1,  4,  9 }, {  1,  5,  6 }, {  1,  5,  8 }, {  1,  5, 10 },
    {  1,  5, 11 }, {  1,  6,  9 }, {  1,  8,  9 }, {  1,  9, 10 },
    {  1,  9, 11 }, {  2,  3,  7 }, {  2,  5,  7 }, {  2,  7,  9 },
    {  3,  4,  5 }, {  3,  4,  9 }, {  3,  5,  6 }, {  3,  5,  8 },
    {  3,  5, 10 }, {  3,  5, 11 }, {  3,  6,  9 }
};

static inline uint8_t get_stage(uint16_t state, int k)
{
    return (state >> (11 - k)) & 1u;
}

static inline uint8_t parity(uint16_t state, const int *taps, int cnt)
{
    uint8_t p = 0;
    for (int i = 0; i < cnt; ++i) p ^= get_stage(state, taps[i]);
    return p;
}

static inline uint8_t lfsr_step(uint16_t *state, const int *fb_taps, int cnt)
{
    uint8_t out = (*state >> 10) & 1u;
    uint8_t fb  = parity(*state, fb_taps, cnt);
    *state = ((*state << 1) & LFSR_MASK) | fb;
    return out;
}

static void cb1i_generate(int prn, uint8_t *dst)
{
    if (prn < 1 || prn > 63) return;

    uint16_t g1 = 0x2AA;
    uint16_t g2 = 0x2AA;
    const int8_t *sel = g2_prn_taps[prn];

    for (int i = 0; i < CODE_LEN + 1; ++i) {
        uint8_t g1_chip = lfsr_step(&g1, g1_taps, g1_tap_cnt);

        uint8_t g2_sel = 0;
        for (int k = 0; k < 3 && sel[k] != -1; ++k)
            g2_sel ^= get_stage(g2, sel[k]);

        if (i < CODE_LEN) dst[i] = g1_chip ^ g2_sel;

        lfsr_step(&g2, g2_taps, g2_tap_cnt);
    }
}

static void init_prn_table(void)
{
    for (int p = 1; p <= 63; ++p)
        cb1i_generate(p, prn_code[p]);
    puts("[bdssim] PRN 表已產生 (1–63)");
}

/* ───────────── 對外 API ───────────── */
bool init_simulator(sim_config_t *cfg)
{
    if(simulator_inited) return true;
    init_prn_table();
    if(read_rinex_nav(cfg->rinex_file)!=0) return false;
    nav_week = (int)(nav_time_min/604800.0);
    simulator_inited = 1;
    return true;
}
void cleanup_simulator(void){ simulator_inited = 0; }

