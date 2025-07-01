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

/* ───────────── Gold 2046 產生 ───────────── */
#define G_INIT 0x2AA
#define MASK11 0x7FF
static inline uint8_t fb_g1(uint16_t r){ return ((r>>10)^(r>>9)^(r>>3)^(r>>2)^(r>>1)^r)&1; }
static inline uint8_t fb_g2(uint16_t r){ return ((r>>10)^(r>>9)^(r>>8)^(r>>7)^(r>>6)^(r>>5)^(r>>2)^(r>>1)^r)&1; }
static inline uint8_t tap(const uint16_t r,uint8_t t){ return (r>>(11-t))&1; }

/* 依 3-tap / 2-tap 表，完整 1-63 */
static const uint8_t tbl[64][3]={
 [ 0]={0},
 [ 1]={ 1, 3, 0}, [ 2]={ 1, 4, 0}, [ 3]={ 1, 5, 0}, [ 4]={ 1, 6, 0},
 [ 5]={ 1, 8, 0}, [ 6]={ 1, 9, 0}, [ 7]={ 1,10, 0}, [ 8]={ 1,11, 0},
 [ 9]={ 2, 7, 0}, [10]={ 3, 4, 0}, [11]={ 3, 5, 0}, [12]={ 3, 6, 0},
 [13]={ 3, 8, 0}, [14]={ 3, 9, 0}, [15]={ 3,10, 0}, [16]={ 3,11, 0},
 [17]={ 4, 5, 0}, [18]={ 4, 6, 0}, [19]={ 4, 8, 0}, [20]={ 4, 9, 0},
 [21]={ 4,10, 0}, [22]={ 4,11, 0}, [23]={ 5, 6, 0}, [24]={ 5, 8, 0},
 [25]={ 5, 9, 0}, [26]={ 5,10, 0}, [27]={ 5,11, 0}, [28]={ 6, 8, 0},
 [29]={ 6, 9, 0}, [30]={ 6,10, 0}, [31]={ 6,11, 0}, [32]={ 8, 9, 0},
 [33]={ 8,10, 0}, [34]={ 8,11, 0}, [35]={ 9,10, 0}, [36]={ 9,11, 0},
 [37]={10,11, 0},
 [38]={ 1, 2, 7}, [39]={ 2, 3, 7}, [40]={ 2, 4, 7}, [41]={ 2, 5, 7},
 [42]={ 2, 6, 7}, [43]={ 2, 8, 7}, [44]={ 2, 9, 7}, [45]={ 2,10, 7},
 [46]={ 2,11, 7}, [47]={ 3, 7, 8}, [48]={ 3, 7, 9}, [49]={ 3, 7,10},
 [50]={ 3, 7,11}, [51]={ 4, 7, 8}, [52]={ 4, 7, 9}, [53]={ 4, 7,10},
 [54]={ 4, 7,11}, [55]={ 5, 7, 8}, [56]={ 5, 7, 9}, [57]={ 5, 7,10},
 [58]={ 5, 7,11}, [59]={ 6, 7, 8}, [60]={ 6, 7, 9}, [61]={ 6, 7,10},
 [62]={ 6, 7,11}, [63]={ 7, 8, 9}
};

static void gen_prn(int prn)
{
    uint16_t g1=G_INIT,g2=G_INIT;
    const uint8_t *tp=tbl[prn];
    for(int i=0;i<CODE_LEN;i++){
        uint8_t out = ((g1>>10)&1) ^ tap(g2,tp[0]) ^ tap(g2,tp[1]) ^ (tp[2]?tap(g2,tp[2]):0);
        prn_code[prn][i]=out;
        g1 = ((g1<<1)|fb_g1(g1)) & MASK11;
        g2 = ((g2<<1)|fb_g2(g2)) & MASK11;
    }
}

static void init_prn_table(void){
    for(int p=1;p<=63;p++) gen_prn(p);
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

