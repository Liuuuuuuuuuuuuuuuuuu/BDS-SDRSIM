#ifndef GLOBALS_H
#define GLOBALS_H
#define DEBUG
#include <stdint.h>
#include "bdssim.h"          /* 已含 MAX_SAT、CODE_LEN、ephemeris_t */

/* === CN0 / NH20 全域設定 === */
#define CN0_TARGET_DBHZ   42.0     /* 目標載波-雜訊比；40~45 均可 */
#define HEADROOM_RATIO    0.8      /* -2 dB 的峰值餘量 */
#define FS_OUTPUT_HZ      6144000  /* 6.144 MHz = 2.046 MHz × 3 */

/* 由 globals.c 定義 */
extern uint8_t      prn_code[MAX_SAT][CODE_LEN];
extern ephemeris_t  eph[MAX_SAT];
extern int          simulator_inited;
extern double       nav_time_min;
extern double       nav_time_max;
extern int          nav_week;       /* ephemeris week reference */
extern double       iono_alpha[4];  /* ionospheric parameters */
extern double       iono_beta[4];
extern int          utc_bdt_diff;   /* UTC->BDT offset seconds */

#endif

