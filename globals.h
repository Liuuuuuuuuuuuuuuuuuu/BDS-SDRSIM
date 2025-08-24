#ifndef GLOBALS_H
#define GLOBALS_H
#define DEBUG
#include <stdint.h>
#include "bdssim.h"          /* 已含 MAX_SAT、CODE_LEN、ephemeris_t */

/* ===== gps-sdr-sim 風格的基本參數 ===== */
#define FS_OUTPUT_HZ      8192000.0   /* 8.192 MHz output sample rate */
#define CN0_TARGET_DBHZ   42.0        /* 40–45 dB-Hz 常用 */
#define HEADROOM_RATIO    0.80        /* 約 -2 dB 頭房，防飽和 */
#define AMP_SMOOTH_TC_MS  1000        /* 振幅平滑時間常數 */
/* Physical constant */
#define CLIGHT            299792458.0 /* Speed of light (m/s) */
/* Per-orbit-class amplitude offsets (dB) */
#define GAIN_MEO_DB   +1.5
#define GAIN_IGSO_DB  +1.5

/* 由 globals.c 定義 */
extern uint8_t      prn_code[MAX_SAT][CODE_LEN];
extern ephemeris_t  eph[MAX_SAT];
extern int          simulator_inited;
extern int          prn_max;        /* 上限 PRN 編號 */
extern double       nav_time_min;
extern double       nav_time_max;
extern int          nav_week;       /* ephemeris week reference */
extern double       iono_alpha[4];  /* ionospheric parameters */
extern double       iono_beta[4];
extern int          utc_bdt_diff;   /* UTC->BDT offset seconds */

/* Global transmit time in BDT seconds */
extern double       g_t_tx;

#endif

