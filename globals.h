/* globals.h */
#ifndef GLOBALS_H
#define GLOBALS_H

#include <stdint.h>

/* Some platforms do not define M_PI in math.h */
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

/* Physical and signal constants */
#define CLIGHT      299792458.0    /* Speed of light (m/s) */
#define F_B1I       1561.098e6     /* B1I carrier frequency (Hz) */
#define CHIPRATE    2.046e6        /* B1I code chipping rate (Hz) */

/* Simulation limits */
#define MAX_CH      8              /* Maximum channels */
#define CODE_LEN    2046           /* PRN code length */
#define MAX_SAT     65             /* PRN slots */

/* Amplitude and sampling parameters */
#define FS_OUTPUT_HZ    5200000.0
#define CN0_TARGET_DBHZ 42.0
#define HEADROOM_RATIO  0.80
#define AMP_SMOOTH_TC_MS 1000
#define GAIN_MEO_DB     +0.5
#define GAIN_IGSO_DB    +1.5

/* G2 tap table for B1I PRN generation */
static const uint8_t g2_taps[64][3] = {
 {0,0,0},
 {1,3,0},{1,4,0},{1,5,0},{1,6,0},{1,8,0},
 {1,9,0},{1,10,0},{1,11,0},{2,7,0},{3,4,0},
 {3,5,0},{3,6,0},{3,8,0},{3,9,0},{3,10,0},{3,11,0},
 {4,5,0},{4,6,0},{4,8,0},{4,9,0},{4,10,0},{4,11,0},
 {5,6,0},{5,8,0},{5,9,0},{5,10,0},{5,11,0},
 {6,8,0},{6,9,0},{6,10,0},{6,11,0},
 {8,9,0},{8,10,0},{8,11,0},
 {9,10,0},{9,11,0},{10,11,0},
 {1,2,7},{1,3,4},{1,3,6},{1,3,8},{1,3,10},
 {1,3,11},{1,4,5},{1,4,9},{1,5,6},{1,5,8},
 {1,5,10},{1,5,11},{1,6,9},{1,8,9},{1,9,10},{1,9,11},
 {2,3,7},{2,5,7},{2,7,9},{3,4,5},{3,4,9},
 {3,5,6},{3,5,8},{3,5,10},{3,5,11},{3,6,9}
};

extern int simulator_inited;
extern int prn_max;
extern double nav_time_min;
extern double nav_time_max;
extern int nav_week;
extern double iono_alpha[4];
extern double iono_beta[4];
extern int utc_bdt_diff;
extern double g_t_tx;
extern double g_target_cn0;
extern int g_enable_iono;

#endif /* GLOBALS_H */

