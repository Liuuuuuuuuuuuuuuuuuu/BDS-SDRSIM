#ifndef GLOBALS_H
#define GLOBALS_H
#define DEBUG
#include <stdint.h>
#include "bdssim.h"          /* 已含 MAX_SAT、CODE_LEN、ephemeris_t */

/* 由 globals.c 定義 */
extern uint8_t      prn_code[MAX_SAT][CODE_LEN];
extern ephemeris_t  eph[MAX_SAT];
extern int          simulator_inited;

#endif

