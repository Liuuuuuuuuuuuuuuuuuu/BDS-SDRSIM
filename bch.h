#ifndef BCH_H
#define BCH_H
#include <stdint.h>
/* 將 11 bit 資訊位編成 15 bit BCH(15,11,1) */
uint16_t bch_encode(uint16_t info);
#endif

