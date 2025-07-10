#ifndef BCH_H
#define BCH_H
#include <stdint.h>
/* 將 11 bit 資訊位編成 15 bit BCH(15,11,1) */
uint16_t bch_encode(uint16_t info);
/* 通用 BCH(15,11,1) for longer payloads (22 or 26 bits). */
uint32_t bch1511(uint32_t in, int payloadBits);
#endif

