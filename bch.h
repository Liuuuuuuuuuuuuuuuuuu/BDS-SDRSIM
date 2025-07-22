#ifndef BCH_H
#define BCH_H
#include <stdint.h>
/* 將 11 bit 資訊位編成 15 bit BCH(15,11,1) */
uint16_t bch_encode(uint16_t info);
/* 將 26 位元資料尾端 11 位做 BCH，輸出 30 位 */
uint32_t bch_encode_26bit(uint32_t payload);
/* 將兩組 11 位元資料分別 BCH 後交錯成 30 位 */
uint32_t bch_interleave_22bit(uint32_t payload);
#endif

