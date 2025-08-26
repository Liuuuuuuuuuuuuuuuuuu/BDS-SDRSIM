#include "bch.h"

uint16_t bch_encode(uint16_t d11)
{
    d11 &= 0x7FF;
    uint16_t reg = d11 << 4;          /* 高 4 位空出校驗 */
    for(int i=14;i>=4;i--){
        if(reg & (1<<i))
            reg ^= BCH_G << (i-4);
    }
    return (d11<<4) | (reg & 0xF);
}

/* --------------------------- 對 26 位元資料的後 11 位做 BCH，並在尾端附上 4 位校驗 ------------------------------*/
uint32_t bch_encode_26bit(uint32_t payload)
{
    payload &= 0x3FFFFFF;                 /* 26 bits */
    uint16_t info = payload & 0x7FF;      /* 取出最後 11 位 */
    uint16_t parity = bch_encode(info) & 0xF;
    return (payload << 4) | parity;       /* 26 資料 + 4 校驗 */
}

/* --------------------------- 將 22 位元資料分成兩組 11 位做 BCH，並交錯輸出 ------------------------------*/
uint32_t bch_interleave_22bit(uint32_t payload)
{
    uint16_t high11 = (payload >> 11) & 0x7FF;
    uint16_t low11  = payload & 0x7FF;
    uint16_t codeA = bch_encode(high11);
    uint16_t codeB = bch_encode(low11);
    uint32_t result = 0;
    for (int i = 14; i >= 0; --i) {
        result = (result << 1) | ((codeA >> i) & 1);
        result = (result << 1) | ((codeB >> i) & 1);
    }
    return result;
}
/* ---------------------------  End  ------------------------------*/
