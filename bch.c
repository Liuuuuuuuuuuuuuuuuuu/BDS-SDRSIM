#include "bch.h"

/* 生成多項式 g(x)=1+x+x⁴ = 0b1_0011 = 0x13 */
static const uint16_t G = 0x13;

uint16_t bch_encode(uint16_t d11)
{
    d11 &= 0x7FF;
    uint16_t reg = d11 << 4;          /* 高 4 位空出校驗 */
    for(int i=14;i>=4;i--){
        if(reg & (1<<i))
            reg ^= G << (i-4);
    }
    return (d11<<4) | (reg & 0xF);
}



/* Generic BCH(15,11,1) encoder returning parity bits for 22 or 26 bit payloads */
uint32_t bch1511(uint32_t in, int payloadBits)
{
    if (payloadBits != 22 && payloadBits != 26)
        return 0;
    int parityBits = (payloadBits == 26) ? 4 : 8;
    uint32_t reg = in << parityBits;  /* space for parity */
    int top = payloadBits + parityBits - 1;
    for (int i = top; i >= parityBits; --i)
        if (reg & (1u << i))
            reg ^= (uint32_t)G << (i - 4);
    return reg & ((1u << parityBits) - 1);
}
