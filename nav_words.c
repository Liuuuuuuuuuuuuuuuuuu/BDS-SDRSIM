#include <stdint.h>
#include "bch.h"
#include "nav_words.h"

/* Build a 30-bit navigation word following the BeiDou ICD.
 * bits must be either 26 (first word with a 15-bit unencoded prefix
 * and 11-bit BCH block) or 22 (two interleaved BCH blocks).
 */
uint32_t build_word(uint32_t payload, int bits)
{
    if (bits == 22) {
        /* split into two 11-bit groups and interleave the BCH codes */
        uint16_t infoA = (payload >> 11) & 0x7FF;
        uint16_t infoB = payload & 0x7FF;
        uint16_t codeA = bch_encode(infoA);
        uint16_t codeB = bch_encode(infoB);
        uint32_t w = 0;
        for (int i = 14; i >= 0; --i) {
            w = (w << 1) | ((codeA >> i) & 1);
            w = (w << 1) | ((codeB >> i) & 1);
        }
        return w;
    } else if (bits == 26) {
        /* 15-bit prefix followed by one BCH(15,11,1) block */
        uint32_t prefix = payload >> 11;         /* upper 15 bits */
        uint16_t info   = payload & 0x7FF;       /* lower 11 bits */
        uint16_t code   = bch_encode(info);
        return (prefix << 15) | code;
    }

    /* fallback: behave like the previous implementation */
    uint32_t word = payload << (30 - bits);
    uint32_t parity = bch1511(payload, bits);
    word |= parity & ((bits==26)?0xF:0xFF);
    return word;
}
