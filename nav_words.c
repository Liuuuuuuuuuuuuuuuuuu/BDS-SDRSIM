#include <stdint.h>
#include "bch.h"
#include "nav_words.h"

/* Build a 30-bit navigation word from payload bits.
 * payload: left aligned (MSB first) within 'bits' width.
 * bits    : number of information bits (22 or 26).
 */
uint32_t build_word(uint32_t payload, int bits)
{
    uint32_t word = payload << (30 - bits);
    uint32_t parity = bch1511(payload, bits);
    word |= parity & ((bits==26)?0xF:0xFF);
    return word;
}
