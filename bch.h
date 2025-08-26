#ifndef BCH_H
#define BCH_H

#include <stdint.h>

/* generator polynomial g(x)=1+x+x^4 */
static const uint16_t BCH_G = 0x13;

uint16_t bch_encode(uint16_t info);
uint32_t bch_encode_26bit(uint32_t payload);
uint32_t bch_interleave_22bit(uint32_t payload);

#endif /* BCH_H */

