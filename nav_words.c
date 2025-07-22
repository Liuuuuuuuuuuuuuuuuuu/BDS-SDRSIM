#include <stdint.h>
#include "bch.h"
#include "nav_words.h"

/* Build a 30-bit navigation word following the BeiDou ICD.
 * bits must be either 26 (first word with a 15-bit unencoded prefix
 * and 11-bit BCH block) or 22 (two interleaved BCH blocks).
 */
uint32_t build_word(uint32_t payload, int bits)
{
    if (bits == 22)
        return bch_interleave_22bit(payload);
    else if (bits == 26)
        return bch_encode_26bit(payload);
    return 0;
}
