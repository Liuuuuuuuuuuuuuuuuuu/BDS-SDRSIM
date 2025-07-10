#ifndef NAV_WORDS_H
#define NAV_WORDS_H
#include <stdint.h>
/* Construct a 30‑bit navigation word following the BeiDou ICD.
 * For bits == 26 the upper 15 bits are transmitted directly and the
 * lower 11 bits are BCH encoded.  For bits == 22 two 11‑bit groups are
 * BCH encoded and interleaved bit by bit.
 */
uint32_t build_word(uint32_t payload, int bits);
#endif
