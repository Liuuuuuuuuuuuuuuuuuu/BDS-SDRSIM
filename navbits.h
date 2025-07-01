#ifndef NAVBITS_H
#define NAVBITS_H

#include <stdint.h>
#include "bdssim.h"

#define SUBFRAME_BITS (300)          /* 10 words × 30 bits */

#define HALF_SUBFRAME_BITS (SUBFRAME_BITS/2)

#define SF_STREAM_LEN (SUBFRAME_BITS)

void navbits_init(void);
void get_subframe_bits(int prn, int sf_id, int week, double sow, uint8_t *out); /* sf_id=1,2,3 */

#endif






