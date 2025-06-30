#ifndef NAVBITS_H
#define NAVBITS_H

#include <stdint.h>
#include "bdssim.h"

#define SUBFRAME_BITS (300)          /* 10 words × 30 bits */
#define SF_STREAM_LEN (SUBFRAME_BITS)

void navbits_init(int bdt_week, double sow);
void get_subframe_bits(int prn, int sf_id, uint8_t *out); /* sf_id=1,2,3 */

#endif






