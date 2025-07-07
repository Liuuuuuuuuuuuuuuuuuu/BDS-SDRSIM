#ifndef TEST_HELPERS_H
#define TEST_HELPERS_H

#include <stdint.h>
#include <stddef.h>

/* Load RINEX navigation file from tests/vectors/ */
int load_rinex(const char *fname);

/* Ensure PRN table is generated once for tests */
void ensure_prn(void);

/* Generate PRN sequence for given PRN id into buf (length CODE_LEN) */
void generate_prn(int prn, uint8_t *buf, size_t len);

/* Compute Hamming(32,26) parity for info bits */
uint32_t hamming_parity(uint32_t info);

#endif
