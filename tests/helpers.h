#ifndef TEST_HELPERS_H
#define TEST_HELPERS_H

#include <stdint.h>
#include <stddef.h>

/* Load minimal RINEX navigation data for tests */
int load_rinex(const char *path);

/* Generate PRN sequence for given PRN id into buf (length CODE_LEN) */
void generate_prn(int prn, uint8_t *buf, size_t len);

/* Compute Hamming(32,26) parity for info bits */
uint32_t hamming_parity(uint32_t info);

#endif
