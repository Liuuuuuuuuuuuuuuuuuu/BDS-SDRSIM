#ifndef RINEX_H
#define RINEX_H

#include <stdint.h>

typedef struct ephemeris {
    int    prn, week;
    double toe;
    double sqrtA, e, i0, omega0, w, M0;
    double deltan, idot, omegadot;
    double cuc, cus, cic, cis, crc, crs;
    double af0, af1, af2;
    int     aode;
    double  tgd1, tgd2;
    int     aodc;
    int     toc;
    uint8_t ura;
    uint8_t health;
} ephemeris_t;

int read_rinex_nav(const char *file, double start_bdt);

#endif /* RINEX_H */


