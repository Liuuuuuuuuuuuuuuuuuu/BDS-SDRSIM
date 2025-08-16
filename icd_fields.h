#ifndef ICD_FIELDS_H
#define ICD_FIELDS_H
#include <stdint.h>

/* BeiDou B1I D1 frame fields extracted from ICD.
 * Each member width is noted in bits for reference. */
typedef struct {
    uint32_t sow;        /* 20 bits */
    uint32_t satH1;      /* 1 bit  */
    uint32_t aodc;       /* 5 bits */
    uint32_t urai;       /* 4 bits */
    uint32_t toc;        /* 17 bits */
    int32_t  a0;         /* 24 bits signed */
    int32_t  a1;         /* 22 bits signed */
    int32_t  a2;         /* 11 bits signed */
    int32_t  tgd1;       /* 10 bits signed */
    int32_t  tgd2;       /* 10 bits signed */
    int32_t  alpha[4];   /* 8 bits each signed */
    int32_t  beta[4];    /* 8 bits each signed */

    /* orbital parameters --------------------------------------- */
    int32_t  delta_n;    /* 16 bits signed */
    int32_t  M0;         /* 32 bits signed */
    uint32_t e;          /* 32 bits unsigned */
    uint32_t sqrtA;      /* 32 bits unsigned */
    int32_t  cuc;        /* 18 bits signed */
    int32_t  cus;        /* 18 bits signed */
    int32_t  crc;        /* 18 bits signed */
    int32_t  crs;        /* 18 bits signed */
    uint32_t toe;        /* 17 bits */
    uint32_t aode;       /* 5 bits */
    int32_t  i0;         /* 32 bits signed */
    int32_t  omega0;     /* 32 bits signed */
    int32_t  omega;      /* 32 bits signed */
    int32_t  omegadot;   /* 24 bits signed */
    int32_t  idot;       /* 14 bits signed */
    int32_t  cic;        /* 18 bits signed */
    int32_t  cis;        /* 18 bits signed */
} B1I_D1_Frame;


#endif /* ICD_FIELDS_H */
