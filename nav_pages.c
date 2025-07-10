#include <stdint.h>
#include <string.h>
#include "icd_fields.h"
#include "nav_words.h"

static void sf1(const B1I_D1_Frame *f, int week, uint32_t sow, uint32_t w[10])
{
    uint32_t p;
    /* Word1: preamble + reserved + sfID + SOW[19:12] */
    p = 0x712;             /* 11 bits */
    p = (p << 4);          /* reserved */
    p = (p << 3) | 1;      /* FraID */
    p = (p << 8) | ((sow >> 12) & 0xFF);
    w[0] = build_word(p, 26);

    /* Word2: SOW[11:0] + week */
    p = ((sow & 0xFFF) << 13) | (week & 0x1FFF);
    w[1] = build_word(p, 26);

    /* Word3: URAI, SatH1, toc, AODC */
    p = (f->urai & 0xF) << 26;
    p |= (f->satH1 & 0x1) << 25;
    p |= (f->toc & 0x1FFFF) << 8;
    p |= (f->aodc & 0x1F);
    w[2] = build_word(p, 26);

    /* Word4: a0 (24 bits signed) */
    p = ((uint32_t)f->a0 & 0xFFFFFF);
    w[3] = build_word(p, 26);

    /* Word5: a1(22) + a2(11) */
    p = ((uint32_t)f->a1 & 0x3FFFFF) << 11;
    p |= ((uint32_t)f->a2 & 0x7FF);
    w[4] = build_word(p, 22);

    /* repeat words */
    for (int i = 0; i < 5; ++i) w[i + 5] = w[i];
}

static void sf2(const B1I_D1_Frame *f, uint32_t w[10])
{
    uint32_t p;
    /* Word1: sync + FraID=2 + toe[16:9] */
    p = 0x712;
    p = (p << 4);
    p = (p << 3) | 2;
    p = (p << 8) | ((f->toe >> 9) & 0xFF);
    w[0] = build_word(p, 26);

    /* Word2: toe[8:0] + sqrtA[31:23] */
    p = ((f->toe & 0x1FF) << 13) | ((uint32_t)f->sqrtA >> 19 & 0x1FFF);
    w[1] = build_word(p, 26);

    /* Word3: sqrtA[22:0] + e[31:21] */
    p = ((uint32_t)f->sqrtA & 0x7FFFFF) << 11 | ((uint32_t)f->e >> 21 & 0x7FF);
    w[2] = build_word(p, 26);

    /* Word4: e[20:0] + delta_n[15:0] + M0[31] */
    p = ((uint32_t)f->e & 0x1FFFFF) << 9 | ((uint32_t)f->delta_n & 0xFFFF) << 3 | ((uint32_t)f->M0 >> 29 & 0x7);
    w[3] = build_word(p, 26);

    /* Word5: M0[30:10] */
    p = ((uint32_t)f->M0 >> 8) & 0x3FFFFF;
    w[4] = build_word(p, 22);

    for (int i = 0; i < 5; ++i) w[i + 5] = w[i];
}

static void sf3(const B1I_D1_Frame *f, uint32_t w[10])
{
    uint32_t p;
    /* Word1: sync + FraID=3 + omega0[31:24] */
    p = 0x712;
    p = (p << 4);
    p = (p << 3) | 3;
    p = (p << 8) | ((uint32_t)f->omega0 >> 24 & 0xFF);
    w[0] = build_word(p, 26);

    /* Word2: omega0[23:10] + i0[31:23] */
    p = ((uint32_t)f->omega0 >> 10 & 0x3FFF) << 13 | ((uint32_t)f->i0 >> 23 & 0x1FFF);
    w[1] = build_word(p, 26);

    /* Word3: i0[22:0] + omega[31:21] */
    p = ((uint32_t)f->i0 & 0x7FFFFF) << 11 | ((uint32_t)f->omega >> 21 & 0x7FF);
    w[2] = build_word(p, 26);

    /* Word4: omega[20:0] + crc[17:0] + idot[13] */
    p = ((uint32_t)f->omega & 0x1FFFFF) << 9 | ((uint32_t)f->crc & 0x3FFFF) << 1 | ((uint32_t)f->idot >> 13 & 0x1);
    w[3] = build_word(p, 26);

    /* Word5: idot[12:0] + omegadot[23:0] */
    p = ((uint32_t)f->idot & 0x1FFF) << 20 | ((uint32_t)f->omegadot & 0xFFFFF);
    w[4] = build_word(p, 22);

    for (int i = 0; i < 5; ++i) w[i + 5] = w[i];
}

void assemble_d1_subframe(int sfid, const B1I_D1_Frame *frm,
                          int week, uint32_t sow, uint32_t words[10])
{
    memset(words, 0, sizeof(uint32_t)*10);
    if (!frm) return;

    switch (sfid) {
    case 1: sf1(frm, week, sow, words); break;
    case 2: sf2(frm, words);           break;
    case 3: sf3(frm, words);           break;
    default: break;
    }
}

void assemble_d1_subframe2(const B1I_D1_Frame *frm, uint32_t words[10])
{
    if (frm) sf2(frm, words); else memset(words,0,sizeof(uint32_t)*10);
}

void assemble_d1_subframe3(const B1I_D1_Frame *frm, uint32_t words[10])
{
    if (frm) sf3(frm, words); else memset(words,0,sizeof(uint32_t)*10);
}
