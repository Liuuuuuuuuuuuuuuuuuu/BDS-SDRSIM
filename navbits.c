#include "navbits.h"
#include "bch.h"
#include "globals.h"
#include "icd_fields.h"
#include <math.h>
#include <string.h>

extern ephemeris_t eph[MAX_SAT];
extern double iono_alpha[4], iono_beta[4];

static inline double rad_to_sc(double v) { return v / M_PI; }
static uint32_t build_word(uint32_t p, int b) {
  return b == 22 ? bch_interleave_22bit(p) : b == 26 ? bch_encode_26bit(p) : 0;
}
static void put_word(uint8_t *b, int pos, uint32_t w) {
  for (int i = 0; i < 30; i++)
    b[pos + i] = (w >> (29 - i)) & 1;
}

static void frame_from_ephemeris(const ephemeris_t *e, B1I_D1_Frame *f) {
  memset(f, 0, sizeof(*f));
  f->toc = (uint32_t)(e->toc / 8);
  f->aodc = e->aodc & 0x1F;
  f->urai = e->ura & 0xF;
  f->satH1 = 0;
  f->aode = e->aode & 0x1F;
  f->tgd1 = (int32_t)llround(e->tgd1 / 1e-10);
  f->tgd2 = (int32_t)llround(e->tgd2 / 1e-10);
  static const double as[4] = {pow(2, -30), pow(2, -27), pow(2, -24),
                               pow(2, -24)};
  static const double bs[4] = {pow(2, 11), pow(2, 14), pow(2, 16), pow(2, 16)};
  for (int i = 0; i < 4; i++) {
    f->alpha[i] = (int32_t)llround(iono_alpha[i] / as[i]);
    f->beta[i] = (int32_t)llround(iono_beta[i] / bs[i]);
  }
  f->a0 = (int32_t)llround(e->af0 / pow(2, -33));
  f->a1 = (int32_t)llround(e->af1 / pow(2, -50));
  f->a2 = (int32_t)llround(e->af2 / pow(2, -66));
  f->toe = (uint32_t)(e->toe / 8);
  f->sqrtA = (uint32_t)llround(e->sqrtA / pow(2, -19));
  f->e = (uint32_t)llround(e->e / pow(2, -33));
  f->delta_n = (int32_t)llround(rad_to_sc(e->deltan) / pow(2, -43));
  f->M0 = (int32_t)llround(rad_to_sc(e->M0) / pow(2, -31));
  f->omega0 = (int32_t)llround(rad_to_sc(e->omega0) / pow(2, -31));
  f->i0 = (int32_t)llround(rad_to_sc(e->i0) / pow(2, -31));
  f->omega = (int32_t)llround(rad_to_sc(e->w) / pow(2, -31));
  f->crc = (int32_t)llround(e->crc / pow(2, -6));
  f->crs = (int32_t)llround(e->crs / pow(2, -6));
  f->cuc = (int32_t)llround(e->cuc / pow(2, -31));
  f->cus = (int32_t)llround(e->cus / pow(2, -31));
  f->cic = (int32_t)llround(e->cic / pow(2, -31));
  f->cis = (int32_t)llround(e->cis / pow(2, -31));
  f->idot = (int32_t)llround(rad_to_sc(e->idot) / pow(2, -43));
  f->omegadot = (int32_t)llround(rad_to_sc(e->omegadot) / pow(2, -43));
}

void navbits_init(void) {}

void get_subframe_bits(int prn, int sf, int week, double sow, double len,
                       uint8_t *out) {
  B1I_D1_Frame f;
  frame_from_ephemeris(&eph[prn], &f);
  double start = floor(sow / len) * len;
  uint32_t sowi = (uint32_t)start;
  memset(out, 0, SF_STREAM_LEN);
  switch (sf) {
  case 1: {
    uint32_t w = 0x712;
    w = (w << 4);
    w = (w << 3) | 1;
    w = (w << 8) | ((sowi >> 12) & 0xFF);
    put_word(out, 0, build_word(w, 26));
    w = ((sowi & 0xFFF) << 10) | ((f.satH1 & 1) << 9) | ((f.aodc & 0x1F) << 4) |
        (f.urai & 0xF);
    put_word(out, 30, build_word(w, 22));
    w = ((week & 0x1FFF) << 9) | ((f.toc >> 8) & 0x1FF);
    put_word(out, 60, build_word(w, 22));
    w = ((f.toc & 0xFF) << 14) | ((f.tgd1 & 0x3FF) << 4) |
        ((f.tgd2 >> 6) & 0xF);
    put_word(out, 90, build_word(w, 22));
    w = ((f.tgd2 & 0x3F) << 16) | ((f.alpha[0] & 0xFF) << 8) |
        (f.alpha[1] & 0xFF);
    put_word(out, 120, build_word(w, 22));
    uint32_t b0 = f.beta[0], b1 = f.beta[1], b2 = f.beta[2], b3 = f.beta[3];
    w = ((f.alpha[2] & 0xFF) << 14) | ((f.alpha[3] & 0xFF) << 6) |
        ((b0 >> 2) & 0x3F);
    put_word(out, 150, build_word(w, 22));
    w = ((b0 & 3) << 20) | ((b1 & 0xFF) << 12) | ((b2 & 0xFF) << 4) |
        ((b3 >> 4) & 0xF);
    put_word(out, 180, build_word(w, 22));
    int32_t a0 = f.a0, a1 = f.a1, a2 = f.a2;
    w = ((b3 & 0xF) << 18) | ((a2 & 0x7FF) << 7) |
        (((uint32_t)a0 >> 17) & 0x7F);
    put_word(out, 210, build_word(w, 22));
    w = (((uint32_t)a0 & 0x1FFFF) << 5) | (((uint32_t)a1 >> 17) & 0x1F);
    put_word(out, 240, build_word(w, 22));
    w = (((uint32_t)a1 & 0x1FFFF) << 5) | (f.aode & 0x1F);
    put_word(out, 270, build_word(w, 22));
    break;
  }
  case 2: {
    uint32_t w = 0x712;
    w = (w << 4);
    w = (w << 3) | 2;
    w = (w << 8) | ((sowi >> 12) & 0xFF);
    put_word(out, 0, build_word(w, 26));
    w = ((sowi & 0xFFF) << 10) | ((uint32_t)f.delta_n >> 6 & 0x3FF);
    put_word(out, 30, build_word(w, 22));
    w = (((uint32_t)f.delta_n & 0x3F) << 16) |
        (((uint32_t)f.cuc >> 2) & 0xFFFF);
    put_word(out, 60, build_word(w, 22));
    w = (((uint32_t)f.cuc & 3) << 20) | (((uint32_t)f.M0 >> 12) & 0xFFFFF);
    put_word(out, 90, build_word(w, 22));
    w = (((uint32_t)f.M0 & 0xFFF) << 10) | ((uint32_t)f.e >> 22 & 0x3FF);
    put_word(out, 120, build_word(w, 22));
    w = ((uint32_t)f.e & 0x3FFFFF);
    put_word(out, 150, build_word(w, 22));
    w = (((uint32_t)f.cus & 0x3FFFF) << 4) | ((uint32_t)f.crc >> 14 & 0xF);
    put_word(out, 180, build_word(w, 22));
    w = (((uint32_t)f.crc & 0x3FFF) << 8) | ((uint32_t)f.crs >> 10 & 0xFF);
    put_word(out, 210, build_word(w, 22));
    w = (((uint32_t)f.crs & 0x3FF) << 12) | ((uint32_t)f.sqrtA >> 20 & 0xFFF);
    put_word(out, 240, build_word(w, 22));
    w = (((uint32_t)f.sqrtA & 0xFFFFF) << 2) | ((f.toe >> 15) & 3);
    put_word(out, 270, build_word(w, 22));
    break;
  }
  case 3: {
    uint32_t w = 0x712;
    w = (w << 4);
    w = (w << 3) | 3;
    w = (w << 8) | ((sowi >> 12) & 0xFF);
    put_word(out, 0, build_word(w, 26));
    w = ((sowi & 0xFFF) << 10) | ((f.toe >> 5) & 0x3FF);
    put_word(out, 30, build_word(w, 22));
    w = (((f.toe & 0x1F) << 17) | ((uint32_t)f.i0 >> 15 & 0x1FFFF));
    put_word(out, 60, build_word(w, 22));
    w = (((uint32_t)f.i0 & 0x7FFF) << 7) | ((uint32_t)f.cic >> 11 & 0x7F);
    put_word(out, 90, build_word(w, 22));
    w = (((uint32_t)f.cic & 0x7FF) << 11) |
        ((uint32_t)f.omegadot >> 13 & 0x7FF);
    put_word(out, 120, build_word(w, 22));
    w = (((uint32_t)f.omegadot & 0x1FFF) << 9) | ((uint32_t)f.cis >> 9 & 0x1FF);
    put_word(out, 150, build_word(w, 22));
    w = (((uint32_t)f.cis & 0x1FF) << 13) | ((uint32_t)f.idot >> 1 & 0x1FFF);
    put_word(out, 180, build_word(w, 22));
    w = (((uint32_t)f.idot & 1) << 21) | ((uint32_t)f.omega0 >> 11 & 0x1FFFFF);
    put_word(out, 210, build_word(w, 22));
    w = (((uint32_t)f.omega0 & 0x7FF) << 11) |
        ((uint32_t)f.omega >> 21 & 0x7FF);
    put_word(out, 240, build_word(w, 22));
    w = (((uint32_t)f.omega & 0x1FFFFF) << 1);
    put_word(out, 270, build_word(w, 22));
    break;
  }
  case 4: {
    uint32_t w = 0x712;
    w = (w << 4);
    w = (w << 3) | 4;
    w = (w << 8) | ((sowi >> 12) & 0xFF);
    put_word(out, 0, build_word(w, 26));

    int frame = (int)floor(sow / (len * 5.0));
    int pnum = frame % 24 + 1;

    w = ((sowi & 0xFFF) << 10) | ((pnum & 0x7F) << 2) | 0x2;
    put_word(out, 30, build_word(w, 22));

    uint32_t filler = 0x2AAAAA;
    for (int i = 2; i < 10; i++)
      put_word(out, i * 30, build_word(filler, 22));
    break;
  }
  case 5: {
    uint32_t w = 0x712;
    w = (w << 4);
    w = (w << 3) | 5;
    w = (w << 8) | ((sowi >> 12) & 0xFF);
    put_word(out, 0, build_word(w, 26));

    int frame = (int)floor(sow / (len * 5.0));
    int pnum = frame % 24 + 1;

    w = ((sowi & 0xFFF) << 10) | ((pnum & 0x7F) << 2) | 0x2;
    put_word(out, 30, build_word(w, 22));

    uint32_t filler = 0x2AAAAA;
    for (int i = 2; i < 10; i++)
      put_word(out, i * 30, build_word(filler, 22));
    break;
  }
  default:
    break;
  }
}
