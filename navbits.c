/* navbits.c : 產生 B1I D1/D2 子帧 (星曆) */

#include <math.h>
#include <string.h>
#include "navbits.h"
#include "bch.h"

/* Construct a 30-bit navigation word. "bits" is either 26 for the
 * first word type or 22 for the remaining words. */
static uint32_t build_word(uint32_t payload, int bits)
{
    if(bits == 22)
        return bch_interleave_22bit(payload);
    else if(bits == 26)
        return bch_encode_26bit(payload);
    return 0;
}

extern ephemeris_t eph[MAX_SAT];
#include "icd_fields.h"

/* Convert ephemeris_t to the simplified B1I_D1_Frame structure.
 * Only the fields needed for navigation message assembly are
 * populated here. */
static void frame_from_ephemeris(const ephemeris_t *e, B1I_D1_Frame *f)
{
    memset(f, 0, sizeof(*f));
    f->toc      = (uint32_t)(e->toc / 8);
    f->aodc     = e->aodc & 0x1F;
    f->urai     = e->ura & 0xF;
    f->satH1    = e->health & 0x1;

    /* TGD values: 0.1 ns resolution (ICD 5.2.4.10) */
    f->tgd1     = (int32_t)llround(e->tgd1 / 1e-10);
    f->tgd2     = (int32_t)llround(e->tgd2 / 1e-10);

    /* Ionospheric parameters not parsed – leave as zero */
    for (int i = 0; i < 4; ++i) {
        f->alpha[i] = 0;
        f->beta[i]  = 0;
    }

    f->a0       = (int32_t)llround(e->af0 / pow(2, -33));
    f->a1       = (int32_t)llround(e->af1 / pow(2, -50));
    f->a2       = (int32_t)llround(e->af2 / pow(2, -66));

    f->toe      = (uint32_t)(e->toe / 8);
    f->sqrtA    = (uint32_t)llround(e->sqrtA / pow(2, -19));
    f->e        = (uint32_t)llround(e->e / pow(2, -33));
    f->delta_n  = (int32_t)llround(e->deltan / pow(2, -43));
    f->M0       = (int32_t)llround(e->M0 / pow(2, -31));

    f->omega0   = (int32_t)llround(e->omega0 / pow(2, -31));
    f->i0       = (int32_t)llround(e->i0 / pow(2, -31));
    f->omega    = (int32_t)llround(e->w / pow(2, -31));
    f->crc      = (int32_t)llround(e->crc / pow(2, -6));
    f->crs      = (int32_t)llround(e->crs / pow(2, -6));
    f->cuc      = (int32_t)llround(e->cuc / pow(2, -31));
    f->cus      = (int32_t)llround(e->cus / pow(2, -31));
    f->cic      = (int32_t)llround(e->cic / pow(2, -31));
    f->cis      = (int32_t)llround(e->cis / pow(2, -31));
    f->idot     = (int32_t)llround(e->idot / pow(2, -43));
    f->omegadot = (int32_t)llround(e->omegadot / pow(2, -43));
}

/* --------------------------------- 宏 & 工具 -------------------------------- */

/* 填 30-bit word 至 bit 流。pos 為 bit index 0..299 (MSB first) */
static void put_word(uint8_t *b, int pos, uint32_t w30)
{
    for(int i=0;i<30;i++)
        b[pos+i] = (w30>>(29-i)) & 1;
}

/* --------------------------------- 子帧 1 ----------------------------------- */
static void build_subframe1_d1(uint8_t *out, const ephemeris_t *e,
                               int week, double sow, double frame_len)
{
    memset(out,0,SF_STREAM_LEN);

    /* word1：帧同步 11 bits（11100010010） + FraID(001) + SOW[19:12] */
    uint32_t info = 0x712;            /* 11-bit preamble */
    info = (info << 4);               /* reserved bits */
    info = (info << 3) | 0x1;         /* FraID=001 */
    /*
     * SOW should reflect the start time of this subframe.
     * According to the B1I ICD section 5.2.4.3 the value
     * corresponds to the rising edge of the frame sync
     * at the beginning of the subframe.
     */
    /* For D2 the message repeats every 0.6 s whereas D1 repeats
       every 6 s. The time-of-week field shall therefore reflect the
       start of the current subframe length. */
    uint32_t sow_int = (uint32_t)(floor(sow/frame_len)*frame_len);
    info = (info<<8) | ((sow_int>>12) & 0xFF);
    put_word(out,0, build_word(info, 26));

    /* word2: SOW[11:0] + SatH1 + AODC + URAI */
    uint32_t w2 = ((sow_int & 0xFFF) << 10) |
                  ((e->health & 1) << 9) |
                  ((e->aodc & 0x1F) << 4) |
                  (e->ura & 0xF);
    put_word(out,30, build_word(w2,22));

    /* word3: WN + toc high 9 bits */
    uint32_t toc_ticks = (uint32_t)(e->toc/8);
    uint32_t w3 = ((week & 0x1FFF) << 9) |
                  ((toc_ticks >> 8) & 0x1FF);
    put_word(out,60, build_word(w3,22));

    /* word4: toc low 8 + TGD1 + TGD2 high 4 */
    int32_t tgd1_i = (int32_t)llround(e->tgd1/1e-10);
    int32_t tgd2_i = (int32_t)llround(e->tgd2/1e-10);
    uint32_t w4 = ((toc_ticks & 0xFF) << 14) |
                  ((tgd1_i & 0x3FF) << 4) |
                  ((tgd2_i >> 6) & 0xF);
    put_word(out,90, build_word(w4,22));

    /* word5: TGD2 low 6 + alpha0 + alpha1 */
    uint32_t w5 = ((tgd2_i & 0x3F) << 16) |
                  ((0 & 0xFF) << 8) |
                  (0 & 0xFF);
    put_word(out,120, build_word(w5,22));

    /* word6: alpha2 + alpha3 + beta0 high 6 */
    uint32_t beta0 = 0;
    uint32_t w6 = ((0 & 0xFF) << 14) |
                  ((0 & 0xFF) << 6) |
                  ((beta0 >> 2) & 0x3F);
    put_word(out,150, build_word(w6,22));

    /* word7: beta0 low 2 + beta1 + beta2 + beta3 high 4 */
    uint32_t beta1 = 0, beta2 = 0, beta3 = 0;
    uint32_t w7 = ((beta0 & 0x3) << 20) |
                  ((beta1 & 0xFF) << 12) |
                  ((beta2 & 0xFF) << 4) |
                  ((beta3 >> 4) & 0xF);
    put_word(out,180, build_word(w7,22));

    /* word8: beta3 low 4 + a2 + a0 high 7 */
    int32_t a0_i = (int32_t)llround(e->af0 / pow(2,-33));
    int32_t a1_i = (int32_t)llround(e->af1 / pow(2,-50));
    int32_t a2_i = (int32_t)llround(e->af2 / pow(2,-66));
    uint32_t w8 = ((beta3 & 0xF) << 18) |
                  ((a2_i & 0x7FF) << 7) |
                  (((uint32_t)a0_i >> 17) & 0x7F);
    put_word(out,210, build_word(w8,22));

    /* word9: a0 low 17 + a1 high 5 */
    uint32_t w9 = ((uint32_t)a0_i & 0x1FFFF) << 5 |
                  (((uint32_t)a1_i >> 17) & 0x1F);
    put_word(out,240, build_word(w9,22));

    /* word10: a1 low 17 + AODE */
    uint32_t w10 = (((uint32_t)a1_i & 0x1FFFF) << 5) |
                   (e->aode & 0x1F);
    put_word(out,270, build_word(w10,22));
}

/* --------------------------------- 子帧 2 ----------------------------------- */
static void build_subframe2_d1(const B1I_D1_Frame *f, uint32_t sow, uint32_t w[10])
{
    uint32_t p;
    /* Word1: preamble + reserved + FraID=2 + SOW[19:12] */
    p = 0x712;
    p = (p << 4);
    p = (p << 3) | 2;
    p = (p << 8) | ((sow >> 12) & 0xFF);
    w[0] = build_word(p, 26);

    /* Word2: SOW[11:0] + delta_n high 10 bits */
    p = ((sow & 0xFFF) << 10) | ((uint32_t)f->delta_n >> 6 & 0x3FF);
    w[1] = build_word(p, 22);

    /* Word3: delta_n low 6 + Cuc high 16 */
    p = ((uint32_t)f->delta_n & 0x3F) << 16 |
        ((uint32_t)f->cuc >> 2 & 0xFFFF);
    w[2] = build_word(p, 22);

    /* Word4: Cuc low 2 + M0 high 20 */
    p = ((uint32_t)f->cuc & 0x3) << 20 |
        ((uint32_t)f->M0 >> 12 & 0xFFFFF);
    w[3] = build_word(p, 22);

    /* Word5: M0 low 12 + e high 10 */
    p = ((uint32_t)f->M0 & 0xFFF) << 10 |
        ((uint32_t)f->e >> 22 & 0x3FF);
    w[4] = build_word(p, 22);

    /* Word6: e low 22 */
    p = (uint32_t)f->e & 0x3FFFFF;
    w[5] = build_word(p, 22);

    /* Word7: Cus + Crc high 4 */
    p = ((uint32_t)f->cus & 0x3FFFF) << 4 |
        ((uint32_t)f->crc >> 14 & 0xF);
    w[6] = build_word(p, 22);

    /* Word8: Crc low 14 + Crs high 8 */
    p = ((uint32_t)f->crc & 0x3FFF) << 8 |
        ((uint32_t)f->crs >> 10 & 0xFF);
    w[7] = build_word(p, 22);

    /* Word9: Crs low 10 + sqrtA high 12 */
    p = ((uint32_t)f->crs & 0x3FF) << 12 |
        ((uint32_t)f->sqrtA >> 20 & 0xFFF);
    w[8] = build_word(p, 22);

    /* Word10: sqrtA low 20 + toe high 2 */
    p = ((uint32_t)f->sqrtA & 0xFFFFF) << 2 |
        ((f->toe >> 15) & 0x3);
    w[9] = build_word(p, 22);
}
/* --------------------------------- 子帧 3 ----------------------------------- */
static void build_subframe3_d1(const B1I_D1_Frame *f, uint32_t sow, uint32_t w[10])
{
    uint32_t p;
    /* Word1: preamble + reserved + FraID=3 + SOW[19:12] */
    p = 0x712;
    p = (p << 4);
    p = (p << 3) | 3;
    p = (p << 8) | ((sow >> 12) & 0xFF);
    w[0] = build_word(p, 26);

    /* Word2: SOW[11:0] + toe bits[14:5] */
    p = ((sow & 0xFFF) << 10) | ((f->toe >> 5) & 0x3FF);
    w[1] = build_word(p, 22);

    /* Word3: toe low5 + i0 high17 */
    p = ((f->toe & 0x1F) << 17) | ((uint32_t)f->i0 >> 15 & 0x1FFFF);
    w[2] = build_word(p, 22);

    /* Word4: i0 low15 + Cic high7 */
    p = ((uint32_t)f->i0 & 0x7FFF) << 7 | ((uint32_t)f->cic >> 11 & 0x7F);
    w[3] = build_word(p, 22);

    /* Word5: Cic low11 + omegadot high11 */
    p = ((uint32_t)f->cic & 0x7FF) << 11 | ((uint32_t)f->omegadot >> 13 & 0x7FF);
    w[4] = build_word(p, 22);

    /* Word6: omegadot low13 + Cis high9 */
    p = ((uint32_t)f->omegadot & 0x1FFF) << 9 | ((uint32_t)f->cis >> 9 & 0x1FF);
    w[5] = build_word(p, 22);

    /* Word7: Cis low9 + IDOT high13 */
    p = ((uint32_t)f->cis & 0x1FF) << 13 | ((uint32_t)f->idot >> 1 & 0x1FFF);
    w[6] = build_word(p, 22);

    /* Word8: IDOT low1 + omega0 high21 */
    p = ((uint32_t)f->idot & 0x1) << 21 |
        ((uint32_t)f->omega0 >> 11 & 0x1FFFFF);
    w[7] = build_word(p, 22);

    /* Word9: omega0 low11 + omega high11 */
    p = ((uint32_t)f->omega0 & 0x7FF) << 11 | ((uint32_t)f->omega >> 21 & 0x7FF);
    w[8] = build_word(p, 22);

    /* Word10: omega low21 + reserved bit */
    p = ((uint32_t)f->omega & 0x1FFFFF) << 1;
    w[9] = build_word(p, 22);
}

/* --------------------------------- 子帧 4 ----------------------------------- */
static void build_subframe4_d1(uint8_t *out, int week, double sow,
                               double frame_len)
{
    (void)week; /* not used */

    /* Subframe 4 carries almanac pages.  Actual parameters are not yet
     * implemented, so we keep the word layout but fill the data portion
     * with the official "10" placeholder pattern (0x2AAAAA). */

    memset(out, 0, SF_STREAM_LEN);

    uint32_t sow_int = (uint32_t)(floor(sow / frame_len) * frame_len);

    /* word1: preamble + reserved + FraID=4 + SOW[19:12] */
    uint32_t w = 0x712;
    w = (w << 4);
    w = (w << 3) | 4;
    w = (w << 8) | ((sow_int >> 12) & 0xFF);
    put_word(out, 0, build_word(w, 26));

    /* word2: SOW[11:0] + reserved + PNUM + spare (1,0) */
    uint32_t w2 = ((sow_int & 0xFFF) << 10) | (0 << 9) | (0 << 2) | 0x2;
    put_word(out, 30, build_word(w2, 22));

    uint32_t filler = 0x2AAAAA; /* 1010... repeated */
    for (int i = 2; i < 10; ++i)
        put_word(out, i * 30, build_word(filler, 22));
}
/* --------------------------------- 子帧 5 ----------------------------------- */
static void build_subframe5_d1(uint8_t *out, int week, double sow,
                               double frame_len)
{
    (void)week;

    /* Identical to subframe 4 but with FraID=5.  Data words use the same
     * placeholder pattern since almanac/UTC parameters are not generated. */

    memset(out, 0, SF_STREAM_LEN);

    uint32_t sow_int = (uint32_t)(floor(sow / frame_len) * frame_len);

    /* word1: preamble + reserved + FraID=5 + SOW[19:12] */
    uint32_t w = 0x712;
    w = (w << 4);
    w = (w << 3) | 5;
    w = (w << 8) | ((sow_int >> 12) & 0xFF);
    put_word(out, 0, build_word(w, 26));

    /* word2: SOW[11:0] + reserved + PNUM + spare (1,0) */
    uint32_t w2 = ((sow_int & 0xFFF) << 10) | (0 << 9) | (0 << 2) | 0x2;
    put_word(out, 30, build_word(w2, 22));

    uint32_t filler = 0x2AAAAA;
    for (int i = 2; i < 10; ++i)
        put_word(out, i * 30, build_word(filler, 22));
}
/* D2 subframe 1 with page number field */
static void build_subframe1_d2(uint8_t *out, const ephemeris_t *e,
                               int week, double sow)
{
    (void)out; (void)e; (void)week; (void)sow;
}
/* D2 subframe 2 with page number field */
static void build_subframe2_d2(uint8_t *out, const ephemeris_t *e, double sow)
{
    (void)out; (void)e; (void)sow;
}

/* ------------------ D1 Subframe Assembly Helpers ------------------------- */




static void assemble_subframe2_d1(const B1I_D1_Frame *frm, uint32_t sow,
                                  uint32_t words[10])
{
    if (frm) build_subframe2_d1(frm, sow, words); else memset(words, 0, sizeof(uint32_t) * 10);
}

static void assemble_subframe3_d1(const B1I_D1_Frame *frm, uint32_t sow,
                                  uint32_t words[10])
{
    if (frm) build_subframe3_d1(frm, sow, words); else memset(words, 0, sizeof(uint32_t) * 10);
}

/* ---------------------------------------------------------------------------- */
static uint8_t sf_static[MAX_SAT][4][SF_STREAM_LEN]; /* subframe 2-5 */

void navbits_init(void)
{
    for(int prn=1;prn<=63;prn++){
        B1I_D1_Frame frm;
        frame_from_ephemeris(&eph[prn], &frm);
        uint32_t words[10];

        assemble_subframe2_d1(&frm, 0, words);
        for(int w=0; w<10; ++w)
            put_word(sf_static[prn][0], w*30, words[w]);

        assemble_subframe3_d1(&frm, 0, words);
        for(int w=0; w<10; ++w)
            put_word(sf_static[prn][1], w*30, words[w]);

        build_subframe4_d1(sf_static[prn][2],0,0,6.0);
        build_subframe5_d1(sf_static[prn][3],0,0,6.0);
    }
}

/* 根據時間取得子帧 bit 流 (300 bits) */
void get_subframe_bits(int prn,int sf_id,int week,double sow,
                       double frame_len,uint8_t *out)
{
    double start = floor(sow/frame_len)*frame_len;
    int is_d2 = (frame_len < 6.0);
    if(sf_id==1){
        if(is_d2) build_subframe1_d2(out,&eph[prn],week,start);
        else      build_subframe1_d1(out,&eph[prn],week,start,frame_len);
    }else if(sf_id==2){
        if(is_d2) build_subframe2_d2(out,&eph[prn],start);
        else {
            B1I_D1_Frame frm;
            frame_from_ephemeris(&eph[prn], &frm);
            uint32_t words[10];
            assemble_subframe2_d1(&frm, (uint32_t)start, words);
            for(int w=0; w<10; ++w)
                put_word(out, w*30, words[w]);
        }
    }else if(sf_id==3){
        if(is_d2) memcpy(out,sf_static[prn][1],SF_STREAM_LEN);
        else {
            B1I_D1_Frame frm;
            frame_from_ephemeris(&eph[prn], &frm);
            uint32_t words[10];
            assemble_subframe3_d1(&frm, (uint32_t)start, words);
            for(int w=0; w<10; ++w)
                put_word(out, w*30, words[w]);
        }
    }else if(sf_id==4){
        build_subframe4_d1(out,week,start,frame_len);
    }else if(sf_id==5){
        build_subframe5_d1(out,week,start,frame_len);

    }else{
        memset(out,0,SF_STREAM_LEN);
    }
}

