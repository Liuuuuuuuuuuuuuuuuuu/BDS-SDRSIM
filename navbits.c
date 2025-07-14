/* navbits.c : 產生 B1I D1/D2 子帧 (星曆) */

#include <math.h>
#include <string.h>
#include "navbits.h"
#include "bch.h"
#include "nav_words.h"

extern ephemeris_t eph[MAX_SAT];
#include "nav_pages.h"
#include "icd_fields.h"

/* Convert ephemeris_t to the simplified B1I_D1_Frame used by
 * nav_pages.c. Only the fields actually consumed by those
 * helpers are populated. */
static void eph_to_frame(const ephemeris_t *e, B1I_D1_Frame *f)
{
    memset(f, 0, sizeof(*f));
    f->toc      = (uint32_t)(e->toc / 8);
    f->aodc     = e->iode & 0x1F;
    f->urai     = e->ura & 0xF;
    f->satH1    = e->health & 0x1;

    f->a0       = (int32_t)llround(e->af0 / pow(2, -33));
    f->a1       = (int32_t)llround(e->af1 / pow(2, -50));
    f->a2       = (int32_t)llround(e->af2 / pow(2, -66));

    f->toe      = (uint32_t)(e->toe / 8);
    f->sqrtA    = (int32_t)llround(e->sqrtA / pow(2, -19));
    f->e        = (int32_t)llround(e->e / pow(2, -33));
    f->delta_n  = (int32_t)llround(e->deltan / pow(2, -43));
    f->M0       = (int32_t)llround(e->M0 / pow(2, -31));

    f->omega0   = (int32_t)llround(e->omega0 / pow(2, -31));
    f->i0       = (int32_t)llround(e->i0 / pow(2, -31));
    f->omega    = (int32_t)llround(e->w / pow(2, -31));
    f->crc      = (int32_t)llround(e->crc / pow(2, -5));
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
static void build_subframe1(uint8_t *out, const ephemeris_t *e,
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

    /* word2: SOW[11:0] (12 bits) + WN (13 bits) + reserved(?) */
    uint32_t info2 = ((sow_int & 0xFFF)<<13) | (week&0x1FFF);
    put_word(out,30, build_word(info2, 26));

    /* word3: URA 4b, health 1b, toc(17b, /8), AODC 5b(=0) */
    uint32_t info3 = (e->ura & 0xF)<<26 | (e->health&1)<<25 | ((uint32_t)(e->toc/8)&0x1FFFF)<<8;
    put_word(out,60, build_word(info3, 26));

    /* word4: af0 (dt=2^-33) 22 bits → twos-comp */
    int32_t a0_i = (int32_t)llround(e->af0 / pow(2,-33));
    put_word(out,90, build_word((uint32_t)(a0_i & 0x3FFFFF), 26));

    /* word5: af1 16 bits(+2^-50)、af2 11 bits(+2^-66) */
    int32_t a1_i = (int32_t)llround(e->af1 / pow(2,-50));
    int32_t a2_i = (int32_t)llround(e->af2 / pow(2,-66));
    uint32_t info5 = ((a1_i & 0xFFFF)<<11) | (a2_i & 0x7FF);
    put_word(out,120, build_word(info5, 22));

    /* 前 5 words 重複一次構成 10 words */
    for(int i=0;i<HALF_SUBFRAME_BITS;i++) out[i+HALF_SUBFRAME_BITS]=out[i];
}

/* D2 subframe 1 with page number field */
static void build_subframe1_d2(uint8_t *out,const ephemeris_t *e,
                               int week,double sow)
{
    memset(out,0,SF_STREAM_LEN);
    uint32_t sow_int = (uint32_t)(floor(sow/0.6)*0.6);
    int frame = (int)fmod(sow_int,360.0)/3; /* 0..119 */
    uint32_t pnum = (frame % 10) + 1;       /* 1..10 */

    uint32_t info = 0x712;                  /* preamble */
    info = (info<<4) | (pnum & 0xF);        /* PNUM1 */
    info = (info<<3) | 0x1;                 /* FraID=001 */
    info = (info<<8) | ((sow_int>>12)&0xFF);
    put_word(out,0, build_word(info,26));

    uint32_t info2 = ((sow_int&0xFFF)<<13) | (week&0x1FFF);
    put_word(out,30, build_word(info2,26));

    uint32_t info3 = (e->ura&0xF)<<26 | (e->health&1)<<25 | ((uint32_t)(e->toc/8)&0x1FFFF)<<8;
    put_word(out,60, build_word(info3,26));

    int32_t a0_i = (int32_t)llround(e->af0/pow(2,-33));
    put_word(out,90, build_word((uint32_t)(a0_i & 0x3FFFFF),26));

    int32_t a1_i = (int32_t)llround(e->af1/pow(2,-50));
    int32_t a2_i = (int32_t)llround(e->af2/pow(2,-66));
    uint32_t info5 = ((a1_i & 0xFFFF)<<11) | (a2_i & 0x7FF);
    put_word(out,120, build_word(info5,22));

    for(int i=0;i<HALF_SUBFRAME_BITS;i++) out[i+HALF_SUBFRAME_BITS]=out[i];
}

/* --------------------------------- 子帧 2 ----------------------------------- */
/*
 * build_subframe2 and build_subframe3 are kept for reference. They are not
 * invoked in the current code path but may be useful for future extensions.
 * Mark them as unused to silence compiler warnings.
 */
__attribute__((unused)) static void build_subframe2(uint8_t *out, const ephemeris_t *e)
{
    memset(out,0,SF_STREAM_LEN);

    /* word1: FrameSync + FraID=010 + toe[15:8] */
    uint32_t info = 0x712;            /* preamble */
    info = (info << 4);               /* reserved */
    info = (info << 3) | 0x2;         /* FraID=010 */
    info = (info << 8) | ((uint32_t)e->toe>>8 &0xFF);
    put_word(out,0, build_word(info, 26));

    /* word2: toe[7:0] + √A[21:13] */
    uint32_t info2 = ((uint32_t)e->toe &0xFF)<<13 |
                     (uint32_t)(llround(e->sqrtA/pow(2,-19))>>13 &0x1FFF);
    put_word(out,30, build_word(info2, 26));

    /* word3: √A[12:0] + e[21:11] */
    uint32_t sqrtA_i = (uint32_t)llround(e->sqrtA/pow(2,-19));
    uint32_t e_i     = (uint32_t)llround(e->e     /pow(2,-33));
    uint32_t info3 = (sqrtA_i &0x1FFF)<<11 | (e_i>>11 &0x7FF);
    put_word(out,60, build_word(info3, 26));

    /* word4: e[10:0] + Δn(16b,+2^-43) + M0[21] */
    uint32_t dn_i = (int32_t)llround(e->deltan/pow(2,-43)) & 0xFFFF;
    uint32_t M0_i = (int32_t)llround(e->M0/pow(2,-31));
    uint32_t info4 = (e_i &0x7FF)<<19 | dn_i<<3 | (M0_i>>21 &0x7);
    put_word(out,90, build_word(info4, 26));

    /* word5: M0[20:0] 首段 */
    put_word(out,120, build_word(M0_i & 0x1FFFFF, 22));

    for(int i=0;i<HALF_SUBFRAME_BITS;i++) out[i+HALF_SUBFRAME_BITS]=out[i];
}

/* D2 subframe 2 with page number field */
static void build_subframe2_d2(uint8_t *out,const ephemeris_t *e,double sow)
{
    memset(out,0,SF_STREAM_LEN);
    uint32_t sow_int = (uint32_t)(floor(sow/0.6)*0.6);
    int frame = (int)fmod(sow_int,360.0)/3; /* 0..119 */
    uint32_t pnum = (frame % 6) + 1;        /* 1..6 */

    uint32_t info = 0x712;
    info = (info<<4) | (pnum & 0xF);        /* PNUM2 */
    info = (info<<3) | 0x2;                 /* FraID=010 */
    info = (info<<8) | ((uint32_t)e->toe>>8 & 0xFF);
    put_word(out,0, build_word(info,26));

    uint32_t info2 = ((uint32_t)e->toe & 0xFF)<<13 |
                     (uint32_t)(llround(e->sqrtA/pow(2,-19))>>13 &0x1FFF);
    put_word(out,30, build_word(info2,26));

    uint32_t sqrtA_i = (uint32_t)llround(e->sqrtA/pow(2,-19));
    uint32_t e_i     = (uint32_t)llround(e->e/pow(2,-33));
    uint32_t info3 = (sqrtA_i &0x1FFF)<<11 | (e_i>>11 &0x7FF);
    put_word(out,60, build_word(info3,26));

    uint32_t dn_i = (int32_t)llround(e->deltan/pow(2,-43)) & 0xFFFF;
    uint32_t M0_i = (int32_t)llround(e->M0/pow(2,-31));
    uint32_t info4 = (e_i &0x7FF)<<19 | dn_i<<3 | (M0_i>>21 &0x7);
    put_word(out,90, build_word(info4,26));

    put_word(out,120, build_word(M0_i & 0x1FFFFF,22));

    for(int i=0;i<HALF_SUBFRAME_BITS;i++) out[i+HALF_SUBFRAME_BITS]=out[i];
}

/* --------------------------------- 子帧 3 ----------------------------------- */
__attribute__((unused)) static void build_subframe3(uint8_t *out, const ephemeris_t *e)
{
    memset(out,0,SF_STREAM_LEN);

    /* word1 : sync + FraID=011 + Ω0[21:14] */
    uint32_t info = 0x712;            /* preamble */
    info = (info << 4);
    info = (info << 3) | 0x3;         /* FraID=011 */
    int32_t O0_i = (int32_t)llround(e->omega0 / pow(2,-31));
    info = (info<<8) | ((O0_i>>14)&0xFF);
    put_word(out,0, build_word(info, 26));

    /* word2: Ω0[13:0] + i0[21:13] */
    int32_t i0_i = (int32_t)llround(e->i0 / pow(2,-31));
    uint32_t info2 = ((O0_i &0x3FFF)<<13) | ((i0_i>>13)&0x1FFF);
    put_word(out,30, build_word(info2, 26));

    /* word3: i0[12:0] + ω[21:11] */
    int32_t w_i = (int32_t)llround(e->w / pow(2,-31));
    uint32_t info3 = ((i0_i &0x1FFF)<<11) | ((w_i>>11)&0x7FF);
    put_word(out,60, build_word(info3, 26));

    /* word4: ω[10:0] + CRC + IDOT[10] */
    int32_t crc_i = (int32_t)llround(e->crc / pow(2,-5));       /* 18 bits */
    int32_t idot_i= (int32_t)llround(e->idot/ pow(2,-43));      /* 14 bits → 11+3 */
    uint32_t info4 = ((w_i &0x7FF)<<19) | ((crc_i&0x3FFFF)<<1) | ((idot_i>>10)&1);
    put_word(out,90, build_word(info4, 26));

    /* word5: IDOT[9:0] + Ω̇ dot[23:?] */
    int32_t odot_i = (int32_t)llround(e->omegadot / pow(2,-43));/* 24 bits */
    uint32_t info5 = ((idot_i &0x3FF)<<20) | (odot_i &0xFFFFF);
    put_word(out,120, build_word(info5, 22));

    for(int i=0;i<HALF_SUBFRAME_BITS;i++) out[i+HALF_SUBFRAME_BITS]=out[i];
}

/* --------------------------------- 子帧 4 ----------------------------------- */

static void build_sf45_template(uint8_t *out)
{
    /*
     * 每個 word 填入官方指定的固定值 0xAAAAAAAA (= 1010... pattern)。
     * 此為子幀 4/5 的預留欄位內容。
     */
    for(int w=2; w<=10; ++w)
        put_word(out,(w-1)*30, 0xAAAAAAAA);
}

static void build_subframe4(uint8_t *out,int week,double sow,double frame_len)

{
    memset(out,0,SF_STREAM_LEN);
    uint32_t info = 0x712;            /* preamble */
    info = (info << 4);
    info = (info << 3) | 0x4;         /* FraID=100 */

    uint32_t sow_int = (uint32_t)(floor(sow/frame_len)*frame_len);
    info = (info<<8) | ((sow_int>>12) & 0xFF);
    put_word(out,0, build_word(info, 26));
    uint32_t info2 = ((sow_int&0xFFF)<<13) | (week&0x1FFF);
    put_word(out,30, build_word(info2, 26));
    build_sf45_template(out);
    for(int i=0;i<HALF_SUBFRAME_BITS;i++) out[i+HALF_SUBFRAME_BITS]=out[i];
}

/* --------------------------------- 子帧 5 ----------------------------------- */

static void build_subframe5(uint8_t *out,int week,double sow,double frame_len)

{
    memset(out,0,SF_STREAM_LEN);
    uint32_t info = 0x712;            /* preamble */
    info = (info << 4);
    info = (info << 3) | 0x5;         /* FraID=101 */

    uint32_t sow_int = (uint32_t)(floor(sow/frame_len)*frame_len);
    info = (info<<8) | ((sow_int>>12) & 0xFF);
    put_word(out,0, build_word(info, 26));
    uint32_t info2 = ((sow_int & 0xFFF)<<13) | (week&0x1FFF);
    put_word(out,30, build_word(info2, 26));
    build_sf45_template(out);
    for(int i=0;i<HALF_SUBFRAME_BITS;i++) out[i+HALF_SUBFRAME_BITS]=out[i];
}

/* ---------------------------------------------------------------------------- */
static uint8_t sf_static[MAX_SAT][4][SF_STREAM_LEN]; /* subframe 2-5 */

void navbits_init(void)
{
    for(int prn=1;prn<=63;prn++){
        B1I_D1_Frame frm;
        eph_to_frame(&eph[prn], &frm);
        uint32_t words[10];

        assemble_d1_subframe2(&frm, words);
        for(int w=0; w<10; ++w)
            put_word(sf_static[prn][0], w*30, words[w]);

        assemble_d1_subframe3(&frm, words);
        for(int w=0; w<10; ++w)
            put_word(sf_static[prn][1], w*30, words[w]);

        build_subframe4(sf_static[prn][2],0,0,6.0);
        build_subframe5(sf_static[prn][3],0,0,6.0);
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
        else      build_subframe1(out,&eph[prn],week,start,frame_len);
    }else if(sf_id==2){
        if(is_d2) build_subframe2_d2(out,&eph[prn],start);
        else memcpy(out,sf_static[prn][0],SF_STREAM_LEN);
    }else if(sf_id==3){
        if(is_d2) memcpy(out,sf_static[prn][1],SF_STREAM_LEN);
        else memcpy(out,sf_static[prn][1],SF_STREAM_LEN);
    }else if(sf_id==4){
        build_subframe4(out,week,start,frame_len);
    }else if(sf_id==5){
        build_subframe5(out,week,start,frame_len);

    }else{
        memset(out,0,SF_STREAM_LEN);
    }
}

