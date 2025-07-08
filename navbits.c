/* navbits.c : 產生 B1I D1 子帧 1–3 (星曆) */

#include <math.h>
#include <string.h>
#include "navbits.h"
#include "bch.h"

extern ephemeris_t eph[MAX_SAT];

/* --------------------------------- 宏 & 工具 -------------------------------- */
static uint32_t make_word30(uint32_t info11)
{
    /* info11(高→低) 先 BCH(15,11,1) 產出 15 bits，再左右交织成 30 bits */
    uint16_t bch = bch_encode(info11);
    uint32_t w = 0;
    for(int i=0;i<15;i++){
        w <<= 2;
        w |= ((info11>>(14-i))&1);      /* 左：信息位 */
        w |= ((bch    >>i)&1)<<1;       /* 右：校验位 */
    }
    return w;
}

/* 填 30-bit word 至 bit 流。pos 為 bit index 0..299 (MSB first) */
static void put_word(uint8_t *b, int pos, uint32_t w30)
{
    for(int i=0;i<30;i++)
        b[pos+i] = (w30>>(29-i)) & 1;
}

/* --------------------------------- 子帧 1 ----------------------------------- */
static void build_subframe1(uint8_t *out, const ephemeris_t *e,
                             int week, double sow)
{
    memset(out,0,SF_STREAM_LEN);

    /* word1：帧同步 11 bits（11100010010） + FraID(001) + SOW[19:12] */
    uint16_t info = 0x712;            /* 11100010010 */
    info = ((info << 3)|0x1) & 0x7FF; /* +FraID=001 */
    uint32_t sow_int = (uint32_t)(floor(sow/6.0)*6.0);
    info = (info<<8) | ((sow_int>>12) & 0xFF);
    put_word(out,0, make_word30(info));

    /* word2: SOW[11:0] (12 bits) + WN (13 bits) + reserved(?) */
    uint32_t info2 = ((sow_int & 0xFFF)<<13) | (week&0x1FFF);
    put_word(out,30, make_word30(info2));

    /* word3: URA 4b, health 1b, toc(17b, /8), AODC 5b(=0) */
    uint32_t info3 = (e->ura & 0xF)<<26 | (e->health&1)<<25 | ((uint32_t)(e->toc/8)&0x1FFFF)<<8;
    put_word(out,60, make_word30(info3));

    /* word4: af0 (dt=2^-33) 22 bits → twos-comp */
    int32_t a0_i = (int32_t)llround(e->af0 / pow(2,-33));
    put_word(out,90, make_word30((uint32_t)(a0_i & 0x3FFFFF)));

    /* word5: af1 16 bits(+2^-50)、af2 11 bits(+2^-66) */
    int32_t a1_i = (int32_t)llround(e->af1 / pow(2,-50));
    int32_t a2_i = (int32_t)llround(e->af2 / pow(2,-66));
    uint32_t info5 = ((a1_i & 0xFFFF)<<11) | (a2_i & 0x7FF);
    put_word(out,120, make_word30(info5));

    /* 前 5 words 重複一次構成 10 words */
    for(int i=0;i<HALF_SUBFRAME_BITS;i++) out[i+HALF_SUBFRAME_BITS]=out[i];
}

/* --------------------------------- 子帧 2 ----------------------------------- */
static void build_subframe2(uint8_t *out, const ephemeris_t *e)
{
    memset(out,0,SF_STREAM_LEN);

    /* word1: FrameSync + FraID=010 + toe[15:8] */
    uint16_t info = 0x712;            /* sync */
    info = ((info<<3)|0x2) & 0x7FF;   /* FraID=010 */
    info = (info<<8) | ((uint32_t)e->toe>>8 &0xFF);
    put_word(out,0, make_word30(info));

    /* word2: toe[7:0] + √A[21:13] */
    uint32_t info2 = ((uint32_t)e->toe &0xFF)<<13 |
                     (uint32_t)(llround(e->sqrtA/pow(2,-19))>>13 &0x1FFF);
    put_word(out,30, make_word30(info2));

    /* word3: √A[12:0] + e[21:11] */
    uint32_t sqrtA_i = (uint32_t)llround(e->sqrtA/pow(2,-19));
    uint32_t e_i     = (uint32_t)llround(e->e     /pow(2,-33));
    uint32_t info3 = (sqrtA_i &0x1FFF)<<11 | (e_i>>11 &0x7FF);
    put_word(out,60, make_word30(info3));

    /* word4: e[10:0] + Δn(16b,+2^-43) + M0[21] */
    uint32_t dn_i = (int32_t)llround(e->deltan/pow(2,-43)) & 0xFFFF;
    uint32_t M0_i = (int32_t)llround(e->M0/pow(2,-31));
    uint32_t info4 = (e_i &0x7FF)<<19 | dn_i<<3 | (M0_i>>21 &0x7);
    put_word(out,90, make_word30(info4));

    /* word5: M0[20:0] 首段 */
    put_word(out,120, make_word30(M0_i & 0x1FFFFF));

    for(int i=0;i<HALF_SUBFRAME_BITS;i++) out[i+HALF_SUBFRAME_BITS]=out[i];
}

/* --------------------------------- 子帧 3 ----------------------------------- */
static void build_subframe3(uint8_t *out, const ephemeris_t *e)
{
    memset(out,0,SF_STREAM_LEN);

    /* word1 : sync + FraID=011 + Ω0[21:14] */
    uint16_t info = 0x712;
    info = ((info<<3)|0x3) & 0x7FF;
    int32_t O0_i = (int32_t)llround(e->omega0 / pow(2,-31));
    info = (info<<8) | ((O0_i>>14)&0xFF);
    put_word(out,0, make_word30(info));

    /* word2: Ω0[13:0] + i0[21:13] */
    int32_t i0_i = (int32_t)llround(e->i0 / pow(2,-31));
    uint32_t info2 = ((O0_i &0x3FFF)<<13) | ((i0_i>>13)&0x1FFF);
    put_word(out,30, make_word30(info2));

    /* word3: i0[12:0] + ω[21:11] */
    int32_t w_i = (int32_t)llround(e->w / pow(2,-31));
    uint32_t info3 = ((i0_i &0x1FFF)<<11) | ((w_i>>11)&0x7FF);
    put_word(out,60, make_word30(info3));

    /* word4: ω[10:0] + CRC + IDOT[10] */
    int32_t crc_i = (int32_t)llround(e->crc / pow(2,-5));       /* 18 bits */
    int32_t idot_i= (int32_t)llround(e->idot/ pow(2,-43));      /* 14 bits → 11+3 */
    uint32_t info4 = ((w_i &0x7FF)<<19) | ((crc_i&0x3FFFF)<<1) | ((idot_i>>10)&1);
    put_word(out,90, make_word30(info4));

    /* word5: IDOT[9:0] + Ω̇ dot[23:?] */
    int32_t odot_i = (int32_t)llround(e->omegadot / pow(2,-43));/* 24 bits */
    uint32_t info5 = ((idot_i &0x3FF)<<20) | (odot_i &0xFFFFF);
    put_word(out,120, make_word30(info5));

    for(int i=0;i<HALF_SUBFRAME_BITS;i++) out[i+HALF_SUBFRAME_BITS]=out[i];
}

/* --------------------------------- 子帧 4 ----------------------------------- */

static void build_sf45_template(uint8_t *out)
{
    /*
     * 每個 word 填入 30-bit 固定值 0x2AAAAAAA (= 1010... pattern)。
     * 依官方文件，此為子幀 4/5 的預留欄位內容。
     */
    for(int w=2; w<=10; ++w)
        put_word(out,(w-1)*30, 0x2AAAAAAA);
}

static void build_subframe4(uint8_t *out,int week,double sow)

{
    memset(out,0,SF_STREAM_LEN);
    uint16_t info = 0x712;            /* sync */
    info = ((info<<3)|0x4) & 0x7FF;   /* FraID=100 */

    uint32_t sow_int = (uint32_t)(floor(sow/6.0)*6.0);
    info = (info<<8) | ((sow_int>>12) & 0xFF);
    put_word(out,0, make_word30(info));
    uint32_t info2 = ((sow_int&0xFFF)<<13) | (week&0x1FFF);
    put_word(out,30, make_word30(info2));
    build_sf45_template(out);
    for(int i=0;i<HALF_SUBFRAME_BITS;i++) out[i+HALF_SUBFRAME_BITS]=out[i];
}

/* --------------------------------- 子帧 5 ----------------------------------- */

static void build_subframe5(uint8_t *out,int week,double sow)

{
    memset(out,0,SF_STREAM_LEN);
    uint16_t info = 0x712;            /* sync */
    info = ((info<<3)|0x5) & 0x7FF;   /* FraID=101 */

    uint32_t sow_int = (uint32_t)(floor(sow/6.0)*6.0);
    info = (info<<8) | ((sow_int>>12) & 0xFF);
    put_word(out,0, make_word30(info));
    uint32_t info2 = ((sow_int & 0xFFF)<<13) | (week&0x1FFF);
    put_word(out,30, make_word30(info2));
    build_sf45_template(out);
    for(int i=0;i<HALF_SUBFRAME_BITS;i++) out[i+HALF_SUBFRAME_BITS]=out[i];
}

/* ---------------------------------------------------------------------------- */
static uint8_t sf_static[MAX_SAT][4][SF_STREAM_LEN]; /* subframe 2-5 */

void navbits_init(void)
{
    for(int prn=1;prn<=63;prn++){
        build_subframe2(sf_static[prn][0], &eph[prn]);
        build_subframe3(sf_static[prn][1], &eph[prn]);
        build_subframe4(sf_static[prn][2],0,0);
        build_subframe5(sf_static[prn][3],0,0);
    }
}

/* 根據時間取得子帧 bit 流 (300 bits) */
void get_subframe_bits(int prn,int sf_id,int week,double sow,uint8_t *out)
{
    if(sf_id==1){
        build_subframe1(out,&eph[prn],week,sow);
    }else if(sf_id==2){
        memcpy(out,sf_static[prn][0],SF_STREAM_LEN);
    }else if(sf_id==3){
        memcpy(out,sf_static[prn][1],SF_STREAM_LEN);
    }else if(sf_id==4){

        build_subframe4(out,week,sow);
    }else if(sf_id==5){
        build_subframe5(out,week,sow);

    }else{
        memset(out,0,SF_STREAM_LEN);
    }
}

