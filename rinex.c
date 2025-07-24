/* rinex.c  –  解析 BDS RINEX3.04 NAV，僅提取 MEO/IGSO D1 需要的欄位
 *            (c) 2025  your-name
 * ------------------------------------------------------------------------ */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <time.h>
#include "bdssim.h"
#include "globals.h"

/* ---------- 小工具 ---------- */
static double fld(const char *ln, int idx, int ind)
{
    char buf[20];
    int pos = ind + idx * 19;
    if ((int)strlen(ln) <= pos) return 0.0;
    memcpy(buf, ln + pos, 19);
    buf[19] = '\0';
    for (int i = 0; i < 19; i++)
        if (buf[i] == 'D' || buf[i] == 'd') buf[i] = 'E';
    char *end;
    double v = strtod(buf, &end);
    return (end == buf) ? 0.0 : v;
}
static int ifld(const char *ln, int idx, int ind)
{
    return (int)fld(ln, idx, ind);
}

/* very-tiny, zone-less timegm (UTC→UNIX)  –  只考慮閏年規則，可跨平台 */

/* ---------- 主要解析 ---------- */
int read_rinex_nav(const char *fname)
{
    FILE *fp = fopen(fname, "r");
    if (!fp) { perror(fname);  return -1; }

    char l[120];

    /* skip header --------------------------------------------------- */
    while (fgets(l, sizeof(l), fp))
        if (strstr(l, "END OF HEADER")) break;

    const int IND1 = 23;   /* 第一行 offset */
    const int INDN = 4;    /* 後 7 行 offset */

    nav_time_min = 1e20;
    nav_time_max = -1e20;

    while (fgets(l, sizeof(l), fp))
    {
        if (l[0] != 'C' || !isdigit(l[1]))          /* 只處理北斗 */
            continue;

        int prn = atoi(l + 1);
        if (prn < 1 || prn > 63) {                  /* out of range：跳過 7 行 */
            char skip[120];
            for (int i = 0; i < 7; ++i)
                if (!fgets(skip, sizeof skip, fp)) break;
            continue;
        }

        ephemeris_t *e = &eph[prn];
        if(e->prn){
            char skip[120];
            for(int i=0;i<7;++i) if(!fgets(skip,sizeof skip,fp)) return -1;
            continue;                       /* keep first record only */
        }
        memset(e, 0, sizeof(*e));
        e->prn = prn;

        /* ── 1st row：af0/1/2 & UTC stamp -------------------------- */
        e->af0 = fld(l, 0, IND1);
        e->af1 = fld(l, 1, IND1);
        e->af2 = fld(l, 2, IND1);

        /* 先略過年月日等欄位，僅取後續週數與 TOE/SOW  */

        /* ── 讀後 7 行 -------------------------------------------- */
        char r[7][120];
        for (int i = 0; i < 7; ++i)
            if (!fgets(r[i], sizeof r[i], fp)) return -1;

        e->aode    = ifld(r[0], 0, INDN);
        e->crs     = fld(r[0], 1, INDN);
        e->deltan  = fld(r[0], 2, INDN);
        e->M0      = fld(r[0], 3, INDN);

        e->cuc     = fld(r[1], 0, INDN);
        e->e       = fld(r[1], 1, INDN);
        e->cus     = fld(r[1], 2, INDN);
        e->sqrtA   = fld(r[1], 3, INDN);

        e->toe     = fld(r[2], 0, INDN);
        e->toc     = (int)e->toe;          /* RINEX 無 toc 欄位，沿用 toe */
        e->cic     = fld(r[2], 1, INDN);
        e->omega0  = fld(r[2], 2, INDN);
        e->cis     = fld(r[2], 3, INDN);

        e->i0      = fld(r[3], 0, INDN);
        e->crc     = fld(r[3], 1, INDN);
        e->w       = fld(r[3], 2, INDN);
        e->omegadot= fld(r[3], 3, INDN);

        e->idot    = fld(r[4], 0, INDN);
        e->freq_num = ifld(r[4], 1, INDN); /* reserved/freq number */
        int week_r = ifld(r[4], 2, INDN);  /* BDS week number */
        e->week = week_r;

        /* line 7: SV accuracy/health and TGD1/2 */
        e->ura     = ifld(r[5], 0, INDN) & 0xF;
        e->health  = ifld(r[5], 1, INDN) & 0x1;
        e->tgd1    = fld(r[5], 2, INDN);
        e->tgd2    = fld(r[5], 3, INDN);

        /* line 8: transmission time of message and AODC */
        e->toe_msg = fld(r[6], 0, INDN);
        e->aodc    = ifld(r[6], 1, INDN);

        double t_bdt = e->week * 604800.0 + e->toe;
        if (t_bdt < nav_time_min) nav_time_min = t_bdt;
        if (t_bdt > nav_time_max) nav_time_max = t_bdt;
    }

    fclose(fp);
    puts("[rinex] 北斗星曆已載入");
    return 0;
}

