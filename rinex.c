/* rinex.c  –  解析 BDS RINEX3.05 NAV，僅提取 MEO/IGSO D1 需要的欄位
 *            (c) 2025  your-name
 * ------------------------------------------------------------------------ */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <time.h>
#include "bdssim.h"

extern ephemeris_t eph[MAX_SAT];

/* ---------- 小工具 ---------- */
static double fld(const char *ln, int idx, int ind)
{
    char buf[20];
    memcpy(buf, &ln[ind + idx * 19], 19);
    buf[19] = '\0';
    for (int i = 0; i < 19; i++)
        if (buf[i] == 'D' || buf[i] == 'd') buf[i] = 'E';
    return strtod(buf, NULL);
}
static int ifld(const char *ln, int idx, int ind)
{
    return (int)fld(ln, idx, ind);
}

/* very-tiny, zone-less timegm (UTC→UNIX)  –  只考慮閏年規則，可跨平台 */
static time_t my_timegm(struct tm *t)
{
    static const int days_in_mon[12] =
        {31,28,31,30,31,30,31,31,30,31,30,31};

    int y = t->tm_year + 1900;
    int m = t->tm_mon;
    int d = t->tm_mday - 1;             /* 0-based */

    /*  days since 1970-01-01 */
    int days = 0;
    for (int yr = 1970; yr < y; ++yr)
        days += 365 + ((yr % 4 == 0 && yr % 100) || (yr % 400 == 0));
    for (int mo = 0; mo < m; ++mo) {
        days += days_in_mon[mo];
        if (mo == 1 && ((y % 4 == 0 && y % 100) || (y % 400 == 0)))
            ++days;                     /* Feb 29 */
    }
    days += d;

    return (time_t)days * 86400
         + t->tm_hour * 3600 + t->tm_min * 60 + t->tm_sec;
}

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
        memset(e, 0, sizeof(*e));
        e->prn = prn;

        /* ── 1st row：af0/1/2 & UTC stamp -------------------------- */
        e->af0 = fld(l, 0, IND1);
        e->af1 = fld(l, 1, IND1);
        e->af2 = fld(l, 2, IND1);

        int Y, M, D, h, m;
        double s;
        sscanf(l + 4, "%4d %2d %2d %2d %2d %lf", &Y, &M, &D, &h, &m, &s);

        struct tm utc = {0};
        utc.tm_year = Y - 1900;  utc.tm_mon  = M - 1;  utc.tm_mday = D;
        utc.tm_hour = h;         utc.tm_min  = m;      utc.tm_sec  = (int)s;

        time_t unix_t = my_timegm(&utc);

        /* BDT: TAI-32.0 (leap) – 既然 RINEX 是 UTC，需加上 UTC-BDT 偏移 */
        const time_t BDT_EPOCH = 1136073600;    /* 2006-01-01 00:00:00 UTC */
        const int    LS_GPS_2025 = 18;          /* leap seconds 到 2025-06 */
        const int    LS_BDT      = 0;           /* BDT 無閏秒 */
        int bdt_off = (LS_GPS_2025 - LS_BDT);   /* 一律 18 秒 (現行) */

        time_t bdt_time = unix_t + bdt_off - BDT_EPOCH;
        e->toc = (int)(bdt_time % 604800);
        if (e->toc < 0) e->toc += 604800;
        e->week = (int)(bdt_time / 604800);

        /* ── 讀後 7 行 -------------------------------------------- */
        char r[7][120];
        for (int i = 0; i < 7; ++i)
            if (!fgets(r[i], sizeof r[i], fp)) return -1;

        e->crs     = fld(r[0], 1, INDN);
        e->deltan  = fld(r[0], 2, INDN);
        e->M0      = fld(r[0], 3, INDN);

        e->cuc     = fld(r[1], 0, INDN);
        e->e       = fld(r[1], 1, INDN);
        e->cus     = fld(r[1], 2, INDN);
        e->sqrtA   = fld(r[1], 3, INDN);

        e->toe     = fld(r[2], 0, INDN);
        e->cic     = fld(r[2], 1, INDN);
        e->omega0  = fld(r[2], 2, INDN);
        e->cis     = fld(r[2], 3, INDN);

        e->i0      = fld(r[3], 0, INDN);
        e->crc     = fld(r[3], 1, INDN);
        e->w       = fld(r[3], 2, INDN);
        e->omegadot= fld(r[3], 3, INDN);

        e->idot    = fld(r[4], 0, INDN);

        /* row 5-field2 = URAI / field3 = SatH1 (ICD) */
        e->ura     = ifld(r[4], 1, INDN) & 0xF;
        e->health  = ifld(r[4], 2, INDN) & 0x1;

        /* 星期計數也寫一次（部分接收機只讀 .week） */
        int week_rnx = ifld(r[4], 2, INDN);
        if (week_rnx > 0) e->week = week_rnx;
    }

    fclose(fp);
    puts("[rinex] 北斗星曆已載入");
    return 0;
}

