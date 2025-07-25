/* rinex.c  –  解析 BDS RINEX3.04 NAV，僅提取 MEO/IGSO D1 需要的欄位
 *            (c) 2025  your-name
 * ------------------------------------------------------------------------ */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <time.h>
#include <math.h>
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
int read_rinex_nav(const char *fname, double start_bdt)
{
    FILE *fp = fopen(fname, "r");
    if (!fp) { perror(fname);  return -1; }

    char l[120];

    /* parse header -------------------------------------------------- */
    while (fgets(l, sizeof(l), fp))
    {
        if(strncmp(l,"BDSA",4)==0){
            sscanf(l+4, "%lf %lf %lf %lf", &iono_alpha[0], &iono_alpha[1],
                   &iono_alpha[2], &iono_alpha[3]);
            continue;
        }
        if(strncmp(l,"BDSB",4)==0){
            sscanf(l+4, "%lf %lf %lf %lf", &iono_beta[0], &iono_beta[1],
                   &iono_beta[2], &iono_beta[3]);
            continue;
        }
        if(strncmp(l,"BDUT",4)==0){
            double a0,a1; int off;
            if(sscanf(l+4, "%lf %lf %d", &a0, &a1, &off)==3)
                utc_bdt_diff = off;
            continue;
        }
        if (strstr(l, "END OF HEADER")) break;
    }

    const int IND1 = 23;   /* 第一行 offset */
    const int INDN = 4;    /* 後 7 行 offset */

    nav_time_min = 1e20;
    nav_time_max = -1e20;
    double best_diff[MAX_SAT];
    for(int i=0;i<MAX_SAT;i++) best_diff[i]=1e20;

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

        ephemeris_t tmp = {0};
        tmp.prn = prn;

        /* ── 1st row：af0/1/2 & UTC stamp -------------------------- */
        tmp.af0 = fld(l, 0, IND1);
        tmp.af1 = fld(l, 1, IND1);
        tmp.af2 = fld(l, 2, IND1);

        /* 先略過年月日等欄位，僅取後續週數與 TOE/SOW  */

        /* ── 讀後 7 行 -------------------------------------------- */
        char r[7][120];
        for (int i = 0; i < 7; ++i)
            if (!fgets(r[i], sizeof r[i], fp)) return -1;

        tmp.aode    = ifld(r[0], 0, INDN);
        tmp.crs     = fld(r[0], 1, INDN);
        tmp.deltan  = fld(r[0], 2, INDN);
        tmp.M0      = fld(r[0], 3, INDN);

        tmp.cuc     = fld(r[1], 0, INDN);
        tmp.e       = fld(r[1], 1, INDN);
        tmp.cus     = fld(r[1], 2, INDN);
        tmp.sqrtA   = fld(r[1], 3, INDN);

        tmp.toe     = fld(r[2], 0, INDN);
        tmp.toc     = (int)tmp.toe;          /* RINEX 無 toc 欄位，沿用 toe */
        tmp.cic     = fld(r[2], 1, INDN);
        tmp.omega0  = fld(r[2], 2, INDN);
        tmp.cis     = fld(r[2], 3, INDN);

        tmp.i0      = fld(r[3], 0, INDN);
        tmp.crc     = fld(r[3], 1, INDN);
        tmp.w       = fld(r[3], 2, INDN);
        tmp.omegadot= fld(r[3], 3, INDN);

        tmp.idot    = fld(r[4], 0, INDN);
        int week_r = ifld(r[4], 2, INDN);  /* BDS week number */
        tmp.week = week_r;

        /* line 7: SV accuracy/health and TGD1/2 */
        tmp.ura     = ifld(r[5], 0, INDN) & 0xF;
        tmp.health  = ifld(r[5], 1, INDN) & 0x1;
        tmp.tgd1    = fld(r[5], 2, INDN);
        tmp.tgd2    = fld(r[5], 3, INDN);

        /* line 8: transmission time of message and AODC */
        tmp.toe_msg = fld(r[6], 0, INDN);
        tmp.aodc    = ifld(r[6], 1, INDN);

        double t_bdt = tmp.week * 604800.0 + tmp.toe;
        if (t_bdt < nav_time_min) nav_time_min = t_bdt;
        if (t_bdt > nav_time_max) nav_time_max = t_bdt;

        double diff = fabs(t_bdt - start_bdt);
        if(diff < best_diff[prn]){
            best_diff[prn] = diff;
            eph[prn] = tmp;
        }
    }

    fclose(fp);
    puts("[rinex] 北斗星曆已載入");
    return 0;
}

