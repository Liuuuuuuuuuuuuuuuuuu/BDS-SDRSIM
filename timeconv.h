/* timeconv.h ------------------------------------------------------------ */
/* 必須在 *任何* 系統標頭之前宣告： */
#ifndef TIMECONV_H
#define TIMECONV_H

#ifndef _GNU_SOURCE
#define _GNU_SOURCE          /* 觸發 strptime() 宣告 (glibc/musl 皆可) */
#endif

#include <time.h>
#include <stddef.h>

/* 若編譯器仍找不到原型，手動補一份 ── 可移植到非 glibc 平台            */
#if !defined(__USE_XOPEN2K) && !defined(__USE_XOPEN) && !defined(__USE_GNU)
char *strptime(const char *restrict, const char *restrict, struct tm *restrict);
#endif

/* yyyy/mm/dd,hh:mm:ss → BDT week & sow；成功傳回 0，失敗 -1 */
static inline int utc_to_bdt(const char *utc_str, int *week, double *sow)
{
    struct tm tm = {0};
    if (!strptime(utc_str, "%Y/%m/%d,%H:%M:%S", &tm))
        return -1;

    time_t t = timegm(&tm);          /* UTC → Unix epoch */
    const time_t bdt0 = 1136073600;  /* 2006-01-01 00:00:00 UTC */
    time_t diff = t - bdt0;

    if (diff < 0) diff += 604800;    /* 粗略容錯（負值時補 1 週） */

    *week = diff / 604800;
    *sow  = (double)(diff % 604800);
    return 0;
}
#endif /* TIMECONV_H */

