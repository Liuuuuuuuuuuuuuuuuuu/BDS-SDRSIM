#ifndef TIMECONV_H
#define TIMECONV_H

#ifndef _GNU_SOURCE
#define _GNU_SOURCE
#endif

#include <time.h>
#include <stddef.h>

#if !defined(__USE_XOPEN2K) && !defined(__USE_XOPEN) && !defined(__USE_GNU)
char *strptime(const char *restrict, const char *restrict, struct tm *restrict);
#endif

static inline int utc_to_bdt(const char *utc_str, int *week, double *sow)
{
    struct tm tm = {0};
    if (!strptime(utc_str, "%Y/%m/%d,%H:%M:%S", &tm))
        return -1;

    time_t t = timegm(&tm);
    const time_t bdt0 = 1136073600;
    extern int utc_bdt_diff;
    const time_t week_sec = 604800;
    time_t diff = t - bdt0 + utc_bdt_diff;

    long w = diff / week_sec;
    long r = diff % week_sec;
    if (r < 0) { r += week_sec; --w; }

    *week = (int)w;
    *sow  = (double)r;
    return 0;
}

#endif /* TIMECONV_H */

