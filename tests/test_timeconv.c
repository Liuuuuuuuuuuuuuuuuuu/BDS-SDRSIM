#include <stdio.h>
#include <time.h>
#include "timeconv.h"

int utc_bdt_diff = 4;

int main(int argc, char *argv[])
{
    if (argc != 2) {
        fprintf(stderr, "Usage: %s YYYY/MM/DD,HH:MM:SS\n", argv[0]);
        return 1;
    }

    int week; double sow;
    if (utc_to_bdt(argv[1], &week, &sow) != 0) {
        fprintf(stderr, "Invalid UTC format\n");
        return 1;
    }

    struct tm tm = {0};
    if (!strptime(argv[1], "%Y/%m/%d,%H:%M:%S", &tm)) {
        fprintf(stderr, "Invalid UTC format\n");
        return 1;
    }
    extern int utc_bdt_diff;
    time_t t = timegm(&tm) + utc_bdt_diff; /* UTC -> BDT */

    struct tm bdt;
    gmtime_r(&t, &bdt);

    printf("BDT calendar: %04d/%02d/%02d,%02d:%02d:%02d\n",
           bdt.tm_year + 1900, bdt.tm_mon + 1, bdt.tm_mday,
           bdt.tm_hour, bdt.tm_min, bdt.tm_sec);
    printf("BDT week/sow: %d %.0f\n", week, sow);
    return 0;
}
