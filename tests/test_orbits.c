#include <stdio.h>
#include "../bdssim.h"
#include "../globals.h"
#include "../orbits.h"
#include <time.h>

static void bdt_to_utc(int week, double sow, char *buf, size_t len)
{
    const time_t BDT_EPOCH = 1136073600; /* 2006-01-01 00:00:00 UTC */
    double sec = week * 604800.0 + sow;  /* BDT seconds since epoch */
    time_t t = (time_t)(BDT_EPOCH - 14 + sec); /* UTC epoch time */
    struct tm tm;
    gmtime_r(&t, &tm);
    strftime(buf, len, "%Y/%m/%d,%H:%M:%S", &tm);
}

int main(void)
{
    sim_config_t cfg = {0};
    snprintf(cfg.rinex_file, sizeof(cfg.rinex_file),
             "BRDM00DLR_S_20251760000_01D_MN.rnx");
    if(!init_simulator(&cfg)){
        fprintf(stderr, "init failed\n");
        return 1;
    }

    double sow = nav_time_min - nav_week*604800.0;
    char utc[32];
    bdt_to_utc(nav_week, sow, utc, sizeof utc);
    for(int prn=1; prn<=63; ++prn){
        if(eph[prn].prn==0) continue;
        double xyz[3];
        calc_sat_position_velocity(prn, nav_week, sow, xyz, NULL);
        printf("PRN%02d @ %s -> X %.3f Y %.3f Z %.3f\n",
               prn, utc, xyz[0], xyz[1], xyz[2]);
    }

    cleanup_simulator();
    return 0;
}
