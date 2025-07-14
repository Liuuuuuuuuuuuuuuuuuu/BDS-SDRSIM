#include <stdio.h>
#include "../bdssim.h"
#include "../globals.h"
#include "../orbits.h"

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
    for(int prn=1; prn<=63; ++prn){
        if(eph[prn].prn==0) continue;
        double xyz[3];
        calc_sat_position_velocity(prn, nav_week, sow, xyz, NULL);
        printf("PRN%02d @ W%d %.1f -> X %.3f Y %.3f Z %.3f\n",
               prn, nav_week, sow, xyz[0], xyz[1], xyz[2]);
    }

    cleanup_simulator();
    return 0;
}
