#include <stdio.h>
#include <stdbool.h>
#include <string.h>
#include "../bdssim.h"
#include "../globals.h"
#include "../orbits.h"
#include "../timeconv.h"


int main(void)
{

    FILE *fp = fopen("GBM0MGXRAP_20251760000_01D_05M_ORB.SP3", "r");
    if(!fp){ perror("SP3"); return 1; }

    char l[256];
    int bw = 0;             /* first epoch BDT week */
    double bsow = 0.0;      /* first epoch BDT second-of-week */
    bool got_epoch = false; /* set once first epoch parsed */

    char current_utc[32] = ""; /* epoch time string */

    while(fgets(l, sizeof l, fp))
    {
        if(l[0] == '*'){
            int Y,M,D,h,m; double s;
            sscanf(l+1, "%d %d %d %d %d %lf", &Y,&M,&D,&h,&m,&s);
            int is = (int)s;
            double frac = s - is;

            /* original SP3 epoch (GPS time) */
            char sp3_utc[32];
            snprintf(sp3_utc, sizeof sp3_utc, "%04d/%02d/%02d,%02d:%02d:%02d",
                     Y, M, D, h, m, is);
            strncpy(current_utc, sp3_utc, sizeof current_utc);

            /* convert (epoch - 18 s) to BDT */
            struct tm tm = {0};
            tm.tm_year = Y - 1900;
            tm.tm_mon  = M - 1;
            tm.tm_mday = D;
            tm.tm_hour = h;
            tm.tm_min  = m;
            tm.tm_sec  = is;

            time_t t = timegm(&tm) - 18;
            struct tm tm2;
            gmtime_r(&t, &tm2);

            char utc[32];
            snprintf(utc, sizeof utc, "%04d/%02d/%02d,%02d:%02d:%02d",
                     tm2.tm_year + 1900, tm2.tm_mon + 1, tm2.tm_mday,
                     tm2.tm_hour, tm2.tm_min, tm2.tm_sec);

            if (utc_to_bdt(utc, &bw, &bsow) != 0) {
                fprintf(stderr, "utc_to_bdt failed\n");
                return 1;
            }
            bsow += frac; /* fractional seconds */

            if(!got_epoch){
                double start_bdt = bw*604800.0 + bsow;
                sim_config_t cfg = {0};
                snprintf(cfg.rinex_file, sizeof(cfg.rinex_file),
                         "BRDM00DLR_S_20251760000_01D_MN.rnx");
                if(!init_simulator(&cfg, start_bdt)){
                    fprintf(stderr, "init failed\n");
                    return 1;
                }
                got_epoch = true;
                continue;
            }else{
                break; /* second epoch -> stop */
            }
        }
        if(!got_epoch)
            continue;            /* skip header before first epoch */

        if(l[0]=='P' && l[1]=='C'){
            int prn; double x,y,z,clk;
            if(sscanf(l+1, "C%2d %lf %lf %lf %lf", &prn,&x,&y,&z,&clk)!=5) continue;
            if(prn<1 || prn>63 || eph[prn].prn==0) continue;
            if(prn <= 5 || (prn >= 59 && prn <= 63))
                continue;        /* skip GEO satellites */

            double xyz[3];
            calc_sat_position_velocity(prn, bw, bsow, xyz, NULL, NULL);

            printf("time=%s PRN%02d calc=(%.3f %.3f %.3f) sp3=(%.3f %.3f %.3f)\n",
                   current_utc, prn, xyz[0], xyz[1], xyz[2],
                   x*1000.0, y*1000.0, z*1000.0);
        }
    }
    fclose(fp);

    return 0;
}
