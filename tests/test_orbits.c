#include <stdio.h>
#include <stdbool.h>
#include <string.h>
#include "../bdssim.h"
#include "../globals.h"
#include "../orbits.h"
#include "../timeconv.h"


int main(void)
{
    sim_config_t cfg = {0};
    snprintf(cfg.rinex_file, sizeof(cfg.rinex_file),
             "BRDM00DLR_S_20251760000_01D_MN.rnx");
    if(!init_simulator(&cfg, 0.0)){
        fprintf(stderr, "init failed\n");
        return 1;
    }

    FILE *fp = fopen("GBM0MGXRAP_20251760000_01D_05M_ORB.SP3", "r");
    if(!fp){ perror("SP3"); return 1; }

    char l[256];
    int bw = 0;             /* current BDT week */
    double bsow = 0.0;      /* current BDT second-of-week */
    bool got_epoch = false; /* true once the first epoch is parsed */

    int count = 0;          /* number of samples compared */
    double err2_sum = 0.0;  /* squared error accumulator  */
    bool printed[MAX_SAT] = {0};

    while(fgets(l, sizeof l, fp))
    {
        if(l[0] == '*'){
            int Y,M,D,h,m; double s;
            sscanf(l+1, "%d %d %d %d %d %lf", &Y,&M,&D,&h,&m,&s);
            char utc[32];
            int is = (int)s;
            snprintf(utc, sizeof utc, "%04d/%02d/%02d,%02d:%02d:%02d", Y,M,D,h,m,is);
            if (utc_to_bdt(utc, &bw, &bsow) != 0) {
                fprintf(stderr, "utc_to_bdt failed\n");
                return 1;
            }
            bsow += s - is; /* add fractional seconds */
            if(got_epoch) break; /* second epoch -> stop */
            got_epoch = true;
            continue;
        }
        if(!got_epoch)
            continue;            /* skip header before first epoch */

        if(l[0]=='P' && l[1]=='C'){
            int prn; double x,y,z,clk;
            if(sscanf(l+1, "C%2d %lf %lf %lf %lf", &prn,&x,&y,&z,&clk)!=5) continue;
            if(prn<1 || prn>63 || eph[prn].prn==0) continue;

            if(!printed[prn]){
                printf("PRN%02d sqrtA=%.3f omegadot=%g\n",
                       prn, eph[prn].sqrtA, eph[prn].omegadot);
                printed[prn] = true;
            }


            double xyz[3];
            calc_sat_position_velocity(prn, bw, bsow, xyz, NULL);

            double dx = xyz[0] - x*1000.0;
            double dy = xyz[1] - y*1000.0;
            double dz = xyz[2] - z*1000.0;
            double err = sqrt(dx*dx + dy*dy + dz*dz);
            char *nl = strchr(l, '\n');
            if(nl) *nl = '\0';
            printf("%s bw=%d bsow=%.3f PRN%02d calc=(%.3f %.3f %.3f) err=%.3f\n",
                   l, bw, bsow, prn, xyz[0], xyz[1], xyz[2], err);
            err2_sum += err*err;
            ++count;
        }
    }
    fclose(fp);

    double rms = (count>0)? sqrt(err2_sum/count) : 0.0;
    printf("Compared %d points, RMS error = %.3f m\n", count, rms);

    cleanup_simulator();
    return (rms < 10.0) ? 0 : 1; /* return error if RMS >10 m */
}
