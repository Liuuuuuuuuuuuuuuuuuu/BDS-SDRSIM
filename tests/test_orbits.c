#include <stdio.h>
#include <stdbool.h>
#include "../bdssim.h"
#include "../globals.h"
#include "../orbits.h"
#include <time.h>

/* very small timegm since we cannot rely on the libc one */
static time_t tiny_timegm(struct tm *t)
{
    static const int dim[12] = {31,28,31,30,31,30,31,31,30,31,30,31};
    int y = t->tm_year + 1900;
    int m = t->tm_mon;
    int d = t->tm_mday - 1;             /* 0-based */

    int days = 0;
    for(int yr=1970; yr<y; ++yr)
        days += 365 + ((yr%4==0 && yr%100) || (yr%400==0));
    for(int mo=0; mo<m; ++mo){
        days += dim[mo];
        if(mo==1 && ((y%4==0 && y%100) || (y%400==0)))
            ++days;                     /* Feb 29 */
    }
    days += d;
    return (time_t)days*86400 + t->tm_hour*3600 + t->tm_min*60 + t->tm_sec;
}

/* convert GPS week/time to BDT week/time (BDT lags GPS by 14 s, epoch 1356w) */
static void gps_to_bdt(int gps_w, double gps_sow, int *bdt_w, double *bdt_sow)
{
    int w = gps_w - 1356;
    double sow = gps_sow - 14.0;        /* BDT is 14 s behind GPS */
    if(sow < 0){ sow += 604800.0; --w; }
    *bdt_w = w; *bdt_sow = sow;
}

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
    int gps_w = 0;          /* current GPS week from SP3 */
    double gps_sow = 0.0;   /* current GPS second-of-week */

    int count = 0;          /* number of samples compared */
    double err2_sum = 0.0;  /* squared error accumulator  */
    bool printed[MAX_SAT] = {0};

    while(fgets(l, sizeof l, fp))
    {
        if(l[0] == '*'){
            int Y,M,D,h,m; double s;
            sscanf(l+1, "%d %d %d %d %d %lf", &Y,&M,&D,&h,&m,&s);
            struct tm utc = {0};
            utc.tm_year=Y-1900; utc.tm_mon=M-1; utc.tm_mday=D;
            utc.tm_hour=h;     utc.tm_min=m;   utc.tm_sec=(int)s;
            time_t t = tiny_timegm(&utc);
            const time_t GPS_EPOCH = 315964800; /* 1980-01-06 */
            double sec = difftime(t, GPS_EPOCH) + s - (int)s; /* include fractional */
            gps_w  = (int)(sec / 604800.0);
            gps_sow= sec - gps_w*604800.0;
            continue;
        }
        if(l[0]=='P' && l[1]=='C'){
            int prn; double x,y,z,clk;
            if(sscanf(l+1, "C%2d %lf %lf %lf %lf", &prn,&x,&y,&z,&clk)!=5) continue;
            if(prn<1 || prn>63 || eph[prn].prn==0) continue;
            if(!((prn>=1 && prn<=5) || (prn>=59 && prn<=63)))
                continue; /* only GEO */

            if(!printed[prn]){
                printf("PRN%02d sqrtA=%.3f omegadot=%g\n",
                       prn, eph[prn].sqrtA, eph[prn].omegadot);
                printed[prn] = true;
            }

            int bw; double bsow;
            gps_to_bdt(gps_w, gps_sow, &bw, &bsow);

            double xyz[3];
            calc_sat_position_velocity(prn, bw, bsow, xyz, NULL);

            double dx = xyz[0] - x*1000.0;
            double dy = xyz[1] - y*1000.0;
            double dz = xyz[2] - z*1000.0;
            double err = sqrt(dx*dx + dy*dy + dz*dz);
            if(count < 3){
                printf("%s bw=%d bsow=%.1f PRN%02d err=%.3f\n",
                       l+1, bw, bsow, prn, err);
            }
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
