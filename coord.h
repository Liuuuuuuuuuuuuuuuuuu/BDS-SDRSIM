#ifndef COORD_H
#define COORD_H
#include <math.h>

/* WGS-84 */
#define WGS_A   6378137.0
#define WGS_F   (1.0/298.257223563)
#define WGS_E2  (WGS_F*(2.0-WGS_F))

typedef struct {
    double llh[3];   /* [deg,deg,m]  使用者原始經緯高 */
    double xyz[3];   /* [m]          對應 ECEF */
    int    week;     /* BDT 周 & 秒，可由主程式填入 */
    double sow;
} coord_t;

/* 介面 ----------------------------------------------------- */
void   llh2xyz(const double llh_deg[3], coord_t *c);
void   xyz2llh(const double xyz[3], coord_t *c);
void   ecef2enu(const coord_t *usr, const double sat_xyz[3], double enu[3]);
double enu_elevation_deg(const double enu[3]);

#endif

