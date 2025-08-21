#ifndef COORD_H
#define COORD_H
#include <math.h>
#include "orbits.h"  /* WGS84 constants */

typedef struct {
    double llh[3];   /* [rad,rad,m]  使用者原始經緯高 */
    double xyz[3];   /* [m]          對應 ECEF */
    int    week;     /* BDT 周 & 秒，可由主程式填入 */
    double sow;
} coord_t;

/* 介面 ----------------------------------------------------- */
void   lla_to_ecef(const double lla[3], coord_t *c);
void   ecef_to_lla(const double xyz[3], coord_t *c);
void   ecef_to_enu(const coord_t *usr, const double sat_xyz[3], double enu[3]);
double enu_elevation_deg(const double enu[3]);
void   static_user_at(int week,double sow,const coord_t *ref,coord_t *out,double vel[3]);

#endif

