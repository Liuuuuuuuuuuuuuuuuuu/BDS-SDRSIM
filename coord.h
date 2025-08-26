#ifndef COORD_H
#define COORD_H

#include <math.h>
#include "orbits.h"

typedef struct {
    double llh[3];
    double xyz[3];
    int    week;
    double sow;
} coord_t;

void lla_to_ecef(const double lla[3], coord_t *c);
void ecef_to_lla(const double xyz[3], coord_t *c);
void ecef_to_enu(const coord_t *usr, const double sat_xyz[3], double enu[3]);
double enu_elevation_deg(const double enu[3]);
void static_user_at(int week, double sow, const coord_t *ref,
                    coord_t *out, double vel[3]);

#endif /* COORD_H */

