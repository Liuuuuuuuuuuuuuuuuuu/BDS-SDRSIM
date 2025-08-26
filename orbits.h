#ifndef ORBITS_H
#define ORBITS_H

#define OMEGA_E  7.2921151467e-5
#define WGS84_A  6378137.0
#define WGS84_E2 6.69437999014e-3

void calc_sat_position_velocity(int prn, int week, double sow,
                                double *xyz, double *vel, double *clk);

#endif /* ORBITS_H */



