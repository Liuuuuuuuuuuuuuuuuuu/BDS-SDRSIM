#ifndef ORBITS_H
#define ORBITS_H
#include "bdssim.h"
void calc_sat_position_velocity(int prn,int week,double sow,
                                double *xyz,double *vel);
void eph_pos_vel_ecef(const ephemeris_t *ep,int week,double sow,
                      double *xyz,double *vel);
#endif



