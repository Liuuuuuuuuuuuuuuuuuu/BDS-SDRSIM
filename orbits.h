#ifndef ORBITS_H
#define ORBITS_H
/* Compute satellite position, velocity and clock bias/drift at given time */
void calc_sat_position_velocity(int prn,int week,double sow,
                                double *xyz,double *vel,double *clk);
#endif



