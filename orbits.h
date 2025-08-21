#ifndef ORBITS_H
#define ORBITS_H

/* WGS-84 and Earth-rotation constants */
#define OMEGA_E   7.2921151467e-5      /* Earth rotation rate (rad/s) */
#define WGS84_A   6378137.0            /* semi-major axis (m) */
#define WGS84_E2  6.69437999014e-3     /* first eccentricity squared */

/* Compute satellite position/velocity and clock bias/drift.
 *
 * All returned coordinates are in the ECEF frame at transmit time t_tx.
 */
void calc_sat_position_velocity(int prn, int week, double sow,
                                double *xyz, double *vel, double *clk);

#endif



