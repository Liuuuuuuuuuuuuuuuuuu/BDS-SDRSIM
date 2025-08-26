#ifndef IONO_H
#define IONO_H

typedef struct {
    double F;        /* slant factor */
    double psi;      /* earth-centered angle */
    double phi_i;    /* sub-ionospheric latitude (semicircle) */
    double lambda_i; /* sub-ionospheric longitude (semicircle) */
    double phi_m;    /* geomagnetic latitude (semicircle) */
    double t_local;  /* local time (s) */
    double I_vert;   /* vertical delay (s) */
    double I_slant;  /* slant delay (s) */
    double delay_m;  /* slant delay (m) */
} iono_res_t;

/* Compute Klobuchar ionospheric delay.  Returns slant-path delay in meters.
 * Add this value to the geometric range for pseudorange; subtract the
 * equivalent phase (delay_m / wavelength) from carrier phase measurements.
 */
double iono_delay(double lat_deg, double lon_deg,
                  double az_deg, double el_deg,
                  double sow, const double alpha[4], const double beta[4],
                  double f, double f_ref, iono_res_t *out);

#endif /* IONO_H */
