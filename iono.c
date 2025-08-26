#include <math.h>
#include <float.h>
#include "globals.h"  /* CLIGHT, M_PI */
#include "iono.h"

/* Compute ionospheric delay using Klobuchar model.
 * Inputs are in degrees and seconds; alpha/beta coefficients follow
 * broadcast convention.  Frequencies are in Hz.
 * Returns slant delay in meters and fills optional debug structure.
 */

double iono_delay(double lat_deg, double lon_deg,
                  double az_deg, double el_deg,
                  double sow, const double alpha[4], const double beta[4],
                  double f, double f_ref, iono_res_t *out)
{
    /* Convert degrees to semicircles */
    double phi_u = lat_deg / 180.0;
    double lam_u = lon_deg / 180.0;
    double az    = az_deg  / 180.0;
    double el    = el_deg  / 180.0;

    /* Clamp elevation to [0,0.9] semicircles (0°–162°) */
    if (el < 0.0) el = 0.0;
    else if (el > 0.9) el = 0.9;

    /* Wrap SOW to [0,86400) */
    sow = fmod(sow, 86400.0);
    if (sow < 0.0) sow += 86400.0;

    /* 1. Slant factor F */
    double tmp = 0.53 - el;
    double F = 1.0 + 16.0 * tmp * tmp * tmp;

    /* 2. Earth-centered angle psi */
    double psi = 0.0137 / (el + 0.11) - 0.022;

    /* 3. Sub-ionospheric latitude phi_i */
    double phi_i = phi_u + psi * cos(M_PI * az);
    if (phi_i > 0.416) phi_i = 0.416;
    if (phi_i < -0.416) phi_i = -0.416;

    /* 4. Sub-ionospheric longitude lambda_i */
    double cos_phi = cos(M_PI * phi_i);
    if (fabs(cos_phi) < 1e-12)
        cos_phi = cos_phi >= 0.0 ? 1e-12 : -1e-12;
    double lambda_i = lam_u + (psi * sin(M_PI * az)) / cos_phi;

    /* Normalize lambda_i to [-1,+1] semicircles */
    if (lambda_i > 1.0 || lambda_i < -1.0)
        lambda_i -= floor((lambda_i + 1.0) / 2.0) * 2.0;

    /* 5. Geomagnetic latitude phi_m */
    double phi_m = phi_i + 0.064 * cos(M_PI * (lambda_i - 1.617));

    /* 6. Local time t_local */
    double t_local = 43200.0 * lambda_i + sow;
    t_local = fmod(t_local, 86400.0);
    if (t_local < 0.0) t_local += 86400.0;

    /* 7. Amplitude A */
    double phi_m2 = phi_m * phi_m;
    double phi_m3 = phi_m2 * phi_m;
    double A = alpha[0] + alpha[1]*phi_m + alpha[2]*phi_m2 + alpha[3]*phi_m3;
    if (A < 0.0) A = 0.0;

    /* 8. Period P */
    double P = beta[0] + beta[1]*phi_m + beta[2]*phi_m2 + beta[3]*phi_m3;
    if (P < 72000.0) P = 72000.0;

    /* 9. Phase x */
    double x = 2.0 * M_PI * (t_local - 50400.0) / P;

    /* 10. Vertical delay I_vert */
    double I_vert;
    if (fabs(x) < (M_PI/2))
        I_vert = 5e-9 + A * (1.0 - x*x/2.0 + (x*x*x*x)/24.0);
    else
        I_vert = 5e-9;

    /* 11. Slant delay seconds */
    double I_slant = F * I_vert;

    /* 12. Frequency mapping */
    if (f > 0.0 && f_ref > 0.0) {
        double ratio = f_ref / f;
        I_slant *= ratio * ratio;
    }

    /* 13. Convert to meters */
    double delay_m = CLIGHT * I_slant;

    if (out) {
        out->F = F;
        out->psi = psi;
        out->phi_i = phi_i;
        out->lambda_i = lambda_i;
        out->phi_m = phi_m;
        out->t_local = t_local;
        out->I_vert = I_vert;
        out->I_slant = I_slant;
        out->delay_m = delay_m;
    }
    return delay_m;
}
