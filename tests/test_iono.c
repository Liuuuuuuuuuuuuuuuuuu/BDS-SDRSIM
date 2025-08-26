#include <stdio.h>
#include <math.h>
#include "iono.h"

int main(void)
{
    double lat = 25.0, lon = 121.0;
    double az = 45.0, el = 30.0;
    double sow = 45000.0;
    double alpha[4] = {2.1420e-08, 5.9605e-08, -6.5565e-07, 1.0729e-06};
    double beta[4]  = {1.1264e5, 9.8304e4, 4.5875e5, -5.8982e5};
    double f = 1561.098e6, f_ref = 1561.098e6;
    iono_res_t out;
    double delay = iono_delay(lat, lon, az, el, sow, alpha, beta, f, f_ref, &out);
    if (fabs(delay - 6.6258123807814675) > 1e-6) {
        fprintf(stderr, "Unexpected delay: %.12f m\n", delay);
        return 1;
    }
    printf("Delay: %.6f m\n", delay);
    return 0;
}
