#include <math.h>
#include "bdssim.h"

/* WGS-84 / GNSS 常數 */
#define GM        3.986004418e14
#define OMEGA_E   7.2921150e-5

/* ---------- Kepler 方程 ---------- */
static double solve_kepler(double M,double e)
{
    double E=M, d;
    do{
        d=(E - e*sin(E) - M)/(1 - e*cos(E));
        E-=d;
    }while(fabs(d)>1e-13);
    return E;
}

/* ---------- 位置（共用） ---------- */
static void sat_pos(const ephemeris_t *ep, double tk, double *xyz)
{
    if(ep->toe==0) tk = 0;           /* 防呆：缺 toe 則視為同步 */
    
    double A   = ep->sqrtA*ep->sqrtA;
    double n0  = sqrt(GM/(A*A*A));
    double n   = n0 + ep->deltan;
    double M   = ep->M0 + n*tk;

    double E = solve_kepler(M, ep->e);
    double sinE = sin(E), cosE = cos(E);
    double v = atan2(sqrt(1-ep->e*ep->e)*sinE, cosE-ep->e);
    double phi = v + ep->w;

    double du = ep->cus*sin(2*phi) + ep->cuc*cos(2*phi);
    double dr = ep->crs*sin(2*phi) + ep->crc*cos(2*phi);
    double di = ep->cis*sin(2*phi) + ep->cic*cos(2*phi);

    double u = phi + du;
    double r = A*(1-ep->e*cosE) + dr;
    double i = ep->i0 + di + ep->idot*tk;

    double x_op = r*cos(u);
    double y_op = r*sin(u);

    double omega = ep->omega0 + (ep->omegadot-OMEGA_E)*tk
                   - OMEGA_E*ep->toe;

    double cosO = cos(omega), sinO = sin(omega);
    double cosi = cos(i),     sini = sin(i);

    xyz[0] = x_op*cosO - y_op*cosi*sinO;
    xyz[1] = x_op*sinO + y_op*cosi*cosO;
    xyz[2] = y_op*sini;
}

/* ---------- 主要 API ---------- */
void calc_sat_position_velocity(int prn,int week,double sow,
                                double *xyz,double *vel)
{
    (void)week;                       /* 週號待用 */
    const ephemeris_t *ep=&eph[prn];
    if(ep->prn==0){xyz[0]=xyz[1]=xyz[2]=0; if(vel) vel[0]=vel[1]=vel[2]=0; return;}

    /* 時間差 tk（含跨週折返） */
    double tk = sow - ep->toe;
    if(tk> 302400) tk-=604800;
    if(tk<-302400) tk+=604800;

    /* 位置 */
    sat_pos(ep, tk, xyz);

    if(!vel) return;

    /* 速度：中央差分 Δt = 0.1 s */
    const double dt = 0.1;
    double p_f[3], p_b[3];
    sat_pos(ep, tk+dt, p_f);
    sat_pos(ep, tk-dt, p_b);

    vel[0] = (p_f[0]-p_b[0])/(2*dt);
    vel[1] = (p_f[1]-p_b[1])/(2*dt);
    vel[2] = (p_f[2]-p_b[2])/(2*dt);
}

