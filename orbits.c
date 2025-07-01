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

/* ---------- 主要 API ---------- */
void calc_sat_position_velocity(int prn,int week,double sow,
                                double *xyz,double *vel)
{
    const ephemeris_t *ep=&eph[prn];
    if(ep->prn==0){
        xyz[0]=xyz[1]=xyz[2]=0;
        if(vel) vel[0]=vel[1]=vel[2]=0;
        return;
    }

    /* ---- 時間差 tk：採用 GPS 絕對秒數基準，避免週數交界問題 ---- */
    double simulated_time      = week*604800.0 + sow;
    double ephemeris_toe_time  = ep->week*604800.0 + ep->toe;

    double tk = simulated_time - ephemeris_toe_time;
    if(tk> 302400) tk-=604800;
    if(tk<-302400) tk+=604800;

    /* ---- 位置與速度 ---- */

    /* ---- 衛星速度（解析式，含地球自轉） ---- */
    /* 基本軌道參數 */
    double A    = ep->sqrtA*ep->sqrtA;
    double n0   = sqrt(GM/(A*A*A));
    double n    = n0 + ep->deltan;
    double M    = ep->M0 + n*tk;

    double E = solve_kepler(M, ep->e);
    double sinE = sin(E), cosE = cos(E);
    double nu  = atan2(sqrt(1-ep->e*ep->e)*sinE, cosE-ep->e);
    double phi = nu + ep->w;

    double du = ep->cus*sin(2*phi) + ep->cuc*cos(2*phi);
    double dr = ep->crs*sin(2*phi) + ep->crc*cos(2*phi);
    double di = ep->cis*sin(2*phi) + ep->cic*cos(2*phi);

    double u  = phi + du;
    double r  = A*(1-ep->e*cosE) + dr;
    double i  = ep->i0 + di + ep->idot*tk;

    double x_op = r*cos(u);        /* in orbital plane */
    double y_op = r*sin(u);

    double Omega0 = ep->omega0;
    double Omega  = Omega0 + ep->omegadot*tk;

    double cosO = cos(Omega), sinO = sin(Omega);
    double cosi = cos(i),     sini = sin(i);

    /* ECI position */
    double x_eci = x_op*cosO - y_op*cosi*sinO;
    double y_eci = x_op*sinO + y_op*cosi*cosO;
    double z_eci = y_op*sini;

    /* ECI velocity */
    double Edot = n/(1.0 - ep->e*cosE);
    double nu_dot  = sqrt(1-ep->e*ep->e)*Edot/(1.0 - ep->e*cosE);
    double phi_dot = nu_dot;
    double u_dot = phi_dot + 2.0*(ep->cus*cos(2*phi) - ep->cuc*sin(2*phi))*phi_dot;
    double r_dot = A*ep->e*sinE*Edot + 2.0*(ep->crs*cos(2*phi) - ep->crc*sin(2*phi))*phi_dot;
    double i_dot = ep->idot + 2.0*(ep->cis*cos(2*phi) - ep->cic*sin(2*phi))*phi_dot;

    double x_op_dot = r_dot*cos(u) - r*u_dot*sin(u);
    double y_op_dot = r_dot*sin(u) + r*u_dot*cos(u);

    double Omega_dot = ep->omegadot;

    double x_eci_dot = (x_op_dot - y_op*cosi*Omega_dot)*cosO
                       - (x_op*Omega_dot + y_op_dot*cosi - y_op*sini*i_dot)*sinO;
    double y_eci_dot = (x_op_dot - y_op*cosi*Omega_dot)*sinO
                       + (x_op*Omega_dot + y_op_dot*cosi - y_op*sini*i_dot)*cosO;
    double z_eci_dot = y_op_dot*sini + y_op*cosi*i_dot;

    /* ECEF transformation at t = ep->toe + tk */
    double t = ep->toe + tk;
    double theta = OMEGA_E*t;
    double cosT = cos(theta), sinT = sin(theta);

    xyz[0] =  cosT*x_eci + sinT*y_eci;
    xyz[1] = -sinT*x_eci + cosT*y_eci;
    xyz[2] =  z_eci;

    if(!vel) return;

    vel[0] =  cosT*x_eci_dot + sinT*y_eci_dot + OMEGA_E*xyz[1];
    vel[1] = -sinT*x_eci_dot + cosT*y_eci_dot - OMEGA_E*xyz[0];
    vel[2] =  z_eci_dot;
}

