/* ---------------------------------------------------------------
 *  orbits.c  -  BeiDou / GPS-like satellite orbit propagation
 *  使用標準 GPS/BDS RAAN 公式：
 *      Ω(t) = Ω0 + (Ω̇ − ΩE)·tk − ΩE·Toe
 *  此式適用於 BDS 的 GEO、IGSO 與 MEO 衛星，無需依軌道類型
 *  切換不同計算方式，直接給出 ECEF 座標。
 *
 *  介面與輸入輸出陣列沿用舊版 calc_sat_position_velocity()
 * --------------------------------------------------------------*/

#include <stdio.h>
#include <math.h>
#include "bdssim.h"          /* 提供 ephemeris_t 及 eph[] 全域陣列 */

#define GM        3.986004418e14       /* WGS-84 μ  (m^3/s^2) */
#define OMEGA_E   7.2921150e-5         /* 地球自轉角速度 (rad/s) */

/* --------------------------------------------------------------*/
/*             小工具：解 Kepler 方程 E − e·sinE = M             */
static double kepler(double M, double e)
{
    double E = M, d;
    do {
        d = (E - e * sin(E) - M) / (1.0 - e * cos(E));
        E -= d;
    } while (fabs(d) > 1e-13);
    return E;
}

/* --------------------------------------------------------------*/
/*          主要 API：回傳 ECEF 位置 xyz 與速度 vel               */
void eph_pos_vel_ecef(const ephemeris_t *ep, int week, double sow,
                      double *xyz, double *vel)
{

    /*          輸出保底值 (PRN 不存在)                             */
    if (ep->prn == 0) {
        xyz[0] = xyz[1] = xyz[2] = 0.0;
        if (vel) vel[0] = vel[1] = vel[2] = 0.0;
        return;
    }


    /* -------- tk: 時間差 (ToE → 模擬時刻)，處理週界溢位 ------- */
    const double t_sim = week * 604800.0 + sow;
    const double t_toe = ep->week * 604800.0 + ep->toe;
    double tk = t_sim - t_toe;
    if (tk >  302400.0) tk -= 604800.0;
    if (tk < -302400.0) tk += 604800.0;

    /* -------- 基本軌道參數 & 平近點角 -------------------------- */
    const double A   = ep->sqrtA * ep->sqrtA;
    const double n0  = sqrt(GM / (A * A * A));
    const double n   = n0 + ep->deltan;
    const double M   = ep->M0 + n * tk;
    const double E   = kepler(M, ep->e);

    /* 真近點角 ν 與引數 φ = ν + ω                               */
    const double sinE = sin(E);
    const double cosE = cos(E);
    const double nu   = atan2(sqrt(1 - ep->e * ep->e) * sinE, cosE - ep->e);
    const double phi  = nu + ep->w;

    /* -------- 攝動修正 ---------------------------------------- */
    const double du = ep->cus * sin(2 * phi) + ep->cuc * cos(2 * phi);
    const double dr = ep->crs * sin(2 * phi) + ep->crc * cos(2 * phi);
    const double di = ep->cis * sin(2 * phi) + ep->cic * cos(2 * phi);

    const double u  = phi + du;
    const double r  = A * (1 - ep->e * cosE) + dr;
    const double i  = ep->i0 + di + ep->idot * tk;

    const double x_op = r * cos(u);
    const double y_op = r * sin(u);

    /* -------- RAAN Ω(t)：GPS/BDS 通用公式 --------------------- */
    const double Omega = ep->omega0 +
                         (ep->omegadot - OMEGA_E) * tk -
                         OMEGA_E * ep->toe;

    /* -------- ECI (WGS-84 inertial) 座標 ----------------------- */
    const double cosO = cos(Omega), sinO = sin(Omega);
    const double cosi = cos(i),     sini = sin(i);

    const double x = x_op * cosO - y_op * cosi * sinO;
    const double y = x_op * sinO + y_op * cosi * cosO;
    const double z = y_op * sini;

    /* ---- 若需要速度，先在 ECI 求導 --------------------------- */
    double x_dot = 0.0, y_dot = 0.0, z_dot = 0.0;
    if (vel) {
        const double Edot   = n / (1 - ep->e * cosE);
        const double nu_dot = sqrt(1 - ep->e * ep->e) * Edot / (1 - ep->e * cosE);
        const double phi_dot = nu_dot;
        const double u_dot   = phi_dot + 2 * (ep->cus * cos(2 * phi)
                                            - ep->cuc * sin(2 * phi)) * phi_dot;
        const double r_dot   = A * ep->e * sinE * Edot
                             + 2 * (ep->crs * cos(2 * phi)
                                  - ep->crc * sin(2 * phi)) * phi_dot;
        const double i_dot   = ep->idot + 2 * (ep->cis * cos(2 * phi)
                                             - ep->cic * sin(2 * phi)) * phi_dot;

        const double Omega_dot = ep->omegadot - OMEGA_E;

        const double x_op_dot = r_dot * cos(u) - r * u_dot * sin(u);
        const double y_op_dot = r_dot * sin(u) + r * u_dot * cos(u);

        x_dot = (x_op_dot - y_op * cosi * Omega_dot) * cosO
              - (x_op * Omega_dot + y_op_dot * cosi - y_op * sini * i_dot) * sinO;
        y_dot = (x_op_dot - y_op * cosi * Omega_dot) * sinO
              + (x_op * Omega_dot + y_op_dot * cosi - y_op * sini * i_dot) * cosO;
        z_dot = y_op_dot * sini + y_op * cosi * i_dot;
    }

    xyz[0] = x;
    xyz[1] = y;
    xyz[2] = z;

    if (vel) {
        vel[0] = x_dot;
        vel[1] = y_dot;
        vel[2] = z_dot;
    }
}

/* --------------------------------------------------------------*/
/* 主要 API：依 PRN 查表計算                                     */
void calc_sat_position_velocity(int prn, int week, double sow,
                                double *xyz, double *vel)
{
    const ephemeris_t *ep = &eph[prn];
    eph_pos_vel_ecef(ep, week, sow, xyz, vel);
}
/* ---------------------------  End  ------------------------------*/
