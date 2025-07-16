/* ---------------------------------------------------------------
 *  orbits.c  -  BeiDou / GPS-like satellite orbit propagation
 *
 *  - 增強：自動辨識 GEO / IGSO / MEO 三種軌道型別
 *  -    GEO/IGSO/MEO 均以慣性 RAAN (Ω0 + Ω̇·tk) 計算，再統一旋轉回
 *         ECEF。此方法可減少公式差異並避免重複扣除地球自轉角。
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
/*         軌道類型列舉 & 自動分類 (GEO / IGSO / MEO)           */
typedef enum { ORB_MEO = 0, ORB_IGSO = 1, ORB_GEO = 2 } orb_t;

static orb_t classify_orbit(const ephemeris_t *ep)
{
    const double A        = ep->sqrtA * ep->sqrtA;       /* semi-major axis (m) */
    const double inc_deg  = ep->i0 * (180.0 / M_PI);     /* inclination (°)     */

    /* 北斗：MEO A≈26.6e6 m；IGSO/GEO A≈42.2e6 m */
    if (A < 4.0e7) return ORB_MEO;                       /* 低於 4e7 為中軌 */

    /* A 落在 GEO/IGSO 區 → 用傾角判別 */
    if (fabs(inc_deg) < 15.0) return ORB_GEO;             /* 傾角 <15° 視為 GEO */
    return ORB_IGSO;                                     /* 其餘為 IGSO */
}

/* --------------------------------------------------------------*/
/*          主要 API：回傳 ECEF 位置 xyz 與速度 vel               */
void calc_sat_position_velocity(int prn, int week, double sow,
                                double *xyz, double *vel)
{
    const ephemeris_t *ep = &eph[prn];        /* 全域星曆陣列 (外部定義) */

    /*          輸出保底值 (PRN 不存在)                             */
    if (ep->prn == 0) {
        xyz[0] = xyz[1] = xyz[2] = 0.0;
        if (vel) vel[0] = vel[1] = vel[2] = 0.0;
        return;
    }

    /* -------- 軌道型別自動判斷 & 一次性列印 ------------------- */
    static int printed[64] = {0};
    const orb_t orb = classify_orbit(ep);
    if (!printed[prn]) {
        const char *tag = (orb == ORB_GEO)  ? "GEO"  :
                          (orb == ORB_IGSO) ? "IGSO" : "MEO";
        printf("[orbits] PRN C%02d classified as %s\n", prn, tag);
        printed[prn] = 1;
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

    /* -------- RAAN Ω(t)：統一以 ECI 公式處理 ------------------ */
    /*
     * 先在慣性座標計算 RAAN：Ω = Ω0 + Ω̇·tk
     * 之後再以地球自轉角 OMEGA_E·(ToE+tk) 旋回到 ECEF。
     * 這可避免重複扣除 OMEGA_E，亦能適用於 MEO/IGSO/GEO。  
     */
    const double Omega = ep->omega0 + ep->omegadot * tk;

    /* -------- ECI (WGS-84 inertial) 座標 ----------------------- */
    const double cosO = cos(Omega), sinO = sin(Omega);
    const double cosi = cos(i),     sini = sin(i);

    const double x_eci = x_op * cosO - y_op * cosi * sinO;
    const double y_eci = x_op * sinO + y_op * cosi * cosO;
    const double z_eci = y_op * sini;

    /* ---- 若需要速度，先在 ECI 求導 --------------------------- */
    double x_eci_dot = 0.0, y_eci_dot = 0.0, z_eci_dot = 0.0;
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

        /* RAAN 以慣性系統計，故其導數直接為 omegadot */
        const double Omega_dot = ep->omegadot;

        const double x_op_dot = r_dot * cos(u) - r * u_dot * sin(u);
        const double y_op_dot = r_dot * sin(u) + r * u_dot * cos(u);

        x_eci_dot = (x_op_dot - y_op * cosi * Omega_dot) * cosO
                  - (x_op * Omega_dot + y_op_dot * cosi - y_op * sini * i_dot) * sinO;
        y_eci_dot = (x_op_dot - y_op * cosi * Omega_dot) * sinO
                  + (x_op * Omega_dot + y_op_dot * cosi - y_op * sini * i_dot) * cosO;
        z_eci_dot = y_op_dot * sini + y_op * cosi * i_dot;
    }

    /* -------- 由 ECI 轉到 ECEF (全部軌道皆旋轉一次) -------- */
    const double theta = OMEGA_E * (ep->toe + tk);
    const double cosT  = cos(theta), sinT = sin(theta);

    {
        /* 全部軌道一律旋轉一次地球自轉角 */
        xyz[0] =  cosT * x_eci + sinT * y_eci;
        xyz[1] = -sinT * x_eci + cosT * y_eci;
        xyz[2] =  z_eci;

        if (vel) {
            vel[0] =  cosT * x_eci_dot + sinT * y_eci_dot + OMEGA_E * xyz[1];
            vel[1] = -sinT * x_eci_dot + cosT * y_eci_dot - OMEGA_E * xyz[0];
            vel[2] =  z_eci_dot;
        }
    }
}
/* ---------------------------  End  ------------------------------*/
