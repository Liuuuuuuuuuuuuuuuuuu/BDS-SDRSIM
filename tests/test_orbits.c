/*
 *  讀取 SP3 星曆，驗證 calc_sat_position_velocity() 計算的
 *  ECEF 位置與速度。為了從 SP3 求得速度，需讀取相鄰五個
 *  epoch，利用五點中央差分估算。比較時略過 GEO 衛星。
 */

#include <stdio.h>
#include <stdbool.h>
#include <string.h>
#include "../bdssim.h"
#include "../globals.h"
#include "../orbits.h"
#include "../timeconv.h"

#define MAX_EPOCHS 5
#define MAX_PRN    64

int main(void)
{
    FILE *fp = fopen("GBM0MGXRAP_20251760000_01D_05M_ORB.SP3", "r");
    if (!fp) { perror("SP3"); return 1; }

    char  l[256];
    int   weeks[MAX_EPOCHS];
    double sows[MAX_EPOCHS];
    double sp3_pos[MAX_EPOCHS][MAX_PRN][3] = {{{0}}};
    int   epoch_idx = -1;           /* -1 ⇒ 尚未遇到任何 epoch */

    char current_utc[32] = "";     /* 取中央 epoch 的原始時間字串 */

    while (fgets(l, sizeof l, fp)) {
        if (l[0] == '*') {
            if (epoch_idx + 1 >= MAX_EPOCHS)
                break;              /* 只需前 5 個 epoch */

            ++epoch_idx;
            int Y, M, D, h, m; double s;
            sscanf(l + 1, "%d %d %d %d %d %lf", &Y, &M, &D, &h, &m, &s);
            int    is   = (int)s;
            double frac = s - is;

            /* 原始 SP3 epoch (GPS time) */
            if (epoch_idx == 2)
                snprintf(current_utc, sizeof current_utc,
                         "%04d/%02d/%02d,%02d:%02d:%02d",
                         Y, M, D, h, m, is);

            /* GPS → UTC → BDT */
            struct tm tm = {0};
            tm.tm_year = Y - 1900;
            tm.tm_mon  = M - 1;
            tm.tm_mday = D;
            tm.tm_hour = h;
            tm.tm_min  = m;
            tm.tm_sec  = is;

            time_t t = timegm(&tm) - 18; /* GPS-UTC = 18s */
            struct tm tm2;
            gmtime_r(&t, &tm2);

            char utc[64];
            snprintf(utc, sizeof utc, "%04d/%02d/%02d,%02d:%02d:%02d",
                     tm2.tm_year + 1900, tm2.tm_mon + 1, tm2.tm_mday,
                     tm2.tm_hour, tm2.tm_min, tm2.tm_sec);

            int bw = 0; double bsow = 0.0;
            if (utc_to_bdt(utc, &bw, &bsow) != 0) {
                fprintf(stderr, "utc_to_bdt failed\n");
                return 1;
            }
            bsow += frac;           /* fractional seconds */

            weeks[epoch_idx] = bw;
            sows [epoch_idx] = bsow;

            if (epoch_idx == 2) {
                /* 使用中央 epoch 初始化模擬器 */
                double start_bdt = bw * 604800.0 + bsow;
                sim_config_t cfg = {0};
                snprintf(cfg.rinex_file, sizeof(cfg.rinex_file),
                         "BRDM00DLR_S_20251760000_01D_MN.rnx");
                if (!init_simulator(&cfg, start_bdt)) {
                    fprintf(stderr, "init failed\n");
                    return 1;
                }
            }
            continue;
        }

        if (epoch_idx == -1)
            continue;               /* skip header before first epoch */

        if (l[0] == 'P' && l[1] == 'C') {
            int prn; double x, y, z, clk;
            if (sscanf(l + 1, "C%2d %lf %lf %lf %lf",
                       &prn, &x, &y, &z, &clk) != 5)
                continue;
            if (prn < 1 || prn > 63) continue;
            if (prn <= 5 || (prn >= 59 && prn <= 63))
                continue;          /* skip GEO satellites */

            if (epoch_idx < MAX_EPOCHS) {
                sp3_pos[epoch_idx][prn][0] = x * 1000.0;
                sp3_pos[epoch_idx][prn][1] = y * 1000.0;
                sp3_pos[epoch_idx][prn][2] = z * 1000.0;
            }
        }
    }
    fclose(fp);

    if (epoch_idx < 4) {
        fprintf(stderr, "SP3 file needs at least 5 epochs\n");
        return 1;
    }

    const double h = (weeks[3] * 604800.0 + sows[3]) -
                     (weeks[2] * 604800.0 + sows[2]);

    for (int prn = 1; prn <= 63; ++prn) {
        if (sp3_pos[2][prn][0] == 0.0 &&
            sp3_pos[2][prn][1] == 0.0 &&
            sp3_pos[2][prn][2] == 0.0)
            continue;               /* 無此衛星資料 */
        if ((sp3_pos[0][prn][0] == 0.0 &&
             sp3_pos[0][prn][1] == 0.0 &&
             sp3_pos[0][prn][2] == 0.0) ||
            (sp3_pos[1][prn][0] == 0.0 &&
             sp3_pos[1][prn][1] == 0.0 &&
             sp3_pos[1][prn][2] == 0.0) ||
            (sp3_pos[3][prn][0] == 0.0 &&
             sp3_pos[3][prn][1] == 0.0 &&
             sp3_pos[3][prn][2] == 0.0) ||
            (sp3_pos[4][prn][0] == 0.0 &&
             sp3_pos[4][prn][1] == 0.0 &&
             sp3_pos[4][prn][2] == 0.0))
            continue;               /* 缺資料無法估算速度 */

        double xyz[3], vel[3];
        calc_sat_position_velocity(prn, weeks[2], sows[2], xyz, vel, NULL);

        double sp3_vel[3];
        sp3_vel[0] = (-sp3_pos[4][prn][0] + 8 * sp3_pos[3][prn][0]
                      - 8 * sp3_pos[1][prn][0] + sp3_pos[0][prn][0]) / (12 * h);
        sp3_vel[1] = (-sp3_pos[4][prn][1] + 8 * sp3_pos[3][prn][1]
                      - 8 * sp3_pos[1][prn][1] + sp3_pos[0][prn][1]) / (12 * h);
        sp3_vel[2] = (-sp3_pos[4][prn][2] + 8 * sp3_pos[3][prn][2]
                      - 8 * sp3_pos[1][prn][2] + sp3_pos[0][prn][2]) / (12 * h);

        printf("time=%s PRN%02d calc=(%.3f %.3f %.3f) sp3=(%.3f %.3f %.3f)\n",
               current_utc, prn, xyz[0], xyz[1], xyz[2],
               sp3_pos[2][prn][0], sp3_pos[2][prn][1], sp3_pos[2][prn][2]);
        printf("            vcalc=(%.6f %.6f %.6f) vsp3=(%.6f %.6f %.6f)\n",
               vel[0], vel[1], vel[2],
               sp3_vel[0], sp3_vel[1], sp3_vel[2]);
    }

    return 0;
}
