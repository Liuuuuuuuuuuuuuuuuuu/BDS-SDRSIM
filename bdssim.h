#ifndef BDS_SIM_H
#define BDS_SIM_H
/* ─────────────────────────────────────────────── */
#include <stdint.h>
#include <stdbool.h>
#include "coord.h"
#include "channel.h"           /* 取 CODE_LEN 宏 */

/* ---------- 常數 ---------- */
#define MAX_SAT      65
#define MAX_PRN_LEN  CODE_LEN
#define MAX_CH       12        /* 最多同時播 12 顆 */

/* ---------- 星曆結構 (僅用到的欄位) ---------- */
typedef struct {
    int    prn, week;
    double toe;
    double sqrtA, e, i0, omega0, w, M0;
    double deltan, idot, omegadot;
    double cuc, cus, cic, cis, crc, crs;
    double af0, af1, af2;

    /* 其他 RINEX 欄位（暫未使用） */
    int     aode;            /* Age of Data, Ephemeris (line 2) */
    double  a0utc, a1utc;    /* BDT‑UTC 偏移 */
    double  toe_msg;         /* Transmission time of message */
    double  tgd1, tgd2;      /* TGD1/TGD2 group delays */
    int     aodc;            /* Age of Data, Clock */
    double  reserved;        /* 末行保留值 */

    /* ─ 子幀 1 需要 ─ */
    int     toc;       /* 秒 (0–604799) */
    uint8_t ura;       /* 0–15 */
    uint8_t health;    /* 0 = 正常；1 = 不可用 */
} ephemeris_t;

/* ---------- 全域表 ---------- */
extern ephemeris_t eph[MAX_SAT];
extern uint8_t     prn_code[MAX_SAT][MAX_PRN_LEN];

/* ---------- CLI 設定 ---------- */
typedef struct {
    char   rinex_file[256];    /* RINEX nav 檔 */
    char   time_start[32];     /* UTC 字串 */
    /* 靜態使用者位置 (deg,deg,m) ------------------ */
    double llh[3];
    char   path_file[256];     /* 動態路徑檔案 */
    int    path_type;          /* 0: static, 1:xyz,2:llh,3:nmea */

    /* 其他選項 */
    uint32_t duration;         /* 模擬秒數 */
    uint32_t step_ms;          /* 幾何更新粒度 (ms) */
    double   gain;             /* 輸出增益 */
    double   target_cn0;       /* 目標 CN0 (dB-Hz) */
    unsigned seed;             /* RNG seed */
    bool     byte_output;      /* 以 8-bit 檔輸出 */
    bool     geo_first;        /* 優先可見 GEO 衛星 */
    bool     no_geo;           /* 排除 GEO 衛星 */
    int      single_prn;       /* 僅模擬此 PRN (0 = 全部) */
} sim_config_t;

/* ---------- 介面 ---------- */
bool  init_simulator(sim_config_t *, double start_bdt);
void  generate_signal(const sim_config_t *cfg);   /* ← 加上 const */
void  cleanup_simulator(void);
/* 讓 main 可以先做一次 select_channels() 取得 PRN 清單 */
int  select_channels(channel_t *ch, int *n_ch, const coord_t *usr,
                     bool geo_first, int single_prn, bool no_geo);
#endif
