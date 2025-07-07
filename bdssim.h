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
    uint32_t sample_rate;      /* 取樣率 (Hz) */
    uint32_t duration;         /* 模擬秒數 */
    uint32_t step_ms;          /* 幾何更新粒度 (ms) */
    double   gain;             /* 輸出增益 */
    double   target_cn0;       /* 目標 CN0 (dB-Hz) */
    double   noise_std;        /* AWGN 標準差 (0 表示無) */
    unsigned noise_seed;       /* AWGN 亂數種子 */
    bool     byte_output;      /* 以 8-bit 檔輸出 */
    bool     enable_d2;        /* 啟用 D2 播放 */
} sim_config_t;

/* ---------- 介面 ---------- */
bool  init_simulator(sim_config_t *);
void  generate_signal(const sim_config_t *cfg);   /* ← 加上 const */
void  cleanup_simulator(void);
/* 讓 main 可以先做一次 select_channels() 取得 PRN 清單 */
int  select_channels(channel_t *ch, int *n_ch, const coord_t *usr);
#endif
