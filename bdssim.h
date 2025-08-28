#ifndef BDSSIM_H
#define BDSSIM_H

#include <stdint.h>
#include <stdbool.h>
#include "coord.h"
#include "channel.h"
#include "globals.h"
#include "rinex.h"

extern ephemeris_t eph[MAX_SAT];
extern uint8_t     prn_code[MAX_SAT][CODE_LEN];

typedef struct {
    char     rinex_file[256];
    char     time_start[32];
    double   llh[3];
    char     path_file[256];
    int      path_type;
    uint32_t duration;
    uint32_t step_ms;
    double   gain;
    double   target_cn0;
    unsigned seed;
    bool     byte_output;
    double   fs;           /* sample rate (Hz) */
    bool     meo_only;
    int      single_prn;
    bool     prn37_only;
    bool     iono_on;
} sim_config_t;

bool init_simulator(sim_config_t *, double start_bdt);
void generate_signal(const sim_config_t *cfg);
void cleanup_simulator(void);
int  select_channels(channel_t *ch, int *n_ch, const coord_t *usr,
                     int single_prn, bool meo_only);

#endif /* BDSSIM_H */
