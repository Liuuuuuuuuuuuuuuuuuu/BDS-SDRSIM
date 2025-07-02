#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <getopt.h>

#include "bdssim.h"
#include "timeconv.h"
#include "globals.h"
#include "path.h"

/* ───────────────────────────── */
static void usage(const char *p)
{
    printf("用法: %s --rinex file --start YYYY/MM/DD,hh:mm:ss [選項]\n", p);
    puts("選項:");
    puts("  --llh lat,lon,h      固定使用者位置");
    puts("  --xyz file           ECEF 路徑檔");
    puts("  --llh-file file      LLH 路徑檔");
    puts("  --nmea file          NMEA 路徑檔");
    puts("  --duration sec       模擬秒數 (1-3600)");
    puts("  --gain amp           輸出增益 (>0)");
    puts("  --help               顯示本說明\n");
}

int main(int argc,char *argv[])
{
    /* 1. 預設參數 ----------------------------------- */
    sim_config_t cfg = {0};
    cfg.sample_rate = 5000000;
    cfg.gain        = 1.0;
    cfg.step_ms     = 1;
    cfg.duration    = 300;                /* 預設 300 秒 */

    /* 2. 解析 CLI ----------------------------------- */
    static struct option longopt[] = {
        {"rinex",    required_argument, 0, 'e'},
        {"start",    required_argument, 0, 't'},
        {"llh",      required_argument, 0, 'l'},
        {"xyz",      required_argument, 0, 'x'},
        {"llh-file", required_argument, 0, 'y'},
        {"nmea",     required_argument, 0, 'n'},
        {"duration", required_argument, 0, 'd'},
        {"gain",     required_argument, 0, 'g'},
        {"help",     no_argument,       0, 'h'},
        {0,0,0,0}
    };

    int c, idx;
    while((c = getopt_long(argc, argv, "e:t:l:x:y:n:d:g:h", longopt, &idx)) != -1){
        switch(c){
        case 'e':
            strncpy(cfg.rinex_file, optarg, sizeof(cfg.rinex_file)-1);
            break;
        case 't':
            strncpy(cfg.time_start, optarg, sizeof(cfg.time_start)-1);
            break;
        case 'l': {                      /* lat,lon,h */
            double lat, lon, h;
            if(sscanf(optarg,"%lf,%lf,%lf",&lat,&lon,&h) != 3){
                fprintf(stderr,"--llh 格式應為 lat,lon,h\n");
                return 1;
            }
            cfg.llh[0] = lat;
            cfg.llh[1] = lon;
            cfg.llh[2] = h;
            break;
        }
        case 'x':
            strncpy(cfg.path_file,optarg,sizeof(cfg.path_file)-1);
            cfg.path_type=1;
            break;
        case 'y':
            strncpy(cfg.path_file,optarg,sizeof(cfg.path_file)-1);
            cfg.path_type=2;
            break;
        case 'n':
            strncpy(cfg.path_file,optarg,sizeof(cfg.path_file)-1);
            cfg.path_type=3;
            break;
        case 'd':
            cfg.duration = (uint32_t)atoi(optarg);
            break;
        case 'g':
            cfg.gain = atof(optarg);
            if(cfg.gain <= 0.0){
                fprintf(stderr,"--gain 必須 >0\n");
                return 1;
            }
            break;
        case 'h':
            usage(argv[0]);
            return 0;
        default:
            usage(argv[0]);
            return 1;
        }
    }

    if(cfg.rinex_file[0]=='\0' || cfg.time_start[0]=='\0'){
        usage(argv[0]); return 1;
    }

    /* --duration 1~3600 sec; --llh lat[-90,90], lon[-180,180], h[-1000,20000] */
    if(cfg.duration==0 || cfg.duration>3600){
        fputs("--duration 範圍 1~3600 秒\n", stderr);
        return 1;
    }
    if(cfg.llh[0]<-90.0 || cfg.llh[0]>90.0 ||
       cfg.llh[1]<-180.0 || cfg.llh[1]>180.0 ||
       cfg.llh[2]<-1000.0 || cfg.llh[2]>20000.0){
        fputs("--llh 超出合理範圍\n", stderr);
        return 1;
    }

    /* 3. 初始化 ------------------------------------- */
    if(!init_simulator(&cfg)){
        fprintf(stderr,"初始化失敗\n");
        return 1;
    }

    /* -- start 時間檢查：必須在星曆區間 h-1 ~ h+1 內 -- */
    int start_week; double start_sow;
    if(utc_to_bdt(cfg.time_start, &start_week, &start_sow)!=0){
        fputs("--start 格式錯誤\n", stderr);
        return 1;
    }
    double start_bdt = start_week*604800.0 + start_sow;
    if(start_bdt < nav_time_min - 3600.0 ||
       start_bdt + cfg.duration > nav_time_max + 3600.0){
        fputs("--start 超出星曆可用區間\n", stderr);
        return 1;
    }

    /* ---- 3-b. 使用者初始座標 ---- */
    coord_t usr;
    path_t path={0};
    if(cfg.path_type==1)       load_path_xyz(cfg.path_file,&path);
    else if(cfg.path_type==2)  load_path_llh(cfg.path_file,&path);
    else if(cfg.path_type==3)  load_path_nmea(cfg.path_file,&path);
    if(cfg.path_type!=0 && path.n==0){
        fputs("路徑檔讀取失敗\n",stderr); return 1;
    }
    if(cfg.path_type==0)      llh2xyz(cfg.llh,&usr);
    else {                    
        interpolate_path(&path,0.0,&usr);
        xyz2llh(usr.xyz,&usr);
    }
    usr.week = start_week;
    usr.sow  = start_sow;

    channel_t ch[MAX_CH];
    int n_ch;
    select_channels(ch,&n_ch,&usr);

    /* 印出確認訊息 (簡潔模式) */
    printf("[cfg] UTC %s  BDT W%d %.3f\n",
           cfg.time_start, usr.week, usr.sow);
    printf("[cfg] LLH %.6f %.6f %.1f\n",
           usr.llh[0], usr.llh[1], usr.llh[2]);
    printf("[cfg] XYZ %.3f %.3f %.3f (m)\n",
           usr.xyz[0], usr.xyz[1], usr.xyz[2]);
    printf("[cfg] PRN:");
    for(int i=0;i<n_ch;i++) printf(" %02d", ch[i].prn);
    printf("  Gain %.2f\n\n", cfg.gain);

    /* 4. 產生基帶 ----------------------------------- */
    generate_signal(&cfg);
    free_path(&path);

    /* 5. 結束 --------------------------------------- */
    cleanup_simulator();
    puts("[done] beidou_b1i.bin 已產生");
    return 0;
}

