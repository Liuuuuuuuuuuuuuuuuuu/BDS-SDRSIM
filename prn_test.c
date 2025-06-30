#include <stdio.h>
#include "bdssim.h"
int main(void){
    sim_config_t c={0}; init_simulator(&c);
    puts("PRN01 前 32 chip:");
    for(int i=0;i<32;i++) printf("%d",prn_code[1][i]);
    puts(""); return 0;
}


