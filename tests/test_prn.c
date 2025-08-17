#include <stdio.h>
#include "channel.h"

int main(void){
    int geo_prn[] = {1,2,3,4,5,59,60,61,62,63};
    for(unsigned i=0;i<sizeof(geo_prn)/sizeof(geo_prn[0]);++i){
        if(!is_geo_prn(geo_prn[i])){
            printf("fail GEO %d\n", geo_prn[i]);
            return 1;
        }
    }
    for(int prn=6; prn<59; ++prn){
        if(is_geo_prn(prn)){
            printf("false GEO %d\n", prn);
            return 1;
        }
    }
    puts("PRN test passed");
    return 0;
}
