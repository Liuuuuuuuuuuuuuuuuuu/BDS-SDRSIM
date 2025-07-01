#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include "path.h"

static void trim(char *s){
    char *e=s+strlen(s); while(e>s && isspace((unsigned char)e[-1])) *--e='\0';
    while(*s && isspace((unsigned char)*s)) memmove(s,s+1,e-s+1);
}

static int count_lines(FILE *fp){
    int n=0; char buf[256];
    while(fgets(buf,sizeof(buf),fp)) if(buf[0] && buf[0]!='#' && buf[0]!='\n') n++; 
    fseek(fp,0,SEEK_SET); return n;
}

static int parse_xyz(const char *buf,double *x,double *y,double *z){
    return sscanf(buf,"%lf%lf%lf",x,y,z)==3?0:-1;
}

int load_path_xyz(const char *file,path_t *path)
{
    FILE *fp=fopen(file,"r"); if(!fp) return -1;
    int n=count_lines(fp);
    path->p = malloc(sizeof(coord_t)*n); path->n=0;
    char buf[256];
    while(fgets(buf,sizeof(buf),fp)){
        trim(buf); if(buf[0]=='#' || buf[0]=='\0') continue;
        double x,y,z; if(parse_xyz(buf,&x,&y,&z)==0){
            coord_t c; memset(&c,0,sizeof(c));
            c.xyz[0]=x; c.xyz[1]=y; c.xyz[2]=z;
            path->p[path->n++]=c;
        }
    }
    fclose(fp);
    return path->n==n?0:-1;
}

int load_path_llh(const char *file,path_t *path)
{
    FILE *fp=fopen(file,"r"); if(!fp) return -1;
    int n=count_lines(fp);
    path->p = malloc(sizeof(coord_t)*n); path->n=0;
    char buf[256];
    while(fgets(buf,sizeof(buf),fp)){
        trim(buf); if(buf[0]=='#' || buf[0]=='\0') continue;
        double lat,lon,h; if(sscanf(buf,"%lf%lf%lf",&lat,&lon,&h)==3){
            coord_t c; llh2xyz((double[3]){lat,lon,h},&c);
            path->p[path->n++]=c;
        }
    }
    fclose(fp); return path->n==n?0:-1;
}

static int parse_nmea_gga(const char *buf,double *lat,double *lon,double *h)
{
    if(strncmp(buf,"$G",2)!=0) return -1;
    const char *p=strchr(buf,','); if(!p) return -1; /* time */
    p++; if(!*p) return -1; /* lat */
    double rawlat=0.0; if(sscanf(p,"%lf",&rawlat)!=1) return -1;
    p=strchr(p,','); if(!p) return -1; char ns='N'; sscanf(p+1,"%c",&ns);
    p=strchr(p+1,','); if(!p) return -1; double rawlon=0.0; if(sscanf(p+1,"%lf",&rawlon)!=1) return -1;
    p=strchr(p+1,','); if(!p) return -1; char ew='E'; sscanf(p+1,"%c",&ew);
    for(int i=0;i<4;i++){ p=strchr(p+1,','); if(!p) return -1; }
    if(sscanf(p+1,"%lf",h)!=1) *h=0.0;
    int d=(int)(rawlat/100); *lat=d+(rawlat-d*100)/60.0; if(ns=='S') *lat=-*lat;
    d=(int)(rawlon/100); *lon=d+(rawlon-d*100)/60.0; if(ew=='W') *lon=-*lon;
    return 0;
}

int load_path_nmea(const char *file,path_t *path)
{
    FILE *fp=fopen(file,"r"); if(!fp) return -1;
    int n=count_lines(fp);
    path->p=malloc(sizeof(coord_t)*n); path->n=0;
    char buf[256];
    while(fgets(buf,sizeof(buf),fp)){
        trim(buf); if(buf[0]=='#' || buf[0]=='\0') continue;
        double lat,lon,h; if(parse_nmea_gga(buf,&lat,&lon,&h)==0){
            coord_t c; llh2xyz((double[3]){lat,lon,h},&c);
            path->p[path->n++]=c;
        }
    }
    fclose(fp); return path->n==n?0:-1;
}

void free_path(path_t *path){ if(path->p){ free(path->p); path->p=NULL; } path->n=0; }

void interpolate_path(const path_t *path,double t,coord_t *out)
{
    if(!path->p || path->n==0){ memset(out,0,sizeof(coord_t)); return; }
    if(t<=0){ *out = path->p[0]; return; }
    int i=(int)t; if(i>=path->n-1){ *out=path->p[path->n-1]; return; }
    double f=t-i; for(int k=0;k<3;k++)
        out->xyz[k]=path->p[i].xyz[k]*(1.0-f)+path->p[i+1].xyz[k]*f;
}

