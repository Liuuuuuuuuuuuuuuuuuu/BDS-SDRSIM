#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include "path.h"

static void trim(char *s)
{
    char *e = s + strlen(s);
    while (e > s && isspace((unsigned char)e[-1]))
        *--e = '\0';
    while (*s && isspace((unsigned char)*s))
        memmove(s, s + 1, e - s + 1);
}

static int count_lines(FILE *fp)
{
    int n = 0;
    char buf[256];
    while (fgets(buf, sizeof(buf), fp))
        if (buf[0] && buf[0] != '#' && buf[0] != '\n')
            n++;
    fseek(fp, 0, SEEK_SET);
    return n;
}
typedef int (*parse_fn)(const char *buf, coord_t *out);

static int parse_xyz(const char *buf, coord_t *out)
{
    double x, y, z;
    if (sscanf(buf, "%lf%lf%lf", &x, &y, &z) != 3)
        return -1;
    memset(out, 0, sizeof(*out));
    out->xyz[0] = x;
    out->xyz[1] = y;
    out->xyz[2] = z;
    return 0;
}

static int parse_llh(const char *buf, coord_t *out)
{
    double lat, lon, h;
    if (sscanf(buf, "%lf%lf%lf", &lat, &lon, &h) != 3)
        return -1;
    llh2xyz((double[3]){lat, lon, h}, out);
    return 0;
}

static int parse_nmea_gga(const char *buf, coord_t *out)
{
    double rawlat = 0.0, rawlon = 0.0, alt = 0.0;
    char ns = 'N', ew = 'E';
    int n = sscanf(buf,
                   "$%*[^,],%*[^,],%lf,%c,%lf,%c,%*[^,],%*[^,],%*[^,],%lf",
                   &rawlat, &ns, &rawlon, &ew, &alt);
    if (n < 5)
        return -1;

    int d = (int)(rawlat / 100);
    double lat = d + (rawlat - d * 100) / 60.0;
    if (ns == 'S')
        lat = -lat;

    d = (int)(rawlon / 100);
    double lon = d + (rawlon - d * 100) / 60.0;
    if (ew == 'W')
        lon = -lon;

    llh2xyz((double[3]){lat, lon, alt}, out);
    return 0;
}

static int load_path_common(const char *file, path_t *path, parse_fn parse)
{
    FILE *fp = fopen(file, "r");
    if (!fp)
        return -1;

    int n = count_lines(fp);
    path->p = malloc(sizeof(*path->p) * n);
    path->n = 0;

    char buf[256];
    while (fgets(buf, sizeof(buf), fp)) {
        trim(buf);
        if (buf[0] == '#' || buf[0] == '\0')
            continue;
        coord_t c;
        if (parse(buf, &c) == 0)
            path->p[path->n++] = c;
    }

    fclose(fp);
    return path->n == n ? 0 : -1;
}

int load_path_xyz(const char *file, path_t *path)
{
    return load_path_common(file, path, parse_xyz);
}

int load_path_llh(const char *file, path_t *path)
{
    return load_path_common(file, path, parse_llh);
}

int load_path_nmea(const char *file, path_t *path)
{
    return load_path_common(file, path, parse_nmea_gga);
}

void free_path(path_t *path)
{
    if (path->p) {
        free(path->p);
        path->p = NULL;
    }
    path->n = 0;
}

void interpolate_path(const path_t *path, double t, coord_t *out)
{
    if (!path->p || path->n == 0) {
        memset(out, 0, sizeof(coord_t));
        return;
    }

    if (t <= 0) {
        *out = path->p[0];
        return;
    }

    int i = (int)t;
    if (i >= path->n - 1) {
        *out = path->p[path->n - 1];
        return;
    }

    double f = t - i;
    for (int k = 0; k < 3; ++k)
        out->xyz[k] = path->p[i].xyz[k] * (1.0 - f) +
                      path->p[i + 1].xyz[k] * f;
}

void interpolate_path_llh(const path_t *path,double t,double llh[3])
{
    if(!path->p || path->n==0){
        llh[0]=llh[1]=llh[2]=0.0;
        return;
    }
    if(t<=0){
        memcpy(llh,path->p[0].llh,sizeof(double)*3);
        return;
    }
    int i=(int)t;
    if(i>=path->n-1){
        memcpy(llh,path->p[path->n-1].llh,sizeof(double)*3);
        return;
    }
    double f=t-i;
    for(int k=0;k<3;++k)
        llh[k]=path->p[i].llh[k]*(1.0-f)+path->p[i+1].llh[k]*f;
}

