#ifndef PATH_H
#define PATH_H
#include "coord.h"

typedef struct {
    int n;           /* number of 1 Hz points */
    coord_t *p;      /* array of coordinates (ECEF) */
} path_t;

int load_path_xyz(const char *file, path_t *path);
int load_path_llh(const char *file, path_t *path);
int load_path_nmea(const char *file, path_t *path);
void free_path(path_t *path);
void interpolate_path(const path_t *path,double t,coord_t *out);

#endif
