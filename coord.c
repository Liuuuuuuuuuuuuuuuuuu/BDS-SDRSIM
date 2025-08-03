#include "coord.h"
#include "globals.h"     /* nav_week */
/* 地球自轉角速度 (rad/s) */
#define OMEGA_E   7.2921150e-5

/* LLH → ECEF，並把 llh/xyz 一併寫回 coord_t -------------- */
void llh2xyz(const double llh_deg[3], coord_t *c)
{
    const double lat = llh_deg[0] * (M_PI/180.0);
    const double lon = llh_deg[1] * (M_PI/180.0);
    const double  h  = llh_deg[2];

    /* 儲存原始 llh (deg) 方便之後 ecef2enu 使用 */
    c->llh[0] = llh_deg[0];
    c->llh[1] = llh_deg[1];
    c->llh[2] = llh_deg[2];

    /* WGS-84 精確轉換 ------------------------------------------------ */
    const double sinp = sin(lat);
    const double N = WGS_A / sqrt(1.0 - WGS_E2 * sinp * sinp);

    c->xyz[0] = (N + h) * cos(lat) * cos(lon);
    c->xyz[1] = (N + h) * cos(lat) * sin(lon);
    c->xyz[2] = (N * (1.0 - WGS_E2) + h) * sinp;
}

/* ECEF → LLH ------------------------------------------------------- */
void xyz2llh(const double xyz[3], coord_t *c)
{
    const double x=xyz[0], y=xyz[1], z=xyz[2];
    const double a=WGS_A, e2=WGS_E2;
    double lon=atan2(y,x);
    double p=hypot(x,y);
    double lat=atan2(z,p*(1.0-e2));
    for(int i=0;i<5;i++){
        double N=a/sqrt(1.0-e2*sin(lat)*sin(lat));
        double h=p/cos(lat)-N;
        lat=atan2(z,p*(1.0-e2*N/(N+h)));
    }
    double N=a/sqrt(1.0-e2*sin(lat)*sin(lat));
    double h=p/cos(lat)-N;
    c->llh[0]=lat*(180.0/M_PI);
    c->llh[1]=lon*(180.0/M_PI);
    c->llh[2]=h;
    c->xyz[0]=x; c->xyz[1]=y; c->xyz[2]=z;
}

/* ECEF → ENU ------------------------------------------------------- */
void ecef2enu(const coord_t *usr, const double sat_xyz[3], double enu[3])
{
    const double lat = usr->llh[0] * (M_PI/180.0);
    const double lon = usr->llh[1] * (M_PI/180.0);

    const double dx = sat_xyz[0] - usr->xyz[0];
    const double dy = sat_xyz[1] - usr->xyz[1];
    const double dz = sat_xyz[2] - usr->xyz[2];

    const double sinLat = sin(lat),  cosLat = cos(lat);
    const double sinLon = sin(lon),  cosLon = cos(lon);

    enu[0] = -sinLon * dx +  cosLon * dy;                           /* E */
    enu[1] = -cosLon*sinLat*dx - sinLon*sinLat*dy + cosLat*dz;      /* N */
    enu[2] =  cosLon*cosLat*dx + sinLon*cosLat*dy + sinLat*dz;      /* U */
}

/* 由 ENU 向量計算仰角（deg） --------------------------------------- */
double enu_elevation_deg(const double enu[3])
{
    const double rho = sqrt(enu[0]*enu[0] + enu[1]*enu[1] + enu[2]*enu[2]);
    if (rho == 0.0) return -90.0;   /* 防除以零 */
    return asin( enu[2] / rho ) * 180.0 / M_PI;
}

/* ECEF → ECI ------------------------------------------------------- */
void ecef_to_eci(const double ecef[3],int week,double sow,double eci[3])
{
    double t = (week - nav_week)*604800.0 + sow;
    double theta = -OMEGA_E * t;
    double cosT = cos(theta), sinT = sin(theta);
    eci[0] = cosT*ecef[0] - sinT*ecef[1];
    eci[1] = sinT*ecef[0] + cosT*ecef[1];
    eci[2] = ecef[2];
}

/* Rotate fixed LLH with Earth rotation ------------------------------- */
void static_user_at(int week,double sow,const coord_t *ref,
                    coord_t *out,double vel[3])
{
    double lat = ref->llh[0]*(M_PI/180.0);
    double lon0= ref->llh[1]*(M_PI/180.0);
    double h   = ref->llh[2];
    double sinp= sin(lat); double cosp=cos(lat);
    double N   = WGS_A/sqrt(1.0-WGS_E2*sinp*sinp);
    double t   = (week-nav_week)*604800.0 + sow;
    double lon = lon0 + OMEGA_E*t;
    double cosL= cos(lon), sinL= sin(lon);
    out->xyz[0]=(N+h)*cosp*cosL;
    out->xyz[1]=(N+h)*cosp*sinL;
    out->xyz[2]=(N*(1.0-WGS_E2)+h)*sinp;
    out->llh[0]=ref->llh[0];
    out->llh[1]=ref->llh[1];
    out->llh[2]=ref->llh[2];
    out->week=week;
    out->sow=sow;
    if(vel){
        vel[0]=-OMEGA_E*out->xyz[1];
        vel[1]= OMEGA_E*out->xyz[0];
        vel[2]=0.0;
    }
}


