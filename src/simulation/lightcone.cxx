#include <cassert>
#include <cmath>
#include <cstdio>

#include <algorithm>

#include "lightcone.h"

static inline double sq(double x) { return x*x; }
static inline double cb(double x) { return x*x*x; }

/* Calculate the quantity d\chi/dz = c/H(z) in Mpc/h */
double LightCone::cOverH(double z, double OmegaM) {
    return 2997.925/sqrt(OmegaM*cb(1+z) + 1-OmegaM);
}

LightCone::LightCone(double OmegaM, double x0_, double y0_, double z0_, double zmax) {
    x0 = x0_;
    y0 = y0_;
    z0 = z0_;
    amin = 1/(1 + zmax);
    anear = afar = amin;        // default to empty shell

    /* Compute comoving distance as a function of redshift */
    const int N = 10000;
    double z[N+1], r[N+1], a[N+1];
    z[0] = r[0] = 0;
    a[0] = 1;
    double dz = zmax/N;;
    double k1, k2, k3, k4;
    for(int i = 1; i <= N; i++) {
        /* Integrate dr/dz = c/H(z) using 4th order Runge-Kutta */
        z[i] = i*dz;
        k1 = cOverH(z[i], OmegaM);
        k2 = cOverH(z[i] + 0.5*dz, OmegaM);
        k3 = cOverH(z[i] + 0.5*dz, OmegaM);
        k4 = cOverH(z[i] + dz, OmegaM);
        r[i] = r[i-1] + (k1 + 2*k2 + 2*k3 + k4)*dz/6;
        a[i] = 1/(1 + z[i]);
    }

    /* Create spline for a(r) */
    aofr = CubicSpline(N+1, r, a);

    /* Create initial spline for r(a) */
    std::reverse(&a[0], &a[N+1]);
    std::reverse(&r[0], &r[N+1]);
    rofa = CubicSpline(N+1, a, r);

    /* Resample r(a) on uniform grid with linear interpolation */
    for(int i = 1; i < N; i++) {
        a[i] = amin + i*(1 - amin)/N;
        r[i] = rofa(a[i]);
    }
    rofa = LinearSpline(N+1, a, r);
}

void LightCone::SetOrigin(double x, double y, double z) {
    x0 = x;
    y0 = y;
    z0 = z;
}

void LightCone::DefineShell(double afar_, double anear_) {
    assert(amin <= afar_ && afar_ <= anear_ && anear_ <= 1);
    afar = afar_;
    anear = anear_;
    r2far = sq(rofa(afar));
    r2near = sq(rofa(anear));
}

bool LightCone::LiesInShell(double x, double y, double z) const {
    double rx = x - x0;
    double ry = y - y0;
    double rz = z - z0;
    double r2 = rx*rx + ry*ry + rz*rz;
    return (r2near <= r2) && (r2 < r2far);
}

void LightCone::Debug() const {
    FILE* fp = fopen("rofa.dat", "w");
    fprintf(fp, "# Comoving distance r(a) in Mpc/h\n");
    for(int i = 0; i <= 200; i++) {
        double a = amin + i*(1 - amin)/200;
        fprintf(fp, "%g %g\n", a, rofa(a));
    }

    fp = fopen("aofr.dat", "w");
    fprintf(fp, "# Inverse a(r)\n");
    double rmin = rofa(0.999);
    double rmax = rofa(amin);
    for(int i = 0; i <= 200; i++) {
        double r = rmin + i*(rmax - rmin)/200;
        fprintf(fp, "%g %g\n", r, aofr(r));
    }
}
