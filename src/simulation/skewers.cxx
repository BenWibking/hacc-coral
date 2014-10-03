#include <cassert>
#include <cerrno>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>

#include <algorithm>

#include "Domain.h"
#include "Partition.h"

#include "skewers.h"

using std::list;
using std::vector;
using std::min;
using std::max;

static inline double sq(double x) { return x*x; }
static inline double cb(double x) { return x*x*x; }


/***** SkewerGeometry *****/

/* XXX Check that this makes sense */
bool SkewerGeometry::Intersects(const Box& b) const {
    /* Box center */
    double cx = (b.xmax + b.xmin)/2;
    double cy = (b.ymax + b.ymin)/2;
    double cz = (b.zmax + b.zmin)/2;

    /* Box extents */
    double ax = (b.xmax - b.xmin)/2;
    double ay = (b.ymax - b.ymin)/2;
    double az = (b.zmax - b.zmin)/2;

    /* Midpoint of skewer relative to box center */
    double mx = sx + Ls*nx/2 - cx;
    double my = sy + Ls*ny/2 - cy;
    double mz = sz + Ls*nz/2 - cz;

    double Lx = -Ls*nx/2;
    double Ly = -Ls*ny/2;
    double Lz = -Ls*nz/2;
    double ex = fabs(Lx);
    double ey = fabs(Ly);
    double ez = fabs(Lz);

    /* Use separating axis test */
    if(fabs(mx) > ax + ex) return false;
    if(fabs(my) > ay + ey) return false;
    if(fabs(mz) > az + ez) return false;
    if(fabs(my*Lz - mz*Ly) > (ay*ez + az*ey)) return false;
    if(fabs(mx*Lz - mz*Lx) > (ax*ez + az*ex)) return false;
    if(fabs(mx*Ly - my*Lx) > (ax*ey + ay*ex)) return false;

    /* No separating axis, the line intersects */
    return true;
}


/***** Skewer *****/

double Skewer::finterp[Skewer::Ninterp];
bool Skewer::static_initialized = false;

/* This function needs to be fast.  To avoid branches, the preconditions
 * u2 >= 0 and u2 < 1 must be guaranteed by the pixel range calculation. */
double Skewer::f(double u2) {
    assert(0 <= u2 && u2 < 1);
    u2 *= (Ninterp-1);
    int j = int(u2);
    double x = u2 - j;
    return (1-x)*finterp[j] + x*finterp[j+1];
}

void Skewer::StaticInitialize() {
    /* Precompute SPH kernel */
    double u, u2;
    for(int j = 0; j < Ninterp; j++) {
        u2 = j/double(Ninterp-1);
        u = sqrt(u2);
        if(u < 0.5)
            finterp[j] = 8/M_PI * (1 - 6*u2*(1-u));
        else
            finterp[j] = 8/M_PI * 2*cb(1-u);
    }

    static_initialized = true;
}


Skewer::Skewer(const SkewerGeometry& g, int Npixels_, double h_) {
    if(!static_initialized)
        StaticInitialize();

    sx = g.sx;
    sy = g.sy;
    sz = g.sz;
    nx = g.nx;
    ny = g.ny;
    nz = g.nz;
    Ls = g.Ls;
    Npixels = Npixels_;
    h = h_;
    rho = (float*)malloc(Npixels*sizeof(float));
    rhov = (float*)malloc(Npixels*sizeof(float));
    memset(rho, 0, Npixels*sizeof(float));
    memset(rhov, 0, Npixels*sizeof(float));
}

Skewer::Skewer(double sx_, double sy_, double sz_, double nx_, double ny_, double nz_, double Ls_, int Npixels_, double h_) {
    if(!static_initialized)
        StaticInitialize();

    sx = sx_;
    sy = sy_;
    sz = sz_;
    nx = nx_;
    ny = ny_;
    nz = nz_;
    Ls = Ls_;
    Npixels = Npixels_;
    h = h_;
    rho = (float*)malloc(Npixels*sizeof(float));
    rhov = (float*)malloc(Npixels*sizeof(float));
    memset(rho, 0, Npixels*sizeof(float));
    memset(rhov, 0, Npixels*sizeof(float));
}

Skewer::~Skewer() {
    free(rho);
    free(rhov);
}

SkewerGeometry Skewer::GetGeometry() const {
    SkewerGeometry g;
    g.sx = sx;
    g.sy = sy;
    g.sz = sz;
    g.nx = nx;
    g.ny = ny;
    g.nz = nz;
    g.Ls = Ls;
    return g;
}

void Skewer::GetPixelPosition(int i, double& xi, double& yi, double& zi) const {
    double dL = Ls/Npixels;
    xi = sx + nx*(i+0.5)*dL;
    yi = sy + ny*(i+0.5)*dL;
    zi = sz + nz*(i+0.5)*dL;
}

void Skewer::Populate(POSVEL_T x, POSVEL_T y, POSVEL_T z, POSVEL_T vx, POSVEL_T vy, POSVEL_T vz) {
    const double dL = Ls/Npixels;

    /* Find projected distance of particle along skewer */
    double rx = x - sx;
    double ry = y - sy;
    double rz = z - sz;
    double rn = nx*rx + ny*ry + nz*rz;

    /* Find closest approach of particle to skewer */
    double h2 = h*h;
    double rp2 = sq(rx) + sq(ry) + sq(rz) - sq(rn);
    if(rp2 < h2) {
        /* If particle lies within distance h of skewer, populate nearby pixels */
        double I = rn/dL - 0.5;
        double dI = sqrt(h2 - rp2)/dL;
        int imin = (int)ceil(I - dI);
        int imax = (int)floor(I + dI);

        double xi, yi, zi;  // position of pixel i
        double r2, m;
        for(int i = imin; i <= imax; i++) {
            xi = sx + nx*(i+0.5)*dL;
            yi = sy + ny*(i+0.5)*dL;
            zi = sz + nz*(i+0.5)*dL;
            r2 = sq(x - xi) + sq(y - yi) + sq(z - zi);
            m = f(r2/h2);
            rho[i] += m;
            rhov[i] += m*(nx*vx + ny*vy + nz*vz);
        }
    }
}

void Skewer::Populate(int nparticles, POSVEL_T* x, POSVEL_T* y, POSVEL_T* z, POSVEL_T* vx, POSVEL_T* vy, POSVEL_T* vz) {
    const double dL = Ls/Npixels;
    double rx, ry, rz;  // particle position relative to skewer origin
    double rn;          // component of $\vec{r}$ along skewer direction
    double rp2;         // (distance of closest approach)^2 of particle to skewer
    double I;           // pixel coordinate of projected particle position
    double dI;
    int imin, imax;     // actual range of pixels to check
    int i;
    double xi, yi, zi;  // position of pixel i
    double r2, h2 = h*h, m;

    for(int n = 0; n < nparticles; n++) {
        /* Find projected distance of particle along skewer */
        rx = x[n] - sx;
        ry = y[n] - sy;
        rz = z[n] - sz;
        rn = nx*rx + ny*ry + nz*rz;

        /* Find closest approach of particle to skewer */
        rp2 = sq(rx) + sq(ry) + sq(rz) - sq(rn);
        if(rp2 < h2) {
            /* If particle lies within distance h of skewer, populate nearby pixels */
            I = rn/dL - 0.5;
            dI = sqrt(h2 - rp2)/dL;
            imin = (int)ceil(I - dI);
            imax = (int)floor(I + dI);
            for(i = imin; i <= imax; i++) {
                xi = sx + nx*(i+0.5)*dL;
                yi = sy + ny*(i+0.5)*dL;
                zi = sz + nz*(i+0.5)*dL;
                r2 = sq(x[n] - xi) + sq(y[n] - yi) + sq(z[n] - zi);
                m = f(r2/h2);
                rho[i] += m;
                rhov[i] += m*(nx*vx[n] + ny*vy[n] + nz*vz[n]);
            }
        }
    }
}

void Skewer::DumpPixels(FILE* fp) const {
    int n;
    assert(fp && !feof(fp) && !ferror(fp));
    n = fwrite(rho, sizeof(float), Npixels, fp);        assert(n == Npixels);
    n = fwrite(rhov, sizeof(float), Npixels, fp);       assert(n == Npixels);
}

void Skewer::ClearPixels() {
    memset(rho, 0, Npixels*sizeof(float));
    memset(rhov, 0, Npixels*sizeof(float));
}


/***** SkewerSet *****/

SkewerSet::SkewerSet() {
    Ntotal = 0;
    geoms = NULL;
    h = 0;
    Nlocal = 0;
    local_skewer_indices = NULL;
    local_skewers = NULL;
}

SkewerSet::~SkewerSet() {
    free(geoms);
    free(local_skewer_indices);
    if(Nlocal > 0 && local_skewers != NULL)
        for(int j = 0; j < Nlocal; j++)
            free(local_skewers[j]);
    free(local_skewers);
}

void SkewerSet::InitializeMaster(int Ntotal_, int Npixels_, double h_, int masterid) {
    int myid = Partition::getMyProc();
    if(myid == masterid) {
        Ntotal = Ntotal_;
        Npixels = Npixels_;
        h = h_;
        geoms = (SkewerGeometry*)malloc(Ntotal*sizeof(SkewerGeometry));
    }
}

void SkewerSet::InitializeSlaves(int masterid) {
    /* MPI broadcast skewer parameters from the master to slaves */
    MPI::Cartcomm comm = Partition::getComm();
    comm.Bcast(&Ntotal, 1, MPI_INT, masterid);
    comm.Bcast(&Npixels, 1, MPI_INT, masterid);
    comm.Bcast(&h, 1, MPI_DOUBLE, masterid);

    /* Allocate memory for skewer geometries on slaves */
    int myid = Partition::getMyProc();
    if(myid != masterid) {
        free(geoms);
        geoms = (SkewerGeometry*)malloc(Ntotal*sizeof(SkewerGeometry));
    }

    /* MPI broadcast geometries from master to slaves */
    assert(sizeof(SkewerGeometry) == 7*sizeof(double));
    comm.Bcast((double*)geoms, 7*Ntotal, MPI_DOUBLE, masterid);

    /* Get a (padded) bounding box for this processor's domain */
    float dpos[3], dlen[3];       // position and length of domain in global coordinates
    Domain::corner_phys_alive(dpos);
    Domain::rL_local_alive(dlen);
    Box mybox;
    mybox.xmin = dpos[0] - h;
    mybox.xmax = dpos[0] + dlen[0] + h;
    mybox.ymin = dpos[1] - h;
    mybox.ymax = dpos[1] + dlen[1] + h;
    mybox.zmin = dpos[2] - h;
    mybox.zmax = dpos[2] + dlen[2] + h;
      
    /* Create a list of skewers that intersect this processor's domain */
    free(local_skewer_indices);
    free(local_skewers);
    local_skewer_indices = (int*)malloc(Ntotal*sizeof(int));
    local_skewers = (Skewer**)malloc(Ntotal*sizeof(Skewer*));
    Nlocal = 0;
    for(int i = 0; i < Ntotal; i++) {
        if(geoms[i].Intersects(mybox)) {
            local_skewer_indices[Nlocal] = i;
            local_skewers[Nlocal] = new Skewer(geoms[i], Npixels, h);
            Nlocal++;
        }
    }
}

void SkewerSet::RandomSkewers1(double L) {
    assert(geoms != NULL);

    for(int i = 0; i < Ntotal; i++) {
        geoms[i].sx = L*drand48();
        geoms[i].sy = L*drand48();
        geoms[i].sz = -h;
        geoms[i].nx = 0;
        geoms[i].ny = 0;
        geoms[i].nz = 1;
        geoms[i].Ls = L+2*h;
    }
}

void SkewerSet::RandomSkewers2(double L, double zbox) {
    assert(geoms != NULL);

    double R, phi;
    double rx, ry, rz, absr;
    double tstart, tend;
    for(int i = 0; i < Ntotal; i++) {
        /* First pick a point within the circle of radius L/2 on the far face of the box */
        R = (0.99*L/2) * drand48();
        phi = (2*M_PI) * drand48();

        /* Compute the normal vector from the origin to this point */
        rx = R*cos(phi);
        ry = R*sin(phi);
        rz = zbox + L/2;
        absr = sqrt(rx*rx + ry*ry + rz*rz);
        geoms[i].nx = rx/absr;
        geoms[i].ny = ry/absr;
        geoms[i].nz = rz/absr;

        /* Then compute the starting point of the skewer, and its length */
        tstart = (zbox - L/2 - h)/geoms[i].nz;
        tend   = (zbox + L/2 + h)/geoms[i].nz;
        geoms[i].sx = tstart * geoms[i].nx + L/2;
        geoms[i].sy = tstart * geoms[i].ny + L/2;
        geoms[i].sz = tstart * geoms[i].nz + L/2 - zbox;
        geoms[i].Ls = tend - tstart;
    }
}

void SkewerSet::PopulateLocal(int Np, POSVEL_T* x, POSVEL_T* y, POSVEL_T* z, POSVEL_T* vx, POSVEL_T* vy, POSVEL_T* vz) {
    for(int j = 0; j < Nlocal; j++)
        local_skewers[j]->Populate(Np, x, y, z, vx, vy, vz);
}

void SkewerSet::WriteLocalSkewers(const char* filename) {
    FILE* fp = fopen(filename, "wb");
    if(fp == NULL) {
        fprintf(stderr, "Error: could not write local skewers for processor %d to %s [%s]\n", Partition::getMyProc(), filename, strerror(errno));
        return;
    }

    int i, n;
    n = fwrite(&Nlocal, sizeof(int), 1, fp);    assert(n == 1);
    n = fwrite(&Npixels, sizeof(int), 1, fp);   assert(n == 1);
    n = fwrite(&Ntotal, sizeof(int), 1, fp);    assert(n == 1);
    for(int j = 0; j < Nlocal; j++) {
        i = local_skewer_indices[j];
        n = fwrite(&i, sizeof(int), 1, fp);     assert(n == 1);
        n = fwrite(&geoms[i], sizeof(SkewerGeometry), 1, fp);   assert(n == 1);
        local_skewers[j]->DumpPixels(fp);
    }

    fclose(fp);
}

void SkewerSet::ClearSkewers() {
    for(int j = 0; j < Nlocal; j++)
        local_skewers[j]->ClearPixels();
}


ParallelSkewerSet::ParallelSkewerSet() {
    L = 0;
}

ParallelSkewerSet::~ParallelSkewerSet() {
}

void ParallelSkewerSet::InitializeMaster(int Ntotal, int Npixels, double h, double L_, int Nmesh_, int masterid) {
    SkewerSet::InitializeMaster(Ntotal, Npixels, h);

    int myid = Partition::getMyProc();
    if(myid == masterid) {
        L = L_;
        Nmesh = Nmesh_;

        /* Make sure the chaining mesh cells are at least as big as h (so that only
         * need to check particles against skewers in immediately adjacent mesh
         * cells) */
        assert(L/Nmesh >= h);

        /* Generate random skewers */
        SkewerSet::RandomSkewers1(L);
    }
}

void ParallelSkewerSet::InitializeSlaves(int masterid) {
    SkewerSet::InitializeSlaves(masterid);

    /* MPI broadcast chaining mesh parameters from the master to slaves */
    MPI::Cartcomm comm = Partition::getComm();
    comm.Bcast(&L, 1, MPI_DOUBLE, masterid);
    comm.Bcast(&Nmesh, 1, MPI_INT, masterid);

    /* Fill chaining mesh with local skewers (on all processors) */
    mesh.resize(Nmesh*Nmesh);
    int ix, iy; // indices for grid cell within chaining mesh
    for(int j = 0; j < Nlocal; j++) {
        int i = local_skewer_indices[j];
        ix = (geoms[i].sx)/L * Nmesh;
        iy = (geoms[i].sy)/L * Nmesh;
        assert(0 <= ix && ix < Nmesh && 0 <= iy && iy < Nmesh);
        mesh[Nmesh*iy + ix].push_back(local_skewers[j]);
    }
}

void ParallelSkewerSet::PopulateLocal(int Np, POSVEL_T* x, POSVEL_T* y, POSVEL_T* z, POSVEL_T* vx, POSVEL_T* vy, POSVEL_T* vz) {
    int ix, iy, ixmin, ixmax, iymin, iymax;
    list<Skewer*>::iterator iter;
    for(int n = 0; n < Np; n++) {
        ix = (x[n]/L) * Nmesh;
        iy = (y[n]/L) * Nmesh;
        ixmin = max(ix-1, 0);
        ixmax = min(ix+1, Nmesh-1);
        iymin = max(iy-1, 0);
        iymax = min(iy+1, Nmesh-1);
        assert(ixmin < ixmax && iymin < iymax);
        for(ix = ixmin; ix <= ixmax; ix++) {
            for(iy = iymin; iy <= iymax; iy++) {
                list<Skewer*>& skewers = mesh[Nmesh*iy + ix];
                for(iter = skewers.begin(); iter != skewers.end(); iter++)
                    (*iter)->Populate(x[n], y[n], z[n], vx[n], vy[n], vz[n]);
            }
        }
    }
}

void ParallelSkewerSet::PopulateLocal(int Np, unsigned int* indices, POSVEL_T* x, POSVEL_T* y, POSVEL_T* z, POSVEL_T* vx, POSVEL_T* vy, POSVEL_T* vz) {
    int n;
    int ix, iy, ixmin, ixmax, iymin, iymax;
    list<Skewer*>::iterator iter;
    for(int m = 0; m < Np; m++) {
        n = indices[m];
        ix = (x[n]/L) * Nmesh;
        iy = (y[n]/L) * Nmesh;
        ixmin = max(ix-1, 0);
        ixmax = min(ix+1, Nmesh-1);
        iymin = max(iy-1, 0);
        iymax = min(iy+1, Nmesh-1);
        assert(ixmin < ixmax && iymin < iymax);
        for(ix = ixmin; ix <= ixmax; ix++) {
            for(iy = iymin; iy <= iymax; iy++) {
                list<Skewer*>& skewers = mesh[Nmesh*iy + ix];
                for(iter = skewers.begin(); iter != skewers.end(); iter++)
                    (*iter)->Populate(x[n], y[n], z[n], vx[n], vy[n], vz[n]);
            }
        }
    }
}


ConvergentSkewerSet::ConvergentSkewerSet() {
    L = 0;
    zbox = 0;
}

ConvergentSkewerSet::~ConvergentSkewerSet() {
}

void ConvergentSkewerSet::InitializeMaster(int Ntotal, int Npixels, double h, double L_, double zbox_, int Nmesh_, int masterid) {
    SkewerSet::InitializeMaster(Ntotal, Npixels, h);

    int myid = Partition::getMyProc();
    if(myid == masterid) {
        L = L_;
        zbox = zbox_;
        Nmesh = Nmesh_;

        /* Make sure the chaining mesh cells are at least as big as h (so that only
         * need to check particles against skewers in immediately adjacent mesh
         * cells) */
        assert(L/Nmesh >= h);

        /* Generate random skewers */
        SkewerSet::RandomSkewers2(L, zbox);
    }
}

void ConvergentSkewerSet::InitializeSlaves(int masterid) {
    SkewerSet::InitializeSlaves(masterid);

    /* MPI broadcast chaining mesh parameters from the master to slaves */
    MPI::Cartcomm comm = Partition::getComm();
    comm.Bcast(&L, 1, MPI_DOUBLE, masterid);
    comm.Bcast(&zbox, 1, MPI_DOUBLE, masterid);
    comm.Bcast(&Nmesh, 1, MPI_INT, masterid);

    /* Fill chaining mesh with local skewers (on all processors) */
    mesh.resize(Nmesh*Nmesh);
    double x, y;        // where skewer intersects z=0 (lower box face)
    double t;           // distance along skewer where this happens
    int ix, iy;         // indices for grid cell within chaining mesh
    int i;
    for(int j = 0; j < Nlocal; j++) {
        i = local_skewer_indices[j];
        t = geoms[i].sz / geoms[i].nz;
        x = geoms[i].sx - t*geoms[i].nx;
        y = geoms[i].sy - t*geoms[i].ny;
        ix = x/L * Nmesh;
        iy = y/L * Nmesh;
        assert(0 <= ix && ix < Nmesh && 0 <= iy && iy < Nmesh);
        mesh[Nmesh*iy + ix].push_back(local_skewers[j]);
    }
}

void ConvergentSkewerSet::PopulateLocal(int Np, POSVEL_T* x, POSVEL_T* y, POSVEL_T* z, POSVEL_T* vx, POSVEL_T* vy, POSVEL_T* vz) {
    int ix, iy, ixmin, ixmax, iymin, iymax;
    double xp, yp;      // position of particle projected onto z=0 plane toward convergence point
    double t;
    list<Skewer*>::iterator iter;
    for(int n = 0; n < Np; n++) {
        t = z[n] / (z[n] + zbox - L/2);
        xp = t*(L/2) + (1-t)*x[n];
        yp = t*(L/2) + (1-t)*y[n];
        ix = (xp/L) * Nmesh;
        iy = (yp/L) * Nmesh;
        ixmin = max(ix-1, 0);
        ixmax = min(ix+1, Nmesh-1);
        iymin = max(iy-1, 0);
        iymax = min(iy+1, Nmesh-1);
        assert(ixmin < ixmax && iymin < iymax);
        for(ix = ixmin; ix <= ixmax; ix++) {
            for(iy = iymin; iy <= iymax; iy++) {
                list<Skewer*>& skewers = mesh[Nmesh*iy + ix];
                for(iter = skewers.begin(); iter != skewers.end(); iter++)
                    (*iter)->Populate(x[n], y[n], z[n], vx[n], vy[n], vz[n]);
            }
        }
    }
}

void ConvergentSkewerSet::PopulateLocal(int Np, unsigned int* indices, POSVEL_T* x, POSVEL_T* y, POSVEL_T* z, POSVEL_T* vx, POSVEL_T* vy, POSVEL_T* vz) {
    int n;
    double xp, yp;      // position of particle projected onto z=0 plane toward convergence point
    double t;
    int ix, iy, ixmin, ixmax, iymin, iymax;
    list<Skewer*>::iterator iter;
    for(int m = 0; m < Np; m++) {
        n = indices[m];
        t = z[n] / (z[n] + zbox - L/2);
        xp = t*(L/2) + (1-t)*x[n];
        yp = t*(L/2) + (1-t)*y[n];
        ix = (xp/L) * Nmesh;
        iy = (yp/L) * Nmesh;
        ix = (x[n]/L) * Nmesh;
        iy = (y[n]/L) * Nmesh;
        ixmin = max(ix-1, 0);
        ixmax = min(ix+1, Nmesh-1);
        iymin = max(iy-1, 0);
        iymax = min(iy+1, Nmesh-1);
        assert(ixmin <= ixmax && iymin <= iymax);
        for(ix = ixmin; ix <= ixmax; ix++) {
            for(iy = iymin; iy <= iymax; iy++) {
                list<Skewer*>& skewers = mesh[Nmesh*iy + ix];
                for(iter = skewers.begin(); iter != skewers.end(); iter++)
                    (*iter)->Populate(x[n], y[n], z[n], vx[n], vy[n], vz[n]);
            }
        }
    }
}
