#ifndef SKEWERS_H
#define SKEWERS_H

#include <list>
#include <vector>

#include "Definition.h"

struct Box {
    double xmin, xmax;
    double ymin, ymax;
    double zmin, zmax;
};

struct SkewerGeometry {
    double sx, sy, sz;  // initial point
    double nx, ny, nz;  // direction
    double Ls;          // length

    /* Return true if this skewer intersects the box b. */
    bool Intersects(const Box& b) const;
};


class Skewer {
public:
    /* Create a skewer with a given skewer geometry. */
    Skewer(const SkewerGeometry& g, int Npixels, double h);
    Skewer(double sx, double sy, double sz, double nx, double ny, double nz, double Ls, int Npixels, double h);

    ~Skewer();

    /* Return the geometry for this skewer. */
    SkewerGeometry GetGeometry() const;

    /* Return the smoothing scale for this skewer. */
    double GetSmoothingScale() { return h; }

    /* Get the location of the ith pixel along this skewer. */
    void GetPixelPosition(int i, double& xi, double& yi, double& zi) const;

    /* Populate this skewer with particles. */
    void Populate(POSVEL_T x, POSVEL_T y, POSVEL_T z, POSVEL_T vx, POSVEL_T vy, POSVEL_T vz);
    void Populate(int nparticles, POSVEL_T* x, POSVEL_T* y, POSVEL_T* z, POSVEL_T* vx, POSVEL_T* vy, POSVEL_T* vz);

    /* Dump raw pixels to file, as a sequence of floats. */
    void DumpPixels(FILE* fp) const;

    /* Reset skewer to zero density */
    void ClearPixels();

protected:
    double sx, sy, sz;
    double nx, ny, nz;
    double Ls;
    int Npixels;
    double h;
    float* rho;
    float* rhov;

    static const int Ninterp = 512;     // number of points for interpolated kernel
    static double finterp[Ninterp];     // data for interpolated kernel
    static double f(double u2);         // SPH kernel (as a function of u2 = r^2/h^2)
    static bool static_initialized;

    /* Pre-calculate universal smoothing kernel. */
    static void StaticInitialize();

    friend class SkewerSet;
};


class SkewerSet {
public:
    SkewerSet();
    virtual ~SkewerSet();

    /* Initialize skewer geometries on master node. */
    virtual void InitializeMaster(int Ntotal, int Npixels, double h, int masterid = 0);

    /* Generate random skewers through the box [0,L]^3 along the z-axis. */
    virtual void RandomSkewers1(double L);

    /* Generate random skewers through the box [0,L]^3 as if the box were
     * centered a distance zbox above the origin, with skewers starting at the
     * origin and ending on the far face of the box. */
    virtual void RandomSkewers2(double L, double zbox);

    /* Broadcast geometries from the master to slaves, and allocate memory for
     * local skewers. */
    virtual void InitializeSlaves(int masterid = 0);

    /* Populate all local skewers with particles */
    virtual void PopulateLocal(int nparticles, POSVEL_T* x, POSVEL_T* y, POSVEL_T* z, POSVEL_T* vx, POSVEL_T* vy, POSVEL_T* vz);

    /* Write local skewers to file, in temporary skewer format.
     * [A temporary skewer file (extension .skt) is a binary file consisting of
     * the following:
     *  (4-byte integer) Nlocal -- number of local skewers described in this file
     *  (4-byte integer) Npixels -- number of pixels per skewer
     *  (4-byte integer) Ntotal -- total number of skewers in this collection
     *  for j = 1 to Nlocal:
     *          (4 byte integer) skewer ID within collection
     *          (28 bytes) SkewerGeometry struct for skewer j
     *          (Npixels 4-byte floats) rho along skewer
     *          (Npixels 4-byte floats) rhov along skewer
     * In total the file size should be 12 + Nlocal*(32 + 8*Npixels).] */
    virtual void WriteLocalSkewers(const char* filename);

    /* Reset all skewers to zero density */
    virtual void ClearSkewers();

protected:
    int Ntotal;                 // total number of skewers
    int Npixels;                // number of pixels per skewer
    double h;                   // smoothing scale
    int Nlocal;                 // number of skewers that intersect this node's domain
    SkewerGeometry* geoms;      // geometries for all skewers
    int* local_skewer_indices;  // indices of skewers that intersect this node's domain
    Skewer** local_skewers;     // full memory allocation for local skewers
};

class ParallelSkewerSet : public SkewerSet {
public:
    ParallelSkewerSet();
    virtual ~ParallelSkewerSet();

    /* Generate Ntotal random skewer geometries through the box [0,L]^3
     * parallel to the z-axis. */
    virtual void InitializeMaster(int Ntotal, int Npixels, double h, double L, int Nmesh = 512, int masterid = 0);
//    virtual void InitializeMaster(const Basedata& indat);

    /* Broadcast geometries from the master to slaves.  Allocate memory for
     * local skewers, and build chaining mesh. */
    virtual void InitializeSlaves(int masterid = 0);

    /* Populate all local skewers with particles */
    virtual void PopulateLocal(int nparticles, POSVEL_T* x, POSVEL_T* y, POSVEL_T* z, POSVEL_T* vx, POSVEL_T* vy, POSVEL_T* vz);
    virtual void PopulateLocal(int nparticles, unsigned int* indices, POSVEL_T* x, POSVEL_T* y, POSVEL_T* z, POSVEL_T* vx, POSVEL_T* vy, POSVEL_T* vz);

protected:
    double L;           // box length
    int Nmesh;          // number of grid cells in chaining mesh
    std::vector<std::list<Skewer*> > mesh;
};

class ConvergentSkewerSet : public SkewerSet {
public:
    ConvergentSkewerSet();
    virtual ~ConvergentSkewerSet();

    /* Generate random skewers through the box [0,L]^3 as if the box were
     * centered a distance zbox above the origin, with skewers starting at the
     * origin and ending on the far face of the box. */
    virtual void InitializeMaster(int Ntotal, int Npixels, double h, double L, double zbox, int Nmesh = 512, int masterid = 0);
//    virtual void InitializeMaster(const Basedata& indat);

    /* Broadcast geometries from the master to slaves.  Allocate memory for
     * local skewers, and build chaining mesh. */
    virtual void InitializeSlaves(int masterid = 0);

    /* Populate all local skewers with particles. */
    virtual void PopulateLocal(int nparticles, POSVEL_T* x, POSVEL_T* y, POSVEL_T* z, POSVEL_T* vx, POSVEL_T* vy, POSVEL_T* vz);
    virtual void PopulateLocal(int nparticles, unsigned int* indices, POSVEL_T* x, POSVEL_T* y, POSVEL_T* z, POSVEL_T* vx, POSVEL_T* vy, POSVEL_T* vz);

protected:
    double L;           // box length
    double zbox;        // distance from convergence point to box center
    int Nmesh;          // number of grid cells in chaining mesh
    std::vector<std::list<Skewer*> > mesh;
};

#endif // SKEWERS_H
