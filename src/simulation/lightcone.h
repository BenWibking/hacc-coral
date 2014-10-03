#ifndef LIGHTCONE_H
#define LIGHTCONE_H

#include "spline.h"

class LightCone {
public:
    /* Construct a light cone around the point (x0,y0,z0) */
    LightCone(double OmegaM, double x0 = 0, double y0 = 0, double z0 = 0, double zmax = 10);

    /* Move the origin of the light cone to the point (x0,y0,z0). */
    void SetOrigin(double x0, double y0, double z0);

    /* Focus attention on a light cone shell (afar < anear). */
    void DefineShell(double afar, double anear);

    /* Returns true if the point (x,y,z) lies within the current light cone shell. */
    bool LiesInShell(double x, double y, double z) const;

    void Debug() const;

    Spline rofa;
    Spline aofr;

protected:
    double x0, y0, z0;
    double amin;
    double afar, anear;
    double r2far, r2near;

    static double cOverH(double a, double OmegaM);
};

#endif // LIGHTCONE_H
