#ifndef SPLINE_H
#define SPLINE_H

#include <vector>
using std::vector;

class Spline;
class SplineImpl;


/***** Spline factory functions *****/

/* Linear interpolation */
Spline LinearSpline(const vector<double>& X, const vector<double>& Y);
Spline LinearSpline(int N, const double* X, const double* Y);

/* Shifted linear interpolation (see Blu, Thevenaz, and Unser, 2004) */
Spline ShiftedLinearSpline(const vector<double>& X, const vector<double>& Y, double tau = 0.2);
Spline ShiftedLinearSpline(int N, const double* X, const double* Y, double tau = 0.2);

/* Natural cubic spline */
Spline CubicSpline(const vector<double>& X, const vector<double>& Y);
Spline CubicSpline(int N, const double* X, const double* Y);


/***************************************************************
 * Spline
 *
 * Generic spline wrapper class.
 ***************************************************************/
class Spline {
public:
    /* Default to cubic spline */
    Spline(const vector<double>& X, const vector<double>& Y);
    Spline(int N, const double* X, const double* Y);

    Spline();
    Spline(SplineImpl* impl, bool clone = false);
    Spline(const Spline& F);
    ~Spline();

    Spline& operator=(const Spline& F);

    double Evaluate(double x) const;
    double EvaluateDerivative(double x) const;

    double operator()(double x) const { return Evaluate(x); }

    /* Find a local maximum or minimum of the interpolated function */
//    double FindMaximum(double xguess, double* ymax = 0);
//    double FindMinimum(double xguess, double& ymin = 0);

protected:
    SplineImpl* impl;   // internal spline implementation
};


/**************************************************
 * SplineImpl
 *
 * Base class for internal spline implementations.
 **************************************************/
struct SplineImpl {
    SplineImpl() {}
    virtual ~SplineImpl() {}

    virtual double y(double x) const = 0;
    virtual double dydx(double x) const = 0;
    virtual SplineImpl* clone() const = 0;

    double xmin, xmax;    // domain
};

#endif // SPLINE_H
