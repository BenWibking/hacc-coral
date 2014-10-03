#include <cassert>
#include <cmath>
#include <cstdlib>
#include <cstring>

#include "spline.h"

/* Preconditions: X[0] <= x <= X[N-1], X[i] < X[i+1] for 0 <= i <= N-2
 * Postcondition: 0 <= i <= N-2 */
static int LookupIndex(double x, int N, const double* X) {
    assert(!isinf(x) && !isnan(x) && X[0] <= x && x <= X[N-1]);

    /* Potentially pathological case */
    if(x == X[N-1])
        return N-2;

    /* First calculate what the index would be for uniformly spaced points */
    double I = (N-1) * (x - X[0])/(X[N-1] - X[0]);
    int i = (int)I;

    /* Check this index and its immediate neighbors, in case the points are
     * actually uniform (or close to it) */
    if(X[i] <= x && x <= X[i+1])
        return i;
    else if(X[i-1] <= x && x < X[i])
        return i-1;
    else if(X[i+1] < x && x <= X[i+2])
        return i+1;
    else {
        /* Finally resort to a binary search starting from the initial guess */
        int ilo = 0;
        int ihi = N-2;
        while(true) {
            if(x < X[i]) {
                ihi = i-1;
                i = (ilo+ihi)/2;
            }
            else if(x > X[i+1]) {
                ilo = i+1;
                i = (ilo+ihi)/2;
            }
            else
                break;
        }
        return i;
    }
}

static double* new_array(int n) {
    return (double*)malloc(n*sizeof(double));
}

static double* new_array(int n, const double* src) {
    double* dest = (double*)malloc(n*sizeof(double));
    memcpy(dest, src, n*sizeof(double));
    return dest;
}


/*****************
 * Linear spline *
 *****************/

struct LinearSplineImpl : public SplineImpl {
    int N;
    double* X;
    double* Y;

    LinearSplineImpl() {
    }

    LinearSplineImpl(int nn, const double* xx, const double* yy) {
        N = nn;
        X = new_array(N, xx);
        Y = new_array(N, yy);

        xmin = X[0];
        xmax = X[N-1];
    }

    virtual ~LinearSplineImpl() {
        free(X);
        free(Y);
    }

    double y(double x) const {
        int i = LookupIndex(x, N, X);
        assert(0 <= i && i <= N-2);
        double x1 = X[i], x2 = X[i+1], y1 = Y[i], y2 = Y[i+1];
        return y1 + (y2 - y1)*(x - x1)/(x2 - x1);
    }

    double dydx(double x) const {
        int i = LookupIndex(x, N, X);
        assert(0 <= i && i <= N-2);
        double x1 = X[i], x2 = X[i+1], y1 = Y[i], y2 = Y[i+1];
        return (y2 - y1)/(x2 - x1);
    }

    double max(double& xmax) const {
        double ymax = -1e100;
        for(int i = 0; i < N; i++) {
            if(Y[i] > ymax) {
                xmax = X[i];
                ymax = Y[i];
            }
        }
        return ymax;
    }

    double min(double& xmin) const {
        double ymin = 1e100;
        for(int i = 0; i < N; i++) {
            if(Y[i] < ymin) {
                xmin = X[i];
                ymin = Y[i];
            }
        }
        return ymin;
    }

    LinearSplineImpl* clone() const {
        return new LinearSplineImpl(N, X, Y);
    }
};

Spline LinearSpline(const vector<double>& X, const vector<double>& Y) {
    assert(X.size() == Y.size());
    return Spline(new LinearSplineImpl(X.size(), &X[0], &Y[0]));
}

Spline LinearSpline(int N, const double* X, const double* Y) {
    return Spline(new LinearSplineImpl(N, X, Y));
}


/**********************************************
 * Shifted linear spline
 * 
 * Uses constant extrapolation at boundaries.
 *********************************************/

struct ShiftedLinearSplineImpl : public LinearSplineImpl {
    double tau;

    ShiftedLinearSplineImpl(int nn, const double* xx, const double* yy, double tau) {
        xmin = xx[0];
        xmax = xx[nn-1];
        double T = xx[1] - xx[0];

#ifdef DEBUG
        /* (Shifted linear interpolation only works for uniformly spaced points.) */
        for(int i = 0; i < nn-1; i++)
            assert(fabs(xx[i+1] - xx[i] - T) < 1e-10);
#endif

        N = nn+1;
        X = new_array(N);
        Y = new_array(N);
        X[0] = xx[0] - (1-tau)*T;
        Y[0] = yy[0];
        for(int i = 0; i < nn; i++) {
            X[i+1] = xx[i] + tau*T;
            Y[i+1] = -tau/(1 - tau) * Y[i] + 1/(1-tau)*yy[i];
        }
    }
};

Spline ShiftedLinearSpline(const vector<double>& X, const vector<double>& Y, double tau) {
    assert(X.size() == Y.size());
    return Spline(new ShiftedLinearSplineImpl(X.size(), &X[0], &Y[0], tau));
}

Spline ShiftedLinearSpline(int N, const double* X, const double* Y, double tau) {
    return Spline(new ShiftedLinearSplineImpl(N, X, Y, tau));
}


/************************************************************
 * Cubic spline
 * 
 * Implementation copied from netlib's fmm/spline.f program.
 ************************************************************/

struct CubicSplineImpl : public SplineImpl {
    /* y(x) = y[i] + b (x-x[i]) + c (x-x[i])^2 + d (x-x[i])^3  for  x[i] <= x < x[i+1] */
    int N;
    double* X;
    double* Y;
    double* b;
    double* c;
    double* d;

    CubicSplineImpl(int nn, const double* xx, const double* yy) {
        N = nn;
        assert(N >= 2);
        X = new_array(N, xx);
        Y = new_array(N, yy);
        b = new_array(N);
        c = new_array(N);
        d = new_array(N);

        xmin = X[0];
        xmax = X[N-1];

        if(N == 2) {
            b[0] = b[1] = (Y[1] - Y[0])/(X[1] - X[0]);
            c[0] = c[1] = 0;
            d[0] = d[1] = 0;
        }
        else {
            /* Set up tridiagonal system: b = diagonal, d = offdiagonal, c = right hand side */
            d[0] = X[1] - X[0];
            c[1] = (Y[1] - Y[0])/d[0];
            for(int i = 1; i < N-1; i++) {
                d[i] = X[i+1] - X[i];
                b[i] = 2*(d[i-1] + d[i]);
                c[i+1] = (Y[i+1] - Y[i])/d[i];
                c[i] = c[i+1] - c[i];
            }

            /* End conditions: third derivatives at X[0] and X[N-1] obtained from divide differences */
            b[0] = -d[0];
            b[N-1] = -d[N-2];
            c[0] = 0;
            c[N-1] = 0;
            if(N > 3) {
                c[0] = c[2]/(X[3] - X[1]) - c[1]/(X[2] - X[0]);
                c[N-1] = c[N-2]/(X[N-1] - X[N-3]) - c[N-3]/(X[N-2] - X[N-4]);
                c[0] = c[0]*d[0]*d[0]/(X[3] - X[0]);
                c[N-1] = -c[N-1]*d[N-1]*d[N-1]/(X[N-1] - X[N-4]);
            }

            /* Forward elimination */
            for(int i = 1; i < N; i++) {
                double t = d[i-1]/b[i-1];
                b[i] -= t*d[i-1];
                c[i] -= t*c[i-1];
            }

            /* Back substitution */
            c[N-1] = c[N-1]/b[N-1];
            for(int j = 0; j < N-1; j++) {
                int i = N - j - 2;
                c[i] = (c[i] - d[i]*c[i+1])/b[i];
            }

            /* Compute polynomial coefficients */
            b[N-1] = (Y[N-1] - Y[N-2])/d[N-2] + d[N-2]*(c[N-2] + 2*c[N-1]);
            for(int i = 0; i < N-1; i++) {
                b[i] = (Y[i+1] - Y[i])/d[i] - d[i]*(c[i+1] + 2*c[i]);
                d[i] = (c[i+1] - c[i])/d[i];
                c[i] = 3*c[i];
            }
            c[N-1] = 3*c[N-1];
            d[N-1] = d[N-2];
        }
    }

    virtual ~CubicSplineImpl() {
        free(X);
        free(Y);
        free(b);
        free(c);
        free(d);
    }

    double y(double x) const {
        int i = LookupIndex(x, N, X);
        assert(0 <= i && i <= N-2);
        double u = x - X[i];
        return Y[i] + b[i]*u + c[i]*u*u + d[i]*u*u*u;
    }

    double dydx(double x) const {
        int i = LookupIndex(x, N, X);
        assert(0 <= i && i <= N-2);
        double u = x - X[i];
        return b[i] + 2*c[i]*u + 3*d[i]*u*u;
    }

//    double max(int i, double& xmax) const {
//        return 0;
//    }

    CubicSplineImpl* clone() const {
        return new CubicSplineImpl(N, X, Y);
    }
};

Spline CubicSpline(const vector<double>& X, const vector<double>& Y) {
    assert(X.size() == Y.size());
    return Spline(new CubicSplineImpl(X.size(), &X[0], &Y[0]));
}

Spline CubicSpline(int N, const double* X, const double* Y) {
    return Spline(new CubicSplineImpl(N, X, Y));
}


/**********
 * Spline *
 **********/

Spline::Spline(const vector<double>& X, const vector<double>& Y) {
    assert(X.size() == Y.size());
    impl = new LinearSplineImpl(X.size(), &X[0], &Y[0]);
}

Spline::Spline(int n, const double* X, const double* Y) {
    impl = new LinearSplineImpl(n, X, Y);
}

Spline::Spline() {
    impl = NULL;
}

Spline::Spline(SplineImpl* otherimpl, bool clone) {
    if(clone && otherimpl != NULL)
        impl = otherimpl->clone();
    else
        impl = otherimpl;
}

Spline::Spline(const Spline& S) {
    if(S.impl)
        impl = S.impl->clone();
    else
        impl = NULL;
}

Spline::~Spline() {
    delete impl;
}

Spline& Spline::operator=(const Spline& S) {
    delete impl;
    impl = (S.impl) ? S.impl->clone() : NULL;
    return *this;
}

double Spline::Evaluate(double x) const {
    if(!impl || x < impl->xmin || x > impl->xmax)
        return 0;
    else
        return impl->y(x);
}

double Spline::EvaluateDerivative(double x) const {
    if(!impl || x < impl->xmin || x > impl->xmax)
        return 0;
    else
        return impl->dydx(x);
}

#if 0
double Spline::FindMaximum(double xguess, double* yret) {
    assert(N >= 0);

    double xmax, ymax;
    if(xguess < X[0]) {
        ymax = yleft;
        xmax = xguess;
    }
    else if(xguess > X[N]) {
        ymax = yright;
        xmax = xguess;
    }
    else {
        /* Now we walk from one segment to the next until we find a maximum */

        int i = LookupIndex(xguess);
        /* If we're right on a boundary the algorithm breaks */
        if(xguess == X[i])
            xguess += (X[i+1]-X[i])/1000;
        else if(xguess == X[i+1])
            xguess -= (X[i+1]-X[i])/1000;

        xmax = GetSegment(i).max(ymax);
        if(xmax == xguess)
            xmax += (X[i+1]-X[i])/100;
        int dir = (xmax > xguess) ? +1 : -1;

        while(true) {
            if(X[i] < xmax && xmax < X[i+1])    // local maximum
                break;
            else if(xmax == X[i]) {     // maximum is the left-most point of the interval
                if(dir == +1)           // if we started off walking right, we're done
                    break;
                if(i == 0)              // if we can't go any further left, we're done
                    break;
                i--;
            }
            else if(xmax == X[i+1]) {   // maximum is the right-most point of the interval
                if(dir == -1)           // if we started off walking left, we're done
                    break;
                if(i == N-1)            // if we can't go any further right, we're done
                    break;
                i++;
            }
            else {
                fprintf(stderr, "InterpolatingFunction::FindMaximum: broken\n");
            }
            xmax = GetSegment(i).max(ymax);
        }
    }

    if(yret != NULL)
        *yret = ymax;
    return xmax;
}
#endif
