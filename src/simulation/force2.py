import numpy as N
import numpy.random as NR
import sys

class ShortRangeForce:
    def __init__(self,
                 rsm=0.01,
                 b=0.72, c=0.01, d=0.27, e=0.0001, 
                 f=360.0, g=100.0, h=0.67, l=17.0,
                 zerocross=3.1163263553866876):
        self.rsm = rsm
        self.b = b
        self.c = c
        self.d = d
        self.e = e
        self.f = f
        self.g = g
        self.h = h
        self.l = l
        self.zerocross = zerocross

    def f_grid(self,r):
        if N.isscalar(r):
            if r > 0:
                return N.tanh(self.b*r)/r**2 \
                    - self.b/r/(N.cosh(self.b*r))**2 \
                    + self.c*r*(1.0+self.d*r**2)*N.exp(-1.0*self.d*r**2) \
                    + self.e*(self.f*r**2+self.g*r**4+self.l*r**6)*N.exp(-1.0*self.h*r**2)
            else:
                return 0.0
        else:
            ret = 0.0*r
            indx = N.where(r>0.0)
            ret[indx] = N.tanh(self.b*r[indx])/r[indx]**2 \
                    - self.b/r[indx]/(N.cosh(self.b*r[indx]))**2 \
                    + self.c*r[indx]*(1.0+self.d*r[indx]**2)*N.exp(-1.0*self.d*r[indx]**2) \
                    + self.e*(self.f*r[indx]**2+self.g*r[indx]**4+self.l*r[indx]**6)*N.exp(-1.0*self.h*r[indx]**2)
            return ret

    def f_gridOR(self,r):
        if N.isscalar(r):
            if r > 0:
                return N.tanh(self.b*r)/r**3 \
                    - self.b/r**2/(N.cosh(self.b*r))**2 \
                    + self.c*(1.0+self.d*r**2)*N.exp(-1.0*self.d*r**2) \
                    + self.e/r*(self.f*r**2+self.g*r**4+self.l*r**6)*N.exp(-1.0*self.h*r**2)
            else:
                return self.c + 2.0/3.0*self.b**3
        else:
            ret = 0.0*r + self.c + 2.0/3.0*self.b**3
            indx = N.where(r>0.0)
            ret[indx] = N.tanh(self.b*r[indx])/r[indx]**3 \
                    - self.b/r[indx]**2/(N.cosh(self.b*r[indx]))**2 \
                    + self.c*(1.0+self.d*r[indx]**2)*N.exp(-1.0*self.d*r[indx]**2) \
                    + self.e/r[indx]*(self.f*r[indx]**2+self.g*r[indx]**4+self.l*r[indx]**6)*N.exp(-1.0*self.h*r[indx]**2)
            return ret

    def f_sr(self,r):
        return 1.0/r**2-self.f_grid(r)

    def one_over_r2(self,r):
        return 1.0/r**2

    def one_over_r3(self,r):
        return 1.0/r**3

    def make_interp_r(self,func,npts=100):
        rlo = 0.0
        rhi = self.zerocross
        step = (rhi - rlo)/(1.0*npts)
        r = N.arange(rlo, rhi+0.5*step, step)
        f = func(r)
        return(r,f)

    def make_interp_r2(self,func,npts=100):
        rlo = 0.0
        rhi = self.zerocross
        rlo2 = rlo**2
        rhi2 = rhi**2
        step = (rhi2 - rlo2)/(1.0*npts)
        r2 = N.arange(rlo2, rhi2+0.5*step, step)
        f = func(N.sqrt(r2))
        return(r2,f)

    def linterp(self, x, y, xEval):
        inRange = (xEval > x[0])*(xEval < x[-1])
        dx = (x[-1] - x[0])/(x.size - 1.0)
        indx = N.cast['int32']( (xEval-x[0])/dx*inRange )
        return (xEval - x[indx])/dx*(y[indx+1] - y[indx]) + y[indx]

    def rand_eval_r(self, rLinterp, fLinterp, npts=1000):
        rEval = NR.rand(npts)*(rLinterp[-1] - rLinterp[0]) + rLinterp[0]
        fEvalLinterp = self.linterp(rLinterp, fLinterp, rEval)
        indx = rEval.argsort()
        return (rEval[indx], fEvalLinterp[indx])

    def rand_err_r(self, rLinterp, fLinterp, fExact, fTotal, npts=1000):
        (rEval,fEvalLinterp) = self.rand_eval_r(rLinterp, fLinterp, npts)
        fEvalExact = fExact(rEval)
        #fEvalErr = N.abs(fEvalExact - fEvalLinterp)
        fEvalErr = fEvalLinterp - fEvalExact
        fEvalFracErr = fEvalErr/fTotal(rEval)
        return (rEval, fEvalFracErr)

    def rand_eval_r2(self, r2Linterp, fLinterp, npts=1000):
        rLinterp = N.sqrt(r2Linterp)
        rEval = NR.rand(npts)*(rLinterp[-1] - rLinterp[0]) + rLinterp[0]
        fEvalLinterp = self.linterp(r2Linterp, fLinterp, rEval**2)
        indx = rEval.argsort()
        return (rEval[indx], fEvalLinterp[indx])

    def rand_err_r2(self, r2Linterp, fLinterp, fExact, fTotal, npts=1000):
        rLinterp = N.sqrt(r2Linterp)
        (rEval,fEvalLinterp) = self.rand_eval_r2(r2Linterp,fLinterp,npts)
        fEvalExact = fExact(rEval)
        #fEvalErr = N.abs(fEvalExact - fEvalLinterp)
        fEvalErr = fEvalLinterp - fEvalExact
        fEvalFracErr = fEvalErr/fTotal(rEval)
        return (rEval, fEvalFracErr)

def main(argv):
    if len(argv) < 5:
        print('')
        print(argv[0] + ' <FG|FGOR> <R|R2> <n> <outName>')
        print('')
        print('FG = grid force, FGOR = grid force over r')
        print('R = linear in r, R2 = linear in r^2')
        print('')
        sys.exit(-1)

    forcetype = argv[1].upper()
    ordinatetype = argv[2].upper()
    n = N.int(argv[3])
    outName = argv[4]

    srf = ShortRangeForce()

    if forcetype == 'FG':
        force = srf.f_grid
    elif forcetype == 'FGOR':
        force = srf.f_gridOR
    else:
        print('ERROR: unknown forcetype ' + forcetype)
        sys.exit(-1)

    if ordinatetype == 'R':
        r,f = srf.make_interp_r(force,n)
    elif ordinatetype == 'R2':
        r,f = srf.make_interp_r2(force,n)
    else:
        print('ERROR: unknown ordinate type ' + ordinatetype)
        sys.exit(-1)

    of = open(outName,'w')
    for i in range(0,n):
        of.write(str(r[i]) + '\t' + str(f[i]) + '\n')
    of.close()

    return

if __name__ == '__main__':
    main(sys.argv)
