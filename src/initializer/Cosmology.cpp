/*
   Initializer:
   Cosmology.cpp

      This file has implementation for routines:

         GetParameters(DataBase)
      which gets cosmology parameters from the
      DataBase (defined in DataBase.h), and stores them 
      with this class as well. 

         Sigma_r(R, AA)
      which returns mass variance on a scale R, of the 
      linear density field with normalization AA.

         TransferFunction(k)
      which returns transfer function for the mode k in 
      Fourier space.

         GrowthFactor(z)
      which returns the linear growth factor at redshift z.

         GrowthDeriv(z)
      which returns the derivative of the linear growth factor 
      at redshift z.

                        Zarija Lukic, February 2009
                           zarija@lanl.gov
*/

#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <iostream>
#include <fstream>

#include "Cosmology.h"

void CosmoClass::SetParameters(GlobalStuff& DataBase, const char *tfName){
  std::ifstream inputFile;
  integer i, ln;
  real tmp, tfcdm, tfbar, norm;
  real akh1, akh2, alpha;
  real file_kmax = M_PI/DataBase.box_size * DataBase.ngrid;
  
  Omega_m = DataBase.Omega_m;
  Omega_bar = DataBase.Omega_bar;
  h = DataBase.Hubble;
  n_s = DataBase.n_s;
  TFFlag = DataBase.TFFlag;
  w_de = DataBase.w_de;
  
  cobe_temp=2.728;  // COBE/FIRAS CMB temperature in K
  tt=cobe_temp/2.7 * cobe_temp/2.7;
  
  if (TFFlag == 0) { // CMBFast transfer function
    //inputFile.open("cmb.tf", std::ios::in);
    inputFile.open(tfName, std::ios::in);
    if (!inputFile){
      printf("Cannot open CMBFast file!");
      abort();
    }
    ln = 0;
    while (! inputFile.eof()) {
      inputFile >> tmp >> tfcdm >> tfbar >> tmp >> tmp;
      //std::cout << tfcdm << " " << tfbar << std::endl;
      ++ln;
    }
    table_size = ln-1; // It also reads newline...
    table_kk = (real *)malloc(table_size*sizeof(real));
    table_tf = (real *)malloc(table_size*sizeof(real));
    inputFile.close();
    
    //inputFile.open("cmb.tf", std::ios::in);
    inputFile.open(tfName, std::ios::in);
    //inputFile.seekg(0, std::ios::beg);
    inputFile >> table_kk[0] >> tfcdm >> tfbar >> tmp >> tmp;
    table_tf[0] = tfbar*Omega_bar/Omega_m + 
      tfcdm*(Omega_m-Omega_bar)/Omega_m;
    norm = table_tf[0];
    for (i=1; i<table_size; ++i){
      inputFile >> table_kk[i] >> tfcdm >> tfbar >> tmp >> tmp;
      table_tf[i] = tfbar*Omega_bar/Omega_m + 
	tfcdm*(Omega_m-Omega_bar)/Omega_m;
      table_tf[i] = table_tf[i]/norm; // Normalization
      //std::cout << table_kk[0] << " " << table_tf[0] << std::endl;
    }
    last_k = table_kk[table_size-1];
    if (last_k < file_kmax) {
      std::cout << "CMBFast file goes only to k = " << last_k << std::endl;
      std::cout << "Aborting" << std::endl << std::flush;
      int rc = 1;
      MPI_Abort(MPI_COMM_WORLD, rc);
    }
    inputFile.close();
    tan_f=&CosmoClass::cmbfast;
  }
  else if (TFFlag == 1) { // Klypin-Holtzmann transfer function
    akh1=pow(46.9*Omega_m*h*h, 0.670)*(1.0+pow(32.1*Omega_m*h*h, -0.532));
    akh2=pow(12.0*Omega_m*h*h, 0.424)*(1.0+pow(45.0*Omega_m*h*h, -0.582));
    alpha=pow(akh1, -1.0*Omega_bar/Omega_m)*pow(akh2, pow(-1.0*Omega_bar/Omega_m, 3.0));
    kh_tmp = Omega_m*h*sqrt(alpha)*pow(1.0-Omega_bar/Omega_m, 0.6);
    tan_f=&CosmoClass::klypin_holtzmann;
    last_k = 10.0;
  }
  else if (TFFlag == 2) { // Hu-Sugiyama transfer function
    akh1=pow(46.9*Omega_m*h*h, 0.670)*(1.0+pow(32.1*Omega_m*h*h, -0.532));
    akh2=pow(12.0*Omega_m*h*h, 0.424)*(1.0+pow(45.0*Omega_m*h*h, -0.582));
    alpha=pow(akh1, -1.0*Omega_bar/Omega_m)*pow(akh2, pow(-1.0*Omega_bar/Omega_m, 3.0));
    kh_tmp = Omega_m*h*sqrt(alpha);
    tan_f = &CosmoClass::hu_sugiyama;
    last_k = 10.0;
  }
  else if (TFFlag == 3) { // Peacock-Dodds transfer function
    kh_tmp = Omega_m*h*exp(-2.0*Omega_bar);
    tan_f = &CosmoClass::peacock_dodds;
    last_k = 10.0;
  }
  else if (TFFlag == 4) { // BBKS transfer function
    kh_tmp = Omega_m*h;
    tan_f = &CosmoClass::bbks;
    last_k = 10.0;
  }
  else{
    printf("Cosmology::SetParameters: TFFlag has to be 0, 1, 2, 3 or 4!\n");
    abort();
  }
  
  return;
}

CosmoClass::CosmoClass() {
  table_kk = NULL;
  table_tf = NULL;
}

CosmoClass::~CosmoClass(){
  if(table_kk)
    free(table_kk);
  if(table_tf)
    free(table_tf);
}

#define pi M_PI

/* Transfer functions */

real CosmoClass::cmbfast(real k){
  real t_f;
  t_f = interpolate(table_kk, table_tf, table_size, k);
  return(t_f);
}

real CosmoClass::klypin_holtzmann(real k){
  real qkh, t_f;
  
  if (k == 0.0) return(0.0);
  qkh = k*tt/kh_tmp;
  /* NOTE: the following line has 0/0 for k=0.
     This was taken care of at the beginning of the routine. */
  t_f = log(1.0 + 2.34*qkh)/(2.34*qkh) * 
    pow(1.0
	+ 13.0*qkh 
	+ pow(10.5*qkh, 2.0)
	+ pow(10.4*qkh, 3.0) 
	+ pow(6.51*qkh, 4.0), 
	-0.25);
  return(t_f);
}

real CosmoClass::hu_sugiyama(real k){
  real qkh, t_f;
      
  if (k == 0.0) return(0.0);
  qkh = k*tt/kh_tmp;
  /* NOTE: the following line has 0/0 for k=0.
     This was taken care of at the beginning of the routine. */
  t_f = log(1.0 + 2.34*qkh)/(2.34*qkh) * 
    pow(1.0
	+ 3.89*qkh 
	+  pow(16.1*qkh, 2.0)
	+ pow(5.46*qkh, 3.0)
	+ pow(6.71*qkh, 4.0), 
	-0.25);
  return(t_f);
}

real CosmoClass::peacock_dodds(real k){
  real qkh, t_f;
  
  if (k == 0.0) return(0.0);
  qkh = k/kh_tmp;
  /* NOTE: the following line has 0/0 for k=0.
     This was taken care of at the beginning of the routine. */
  t_f = log(1.0 + 2.34*qkh)/(2.34*qkh) * 
    pow(1.0
	+ 3.89*qkh
	+ pow(16.1*qkh, 2.0)
	+ pow(5.46*qkh, 3.0)
	+ pow(6.71*qkh, 4.0), 
	-0.25);
  return(t_f);
}

real CosmoClass::bbks(real k){
  real qkh, t_f;
  
  if (k == 0.0) return(0.0);
  qkh = k/kh_tmp;
  /* NOTE: the following line has 0/0 for k=0.
     This was taken care of at the beginning of the routine. */
  t_f = log(1.0 + 2.34*qkh)/(2.34*qkh) * 
    pow(1.0
	+ 3.89*qkh
	+ pow(16.1*qkh, 2.0)
	+ pow(5.46*qkh, 3.0)
	+ pow(6.71*qkh, 4.0), 
	-0.25);
  return(t_f);
}

real CosmoClass::TransferFunction(real k){
  return((this->*tan_f)(k));
}

/* Sigma(R) routines */

real CosmoClass::sigma2(real k){
  real prim_ps, w_f, t_f, s2;
  
  prim_ps = Pk_norm*pow(k, n_s);
  w_f = 3.0*(sin(k*R_M) - k*R_M*cos(k*R_M))/pow(k*R_M, 3.0);
  t_f = TransferFunction(k);
  s2 = k*k * t_f*t_f * w_f*w_f * prim_ps;
  return(s2);
}

real CosmoClass::Sigma_r(real R, real AA){
  real sigma;
  const real k_min=0.0, k_max=last_k;
  
  R_M = R;
  Pk_norm = AA;
  
  sigma = 1.0/(2.0*pi*pi) * integrate(&CosmoClass::sigma2, k_min, k_max);
  sigma = sqrt(sigma);
  return(sigma);
}

#undef pi

/* Growth factor for flat wCDM cosmologies: */

void CosmoClass::GrowthFactor(real z, real *gf, real *g_dot){
  real x1, x2, dplus, ddot;
  const float redshift=200.0;
  real *ystart;
  
  x1 = 1.0/(1.0+100000.0);
  x2 = 1.0/(1.0+z);
  ystart = (real *)malloc(2*sizeof(real));
  ystart[0] = x1;
  ystart[1] = 0.0;
  
  odesolve(ystart, 2, x1, x2, 1.0e-6, 1.0e-6, &CosmoClass::growths, false);
  //printf("Dplus = %f;  Ddot = %f \n", ystart[0], ystart[1]);	
  
  dplus = ystart[0];
  ddot  = ystart[1];
  x1 = 1.0/(1.0+100000.0);
  x2 = 1.0;
  ystart[0] = x1;
  ystart[1] = 0.0;
  
  odesolve(ystart, 2, x1, x2, 1.0e-6, 1.0e-6, &CosmoClass::growths, false);	
  //printf("Dplus = %f;  Ddot = %f \n", ystart[0], ystart[1]);
  
  *gf    = dplus/ystart[0];
  *g_dot = ddot/ystart[0];
  //printf("\nGrowth factor = %f;  Derivative = %f \n", dplus/ystart[0], ddot/ystart[0]);
  free(ystart);
  
  return;
}

void CosmoClass::growths(real a, real y[], real dydx[]){
  real H;
  H = sqrt(Omega_m/(a*a*a) + (1.0-Omega_m)*pow(a, -3.0*(1.0+w_de)));
  dydx[0] = y[1]/(a*H);
  dydx[1] = -2.0*y[1]/a + 1.5*Omega_m*y[0]/(H*pow(a, 4.0));
  return;
}

#define MAXSTP 10000
#define TINY 1.0e-30
#define SIGN(a,b) ((b) >= 0.0 ? fabs(a) : -fabs(a))
void CosmoClass::odesolve(real ystart[], int nvar, real x1, real x2, 
			  real eps, real h1,
                          void (CosmoClass::*derivs)(real, real [], real []), 
			  bool print_stat)
{	
  int i, nstp, nok, nbad, feval;
  real x,hnext,hdid,h;
  real *yscal,*y,*dydx;
  const real hmin=0.0;
  
  feval = 0;
  yscal= (real *)malloc(nvar*sizeof(real));
  y= (real *)malloc(nvar*sizeof(real));
  dydx= (real *)malloc(nvar*sizeof(real));
  x=x1;
  h=SIGN(h1,x2-x1);
  nok = nbad = 0;
  for (i=0; i<nvar; ++i) {y[i]=ystart[i];}
  
  for (nstp=0; nstp<MAXSTP; ++nstp) {
    (this->*derivs)(x, y, dydx);
    ++feval;
    for (i=0; i<nvar; ++i)
    {yscal[i]=fabs(y[i])+fabs(dydx[i]*h)+TINY;}
    if ((x+h-x2)*(x+h-x1) > 0.0) h=x2-x;
    rkqs(y,dydx,nvar,&x,h,eps,yscal,&hdid,&hnext,&feval,derivs);
    if (hdid == h) ++nok; else ++nbad;
    if ((x-x2)*(x2-x1) >= 0.0) {
      for (i=0; i<nvar; ++i) {ystart[i]=y[i];}
      free(dydx);
      free(y);
      free(yscal);
      if (print_stat){
	printf("ODEsolve:\n");
	printf(" Evolved from x = %f to x = %f\n", x1, x2);
	printf(" successful steps: %d\n", nok);
	printf(" bad steps: %d\n", nbad);
	printf(" function evaluations: %d\n", feval);
      }				
      return;
    }
    if (fabs(hnext) <= hmin) {
      printf("Step size too small in ODEsolve");
      exit(1);
    }
    h=hnext;
  }
  printf("Too many steps in ODEsolve");
  exit(1);
}
#undef MAXSTP
#undef TINY
#undef SIGN

#define SAFETY 0.9
#define PGROW -0.2
#define PSHRNK -0.25
#define ERRCON 1.89e-4
static real maxarg1,maxarg2, minarg1, minarg2;
#define FMAX(a,b) (maxarg1=(a),maxarg2=(b),(maxarg1) > (maxarg2) ? (maxarg1) : (maxarg2))
#define FMIN(a,b) (minarg1=(a),minarg2=(b),(minarg1) < (minarg2) ? (minarg1) : (minarg2))
void CosmoClass::rkqs(real y[], real dydx[], int n, real *x, real htry, 
		      real eps, 
                      real yscal[], real *hdid, real *hnext, int *feval, 
                      void (CosmoClass::*derivs)(real, real [], real []))
{
  int i;
  real errmax,h,htemp,xnew,*yerr,*ytemp;
  
  yerr= (real *)malloc(n*sizeof(real));
  ytemp= (real *)malloc(n*sizeof(real));
  h=htry;
  
  for (;;) {
    rkck(y,dydx,n,*x,h,ytemp,yerr,derivs);
    *feval += 5;
    errmax=0.0;
    for (i=0; i<n; ++i) {errmax=FMAX(errmax,fabs(yerr[i]/yscal[i]));}
    errmax /= eps;
    if (errmax <= 1.0) break;
    htemp=SAFETY*h*pow(errmax,PSHRNK);
    h=(h >= 0.0 ? FMAX(htemp,0.1*h) : FMIN(htemp,0.1*h));
    xnew=(*x)+h;
    if (xnew == *x) {
      printf("Stepsize underflow in ODEsolve rkqs");
      exit(1);
    }
  }
  if (errmax > ERRCON) *hnext=SAFETY*h*pow(errmax,PGROW);
  else *hnext=5.0*h;
  *x += (*hdid=h);
  for (i=0; i<n; ++i) {y[i]=ytemp[i];}
  free(ytemp);
  free(yerr);
}
#undef SAFETY
#undef PGROW
#undef PSHRNK
#undef ERRCON
#undef FMAX
#undef FMIN

/* Cash-Karp Runge-Kutta Step, based on the 
work of Fehlberg, who ﬁgured out that six function evaluations could 
be used to determine two Runge-Kutta steps, one fourth-order and one 
ﬁfth-order. The diﬀerence between these estimates can be used as a 
truncation error for adjusting the stepsize. */
void CosmoClass::rkck(real y[], real dydx[], int n, real x, real h, 
		      real yout[], real yerr[], 
		      void (CosmoClass::*derivs)(real, real [], real []))
{
  int i;
  static real a2=0.2,a3=0.3,a4=0.6,a5=1.0,a6=0.875,b21=0.2,
    b31=3.0/40.0,b32=9.0/40.0,b41=0.3,b42 = -0.9,b43=1.2,
    b51 = -11.0/54.0, b52=2.5,b53 = -70.0/27.0,b54=35.0/27.0,
    b61=1631.0/55296.0,b62=175.0/512.0,b63=575.0/13824.0,
    b64=44275.0/110592.0,b65=253.0/4096.0,c1=37.0/378.0,
    c3=250.0/621.0,c4=125.0/594.0,c6=512.0/1771.0,
    dc5 = -277.00/14336.0;
  real dc1=c1-2825.0/27648.0,dc3=c3-18575.0/48384.0,
    dc4=c4-13525.0/55296.0,dc6=c6-0.25;
  real *ak2,*ak3,*ak4,*ak5,*ak6,*ytemp;
  
  ak2= (real *)malloc(n*sizeof(real));
  ak3= (real *)malloc(n*sizeof(real));
  ak4= (real *)malloc(n*sizeof(real));
  ak5= (real *)malloc(n*sizeof(real));
  ak6= (real *)malloc(n*sizeof(real));
  ytemp= (real *)malloc(n*sizeof(real));
  
  for (i=0; i<n; ++i)
    ytemp[i]=y[i]+b21*h*dydx[i];
  (this->*derivs)(x+a2*h,ytemp,ak2);
  for (i=0; i<n; ++i)
    ytemp[i]=y[i]+h*(b31*dydx[i]+b32*ak2[i]);
  (this->*derivs)(x+a3*h,ytemp,ak3);
  for (i=0; i<n; ++i)
    ytemp[i]=y[i]+h*(b41*dydx[i]+b42*ak2[i]+b43*ak3[i]);
  (this->*derivs)(x+a4*h,ytemp,ak4);
  for (i=0; i<n; ++i)
    ytemp[i]=y[i]+h*(b51*dydx[i]+b52*ak2[i]+b53*ak3[i]+b54*ak4[i]);
  (this->*derivs)(x+a5*h,ytemp,ak5);
  for (i=0; i<n; ++i)
    ytemp[i]=y[i]+h*(b61*dydx[i]+b62*ak2[i]+b63*ak3[i]+b64*ak4[i]+b65*ak5[i]);
  (this->*derivs)(x+a6*h,ytemp,ak6);
  for (i=0; i<n; ++i)
    yout[i]=y[i]+h*(c1*dydx[i]+c3*ak3[i]+c4*ak4[i]+c6*ak6[i]);
  for (i=0; i<n; ++i)
    yerr[i]=h*(dc1*dydx[i]+dc3*ak3[i]+dc4*ak4[i]+dc5*ak5[i]+dc6*ak6[i]);
  
  free(ytemp);
  free(ak6);
  free(ak5);
  free(ak4);
  free(ak3);
  free(ak2);
}


/* Numerical integration routines from Numerical Recipes 
   calling:
   x = integrate(&CosmoClass::func, a, b)
   
   where
      func -- is the (real) function whose integral is calculated, 
                                 and has to be member of the CosmoClass
      a    -- is the (real) lower boundary for integration
      b    -- is the (real) upper boundary for integration
*/

#define FUNC(x) ((this->*func)(x))
real CosmoClass::midpoint(real (CosmoClass::*func)(real), 
                          real a, real b, integer n)
{
  real x,tnm,sum,del,ddel;
  static real s;
  integer it,j;
  
  if (n == 1) {
    return (s=(b-a)*FUNC(0.5*(a+b)));
  } 
  else {
    it = 1;
    for(j=1; j<n-1; ++j) it *= 3;
    tnm=it;
    del=(b-a)/(3.0*tnm);
    ddel=del+del;
    x=a+0.5*del;
    sum=0.0;
    for (j=1; j<=it; ++j){
      sum += FUNC(x);
      x += ddel;
      sum += FUNC(x);
      x += del;
    }
    s=(s+(b-a)*sum/tnm)/3.0;
    return (s);
  }
}
#undef FUNC

#define EPS 1.0e-4
#define JMAX 20
#define JMIN 5
real CosmoClass::integrate(real (CosmoClass::*func)(real), real a, real b){
  real sol, old;
  integer j;
  
  sol = 0.0;
  old = -1.0e26;
  for (j=1; j<=JMAX; ++j) {
    sol=midpoint(func, a, b, j);
    if (j > JMIN){
      if (fabs(sol-old) < EPS*fabs(old) || 
	  (sol==0.0 && old==0.0)) return(sol);
    }
    old = sol;
  }
  printf("integrate: no convergence!");
  abort();
  return 0;
}
#undef EPS
#undef JMAX
#undef JMIN

/* Linear interpolation routine
   calling:
   y = interpolate(xx, yy, n, x)
   
   where
   xx -- is the (real) x array of ordered data
   yy -- is the (real) y array
   n  -- size of above arrays
   x  -- is the (real) point whose y value should be interpolated
*/

real CosmoClass::interpolate(real xx[], real yy[], unsigned long n, real x){
  real y, dx, dy;
  unsigned long jlo;
  
  locate(xx, n, x, &jlo);
  // Linear interpolation:
  dx = xx[jlo] - xx[jlo+1];
  dy = yy[jlo] - yy[jlo+1];
  y = dy/dx*(x-xx[jlo]) + yy[jlo];
  
  return(y);
}

/* Routines for locating a value in an ordered table 
   from Numerical Recipes */

void CosmoClass::locate(real xx[], unsigned long n, real x, unsigned long *j){
  unsigned long ju,jm,jl;
  int ascnd;
  
  jl=0;
  ju=n-1;
  ascnd=(xx[n-1] >= xx[0]);
  while (ju-jl > 1) {
    jm=(ju+jl)/2;
    if ((x >= xx[jm]) == ascnd)
      jl=jm;
    else
      ju=jm;
  }
  if (x == xx[0]) *j=0;
  else if(x == xx[n-1]) *j=n-2;
  else *j=jl;
  return;
}

void CosmoClass::hunt(real xx[], unsigned long n, real x, unsigned long *jlo){
  unsigned long jm,jhi,inc;
  int ascnd;
  
  ascnd=(xx[n-1] >= xx[0]);
  if (*jlo <= 0 || *jlo > n) {
    *jlo=0;
    jhi=n+1;
  } else {
    inc=1;
    if (x >= xx[*jlo] == ascnd) {
      if (*jlo == n) return;
      jhi=(*jlo)+1;
      while (x >= xx[jhi] == ascnd) {
	*jlo=jhi;
	inc += inc;
	jhi=(*jlo)+inc;
	if (jhi > n) {
	  jhi=n+1;
	  break;
	}
      }
    } else {
      if (*jlo == 1) {
	*jlo=0;
	return;
      }
      jhi=(*jlo)--;
      while (x < xx[*jlo] == ascnd) {
	jhi=(*jlo);
	inc <<= 1;
	if (inc >= jhi) {
	  *jlo=0;
	  break;
	}
	else *jlo=jhi-inc;
      }
    }
  }
  while (jhi-(*jlo) != 1) {
    jm=(jhi+(*jlo)) >> 1;
    if (x >= xx[jm] == ascnd)
      *jlo=jm;
    else
      jhi=jm;
  }
  if (x == xx[n]) *jlo=n-1;
  if (x == xx[1]) *jlo=1;
}
