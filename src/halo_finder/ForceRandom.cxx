#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

#include <time.h>

#include "ForceLaw.h"

using namespace std;

#define RSM 0.1
#define NINTERP 1024
#define FUDGE 1.2

int main(int argc, char *argv[]) {

  time_t seconds = time(NULL);

  if(argc < 2) {
    fprintf(stderr,"USAGE: %s <n> [seed]\n",argv[0]);
    exit(-1);
  }

  int n = atoi(argv[1]);

  //figure out how to do random number generator seeding
  int seed;
  if(argc > 2)
    seed = atoi(argv[2]);
  else
    seed = (int)seconds;
  srand48(seed);


  //float tmp = atof(argv[1]);

  //set up fitting formula and interpolated force laws
  FGrid *fg = new FGrid();
  float RMAX = fg->rmax();

  FGridEvalFit *fgoref = new FGridEvalFit(fg);
  ForceLawSR *flf = new ForceLawSR(fgoref, RSM);

  /*
  float *r2arr, *farr;
  int nInterp;
  nInterp = fgorei->nInterp();
  r2arr = fgorei->r2();
  farr = fgorei->f();
  for(int i=0; i<nInterp; i++) {
    printf("%f\t%f\n",r2arr[i],farr[i]);
  }
  */

  FGridEvalInterp *fgorei = new FGridEvalInterp(fg, NINTERP);
  ForceLawSR *fli = new ForceLawSR(fgorei, RSM);

  FGridEvalPoly *fgorep = new FGridEvalPoly(fg);
  ForceLawSR *flp = new ForceLawSR(fgorep, RSM);

  float r, r2, f_over_r_fit, f_over_r_interp, f_over_r_poly;

  for(int i=0; i<n; i++) {
    r = RMAX*FUDGE*drand48();
    r2 = r*r;
    f_over_r_fit = flf->f_over_r(r2);
    f_over_r_interp = fli->f_over_r(r2);
    f_over_r_poly = flp->f_over_r(r2);
    printf("%f\t%f\t%f\t%f\n", r, f_over_r_fit, f_over_r_interp, f_over_r_poly);
  }

  //printf("%f\t%f\t%f\n",tmp,flf->f_over_r(tmp),fli->f_over_r(tmp));

  return 0;
}
