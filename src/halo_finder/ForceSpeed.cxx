#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include <math.h>

#include <time.h>
#include <vector>

#include "ForceLaw.h"

using namespace std;

#define RSM 0.1
#define NINTERP 1024

//produces a (mostly) random permutation array
vector<int>* permutation(int n) {
  vector<int> *perm = new vector<int>;
  perm->reserve(n);
  perm->push_back(0);
  for(int i=1; i<n; i++) {
    int j = lrand48() % (i+1);
    perm->push_back( (*perm)[j]);
    (*perm)[j] = i;
  }
  return perm;
}

//figure out elapsed time in seconds
double elapsed(clock_t start, clock_t end) {
  return 1.0*(end-start)/CLOCKS_PER_SEC;
}

int main(int argc, char *argv[]) {

  clock_t start, end;
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

  printf("1/%d = %f second time resolution\n",
	 CLOCKS_PER_SEC, 
	 1.0/CLOCKS_PER_SEC);
  printf("%d force calculations\n",n);
  //printf("\n");

  //set up fitting formula and interpolated force laws
  FGrid *fg = new FGrid();
  float RMAX = fg->rmax();

  FGridEvalFit *fgoref = new FGridEvalFit(fg);
  ForceLawSR *flf = new ForceLawSR(fgoref, RSM);

  FGridEvalInterp *fgorei = new FGridEvalInterp(fg, NINTERP);
  ForceLawSR *fli = new ForceLawSR(fgorei, RSM);

  //pre-calculate random differential positions
  vector<float> *dx = new vector<float>;
  vector<float> *dy = new vector<float>;
  vector<float> *dz = new vector<float>;

  dx->reserve(n);
  dy->reserve(n);
  dz->reserve(n);

  for(int i=0; i<n; i++) {
    dx->push_back(2.0*RMAX*(drand48()-0.5));
    dy->push_back(2.0*RMAX*(drand48()-0.5));
    dz->push_back(2.0*RMAX*(drand48()-0.5));
  }

  //variables for force calculations
  float r2, f_over_r;
  float vx, vy, vz;

  //fitting formula, in-order
  start = clock();
  vx = vy = vz = 0.0;
  for(int i=0; i<n; i++) {
    r2 = (*dx)[i]*(*dx)[i] + (*dy)[i]*(*dy)[i] + (*dz)[i]*(*dz)[i];
    f_over_r = flf->f_over_r(r2);
    vx += f_over_r*(*dx)[i];
    vy += f_over_r*(*dy)[i];
    vz += f_over_r*(*dz)[i];
  }
  end = clock();
  printf("fit in-order: %f seconds\n",elapsed(start,end));

  //interpolation, in-order
  start = clock();
  vx = vy = vz = 0.0;
  for(int i=0; i<n; i++) {
    r2 = (*dx)[i]*(*dx)[i] + (*dy)[i]*(*dy)[i] + (*dz)[i]*(*dz)[i];
    f_over_r = fli->f_over_r(r2);
    vx += f_over_r*(*dx)[i];
    vy += f_over_r*(*dy)[i];
    vz += f_over_r*(*dz)[i];
  }
  end = clock();
  printf("int in-order: %f seconds\n",elapsed(start,end));

  //generate random permutation
  vector<int> *perm = permutation(n);

  //fitting formula, permuted
  start = clock();
  vx = vy = vz = 0.0;
  for(int j=0; j<n; j++) {
    int i = (*perm)[j];
    r2 = (*dx)[i]*(*dx)[i] + (*dy)[i]*(*dy)[i] + (*dz)[i]*(*dz)[i];
    f_over_r = flf->f_over_r(r2);
    vx += f_over_r*(*dx)[i];
    vy += f_over_r*(*dy)[i];
    vz += f_over_r*(*dz)[i];
  }
  end = clock();
  printf("fit permuted: %f seconds\n",elapsed(start,end));

  //interpolation, permuted
  start = clock();
  vx = vy = vz = 0.0;
  for(int j=0; j<n; j++) {
    int i = (*perm)[j];
    r2 = (*dx)[i]*(*dx)[i] + (*dy)[i]*(*dy)[i] + (*dz)[i]*(*dz)[i];
    f_over_r = fli->f_over_r(r2);
    vx += f_over_r*(*dx)[i];
    vy += f_over_r*(*dy)[i];
    vz += f_over_r*(*dz)[i];
  }
  end = clock();
  printf("int permuted: %f seconds\n",elapsed(start,end));

  return 0;
}
