#ifndef MC3_H
#define MC3_H

#include <iostream>
#include <string>
#include <cassert>
#include <algorithm>
#include <vector>

#include <stdlib.h>
#include <math.h>
#include <unistd.h>
#include <string.h>

#include <rru_mpi.h>

#include "BasicDefinition.h"
#include "Partition.h"
#include "ParticleExchange.h"
#include "GridExchange.h"
#include "ParticleDistribute.h"
#include "CosmoHaloFinderP.h"
#include "FOFHaloProperties.h"

#include "mc3types.h"

#include "Domain.h"
#include "Basedata.h"
#include "Particles.h"
#include "Halodata.h"
#include "skewers.h"
#include "lightcone.h"
#include "TimeStepper.h"

#include "SimpleTimings.h"
#include "MC3Options.h"
#include "MC3Extras.h"

#include "dims.h"
//#include "distribution.hpp"
#include <complex-type.h>
#include "solver.hpp"
#ifdef FFTW3
#ifdef ESSL_FFTW
#include <fftw3_essl.h>
#else
#include <fftw3-mpi.h>
#endif
#else
#include <fftw_mpi.h>
#endif

#include "mc3defines.h"

using namespace std;


//MISC
void initial_all_to_all(bool flag);


//MISC (COPY TO OPTIONS)
string create_outName(string outBase, int rank);


//MISC (MOVE TO OPTIONS)
/*
vector<int>* readInts(string inName);
int intInVec(int n, vector<int> *iv);
*/


//GRID
void output_array_alive(GRID_T *arr, string outName);
void output_array_total(GRID_T *arr, string outName);
double sum_rho_alive(GRID_T *arr);
void fft_copy_arr_forward(GRID_T *arr, COMPLEX_T *fft_arr);
void fft_copy_arr_backward(COMPLEX_T *fft_arr, GRID_T *arr);
void poisson_alloc(COMPLEX_T **fft_rho_arr,
		   COMPLEX_T **fft_grad_phi_arr,
		   COMPLEX_T **fft_phi_arr);
void poisson_free(COMPLEX_T **fft_rho_arr,
		  COMPLEX_T **fft_grad_phi_arr,
		  COMPLEX_T **fft_phi_arr);
void writePk(SolverBase *solver, string outName);
void map2_poisson_forward(Particles & particles,
			  SolverBase *solver,
			  GRID_T *rho_arr,
			  COMPLEX_T *fft_rho_arr);
void map2_poisson_backward_gradient(Particles & particles,
				    SolverBase *solver,
				    GRID_T *grad_phi_arr,
				    COMPLEX_T *fft_grad_phi_arr,
				    GridExchange & gexchange,
				    TimeStepper & ts,
				    TS_FLOAT stepFraction);
void map2_poisson_backward_potential(Particles & particles,
				     SolverBase *solver,
				     GRID_T *phi_arr,
				     COMPLEX_T *fft_phi_arr,
				     GridExchange & gexchange);


//PARTICLES
void loadParticles(Basedata & indat,
		   Particles & particles,
		   string inBase,
		   string dataType,
		   string distributeType,
		   int cosmoFormatQ,
		   MC3Options & options);
void coords_Mpc2Mpch(vector<POSVEL_T> xVec,
		     vector<POSVEL_T> yVec,
		     vector<POSVEL_T> zVec,
		     float hubble);
void coords_Mpch2Mpc(vector<POSVEL_T> xVec,
		     vector<POSVEL_T> yVec,
		     vector<POSVEL_T> zVec,
		     float hubble);


//PARTICLES (COPY TO OPTIONS)

enum {
  ReserveXV    = 1<<0,
  ReserveNonXV = 1<<1
};

void allocReserveVectors(vector<POSVEL_T> **xx,
			 vector<POSVEL_T> **yy,
			 vector<POSVEL_T> **zz,
			 vector<POSVEL_T> **vx,
			 vector<POSVEL_T> **vy,
			 vector<POSVEL_T> **vz,
			 vector<POSVEL_T> **mass,
			 vector<POSVEL_T> **phi,
			 vector<ID_T> **tag,
			 vector<MASK_T> **mask,
			 vector<STATUS_T> **status,
			 int Np,
                         unsigned mode = ReserveXV | ReserveNonXV);
void reserveVectors(vector<POSVEL_T> **xx,
		    vector<POSVEL_T> **yy,
		    vector<POSVEL_T> **zz,
		    vector<POSVEL_T> **vx,
		    vector<POSVEL_T> **vy,
		    vector<POSVEL_T> **vz,
		    vector<POSVEL_T> **mass,
		    vector<POSVEL_T> **phi,
		    vector<ID_T> **tag,
		    vector<MASK_T> **mask,
		    vector<STATUS_T> **status,
		    int Np,
                    unsigned mode = ReserveXV | ReserveNonXV);
void deleteVectors(vector<POSVEL_T> **xx,
		   vector<POSVEL_T> **yy,
		   vector<POSVEL_T> **zz,
		   vector<POSVEL_T> **vx,
		   vector<POSVEL_T> **vy,
		   vector<POSVEL_T> **vz,
		   vector<POSVEL_T> **mass,
		   vector<POSVEL_T> **phi,
		   vector<ID_T> **id,
		   vector<MASK_T> **mask,
		   vector<STATUS_T> **status);


//PARTICLES (MOVE TO OPTIONS)
/*
void findHalos(Particles & particles,
	       Basedata & indat,
	       Halodata & halodat,
	       string outBase,
	       float anow,
	       int step,
	       LightCone *lc,
	       int applyLC);
void staticSkewers(Particles & particles,
		   ParallelSkewerSet & zskewers, 
		   string outBase, 
		   float anow, 
		   int step);
void lightconeSkewers(Particles & particles,
		      ConvergentSkewerSet & lcskewers, 
		      string outBase, 
		      float anow, 
		      int step,
		      LightCone lc);
void refreshParticles(Particles & particles, float anow);
*/


//void timingStats(double t, const char *comment);

#endif
