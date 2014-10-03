/* Cosmic initial conditions, roughly following MC2 initialization routines.

   Final output (particle positions and velocities) are in grid units. 

                        Zarija Lukic, May 2009
                           zarija@lanl.gov
*/

#include <mpi.h>
#include <iostream>
#include <fstream>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "TypesAndDefs.h"
#include "Parallelization.h"
#include "DataBase.h"
#include "MT_Random.h"
#include "Cosmology.h"
#include "PerfMon.h"
#include "distribution.h"
#include "Basedata.h"

#define PKTESTS 0
#if PKTESTS == 1
#include "solver.hpp"
#define MAX_STR_LEN 1024
#endif

#ifndef FFTW_ADDR
#define FFTW_ADDR(X) reinterpret_cast<fftw_complex*>(&(X)[0])
#endif

#ifdef FFTW2
#include "fftw_mpi.h"
#endif
#ifdef FFTW3
#if PENCIL
#ifdef ESSL_FFTW
#include <fftw3_essl.h>
#else
#include <fftw3.h>
#endif
#else
#include <fftw3-mpi.h>
#endif
#endif

#define pi M_PI
#define Ndim 3  /* Dimension of the problem */

//#define NOT_REALLY_RANDOM
#ifdef NOT_REALLY_RANDOM
#define genrand_real() (0.5f)
#endif

class Initializer {

 public:

  void initParticles(::real* pos_x, ::real* pos_y, ::real* pos_z, 
		     ::real* vel_x, ::real* vel_y, ::real* vel_z, 
		     Basedata& bdata, const char *tfName, int useWN);
  
  Initializer();
  ~Initializer();

  void CreateMPI_FFTW_COMPLEX(MPI_Datatype *MPI_FFTW_COMPLEX);
  void initspec();
  void indens1D();
  void indens2D();
  void indens3D();
  void indens_single();
  void indens_whitenoise();
  void test_reality();
  void solve_gravity(integer axis);
  void set_particles(real z_in, real d_z, real ddot, integer axis);
  void output(integer axis, real* pos, real* vel);

 protected:

  ParallelTools *Parallel;
  GlobalStuff *DataBase;
  CosmoClass *Cosmology;

  real *Pk;
  my_fftw_complex *rho;
  my_fftw_complex *buf1;
  my_fftw_complex *buf3;
  long My_Ng;  // lenght of the above arrays = number of my grid points
  
  integer nx1, nx2, ny1, ny2, nz1, nz2; // end points of my domain
  integer ngx, ngy, ngz;                // # of zones in my domain
  
  integer MyPE, MasterPE, NumPEs;   // MPI stuff
  
#ifdef FFTW2
  fftwnd_mpi_plan plan_f, plan_b;
  int lnx, lxs, lnyt, lyst, lsize;
#endif

#ifdef FFTW3
#if PENCIL
  fftw_plan plan_f_x, plan_f_y, plan_f_z;
  fftw_plan plan_b_x, plan_b_y, plan_b_z;
#else
  fftw_plan plan_f, plan_b;
#endif
#endif

};


Initializer::Initializer() {
  Parallel = new ParallelTools();
  DataBase = new GlobalStuff();
  Cosmology = new CosmoClass();

  Pk = NULL;
  rho = NULL;
  buf1 = NULL;
  buf3 = NULL;
}


Initializer::~Initializer() {
  delete Parallel;
  delete DataBase;
  delete Cosmology;
}


void Initializer::CreateMPI_FFTW_COMPLEX(MPI_Datatype *MPI_FFTW_COMPLEX){
  MPI_Datatype typelist[3] = {MPI_FLOAT, MPI_FLOAT, MPI_UB};
  MPI_Aint displ[3];
  int blocklen[3] = {1, 1, 1};
  int i, base;
  
  /* compute displacements of structure components */	
  MPI_Address(&rho[0].re, displ);
  MPI_Address(&rho[0].im, displ+1);
  MPI_Address(rho+1, displ+2);
  base = displ[0];
  for (i=0; i<3; ++i) displ[i] -= base;
  
  /* build datatype describing structure */
  MPI_Type_struct(3, blocklen, displ, typelist, MPI_FFTW_COMPLEX);
  //MPI_Type_contiguous(2, MPI_FLOAT, MPI_FFTW_COMPLEX);
  MPI_Type_commit(MPI_FFTW_COMPLEX);
  return;
}


void Initializer::initspec(){
  const real tpi=2.0*pi;
  const real n_s = DataBase->n_s;
  const real sigma8=DataBase->Sigma_8;
  const integer ngrid=DataBase->ngrid;
  const integer nq=ngrid/2;
  const real tpiL=tpi/DataBase->box_size; // k0, physical units
  long i, j, k, index, k_i, k_j, k_k;
  real kk, trans_f, s8, norm;
  
  /* Set P(k)=T^2*k^n array.
     Linear array, taking care of the mirror symmetry
     (reality condition for the density field). */
  for (i=0; i<ngx; ++i){
    k_i = i+nx1;
    if (k_i >= nq) {k_i = -MOD(ngrid-k_i,ngrid);}
    for (j=0; j<ngy; ++j){
      k_j = j+ny1;
      if (k_j >= nq) {k_j = -MOD(ngrid-k_j,ngrid);}
      for (k=0; k<ngz; ++k){
	k_k = k+nz1;
	if (k_k >= nq) {k_k = -MOD(ngrid-k_k,ngrid);}
	index = (i*ngy+j)*ngz+k;
	kk = tpiL*sqrt(k_i*k_i+k_j*k_j+k_k*k_k);
	trans_f = Cosmology->TransferFunction(kk);
	Pk[index] = trans_f*trans_f*pow(kk, n_s);
      }
    }
  }
  
  // Sigma_8 normalization:
  s8 = Cosmology->Sigma_r(8.0, 1.0);
  norm = sigma8*sigma8/(s8*s8);
  s8 = Cosmology->Sigma_r(8.0, norm);
  if (MyPE == MasterPE)
    printf("\nsigma_8 = %f, target was %f\n", s8, sigma8);
  
  // Finally, P(k)=A*T^2*k^n:
  for (i=0; i<My_Ng; ++i) Pk[i] *= norm;
  
  return;
}


void Initializer::indens1D(){

#if PENCIL
  Parallel->ParallelError("indens: 1D not supported with PENCIL!", "stop");
#else

  const real tpi=2.0*pi;
  const integer ngrid=DataBase->ngrid;
  const integer nq=ngrid/2;
  long i, j, k, index, index_cc;
  int proc;
  integer is_conj, rn_size;
  real amp, rn;
  real *am_plane, *ph_plane;
  
  is_conj=Parallel->GetMirrorFlag();
  
  if (is_conj == 0){
    proc = Parallel->GetMyPE();
    init_random(DataBase->seed, proc);
    for (i=0; i<My_Ng; ++i){
      rn  = genrand_real();  // Amplitude
      amp = sqrt(-1.0*log(rn)*Pk[i]);
      rn  = genrand_real();  // Phase
      rho[i].re = amp*cos(tpi*rn);
      rho[i].im = amp*sin(tpi*rn);
    }
  }
  else if (is_conj == 1) {
    rn_size = ngy*ngz;
    am_plane = (real *)malloc(rn_size*sizeof(real));
    ph_plane = (real *)malloc(rn_size*sizeof(real));
    
    // First plane (i=0):
    proc = Parallel->GetMirror1PE();
    init_random(DataBase->seed, proc);
    for (i=0; i<rn_size; ++i){
      am_plane[i] = genrand_real();
      ph_plane[i] = genrand_real();
    }
    for (j=0; j<ngy; ++j){
      for (k=0; k<ngz; ++k){
	index = j*ngz+k;
	index_cc = MOD(ngrid-j,ngrid)*ngz + MOD(ngrid-k,ngrid);
	amp = sqrt(-1.0*log(am_plane[index_cc])*Pk[index]);
	rn = ph_plane[index_cc];
	rho[index].re = amp*cos(tpi*rn);
	rho[index].im = -1.0*amp*sin(tpi*rn);
      }
    }
    
    // Other planes:
    proc = Parallel->GetMirrorPE();
    init_random(DataBase->seed, proc);
    for (i=0; i<rn_size; ++i){ // Skip the first plane
      am_plane[i] = genrand_real();
      ph_plane[i] = genrand_real();
    }
    for (i=ngx-1; i>0; --i){ // Go through other planes
      for (k=0; k<rn_size; ++k){
	am_plane[k] = genrand_real();
	ph_plane[k] = genrand_real();
      }
      for (j=0; j<ngy; ++j){
	for (k=0; k<ngz; ++k){
	  index = (i*ngy+j)*ngz+k;
	  index_cc = MOD(ngrid-j,ngrid)*ngz + MOD(ngrid-k,ngrid);
	  amp = sqrt(-1.0*log(am_plane[index_cc])*Pk[index]);
	  rn = ph_plane[index_cc];
	  rho[index].re = amp*cos(tpi*rn);
	  rho[index].im = -1.0*amp*sin(tpi*rn);
	}
      }
    }
    free(am_plane);
    free(ph_plane);
  }
  else {
    Parallel->ParallelError("indens: wrong mirror flag!", "stop");
  }
  
  // 0 and nq planes:
  if (nx1 == 0 || nx1 == nq){
    i = 0;
    for (j=0; j<ngy; ++j){
      for (k=0; k<ngz; ++k){
	index = (i*ngy + j)*ngz + k;
	if (j < nq){
	  rn  = genrand_real();  // Amplitude
	  amp = sqrt(-1.0*log(rn)*Pk[index]);
	  rn  = genrand_real();  // Phase
	  rho[index].re = amp*cos(tpi*rn);
	  rho[index].im = amp*sin(tpi*rn);
	}
	else {
	  index_cc = (i*ngy + MOD(ngrid-j,ngrid))*ngz + MOD(ngrid-k,ngrid);
	  rho[index].re = rho[index_cc].re;
	  rho[index].im = -1.0*rho[index_cc].im;
	}
      }
    }
    
    for (j=0; j<2; ++j){
      for (k=0; k<ngz; ++k){
	index = (i*ngy + j*nq)*ngz + k;
	if (k < nq){
	  rn  = genrand_real();  // Amplitude
	  amp = sqrt(-1.0*log(rn)*Pk[index]);
	  rn  = genrand_real();  // Phase
	  rho[index].re = amp*cos(tpi*rn);
	  rho[index].im = amp*sin(tpi*rn);
	}
	else{
	  index_cc = (i*ngy + MOD(ngrid-j*nq,ngrid))*ngz + MOD(ngrid-k,ngrid);
	  rho[index].re = rho[index_cc].re;
	  rho[index].im = -1.0*rho[index_cc].im;
	}
      }
    }
    
    // Real values:
    index = (i*ngy + nq)*ngz + 0;
    rho[index].im = 0.0;
    index = (i*ngy + 0)*ngz + 0;
    rho[index].im = 0.0;
    index = (i*ngy + 0)*ngz + nq;
    rho[index].im = 0.0;
    index = (i*ngy + nq)*ngz + nq;
    rho[index].im = 0.0;
    if (nx1 == 0) rho[0].re = 0.0;
  }
  
  return;
#endif
}


void Initializer::indens2D(){

#if !PENCIL      
  Parallel->ParallelError("indens: 2D not supported without PENCIL!", "stop");
#else

  const real tpi=2.0*pi;
  const integer ngrid=DataBase->ngrid;
  const integer nq=ngrid/2;
  long i, j, jc, k, index, index_cc;
  int proc;
  integer is_conj, rn_size;
  real amp, rn;
  real *am_plane, *ph_plane;
  
  is_conj=Parallel->GetMirrorFlag();
  
  if (is_conj == 0){
    proc = Parallel->GetMyPE();
    init_random(DataBase->seed, proc);
    for (i=0; i<My_Ng; ++i){
      rn  = genrand_real();  // Amplitude
      amp = sqrt(-1.0*log(rn)*Pk[i]);
      rn  = genrand_real();  // Phase
      rho[i].re = amp*cos(tpi*rn);
      rho[i].im = amp*sin(tpi*rn);
    }
  }
  else if (is_conj == 1) {
    rn_size = ngy*ngz;
    
    if (nx1 == nq) {
      // We need to generate values for this plane so that it can
      // be symmetrized later.
      
      proc = Parallel->GetMyPE();
      init_random(DataBase->seed, proc);
      for (i=0; i<rn_size; ++i){
	rn  = genrand_real();  // Amplitude
	amp = sqrt(-1.0*log(rn)*Pk[i]);
	rn  = genrand_real();  // Phase
	rho[i].re = amp*cos(tpi*rn);
	rho[i].im = amp*sin(tpi*rn);
      }
    }
    
    am_plane = (real *)malloc(rn_size*sizeof(real));
    ph_plane = (real *)malloc(rn_size*sizeof(real));
    
    if (nx1 != nq) {
      // First plane (i=0), first y slice (j=0):
      // (comes from the first y slice of the first x plane on the conj. rank)
      proc = Parallel->GetMirrorPE(3);
      init_random(DataBase->seed, proc);
      for (i=0; i<rn_size; ++i){
	am_plane[i] = genrand_real();
	ph_plane[i] = genrand_real();
      }
      for (k=0; k<ngz; ++k){
	index = k;
	index_cc = MOD(ngrid-k,ngrid);
	amp = sqrt(-1.0*log(am_plane[index_cc])*Pk[index]);
	rn = ph_plane[index_cc];
	rho[index].re = amp*cos(tpi*rn);
	rho[index].im = -1.0*amp*sin(tpi*rn);
      }
      
      // First plane (i=0), other y slices (j>0):
      // (comes from the first x plane on the conj. rank)
      proc = Parallel->GetMirrorPE(1);
      init_random(DataBase->seed, proc);
      for (i=0; i<rn_size; ++i){
	am_plane[i] = genrand_real();
	ph_plane[i] = genrand_real();
      }
      for (j=ngy-1, jc = 1; j>0; --j, ++jc) {
	for (k=0; k<ngz; ++k){
	  index = j*ngz+k;
	  index_cc = jc*ngz + MOD(ngrid-k,ngrid);
	  amp = sqrt(-1.0*log(am_plane[index_cc])*Pk[index]);
	  rn = ph_plane[index_cc];
	  rho[index].re = amp*cos(tpi*rn);
	  rho[index].im = -1.0*amp*sin(tpi*rn);
	}
      }
    }
    
    // Other planes (i>0), first y slice (j=0)
    // (comes from the first y slice of the conj. plane on the conj. rank)
    proc = Parallel->GetMirrorPE(2);
    init_random(DataBase->seed, proc);
    for (i=0; i<rn_size; ++i){ // Skip the first plane
      am_plane[i] = genrand_real();
      ph_plane[i] = genrand_real();
    }
    for (i=ngx-1; i>0; --i){ // Go through other planes
      for (k=0; k<rn_size; ++k){
	am_plane[k] = genrand_real();
	ph_plane[k] = genrand_real();
      }
      for (k=0; k<ngz; ++k){
	index = i*ngy*ngz+k;
	index_cc = MOD(ngrid-k,ngrid);
	amp = sqrt(-1.0*log(am_plane[index_cc])*Pk[index]);
	rn = ph_plane[index_cc];
	rho[index].re = amp*cos(tpi*rn);
	rho[index].im = -1.0*amp*sin(tpi*rn);
      }
    }
    
    // Other planes (i>0), other y slices (j>0):
    proc = Parallel->GetMirrorPE(0);
    init_random(DataBase->seed, proc);
    for (i=0; i<rn_size; ++i){ // Skip the first plane
      am_plane[i] = genrand_real();
      ph_plane[i] = genrand_real();
    }
    for (i=ngx-1; i>0; --i){ // Go through other planes
      for (k=0; k<rn_size; ++k){
	am_plane[k] = genrand_real();
	ph_plane[k] = genrand_real();
      }
      for (j=ngy-1, jc = 1; j>0; --j, ++jc) { // Skip the first y slice
	for (k=0; k<ngz; ++k){
	  index = (i*ngy+j)*ngz+k;
	  index_cc = jc*ngz + MOD(ngrid-k,ngrid);
	  amp = sqrt(-1.0*log(am_plane[index_cc])*Pk[index]);
	  rn = ph_plane[index_cc];
	  rho[index].re = amp*cos(tpi*rn);
	  rho[index].im = -1.0*amp*sin(tpi*rn);
	}
      }
    }
    free(am_plane);
    free(ph_plane);
  }
  else {
    Parallel->ParallelError("indens: wrong mirror flag!", "stop");
  }
  
  // 0 and nq planes (unlike in the 1-D case, the conj. values are
  // on different nodes). For both the k1 = 0 and k1 = N/2 planes,
  // these are both the first plane on the respective node.
  if (nx1 == 0 || nx1 == nq){
    rn_size = ngy*ngz;
    am_plane = (real *)malloc(rn_size*sizeof(real));
    ph_plane = (real *)malloc(rn_size*sizeof(real));
    
    i = 0;
    
    // F[i,j,k] = F*[i,N-j,N-k]
    // first y slice (which will also be the first y slice
    // on the conj. rank)
    if (ny1 >= nq) {
      proc = Parallel->GetMirrorPE(5);
      init_random(DataBase->seed, proc);
      for (i=0; i<rn_size; ++i){
	am_plane[i] = genrand_real();
	ph_plane[i] = genrand_real();
      }
      for (k=0; k<ngz; ++k){
	index = k;
	index_cc = MOD(ngrid-k,ngrid);
	amp = sqrt(-1.0*log(am_plane[index_cc])*Pk[index]);
	rn = ph_plane[index_cc];
	rho[index].re = amp*cos(tpi*rn);
	rho[index].im = -1.0*amp*sin(tpi*rn);
      }
      
      // other y slices...
      proc = Parallel->GetMirrorPE(4);
      init_random(DataBase->seed, proc);
      for (i=0; i<rn_size; ++i){
	am_plane[i] = genrand_real();
	ph_plane[i] = genrand_real();
      }
      for (j=ngy-1, jc = 1; j>0; --j, ++jc) {
	for (k=0; k<ngz; ++k){
	  index = j*ngz+k;
	  index_cc = jc*ngz + MOD(ngrid-k,ngrid);
	  amp = sqrt(-1.0*log(am_plane[index_cc])*Pk[index]);
	  rn = ph_plane[index_cc];
	  rho[index].re = amp*cos(tpi*rn);
	  rho[index].im = -1.0*amp*sin(tpi*rn);
	}
      }
    }

    // F[i,j,k] = F*[i,j,N-k]
    // these values are already on our node (in our first y slice)
    if (ny1 == 0 || ny1 == nq) {
      for (k=0; k<ngz; ++k){
	index = k;
	if (k >= nq){
	  index_cc = MOD(ngrid-k,ngrid);
	  rho[index].re = rho[index_cc].re;
	  rho[index].im = -1.0*rho[index_cc].im;
	}
      }
    }
    
    // Real values: F[{0,N/2},{0,N/2},{0,N/2}]
    if (ny1 == 0 || ny1 == nq) {
      rho[0].im = 0.0;
      rho[nq].im = 0.0;
    }
    
    // Set the average to 0
    if (nx1 == 0 && ny1 == 0) rho[0].re = 0.0;
    
    free(am_plane);
    free(ph_plane);
  }

//#define INIT_DEBUG_2D_GRF  
#ifdef INIT_DEBUG_2D_GRF
  for (int p = 0; p < Parallel->GetNumPEs(); ++p) { 
    if (p == Parallel->GetMyPE()) {
      std::cout << "p = " << p << std::endl;
      std::cout << "x = " << nx1 << " : " << nx2 << " (" << ngx << ")" << std::endl;
      std::cout << "y = " << ny1 << " : " << ny2 << " (" << ngy << ")" << std::endl;
      std::cout << "z = " << nz1 << " : " << nz2 << " (" << ngz << ")" << std::endl;
      for (int j = 0; j < 6; ++j)
	std::cout << "mirror " << j << " : " << Parallel->GetMirrorPE(j) << std::endl;
      std::cout << std::flush;
    }
    MPI_Barrier(MPI_COMM_WORLD);
  }
  
  MPI_Barrier(MPI_COMM_WORLD);
  
  for (int p = 0; p < Parallel->GetNumPEs(); ++p) { 
    if (p == Parallel->GetMyPE())
      for (int i = 0; i < ngx; ++i)
	for (int j = 0; j < ngy; ++j)
	  for (int k = 0; k < ngz; ++k) {
	    int idx = (i*ngy + j)*ngz + k;
	    std::cout << nx1 + i << " " << ny1 + j << " " << nz1 + k << " (" << idx << "): " << rho[idx].re << ", " << rho[idx].im << std::endl;
	  }
    std::cout << std::flush;
    MPI_Barrier(MPI_COMM_WORLD);
  }
#endif // INIT_DEBUG_2D_GRF
#endif
  return;
}


void Initializer::indens3D(){
  const real tpi=2.0*pi;
  const integer ng=DataBase->ngrid;
  const integer nq=ng/2;
  long i, j, k, index, index_cc;
  int proc;
  
  Parallel->ParallelError("indens: 3D not implemented yet!", "stop");
  
  return;
}


void Initializer::indens_single(){
  const real tpi=2.0*pi;
  const integer ng=DataBase->ngrid;
  const integer nq=ng/2;
  long i, j, k, index, index_cc;
  long tmp_i, tmp_j;
  int proc;
  real amp, rn;
  
  proc = Parallel->GetMyPE();
  init_random(DataBase->seed, proc);
  
  // Set two halves of the domain:
  for (i=0; i<ng; ++i){
    for (j=0; j<ng; ++j){
      for (k=0; k<ng; ++k){
	index = (i*ng + j)*ng + k;
	if (i < nq){
	  rn  = genrand_real();  // Amplitude
	  amp = sqrt(-1.0*log(rn)*Pk[index]);
	  rn  = genrand_real();  // Phase
	  rho[index].re = amp*cos(tpi*rn);
	  rho[index].im = amp*sin(tpi*rn);
	}
	else {
	  index_cc = (MOD(ng-i,ng)*ng + MOD(ng-j,ng))*ng + MOD(ng-k,ng);
	  rho[index].re = rho[index_cc].re;
	  rho[index].im = -1.0*rho[index_cc].im;
	}
      }
    }
  }
  
  // 0 & nq planes:
  for (i=0; i<2; ++i) {
    tmp_i = i*nq;
    for (j=0; j<ng; ++j){
      for (k=0; k<ng; ++k){
	index = (tmp_i*ng + j)*ng + k;
	if (j < nq){
	  rn  = genrand_real();  // Amplitude
	  amp = sqrt(-1.0*log(rn)*Pk[index]);
	  rn  = genrand_real();  // Phase
	  rho[index].re = amp*cos(tpi*rn);
	  rho[index].im = amp*sin(tpi*rn);
	}
	else {
	  index_cc = (tmp_i*ng + MOD(ng-j,ng))*ng + MOD(ng-k,ng);
	  rho[index].re = rho[index_cc].re;
	  rho[index].im = -1.0*rho[index_cc].im;
	}
      }
    }
    for (j=0; j<2; ++j){
      tmp_j = j*nq;
      for (k=0; k<ng; ++k){
	index = (tmp_i*ng + tmp_j)*ng + k;
	if (k <= nq){ /* set re part of rho[nq,nq.nq] here too */
	  rn  = genrand_real();  // Amplitude
	  amp = sqrt(-1.0*log(rn)*Pk[index]);
	  rn  = genrand_real();  // Phase
	  rho[index].re = amp*cos(tpi*rn);
	  rho[index].im = amp*sin(tpi*rn);
	}
	else{
	  index_cc = (tmp_i*ng + tmp_j)*ng + MOD(ng-k,ng);
	  rho[index].re = rho[index_cc].re;
	  rho[index].im = -1.0*rho[index_cc].im;
	}
      }
    }
  }
  
  // Real values:
  rho[0].re = 0.0;
  rho[0].im = 0.0;
  index = (nq*ng + 0)*ng + 0;
  rho[index].im = 0.0;
  index = (0*ng + nq)*ng + 0;
  rho[index].im = 0.0;
  index = (0*ng + 0)*ng + nq;
  rho[index].im = 0.0;
  index = (nq*ng + nq)*ng + 0;
  rho[index].im = 0.0;
  index = (0*ng + nq)*ng + nq;
  rho[index].im = 0.0;
  index = (nq*ng + 0)*ng + nq;
  rho[index].im = 0.0;
  index = (nq*ng + nq)*ng + nq;
  rho[index].im = 0.0;
  
  return;
}


void Initializer::indens_whitenoise(){
  const integer ng=DataBase->ngrid;
  long i;
  int proc;
  double rn, rn1, rn2, scal;

  // Generate Gaussian white noise field:
  proc = Parallel->GetMyPE();
  init_random(DataBase->seed, proc);
  for (i=0; i<My_Ng; ++i){
    do {
      rn1 = -1 + 2*genrand_real();
      rn2 = -1 + 2*genrand_real();
      rn = rn1*rn1 + rn2*rn2;
    } while (rn > 1.0 || rn == 0);
    rho[i].re = rn2 * sqrt(-2.0*log(rn)/rn);
    rho[i].im = 0.0;
  }

  // Go to the k-space:
#if defined (FFTW2) || defined (FFTW3)
  MPI_Barrier(MPI_COMM_WORLD);
#ifdef FFTW2
  fftwnd_mpi(plan_f, 1, (fftw_complex *)rho, NULL, FFTW_TRANSPOSED_ORDER);
#endif
#ifdef FFTW3
#if PENCIL
  Distribution &d = Parallel->d();
  d.redistribute_3_to_2((const complex_t *) rho,  (complex_t *) buf1, 0); // rho  --> buf1
  fftw_execute(plan_f_x);                                                 // buf1 --> buf3
  d.redistribute_2_to_3((const complex_t *) buf3, (complex_t *) buf1, 0); // buf3 --> buf1
  d.redistribute_3_to_2((const complex_t *) buf1, (complex_t *) buf3, 1); // buf1 --> buf3
  fftw_execute(plan_f_y);                                                 // buf3 --> buf1
  d.redistribute_2_to_3((const complex_t *) buf1, (complex_t *) buf3, 1); // buf1 --> buf3
  d.redistribute_3_to_2((const complex_t *) buf3, (complex_t *) buf1, 2); // buf3 --> buf1
  fftw_execute(plan_f_z);                                                 // buf1 --> rho
#else
  fftw_execute(plan_f);
#endif
#endif
  if (MyPE == MasterPE) std::cout << std::endl << "FFT DONE!" << std::endl;
#else
  Parallel->ParallelError("Something bad happened in indens_whitenoise", "no stop");
#endif

  // Multiply by the power spectrum:
  scal = pow(ng, 1.5);
  for (i=0; i<My_Ng; ++i){
    rho[i].re *= sqrt(Pk[i])/scal;
    rho[i].im *= sqrt(Pk[i])/scal;
  }

  // Fix k=0 point:
  if (nx1 == 0 && ny1 == 0 && nz1 == 0) rho[0].re = rho[0].im = 0.0;

  return;
}


void Initializer::test_reality(){
  long i;
  const integer ngrid=DataBase->ngrid;
  
  float test_zero, max_rho, min_rho, ave_rho, scal;

  min_rho = 1.0e30;
  max_rho = -1.0e30;
  ave_rho = 0.0;
  for (i=0; i<My_Ng; ++i) {
    ave_rho += (real)rho[i].re;
    if ((real)rho[i].re > max_rho) max_rho = (real)rho[i].re;
    if ((real)rho[i].re < min_rho) min_rho = (real)rho[i].re;
  }
  ave_rho   = ave_rho/My_Ng;
  
  MPI_Allreduce(MPI_IN_PLACE, &max_rho, 1, MPI_FLOAT, MPI_MAX, MPI_COMM_WORLD);
  MPI_Allreduce(MPI_IN_PLACE, &min_rho, 1, MPI_FLOAT, MPI_MIN, MPI_COMM_WORLD);
  MPI_Allreduce(MPI_IN_PLACE, &ave_rho, 1, MPI_FLOAT, MPI_SUM, MPI_COMM_WORLD);
  if (MyPE == MasterPE) {
    std::cout << std::endl << "Min and max value of density in k space: "
	      << min_rho << " " << max_rho << std::endl;
    ave_rho = ave_rho/NumPEs;
    std::cout << "Average value of density in k space: "
	      << ave_rho << std::endl << std::flush;
  }
  
#if defined (FFTW2) || defined (FFTW3)
  MPI_Barrier(MPI_COMM_WORLD);
#ifdef FFTW2
  fftwnd_mpi(plan_b, 1, (fftw_complex *)rho, NULL, FFTW_TRANSPOSED_ORDER);
#endif
#ifdef FFTW3
#if PENCIL
  Distribution &d = Parallel->d();
  fftw_execute(plan_b_z);                                                 // rho  --> buf3
  d.redistribute_2_to_3((const complex_t *) buf3, (complex_t *) buf1, 2); // buf3 --> buf1
  d.redistribute_3_to_2((const complex_t *) buf1, (complex_t *) buf3, 1); // buf1 --> buf3
  fftw_execute(plan_b_y);                                                 // buf3 --> buf1
  d.redistribute_2_to_3((const complex_t *) buf1, (complex_t *) buf3, 1); // buf1 --> buf3
  d.redistribute_3_to_2((const complex_t *) buf3, (complex_t *) buf1, 0); // buf3 --> buf1
  fftw_execute(plan_b_x);                                                 // buf1 --> buf3
  d.redistribute_2_to_3((const complex_t *) buf3, (complex_t *)  rho, 0); // buf3 --> rho
#else
  fftw_execute(plan_b);
#endif
#endif
  if (MyPE == MasterPE) std::cout << std::endl << "FFT DONE!" << std::endl;
  
  scal = pow(DataBase->box_size, 1.5); // inverse FFT scaling
  for (i=0; i< My_Ng; ++i) rho[i].re = rho[i].re/scal;

 

// start P(k) test
#if PKTESTS == 1
#if !PENCIL
  if (NumPEs > 1) {
    my_fftw_complex *data = (my_fftw_complex *)malloc(My_Ng*sizeof(my_fftw_complex));
    Parallel->d().redistribute_1_to_3((const complex_t *) rho, (complex_t *) data);
    for (i=0; i<My_Ng; ++i){rho[i] = data[i];}
    free(data);
  }
#endif

  SolverDiscrete *solver = new SolverDiscrete(MPI_COMM_WORLD, ngrid);
  solver->forward_solve((complex_t *)rho);
  std::vector<double> power;
  solver->power_spectrum(power);
  
  if(MyPE == MasterPE) {
    float rL = DataBase->box_size;
    float tpi = 2.0*atan(1.0)*4.0;
    float kcoeff = tpi/rL;
    float pkcoeff = powf(rL/ngrid,3.0);
    
    FILE *outFile = fopen("density.pk","w");
    for(int i=0; i<power.size(); i++) {
      fprintf(outFile,"%e\t%e\n",kcoeff*i, pkcoeff*power[i]);
    }
    fclose(outFile);
  }
  MPI_Barrier(MPI_COMM_WORLD);

#if !PENCIL
  if (NumPEs > 1) {
    my_fftw_complex *data = (my_fftw_complex *)malloc(My_Ng*sizeof(my_fftw_complex));
    Parallel->d().redistribute_3_to_1((const complex_t *) rho, (complex_t *) data);
    for (i=0; i<My_Ng; ++i){rho[i] = data[i];}
    free(data);
  }
#endif
#endif
// end P(k) test


  
  min_rho = 1.0e30;
  max_rho = -1.0e30;
  test_zero = 0.0;
  ave_rho = 0.0;
  for (i=0; i<My_Ng; ++i) {
    if (test_zero < fabs((real)rho[i].im))
      test_zero = fabs((real)rho[i].im);
    ave_rho += (real)rho[i].re;
    if ((real)rho[i].re > max_rho) max_rho = (real)rho[i].re;
    if ((real)rho[i].re < min_rho) min_rho = (real)rho[i].re;
  }
  ave_rho = ave_rho/My_Ng;
  
  MPI_Allreduce(MPI_IN_PLACE, &test_zero, 1, MPI_FLOAT, MPI_MAX, MPI_COMM_WORLD);	
  MPI_Allreduce(MPI_IN_PLACE, &max_rho, 1, MPI_FLOAT, MPI_MAX, MPI_COMM_WORLD);
  MPI_Allreduce(MPI_IN_PLACE, &min_rho, 1, MPI_FLOAT, MPI_MIN, MPI_COMM_WORLD);
  MPI_Allreduce(MPI_IN_PLACE, &ave_rho, 1, MPI_FLOAT, MPI_SUM, MPI_COMM_WORLD);
  if (MyPE == MasterPE) {
    std::cout << std::endl << "Max value of the imaginary part of density is "
	      << test_zero << std::endl;
    std::cout << "Min and max value of density in code units: "
	      << min_rho << " " << max_rho << std::endl;
    ave_rho = ave_rho/NumPEs;
    std::cout << "Average value of density in code units: "
	      << ave_rho << std::endl << std::flush;
  }
  
#else
  Parallel->ParallelError("Cannot test rho(r) without FFT", "no stop");
#endif
  
  return;
}


void Initializer::solve_gravity(integer axis){ 
// 0=x, 1=y, 2=z
  long i, j, k, k_i, k_j, k_k, index;
  const real tpi=2.0*pi;
  const integer ngrid=DataBase->ngrid;
  const integer nq=ngrid/2;
  real kk, Grad, Green, scal;
  real xx, yy, zz;
  float min_F, max_F, ave_F, test_zero;
  my_fftw_complex temp;
  
  /* Multiply by the Green's function (-1/k^2), 
     and take the gradient in k-space (-ik_i) */
  switch (axis) {
    case 0:
      xx = 1.0;
      yy = 0.0;
      zz = 0.0;
      break;
    case 1:
      xx = 0.0;
      yy = 1.0;
      zz = 0.0;
      break;
    case 2:
      xx = 0.0;
      yy = 0.0;
      zz = 1.0;
      break;
  }
  for (i=0; i<ngx; ++i){
    k_i = i+nx1;
    if (k_i >= nq) {k_i = -MOD(ngrid-k_i,ngrid);}
    for (j=0; j<ngy; ++j){
      k_j = j+ny1;
      if (k_j >= nq) {k_j = -MOD(ngrid-k_j,ngrid);}
      for (k=0; k<ngz; ++k){
	k_k = k+nz1;
	if (k_k >= nq) {k_k = -MOD(ngrid-k_k,ngrid);}
	Grad = -tpi/ngrid*(k_i*xx+k_j*yy+k_k*zz);
	index = (i*ngy+j)*ngz+k;
	kk = tpi/ngrid*sqrt(k_i*k_i+k_j*k_j+k_k*k_k);
	Green = -1.0/(kk*kk);
	if (kk == 0.0) Green = 0.0;
	temp = rho[index];
	rho[index].re = -Grad * Green * temp.im;
	rho[index].im = Grad * Green * temp.re;
      }
    }
  }
  // Fix k=0 point:
  if (nx1 == 0 && ny1 == 0 && nz1 == 0) rho[0].re = rho[0].im = 0.0;
 
#if defined (FFTW2)	|| defined (FFTW3)
  MPI_Barrier(MPI_COMM_WORLD);
#ifdef FFTW2	
  fftwnd_mpi(plan_b, 1, (fftw_complex *)rho, NULL, FFTW_NORMAL_ORDER);
#endif
#ifdef FFTW3
#if PENCIL
  Distribution &d = Parallel->d();
  fftw_execute(plan_b_z);                                                 // rho  --> buf3
  d.redistribute_2_to_3((const complex_t *) buf3, (complex_t *) buf1, 2); // buf3 --> buf1
  d.redistribute_3_to_2((const complex_t *) buf1, (complex_t *) buf3, 1); // buf1 --> buf3
  fftw_execute(plan_b_y);                                                 // buf3 --> buf1
  d.redistribute_2_to_3((const complex_t *) buf1, (complex_t *) buf3, 1); // buf1 --> buf3
  d.redistribute_3_to_2((const complex_t *) buf3, (complex_t *) buf1, 0); // buf3 --> buf1
  fftw_execute(plan_b_x);                                                 // buf1 --> buf3
  d.redistribute_2_to_3((const complex_t *) buf3, (complex_t *)  rho, 0); // buf3 --> rho
#else
  fftw_execute(plan_b);
#endif
#endif
  if (MyPE == MasterPE) std::cout << std::endl << "FFT DONE!" << std::endl;
  
  /* compensation for inverse FFT scaling */
  scal = pow(DataBase->box_size, 1.5);
  for (i=0; i<My_Ng; ++i) rho[i].re = rho[i].re/scal;
#else
  Parallel->ParallelError("FFT is not done", "no stop");
#endif



// start P(k) test
#if PKTESTS == 1
#if !PENCIL
  if (NumPEs > 1) {
    my_fftw_complex *data = (my_fftw_complex *)malloc(My_Ng*sizeof(my_fftw_complex));
    Parallel->d().redistribute_1_to_3((const complex_t *) rho, (complex_t *) data);
    for (i=0; i<My_Ng; ++i){rho[i] = data[i];}
    free(data);
  }
#endif

  SolverDiscrete *solver = new SolverDiscrete(MPI_COMM_WORLD, ngrid);
  solver->forward_solve((complex_t *)rho);
  std::vector<double> power;
  solver->power_spectrum(power);
  
  if(MyPE == MasterPE) {
    float rL = DataBase->box_size;
    float tpi = 2.0*atan(1.0)*4.0;
    float kcoeff = tpi/rL;
    float pkcoeff = powf(rL/ngrid,3.0);
    
    char outName[MAX_STR_LEN];
    sprintf(outName,"force_%d.pk",axis);
    FILE *outFile = fopen(outName,"w");
    for(int i=0; i<power.size(); i++) {
      fprintf(outFile,"%e\t%e\n",kcoeff*i, pkcoeff*power[i]);
    }
    fclose(outFile);
  }
  MPI_Barrier(MPI_COMM_WORLD);

#if !PENCIL
  if (NumPEs > 1) {
    my_fftw_complex *data = (my_fftw_complex *)malloc(My_Ng*sizeof(my_fftw_complex));
    Parallel->d().redistribute_3_to_1((const complex_t *) rho, (complex_t *) data);
    for (i=0; i<My_Ng; ++i){rho[i] = data[i];}
    free(data);
  }
#endif
#endif
// end P(k) test


 
#ifdef TESTING
  OnClock("Tests");
  min_F = 1.0e30;
  max_F = -1.0e30;
  test_zero = 0.0;
  ave_F = 0.0;
  for (i=0; i<My_Ng; ++i) { 
    if (test_zero < fabs((float)rho[i].im))
      test_zero = fabs((float)rho[i].im);  
    ave_F += (float)rho[i].re;
    if ((float)rho[i].re > max_F) max_F = (float)rho[i].re;
    if ((float)rho[i].re < min_F) min_F = (float)rho[i].re;
  }
  ave_F = ave_F/My_Ng;
  
  MPI_Allreduce(MPI_IN_PLACE, &test_zero, 1, MPI_FLOAT, MPI_MAX, MPI_COMM_WORLD);
  MPI_Allreduce(MPI_IN_PLACE, &max_F, 1, MPI_FLOAT, MPI_MAX, MPI_COMM_WORLD);
  MPI_Allreduce(MPI_IN_PLACE, &min_F, 1, MPI_FLOAT, MPI_MIN, MPI_COMM_WORLD);
  MPI_Allreduce(MPI_IN_PLACE, &ave_F, 1, MPI_FLOAT, MPI_SUM, MPI_COMM_WORLD);
  if (MyPE == MasterPE) {
    std::cout << std::endl << "Max value of the imaginary part of the force is "
	      << test_zero << std::endl;
    std::cout << "Min and max value of force in code units: "
	      << min_F << " " << max_F << std::endl;
    ave_F = ave_F/NumPEs;
    std::cout << "Average value of force in code units: "
	      << ave_F << std::endl << std::flush;
  }
  OffClock("Tests");
#endif
  
  return;
}


void Initializer::set_particles(real z_in, real d_z, real ddot, integer axis){ 
// 0=x, 1=y, 2=z
  const integer ngrid=DataBase->ngrid;
  long i, j, k, index;
  real pos_0, xx, yy, zz;
  float  move, max_move, ave_move;
  my_fftw_complex *F_i;
  
  F_i = rho; /* After gravity solve, this complex array 
		holds F_i component of the force, not density */
  
  switch (axis) {
    case 0:
      xx = 1.0;
      yy = 0.0;
      zz = 0.0;
      break;
    case 1:
      xx = 0.0;
      yy = 1.0;
      zz = 0.0;
      break;
    case 2:
      xx = 0.0;
      yy = 0.0;
      zz = 1.0;
      break;
  }
  
  /* Move particles according to Zeldovich approximation 
     particles will be put in rho array; real array will 
     hold positions, imaginary will hold velocities */
#ifdef TESTING
  OnClock("Tests");
  max_move = 0.0;
  ave_move = 0.0;
  for (i=0; i<My_Ng; ++i) {
    move = fabs(d_z*(float)F_i[i].re);
    ave_move += move;
    if (move > max_move) max_move = move;
  }
  ave_move = ave_move/My_Ng;
  
  MPI_Allreduce(MPI_IN_PLACE, &max_move, 1, MPI_FLOAT, MPI_MAX, MPI_COMM_WORLD);
  MPI_Allreduce(MPI_IN_PLACE, &ave_move, 1, MPI_FLOAT, MPI_SUM, MPI_COMM_WORLD);
  if (MyPE == MasterPE) {
    if (axis == 0) {std::cout << "Max move in X: " << max_move << std::endl;}
    else if (axis == 1) {std::cout << "Max move in Y: " << max_move << std::endl;}
    else {std::cout << "Max move in Z: " << max_move << std::endl;}
    ave_move = ave_move/NumPEs;
    if (axis == 0) {std::cout << "Average move in X: " << ave_move << std::endl;}
    else if (axis == 1) {std::cout << "Average move in Y: " << ave_move << std::endl;}
    else {std::cout << "Average move in Z: " << ave_move << std::endl;}
  }
  OffClock("Tests");
#endif

  // WAIT! 
  // slab: was in slabs in k space
  // slab: now in slabs in r space
  // slab can use common indexing in k and r
  // pencil: was in pencils in k space
  // pencil: now in cubes in r space
  // pencil must use different indexing in k and r

  int tnx1, tny1, tnz1;
  int tngx, tngy, tngz;
  
#if PENCIL
  Distribution &d = Parallel->d();

  tngx = d.local_ng_3d(0);
  tnx1 = d.self_3d(0)*tngx;

  tngy = d.local_ng_3d(1);
  tny1 = d.self_3d(1)*tngy;

  tngz = d.local_ng_3d(2);
  tnz1 = d.self_3d(2)*tngz;

#else
  tnx1 = nx1;
  tny1 = ny1;
  tnz1 = nz1;
  tngx = ngx;
  tngy = ngy;
  tngz = ngz;
#endif

  for (i=0; i<tngx; ++i){
    for (j=0; j<tngy; ++j){
      for (k=0; k<tngz; ++k){
	index = (i*tngy+j)*tngz+k;
	pos_0 = (i+tnx1)*xx + (j+tny1)*yy + (k+tnz1)*zz;
	F_i[index].im = -ddot*F_i[index].re;
	F_i[index].re = pos_0 - d_z*F_i[index].re;
//      if (F_i[index].re < 0.0) {F_i[index].re += (double)ngrid;}
//      else if (F_i[index].re >= ngrid) {F_i[index].re -= (double)ngrid;}
      }
    }
  }



// start P(k) test
#if PKTESTS == 1
#if !PENCIL
  if (NumPEs > 1) {
    my_fftw_complex *data = (my_fftw_complex *)malloc(My_Ng*sizeof(my_fftw_complex));
    Parallel->d().redistribute_1_to_3((const complex_t *) rho, (complex_t *) data);
    for (i=0; i<My_Ng; ++i){rho[i] = data[i];}
    free(data);
  }
#endif

  SolverDiscrete *solver = new SolverDiscrete(MPI_COMM_WORLD, ngrid);
  solver->forward_solve((complex_t *)rho);
  std::vector<double> power;
  solver->power_spectrum(power);
  
  if(MyPE == MasterPE) {
    float rL = DataBase->box_size;
    float tpi = 2.0*atan(1.0)*4.0;
    float kcoeff = tpi/rL;
    float pkcoeff = powf(rL/ngrid,3.0);
    
    char outName[MAX_STR_LEN];
    sprintf(outName,"part_%d.pk",axis);
    FILE *outFile = fopen(outName,"w");
    for(int i=0; i<power.size(); i++) {
      fprintf(outFile,"%e\t%e\n",kcoeff*i, pkcoeff*power[i]);
    }
    fclose(outFile);
  }
  MPI_Barrier(MPI_COMM_WORLD);

#if !PENCIL
  if (NumPEs > 1) {
    my_fftw_complex *data = (my_fftw_complex *)malloc(My_Ng*sizeof(my_fftw_complex));
    Parallel->d().redistribute_3_to_1((const complex_t *) rho, (complex_t *) data);
    for (i=0; i<My_Ng; ++i){rho[i] = data[i];}
    free(data);
  }
#endif
#endif
// end P(k) test


  
  return;
}


void Initializer::output(integer axis, real* pos, real* vel){ 
// 0=x, 1=y, 2=z
  const integer ngrid=DataBase->ngrid;
  long i, j;
  const char *fname;
  FILE *pfile;
  std::ofstream OutFile;
  MPI_Status status;
  MPI_Datatype MPI_FFTW_COMPLEX;
  
  // David's stuff:
#if !PENCIL
  if (NumPEs > 1) {
    my_fftw_complex *data = (my_fftw_complex *)malloc(My_Ng*sizeof(my_fftw_complex));
    Parallel->d().redistribute_1_to_3((const complex_t *) rho, (complex_t *) data);
    for (i=0; i<My_Ng; ++i){rho[i] = data[i];}
    free(data);
  }
#endif
  // End of David's stuff.

  if (DataBase->PrintFlag == 0){
    // No printing, particles are redirected somewhere
    for (i=0; i<My_Ng; ++i){
      pos[i] = (real)rho[i].re;
      vel[i] = (real)rho[i].im;
      if (pos[i] < 0.0) {pos[i] += (real)ngrid;}
      else if (pos[i] >= ngrid) {pos[i] -= (real)ngrid;}
    }
  }
  else if (DataBase->PrintFlag == 1){
    // Serial print through the master processor into an ASCII file
    CreateMPI_FFTW_COMPLEX(&MPI_FFTW_COMPLEX);
    if (MyPE == MasterPE){
      fname = "part.dat";
      OutFile.open(fname, std::ios::out | std::ios::app);
      for (i=0; i<My_Ng; ++i) 
      {OutFile << (float)rho[i].re << " " << (float)rho[i].im << std::endl;}
      for (j=0; j<NumPEs; ++j){  // Get data from processor j:
	if (j == MasterPE) continue;
	MPI_Recv(rho, My_Ng, MPI_FFTW_COMPLEX, j, 101, MPI_COMM_WORLD, &status);
	for (i=0; i<My_Ng; ++i) 
	  OutFile << (float)rho[i].re << " " << (float)rho[i].im << std::endl;
      }
      OutFile.close();
    }
    else{
      MPI_Send(rho, My_Ng, MPI_FFTW_COMPLEX, MasterPE, 101, MPI_COMM_WORLD);
    }
  }
  else if (DataBase->PrintFlag == 2){
    // Serial print through the master processor into a binary file
    CreateMPI_FFTW_COMPLEX(&MPI_FFTW_COMPLEX);
    if (MyPE == MasterPE) {
      fname = "part.bin";
      pfile = fopen(fname, "aw");
      fwrite(rho, sizeof(rho[0]), My_Ng, pfile);
      for (j=0; j<NumPEs; ++j){  // Get data from processor j:
	if (j == MasterPE) continue;
	MPI_Recv(rho, My_Ng, MPI_FFTW_COMPLEX, j, 101, MPI_COMM_WORLD, &status);
	for (i=0; i<My_Ng; ++i){fwrite(rho, sizeof(rho[0]), My_Ng, pfile);}
      }
      fclose(pfile);
    }
    else{
      MPI_Send(rho, My_Ng, MPI_FFTW_COMPLEX, MasterPE, 101, MPI_COMM_WORLD);
    }
  }
  else if (DataBase->PrintFlag == 3){
    // Parallel print -- each processor dumps a (binary) file
    
  }
  
  return;
}


void Initializer::initParticles(real* pos_x, real* pos_y, real* pos_z, 
				real* vel_x, real* vel_y, real* vel_z, 
				Basedata& bdata, const char *tfName, int useWN) {
  integer i;
  real d_z, ddot;
  double t1, t2;
  //void (*indens)();  /* pointer to indens routine */
  void (Initializer::*indens)();  /* pointer to indens routine */
  
#if PENCIL
  integer dim = 2;
#else
  integer dim = 1;
#endif
  
  // Preliminaries:
  StartMonitor();
  OnClock("Initialization");
  Parallel->InitMPI(&t1);
  MyPE = Parallel->GetMyPE();
  MasterPE = Parallel->GetMasterPE();
  NumPEs = Parallel->GetNumPEs();
  DataBase->GetParameters(bdata, *Parallel, dim);
  Cosmology->SetParameters(*DataBase, tfName);
  Parallel->DomainDecompose(DataBase->dim, DataBase->ngrid);
  if (NumPEs == 1){ indens = &Initializer::indens_single; }
  else if(useWN == 1) {
    indens = &Initializer::indens_whitenoise;
    if(MyPE == MasterPE)
      printf("Using white noise initializer\n");
  }
  else {
    if (DataBase->dim == 1) { indens = &Initializer::indens1D; }
    else if (DataBase->dim == 2) { indens = &Initializer::indens2D; }
    else if (DataBase->dim == 3) { indens = &Initializer::indens3D; }
    else {Parallel->ParallelError("Main: Decomposition has to be 1-3", "stop");}
  }
  
  // Find what part of the global grid i got:
  nx1 = Parallel->GetNx1();
  nx2 = Parallel->GetNx2();
  ny1 = Parallel->GetNy1();
  ny2 = Parallel->GetNy2();
  nz1 = Parallel->GetNz1();
  nz2 = Parallel->GetNz2();
  ngx = nx2-nx1+1;
  ngy = ny2-ny1+1;
  ngz = nz2-nz1+1;
  if (useWN == 0 && DataBase->dim==3 && (ngz != ngy || ngx != ngz)) 
    Parallel->AbortMPI("Main: 3D Domain incorrectly decomposed");
  if (useWN == 0 && DataBase->dim==2 && (ngx != ngy || ngz != DataBase->ngrid))
    Parallel->AbortMPI("Main: 2D Domain incorrectly decomposed");
  if (useWN == 0 && DataBase->dim==1 && (ngz != ngy || ngz != DataBase->ngrid))
    Parallel->AbortMPI("Main: 1D Domain incorrectly decomposed");
  My_Ng = ngx*ngy*ngz;
  
  // Allocate arays:
  Pk   = (real *)malloc(My_Ng*sizeof(real));
  rho  = (my_fftw_complex *)malloc(My_Ng*sizeof(my_fftw_complex));
  
  // Initialize FFTW:
#ifdef FFTW2
  const integer ngrid=DataBase.ngrid;
  plan_b = fftw3d_mpi_create_plan(MPI_COMM_WORLD, ngrid, ngrid, ngrid, 
				  FFTW_BACKWARD, FFTW_ESTIMATE);
  plan_f = fftw3d_mpi_create_plan(MPI_COMM_WORLD, ngrid, ngrid, ngrid, 
				  FFTW_FORWARD, FFTW_ESTIMATE);
  fftwnd_mpi_local_sizes(plan, &lnx, &lxs, &lnyt, &lyst, &lsize);
  if (lnx != ngx){
    std::cout << "FFTW has different decomposition!";
    std::cout << "  Mine: " << ngx << "  FFTW's: " << lnx << std::endl;
    Parallel->AbortMPI(" ");
  }
#endif
#ifdef FFTW3
  const integer ngrid=DataBase->ngrid;
#if PENCIL
  buf1  = (my_fftw_complex *)malloc(My_Ng*sizeof(my_fftw_complex));
  buf3  = (my_fftw_complex *)malloc(My_Ng*sizeof(my_fftw_complex));
  
  // from Solver::initialize
  Distribution &d = Parallel->d();
  int n;

  //forward, r to k
  n = d.local_ng_2d_x(0);
  plan_f_x = fftw_plan_many_dft(1, &n, d.local_ng_2d_x(1)*d.local_ng_2d_x(2), 
				FFTW_ADDR(buf1), NULL,
				1, n, FFTW_ADDR(buf3), NULL, 1, n, 
				FFTW_FORWARD, 0);
  n = d.local_ng_2d_y(1);
  plan_f_y = fftw_plan_many_dft(1, &n, d.local_ng_2d_y(0)*d.local_ng_2d_y(2), 
				FFTW_ADDR(buf3), NULL,
				1, n, FFTW_ADDR(buf1), NULL, 1, n, 
				FFTW_FORWARD, 0);
  n = d.local_ng_2d_z(2);
  plan_f_z = fftw_plan_many_dft(1, &n, d.local_ng_2d_z(1)*d.local_ng_2d_z(0), 
				FFTW_ADDR(buf1), NULL,
				1, n, FFTW_ADDR(rho), NULL, 1, n, 
				FFTW_FORWARD, 0);

  //backward, k to r
  n = d.local_ng_2d_z(2);
  plan_b_z = fftw_plan_many_dft(1, &n, d.local_ng_2d_z(1)*d.local_ng_2d_z(0), 
				FFTW_ADDR(rho), NULL,
				1, n, FFTW_ADDR(buf3), NULL, 1, n, 
				FFTW_BACKWARD, 0);
  n = d.local_ng_2d_y(1);
  plan_b_y = fftw_plan_many_dft(1, &n, d.local_ng_2d_y(0)*d.local_ng_2d_y(2), 
				FFTW_ADDR(buf3), NULL,
				1, n, FFTW_ADDR(buf1), NULL, 1, n, 
				FFTW_BACKWARD, 0);
  n = d.local_ng_2d_x(0);
  plan_b_x = fftw_plan_many_dft(1, &n, d.local_ng_2d_x(1)*d.local_ng_2d_x(2), 
				FFTW_ADDR(buf1), NULL,
				1, n, FFTW_ADDR(buf3), NULL, 1, n, 
				FFTW_BACKWARD, 0);

#else
  fftw_mpi_init();
  plan_b = fftw_mpi_plan_dft_3d(ngrid, ngrid, ngrid, (fftw_complex *)rho, 
				(fftw_complex *)rho, MPI_COMM_WORLD, 
				FFTW_BACKWARD, FFTW_ESTIMATE);
  plan_f = fftw_mpi_plan_dft_3d(ngrid, ngrid, ngrid, (fftw_complex *)rho, 
				(fftw_complex *)rho, MPI_COMM_WORLD, 
				FFTW_FORWARD, FFTW_ESTIMATE);
#endif
#endif
  OffClock("Initialization");
  
  // Set P(k) array:
  OnClock("Creating Pk(x,y,z)");
  initspec();
  OffClock("Creating Pk(x,y,z)");
  
  // Growth factor for the initial redshift:
  Cosmology->GrowthFactor(DataBase->z_in, &d_z, &ddot);
  if (MyPE == MasterPE){
    printf("redshift: %f; growth factor = %f; ", DataBase->z_in, d_z);
    printf("derivative = %f\n\n", ddot);
  }
  
#ifdef TESTING
  // Check if density field is real:
  OnClock("Tests");

  //I seriously don't know if this is a good idea
  (this->*indens)();

  test_reality();

  OffClock("Tests");
#endif

  // Set particles:
  for (i=0; i<Ndim; ++i){
    OnClock("Creating rho(x,y,z)");
    //I seriously don't know if this is a good idea
    (this->*indens)();          // Set density in k-space:
    OffClock("Creating rho(x,y,z)");

    OnClock("Poisson solve");
    solve_gravity(i);  // Calculate i-th component of the force
    OffClock("Poisson solve");

    OnClock("Particle move");
    set_particles(DataBase->z_in, d_z, ddot, i); // Zeldovich move
    OffClock("Particle move");
    
    OnClock("Output");
    if (i == 0) {output(i, pos_x, vel_x);}   // Put particle data
    else if (i == 1) {output(i, pos_y, vel_y);}
    else {output(i, pos_z, vel_z);} 
    OffClock("Output");
  }
  
  // Clean all the allocations:
  if(Pk) free(Pk);
  if(rho) free(rho);
  if(buf1) free(buf1);
  if(buf3) free(buf3);

  //Seriously? moved to destructor
  //Cosmology.~CosmoClass();
  //DataBase.~GlobalStuff();
  //Parallel.~ParallelTools();

#ifdef FFTW2
  fftwnd_mpi_destroy_plan(plan_f);
  fftwnd_mpi_destroy_plan(plan_b);
#endif
#ifdef FFTW3
#if PENCIL
  fftw_destroy_plan(plan_f_x);
  fftw_destroy_plan(plan_f_y);
  fftw_destroy_plan(plan_f_z);
  fftw_destroy_plan(plan_b_x);
  fftw_destroy_plan(plan_b_y);
  fftw_destroy_plan(plan_b_z);
#else
  fftw_destroy_plan(plan_f);
  fftw_destroy_plan(plan_b);
#endif
#endif
  
  MPI_Barrier(MPI_COMM_WORLD);
  if (MyPE == MasterPE) PrintClockSummary(stdout);
  Parallel->FinalMPI(&t2);	
  
  return;
}


void init_particles(::real* pos_x, ::real* pos_y, ::real* pos_z, 
		    ::real* vel_x, ::real* vel_y, ::real* vel_z, 
		    Basedata& bdata, const char *tfName, int useWN) {
  Initializer *init = new Initializer();
  init->initParticles(pos_x, pos_y, pos_z, vel_x, vel_y, vel_z, bdata, tfName, useWN);
  delete init;

  return;
}


#undef Ndim
#undef pi
