//-*-C++-*-//////////////////////////////////////////////////////////////////
// $Id: Particles.h,v 1.5 2010/08/06 20:37:28 pope Exp $
/////////////////////////////////////////////////////////////////////////////

#ifndef PARTICLES_H
#define PARTICLES_H

#include <cassert>
#include <cmath>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <sys/time.h>

#include <stdlib.h>
#include <string.h>

#ifdef _OPENMP
#include <omp.h>
#endif

#include "Definition.h"
#include "mc3types.h"
#include "Basedata.h"
#include "array_index.h"
#include "Domain.h"
#include "TimeStepper.h"
#include "SimpleTimings.h"
#include "MC3Options.h"
#include "ForceLaw.h"
#include "BHForceTree.h"
#include "RestartIO.h"
#include "Partition.h"



#define MEM_ALIGN 64
#define N_POSVEL_T 8



class CMLite
{
public:

  CMLite(int ngx, int ngy, int ngz);
  ~CMLite();

  int ng[DIMENSION];
  int Ng;
  int *indxlo;
  int *indxhi;
};



/*
class FGrid {
public:
  FGrid();
  ~FGrid(){};

  float fgor(float r);
  void fgor_r2_interp(int nInterp, float **r2, float **f);

  float m_b, m_c, m_d, m_e, m_f, m_g, m_h, m_l, m_rmax;
};



class FGridEval
{
public:
  FGridEval() {};
  virtual ~FGridEval() {};
  virtual float eval(float) {};
};



class FGridEvalFit : public FGridEval
{
public:
  FGridEvalFit(FGrid *fg);
  ~FGridEvalFit() {};
  float eval(float);

  FGrid *m_fg;
};



class FGridEvalInterp : public FGridEval
{
public:
  FGridEvalInterp(FGrid *fg, int nInterp);
  ~FGridEvalInterp();
  float eval(float);

  float *m_r2;
  float *m_f;
  float m_r2min;
  float m_r2max;
  float m_dr2;
  int m_nInterp;  
};
*/



class Particles
{
public:

  Particles(const Basedata &, MC3Options & options);
  ~Particles();

  GRID_T *field() {return m_field; }
  int Np_local_total() { return m_Np_local_total; }
  int Np_local_alive();

  vector<unsigned int> *aliveIndices();

  void copyAllFromVectors(vector<POSVEL_T> **xx,
			  vector<POSVEL_T> **yy,
			  vector<POSVEL_T> **zz,
			  vector<POSVEL_T> **vx,
			  vector<POSVEL_T> **vy,
			  vector<POSVEL_T> **vz,
			  vector<POSVEL_T> **phi,
			  vector<ID_T> **id,
			  vector<MASK_T> **mask,
			  float anow,
			  int deleteVectorsQ=0);

  void copyAliveIntoVectors(vector<POSVEL_T> *xx,
			    vector<POSVEL_T> *yy,
			    vector<POSVEL_T> *zz,
			    vector<POSVEL_T> *vx,
			    vector<POSVEL_T> *vy,
			    vector<POSVEL_T> *vz,
			    vector<POSVEL_T> *phi,
			    vector<ID_T> *id,
			    vector<MASK_T> *mask,
			    float anow,
			    int npr=0,
			    int dropParticlesQ=0);


  void shoveParticles() {return;};
  void grabParticles() {return;};
  void dropParticles();
  void resortParticles();
  void allocField();

  void writeAliveHCosmo( const char *outName, float anow );
  void writeRestart(const char *outName);
  void readRestart(const char *inName);

  void map1(float pp, float tau, float adot);
  void cic();
  void inverse_cic(float tau, float fscal, int comp);
  void inverse_cic_potential();

  FGrid *f_grid_over_r;
  void loadInterp(int n, float *x, float *y);
  void forceInterp(uint32_t nInterp);

  void subCycle(TimeStepper *gts);
  void map2(TimeStepper *ts, TS_FLOAT stepFraction=1.0);

  void buildChainingMesh();
  void subCycleCM(TimeStepper *gts);
  void map2CM(TimeStepper *ts, TS_FLOAT stepFraction=1.0);

  void writeRawAscii(const char * fname );

private:

  Particles();
  Particles( const Particles& );
  Particles& operator = ( const Particles& );

  void coords_global2local(float anow);
  void coords_local2global(float anow);

  void allocParticles(int Np);
  void updatePointers();
  void updatePointers2();

  int array_index(int xx, int yy, int zz, 
		  int ng[], int lo[], int hi[],
		  int safe);

  POSVEL_T *m_xArr;
  POSVEL_T *m_yArr;
  POSVEL_T *m_zArr;
  POSVEL_T *m_vxArr;
  POSVEL_T *m_vyArr;
  POSVEL_T *m_vzArr;
  POSVEL_T *m_phiArr;
  POSVEL_T *m_massArr;
  POSVEL_T *m_pvData[N_POSVEL_T];
  ID_T *m_idArr;
  MASK_T *m_maskArr;

  GRID_T *m_field;

  int m_Np_local_total;
  int m_Np_last;

  float m_alpha;
  float m_gpscal;

  int m_coords_localQ;

  int m_nsub;
  float m_edge;
  float m_rsm;
  float m_fsrrmax;
  float m_cmsize;
  float m_openAngle;

  CMLite *m_cm;

  FGrid *m_fg;
  FGridEval *m_fgore;

  ForceLaw *m_fl;

  int m_skipStreamQ;
  int m_skipKickSRQ;
  int m_useFastTreeEval;
  int m_useSymmTree;
  int m_useRCBTree;
  int m_rcbTreeExtraLevels;
  int m_rcbTreePPN;
  int m_rcbTreeTaskPartMin;

  int m_mpiio;
};



#endif

