#ifndef MC3EXTRAS_H
#define MC3EXTRAS_H

#include <unistd.h>
#include <stdlib.h>

#include <rru_mpi.h>

#include "BasicDefinition.h"
#include "Partition.h"
#include "ParticleExchange.h"
#include "CosmoHaloFinderP.h"
#include "FOFHaloProperties.h"

#include "mc3types.h"
#include "mc3defines.h"

#include "Basedata.h"
#include "Halodata.h"
#include "Domain.h"
#include "lightcone.h"
#include "skewers.h"
#include "Particles.h"
#include "MC3Options.h"



class MC3Extras {
 public:
  MC3Extras(MC3Options & options, Basedata & indat);
  ~MC3Extras() {}

  void setStep(int step, float aa);
  void particleExtras(Particles & particles,
		      Basedata & indat,
		      string outBase,
		      vector<int> *Nplav);

  int gridQ() { return m_gridQ; }
  int restartDumpQ() { return m_restartDumpQ; }
  string restartDumpName() { return m_restartDumpName; }
  int pkDumpQ() { return m_pkDumpQ; }
  string pkDumpName() { return m_pkDumpName; }
  int aliveDumpQ() { return m_aliveDumpQ; }
  string aliveDumpName() { return m_aliveDumpName; }
  int refreshQ() { return m_refreshQ; }
  string refreshName() { return m_refreshName; }
  int analysisQ() { return m_analysisQ; }
  string analysisdatName() { return m_analysisdatName; }
  int staticQ() { return m_staticQ; }
  string staticDumpName() { return m_staticDumpName; }
  int lightconeQ() { return m_lightconeQ; }
  string LCUpdateName() { return m_LCUpdateName; }
  int skewerQ() { return m_skewerQ; }
  int haloQ() { return m_haloQ; }
  bool initialAlltoallQ() { return m_initialAlltoallQ; }
  int mpiio() { return m_mpiio; }

  vector<int> *restartDumps() { return m_restartDumps; }
  vector<int> *pkDumps() { return m_pkDumps; }
  vector<int> *aliveDumps() { return m_aliveDumps; }
  vector<int> *refreshSteps() { return m_refreshSteps; }
  vector<int> *staticDumps() { return m_staticDumps; }
  vector<int> *LCUpdates() { return m_LCUpdates; }

  Halodata *analysisdat() { return m_analysisdat; }

  ParallelSkewerSet *zskewers() { return &m_zskewers; }
  ConvergentSkewerSet *lcskewers() { return &m_lcskewers; }

  LightCone *lc() { return m_lc; }
  float aa_last() { return m_aa_last; }
  int lc_started() { return m_lc_started; }

  int staticStep() { return m_staticStep; }
  int lcStep() { return m_lcStep; }
  int aliveStep() { return m_aliveStep; }
  int restartStep() { return m_restartStep; }
  int pkStep() { return m_pkStep; }
  int refreshStep() { return m_refreshStep; }

  int extrasStep() { return m_extrasStep; }
  int fftfStep() { return m_fftfStep; }
  int fftbpotStep() { return m_fftbpotStep; }
  int particleStep() { return m_particleStep; }

 private:

  void initializeSkewers(Basedata & indat);
  vector<int>* readInts(string inName);
  int intInVec(int n, vector<int> *iv);
  string create_outName(string outBase, int rank);
  void allocReserveVectors(int Np);
  void reserveVectors(int Np);
  void deleteVectors();
  int calcNReserve(vector<int> *Nplav);
  void findHalos(Basedata & indat, string outBase);
  void refreshParticles(Particles & particles);
  /*
  void staticSkewers(Particles & particles, string outBase);
  void lightconeSkewers(Particles & particles);
  */
  void staticSkewers(string outBase);
  void lightconeSkewers();

  int m_gridQ;
  int m_restartDumpQ;
  string m_restartDumpName;
  int m_pkDumpQ;
  string m_pkDumpName;
  int m_aliveDumpQ;
  string m_aliveDumpName;
  int m_refreshQ;
  string m_refreshName;
  int m_analysisQ;
  string m_analysisdatName;
  int m_staticQ;
  string m_staticDumpName;
  int m_lightconeQ;
  string m_LCUpdateName;
  int m_skewerQ;
  int m_haloQ;
  bool m_initialAlltoallQ;
  int m_mpiio;

  vector<int> *m_restartDumps;
  vector<int> *m_pkDumps;
  vector<int> *m_aliveDumps;
  vector<int> *m_refreshSteps;
  vector<int> *m_staticDumps;
  vector<int> *m_LCUpdates;

  Halodata *m_analysisdat;

  ParallelSkewerSet m_zskewers;
  ConvergentSkewerSet m_lcskewers;
  LightCone *m_lc;
  float m_aa_last;
  int m_lc_started;

  int m_staticStep;
  int m_lcStep;
  int m_aliveStep;
  int m_restartStep;
  int m_pkStep;
  int m_refreshStep;

  int m_extrasStep;
  int m_fftfStep;
  int m_fftbpotStep;
  int m_particleStep;

  vector<POSVEL_T> *m_xx;
  vector<POSVEL_T> *m_yy;
  vector<POSVEL_T> *m_zz;
  vector<POSVEL_T> *m_vx;
  vector<POSVEL_T> *m_vy;
  vector<POSVEL_T> *m_vz;
  vector<POSVEL_T> *m_phi;
  vector<ID_T> *m_id;
  vector<MASK_T> *m_mask;

  vector<POSVEL_T> *m_mass;
  vector<STATUS_T> *m_status;

  int m_step;
  float m_aa;
};

#endif
