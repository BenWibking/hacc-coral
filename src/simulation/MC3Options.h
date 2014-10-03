#ifndef MC3OPTIONS_H
#define MC3OPTIONS_H

#include <unistd.h>
#include <stdlib.h>

#include <rru_mpi.h>

#include "BasicDefinition.h"
#include "mc3types.h"
#include "mc3defines.h"

#include "Partition.h"

class MC3Options {
 public:
  MC3Options(int argc, char *argv[]);
  ~MC3Options() {}

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
  int interpQ() { return m_interpQ; }
  int nInterp() { return m_nInterp; }
  int polyQ() { return m_polyQ; }
  string exeName() { return m_exeName; }
  int bladeQ() { return m_bladeQ; }
  int cmQ() { return m_cmQ; }
  int everyCMQ() { return m_everyCMQ; }
  int skipStreamQ() { return m_skipStreamQ; }
  int skipKickLRQ() { return m_skipKickLRQ; }
  int skipKickSRQ() { return m_skipKickSRQ; }
  int dontDropMemory() { return m_dontDropMemory; }
  int dontUseBigchunk() { return m_dontUseBigchunk; }
  int useFastTreeEval() { return m_useFastTreeEval; }
  int useRCBTree() { return m_useRCBTree; }
  int rcbTreeExtraLevels() { return m_rcbTreeExtraLevels; }
  int rcbTreePPN() { return m_rcbTreePPN; }
  int rcbTreeTaskPartMin() { return m_rcbTreeTaskPartMin; }
  int whiteNoiseInit() { return m_whiteNoiseInit; }
  int mpiio() { return m_mpiio; }
  int topologyQ() { return m_topologyQ; }
  string topologyString() { return m_topologyString; }

 private:

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
  int m_interpQ;
  int m_nInterp;
  int m_polyQ;
  string m_exeName;
  int m_bladeQ;
  int m_cmQ;
  int m_everyCMQ;
  int m_skipStreamQ;
  int m_skipKickLRQ;
  int m_skipKickSRQ;
  int m_dontDropMemory;
  int m_dontUseBigchunk;
  int m_useFastTreeEval;
  int m_useRCBTree;
  int m_rcbTreeExtraLevels;
  int m_rcbTreePPN;
  int m_rcbTreeTaskPartMin;
  int m_whiteNoiseInit;
  int m_mpiio;
  int m_topologyQ;
  string m_topologyString;
};

#endif
