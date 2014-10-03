/*
   Initializer:
   Parallelization.h

      Defines class ParallelTools. 
         The idea is that this class holds parallelization data 
         needed by almost all other routines. 
                  Only class member functions can set the data; other routines 
         that use this class can acces the data via public Get 
         interface.

         There are 5 methods; two are used to set globally visible 
         data which can roughly be divided into 2 categories: 
         mpi stuff (done by InitMPI), and data decomposition stuff 
         (done by DomainDecompose -- see Parallel.cpp for explanation 
         on how it is done). 
         Other three methods provide some basic interface for parallel 
         runs. In particular, for error (or warning) messages one 
         should use 
               ParallelError("Message", "flag")
         if flag is stop, the code will be aborted; if flag is anything 
         else nothing more (but printing) will be done. 
         
                        Zarija Lukic, February 2009
                              zarija@lanl.gov
*/

#ifndef Parallelization_Header_Included
#define Parallelization_Header_Included

#include "TypesAndDefs.h"
#include "distribution.hpp"

class ParallelTools {
 public :
  /* Methods: */
  ParallelTools();
  ~ParallelTools();
  
  void InitMPI(double*);
  void FinalMPI(double*);
  void AbortMPI(const char*);
  void ParallelError(const char*, const char*);
  void DomainDecompose(integer, integer);
      
  /* Public interface to data: */
  /* MPI parameters */
  int GetMyPE() {return MyProc;}
  int GetNumPEs() {return NumProcs;}
  int GetMasterPE() {return MasterProc;}
#if !PENCIL
  int GetMirrorPE() {return MirrorProc;}
  int GetMirror1PE() {return Mirror1Proc;}
#else
  int GetMirrorPE(int i) {return MirrorProc[i];}
#endif
  int GetMirrorFlag() {return MirrorFlag;}	
  
  /* Domain decomposition -- particles/grid */
  integer GetNx1() {return nx1;}
  integer GetNx2() {return nx2;}
  integer GetNy1() {return ny1;}
  integer GetNy2() {return ny2;}
  integer GetNz1() {return nz1;}
  integer GetNz2() {return nz2;}
  
 private:
  int NumProcs, MyProc, MasterProc;
#if !PENCIL
  int MirrorProc, Mirror1Proc;
#else
  int MirrorProc[6];
#endif
  int MirrorFlag;
  integer nx1, nx2, ny1, ny2, nz1, nz2;
  
  void SetMyPE(int i) {MyProc = i;}
  void SetNumPEs(int i) {NumProcs = i;}
  void SetMasterPE(int i) {MasterProc = i;}
#if !PENCIL
  void SetMirrorPE(int i) {MirrorProc = i;}
  void SetMirror1PE(int i) {Mirror1Proc = i;}
#else
  void SetMirrorPE(int i, int j) {MirrorProc[i] = j;}
#endif
  void SetMirrorFlag(int i) {
    if (i == 0 or i == 1) {MirrorFlag = i;}
    else {AbortMPI("Mirror flag can only be 0 or 1!");}
  }
  
  void SetNx1(integer i) {nx1 = i;}
  void SetNx2(integer i) {nx2 = i;}
  void SetNy1(integer i) {ny1 = i;}
  void SetNy2(integer i) {ny2 = i;}
  void SetNz1(integer i) {nz1 = i;}
  void SetNz2(integer i) {nz2 = i;}
 protected:
  Distribution *dist;
  
 public:
  Distribution &d() { return *dist; }
};

#endif
