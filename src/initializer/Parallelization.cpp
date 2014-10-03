/*  
   Initializer:
   Parallelization.cpp

      Contains routines which handle mpi parallelization (see 
         Paralelization.h for the interface). Important ones 
         are:
               InitMPI(double *w_time) 
                  which should be called at the start of the program; it 
         starts mpi, sets process number, master process number, 
         and number of processes. It also returns wall clock time 
         via w_time variable. Similar to that, 
               FinalMPI(double *wtime)
         closes mpi session, and returns end time. 

         The largest and most important routine is
               DomainDecompose(integer dim, integer ngrid)
         where dim defines dimension of decomposition:
               dim = 1 is slab decomposition (x axis is sliced)
               dim = 2 is rod decomposition (rods go along z axis)
               dim = 3 is cube decomposition
         and ngrid is number of grid zones per dimension. 
         Routine sets indicies for the arrays, but also determines 
         which processor hold mirror-symmetric data. This is needed 
         for enforcing reality condition of the density field.

         For example, full array: 
               rho(0:ngrid-1, 0:ngrid-1, 0:ngrid-1)
         is divided such that each processor holds:
               rho(nx1:nx2, ny1:ny2, nz1:nz2)

         Serial job (running on 1 processor) is also handled, and 
                  dim parameters is obviously not important there. 

                           Zarija Lukic, February 2009
                                 zarija@lanl.gov
*/

#include <mpi.h>
#include <iostream>
#include <math.h>
#include "Parallelization.h"


/* Routine should be used for errors which *ALL* processes encounter */
void ParallelTools::ParallelError(const char *message, const char *stop){
  int MyProc, MasterProc, rc = 1;
  MyProc = GetMyPE();
  MasterProc = GetMasterPE();
  if (MyProc == MasterProc)
    std::cout << message << std::endl << std::flush;
  MPI_Barrier(MPI_COMM_WORLD);
  if (stop == "stop" or stop == "Stop" or stop == "STOP"){
    std::cout << std::endl;
    MPI_Abort(MPI_COMM_WORLD, rc);
  }
  return;
}


/* Handling errors which *SOME* processes encounter */
void ParallelTools::AbortMPI(const char *message){
  int rc = 1;
  std::cout << message << std::endl << std::endl << std::flush;
  MPI_Abort(MPI_COMM_WORLD, rc);
  return;
}


void ParallelTools::InitMPI(double *w_time){
  char **argv;
  int argc, rc;
  int MasterProc, MyProc, NumProcs;
//	rc = MPI_Init(&argc, &argv);
//	if (rc != MPI_SUCCESS)
//		ParallelError("Could not start MPI.", "stop");
  *w_time = MPI_Wtime();
  MPI_Comm_size(MPI_COMM_WORLD, &NumProcs);
  MPI_Comm_rank(MPI_COMM_WORLD, &MyProc);
  MasterProc = 0;
  SetMyPE(MyProc);
  SetMasterPE(MasterProc);
  SetNumPEs(NumProcs);
  return;
}


void ParallelTools::FinalMPI(double *w_time){
  *w_time = MPI_Wtime();
//	MPI_Finalize();
  return;
}


ParallelTools::ParallelTools() {
  dist = NULL;
}


ParallelTools::~ParallelTools() {
  if(dist)
    delete dist;
}


void ParallelTools::DomainDecompose(integer dim, integer ngrid){
  int MyProc, NumProcs, MasterProc, MirrorProc, Mirror1Proc, MirrorFlag;
  integer ndiv;
  integer nx1, nx2, ny1, ny2, nz1, nz2;
  
  int nga[3] = { ngrid, ngrid, ngrid };

  //dist = new Distribution();
  //dist->initialize(MPI_COMM_WORLD, nga, nga);
  dist = new Distribution(MPI_COMM_WORLD, nga);

#if !PENCIL
  MPI_Comm_rank(dist->cart_1d(), &MyProc);
#else
  MPI_Comm_rank(dist->cart_2d_z(), &MyProc);
#endif
  NumProcs = GetNumPEs();
  MasterProc = GetMasterPE();
  if (NumProcs == 1) {         // Serial job:
    std::cout << std::endl;
    std::cout << "Running as a serial job." << std::endl;
    std::cout << ngrid << "^3 grid" << std::endl;
    std::cout << "No decomposition" << std::flush;
    SetNx1(0);
    SetNx2(ngrid-1);
    SetNy1(0);
    SetNy2(ngrid-1);
    SetNz1(0);
    SetNz2(ngrid-1);
#if !PENCIL
    SetMirrorPE(0);
#else
    for (int i = 0; i < 6; ++i) {
      SetMirrorPE(i, 0);
    }
#endif
    SetMirrorFlag(0);
    std::cout << "...............done" << std::endl << std::flush;
    return;
  }
  
  switch (dim) {
#if !PENCIL
    case 1:
      if (MyProc == MasterProc){
	std::cout << std::endl;
	std::cout << "Initializer will use " << NumProcs << " processors." << std::endl;
	std::cout << ngrid << "^3 grid" << std::endl;
	std::cout << "Decomposing into slabs" << std::flush;
      }
      SetNumPEs(NumProcs);
      nx1 = dist->self_1d(0)*dist->local_ng_1d(0);
      nx2 = nx1 + dist->local_ng_1d(0) - 1;
      ny1 = dist->self_1d(1)*dist->local_ng_1d(1);
      ny2 = ny1 + dist->local_ng_1d(1) - 1;
      nz1 = dist->self_1d(2)*dist->local_ng_1d(2);
      nz2 = nz1 + dist->local_ng_1d(2) - 1;
      
      SetMyPE(dist->self());
      MyProc = GetMyPE();
      break;
#else
    case 2:
      if (MyProc == MasterProc){
	std::cout << std::endl;
	std::cout << "Initializer will use " << NumProcs << " processors." << std::endl;
	std::cout << ngrid << "^3 grid" << std::endl;
	std::cout << "Decomposing into pencils" << std::flush;
      }
      SetNumPEs(NumProcs);
      nx1 = dist->self_2d_z(0)*dist->local_ng_2d_z(0);
      nx2 = nx1 + dist->local_ng_2d_z(0) - 1;
      ny1 = dist->self_2d_z(1)*dist->local_ng_2d_z(1);
      ny2 = ny1 + dist->local_ng_2d_z(1) - 1;
      nz1 = dist->self_2d_z(2)*dist->local_ng_2d_z(2);
      nz2 = nz1 + dist->local_ng_2d_z(2) - 1;
      
      MPI_Comm_rank(dist->cart_2d_z(), &MyProc);
      SetMyPE(MyProc);
      break;
#endif
    default:
      ParallelError("DomainDecompose: dim has to be 1 or 2 and must match PENCIL", "stop");
      break;
  }
  /* Set domain */
  SetNx1(nx1);
  SetNx2(nx2);
  SetNy1(ny1);
  SetNy2(ny2);
  SetNz1(nz1);
  SetNz2(nz2);
  
#if !PENCIL
  if (dist->self() < NumProcs/2) {
    MirrorFlag = 0;
  }
  else {
    MirrorFlag = 1;
  }
  
  int xmcoord = NumProcs - 1 - dist->self();
  MPI_Cart_rank(dist->cart_1d(), &xmcoord, &MirrorProc);
  
  ++xmcoord;
  MPI_Cart_rank(dist->cart_1d(), &xmcoord, &Mirror1Proc);
  
  SetMirrorPE(MirrorProc);
  SetMirror1PE(Mirror1Proc);
  SetMirrorFlag(MirrorFlag);
#else
  /* Set mirror symmetric processor data */
  /* Force mirror symmetry across the X axis, because F[k1,k2,k3] = F*[N-k1,N-k2,N-k3],
     if k1 == N/2 then F[N/2,k2,k3] == F*[N/2,k2,k3] and if k1 > N/2 ... */
  
  /* 
     So with pencils, this is slightly more complicated than with slabs: 
     With slabs you only had to know which rank held any given x slice. 
     Because each slice was mapped to only one node, knowing the conj. slice 
     have you the conj. node.
     Here, we still symmetry-reduce about the x axis, but each x slice is 
     mapped to multiple nodes. 
     We need to know all of the ranks which hold conj. values to our x slices.
     If we assume that the pencils have square bases, then our conj. values 
     should be on only 4 other ranks. 
  */
  
  int c[6][3];
  for (int j = 0; j < 4; ++j)
    for (int i = 0; i < 3; ++i) {
      c[j][i] = dist->nproc_2d_z(i) - 1 - dist->self_2d_z(i);
    }
  for (int j = 4; j < 6; ++j) {
    c[j][0] = dist->self_2d_z(0);
    for (int i = 1; i < 3; ++i) {
      c[j][i] = dist->nproc_2d_z(i) - 1 - dist->self_2d_z(i);
    }
  }
  
  c[1][0] += 1;
  c[2][1] += 1;
  c[3][0] += 1;
  c[3][1] += 1;
  
  c[5][1] += 1;
  
  for (int j = 0; j < 6; ++j) {
    for (int i = 0; i < 3; ++i)
      if (c[j][i] == dist->nproc_2d_z(i)) c[j][i] = 0;
    
    SetMirrorPE(j, dist->rank_2d_z(c[j]));
  }
  
  if (dist->self_2d_z(0) < dist->nproc_2d_z(0)/2) {
    SetMirrorFlag(0);
  }
  else {
    SetMirrorFlag(1);
  }
#endif
  
  if (MyProc == MasterProc)
    std::cout << "...............done" << std::endl << std::flush;
  
  return;
}
