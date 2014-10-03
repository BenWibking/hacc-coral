//------------------------------------------------------------------------------
//
// parallel io class
//
// Tom Peterka
// Argonne National Laboratory
// 9700 S. Cass Ave.
// Argonne, IL 60439
// tpeterka@mcs.anl.gov
//
// (C) 2012 by Argonne National Laboratory.
// See COPYRIGHT in top-level directory.
//
//--------------------------------------------------------------------------

#ifndef _RESTART_IO
#define _RESTART_IO

#include <stddef.h>
#include <stdint.h>
#include <assert.h>
#include <stdlib.h>
#include <string.h>
#include <vector>

#include "mpi.h"

using namespace std;

// I/O modes
enum {
  IO_WRITE_RESTART,
  IO_READ_RESTART
};

// displacement types
#define OFST 0
#define ADDR 1

// typemap block for creating custom datatypes
struct map_block_t {
  MPI_Datatype base_type; // existing datatype used to create this one
  int disp_type; // diplacement is relative OFST or absolute ADDR
  int count;  // count of each element
  MPI_Aint disp; // displacement of each element in bytes
		 // OFSTs are from the start of the type and ADDRs are from 0x
};

class RestartIO {

 public:

  RestartIO(int mode, char *filename, MPI_Comm comm);
  ~RestartIO(){};

  // restart I/O
  void WriteRestart(int &num_particles, float *xx, float *yy, float *zz,
		    float *vx, float *vy, float *vz, float *phi, int64_t *pid,
		    uint16_t *mask);

  int ReadRestart(float *&xx, float *&yy, float *&zz,
		  float *&vx, float *&vy, float *&vz, float *&phi, int64_t *&pid,
		  uint16_t *&mask);

 private:

  void handle_mpio_error(int errcode, char *str);
  void WriteFooter(MPI_File fd, const int *blk_sizes,
			   int nb_out);
  void ReadFooter(MPI_File fd, int *&ftr, int &tb);
  void CreateRestartWriteDtype(int &num_particles, float *xx, float *yy, 
			       float *zz, float *vx, float *vy, float *vz, 
			       float *phi, int64_t *pid, uint16_t *mask,
			       MPI_Datatype *dtype);
  void CreateRestartReadDtype(int num_particles, float *xx, float *yy, 
			      float *zz, float *vx, float *vy, float *vz, 
			      float *phi, int64_t *pid, uint16_t *mask,
			      MPI_Datatype *dtype);
  void CreateRestartDtype(int num_particles, float *xx, float *yy, 
			  float *zz, float *vx, float *vy, float *vz, 
			  float *phi, int64_t *pid, uint16_t *mask,
			  MPI_Datatype *dtype);
  void CreateDtype(MPI_Aint addr, int num_map_blocks, map_block_t *typemap, 
		   MPI_Datatype *type);
  MPI_Aint mpi_addr(void *addr);

  MPI_Comm comm; // mpi communicator
  int rank; // MPI process rank
  int groupsize; // number of MPI processes
  MPI_File fd_in, fd_out; // input, output restart file descriptor
  bool restart_write_initialized; // whether restart writing was initialized
  bool restart_read_initialized; // whether restart reading was initialized

};

#endif
