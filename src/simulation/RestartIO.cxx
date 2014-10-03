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
// parallel netcdf file writing by Saba Sehrish, Seung Woo Son, Wei-keng Liao
// Northwestern University
//
// (C) 2012 by Argonne National Laboratory.
// See COPYRIGHT in top-level directory.
//
//--------------------------------------------------------------------------

#include "RestartIO.h"

//----------------------------------------------------------------------------
//
// Constructor
// opens the file for collective writing or reading of either analysis or
//  restart data, all processes need to call it
//
// mode: I/O mode
//   IO_WRITE_ANALYSIS, IO_READ_ANALYSIS, IO_WRITE_RESTART, IO_READ_RESTART
// filename: filename
// comm: MPI communicator
//
RestartIO::RestartIO(int mode, char *filename, MPI_Comm comm) {

  this->comm = comm;
  MPI_Comm_rank(comm, &rank);
  MPI_Comm_size(comm, &groupsize);
  restart_write_initialized = false;
  restart_read_initialized = false;

  if (mode == IO_WRITE_RESTART) {
    assert(MPI_File_open(comm, (char *)filename, 
			 MPI_MODE_WRONLY | MPI_MODE_CREATE,
			 MPI_INFO_NULL, &fd_out) == MPI_SUCCESS);
    MPI_File_set_size(fd_out, 0); // start with an empty file every time
    restart_write_initialized = true;
  }

  else if (mode == IO_READ_RESTART) {
    assert(MPI_File_open(comm, (char *)filename, MPI_MODE_RDONLY,
			 MPI_INFO_NULL, &fd_in) == MPI_SUCCESS);
    restart_read_initialized = true;
  }

  else 
    fprintf(stderr, "Unrecognized I/O option, nothing done.\n");

}
//----------------------------------------------------------------------------
//
// Writes checkpoint restart file
//
// num_particles: number of particles in this block
// xx, yy, zz: particle positions
// vx, vy, vz: particle velocities
// phi: particle potentials
// pid: particle IDs
// mask: unused
//
void RestartIO::WriteRestart(int &num_particles, float *xx, float *yy, float *zz,
			     float *vx, float *vy, float *vz, float *phi, int64_t *pid,
			     uint16_t *mask) {

  MPI_Datatype dtype;
  MPI_Status status;
  int errcode;

  // verify that writing was initialized
  assert(restart_write_initialized);

  // create data type
  CreateRestartWriteDtype(num_particles, xx, yy, zz, vx, vy, vz, phi, pid, 
			  mask, &dtype);

  // write block
  int size; // datatype size
  int scan_size; // exclusive scan of sizes before me
  MPI_Offset ofst = 0; // file pointer
  MPI_Type_size(dtype, &size);
  MPI_Exscan(&size, &scan_size, 1, MPI_INT, MPI_SUM, comm);
  if (rank > 0)
    ofst += scan_size;
  errcode = MPI_File_write_at_all(fd_out, ofst, MPI_BOTTOM, 1, dtype, &status);
  if (errcode != MPI_SUCCESS)
    handle_mpio_error(errcode, (char *)"MPI_File_write_all datatype");

  // finalize writing
  int *all_sizes;
  if (rank == 0)
    all_sizes = new int[groupsize];
  MPI_Gather(&size, 1, MPI_INT, all_sizes, 1, MPI_INT, 0, comm);
  if (rank == 0) { // only root writes the footer
    WriteFooter(fd_out, all_sizes, groupsize);
    delete[] all_sizes;
  }

  // cleanup
  MPI_File_close(&fd_out);
  MPI_Type_free(&dtype);

}
//----------------------------------------------------------------------------
//
// Reads checkpoint restart file
//
// num_particles: number of particles in this block
//
//  following arrays are allocated by this function, not the caller
// xx, yy, zz: particle positions
// vx, vy, vz: particle velocities
// phi: particle potentials
// pid: particle IDs
// mask: unused
//
// side effects: allocates above arrays
// returns: local number of particles
//
int RestartIO::ReadRestart(float *&xx, float *&yy, float *&zz,
			   float *&vx, float *&vy, float *&vz, float *&phi, 
			   int64_t *&pid, uint16_t *&mask) {

  int *ftr; // footer data, allocated by ReadFooter()
  int unused;
  int block_start; // start of my data
  MPI_Status status;
  int errcode;

  // verify that writing was initialized
  assert(restart_read_initialized);

  // root reads footer and distributes offset to each process
  if (rank == 0)
    ReadFooter(fd_in, ftr, unused);
  MPI_Scatter(ftr, 1, MPI_INT, &block_start, 1, MPI_INT, 0, comm);
  if (rank == 0)
    delete[] ftr;

  // read number of particles
  int count;
  int num_particles;
  MPI_Offset ofst = block_start;
  errcode = MPI_File_read_at_all(fd_in, ofst, &num_particles, sizeof(int),
				 MPI_BYTE, &status);
  if (errcode != MPI_SUCCESS)
    handle_mpio_error(errcode, (char *)"MPI_File_read_at_all num_particles");
  MPI_Get_count(&status, MPI_BYTE, &count);
  assert(count == sizeof(int));
  ofst += sizeof(int);

  // allocate data arrrays and create data type
  xx = new float[num_particles];
  yy = new float[num_particles];
  zz = new float[num_particles];
  vx = new float[num_particles];
  vy = new float[num_particles];
  vz = new float[num_particles];
  phi = new float[num_particles];
  pid = new int64_t[num_particles];
  mask = new uint16_t[num_particles];
  MPI_Datatype dtype;
  CreateRestartReadDtype(num_particles, xx, yy, zz, vx, vy, vz, phi, pid, 
			 mask, &dtype);

  // read the block
  errcode = MPI_File_read_at_all(fd_in, ofst, MPI_BOTTOM, 1, dtype, &status);
  if (errcode != MPI_SUCCESS)
    handle_mpio_error(errcode, (char *)"MPI_File_read_all datatype");

  // cleanup
  MPI_File_close(&fd_in);
  MPI_Type_free(&dtype);

  return num_particles;

}
//----------------------------------------------------------------------------
//
// independently writes the file footer
//
// fd: open MPI file handle
// blk_sizes: block sizes
// nb_out: total number of output blocks
//
void RestartIO::WriteFooter(MPI_File fd, const int *blk_sizes, int nb_out) {

  MPI_Status status;
  int errcode;
  int *footer = new int[nb_out + 1]; // footer ready for writing

  // populate the footer
  footer[0] = 0;
  for (int i = 1; i < (int)nb_out; i++)
    footer[i] = footer[i - 1] + blk_sizes[i - 1];
  footer[nb_out] = nb_out;

  // write the footer to disk
  MPI_File_seek(fd, 0, MPI_SEEK_END);
  errcode = MPI_File_write(fd, footer, nb_out + 1, MPI_INT, &status);
  int count;
  MPI_Get_count(&status, MPI_INT, &count);
  assert(count == (int)(nb_out + 1));

  // print the file size
  MPI_Offset ofst;
  MPI_Offset disp;
  MPI_File_get_position(fd, &ofst);
  MPI_File_get_byte_offset(fd, ofst, &disp);
  int size = disp / 1048576;
  fprintf(stderr, "Output file size is %d MB\n", size);

  delete[] footer;

}
//----------------------------------------------------------------------------
//
// independently reads the file footer
// footer in file is always ordered by global block id
// output footer is in the same order
// allocates ftr to be the correct size (caller should not allocate it)
//
// fd: open MPI file handle
// ftr: footer data (output)
// tb: total number of blocks (output)
//
void RestartIO::ReadFooter(MPI_File fd, int *&ftr, int &tb) {

  MPI_Status status;
  int errcode;
  int count;

  MPI_File_seek(fd, -sizeof(int), MPI_SEEK_END);

  // read the total number of blocks
  errcode = MPI_File_read(fd, &tb, 1, MPI_INT, &status);
  MPI_Get_count(&status, MPI_INT, &count);
  assert(count == 1);
  if (errcode != MPI_SUCCESS)
    handle_mpio_error(errcode, (char *)"MPI_File_read footer number of blocks");

  // read footer
  if (tb > 0) {
    ftr = new int[tb];
    MPI_File_seek(fd, -(tb + 1) * sizeof(int), MPI_SEEK_CUR);
    errcode = MPI_File_read(fd, ftr, tb, MPI_INT, &status);
    if (errcode != MPI_SUCCESS)
      handle_mpio_error(errcode, (char *)"MPI_File_read footer block start");
    MPI_Get_count(&status, MPI_INT, &count);
    assert(count == tb);
  }

}
//----------------------------------------------------------------------------
//
// creates an MPI struct datatype from a typemap
// use the datatype with the address MPI_BOTTOM
//
// addr: base address of the start of the typemap
// num_map_blocks: number of elements in the typemap
// typemap: typemap
// type: MPI datatype (output)
//
void RestartIO::CreateDtype(MPI_Aint addr, int num_map_blocks, map_block_t *typemap, 
			    MPI_Datatype *type) {

  // vector version of map, skipping map blocks with count 0
  vector<map_block_t>map(typemap, typemap + num_map_blocks);
  for (int i = 0; i < (int)map.size(); i++) {
    if (map[i].count <= 0) {
      map.erase(map.begin() + i);
      i--; // next time, test the new element that got moved up to this slot
    }
  }

  // form the required vectors
  vector <MPI_Aint> addrs(map.size(), addr);
  vector <int> counts(map.size(), 0);
  vector <MPI_Datatype> base_types(map.size(), 0);
  for(int i = 0; i < (int)map.size(); i++) {
    if (map[i].disp_type == OFST)
      addrs[i] = addrs[i] + map[i].disp;
    else
      addrs[i] = map[i].disp;
    counts[i] = map[i].count;
    base_types[i] = map[i].base_type;
  }

  // create the datatype
  MPI_Type_create_struct(map.size(), &counts[0], 
			 (MPI_Aint *)&addrs[0], &base_types[0], type);
  MPI_Type_commit(type);

}
//----------------------------------------------------------------------------
//
// creates datatype for writing restart data
// use with base address MPI_BOTTOM
//
// num_particles: number of particles in this block
// xx, yy, zz: particle positions
// vx, vy, vz: particle velocities
// phi: particle potentials
// pid: particle IDs
// mask: currently unused
// dtype: output datatype, allocated by caller
//
void RestartIO::CreateRestartWriteDtype(int &num_particles, float *xx, float *yy, 
					float *zz, float *vx, float *vy, float *vz, 
					float *phi, int64_t *pid, uint16_t *mask,
					MPI_Datatype *dtype) {
  MPI_Datatype base_type;
  CreateRestartDtype(num_particles, xx, yy, zz, vx, vy, vz, phi, pid, mask, 
		     &base_type);

  struct map_block_t map[] = {

    { MPI_INT,            ADDR, 1, mpi_addr(&num_particles) },
    { base_type,          ADDR, 1, mpi_addr(MPI_BOTTOM)     },

  };

  CreateDtype(NULL, 2, map, dtype);
  MPI_Type_free(&base_type);

}
//----------------------------------------------------------------------------
//
// creates datatype for reading restart data
// use with base address MPI_BOTTOM
//
// num_particles: number of particles in this block
// xx, yy, zz: particle positions
// vx, vy, vz: particle velocities
// phi: particle potentials
// pid: particle IDs
// mask: currently unused
// dtype: output datatype, allocated by caller
//
void RestartIO::CreateRestartReadDtype(int num_particles, float *xx, float *yy, 
				       float *zz, float *vx, float *vy, float *vz, 
				       float *phi, int64_t *pid, uint16_t *mask,
				       MPI_Datatype *dtype) {

  CreateRestartDtype(num_particles, xx, yy, zz, vx, vy, vz, phi, pid, mask, 
		     dtype);

}
//----------------------------------------------------------------------------
//
// creates datatype for reading and writing restart data
// use with base address MPI_BOTTOM
//
// num_particles: number of particles in this block
// xx, yy, zz: particle positions
// vx, vy, vz: particle velocities
// phi: particle potentials
// pid: particle IDs
// mask: currently unused
// dtype: output datatype, allocated by caller
//
void RestartIO::CreateRestartDtype(int num_particles, float *xx, float *yy, 
				   float *zz, float *vx, float *vy, float *vz, 
				   float *phi, int64_t *pid, uint16_t *mask,
				   MPI_Datatype *dtype) {
  
  struct map_block_t map[] = {

    { MPI_FLOAT,          ADDR, num_particles, mpi_addr(xx)             },
    { MPI_FLOAT,          ADDR, num_particles, mpi_addr(yy)             },
    { MPI_FLOAT,          ADDR, num_particles, mpi_addr(zz)             },
    { MPI_FLOAT,          ADDR, num_particles, mpi_addr(vx)             },
    { MPI_FLOAT,          ADDR, num_particles, mpi_addr(vy)             },
    { MPI_FLOAT,          ADDR, num_particles, mpi_addr(vz)             },
    { MPI_FLOAT,          ADDR, num_particles, mpi_addr(phi)            },
    { MPI_LONG_LONG_INT,  ADDR, num_particles, mpi_addr(pid)            },
    { MPI_UNSIGNED_SHORT, ADDR, num_particles, mpi_addr(mask)           },

  };

  CreateDtype(NULL, 9, map, dtype);

}
//----------------------------------------------------------------------------
//
// returns an MPI_Aint address given a pointer or address
//
// addr: pointer or address
//
// returns: MPI address
//
MPI_Aint RestartIO::mpi_addr(void *addr) {

  MPI_Aint p;
  MPI_Get_address(addr, &p);
  return p;

}
//----------------------------------------------------------------------------
//
// MPI error handler
//
void RestartIO::handle_mpio_error(int errcode, char *str) {

  char msg[MPI_MAX_ERROR_STRING];
  int resultlen;
  MPI_Error_string(errcode, msg, &resultlen);
  fprintf(stderr, "%s: %s\n", str, msg);
  MPI_Abort(comm, 1);

}
