/* Example driver for initializer */

#include <mpi.h>
#include <stdlib.h>
#include <iostream>
#include <string>
#include "TypesAndDefs.h"
#include "Basedata.h"
#include "Initializer.h"

#include "Definition.h"
#include "Partition.h"
#include "InitialExchange.h"
#include "distribution.h"

using namespace std;

string create_outName(string outBase, int rank) {
  ostringstream outName;
  outName << outBase << "." << rank;
  return outName.str();
}

void grid2phys(::real *pos_x, ::real *pos_y, ::real *pos_z,
	       ::real *vel_x, ::real *vel_y, ::real *vel_z,
	       long Npart, int np, float rL) {
  long i;
  
  float grid2phys_pos = 1.0*rL/np;
  float grid2phys_vel = 100.0*rL/np;
  
  for(i=0; i<Npart; i++) {
    pos_x[i] *= grid2phys_pos;
    pos_x[i] -= (pos_x[i] >= rL)*rL;

    pos_y[i] *= grid2phys_pos;
    pos_y[i] -= (pos_y[i] >= rL)*rL;

    pos_z[i] *= grid2phys_pos;
    pos_z[i] -= (pos_z[i] >= rL)*rL;

    vel_x[i] *= grid2phys_vel;
    vel_y[i] *= grid2phys_vel;
    vel_z[i] *= grid2phys_vel;
  }
  
  return;
}

void write_hcosmo(::real *pos_x, ::real *pos_y, ::real *pos_z,
		  ::real *vel_x, ::real *vel_y, ::real *vel_z,
		  ID_T *id, MASK_T *maskArr, MASK_T maskBit,
		  long Npart, const char *outName) {
  FILE *outFile;
  long i;
  ::real status=0.0;
  
  outFile = fopen(outName, "wb");
  for(i=0; i<Npart; i++) {
    if( maskArr[i] & maskBit ) {
      fwrite(&pos_x[i], sizeof(::real), 1, outFile);
      fwrite(&vel_x[i], sizeof(::real), 1, outFile);
      fwrite(&pos_y[i], sizeof(::real), 1, outFile);
      fwrite(&vel_y[i], sizeof(::real), 1, outFile);
      fwrite(&pos_z[i], sizeof(::real), 1, outFile);
      fwrite(&vel_z[i], sizeof(::real), 1, outFile);
      fwrite(&status, sizeof(::real), 1, outFile);
      fwrite(&id[i], sizeof(ID_T), 1, outFile);
    }
  }
  fclose(outFile);
  
  return;
}

vector<int>* readSubSampleFile(string subSampleName) {
  vector<int> *subSamples = new vector<int>;
  int tmp;
  ifstream inStream;
  inStream.open(subSampleName.c_str());
  while(1) {
    inStream >> tmp;
    if(!inStream.eof())
      subSamples->push_back(tmp);
    else
      break;
  }
  return subSamples;
}

int main(int argc, char *argv[]) {
  
  if(argc < 6) {
    fprintf(stderr,"USAGE: init2 <indatName> <tfName> <outBase> <ROUND_ROBIN|ONE_TO_ONE> <K|R> [subSampleName]\n");
    fprintf(stderr,"<K|R> = start with density in K space or with white noise in R space\n");
    exit(-1);
  }
  string indatName = argv[1];
  string tfName = argv[2];
  string outBase = argv[3];
  string distributeType = argv[4];
  string KorR = argv[5];

  int useWN = 0;
  if(KorR == "R")
    useWN = 1;
  else if (KorR != "K") {
    fprintf(stderr,"ERROR: KorR = %s\n",KorR);
    exit(-1);
  }

  int oneToOneQ = 0;
  if(distributeType == "ONE_TO_ONE")
    oneToOneQ = 1;

  int subSampleQ = 0;
  string subSampleName;
  if(argc > 6) {
    subSampleName = argv[6];
    subSampleQ = 1;
  }

  Basedata bdata(indatName.c_str());

  ::real *pos_x, *pos_y, *pos_z, *vel_x, *vel_y, *vel_z;
  int64_t Npart, np, NumProcs64;
  int NumProcs, MyPE;
  
  MPI::Init();

  MPI_Comm_size(MPI_COMM_WORLD, &NumProcs);
  MPI_Comm_rank(MPI_COMM_WORLD, &MyPE);

  np = bdata.np();
  NumProcs64 = NumProcs;

  Npart = (np*np*np)/NumProcs;

  pos_x = (::real *)malloc(Npart*sizeof(::real));
  pos_y = (::real *)malloc(Npart*sizeof(::real));
  pos_z = (::real *)malloc(Npart*sizeof(::real));
  vel_x = (::real *)malloc(Npart*sizeof(::real));
  vel_y = (::real *)malloc(Npart*sizeof(::real));
  vel_z = (::real *)malloc(Npart*sizeof(::real));

  init_particles(pos_x, pos_y, pos_z, vel_x, vel_y, vel_z, 
		 bdata, tfName.c_str(), useWN);

  grid2phys(pos_x, pos_y, pos_z, vel_x, vel_y, vel_z,
            Npart, bdata.np(), bdata.rL() );

  ID_T *id, irank;
  MASK_T *maskArr;
  MASK_T one = 1;
  MASK_T nbits = 8*sizeof(MASK_T);
  MASK_T maskBit = one << (nbits-1);

  id = (ID_T *)malloc(Npart*sizeof(ID_T));
  maskArr = (MASK_T *)malloc(Npart*sizeof(MASK_T));
  for(ID_T i=0; i<Npart; i++) {
     id[i] = Npart*MyPE + i;
     maskArr[i] = maskBit;
  }

  //mark subsamples with mask bits
  int nss=0;
  vector<int> *subSamples;
  if(subSampleQ) {
    distribution_t dist;
    int n_dist[3];
    int padding[3] = {0,0,0};
    int nproc[3], self[3];
    int n_local[3], lo_global[3], hi_global[3];
    for(int i=0; i<3; i++)
      n_dist[i] = bdata.np();
    distribution_init(MPI_COMM_WORLD, n_dist, padding, &dist, false);
    distribution_assert_commensurate(&dist);
    for(int i=0; i<3; i++) {
      nproc[i] = distribution_get_nproc_3d(&dist, i);
      self[i] = distribution_get_self_3d(&dist, i);
      n_local[i] = dist.process_topology_3.n[i];
      
      lo_global[i] = n_local[i]*self[i];
      hi_global[i] = lo_global[i] + n_local[i];
    }
    /*
    printf("%d/%d: (%d, %d, %d) (%d, %d, %d)\n",rank, NumProcs,
           lo_global[0], lo_global[1], lo_global[2],
           hi_global[0], hi_global[1], hi_global[2]);
    fflush(stdout);
    */
    distribution_fini(&dist);

    subSamples = readSubSampleFile(subSampleName);
    if(subSamples->size() > nbits-1) {
      fprintf(stderr,"WARNING: only doing first %d subsamples.\n",nbits-1);
      fflush(stderr);
      subSamples->resize(nbits-1);
    }
    nss = subSamples->size();

    for(int ssn=0; ssn<nss; ssn++) {
      int ss = (*subSamples)[ssn];
      maskBit = one << ssn;
      long cntr = 0;
      for(int i=0; i<n_local[0]; i++)
	for(int j=0; j<n_local[1]; j++)
	  for(int k=0; k<n_local[2]; k++) {
	    MASK_T tmp = maskBit;
	    tmp *= !( (i+lo_global[0])%ss );
	    tmp *= !( (j+lo_global[1])%ss );
	    tmp *= !( (k+lo_global[2])%ss );
	    maskArr[cntr] = maskArr[cntr] | tmp;
	    cntr++;
	  }
    }
  }

  //do initialexchange if one to one
  vector< ::real > *vpos_x, *vpos_y, *vpos_z, *vvel_x, *vvel_y, *vvel_z, *vpot;
  vector<MASK_T> *vmaskArr;
  vector<ID_T> *vid;
  vector<STATUS_T> *vstatus;

  if(oneToOneQ && NumProcs > 1) {
    //alloc vectors
    vpos_x = new vector< ::real >;
    vpos_y = new vector< ::real >;
    vpos_z = new vector< ::real >;
    vvel_x = new vector< ::real >;
    vvel_y = new vector< ::real >;
    vvel_z = new vector< ::real >;
    vpot = new vector< ::real >;
    vmaskArr = new vector<MASK_T>;
    vid = new vector<ID_T>;
    vstatus = new vector<STATUS_T>;

    //reserve vectors
    float Npartf = static_cast<float>(Npart);
    float Nresf = Npartf + 3.0*pow((double) Npartf,2.0/3.0);
    long Nres = static_cast<long>(Nresf);
    vpos_x->reserve(Nres);
    vpos_y->reserve(Nres);
    vpos_z->reserve(Nres);
    vvel_x->reserve(Nres);
    vvel_y->reserve(Nres);
    vvel_z->reserve(Nres);
    vpot->reserve(Nres);
    vmaskArr->reserve(Nres);
    vid->reserve(Nres);

    //fake potential array
    ::real *pot = (::real *)malloc(Npart*sizeof(::real));
    for(long i=0; i<Npart; i++)
      pot[i] = 0.0;

    //set up exchanger
    MPI_Barrier(MPI_COMM_WORLD);
    Partition::initialize();
    InitialExchange exchange;

    //
    exchange.setParameters(bdata.rL(), INITIAL_EXCHANGE_FUDGE*bdata.rL()/bdata.np());

    exchange.initialize();
    exchange.setParticleArrays(Npart,
			       pos_x, pos_y, pos_z,
			       vel_x, vel_y, vel_z, 
			       pot, id, maskArr);
    exchange.setParticleVectors(vpos_x, vpos_y, vpos_z,
				vvel_x, vvel_y, vvel_z,
				vpot, vid, vmaskArr, vstatus);
    exchange.exchangeParticles();

    MPI_Barrier(MPI_COMM_WORLD);

    //free up arrays
    free(pos_x);
    free(pos_y);
    free(pos_z);
    free(vel_x);
    free(vel_y);
    free(vel_z);
    free(pot);
    free(id);
    free(maskArr);

    //re-point arrays
    pos_x = &(*vpos_x)[0];
    pos_y = &(*vpos_y)[0];
    pos_z = &(*vpos_z)[0];
    vel_x = &(*vvel_x)[0];
    vel_y = &(*vvel_y)[0];
    vel_z = &(*vvel_z)[0];
    id = &(*vid)[0];
    maskArr = &(*vmaskArr)[0];

    //reset Npart
    Npart = vpos_x->size();

    //get rid of fake potential
    delete vpot;

    delete vstatus;
  }

  MPI_Barrier(MPI_COMM_WORLD);

  //output
  outBase += ".hcosmo";
  if(!subSampleQ) {
    maskBit = one << (nbits-1);
    write_hcosmo(pos_x, pos_y, pos_z, vel_x, vel_y, vel_z,
		 id, maskArr, maskBit, Npart,
		 create_outName(outBase,MyPE).c_str() );
  } else {
    int np = bdata.np();
    maskBit = one << (nbits-1);
    write_hcosmo(pos_x, pos_y, pos_z, vel_x, vel_y, vel_z,
		 id, maskArr, maskBit, Npart,
		 create_outName(create_outName(outBase,np),MyPE).c_str() );

    for(int ssn=0; ssn<nss; ssn++) {
      int ss = (*subSamples)[ssn];
      np = bdata.np()/ss;
      maskBit = one << ssn;
      write_hcosmo(pos_x, pos_y, pos_z, vel_x, vel_y, vel_z,
		   id, maskArr, maskBit, Npart,
		   create_outName(create_outName(outBase,np),MyPE).c_str());
    }
  }

  MPI_Barrier(MPI_COMM_WORLD);

  //cleanup memory
  if(oneToOneQ) {
    delete vpos_x;
    delete vpos_y;
    delete vpos_z;
    delete vvel_x;
    delete vvel_y;
    delete vvel_z;
    delete vid;
    delete vmaskArr;
    Partition::finalize();
  } else {
    free(pos_x);
    free(pos_y);
    free(pos_z);
    free(vel_x);
    free(vel_y);
    free(vel_z);
    free(id);
    free(maskArr);
  }

  MPI_Barrier(MPI_COMM_WORLD);
  MPI_Finalize();

  return(0);
}
