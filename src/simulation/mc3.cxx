#include "mc3.h"
#include "Initializer.h"
#include "InitialExchange.h"

#include "bigchunk.h"

#include <algorithm>

//MISC



// Functor for assert(x == 1)
struct AssertEqualsOne {
  void operator() (int x) { assert(x == 1 && "error in recvbuf"); }
};

// Conditionally execute a trivial all to all.
void initial_all_to_all(bool flag)
{
  bool const debug = true;

  if (flag) {
    int rank;
    int size;
    std::vector<int> sendbuf;
    std::vector<int> recvbuf;

    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    if (rank == 0) {
      std::cout << "Executing initial all-to-all...";
    }

    sendbuf.resize(size);
    recvbuf.resize(size);
    if (debug) {
      sendbuf.assign(size, 1);
      recvbuf.assign(size, 0);
    }

    MPI_Alltoall(&sendbuf[0], 1, MPI_INT,
		 &recvbuf[0], 1, MPI_INT,
		 MPI_COMM_WORLD);

    if (debug) {
      std::for_each(recvbuf.begin(), recvbuf.end(), AssertEqualsOne());
    }
    if (rank == 0) {
      std::cout << "done." << std::endl;
    }
  }
}



//MISC (COPY TO OPTIONS)



string create_outName(string outBase, int rank) {
  ostringstream outName;
  outName << outBase << "." << rank;
  return outName.str();
}



//MISC (MOVE TO OPTIONS)



/*
vector<int>* readInts(string inName) {
  vector<int> *ret = new vector<int>;
  int tmp;
  ifstream inStream;
  inStream.open(inName.c_str());
  while(1) {
    inStream >> tmp;
    if(!inStream.eof())
      ret->push_back(tmp);
    else
      break;
  }
  return ret;
}



int intInVec(int n, vector<int> *iv) {
  int ret=0;

  if(iv)
    for(int i=0; i<iv->size(); i++)
      if( (*iv)[i] == n )
	ret = 1;

  return ret;
}
*/



//GRID



void output_array_alive(GRID_T *arr, string outName) {
  FILE *outFile;
  int i, j, k, ca[DIMENSION], ngla[DIMENSION], nglt[DIMENSION], ngo;
  int64_t gindx, ng, ig, jg, kg;
  GRID_T val;

  Domain::corner_grid_alive(ca);
  ng = Domain::ng();
  Domain::ng_local_alive(ngla);
  Domain::ng_local_total(nglt);
  ngo = Domain::ng_overload();

  outFile = fopen(outName.c_str(),"w");

  for(i=0; i<ngla[0]; i++) {
    for(j=0; j<ngla[1]; j++) {
      for(k=0; k<ngla[2]; k++) {

        val = arr[_l_array_index(i+ngo, j+ngo, k+ngo, nglt[1], nglt[2])];
	ig = i + ca[0];
	jg = j + ca[1];
	kg = k + ca[2];
        gindx = _g_array_index(ig, jg, kg, ng, ng);
        fprintf(outFile,"%d %d %d %d %f\n",
                gindx,ig,jg,kg,val);

      }
    }
  }

  fclose(outFile);

  return;
}



void output_array_total(GRID_T *arr, string outName) {
  FILE *outFile;
  int i, j, k, ct[DIMENSION], nglt[DIMENSION];
  int64_t gindx, ng, ig, jg, kg;
  GRID_T val;

  Domain::corner_grid_total(ct);
  ng = Domain::ng();
  Domain::ng_local_total(nglt);

  outFile = fopen(outName.c_str(),"w");
  for(i=0; i<nglt[0]; i++) {
    for(j=0; j<nglt[1]; j++) {
      for(k=0; k<nglt[2]; k++) {
	val = arr[_l_array_index(i, j, k, nglt[1], nglt[2])];
	ig = i + ct[0];
	jg = j + ct[1];
	kg = k + ct[2];
	if(ig < 0) ig += ng;
	if(jg < 0) jg += ng;
	if(kg < 0) kg += ng;
	if(ig >= ng) ig -= ng;
	if(jg >= ng) jg -= ng;
	if(kg >= ng) kg -= ng;
	gindx = _g_array_index(ig,jg,kg,ng,ng);
	fprintf(outFile,"%d %d %d %d %f\n",
		gindx,ig,jg,kg,val);
      }
    }
  }
  fclose(outFile);
  
  return;
}



double sum_rho_alive(GRID_T *arr) {
  int i, j, k, ca[DIMENSION], ngla[DIMENSION], nglt[DIMENSION], ngo;
  double val=0.0;

  Domain::corner_grid_alive(ca);
  Domain::ng_local_alive(ngla);
  Domain::ng_local_total(nglt);
  ngo = Domain::ng_overload();

  for(i=0; i<ngla[0]; i++) {
    for(j=0; j<ngla[1]; j++) {
      for(k=0; k<ngla[2]; k++) {
        val += arr[_l_array_index(i+ngo, j+ngo, k+ngo, nglt[1], nglt[2])];
      }
    }
  }

  return val;
}



void fft_copy_arr_forward(GRID_T *arr, COMPLEX_T *fft_arr) {
  int i, j, k;
  int ngla[DIMENSION], nglt[DIMENSION], ngo;

  Domain::ng_local_alive(ngla);
  Domain::ng_local_total(nglt);
  ngo = Domain::ng_overload();

  for(i=0; i<ngla[0]; i++) {
    for(j=0; j<ngla[1]; j++) {
      for(k=0; k<ngla[2]; k++) {
	fft_arr[_l_array_index(i,j,k,ngla[1],ngla[2])] =
	  arr[_l_array_index(i+ngo, j+ngo, k+ngo, nglt[1], nglt[2])];
      }
    }
  }
 
  return;
}



void fft_copy_arr_backward(COMPLEX_T *fft_arr, GRID_T *arr) {
  int i, j, k;
  int ngla[DIMENSION], nglt[DIMENSION], ngo;

  Domain::ng_local_alive(ngla);
  Domain::ng_local_total(nglt);
  ngo = Domain::ng_overload();

  for(i=0; i<ngla[0]; i++) {
    for(j=0; j<ngla[1]; j++) {
      for(k=0; k<ngla[2]; k++) {
	arr[_l_array_index(i+ngo, j+ngo, k+ngo, nglt[1], nglt[2])] =
	  std::real(fft_arr[_l_array_index(i, j, k, ngla[1], ngla[2])]);
      }
    }
  }

  return;
}



void poisson_alloc(COMPLEX_T **fft_rho_arr,
		   COMPLEX_T **fft_grad_phi_arr,
		   COMPLEX_T **fft_phi_arr) {
  int Ngla = Domain::Ng_local_alive();
  int ng = Domain::ng();

  //maybe this should be fftw_malloc'd
  *fft_rho_arr = (COMPLEX_T*) bigchunk_malloc(Ngla*sizeof(COMPLEX_T));
  *fft_grad_phi_arr = *fft_rho_arr;
  *fft_phi_arr = *fft_rho_arr;

  return;
}



void poisson_free(COMPLEX_T **fft_rho_arr,
		  COMPLEX_T **fft_grad_phi_arr,
		  COMPLEX_T **fft_phi_arr) {
  bigchunk_free(*fft_rho_arr);
  *fft_rho_arr = *fft_grad_phi_arr = *fft_phi_arr = NULL;

  return;
}



void writePk(SolverBase *solver, string outName) {
  vector<double> power;

  solver->power_spectrum(power);
  //MPI_Barrier(MPI_COMM_WORLD);

  int rank = Partition::getMyProc();
  if(rank==0) {
    float rL = Domain::rL();
    float ng = 1.0*Domain::ng();
    float pi = 4.0*atan(1.0);
    float kcoeff = 2.0*pi/rL;
    float pkcoeff = powf(rL/ng,3.0);

    FILE *outFile = fopen(outName.c_str(),"w");

    for(int i=0; i<power.size(); i++)
      fprintf(outFile, "%e\t%e\n", kcoeff*i, pkcoeff*power[i]);

    fclose(outFile);
  }

  return;
}



void map2_poisson_forward(Particles & particles,
			  SolverBase *solver,
			  GRID_T *rho_arr,
			  COMPLEX_T *fft_rho_arr) {

  // SimpleTimings::TimerRef t_cic = SimpleTimings::getTimer("cic");
  // SimpleTimings::TimerRef t_copyf = SimpleTimings::getTimer("copyf");
  // SimpleTimings::TimerRef t_fftf = SimpleTimings::getTimer("fftf");

  //CIC
  // SimpleTimings::startTimer(t_cic);
  particles.cic();
  // SimpleTimings::stopTimerStats(t_cic);
  //MPI_Barrier(MPI_COMM_WORLD);
  
  //COPY RHO FOR POISSON SOLVE
  // SimpleTimings::startTimer(t_copyf);
  fft_copy_arr_forward(rho_arr, fft_rho_arr);
  // SimpleTimings::stopTimerStats(t_copyf);
  //MPI_Barrier(MPI_COMM_WORLD);
  
  //FORWARD POISSON SOLVE
  // SimpleTimings::startTimer(t_fftf);
  solver->forward_solve(fft_rho_arr);
  // SimpleTimings::stopTimerStats(t_fftf);
  //MPI_Barrier(MPI_COMM_WORLD);

  return;
}



void map2_poisson_backward_gradient(Particles & particles,
				    SolverBase *solver,
				    GRID_T *grad_phi_arr,
				    COMPLEX_T *fft_grad_phi_arr,
				    GridExchange & gexchange,
				    TimeStepper & ts,
				    TS_FLOAT stepFraction) {

  // SimpleTimings::TimerRef t_fftb = SimpleTimings::getTimer("fftb");
  // SimpleTimings::TimerRef t_copyb = SimpleTimings::getTimer("copyb");
  // SimpleTimings::TimerRef t_xch = SimpleTimings::getTimer("xch");
  // SimpleTimings::TimerRef t_cici = SimpleTimings::getTimer("cici");

  //FOR EACH COMPONENT OF THE GRADIENT OF PHI
  for(int d=0; d<DIMENSION; d++) {
    
    //BACKWARD POISSON SOLVE
    // SimpleTimings::startTimer(t_fftb);
    solver->backward_solve_gradient(d, fft_grad_phi_arr);
    // SimpleTimings::stopTimerStats(t_fftb);
    //MPI_Barrier(MPI_COMM_WORLD);
    
    //COPY GRAD_PHI FROM POISSON SOLVE
    // SimpleTimings::startTimer(t_copyb);
    fft_copy_arr_backward(fft_grad_phi_arr, grad_phi_arr);
    // SimpleTimings::stopTimerStats(t_copyb);
    //MPI_Barrier(MPI_COMM_WORLD);
    
    //EXCHANGE GRID FOR OVERLOADING CIC^-1
    // SimpleTimings::startTimer(t_xch);
    gexchange.exchangeGrid(grad_phi_arr);
    // SimpleTimings::stopTimerStats(t_xch);
    //MPI_Barrier(MPI_COMM_WORLD);
    
    //CIC^-1
    // SimpleTimings::startTimer(t_cici);
    particles.inverse_cic(ts.tau()*stepFraction, ts.fscal(), d);
    // SimpleTimings::stopTimerStats(t_cici);
    //MPI_Barrier(MPI_COMM_WORLD);
  }
  //MPI_Barrier(MPI_COMM_WORLD);

  return;
}



void map2_poisson_backward_potential(Particles & particles,
				     SolverBase *solver,
				     GRID_T *phi_arr,
				     COMPLEX_T *fft_phi_arr,
				     GridExchange & gexchange) {

  // SimpleTimings::TimerRef t_fftb = SimpleTimings::getTimer("fftb");
  // SimpleTimings::TimerRef t_copyb = SimpleTimings::getTimer("copyb");
  // SimpleTimings::TimerRef t_xch = SimpleTimings::getTimer("xch");
  // SimpleTimings::TimerRef t_cici = SimpleTimings::getTimer("cici");

  // SimpleTimings::startTimer(t_fftb);
  solver->backward_solve(fft_phi_arr);
  // SimpleTimings::stopTimerStats(t_fftb);
  //MPI_Barrier(MPI_COMM_WORLD);
  
  // SimpleTimings::startTimer(t_copyb);
  fft_copy_arr_backward(fft_phi_arr, phi_arr);
  // SimpleTimings::stopTimerStats(t_copyb);
  //MPI_Barrier(MPI_COMM_WORLD);
  
  // SimpleTimings::startTimer(t_xch);
  gexchange.exchangeGrid(phi_arr);
  // SimpleTimings::stopTimerStats(t_xch);
  //MPI_Barrier(MPI_COMM_WORLD);
  
  // SimpleTimings::startTimer(t_cici);
  particles.inverse_cic_potential();
  // SimpleTimings::stopTimerStats(t_cici);
  //MPI_Barrier(MPI_COMM_WORLD);

  return;
}



//PARTICLES

void grid2phys(POSVEL_T *pos_x, POSVEL_T *pos_y, POSVEL_T *pos_z,
	       POSVEL_T *vel_x, POSVEL_T *vel_y, POSVEL_T *vel_z,
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

void loadParticles(Basedata & indat,
		   Particles & particles,
		   string inBase,
		   string dataType,
		   string distributeType,
		   int cosmoFormatQ,
		   MC3Options & options)
{
  // SimpleTimings::TimerRef t_sort = SimpleTimings::getTimer("sort");
  int numranks = Partition::getNumProc();
  int rank = Partition::getMyProc();

  //SET UP LOCAL PARTICLE VECTORS
  vector<POSVEL_T> *xx, *yy, *zz, *vx, *vy, *vz, *phi, *mass;
  vector<ID_T> *tag;
  vector<MASK_T> *mask;
  vector<STATUS_T> *status;
  
  allocReserveVectors(&xx, &yy, &zz, &vx, &vy, &vz, &mass,
		      &phi, &tag, &mask, &status,
                      Domain::maxNpLocal(), ReserveXV);
  //MPI_Barrier(MPI_COMM_WORLD);
  
  
  //READ AND DISTRIBUTE PARTICLES
  ParticleDistribute *distribute = new ParticleDistribute;
  int tmpNpla;
  if (dataType == "INIT") {

    uint64_t Npart, np, NumProcs64;
    np = indat.np();
    NumProcs64 = numranks;
    Npart = (np*np*np)/NumProcs64;

    string tfName = inBase;

    xx->resize(Npart);
    yy->resize(Npart);
    zz->resize(Npart);

    vx->resize(Npart);
    vy->resize(Npart);
    vz->resize(Npart);

    init_particles(&xx->at(0), &yy->at(0), &zz->at(0),
                   &vx->at(0), &vy->at(0), &vz->at(0),
                   indat, tfName.c_str(), options.whiteNoiseInit());

    grid2phys(&xx->at(0), &yy->at(0), &zz->at(0),
              &vx->at(0), &vy->at(0), &vz->at(0),
              Npart, indat.np(), indat.rL());

    reserveVectors(&xx, &yy, &zz, &vx, &vy, &vz, &mass,
  		        &phi, &tag, &mask, &status,
                        Domain::maxNpLocal(), ReserveNonXV);

    MASK_T one = 1;
    MASK_T nbits = 8*sizeof(MASK_T);
    MASK_T maskBit = one << (nbits-1);
    mask->resize(Npart);
    std::fill(mask->begin(), mask->end(), maskBit);

    tag->resize(Npart);
    for (ID_T i = 0; i < Npart; ++i) {
        tag->at(i) = Npart*rank + i;
	mask->at(i) = maskBit;
    }

    // fake potential array
    phi->resize(Npart);
    std::fill(phi->begin(), phi->end(), 0.0);
 
    // now we need to redistribute the particles...
    // make copies to be the source for initial exchange
    vector<POSVEL_T> xx_old(*xx), yy_old(*yy), zz_old(*zz),
                     vx_old(*vx), vy_old(*vy), vz_old(*vz);
    vector<POSVEL_T> phi_old(*phi);
    vector<ID_T>     tag_old(*tag);
    vector<MASK_T>   mask_old(*mask);

    vector<STATUS_T> status;

    // clear out destination vectors
    xx->clear();
    yy->clear();
    zz->clear();
    vx->clear();
    vy->clear();
    vz->clear();
    phi->clear();
    tag->clear();
    mask->clear();

    // increase the memory reservation for the destination
    long Nres = (long) ((double) Npart + 3.0*pow((double) Npart,2.0/3.0));
    xx->reserve(Nres);
    yy->reserve(Nres);
    zz->reserve(Nres);
    vx->reserve(Nres);
    vy->reserve(Nres);
    vz->reserve(Nres);
    phi->reserve(Nres);
    tag->reserve(Nres);
    mask->reserve(Nres);
   
    //MPI_Barrier(MPI_COMM_WORLD);

    Partition::initialize();
    InitialExchange exchange;

    //
    exchange.setParameters(indat.rL(), INITIAL_EXCHANGE_FUDGE*indat.rL()/indat.np());

    exchange.initialize();
    exchange.setParticleArrays(Npart, &xx_old[0], &yy_old[0], &zz_old[0],
                               &vx_old[0], &vy_old[0], &vz_old[0], &phi_old[0],
                               &tag_old[0], &mask_old[0]);
    exchange.setParticleVectors(xx, yy, zz, vx, vy, vz,
                                phi, tag, mask, &status);
    exchange.exchangeParticles();

    //MPI_Barrier(MPI_COMM_WORLD);
  }
  else {
    reserveVectors(&xx, &yy, &zz, &vx, &vy, &vz, &mass,
  		        &phi, &tag, &mask, &status,
                        Domain::maxNpLocal(), ReserveNonXV);

    float rLdistribute = indat.rL();
    if(cosmoFormatQ)
      rLdistribute *= 1.0/indat.hubble();
    distribute->setParameters(inBase, rLdistribute, dataType);
    distribute->initialize();
    distribute->setParticles(xx, yy, zz, vx, vy, vz, mass, tag);
  
    if (distributeType == "ROUND_ROBIN")
      distribute->readParticlesRoundRobin(0);
    else if (distributeType == "ALL_TO_ALL")
      distribute->readParticlesAllToAll(0);
    else if (distributeType == "ONE_TO_ONE")
      distribute->readParticlesOneToOne(0);
  
    tmpNpla = xx->size();
  }

  mass->clear();
  phi->clear();
  mask->clear();
 
  //FAKE POTENTIAL VECTOR
  for(int i=0; i<xx->size(); i++) {
    mass->push_back(1.0);
    phi->push_back(0.0);
    mask->push_back((MASK_T)0);
  }
  
  //MPI_Barrier(MPI_COMM_WORLD);
  
  
  //CONVERT Mpc to Mpc/h
  if(cosmoFormatQ)
    coords_Mpc2Mpch(*xx, *yy, *zz, indat.hubble());
  //MPI_Barrier(MPI_COMM_WORLD);
  
  //EXCHANGE PARTICLES FOR OVERLOADING
  float rLexchange = indat.rL();
  float oLexchange = Domain::oL();
  ParticleExchange *pexchange = new ParticleExchange;
  pexchange->setParameters(rLexchange, oLexchange);
  pexchange->initialize();
  pexchange->setParticles(xx, yy, zz, vx, vy, vz, mass, phi, 
			  tag, mask, status);
  pexchange->exchangeParticles();
  //MPI_Barrier(MPI_COMM_WORLD);
  
  
#if 0
  printf("INITIAL PXCH rank = %d Npla = %d nReserve = %d size = %d\n",
	 Partition::getMyProc(), tmpNpla, Domain::maxNpLocal(), xx->size());
  fflush(stdout);
  //MPI_Barrier(MPI_COMM_WORLD);
#endif
  

  //COPY FROM LOCAL VECTORS
  particles.copyAllFromVectors(&xx, &yy, &zz, &vx, &vy, &vz, &phi, &tag, 
			       &mask, 1.0/(1.0+indat.zin()), 1);

  
  //DROP LOCAL VECTORS
  deleteVectors(&xx,&yy,&zz,&vx,&vy,&vz,&mass,&phi,&tag,&mask,&status);
  delete distribute;
  delete pexchange;
  //MPI_Barrier(MPI_COMM_WORLD);
  
  // If not already done, allocate the particle field.
  particles.allocField();
  
  //SORT PARTICLES
  // SimpleTimings::startTimer(t_sort);
  particles.resortParticles();
  // SimpleTimings::stopTimerStats(t_sort);
  if(rank==0) printf("\n");
  //MPI_Barrier(MPI_COMM_WORLD);

  return;
}



void coords_Mpc2Mpch(vector<POSVEL_T> xVec,
		     vector<POSVEL_T> yVec,
		     vector<POSVEL_T> zVec,
		     float hubble) {
  int i;
  int Np = xVec.size();

  for(i=0; i<Np; i++) {
    xVec[i] *= hubble;
    yVec[i] *= hubble;
    zVec[i] *= hubble;
  }

  return;
}



void coords_Mpch2Mpc(vector<POSVEL_T> xVec,
		     vector<POSVEL_T> yVec,
		     vector<POSVEL_T> zVec,
		     float hubble) {
  int i;
  float hubblei = 1.0/hubble;
  int Np = xVec.size();

  for(i=0; i<Np; i++) {
    xVec[i] *= hubblei;
    yVec[i] *= hubblei;
    zVec[i] *= hubblei;
  }

  return;
}



//PARTICLES (COPY TO OPTIONS)



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
			 int Np, unsigned mode) {
  *xx = new vector<POSVEL_T>;
  *yy = new vector<POSVEL_T>;
  *zz = new vector<POSVEL_T>;
  *vx = new vector<POSVEL_T>;
  *vy = new vector<POSVEL_T>;
  *vz = new vector<POSVEL_T>;
  *mass = new vector<POSVEL_T>;
  *phi = new vector<POSVEL_T>;
  *tag = new vector<ID_T>;
  *mask = new vector<MASK_T>;
  *status = new vector<STATUS_T>;

  reserveVectors(xx, yy, zz, vx, vy, vz,
                 mass, phi, tag, mask, status,
                 Np, mode);

  return;
}



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
		    int Np, unsigned mode) {
  bool reserveXV = mode & ReserveXV;
  bool reserveNonXV = mode & ReserveNonXV;

  if (reserveXV) {
    (*xx)->reserve(Np);
    (*yy)->reserve(Np);
    (*zz)->reserve(Np);
    (*vx)->reserve(Np);
    (*vy)->reserve(Np);
    (*vz)->reserve(Np);
  }

  if (reserveNonXV) {
    (*mass)->reserve(Np);
    (*phi)->reserve(Np);
    (*tag)->reserve(Np);
    (*mask)->reserve(Np);
    (*status)->reserve(Np);
  }

  return;
}



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
		   vector<STATUS_T> **status) {
  if(*xx) delete *xx;
  *xx = NULL;

  if(*yy) delete *yy;
  *yy = NULL;

  if(*zz) delete *zz;
  *zz = NULL;

  if(*vx) delete *vx;
  *vx = NULL;

  if(*vy) delete *vy;
  *vy = NULL;

  if(*vz) delete *vz;
  *vz = NULL;

  if(*mass) delete *mass;
  *mass = NULL;

  if(*phi) delete *phi;
  *phi = NULL;

  if(*id) delete *id;
  *id = NULL;

  if(*mask) delete *mask;
  *mask = NULL;

  if(*status) delete *status;
  *status = NULL;

  return;
}



//PARTICLES (MOVE TO OPTIONS)



/*
void findHalos(Particles & particles,
	       Basedata & indat,
	       Halodata & halodat,
	       string outBase,
	       float anow,
	       int step,
	       LightCone *lc,
	       int applyLC)
{
  //INITIALIZE EXCHANGE
  float rL = indat.rL();
  float oL = halodat.oL();
  ParticleExchange *hexchange = new ParticleExchange;
  hexchange->setParameters(rL, oL);
  hexchange->initialize();


  //FIGURE OUT HOW MUCH TO RESERVE
  float rla[DIMENSION];
  Domain::rL_local_alive(rla);
  float rlt[DIMENSION];
  for(int i=0; i<DIMENSION; i++)
    rlt[i] = rla[i] + 2.0*oL;
  int Npla = particles.Np_local_alive();
  int nReserve = static_cast<int>( ceil( FUDGE_HALO_OL*(rlt[0]/rla[0])*(rlt[1]/rla[1])*(rlt[2]/rla[2])*Npla ) );


  //GRAB ALIVE PARTICLES
  vector<POSVEL_T> *xx, *yy, *zz, *vx, *vy, *vz, *phi;
  vector<POSVEL_T> *mass;
  vector<ID_T> *tag;
  vector<MASK_T> *mask;
  vector<STATUS_T> *status;
  allocReserveVectors(&xx, &yy, &zz, &vx, &vy, &vz, &mass, &phi, 
		      &tag, &mask, &status, nReserve);
  particles.copyAliveIntoVectors(xx,yy,zz,vx,vy,vz,phi,tag,mask,anow);
  particles.shoveParticles();

  //FIGURE OUT MASS OF PARTICLE
  float npf = static_cast<float>(indat.np());
  float omegadm = indat.omegadm();
  float hubble = indat.hubble();
  float deut = indat.deut();
  float omegatot = omegadm + deut/hubble/hubble;
  float mp = 2.77536627e11*rL*rL*rL*omegatot/npf/npf/npf;

  for(int i=0; i<xx->size(); i++)
    mass->push_back(mp);

  hexchange->setParticles(xx, yy, zz, vx, vy, vz, mass, phi, tag, mask, status);
  //MPI_Barrier(MPI_COMM_WORLD);


  //EXCHANGE
  hexchange->exchangeParticles();
  //MPI_Barrier(MPI_COMM_WORLD);


  printf("HALO    PXCH rank = %d Npla = %d nReserve = %d size = %d\n",
	 Partition::getMyProc(), Npla, nReserve, xx->size());
  fflush(stdout);
  //MPI_Barrier(MPI_COMM_WORLD);


  //FIGURE OTHER FOF PARAMETERS
  float bb = halodat.bb();
  float minmass = halodat.minmass();
  int nmin = static_cast<int>(floor(minmass/mp));
  if(nmin < halodat.nmin())
    nmin = halodat.nmin();


  //FIND HALOS
  CosmoHaloFinderP *haloFinder = new CosmoHaloFinderP;
  int np = indat.np();

  string haloParticleSuffix;
  if(applyLC)
    haloParticleSuffix = LC_HALO_PARTICLE_SUFFIX;
  else
    haloParticleSuffix = STATIC_HALO_PARTICLE_SUFFIX;
  string hfOutBase = create_outName(outBase+"."+haloParticleSuffix,step);

  haloFinder->setParameters(hfOutBase, rL, oL, np, nmin, bb);
  haloFinder->setParticles(xx, yy, zz, vx, vy, vz, phi, tag, mask, status);
  //MPI_Barrier(MPI_COMM_WORLD);
  haloFinder->executeHaloFinder();
  //MPI_Barrier(MPI_COMM_WORLD);
  haloFinder->collectHalos();
  //MPI_Barrier(MPI_COMM_WORLD);
  haloFinder->mergeHalos();
  //MPI_Barrier(MPI_COMM_WORLD);


  //FIND HALO PROPERTIES
  int numberOfHalos = haloFinder->getNumberOfHalos();
  int *halos = haloFinder->getHalos();
  int *haloCount = haloFinder->getHaloCount();
  int *haloList = haloFinder->getHaloList();

  FOFHaloProperties *property = new FOFHaloProperties;
  property->setHalos(numberOfHalos, halos, haloCount, haloList);

  string haloSuffix;
  if(applyLC)
    haloSuffix = LC_HALO_SUFFIX;
  else
    haloSuffix = STATIC_HALO_SUFFIX;
  string hpOutBase = create_outName(outBase+"."+haloSuffix, step);

  property->setParameters(hpOutBase, rL, oL, mp, bb);
  property->setParticles(xx, yy, zz, vx, vy, vz, mass, phi, tag, mask, status);

  vector<int> *haloCenter = new vector<int>;
  property->FOFHaloCenterMinimumPotential(haloCenter);

  vector<POSVEL_T> *haloMass = new vector<POSVEL_T>;
  property->FOFHaloMass(haloMass);

  vector<POSVEL_T> *hvx = new vector<POSVEL_T>;
  vector<POSVEL_T> *hvy = new vector<POSVEL_T>;
  vector<POSVEL_T> *hvz = new vector<POSVEL_T>;
  property->FOFVelocity(hvx, hvy, hvz);
  //MPI_Barrier(MPI_COMM_WORLD);


  //OUTPUT
  FILE *hpOutFile = fopen(create_outName(hpOutBase,Partition::getMyProc()).c_str(),"wb");
  int inLC=1;
  for(int i=0; i<numberOfHalos; i++) {
    inLC = 1;
    if(applyLC)
      inLC = lc->LiesInShell( (*xx)[(*haloCenter)[i]],
			     (*yy)[(*haloCenter)[i]],
			     (*zz)[(*haloCenter)[i]] );
    if(inLC) {
      fwrite(&(*xx)[(*haloCenter)[i]], sizeof(POSVEL_T), 1, hpOutFile);
      fwrite(&(*hvx)[i], sizeof(POSVEL_T), 1, hpOutFile);
      
      fwrite(&(*yy)[(*haloCenter)[i]], sizeof(POSVEL_T), 1, hpOutFile);
      fwrite(&(*hvy)[i], sizeof(POSVEL_T), 1, hpOutFile);
      
      fwrite(&(*zz)[(*haloCenter)[i]], sizeof(POSVEL_T), 1, hpOutFile);
      fwrite(&(*hvz)[i], sizeof(POSVEL_T), 1, hpOutFile);
      
      fwrite(&(*haloMass)[i], sizeof(POSVEL_T), 1, hpOutFile);
      
      fwrite(&(*tag)[halos[i]], sizeof(ID_T), 1, hpOutFile);
    }
  }
  fclose(hpOutFile);
  //MPI_Barrier(MPI_COMM_WORLD);


  //CLEANUP
  delete hvx;
  delete hvy;
  delete hvz;
  delete haloMass;
  delete haloCenter;
  delete property;
  delete haloFinder;
  deleteVectors(&xx,&yy,&zz,&vx,&vy,&vz,&mass,&phi,&tag,&mask,&status);
  delete hexchange;

  //MPI_Barrier(MPI_COMM_WORLD);

  return;
}



void staticSkewers(Particles & particles,
		   ParallelSkewerSet & zskewers, 
		   string outBase, 
		   float anow, 
		   int step) {
  int Npla = particles.Np_local_alive();
  vector<POSVEL_T> *xx, *yy, *zz, *vx, *vy, *vz, *phi;
  vector<POSVEL_T> *mass;
  vector<ID_T> *tag;
  vector<MASK_T> *mask;
  vector<STATUS_T> *status;
  allocReserveVectors(&xx, &yy, &zz, &vx, &vy, &vz, &mass, &phi, 
		      &tag, &mask, &status, Npla);
  particles.copyAliveIntoVectors(xx, yy, zz, vx, vy, vz, phi, tag, mask, anow);
  int np = xx->size();
  POSVEL_T *xxa = &(*xx)[0];
  POSVEL_T *yya = &(*yy)[0];
  POSVEL_T *zza = &(*zz)[0];
  POSVEL_T *vxa = &(*vx)[0];
  POSVEL_T *vya = &(*vy)[0];
  POSVEL_T *vza = &(*vz)[0];

  zskewers.PopulateLocal(np, xxa, yya, zza, vxa, vya, vza);

  deleteVectors(&xx,&yy,&zz,&vx,&vy,&vz,&mass,&phi,&tag,&mask,&status);

  string outName = create_outName(create_outName(outBase+"."+STATIC_SKEWER_SUFFIX, step), Partition::getMyProc());
  zskewers.WriteLocalSkewers(outName.c_str());
  zskewers.ClearSkewers();

  return;
}



void lightconeSkewers(Particles & particles,
		      ConvergentSkewerSet & lcskewers, 
		      string outBase, 
		      float anow, 
		      int step,
		      LightCone lc) {
  int Npla = particles.Np_local_alive();
  vector<POSVEL_T> *xx, *yy, *zz, *vx, *vy, *vz, *phi;
  vector<POSVEL_T> *mass;
  vector<ID_T> *tag;
  vector<MASK_T> *mask;
  vector<STATUS_T> *status;
  allocReserveVectors(&xx, &yy, &zz, &vx, &vy, &vz, &mass,
		      &phi, &tag, &mask, &status, Npla);
  particles.copyAliveIntoVectorsLC(&lc, xx, yy, zz, vx, vy, vz, phi, tag, mask, 
				   anow);
  int np = xx->size();
  POSVEL_T *xxa = &(*xx)[0];
  POSVEL_T *yya = &(*yy)[0];
  POSVEL_T *zza = &(*zz)[0];
  POSVEL_T *vxa = &(*vx)[0];
  POSVEL_T *vya = &(*vy)[0];
  POSVEL_T *vza = &(*vz)[0];

  lcskewers.PopulateLocal(np, xxa, yya, zza, vxa, vya, vza);

  deleteVectors(&xx,&yy,&zz,&vx,&vy,&vz,&mass,&phi,&tag,&mask,&status);

  //string outName = create_outName(create_outName(outBase+"."+LC_SKEWER_SUFFIX, step), Partition::getMyProc());
  //lcskewers.WriteLocalSkewers(outName.c_str());
  //lcskewers.ClearSkewers();

  return;
}



void refreshParticles(Particles & particles, float anow) {

  //SET UP EXCHANGE
  float oL = Domain::oL();
  float rL = Domain::rL();
  ParticleExchange *rexchange = new ParticleExchange;
  rexchange->setParameters(rL, oL);
  rexchange->initialize();


  //FIGURE OUT HOW MUCH TO RESERVE
  float rla[DIMENSION];
  Domain::rL_local_alive(rla);
  float rlt[DIMENSION];
  for(int i=0; i<DIMENSION; i++)
    rlt[i] = rla[i] + 2.0*oL;
  int Npla = particles.Np_local_alive();
  int nReserve = static_cast<int>( ceil( FUDGE_REFRESH_OL*(rlt[0]/rla[0])*(rlt[1]/rla[1])*(rlt[2]/rla[2])*Npla ) );


  //GRAB ALIVE PARTICLES
  vector<POSVEL_T> *xx, *yy, *zz, *vx, *vy, *vz, *phi;
  vector<POSVEL_T> *mass;
  vector<ID_T> *tag;
  vector<MASK_T> *mask;
  vector<STATUS_T> *status;
  allocReserveVectors(&xx, &yy, &zz, &vx, &vy, &vz, &mass, &phi, 
		      //&tag, &mask, &status, nReserve);
		      &tag, &mask, &status, 0);
  mass->reserve(nReserve);
  status->reserve(nReserve);
  particles.copyAliveIntoVectors(xx, yy, zz, vx, vy, vz, phi, tag, mask, anow,
				 nReserve, 1);
  for(int i=0; i<xx->size(); i++)
    mass->push_back(1.0);
  rexchange->setParticles(xx, yy, zz, vx, vy, vz, mass, phi, tag, mask, status);


  //DROP PARTICLES CLASS MEMORY
  particles.dropParticles();
  //MPI_Barrier(MPI_COMM_WORLD);


  //EXCHANGE
  rexchange->exchangeParticles();
  //MPI_Barrier(MPI_COMM_WORLD);


  printf("REFRESH PXCH rank = %d Npla = %d nReserve = %d size = %d\n",
	 Partition::getMyProc(), Npla, nReserve, xx->size());
  fflush(stdout);
  //MPI_Barrier(MPI_COMM_WORLD);


  //COPY REFRESHED PARTICLES BACK INTO CLASS
  particles.copyAllFromVectors(&xx, &yy, &zz, &vx, &vy, &vz, &phi, &tag, 
			       &mask, anow, 1);


  //CLEANUP MEMORY
  deleteVectors(&xx,&yy,&zz,&vx,&vy,&vz,&mass,&phi,&tag,&mask,&status);
  delete rexchange;


  //RESORT PARTICLES
  particles.resortParticles();
  //MPI_Barrier(MPI_COMM_WORLD);


  return;
}
*/



/*
void timingStats(double t, const char *comment) {
  double tmax, tmin, tavg;
  int rank;
  double numranks;

  rank = Partition::getMyProc();
  numranks = 1.0*Partition::getNumProc();

  MPI_Allreduce(&t, &tmax, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
  MPI_Allreduce(&t, &tmin, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
  MPI_Allreduce(&t, &tavg, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  tavg /= numranks;
  
  if(rank==0) {
    printf("%s max %.3es avg %.3es min %.3es\n",
	   comment, tmax, tavg, tmin);
    fflush(stdout);
  }

  return;
}
*/
