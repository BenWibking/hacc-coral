#include "mc3.h"



//MAIN
int main(int argc, char* argv[])
{
  int step0 = 0;


  if(argc < 6) {
    fprintf(stderr,"USAGE: mc3 <indat> <inBase|tfName> <outBase> <INIT|RECORD|BLOCK|COSMO|RESTART> <ROUND_ROBIN|ALL_TO_ALL|ONE_TO_ONE|restart_step>\n");
    fprintf(stderr,"-a <aliveDumpName>   : alive particle dumps\n");
    fprintf(stderr,"-r <restartDumpName> : restart particle dumps\n");
    fprintf(stderr,"-f <refreshStepName> : steps for particle refresh\n");
    fprintf(stderr,"-o <analysisdat>     : config file for analysis\n");
    fprintf(stderr,"-s <staticDumpName>  : static time analysis dumps\n");
    fprintf(stderr,"-l <LCUpdateName>    : lightcone time updates\n");
    fprintf(stderr,"-h                   : halo outputs\n");
    fprintf(stderr,"-z                   : skewerz\n");
    fprintf(stderr,"-g                   : final grid output\n");
    fprintf(stderr,"-m                   : initialize MPI_Alltoall\n");
    fprintf(stderr,"-p <pkDumpName>      : P(k) dumps (+ initial, final, restarts)\n");
    fprintf(stderr,"-t <NXxNYxNZ>        : use 3D topology NX by NY by NZ\n");
    exit(-1);
  }


  //sort command line options
  MC3Options options(argc, argv);


  int argvi = optind;
  string indatName = argv[argvi++];
  string inBase = argv[argvi++];
  string outBase = argv[argvi++];
  string dataType = argv[argvi++];
  string distributeType = argv[argvi++];


  //starting from restart file
  int restartQ = 0;
  if( dataType == "RESTART") {
    restartQ = 1;
    step0 = atoi(distributeType.c_str());
  }


  //need to convert Mpc to Mpc/h for true cosmo format
  int cosmoFormatQ=0;
  if( dataType == "COSMO") {
    cosmoFormatQ=1;
    dataType = "RECORD";
  }


  //INITIALIZE MPI
  MPI_Init(&argc, &argv);


  //INITIALIZE TOPOLOGY
  int tmpDims[DIMENSION];
  if(options.topologyQ()) {
    char *tmpDimsStr = (char *)options.topologyString().c_str();
    char *tmpDimsTok;
    tmpDimsTok = strtok(tmpDimsStr,"x");
    tmpDims[0] = atoi(tmpDimsTok);
    tmpDimsTok = strtok(NULL,"x");
    tmpDims[1] = atoi(tmpDimsTok);
    tmpDimsTok = strtok(NULL,"x");
    tmpDims[2] = atoi(tmpDimsTok);
    int nnodes;
    MPI_Comm_size(MPI_COMM_WORLD, &nnodes);
    MY_Dims_init_3D(nnodes, DIMENSION, tmpDims);
  }
  Partition::initialize();
  int numranks = Partition::getNumProc();
  int rank = Partition::getMyProc();


  //READ INDAT FILE
  Basedata indat( indatName.c_str() );


  //INITIALIZE GEOMETRY
  Domain::initialize(indat);


  //
  MC3Extras *extras = new MC3Extras(options, indat);


  //start some timers
  SimpleTimings::TimerRef t_total = SimpleTimings::getTimer("total");
  SimpleTimings::startTimer(t_total);


  //(OPTIONALLY) INITIALIZE MPI_ALLTOALL
  if(options.initialAlltoallQ()) {
    if(rank==0) {
      printf("\nStarting MPI_Alltoall initialization\n");
      fflush(stdout);
    }
    initial_all_to_all(options.initialAlltoallQ());
    MPI_Barrier(MPI_COMM_WORLD);
  }


  //INSTANTIATE PARTICLES
  Particles particles(indat, options);

    
  if(!restartQ) {
    //START FROM THE BEGINNING
    loadParticles(indat, particles, 
		  inBase, dataType, distributeType, 
		  cosmoFormatQ, options);
  } else {
    //RESTARTING
    particles.readRestart( create_outName( create_outName( inBase + "." + RESTART_SUFFIX, step0), rank).c_str() );
    step0++;
  }
  MPI_Barrier(MPI_COMM_WORLD);


  //MOVE PARTICLES TO CELL BEFORE ALLOCATING FFT MEMORY
  particles.shoveParticles();
  MPI_Barrier(MPI_COMM_WORLD);


  //LOCAL GRID
  GRID_T *rho_arr = particles.field();
  GRID_T *grad_phi_arr = rho_arr;
  GRID_T *phi_arr = rho_arr;


  //LOCAL COPIES OF GEOMETRIC INFO
  int ngla[DIMENSION], nglt[DIMENSION];
  Domain::ng_local_alive(ngla);
  Domain::ng_local_total(nglt);
  int Ngla = Domain::Ng_local_alive();
  int Nglt = Domain::Ng_local_total();
  int ngo = Domain::ng_overload();
  int ng = Domain::ng();


  //INITIALIZE GRID EXCHANGE
  GridExchange gexchange(nglt, ngo, ngo+1);


  //ALLOC POISSON SOLVER AND BUFFERS
  SolverDiscrete *solver = new SolverDiscrete(MPI_COMM_WORLD, ng);
  COMPLEX_T *fft_rho_arr, *fft_grad_phi_arr, *fft_phi_arr;
  poisson_alloc(&fft_rho_arr, &fft_grad_phi_arr, &fft_phi_arr);
  MPI_Barrier(MPI_COMM_WORLD);


  //TIMESTEPPER VARIABLES
  TimeStepper ts(indat.alpha(), indat.ain(), indat.afin(),
		 indat.nsteps(), indat.omegatot() );
  int64_t Np_local_alive, Np_global_alive;
  double rho_local_alive, rho_global_alive;


  MPI_Barrier(MPI_COMM_WORLD);


  //P(k) INITIAL
  if(!restartQ) {
    if(rank==0) {
      printf("P(k) initial\n");
      fflush(stdout);
    }
    MPI_Barrier(MPI_COMM_WORLD);

    map2_poisson_forward(particles, solver, rho_arr, fft_rho_arr);
    MPI_Barrier(MPI_COMM_WORLD);

    //writePk(solver, outBase + "." + PK_SUFFIX + ".ini");
    writePk(solver, outBase);
    MPI_Barrier(MPI_COMM_WORLD);

    if(rank==0) printf("\n");
    MPI_Barrier(MPI_COMM_WORLD);
  }

  delete solver;
  poisson_free(&fft_rho_arr, &fft_grad_phi_arr, &fft_phi_arr);
  MPI_Barrier(MPI_COMM_WORLD);


  SimpleTimings::stopTimer(t_total);
  SimpleTimings::accumStats();


  // Shut down MPI
  Partition::finalize();
  MPI_Finalize();

  return 0;
}
