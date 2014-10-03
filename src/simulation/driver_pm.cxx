#include "mc3.h"



//MAIN
int main(int argc, char* argv[])
{
  int step0 = 0;

  SimpleTimings::TimerRef t_total = SimpleTimings::getTimer("total");
  SimpleTimings::TimerRef t_init = SimpleTimings::getTimer("init");
  SimpleTimings::TimerRef t_stepr = SimpleTimings::getTimer("stepr");
  SimpleTimings::TimerRef t_step = SimpleTimings::getTimer("step");
  SimpleTimings::TimerRef t_xtra = SimpleTimings::getTimer("xtra");
  SimpleTimings::TimerRef t_sort = SimpleTimings::getTimer("sort");
  SimpleTimings::TimerRef t_map1 = SimpleTimings::getTimer("map1");

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
    fprintf(stderr,"-w                   : use white noise initializer\n");
    fprintf(stderr,"-I                   : use MPI IO for restart files and pseudo-outputs\n");
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
  SimpleTimings::startTimer(t_total);
  SimpleTimings::startTimer(t_init);


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
    if(options.mpiio()) {
      particles.readRestart( create_outName( inBase + "." + MPI_RESTART_SUFFIX, step0).c_str() );
    } else {
      particles.readRestart( create_outName( create_outName( inBase + "." + RESTART_SUFFIX, step0), rank).c_str() );
    }
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

    writePk(solver, outBase + "." + PK_SUFFIX + ".ini");
    MPI_Barrier(MPI_COMM_WORLD);

    if(rank==0) printf("\n");
    MPI_Barrier(MPI_COMM_WORLD);
  }


  //DONE WITH INITIALIZATION
  SimpleTimings::stopTimerStats(t_init);
  if(rank==0) printf("\n");


  vector<int> Nplav;
  Nplav.reserve(ts.nsteps()+1);
  Nplav.push_back(particles.Np_local_alive());


  //TIMESTEPPER
  SimpleTimings::startTimer(t_stepr);

  //if restart, get timestepper variables up to speed
  for(int step = 0; step < step0; step++) {
    ts.advanceFullStep();
    extras->setStep(step, ts.aa());
  }
  MPI_Barrier(MPI_COMM_WORLD);


  //actual timestepping
  for(int step = step0; step < ts.nsteps(); step++) {
    SimpleTimings::startTimer(t_step);

    if(rank==0) {
      printf("STEP %d, pp = %f, a = %f, z = %f\n",
	     step, ts.pp(), ts.aa(), ts.zz());
      fflush(stdout);
    }
    MPI_Barrier(MPI_COMM_WORLD);


    //HALF STREAM
    SimpleTimings::startTimer(t_map1);
    particles.map1(ts.pp(), ts.tau2(), ts.adot());
    SimpleTimings::stopTimerStats(t_map1);
    MPI_Barrier(MPI_COMM_WORLD);


    //UPDATE TIME
    ts.advanceHalfStep();
    MPI_Barrier(MPI_COMM_WORLD);


    //POISSON FORWARD
    map2_poisson_forward(particles, solver, rho_arr, fft_rho_arr);
    MPI_Barrier(MPI_COMM_WORLD);


    //CHECKSUM DENSITY GRID
    rho_local_alive = sum_rho_alive(rho_arr);
    MPI_Allreduce(&rho_local_alive, &rho_global_alive, 1, 
		  MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    MPI_Barrier(MPI_COMM_WORLD);


    //POISSON BACKWARD GRADIENT
    map2_poisson_backward_gradient(particles, solver,
				   grad_phi_arr, fft_grad_phi_arr,
				   gexchange, ts, 1.0);
    MPI_Barrier(MPI_COMM_WORLD);


    //UPDATE TIME
    ts.advanceHalfStep();
    MPI_Barrier(MPI_COMM_WORLD);


    //HALF STREAM
    SimpleTimings::startTimer(t_map1);
    particles.map1(ts.pp(), ts.tau2(), ts.adot());
    SimpleTimings::stopTimerStats(t_map1);
    MPI_Barrier(MPI_COMM_WORLD);


    //CHECKSUM PARTICLES
    Np_local_alive = particles.Np_local_alive();
    Nplav.push_back(Np_local_alive);
    MPI_Allreduce(&Np_local_alive, &Np_global_alive, 1, 
		  MPI_LONG_LONG, MPI_SUM, MPI_COMM_WORLD);
    if(rank==0) {
      printf( "total alive density   = %f\n",rho_global_alive);
      cout << "total alive particles = " << Np_global_alive << endl;
      cout.flush();
    }
    MPI_Barrier(MPI_COMM_WORLD);




    //EXTRA STUFF
    SimpleTimings::startTimer(t_xtra);
    extras->setStep(step, ts.aa());


    if(extras->extrasStep()) {
      if(rank==0) {
	printf("EXTRAS: ");
	if(extras->staticStep())printf("(static output) ");
	if(extras->lcStep())printf("(light cone update) ");
	if(extras->aliveStep())printf("(alive particle output) ");
	if(extras->restartStep())printf("(restart dump) ");
	if(extras->pkStep())printf("(power spectrum) ");
	if(extras->refreshStep())printf("(overload particle refresh)");
	printf("\n");
	fflush(stdout);
      }
    }
    MPI_Barrier(MPI_COMM_WORLD);


    //FORWARD FFT SOLVE
    if(extras->fftfStep())
      map2_poisson_forward(particles, solver, rho_arr, fft_rho_arr);
    MPI_Barrier(MPI_COMM_WORLD);


    //P(k) DUMP
    if(extras->pkStep()) {
      writePk(solver, create_outName(outBase + "." + PK_SUFFIX, step) );
      MPI_Barrier(MPI_COMM_WORLD);
    }


    //BACKWARD POTENTIAL CALCULATION
    if(extras->fftbpotStep())
      map2_poisson_backward_potential(particles, solver, 
				      phi_arr, fft_phi_arr,
				      gexchange);
    MPI_Barrier(MPI_COMM_WORLD);


    //DROP FFT MEMORY
    int fftMemDropped = 0;
    if(extras->particleStep()) {
      fftMemDropped = 1;
      delete solver;
      poisson_free(&fft_rho_arr, &fft_grad_phi_arr, &fft_phi_arr);
      MPI_Barrier(MPI_COMM_WORLD);
    }


    if(extras->particleStep())
      extras->particleExtras(particles, indat, outBase, &Nplav);


    //RE-ALLOC FFT MEMORY
    if(fftMemDropped) {
      particles.shoveParticles();
      solver = new SolverDiscrete(MPI_COMM_WORLD, ng);
      poisson_alloc(&fft_rho_arr, &fft_grad_phi_arr, &fft_phi_arr);
      MPI_Barrier(MPI_COMM_WORLD);
    }

    SimpleTimings::stopTimerStats(t_xtra);
    SimpleTimings::stopTimerStats(t_step);

    if(rank==0) {
      printf("\n");
      fflush(stdout);
    }

  } // end timestepper 


  SimpleTimings::stopTimerStats(t_stepr);


  //OUTPUT REST OF LIGHT CONE SKEWERS ACCUMULATED
  if(extras->lightconeQ() && extras->skewerQ()) {
    string outName = create_outName(create_outName(outBase+"."+LC_SKEWER_SUFFIX, ts.nsteps()), Partition::getMyProc());
    (extras->lcskewers())->WriteLocalSkewers(outName.c_str());
    (extras->lcskewers())->ClearSkewers();
  }


  //P(k) FINAL
  if(rank==0) {
    printf("P(k) final\n");
    fflush(stdout);
  }
  MPI_Barrier(MPI_COMM_WORLD);
  map2_poisson_forward(particles, solver, rho_arr, fft_rho_arr);
  MPI_Barrier(MPI_COMM_WORLD);
  writePk(solver, outBase + "." + PK_SUFFIX + ".fin");  
  MPI_Barrier(MPI_COMM_WORLD);
  if(rank==0)
    printf("\n");
  MPI_Barrier(MPI_COMM_WORLD);

  
  //GRID OUTPUT
  if(extras->gridQ()) {
    //CIC ALREADY DONE FOR P(k)
    output_array_alive(rho_arr,create_outName(create_outName(outBase+"."+GRID_SUFFIX,ts.nsteps()),rank).c_str());
  }
  MPI_Barrier(MPI_COMM_WORLD);


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
