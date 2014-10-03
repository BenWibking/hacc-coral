#include "MC3Extras.h"

MC3Extras::MC3Extras(MC3Options & options, Basedata & indat) :
  m_gridQ(0),
  m_restartDumpQ(0),
  m_pkDumpQ(0),
  m_aliveDumpQ(0),
  m_refreshQ(0),
  m_analysisQ(0),
  m_staticQ(0),
  m_lightconeQ(0),
  m_initialAlltoallQ(false),
  m_mpiio(0),
  m_skewerQ(0),
  m_haloQ(0),
  m_restartDumps(NULL),
  m_pkDumps(NULL),
  m_aliveDumps(NULL),
  m_refreshSteps(NULL),
  m_staticDumps(NULL),
  m_LCUpdates(NULL),
  m_analysisdat(NULL),
  m_lc(NULL),
  m_aa_last(0.0),
  m_lc_started(0),
  m_staticStep(0),
  m_lcStep(0),
  m_aliveStep(0),
  m_restartStep(0),
  m_pkStep(0),
  m_refreshStep(0),
  m_step(0),
  m_aa(0.0)
{
  m_gridQ = options.gridQ();
  m_initialAlltoallQ = options.initialAlltoallQ();
  m_haloQ = options.haloQ();
  m_skewerQ = options.skewerQ();
  m_restartDumpQ = options.restartDumpQ();
  m_restartDumpName = options.restartDumpName();
  m_pkDumpQ = options.pkDumpQ();
  m_pkDumpName = options.pkDumpName();
  m_aliveDumpQ = options.aliveDumpQ();
  m_aliveDumpName = options.aliveDumpName();
  m_analysisQ = options.analysisQ();
  m_analysisdatName = options.analysisdatName();
  m_staticQ = options.staticQ();
  m_staticDumpName = options.staticDumpName();
  m_lightconeQ = options.lightconeQ();
  m_LCUpdateName = options.LCUpdateName();
  m_refreshQ = options.refreshQ();
  m_refreshName = options.refreshName();
  m_mpiio = options.mpiio();
  
  if(m_restartDumpQ)
    m_restartDumps = readInts(m_restartDumpName);

  if(m_pkDumpQ)
    m_pkDumps = readInts(m_pkDumpName);

  if(m_aliveDumpQ)
    m_aliveDumps = readInts(m_aliveDumpName);

  if(m_refreshQ)
    m_refreshSteps = readInts(m_refreshName);

  if(m_analysisQ)
    m_analysisdat = new Halodata(m_analysisdatName.c_str());

  if(m_staticQ) {
    m_staticDumps = readInts(m_staticDumpName);
    if(!m_analysisQ) {
      fprintf(stderr,"ERROR: static output needs an analysis config file.\n");
      exit(-1);
    }
  }

  if(m_lightconeQ) {
    m_LCUpdates = readInts(m_LCUpdateName);
    if(!m_analysisQ) {
      fprintf(stderr,"ERROR: lighcone output needs an analysis config file.\n");
      exit(-1);
    }
  }
  
  initializeSkewers(indat);
}



vector<int>* MC3Extras::readInts(string inName) {
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



int MC3Extras::intInVec(int n, vector<int> *iv) {
  int ret=0;

  if(iv)
    for(int i=0; i<iv->size(); i++)
      if( (*iv)[i] == n )
	ret = 1;

  return ret;
}



void MC3Extras::initializeSkewers(Basedata & indat) {
  //SET UP SKEWERS
  srand48(indat.iseed());


  //STATIC TIME SKEWERS
  if(m_staticQ) {
    m_zskewers.InitializeMaster(m_analysisdat->Nskewers(),
				m_analysisdat->Npixels(),
				m_analysisdat->h(),
				Domain::rL(),
				m_analysisdat->Nmesh(),
				0);
    m_zskewers.InitializeSlaves(0);
  }


  //LIGHTCONE SKEWERS
  m_lc = new LightCone(indat.omegatot(), 0, 0, 0, FUDGE_Z_IN*indat.zin());
  float rL = Domain::rL();
  float zbox = m_lc->rofa(indat.afin()) + 0.5*rL;
  m_lc->SetOrigin(0.5*rL, 0.5*rL, 0.5*rL-zbox);

  m_aa_last = indat.ain();
  m_lc_started = 0;

  if(m_lightconeQ) {
    m_lcskewers.InitializeMaster(m_analysisdat->Nskewers(),
				 m_analysisdat->Npixels(),
				 m_analysisdat->h(),
				 rL,
				 zbox,
				 m_analysisdat->Nmesh(),
				 0);
    m_lcskewers.InitializeSlaves(0);
  }

  return;
}



void MC3Extras::setStep(int step, float aa) {
  m_step = step;
  m_aa = aa;

  int SorH = m_skewerQ || m_haloQ;

  m_staticStep = (m_staticQ && intInVec(step, m_staticDumps) && SorH);
  m_lcStep = (m_lightconeQ && intInVec(step, m_LCUpdates) && SorH);
  m_aliveStep = (m_aliveDumpQ && intInVec(step, m_aliveDumps));
  m_restartStep = (m_restartDumpQ && intInVec(step, m_restartDumps));
  m_pkStep = (m_pkDumpQ && intInVec(step, m_pkDumps));
  m_refreshStep = (m_refreshQ && intInVec(step,m_refreshSteps));

  //ALL STEPS WHERE WE OUTPUT P(k)
  m_pkStep = (m_pkStep || m_restartStep);

  //ALL STEPS WHERE WE XFER PARTICLES TO EXTERNAL VECTORS
  m_refreshStep = (m_refreshStep || 
		   (SorH && (m_staticStep || m_lcStep)));

  if(m_lcStep) {
    m_lc->DefineShell(m_aa_last, aa);
    m_aa_last = aa;
    m_lc_started = 1;
  }

  m_extrasStep = (m_staticStep || m_lcStep || m_aliveStep || m_restartStep || 
		  m_pkStep || m_refreshStep);
  m_fftfStep = (m_staticStep || m_lcStep || m_aliveStep || m_pkStep);
  m_fftbpotStep = (m_staticStep || m_lcStep || m_aliveStep);
  m_particleStep = (m_staticStep || m_lcStep || m_aliveStep || m_restartStep ||
		    m_refreshStep);

  return;
}



string MC3Extras::create_outName(string outBase, int rank) {
  ostringstream outName;
  outName << outBase << "." << rank;
  return outName.str();
}



void MC3Extras::allocReserveVectors(int Np) {
  m_xx = new vector<POSVEL_T>;
  m_yy = new vector<POSVEL_T>;
  m_zz = new vector<POSVEL_T>;
  m_vx = new vector<POSVEL_T>;
  m_vy = new vector<POSVEL_T>;
  m_vz = new vector<POSVEL_T>;
  m_phi = new vector<POSVEL_T>;
  m_id = new vector<ID_T>;
  m_mask = new vector<MASK_T>;

  m_mass = new vector<POSVEL_T>;
  m_status = new vector<STATUS_T>;

  reserveVectors(Np);

  return;
}



void MC3Extras::reserveVectors(int Np) {
  m_xx->reserve(Np);
  m_yy->reserve(Np);
  m_zz->reserve(Np);
  m_vx->reserve(Np);
  m_vy->reserve(Np);
  m_vz->reserve(Np);
  m_phi->reserve(Np);
  m_id->reserve(Np);
  m_mask->reserve(Np);

  m_mass->reserve(Np);
  m_status->reserve(Np);

  return;
}



void MC3Extras::deleteVectors() {
  if(m_xx) delete m_xx;
  if(m_yy) delete m_yy;
  if(m_zz) delete m_zz;
  if(m_vx) delete m_vx;
  if(m_vy) delete m_vy;
  if(m_vz) delete m_vz;
  if(m_phi) delete m_phi;
  if(m_id) delete m_id;
  if(m_mask) delete m_mask;

  if(m_mass) delete m_mass;
  if(m_status) delete m_status;

  m_xx = NULL;
  m_yy = NULL;
  m_zz = NULL;
  m_vx = NULL;
  m_vy = NULL;
  m_vz = NULL;
  m_phi = NULL;
  m_id = NULL;
  m_mask = NULL;

  m_mass = NULL;
  m_status = NULL;

  return;
}



#define MAX(a, b) (a > b ? a : b)

int MC3Extras::calcNReserve(vector<int> *Nplav) {
  //FIGURE OUT HOW MUCH TO RESERVE

  int nReserve;

  float oL;
  if(m_haloQ && (m_staticStep || m_lcStep))
    oL = MAX(m_analysisdat->oL(), Domain::oL());
  else
    oL = Domain::oL();

  float rla[DIMENSION];
  Domain::rL_local_alive(rla);
  float rlt[DIMENSION];
  for(int i=0; i<DIMENSION; i++)
    rlt[i] = rla[i] + 2.0*oL;

  int nNpla = Nplav->size();
  int Npla = (*Nplav)[nNpla-1];

  nReserve = static_cast<int>( ceil( FUDGE_OL_RESERVE*(rlt[0]/rla[0])*(rlt[1]/rla[1])*(rlt[2]/rla[2])*Npla ) );

  return nReserve;
}



void MC3Extras::particleExtras(Particles & particles,
				Basedata & indat,
				string outBase,
				vector<int> *Nplav) {
  int numranks = Partition::getNumProc();
  int rank = Partition::getMyProc();
  int Npla = particles.Np_local_alive();
  int step = m_step;
  float aa = m_aa;

  //ALIVE PARTICLE DUMP
  if(m_aliveStep) {
    if(m_mpiio) {
      particles.writeAliveHCosmo( create_outName( outBase + "." + MPI_ALIVE_SUFFIX, step).c_str(), aa);
    } else {
      particles.writeAliveHCosmo( create_outName( create_outName( outBase + "." + ALIVE_SUFFIX, step), rank).c_str(), aa);
    }
    //MPI_Barrier(MPI_COMM_WORLD);
  }
  
  //MOVE PARTICLE INFO INTO VECTORS FOR SKEWER, HALO, REFRESH
  if(m_refreshStep) {
    allocReserveVectors(0);
    int nReserve = calcNReserve(Nplav);
    particles.copyAliveIntoVectors(m_xx, m_yy, m_zz, m_vx, m_vy, m_vz, m_phi,
				   m_id, m_mask, m_aa, nReserve, 1);
    //particles is now empty
    //vectors have alive particles
    //vectors have space reserved for halo and refresh overloading

    //SKEWERS
    if(m_skewerQ) {
      
      if(m_staticStep) {
	staticSkewers(outBase);
	//MPI_Barrier(MPI_COMM_WORLD);
      }
    
      if(m_lcStep) {
	lightconeSkewers();
	//MPI_Barrier(MPI_COMM_WORLD);
      }
    }
    
    //HALOS
    if(m_haloQ && (m_staticStep || m_lcStep)) {
      findHalos(indat, outBase);
      //MPI_Barrier(MPI_COMM_WORLD);
    }

    //REFRESH
    if(m_refreshStep) {
      refreshParticles(particles);
      //MPI_Barrier(MPI_COMM_WORLD);
    }    
  }

  //RESTART PARTICLE DUMP
  if(m_restartStep) {
    if(m_mpiio) {
      particles.writeRestart( create_outName( outBase + "." + MPI_RESTART_SUFFIX, step).c_str() );
    } else {
      particles.writeRestart( create_outName( create_outName( outBase + "." + RESTART_SUFFIX, step), rank).c_str() );
    }
    //MPI_Barrier(MPI_COMM_WORLD);

    //OUTPUT LIGHT CONE SKEWERS ACCUMULATED SO FAR
    if(m_lightconeQ && m_skewerQ && m_lc_started) {
      string outName = create_outName(create_outName(outBase+"."+LC_SKEWER_SUFFIX, step), Partition::getMyProc());
      m_lcskewers.WriteLocalSkewers(outName.c_str());
      m_lcskewers.ClearSkewers();
      //MPI_Barrier(MPI_COMM_WORLD);
    }
  }

  return;
}



void MC3Extras::staticSkewers(string outBase) {
  int np = m_xx->size();
  POSVEL_T *xxa = &(*m_xx)[0];
  POSVEL_T *yya = &(*m_yy)[0];
  POSVEL_T *zza = &(*m_zz)[0];
  POSVEL_T *vxa = &(*m_vx)[0];
  POSVEL_T *vya = &(*m_vy)[0];
  POSVEL_T *vza = &(*m_vz)[0];

  m_zskewers.PopulateLocal(np, xxa, yya, zza, vxa, vya, vza);

  string outName = create_outName(create_outName(outBase+"."+STATIC_SKEWER_SUFFIX, m_step), Partition::getMyProc());
  m_zskewers.WriteLocalSkewers(outName.c_str());
  m_zskewers.ClearSkewers();

  return;
}



void MC3Extras::lightconeSkewers() {
  POSVEL_T *xxa = &(*m_xx)[0];
  POSVEL_T *yya = &(*m_yy)[0];
  POSVEL_T *zza = &(*m_zz)[0];
  POSVEL_T *vxa = &(*m_vx)[0];
  POSVEL_T *vya = &(*m_vy)[0];
  POSVEL_T *vza = &(*m_vz)[0];

  int Npla = m_xx->size();

  vector<unsigned int> *lci = new vector<unsigned int>;
  lci->reserve(Npla);
  for(int i=0; i<Npla; i++)
    if(m_lc->LiesInShell((*m_xx)[i], (*m_yy)[i], (*m_zz)[i]))
      lci->push_back(i);
  int np = lci->size();

  m_lcskewers.PopulateLocal(np, &(*lci)[0], xxa, yya, zza, vxa, vya, vza);

  delete lci;

  return;
}



/*
void MC3Extras::staticSkewers(Particles & particles, string outBase) {
  POSVEL_T *xxa = particles.xArr();
  POSVEL_T *yya = particles.yArr();
  POSVEL_T *zza = particles.zArr();
  POSVEL_T *vxa = particles.vxArr();
  POSVEL_T *vya = particles.vyArr();
  POSVEL_T *vza = particles.vzArr();

  vector<unsigned int> *ai = particles.aliveIndices();
  int Npla = ai->size();

  particles.coords_local2global(m_aa);
  m_zskewers.PopulateLocal(Npla, &(*ai)[0], xxa, yya, zza, vxa, vya, vza);
  particles.coords_global2local(m_aa);

  string outName = create_outName(create_outName(outBase+"."+STATIC_SKEWER_SUFFIX, m_step), Partition::getMyProc());
  m_zskewers.WriteLocalSkewers(outName.c_str());
  m_zskewers.ClearSkewers();

  delete ai;

  return;
}



void MC3Extras::lightconeSkewers(Particles & particles) {
  POSVEL_T *xxa = particles.xArr();
  POSVEL_T *yya = particles.yArr();
  POSVEL_T *zza = particles.zArr();
  POSVEL_T *vxa = particles.vxArr();
  POSVEL_T *vya = particles.vyArr();
  POSVEL_T *vza = particles.vzArr();

  vector<unsigned int> *ai = particles.aliveIndices();
  int Npla = ai->size();

  vector<unsigned int> *lci = new vector<unsigned int>;
  lci->reserve(Npla);
  for(int i=0; i<Npla; i++)
    if(m_lc->LiesInShell( xxa[(*ai)[i]], yya[(*ai)[i]], zza[(*ai)[i]] ))
      lci->push_back((*ai)[i]);
  int np = lci->size();

  particles.coords_local2global(m_aa);
  m_lcskewers.PopulateLocal(np, &(*lci)[0], xxa, yya, zza, vxa, vya, vza);
  particles.coords_global2local(m_aa);

  delete ai;
  delete lci;

  return;
}
*/



void MC3Extras::findHalos(Basedata & indat, string outBase) {
  int Npla = m_xx->size();
  int nReserve = m_xx->capacity();
  int step = m_step;


  //INITIALIZE EXCHANGE
  float rL = Domain::rL();
  float oL = m_analysisdat->oL();
  ParticleExchange *hexchange = new ParticleExchange;
  hexchange->setParameters(rL, oL);
  hexchange->initialize();


  //FIGURE OUT MASS OF PARTICLE
  float npf = static_cast<float>(indat.np());
  float omegadm = indat.omegadm();
  float hubble = indat.hubble();
  float deut = indat.deut();
  float omegatot = omegadm + deut/hubble/hubble;
  float mp = 2.77536627e11*rL*rL*rL*omegatot/npf/npf/npf;


  m_mass->reserve(nReserve);
  m_status->reserve(nReserve);
  for(int i=0; i<Npla; i++)
    m_mass->push_back(mp);


  hexchange->setParticles(m_xx, m_yy, m_zz, m_vx, m_vy, m_vz, 
			  m_mass, m_phi, m_id, m_mask, m_status);
  //MPI_Barrier(MPI_COMM_WORLD);


  //EXCHANGE
  hexchange->exchangeParticles();
  //MPI_Barrier(MPI_COMM_WORLD);


#if 0
  printf("HALO    PXCH rank = %d Npla = %d nReserve = %d size = %d\n",
	 Partition::getMyProc(), Npla, nReserve, m_xx->size());
  fflush(stdout);
  //MPI_Barrier(MPI_COMM_WORLD);
#endif


  //FIGURE OTHER FOF PARAMETERS
  float bb = m_analysisdat->bb();
  float minmass = m_analysisdat->minmass();
  int nmin = static_cast<int>(floor(minmass/mp));
  if(nmin < m_analysisdat->nmin())
    nmin = m_analysisdat->nmin();


  //FIND HALOS
  CosmoHaloFinderP *haloFinder = new CosmoHaloFinderP;
  int np = indat.np();
  string hfOutBase = create_outName(outBase+"."+HALO_PARTICLE_SUFFIX,step);
  haloFinder->setParameters(hfOutBase, rL, oL, np, nmin, bb);
  haloFinder->setParticles(m_xx, m_yy, m_zz, m_vx, m_vy, m_vz, 
			   m_phi, m_id, m_mask, m_status);
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
  string hpOutBase = create_outName(outBase+"."+HALO_SUFFIX, step);
  property->setParameters(hpOutBase, rL, oL, bb);
  property->setParticles(m_xx, m_yy, m_zz, m_vx, m_vy, m_vz, 
			 m_mass, m_phi, m_id, m_mask, m_status);
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
  if(m_lcStep) {
    string hpOutBase = create_outName(outBase+"."+LC_HALO_SUFFIX, step);
    FILE *hpOutFile = fopen(create_outName(hpOutBase,Partition::getMyProc()).c_str(),"wb");

    int inLC;
    for(int i=0; i<numberOfHalos; i++) {
      inLC = m_lc->LiesInShell( (*m_xx)[(*haloCenter)[i]],
				(*m_yy)[(*haloCenter)[i]],
				(*m_zz)[(*haloCenter)[i]] );
      if(inLC) {
	fwrite(&(*m_xx)[(*haloCenter)[i]], sizeof(POSVEL_T), 1, hpOutFile);
	fwrite(&(*hvx)[i], sizeof(POSVEL_T), 1, hpOutFile);
	fwrite(&(*m_yy)[(*haloCenter)[i]], sizeof(POSVEL_T), 1, hpOutFile);
	fwrite(&(*hvy)[i], sizeof(POSVEL_T), 1, hpOutFile);
	fwrite(&(*m_zz)[(*haloCenter)[i]], sizeof(POSVEL_T), 1, hpOutFile);
	fwrite(&(*hvz)[i], sizeof(POSVEL_T), 1, hpOutFile);
	fwrite(&(*haloMass)[i], sizeof(POSVEL_T), 1, hpOutFile);
	fwrite(&(*m_id)[halos[i]], sizeof(ID_T), 1, hpOutFile);
      }
    }
    fclose(hpOutFile);
    //MPI_Barrier(MPI_COMM_WORLD);
  }

  if(m_staticStep) {
    string hpOutBase = create_outName(outBase+"."+STATIC_HALO_SUFFIX, step);
    FILE *hpOutFile = fopen(create_outName(hpOutBase,Partition::getMyProc()).c_str(),"wb");

    for(int i=0; i<numberOfHalos; i++) {
      fwrite(&(*m_xx)[(*haloCenter)[i]], sizeof(POSVEL_T), 1, hpOutFile);
      fwrite(&(*hvx)[i], sizeof(POSVEL_T), 1, hpOutFile);
      fwrite(&(*m_yy)[(*haloCenter)[i]], sizeof(POSVEL_T), 1, hpOutFile);
      fwrite(&(*hvy)[i], sizeof(POSVEL_T), 1, hpOutFile);
      fwrite(&(*m_zz)[(*haloCenter)[i]], sizeof(POSVEL_T), 1, hpOutFile);
      fwrite(&(*hvz)[i], sizeof(POSVEL_T), 1, hpOutFile);
      fwrite(&(*haloMass)[i], sizeof(POSVEL_T), 1, hpOutFile);
      fwrite(&(*m_id)[halos[i]], sizeof(ID_T), 1, hpOutFile);
    }
    fclose(hpOutFile);
    //MPI_Barrier(MPI_COMM_WORLD);
  }

  //CLEANUP
  delete hvx;
  delete hvy;
  delete hvz;
  delete haloMass;
  delete haloCenter;
  delete property;
  delete haloFinder;
  delete hexchange;

  m_xx->resize(Npla);
  m_yy->resize(Npla);
  m_zz->resize(Npla);
  m_vx->resize(Npla);
  m_vy->resize(Npla);
  m_vz->resize(Npla);
  m_phi->resize(Npla);
  m_id->resize(Npla);
  m_mask->resize(Npla);

  m_mass->resize(0);
  m_status->resize(0);

  //MPI_Barrier(MPI_COMM_WORLD);

  return;
}



void MC3Extras::refreshParticles(Particles & particles) {
  int Npla = m_xx->size();
  int nReserve = m_xx->capacity();


  //SET UP EXCHANGE
  float oL = Domain::oL();
  float rL = Domain::rL();
  ParticleExchange *rexchange = new ParticleExchange;
  rexchange->setParameters(rL, oL);
  rexchange->initialize();
  m_mass->reserve(nReserve);
  m_status->reserve(nReserve);
  for(int i=0; i<Npla; i++)
    m_mass->push_back(1.0);
  rexchange->setParticles(m_xx, m_yy, m_zz, m_vx, m_vy, m_vz, 
			  m_mass, m_phi, m_id, m_mask, m_status);
  //MPI_Barrier(MPI_COMM_WORLD);


  //EXCHANGE
  rexchange->exchangeParticles();
  //MPI_Barrier(MPI_COMM_WORLD);

#if 0
  printf("REFRESH PXCH rank = %d Npla = %d nReserve = %d size = %d\n",
	 Partition::getMyProc(), Npla, nReserve, xx->size());
  fflush(stdout);
  //MPI_Barrier(MPI_COMM_WORLD);
#endif


  //COPY REFRESHED PARTICLES BACK INTO CLASS
  particles.copyAllFromVectors(&m_xx, &m_yy, &m_zz, &m_vx, &m_vy, &m_vz, 
			       &m_phi, &m_id, &m_mask, m_aa, 1);


  //CLEANUP MEMORY
  deleteVectors();
  delete rexchange;


  //RESORT PARTICLES
  particles.resortParticles();
  //MPI_Barrier(MPI_COMM_WORLD);


  return;
}
