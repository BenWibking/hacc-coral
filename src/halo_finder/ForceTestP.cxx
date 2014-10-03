/*=========================================================================
                                                                                
Copyright (c) 2007, Los Alamos National Security, LLC

All rights reserved.

Copyright 2007. Los Alamos National Security, LLC. 
This software was produced under U.S. Government contract DE-AC52-06NA25396 
for Los Alamos National Laboratory (LANL), which is operated by 
Los Alamos National Security, LLC for the U.S. Department of Energy. 
The U.S. Government has rights to use, reproduce, and distribute this software. 
NEITHER THE GOVERNMENT NOR LOS ALAMOS NATIONAL SECURITY, LLC MAKES ANY WARRANTY,
EXPRESS OR IMPLIED, OR ASSUMES ANY LIABILITY FOR THE USE OF THIS SOFTWARE.  
If software is modified to produce derivative works, such modified software 
should be clearly marked, so as not to confuse it with the version available 
from LANL.
 
Additionally, redistribution and use in source and binary forms, with or 
without modification, are permitted provided that the following conditions 
are met:
-   Redistributions of source code must retain the above copyright notice, 
    this list of conditions and the following disclaimer. 
-   Redistributions in binary form must reproduce the above copyright notice,
    this list of conditions and the following disclaimer in the documentation
    and/or other materials provided with the distribution. 
-   Neither the name of Los Alamos National Security, LLC, Los Alamos National
    Laboratory, LANL, the U.S. Government, nor the names of its contributors
    may be used to endorse or promote products derived from this software 
    without specific prior written permission. 

THIS SOFTWARE IS PROVIDED BY LOS ALAMOS NATIONAL SECURITY, LLC AND CONTRIBUTORS
"AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, 
THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE 
ARE DISCLAIMED. IN NO EVENT SHALL LOS ALAMOS NATIONAL SECURITY, LLC OR 
CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, 
EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, 
PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; 
OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, 
WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR 
OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF 
ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
                                                                                
=========================================================================*/

// .NAME ForceTestP - drive the tree force calculation on particles
//
// .SECTION Description
// Reads the input description file
// Uses ParticleDistribute to read data from files
// Uses ParticleExchange to place overloaded particles on neighbor processors
// Uses ForceTree to build the BHForceTree and do force calculation
//

#include "HaloFinderInput.h"
#include "Partition.h"
#include "Timings.h"

#include "ParticleDistribute.h"
#include "ParticleExchange.h"

#include "ForceTree.h"

#include <iostream>
#include <fstream>
#include <sstream>
#include <rru_mpi.h>

using namespace std;

/////////////////////////////////////////////////////////////////////////////
//
// Class for testing cosmology code
//
/////////////////////////////////////////////////////////////////////////////

class ForceTestP {
public:
  ForceTestP(int argc, char** argv);
  ~ForceTestP();

  // Reads particles from file, distributes to processor
  void DistributeParticles();

  // Do the force calculation on a tree
  void ForceTreeCalculation();
  void ForceTreeCalculationDev();

  int numProc;			// Number of processors
  int myProc;			// Rank of this processor

private:
  HaloFinderInput    haloIn;	// Read input file to direct computation
  ParticleDistribute distribute;// Distributes particles to processors
  ParticleExchange   exchange;	// Exchanges ghost particles

  string haloInputFile;		// Name of driver input file
  string inFile;		// Base name of input particle files
  string outFile;		// Base name of output particle files
  string dataType;		// BLOCK or RECORD structure input
  string distributeType;	// ROUND_ROBIN every proc looks at all files
				// ONE_TO_ONE files match to processors

  float massConvertFactor;	// Multiply every mass read by factor
  float distConvertFactor;	// Multiply every pos read by factor
  float rhocConvertFactor;	// RHO_C based on Mpc/h, convert to match units
  float sodMassConvertFactor;	// SOD_MASS based on Msun/h, convert units
				// If units are Msun and Mpc and user wants to
                                // keep them, rhoc * hubble * hubble
                                // If units are Msun and Mpc and user wants to
                                // convert massFactor = hubble, 
                                // distFactor = hubble, rhocFactor = 1.0

  POSVEL_T rL;			// Physical coordinate box size
  POSVEL_T deadSize;		// Physical coordinate dead or ghost zone
  POSVEL_T bb;			// Distance between particles in a halo
  POSVEL_T omegadm;		// Matter density contribution dark matter
  POSVEL_T deut;		// Matter density contribution baryonic matter
  POSVEL_T omegatot;		// Total matter density
  POSVEL_T hubble;		// Hubble constant

  POSVEL_T RHOC;		// Critical density of universe, Mpc/h, Msun/h
  POSVEL_T SODMASS;		// Factor used to find initial radius for SOD
  POSVEL_T particleMass;	// Physical mass of one particle
  int pmin;			// Minimum number of particles to make a halo
  int np;			// Number of particles in problem cubed

  POSVEL_T alphaFactor;		// Subhalo finding cut/grow factor
  POSVEL_T betaFactor;		// Subhalo finding Poisson noise significance
  int minCandidateSize;		// Smallest particle count for subhalo
  int numSPHNeighbors;		// For calculating smoothing length for subhalo
  int numNeighbors;		// Number of neighbors for subhalo build

  int numberOfParticles;	// Total particles on this processor
  int numberOfFOFHalos;		// Total FOF halos on this processor

  vector<POSVEL_T>* xx;		// Locations of particles on this processor
  vector<POSVEL_T>* yy;
  vector<POSVEL_T>* zz;
  vector<POSVEL_T>* vx;		// Velocities of particles on this processor
  vector<POSVEL_T>* vy;
  vector<POSVEL_T>* vz;
  vector<POSVEL_T>* mass;	// Mass of particles on this processor
  vector<ID_T>* tag;		// Id within entire problem of particle
  vector<STATUS_T>* status;	// ALIVE, DEAD and other information
  vector<POTENTIAL_T>* potential;
  vector<MASK_T>* mask;
};


/////////////////////////////////////////////////////////////////////////////
//
// Constructor
//
/////////////////////////////////////////////////////////////////////////////

ForceTestP::ForceTestP(int argc, char* argv[])
{
  this->numProc = Partition::getNumProc();
  this->myProc = Partition::getMyProc();

  if (argc != 2) {
    cout << "Usage: mpirun -np # HaloTestP halo_finder_input_file" << endl;
  }

  this->haloInputFile = argv[1];
  this->haloIn.initialize(this->haloInputFile);

  // Base file name (Actual file names are basename.proc#
  this->inFile = this->haloIn.getInputBaseName();
  this->outFile = this->haloIn.getOutputBaseName();

  // Mass and distance units
  this->massConvertFactor = this->haloIn.getMassConvertFactor();
  this->distConvertFactor = this->haloIn.getDistConvertFactor();
  this->rhocConvertFactor = this->haloIn.getRHOCConvertFactor();
  this->sodMassConvertFactor = this->haloIn.getSODMassConvertFactor();

  // Physical coordinate box size
  this->rL = (POSVEL_T) this->haloIn.getBoxSize();
  this->rL *= this->distConvertFactor;

  // Physical coordinate dead zone area
  this->deadSize = (POSVEL_T) this->haloIn.getOverloadSize();
  this->deadSize *= this->distConvertFactor;

  // Superimposed grid on physical box used to determine wraparound
  this->np = this->haloIn.getNumberOfParticles();

  // BB parameter for distance between particles comprising a halo
  // This is in grid units for positions normalized on np x np x np grid
  this->bb = (POSVEL_T) this->haloIn.getMinParticleDistance();

  // Minimum number of particles to make a halo
  this->pmin = this->haloIn.getMinParticlesPerHalo();

  // Omegadm
  this->omegadm = (POSVEL_T) this->haloIn.getOmegadm();

  // Hubble constant
  this->hubble = (POSVEL_T) this->haloIn.getHubbleConstant();

  // Deut
  this->deut = (POSVEL_T) this->haloIn.getDeut();

  // Critical density of the universe
  this->RHOC = RHO_C * this->rhocConvertFactor;

  // Factor used to determine initial radius guess for SOD halos from FOF halos
  this->SODMASS = SOD_MASS * this->sodMassConvertFactor;

  // Input is BLOCK or RECORD structured data
  this->dataType = this->haloIn.getInputType();

  // Distribution mechanism is ROUND_ROBIN which selects alive and dead
  // as the data files are read and passed around
  // or EXCHANGE where the alive data is immediately read on a processor
  // and the dead particles must be bundled and shared
  this->distributeType = this->haloIn.getDistributeType();

  // Mass of one particle based on size of problem
  // test4:
  //   rL = 90.1408, hubble = 0.71, omegadm = .27, deut = 0.02218, np = 256
  //        particleMass = 1.917756e09
  // sb256:
  //   rL = 64.0, hubble = 0.5, omegadm = 1.0, deut = 0.0, np = 256
  //        particleMass = 1.08413e09
  // Both of these masses appear in the files and so are correct
  //
  this->omegatot = omegadm + deut/hubble/hubble;
  this->particleMass = this->RHOC * rL * rL * rL * omegatot / np / np / np;
  this->particleMass /= this->massConvertFactor;
  if (this->myProc == 0)
    cout << "Particle mass calculated: " << this->particleMass << endl;

  // Alpha factor controls the CUT/GROW of candidates
  // 1.0 / alphaFactor is the number of times larger a candidate must be
  // in order for the smaller to be CUT rather than allowed to GROW
  // Set to 1.0 and it always cuts, set to 0.01 to aggressively grow
  // small subhalos
  this->alphaFactor = (POSVEL_T) this->haloIn.getAlphaSubhalo();
  
  // Beta factor controls the Poisson noise significance of candidates
  // Original SUBFIND algorithm would have beta = 0.0 meaning that all
  // candidates are allowed to remain as separate and not be COMBINED
  // immediately into the saddle point partner.  If beta is larger it will
  // allow the identification of very small substructures such as tails
  this->betaFactor = (POSVEL_T) this->haloIn.getBetaSubhalo();
  
  // Minimum size of a subhalo candidate
  this->minCandidateSize = this->haloIn.getMinSubhaloSize();

  // BHTree, SPH and density parameters
  this->numSPHNeighbors = this->haloIn.getNumSPHDensity();
  this->numNeighbors = this->haloIn.getNumSubhaloNeighbors();

  if (this->myProc == 0 && this->haloIn.getOutputSubhaloProperties() == 1) {
    cout << "Particle mass: " << this->particleMass << endl;
    cout << "Gravitational constant: " << GRAVITY_C << endl;
    cout << "Potential energy factor: "
         << (this->particleMass * GRAVITY_C) << endl;
    cout << "Cut/Grow factor: " << this->alphaFactor << endl;
    cout << "Poisson noise factor: " << this->betaFactor << endl;
    cout << "Minimum candidate size: " << this->minCandidateSize << endl; 
    cout << "Number of neighbors for SPH: " << this->numSPHNeighbors << endl;
    cout << "Number of neighbors for subgroups: " << this->numNeighbors << endl;
  }

  this->mass = 0;
  this->tag = 0;
  this->status = 0;
  this->potential = 0;
  this->mask = 0;
}

/////////////////////////////////////////////////////////////////////////////
//
// Destructor
//
/////////////////////////////////////////////////////////////////////////////

ForceTestP::~ForceTestP()
{
  if (this->xx != 0) delete this->xx;
  if (this->yy != 0) delete this->yy;
  if (this->zz != 0) delete this->zz;
  if (this->vx != 0) delete this->vx;
  if (this->vy != 0) delete this->vy;
  if (this->vz != 0) delete this->vz;

  if (this->mass != 0) delete this->mass;
  if (this->tag != 0) delete this->tag;
  if (this->status != 0) delete this->status;
  if (this->potential != 0) delete this->potential;
  if (this->mask != 0) delete this->mask;
}

/////////////////////////////////////////////////////////////////////////////
//
// ParticleDistribute reads particles from files and distributes to the
// appropriate processor based on decomposition.
//
// ParticleExchange causes particles on the edges to be shared with
// neighboring processsors making some particles ALIVE and some DEAD
//
/////////////////////////////////////////////////////////////////////////////

void ForceTestP::DistributeParticles()
{
  static Timings::TimerRef dtimer = Timings::getTimer("Distribute Particles");
  Timings::startTimer(dtimer);

  // Initialize classes for reading, exchanging and calculating
  this->distribute.setParameters(this->inFile, this->rL, this->dataType);
  this->distribute.setConvertParameters(this->massConvertFactor,
                                        this->distConvertFactor);
  this->exchange.setParameters(this->rL, this->deadSize);

  this->distribute.initialize();
  this->exchange.initialize();

  // Read alive particles only from files
  // In ROUND_ROBIN all files are read and particles are passed round robin
  // to every other processor so that every processor chooses its own
  // In ONE_TO_ONE every processor reads its own processor in the topology
  // which has already been populated with the correct alive particles

  this->xx = new vector<POSVEL_T>;
  this->yy = new vector<POSVEL_T>;
  this->zz = new vector<POSVEL_T>;
  this->vx = new vector<POSVEL_T>;
  this->vy = new vector<POSVEL_T>;
  this->vz = new vector<POSVEL_T>;
  this->mass = new vector<POSVEL_T>;
  this->tag = new vector<ID_T>;
  this->status = new vector<STATUS_T>;

  this->distribute.setParticles(this->xx, this->yy, this->zz, 
                                this->vx, this->vy, this->vz,
                                this->mass, this->tag);
  if (this->distributeType == "ROUND_ROBIN")
    this->distribute.readParticlesRoundRobin();
  else if (this->distributeType == "ONE_TO_ONE")
    this->distribute.readParticlesOneToOne();

  // Create the mask and potential vectors which will be filled in elsewhere
  this->numberOfParticles = this->xx->size();
  this->potential = new vector<POTENTIAL_T>(this->numberOfParticles);
  this->mask = new vector<MASK_T>(this->numberOfParticles);

  // Exchange particles adds dead particles to all the vectors
  this->exchange.setParticles(this->xx, this->yy, this->zz, 
                              this->vx, this->vy, this->vz, this->mass,
                              this->potential, this->tag, 
                              this->mask, this->status);
  this->exchange.exchangeParticles();

  // If mass is a constant 1.0 must reset to particleMass
  for (int i = 0; i < this->numberOfParticles; i++)
    if ((*this->mass)[i] == 1.0)
      (*this->mass)[i] = this->particleMass;

  Timings::stopTimer(dtimer);
}

/////////////////////////////////////////////////////////////////////////////
//
// Build the BHForceTree and run a force calculation on it
//
/////////////////////////////////////////////////////////////////////////////

void ForceTestP::ForceTreeCalculation()
{
  static Timings::TimerRef shtimer = Timings::getTimer("Tree Force");
  Timings::startTimer(shtimer);

  // Normalize the particle locations to be placed on a grid
  POSVEL_T normalizeFactor = (POSVEL_T) ((1.0 * this->np) / this->rL);
  for (int p = 0; p < this->numberOfParticles; p++) {
    (*this->xx)[p] *= normalizeFactor;
    (*this->yy)[p] *= normalizeFactor;
    (*this->zz)[p] *= normalizeFactor;
  }

  // Create the force tree
  ForceTree* forceTree = new ForceTree();

  // Overload size in normalized grid units
  POSVEL_T overLoadSize = this->deadSize * normalizeFactor;

  // Set the range of the force tree for this processor
  // Must include the overload zone
  POSVEL_T minPosition[DIMENSION];
  POSVEL_T maxPosition[DIMENSION];
  for (int dim = 0; dim < DIMENSION; dim++) {
    minPosition[dim] = -overLoadSize;
    maxPosition[dim] = (POSVEL_T) this->np + overLoadSize;
  }

  // Criteria for opening a node instead of accepting as a collection
  POSVEL_T openingAngle = 0.5;

  // Distance from a particle beyond which particles and nodes are ignored
  POSVEL_T criticalRadius = 3.12;

  // Subgroup method sometimes produces small numbers of particles
  // Less than minimum, run the top down method instead
  // Greater than maximum and continue to open to find smaller subgroups
  int minimumGroup = 16;
  int maximumGroup = 256;

  forceTree->setParameters(minPosition, maxPosition, 
                           openingAngle,
                           criticalRadius,
                           minimumGroup,
                           maximumGroup,
                           this->particleMass);

  forceTree->setParticles(this->xx, this->yy, this->zz, 
                          this->vx, this->vy, this->vz,
                          this->mass, this->potential);

  // Construct the BH force tree from normalized particles
  // Threading also calculates particle radius for each node
  static Timings::TimerRef bhtimer = Timings::getTimer("Build/Thread BHTree");
  Timings::startTimer(bhtimer);
  forceTree->buildForceTree();
  Timings::stopTimer(bhtimer);

/*
  static Timings::TimerRef n2timer = Timings::getTimer("N^2 Force");
  Timings::startTimer(n2timer);
  if (this->myProc == 0)
    cout << "GADGET N^2 METHOD" << endl;
  if (this->myProc == 0)
  forceTree->forceCalculationN2();
  if (this->myProc == 0)
  forceTree->printForceValues();
  Timings::stopTimer(n2timer);
  MPI_Barrier(MPI_COMM_WORLD);
*/

/*
  static Timings::TimerRef gbtimer = Timings::getTimer("Gadget Bottom Up Walk");
  Timings::startTimer(gbtimer);
  if (this->myProc == 0)
    cout << "GADGET BOTTOM UP METHOD" << endl;
//  if (this->myProc == 0)
  forceTree->forceCalculationGadgetBottomUp();
//  if (this->myProc == 0)
//  forceTree->printForceValues();
  Timings::stopTimer(gbtimer);
  MPI_Barrier(MPI_COMM_WORLD);
*/

/*
  static Timings::TimerRef gtimer = Timings::getTimer("Gadget Top Down Walk");
  Timings::startTimer(gtimer);
  if (this->myProc == 0)
    cout << "GADGET TOP DOWN METHOD" << endl;
//  if (this->myProc == 0)
  forceTree->forceCalculationGadgetTopDown();
//  if (this->myProc == 0)
//  forceTree->printForceValues();
  Timings::stopTimer(gtimer);
  MPI_Barrier(MPI_COMM_WORLD);
*/

  static Timings::TimerRef grtimer = Timings::getTimer("Group Walk");
  Timings::startTimer(grtimer);
  if (this->myProc == 0)
    cout << "GROUP METHOD " << maximumGroup << endl;
//  if (this->myProc == 0)
  forceTree->forceCalculationGroup();
//  if (this->myProc == 0)
//  forceTree->printForceValues();
  Timings::stopTimer(grtimer);
  MPI_Barrier(MPI_COMM_WORLD);

/*
  static Timings::TimerRef batimer = Timings::getTimer("Barnes Adjusted Walk");
  Timings::startTimer(batimer);
  if (this->myProc == 0)
    cout << "BARNES ADJUST METHOD" << endl;
//  if (this->myProc == 0)
  forceTree->forceCalculationBarnesAdjust();
//  if (this->myProc == 0)
//  forceTree->printForceValues();
  Timings::stopTimer(batimer);
  MPI_Barrier(MPI_COMM_WORLD);
*/

/*
  static Timings::TimerRef bqtimer = Timings::getTimer("Barnes Quick Walk");
  Timings::startTimer(bqtimer);
  if (this->myProc == 0)
    cout << "BARNES QUICK METHOD" << endl;
//  if (this->myProc == 0)
  forceTree->forceCalculationBarnesQuick();
//  if (this->myProc == 0)
//  forceTree->printForceValues();
  Timings::stopTimer(bqtimer);
  MPI_Barrier(MPI_COMM_WORLD);
*/

  delete forceTree;
  Timings::stopTimer(shtimer);
}

/////////////////////////////////////////////////////////////////////////////
//
// Test driver
//
/////////////////////////////////////////////////////////////////////////////

int main(int argc, char* argv[])
{
  // Initialize MPI and set the decomposition in the Partition class
  MPI_Init(&argc, &argv);
  Partition::initialize();

  // Construct the tester class
  ForceTestP forceTestP(argc, argv);

  // Read, distribute and exchange particles
  forceTestP.DistributeParticles();
  MPI_Barrier(MPI_COMM_WORLD);

  // Find force calculation
  forceTestP.ForceTreeCalculation();
  cout << "Rank " << Partition::getMyProc() << " FINISHED " << endl;

  MPI_Barrier(MPI_COMM_WORLD);
  Timings::print();

  // Shut down MPI
  Partition::finalize();
  MPI_Finalize();

  return 0;
}
