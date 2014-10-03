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

// .NAME SubHaloTest - drive the subhalo finding from input to output
//
// .SECTION Description
//

#include "HaloFinderInput.h"
#include "Timings.h"

#include "CosmoHaloFinder.h"

#include "FOFHaloProperties.h"
#include "HaloCenterFinder.h"
#include "SubHaloFinderDev.h"

#include "SODHalo.h"
#include "ChainingMesh.h"

#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <rru_mpi.h>

using namespace std;

/////////////////////////////////////////////////////////////////////////////
//
// Class for testing cosmology code
//
/////////////////////////////////////////////////////////////////////////////

class SubHaloTest {
public:
  SubHaloTest(int argc, char** argv);
  ~SubHaloTest();

  // Reads particles from file
  void ReadFromRecordFile();

  // Finds FOF halos uniquely on each processor
  void FOFSubHalo();

  // Subhalo finding requires extraction of halo particles from FOF halo list
  // and building of a Barnes Hut tree of those particles
  void FindSubHalos();

  // Write a .cosmo file of halos of this size or greater
  void WriteSubhaloCosmoFile(string name);

private:
  HaloFinderInput    haloIn;	// Read input file to direct computation
  CosmoHaloFinder    haloFinder;// Serial FOF halo finder

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
  POSVEL_T alphaFactor;		// Factor for cut/grow
  POSVEL_T betaFactor;		// Factor for Poisson noise significance
  int minCandidateSize;         // Smallest particle count for subhalo

  int pmin;			// Minimum number of particles to make a halo
  int np;			// Number of particles in problem cubed

  ID_T numberOfParticles;	// Total particles on this processor
  int numberOfFOFHalos;		// Total FOF halos on this processor

  POSVEL_T* xx;			// Locations of particles on this processor
  POSVEL_T* yy;
  POSVEL_T* zz;
  POSVEL_T* vx;			// Velocities of particles on this processor
  POSVEL_T* vy;
  POSVEL_T* vz;
  POSVEL_T* mass;		// Mass of particles on this processor
// PKF
  POSVEL_T* smoothingLength;
  POSVEL_T* density;	

  ID_T* tag;			// Id within entire problem of particle

  POSVEL_T normalizeFactor;     // Convert physical location to grid location
  POSVEL_T** haloData;          // Normalized data for serial halo finder
  int* haloTag;                 // From serial halo finder, the index of the
                                // first particle in a halo
  int* haloSize;                // From serial halo finder, the size of a halo
                                // where the first particle has the actual size
                                // and other member particles have size=0

  vector<int> halos;            // First particle index into haloList
  vector<int> haloCount;        // Size of each halo 

  int* haloList;                // Indices of next particle in halo
  int* haloStart;               // Index of first particle in halo
                                // Chain is built backwards but using these two
                                // arrays, all particle indices for a halo
                                // can be found
  int numberOfHaloParticles;    // Number of particles in all halos
  int numberOfHalos;            // Number of halos

  int numSPHNeighbors;		// Number of neighbor for local density
  int numNeighbors;		// Number of neighbors for algorithm
};


/////////////////////////////////////////////////////////////////////////////
//
// Constructor
//
/////////////////////////////////////////////////////////////////////////////

SubHaloTest::SubHaloTest(int argc, char* argv[])
{
  if (argc != 2) {
    cout << "Usage: mpirun -np # SubHaloTest halo_finder_input_file" << endl;
  }

  this->haloInputFile = argv[1];
  this->haloIn.initialize(this->haloInputFile);

  // Base file name (Actual file names are basename.proc#
  this->inFile = haloIn.getInputBaseName();
  this->outFile = haloIn.getOutputBaseName();

  // Mass and distance units
  this->massConvertFactor = haloIn.getMassConvertFactor();
  this->distConvertFactor = haloIn.getDistConvertFactor();
  this->rhocConvertFactor = haloIn.getRHOCConvertFactor();
  this->sodMassConvertFactor = haloIn.getSODMassConvertFactor();

  // Physical coordinate box size
  this->rL = (POSVEL_T) haloIn.getBoxSize();
  this->rL *= this->distConvertFactor;

  // Physical coordinate dead zone area
  this->deadSize = (POSVEL_T) haloIn.getOverloadSize();
  this->deadSize *= this->distConvertFactor;

  // Superimposed grid on physical box used to determine wraparound
  this->np = haloIn.getNumberOfParticles();

  // BB parameter for distance between particles comprising a halo
  // This is in grid units for positions normalized on np x np x np grid
  this->bb = (POSVEL_T) haloIn.getMinParticleDistance();

  // Minimum number of particles to make a halo
  this->pmin = haloIn.getMinParticlesPerHalo();

  // Omegadm
  this->omegadm = (POSVEL_T) (POSVEL_T) haloIn.getOmegadm();

  // Hubble constant
  this->hubble = (POSVEL_T) haloIn.getHubbleConstant();

  // Deut
  this->deut = (POSVEL_T) haloIn.getDeut();

  // Critical density of the universe
  this->RHOC = RHO_C * this->rhocConvertFactor;

  // Factor used to determine initial radius guess for SOD halos from FOF halos
  this->SODMASS = SOD_MASS * this->sodMassConvertFactor;

  // Input is BLOCK or RECORD structured data
  this->dataType = haloIn.getInputType();

  // Distribution mechanism is ROUND_ROBIN which selects alive and dead
  // as the data files are read and passed around
  // or EXCHANGE where the alive data is immediately read on a processor
  // and the dead particles must be bundled and shared
  this->distributeType = haloIn.getDistributeType();

  // For now enter the particle mass, later I can get from input file
  this->particleMass = 0.0860657e10;

  // Alpha factor controls the CUT/GROW of candidates
  // 1.0 / alphaFactor is the number of times larger a candidate must be
  // in order for the smaller to be CUT rather than allowed to GROW
  // Set to 1.0 means always cut
  this->alphaFactor = ALPHA_SUBHALO;

  // Beta factor controls the Poisson noise significance of candidates
  // Original SUBFIND algorithm would have beta = 0.0 meaning that all
  // candidates are allowed to remain as separate and not be COMBINED
  // immediately into the saddle point partner.  If beta is larger it will
  // allow the identification of very small substructures such as tails
  this->betaFactor = BETA_SUBHALO;

  // Minimum size of a subhalo candidate
  this->minCandidateSize = MIN_SUBHALO_SIZE;

  // BHTree, SPH and density parameters
  this->numSPHNeighbors = NUM_SPH_DENSITY;
  this->numNeighbors = NUM_SUBHALO_NEIGHBOR;

  cout << "Particle mass: " << this->particleMass << endl;
  cout << "Gravitational constant: " << GRAVITY_C << endl;
  cout << "Potential energy factor: " 
       << (this->particleMass * GRAVITY_C) << endl;
  cout << "Cut/Grow factor: " << this->alphaFactor << endl;
  cout << "Poisson noise factor: " << this->betaFactor << endl;
  cout << "Minimum candidate size: " << this->minCandidateSize << endl;
  cout << "Number of neighbors for SPH: " << this->numSPHNeighbors << endl;
  cout << "Number of neighbors for subgroups: " << this->numNeighbors << endl;

  this->xx = 0;
  this->yy = 0;
  this->zz = 0;
  this->vx = 0;
  this->vy = 0;
  this->vz = 0;
  this->mass = 0;
  this->tag = 0;
}

/////////////////////////////////////////////////////////////////////////////
//
// Destructor
//
/////////////////////////////////////////////////////////////////////////////

SubHaloTest::~SubHaloTest()
{
  if (this->xx != 0) delete [] this->xx;
  if (this->yy != 0) delete [] this->yy;
  if (this->zz != 0) delete [] this->zz;
  if (this->vx != 0) delete [] this->vx;
  if (this->vy != 0) delete [] this->vy;
  if (this->vz != 0) delete [] this->vz;

  if (this->mass != 0) delete [] this->mass;
  if (this->tag != 0) delete [] this->tag;
}

/////////////////////////////////////////////////////////////////////////////
//
// Read particles for one halo in directly to one processor
//
/////////////////////////////////////////////////////////////////////////////

void SubHaloTest::ReadFromRecordFile()
{
  // Open halo particle file
  ifstream inStream(this->inFile.c_str(), ios::in);

  // Compute the number of particles from file size
  inStream.seekg(0L, ios::end);
// PKF
int pkfRecordSize = sizeof(POSVEL_T) * (COSMO_FLOAT+2) + 
                          sizeof(ID_T) * COSMO_INT;
  //this->numberOfParticles = inStream.tellg() / RECORD_SIZE;
  this->numberOfParticles = inStream.tellg() / pkfRecordSize;
  inStream.seekg(0L, ios::beg);

  cout << "Open file " << this->inFile
       << " with " << this->numberOfParticles << " particles" << endl;

  this->xx = new POSVEL_T[this->numberOfParticles];
  this->yy = new POSVEL_T[this->numberOfParticles];
  this->zz = new POSVEL_T[this->numberOfParticles];
  this->vx = new POSVEL_T[this->numberOfParticles];
  this->vy = new POSVEL_T[this->numberOfParticles];
  this->vz = new POSVEL_T[this->numberOfParticles];
  this->mass = new POSVEL_T[this->numberOfParticles];
// PKF
  this->smoothingLength = new POSVEL_T[this->numberOfParticles];
  this->density = new POSVEL_T[this->numberOfParticles];

  this->tag = new ID_T[this->numberOfParticles];

//PKF
  //POSVEL_T* fBlock = new POSVEL_T[COSMO_FLOAT];
  POSVEL_T* fBlock = new POSVEL_T[COSMO_FLOAT + 2];
  ID_T* iBlock = new ID_T[COSMO_INT];

  // Store each particle location, velocity and tag
  for (int p = 0; p < this->numberOfParticles; p++) {

    // Set file pointer to the requested particle
    inStream.read(reinterpret_cast<char*>(fBlock),
// PKF
                   //COSMO_FLOAT * sizeof(POSVEL_T));
                   (COSMO_FLOAT+2) * sizeof(POSVEL_T));

// PKF
    //if (inStream.gcount() != COSMO_FLOAT * sizeof(POSVEL_T)) {
    if (inStream.gcount() != (COSMO_FLOAT+2) * sizeof(POSVEL_T)) {
      cout << "Premature end-of-file" << endl;
      exit (-1);
    }

    // Convert units if requested
    fBlock[0] *= this->distConvertFactor;
    fBlock[2] *= this->distConvertFactor;
    fBlock[4] *= this->distConvertFactor;
    fBlock[6] *= this->massConvertFactor;

    inStream.read(reinterpret_cast<char*>(iBlock),
                   COSMO_INT * sizeof(ID_T));

    if (inStream.gcount() != COSMO_INT * sizeof(ID_T)) {
      cout << "Premature end-of-file" << endl;
      exit (-1);
    }

    this->xx[p] = fBlock[0];
    this->vx[p] = fBlock[1];
    this->yy[p] = fBlock[2];
    this->vy[p] = fBlock[3];
    this->zz[p] = fBlock[4];
    this->vz[p] = fBlock[5];
    this->mass[p] = fBlock[6];
// PKF
    this->smoothingLength[p] = fBlock[7];
    this->density[p] = fBlock[8];

    this->tag[p] = iBlock[0];
//cout << "Read particle " << p << " loc " << xx[p] << "," << yy[p] << "," << zz[p] << " vel " << vx[p] << "," << vy[p] << "," << vz[p] << " mass " << mass[p] << " smooth " << smoothingLength[p] << " density " << density[p] << " tag " << tag[p] << endl;
//cout << "Read particle " << p << " loc " << xx[p] << "," << yy[p] << "," << zz[p] << " vel " << vx[p] << "," << vy[p] << "," << vz[p] << " mass " << mass[p] << " tag " << tag[p] << endl;
  }

  inStream.close();
  delete [] fBlock;
  delete [] iBlock;
}

/////////////////////////////////////////////////////////////////////////////
//
// Find the FOF Halos
//
/////////////////////////////////////////////////////////////////////////////

void SubHaloTest::FOFSubHalo()
{
  static Timings::TimerRef h1timer = Timings::getTimer("FOF Halo Finder");
  Timings::startTimer(h1timer);

  // Serial halo finder
  this->haloFinder.np = this->np;
  this->haloFinder.pmin = this->pmin;
  this->haloFinder.bb = this->bb;
  this->haloFinder.rL = rL;
  this->haloFinder.periodic = false;
  this->haloFinder.textmode = "ascii";

  this->normalizeFactor = (POSVEL_T)((1.0 * this->np) / this->rL);

  // Execute the serial halo finder
  this->haloData = new POSVEL_T*[DIMENSION];
  for (int dim = 0; dim < DIMENSION; dim++)
    this->haloData[dim] = new POSVEL_T[this->numberOfParticles];

  // Fill it with normalized x,y,z of all particles on this processor
  for (int p = 0; p < this->numberOfParticles; p++) {
    this->haloData[0][p] = this->xx[p] * this->normalizeFactor;
    this->haloData[1][p] = this->yy[p] * this->normalizeFactor;
    this->haloData[2][p] = this->zz[p] * this->normalizeFactor;
  }

  this->haloFinder.setParticleLocations(haloData);
  this->haloFinder.setNumberOfParticles(this->numberOfParticles);
  this->haloFinder.setOutFile(this->outFile);

  cout << "RUNNING SERIAL HALO FINDER on "
       << this->numberOfParticles << " particles" << endl;

  if (this->numberOfParticles > 0)
    this->haloFinder.Finding();

  // Collect the serial halo finder results into a usable structure
  // Halo tag returned from the serial halo finder is actually the index
  // of the particle on this processor.  Must map to get to actual tag
  // which is common information between all processors.
  this->haloTag = haloFinder.getHaloTag();

  // Record the halo size of each particle on this processor
  // Create a list of particles in any halo by recording the index of the
  // first particle and having that index give the index of the next particle
  // Last particle index reports a -1
  // List is built by iterating on the tags and storing in opposite order so
  this->haloSize = new int[this->numberOfParticles];
  this->haloList = new int[this->numberOfParticles];
  this->haloStart = new int[this->numberOfParticles];

  for (int p = 0; p < this->numberOfParticles; p++) {
    this->haloList[p] = -1;
    this->haloStart[p] = p;
    this->haloSize[p] = 0;
  }

  // Link all particles in the same halo together
  for (int p = 0; p < this->numberOfParticles; p++) {
    if (this->haloTag[p] != p) {
      this->haloList[p] = haloStart[this->haloTag[p]];
      this->haloStart[this->haloTag[p]] = p;
    }
    this->haloSize[this->haloTag[p]]++;
  }

  // Only the first particle id for a halo records the size
  // Succeeding particles which are members of a halo have a size of 0
  // Record the start index of any legal halo which will allow the
  // following of the chaining mesh to identify all particles in a halo
  this->numberOfHaloParticles = 0;
  this->numberOfHalos = 0;

  for (ID_T p = 0; p < this->numberOfParticles; p++) {

    if (this->haloSize[p] >= this->pmin) {
cout << "Particle " << p << " haloSize " << this->haloSize[p] << endl;
      this->numberOfHalos++;
      this->numberOfHaloParticles += this->haloSize[p];
      this->halos.push_back(this->haloStart[p]);
      this->haloCount.push_back(this->haloSize[p]);

/*
      // Color all particle masses for this halo by the halo index for viz
      int q = this->haloStart[p];
      while (q != -1) {
        this->mass[q] = (POSVEL_T) this->numberOfHalos;
        q = this->haloList[q];
      }
*/
    }
  }

  Timings::stopTimer(h1timer);
}

/////////////////////////////////////////////////////////////////////////////
//
// Subhalo finder requires special data structures so extract the locations 
// for each individual FOF halo into arrays for passing to SubHaloFinder
// An advantage to passing arrays to these classes is that they can
// be used with any set of particles a user wants processed
// because they don't iterate over halos and haloList
//
/////////////////////////////////////////////////////////////////////////////

void SubHaloTest::FindSubHalos()
{
  static Timings::TimerRef shtimer = Timings::getTimer("Find SubHalos");
  Timings::startTimer(shtimer);

  // Look for subhalos within the FOF halo using extract location arrays
  SubHaloFinder subFinder;
  subFinder.setParameters(this->particleMass, GRAVITY_C, 
                          this->alphaFactor, this->betaFactor,
                          this->minCandidateSize,
                          this->numSPHNeighbors, this->numNeighbors);
  subFinder.setParticles(this->numberOfParticles, 
                        &this->xx[0], &this->yy[0], &this->zz[0], 
                        &this->vx[0], &this->vy[0], &this->vz[0], 
// PKF
                        //&this->mass[0], &this->tag[0]);
                        &this->mass[0], &this->smoothingLength[0], &this->density[0], &this->tag[0]);
  subFinder.findSubHalos();

  Timings::stopTimer(shtimer);
}

/////////////////////////////////////////////////////////////////////////////
//
// Write a .cosmo file of the requested halo
//
/////////////////////////////////////////////////////////////////////////////
  
void SubHaloTest::WriteSubhaloCosmoFile(string str)
{       
  ostringstream name;
  name << str << this->numberOfParticles << ".cosmo";
  ofstream cStream(name.str().c_str(), ios::out|ios::binary);

  float fBlock[COSMO_FLOAT];
  int iBlock[COSMO_INT];

  for (int p = 0; p < this->numberOfParticles; p++) {
    fBlock[0] = this->xx[p];
    fBlock[1] = this->vx[p];
    fBlock[2] = this->yy[p];
    fBlock[3] = this->vy[p];
    fBlock[4] = this->zz[p];
    fBlock[5] = this->vz[p];
    fBlock[6] = this->mass[p];
    cStream.write(reinterpret_cast<char*>(fBlock), 
                  COSMO_FLOAT * sizeof(POSVEL_T));
    iBlock[0] = this->tag[p];
    cStream.write(reinterpret_cast<char*>(iBlock), 
                  COSMO_INT * sizeof(ID_T));
  }
  cStream.close();
}   

/////////////////////////////////////////////////////////////////////////////
//
// Test driver
//
/////////////////////////////////////////////////////////////////////////////

int main(int argc, char* argv[])
{
  // Construct the tester class
  SubHaloTest subhaloTest(argc, argv);

  // Read particles
  subhaloTest.ReadFromRecordFile();

  // Recursive FOF method
  // Don't run this if we start with particles as in Millennium data
  //subhaloTest.FOFSubHalo();

  // Write the halo catalog and a summary of FOF properties
  string fofStr = "FOFSubHalo_";
  //subhaloTest.WriteSubhaloCosmoFile(fofStr);

  // Find subhalos method
  subhaloTest.FindSubHalos();

  // Write the halo catalog and a summary of subhalo properties
  string subfindStr = "SUBHalo_";
  subhaloTest.WriteSubhaloCosmoFile(subfindStr);

  return 0;
}
