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

#include <stdlib.h>
#include <fstream>
#include <iostream>
#include <string.h>

#include "BasicDefinition.h"

using namespace std;

#define DEBUG    1

/////////////////////////////////////////////////////////////////////////
//
// Read in the requested number of characters
//
/////////////////////////////////////////////////////////////////////////

string readString(ifstream* inStr, int size)
{
   char* buffer = new char[size + 1];
   inStr->read(buffer, size);
   buffer[size] = '\0';

   // Make sure string has legal values
   if (isalnum(buffer[0]) == 0)
      buffer[0] = '\0';
   for (int i = 1; i < size; i++)
      if (isprint(buffer[i]) == 0)
         buffer[i] = '\0';

   string retString = buffer;
   delete [] buffer;
   return retString;
}

/////////////////////////////////////////////////////////////////////////
//
// Read in the number of items from the file pointer and
// byte swap if necessary
//
/////////////////////////////////////////////////////////////////////////

void readData(
        bool swap,
        void* data,
        unsigned long dataSize,
        unsigned long dataCount,
        ifstream* inStr)
{
   // Read all the data from the file
   inStr->read(reinterpret_cast<char*>(data), dataSize*dataCount);

   if (swap == true) {

      // Byte swap each integer
      char* dataPtr = (char*) data;
      char temp;
      for (unsigned long item = 0; item < dataCount; item++) {

         // Do a byte-by-byte swap, reversing the order.
         for (unsigned int i = 0; i < dataSize / 2; i++) {
            temp = dataPtr[i];
            dataPtr[i] = dataPtr[dataSize - 1 - i];
            dataPtr[dataSize - 1 - i] = temp;
         }
         dataPtr += dataSize;
      }
   }
}

/////////////////////////////////////////////////////////////////////////////
//
// First command line parameter is the Gadget file name
// Second command line parameter is the Cosmo file name
// Third command line parameter is conversion factor on mass
// Fourth command line parameter is conversion factor on position
//
// NOTE: Compile this program using the definitions in include.mk for
//       data type sizes.
//
//       HF_TYPE_FLAGS := -DID_32 -DPOSVEL_32 -DGRID_32 (32 bit data and tag)
//       HF_TYPE_FLAGS := -DID_64 -DPOSVEL_32 -DGRID_32 (32 bit data
//                                                       64 bit tag)
//       The latter is needed for data sets larger than 512^3 particles
//
/////////////////////////////////////////////////////////////////////////////
//
// Gadget format (BLOCK):
//    SKIP_GADGET_2 has extra 16 bytes
//    SKIP_H 4 bytes (size of header)
//    Header (256 bytes)
//    SKIP_H 4 bytes (size of header)
//
//    SKIP_GADGET_2 has extra 16 bytes
//    SKIP_L 4 bytes (size of location block in bytes)
//    Block of location data where each particle's x,y,z is stored together
//    SKIP_L 4 bytes (size of location block in bytes)
//
//    SKIP_GADGET_2 has extra 16 bytes
//    SKIP_V 4 bytes (size of velocity block in bytes)
//    Block of velocity data where each particle's xv,yv,zv is stored together
//    SKIP_V 4 bytes (size of velocity block in bytes)
//
//    SKIP_GADGET_2 has extra 16 bytes
//    SKIP_T 4 bytes (size of tag block in bytes)
//    Block of tag data
//    SKIP_T 4 bytes (size of tag block in bytes)
//
//    Header file npart[6] array indicates the number of particles of each
//    type stored in the file.  The types are:
//
//       0 Gas
//       1 Halo
//       2 Disk
//       3 Bulge
//       4 Stars
//       5 Boundary
//
//    If a mass is given as 0, but there are particles, then there is a
//    variable mass block
//
//    Gas particles have extra information blocks after the tags
//    Internal energy, density and smoothing length
//
/////////////////////////////////////////////////////////////////////////////
//
// Cosmo format (RECORD):
//    X location
//    X velocity
//    Y location
//    Y velocity
//    Z location
//    Z velocity
//    Mass
//    Tag
//
/////////////////////////////////////////////////////////////////////////////

int main(int argc, char** argv)
{
  if (argc != 5) {
    cout << "Usage: Gadget2ToCosmo " << endl;
    cout << "         gadget-file (input file) " << endl;
    cout << "         cosmo-file (output file " << endl;
    cout << "         mass-factor (to make M_sun/h) " << endl;
    cout << "         loc-factor (to make Mpc/h) " << endl;
    exit(-1);
  }

  int cnt = 1;
  string inFile = argv[cnt++];
  string outFile = argv[cnt++];
  float massFactor = atof(argv[cnt++]);
  float posFactor = atof(argv[cnt++]);

  int blockSize, blockSize2;
  string gadget2;

  // Initialize for Gadget-1 and no swap
  string gadgetType[NUM_GADGET_TYPES] = {
    "Gas", "Halo", "Disk", "Bulge", "Stars", "Boundary"};
  struct GadgetHeader header;
  int format = GADGET_1;
  bool swap = false;

  // Open input Gadget-2 file
  ifstream *gStr = new ifstream(inFile.c_str(), ios::in | ios::binary);
  if (gStr->fail()) {
    cout << "File: " << inFile << " cannot be opened" << endl;
    exit(-1);
  }

  // Set the gadget format type by reading the first 4 byte integer
  // If it is not "256" or "65536" then gadget-2 format with 16 bytes in front
  readData(swap, (void*) &blockSize, GADGET_SKIP, 1, gStr);
  if (blockSize != GADGET_HEADER_SIZE && blockSize != GADGET_HEADER_SIZE_SWP) {
    format = GADGET_2;
    gadget2 = readString(gStr, GADGET_2_SKIP - GADGET_SKIP);
    readData(swap, (void*) &blockSize, GADGET_SKIP, 1, gStr);
  }

  // Set the swap type
  if (blockSize != GADGET_HEADER_SIZE) {
    swap = true;
    blockSize = GADGET_HEADER_SIZE;
  }

  // Read the Gadget header
  readData(swap, (void*) &header.npart[0], 
                         sizeof(int), NUM_GADGET_TYPES, gStr);
  readData(swap, (void*) &header.mass[0], 
                         sizeof(double), NUM_GADGET_TYPES, gStr);
  readData(swap, (void*) &header.time, sizeof(double), 1, gStr);
  readData(swap, (void*) &header.redshift, sizeof(double), 1, gStr);
  readData(swap, (void*) &header.flag_sfr, sizeof(int), 1, gStr);
  readData(swap, (void*) &header.flag_feedback, sizeof(int), 1, gStr);
  readData(swap, (void*) &header.npartTotal[0], 
                         sizeof(int), NUM_GADGET_TYPES, gStr);
  readData(swap, (void*) &header.flag_cooling, sizeof(int), 1, gStr);
  readData(swap, (void*) &header.num_files, sizeof(int), 1, gStr);
  readData(swap, (void*) &header.BoxSize, sizeof(double), 1, gStr);
  readData(swap, (void*) &header.Omega0, sizeof(double), 1, gStr);
  readData(swap, (void*) &header.OmegaLambda, sizeof(double), 1, gStr);
  readData(swap, (void*) &header.HubbleParam, sizeof(double), 1, gStr);
  readData(swap, (void*) &header.flag_stellarage, sizeof(int), 1, gStr);
  readData(swap, (void*) &header.flag_metals, sizeof(int), 1, gStr);
  readData(swap, (void*) &header.HighWord[0], 
                         sizeof(int), NUM_GADGET_TYPES, gStr);
  readData(swap, (void*) &header.flag_entropy, sizeof(int), 1, gStr);
  string fill = readString(gStr, 60);
  strcpy(&header.fill[0], fill.c_str());

  // Read the Gadget header size to verify block
  readData(swap, (void*) &blockSize2, GADGET_SKIP, 1, gStr);
  if (blockSize != blockSize2)
    cout << "Error reading header: end position is wrong" << endl;

  // Every type particle will have location, velocity and tag so sum up
  long int particleCount = 0;
  for (int i = 0; i < NUM_GADGET_TYPES; i++)
    particleCount += header.npart[i];

#ifdef DEBUG
  for (int i = 0; i < NUM_GADGET_TYPES; i++) {
    if (header.npart[i] > 0) {
      cout << endl;
      cout << "Type of particle: " << gadgetType[i] << endl;
      cout << "Number of particles: " << header.npart[i] << endl;
      cout << "Particle mass: " << header.mass[i] << endl << endl;
    }
  }
  cout << "Total particles this file: " << particleCount << endl;
  cout << "Number of files: " << header.num_files << endl;
  cout << "BoxSize: " << header.BoxSize << endl;
  cout << "Omega0: " << header.Omega0 << endl;
  cout << "OmegaLambda: " << header.OmegaLambda << endl;
  cout << "HubbleParam: " << header.HubbleParam << endl;
#endif
  
  // Allocate storage for all blocks because Cosmo format requires records
  POSVEL_T* location = new POSVEL_T[particleCount * DIMENSION];
  POSVEL_T* velocity = new POSVEL_T[particleCount * DIMENSION];
  ID_T* tag = new ID_T[particleCount];

  // Only gas particles which are SPH have extra blocks of data
  POSVEL_T* u = new POSVEL_T[header.npart[GADGET_GAS]];
  POSVEL_T* rho = new POSVEL_T[header.npart[GADGET_GAS]];
  POSVEL_T* hsml = new POSVEL_T[header.npart[GADGET_GAS]];

  // Skip for Gadget-2
  if (format == GADGET_2)
    gadget2 = readString(gStr, GADGET_2_SKIP);

  //
  // Read locations
  //
  readData(swap, (void*) &blockSize, GADGET_SKIP, 1, gStr);
  readData(swap, (void*) &location[0], 
                         sizeof(POSVEL_T), particleCount*DIMENSION, gStr);
  readData(swap, (void*) &blockSize2, GADGET_SKIP, 1, gStr);
  if (blockSize != blockSize2)
    cout << "Error reading locations: end position is wrong "
         << blockSize << " vs " << blockSize2 << endl;

  // Skip for Gadget-2
  if (format == GADGET_2)
    gadget2 = readString(gStr, GADGET_2_SKIP);

  //
  // Read velocities
  //
  readData(swap, (void*) &blockSize, GADGET_SKIP, 1, gStr);
  readData(swap, (void*) &velocity[0], 
                         sizeof(POSVEL_T), particleCount*DIMENSION, gStr);
  readData(swap, (void*) &blockSize2, GADGET_SKIP, 1, gStr);
  if (blockSize != blockSize2)
    cout << "Error reading velocities: end position is wrong "
         << blockSize << " vs " << blockSize2 << endl;

  // Skip for Gadget-2
  if (format == GADGET_2)
    gadget2 = readString(gStr, GADGET_2_SKIP);

  //
  // Read tags
  //
  readData(swap, (void*) &blockSize, GADGET_SKIP, 1, gStr);
  readData(swap, (void*) &tag[0], 
                         sizeof(ID_T), particleCount, gStr);
  readData(swap, (void*) &blockSize2, GADGET_SKIP, 1, gStr);
  if (blockSize != blockSize2)
    cout << "Error reading tags: end position is wrong "
         << blockSize << " vs " << blockSize2 << endl;

  //
  // Gas particles contain extra information
  // How do I use this in .cosmo format?
  //
  if (header.npart[GADGET_GAS] > 0) {

    // Skip for Gadget-2
    if (format == 2)
      gadget2 = readString(gStr, GADGET_2_SKIP);

    //
    // Read internal energy (u)
    //
    readData(swap, (void*) &blockSize, GADGET_SKIP, 1, gStr);
    readData(swap, (void*) &u[0], 
                           sizeof(POSVEL_T), header.npart[GADGET_GAS], gStr);
    readData(swap, (void*) &blockSize2, GADGET_SKIP, 1, gStr);
    if (blockSize != blockSize2)
      cout << "Error reading internal energy: end position is wrong "
           << blockSize << " vs " << blockSize2 << endl;

    // Skip for Gadget-2
    if (format == 2)
      gadget2 = readString(gStr, GADGET_2_SKIP);

    //
    // Read density (rho)
    //
    readData(swap, (void*) &blockSize, GADGET_SKIP, 1, gStr);
    readData(swap, (void*) &rho[0], 
                           sizeof(POSVEL_T), header.npart[GADGET_GAS], gStr);
    readData(swap, (void*) &blockSize2, GADGET_SKIP, 1, gStr);
    if (blockSize != blockSize2)
      cout << "Error reading density: end position is wrong "
           << blockSize << " vs " << blockSize2 << endl;

    // Skip for Gadget-2
    if (format == 2)
      gadget2 = readString(gStr, GADGET_2_SKIP);

    //
    // Read smoothing length (hsml)
    //
    readData(swap, (void*) &blockSize, GADGET_SKIP, 1, gStr);
    readData(swap, (void*) &hsml[0], 
                           sizeof(POSVEL_T), header.npart[GADGET_GAS], gStr);
    readData(swap, (void*) &blockSize2, GADGET_SKIP, 1, gStr);
    if (blockSize != blockSize2)
      cout << "Error reading smoothing length: end position is wrong "
           << blockSize << " vs " << blockSize2 << endl;
  }
  gStr->close();

  // Open the output Cosmo file
  ofstream *cStr = new ofstream(outFile.c_str(), ios::out | ios::binary);
  if (cStr->fail()) {
    cout << "File: " << outFile << " cannot be opened" << endl;
    exit(-1);
  }

#ifdef DEBUG
  POSVEL_T minLoc[DIMENSION], maxLoc[DIMENSION];
  POSVEL_T minVel[DIMENSION], maxVel[DIMENSION];
  ID_T minTag = particleCount;
  ID_T maxTag = -1;
  for (int dim = 0; dim < DIMENSION; dim++) {
    minLoc[dim] = MAX_FLOAT;
    maxLoc[dim] = MIN_FLOAT;
    minVel[dim] = MAX_FLOAT;
    maxVel[dim] = MIN_FLOAT;
  }
#endif

  ID_T indx = 0;
  ID_T tagindx = 0;
  for (int type = 0; type < NUM_GADGET_TYPES; type++) {

    // Multiply the Gadget mass by a factor to put in .cosmo format (M_sun/h)
    POSVEL_T mass = (POSVEL_T) header.mass[type] * massFactor;

    for (int i = 0; i < header.npart[type]; i++) {

      // Multiply Gadget location by a factor to put in .cosmo format  Mpc/h)
      for (int dim = 0; dim < DIMENSION; dim++)
        location[indx+dim] *= posFactor;

#ifdef DEBUG
      // Collect ranges on this file
      if (minTag > tag[tagindx]) minTag = tag[tagindx];
      if (maxTag < tag[tagindx]) maxTag = tag[tagindx];

      for (int dim = 0; dim < DIMENSION; dim++) {
        if (minLoc[dim] > location[indx+dim])
          minLoc[dim] = location[indx+dim];
        if (maxLoc[dim] < location[indx+dim])
          maxLoc[dim] = location[indx+dim];
        if (minVel[dim] > velocity[indx+dim])
          minVel[dim] = velocity[indx+dim];
        if (maxVel[dim] < velocity[indx+dim])
          maxVel[dim] = velocity[indx+dim];
      }
#endif

      cStr->write(reinterpret_cast<char*>(&location[indx]), sizeof(POSVEL_T));
      cStr->write(reinterpret_cast<char*>(&velocity[indx]), sizeof(POSVEL_T));
      cStr->write(reinterpret_cast<char*>(&location[indx+1]), sizeof(POSVEL_T));
      cStr->write(reinterpret_cast<char*>(&velocity[indx+1]), sizeof(POSVEL_T));
      cStr->write(reinterpret_cast<char*>(&location[indx+2]), sizeof(POSVEL_T));
      cStr->write(reinterpret_cast<char*>(&velocity[indx+2]), sizeof(POSVEL_T));
      cStr->write(reinterpret_cast<char*>(&mass), sizeof(POSVEL_T));
      cStr->write(reinterpret_cast<char*>(&tag[tagindx]), sizeof(ID_T));

      indx += DIMENSION;
      tagindx++;
    }
  }
  cStr->close();

#ifdef DEBUG
  // Ranges of location and velocity in file
  cout << endl;
  cout << "Number of particles: " << particleCount << endl;
  cout << "Location: ";
  for (int dim = 0; dim < DIMENSION; dim++)
    cout << " [" << minLoc[dim] << ":" << maxLoc[dim] << "] ";
  cout << endl;
  cout << "Velocity: ";
  for (int dim = 0; dim < DIMENSION; dim++)
    cout << " [" << minVel[dim] << ":" << maxVel[dim] << "] ";
  cout << endl;
  cout << "Tag:       [" << minTag << ":" << maxTag << "]" << endl << endl;
#endif
}
