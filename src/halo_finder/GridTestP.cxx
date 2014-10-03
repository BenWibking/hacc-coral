/*=========================================================================

  Program:   Grid exchange driver program
  Module:    $RCSfile: GridTestP.cxx,v $

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/
//
// Test the parallel exchange of grids with neighbors which is not quite
// the standard ghost cell exchange.  In this case each processor has
// a contiguous array with an alive portion in the middle surrounded by
// a dead portion where the size on the front part of a dimension is not
// necessarily the size on the back.  When data is sent to a neighbor it
// comes from the alive part of the array, but when it is received it is
// unpacked into the dead part.
//

#include <iostream>
#include <iomanip>
#include "Partition.h"
#include "GridExchange.h"

#include <stdlib.h>
#include <rru_mpi.h>

using namespace std;

int main(int argc, char* argv[])
{
  if (argc != 5) {
    cout << "Usage: mpirun -np # GridTestP totalSize dead0 dead1 proc" << endl;
  }

  // Size of the data array on a side
  int cnt = 1;
  int size = atoi(argv[cnt++]);

  // Size of dead regions on the front and back of each dimension
  int dead0 = atoi(argv[cnt++]);
  int dead1 = atoi(argv[cnt++]);
  int printProc = atoi(argv[cnt++]);

  // Initialize the partitioner which uses MPI Cartesian Topology
  //Partition::initialize(argc, argv);
  MPI_Init(&argc, &argv);
  Partition::initialize();

  // Position of each rank within the decomposition
  int rank = Partition::getMyProc();
  int layoutPos[3], layoutSize[3];
  Partition::getMyPosition(layoutPos);
  Partition::getDecompSize(layoutSize);
  cout << "Rank " << rank << " Pos " << layoutPos[0] << "," 
       << layoutPos[1] << "," << layoutPos[2] << endl;
  MPI_Barrier(MPI_COMM_WORLD);

  // Allocate data to share
  int totalSize[DIMENSION];
  int dataSize = 1;

  for (int dim = 0; dim < DIMENSION; dim++) {
    totalSize[dim] = size / layoutSize[dim];
    dataSize *= totalSize[dim];
  }

  GRID_T* data = new GRID_T[dataSize];
  for (int i = 0; i < dataSize; i++)
    data[i] = -9;

  int planeSize = totalSize[1] * totalSize[2];
  int rowSize = totalSize[2];

  for (int i = 0; i < totalSize[0]; i++) {
    for (int j = 0; j < totalSize[1]; j++) {
      for (int k = 0; k < totalSize[2]; k++) {
        if (k >= dead0 && k < (totalSize[2] - dead1) &&
            j >= dead0 && j < (totalSize[1] - dead1) &&
            i >= dead0 && i < (totalSize[0] - dead1)) {
              int index = (i * planeSize) + (j * rowSize) + k;
              data[index] = rank + (index * .001);
        }
      }
    }
  }

  if (rank == printProc) {
    cout << "BEFORE (Z Planes front to back)" << endl;
    for (int i = 0; i < totalSize[0]; i++) {
      cout << endl << " Plane " << i << endl << endl;
      for (int j = totalSize[1] - 1; j >= 0; j--) {
        for (int k = 0; k < totalSize[2]; k++) {
          int index = (i * planeSize) + (j * rowSize) + k;
          cout << "  " << setprecision(5) << setw(6) << data[index];
        }
        cout << endl;
      }
      cout << endl << endl << endl;
    }
  }

  if (rank == printProc) {
    cout << "BEFORE (X Planes front to back)" << endl;
    for (int k = 0; k < totalSize[2]; k++) {
      cout << endl << " Plane " << k << endl << endl;
      for (int j = totalSize[1] - 1; j >= 0; j--) {
        for (int i = 0; i < totalSize[0]; i++) {
          int index = (i * planeSize) + (j * rowSize) + k;
          cout << "  " << setprecision(5) << setw(6) << data[index];
        }
        cout << endl;
      }
      cout << endl << endl << endl;
    }
  }

  // Construct the grid exchanger
  GridExchange exchange(totalSize, dead0, dead1);

  // Exchange grid
  exchange.exchangeGrid(data);

  if (rank == printProc) {
    cout << "AFTER (Z Planes front to back)" << endl;
    for (int i = 0; i < totalSize[0]; i++) {
      cout << endl << " Plane " << i << endl << endl;
      for (int j = totalSize[1] - 1; j >= 0; j--) {
        for (int k = 0; k < totalSize[2]; k++) {
          int index = (i * planeSize) + (j * rowSize) + k;
          cout << "  " << setprecision(5) << setw(6) << data[index];
        }
        cout << endl;
      }
      cout << endl << endl << endl;
    }
  }

  if (rank == printProc) {
    cout << "AFTER (X Planes front to back)" << endl;
    for (int k = 0; k < totalSize[2]; k++) {
      cout << endl << " Plane " << k << endl << endl;
      for (int j = totalSize[1] - 1; j >= 0; j--) {
        for (int i = 0; i < totalSize[0]; i++) {
          int index = (i * planeSize) + (j * rowSize) + k;
          cout << "  " << setprecision(5) << setw(6) << data[index];
        }
        cout << endl;
      }
      cout << endl << endl << endl;
    }
  }

  // Shut down MPI
  Partition::finalize();
  MPI_Finalize();

  return 0;
}
