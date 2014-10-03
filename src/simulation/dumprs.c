#include <stdlib.h>
#include <stdio.h>

typedef float POSVEL_T;
typedef int ID_T;
typedef unsigned short MASK_T;

int Np = 0;
POSVEL_T *xArr, *vxArr, *yArr, *vyArr, *zArr, *vzArr, *phiArr;
ID_T *idArr;
MASK_T *maskArr;

int allocParticles() {
  xArr = (POSVEL_T*) calloc(sizeof(POSVEL_T), Np);
  vxArr = (POSVEL_T*) calloc(sizeof(POSVEL_T), Np);

  yArr = (POSVEL_T*) calloc(sizeof(POSVEL_T), Np);
  vyArr = (POSVEL_T*) calloc(sizeof(POSVEL_T), Np);

  zArr = (POSVEL_T*) calloc(sizeof(POSVEL_T), Np);
  vzArr = (POSVEL_T*) calloc(sizeof(POSVEL_T), Np);

  phiArr = (POSVEL_T*) calloc(sizeof(POSVEL_T), Np);

  idArr = (ID_T*) calloc(sizeof(ID_T), Np);

  maskArr = (MASK_T*) calloc(sizeof(MASK_T), Np);

  if (!xArr || !vxArr || !yArr || !vyArr || !zArr || !vzArr ||
      !phiArr || !idArr || !maskArr) {
    fprintf(stderr, "malloc failed!\n");
    return 1;
  }

  return 0;
}

int readRestart( const char *inName ) {
  FILE *inFile = fopen(inName, "rb");
  if (!inFile) {
    fprintf(stderr, "fopen failed on %s!\n", inName);
    return 1;
  }

  fread(&Np, sizeof(int), 1, inFile);
  if (allocParticles() != 0) {
    return 1;
  }

  fread(&xArr[0], sizeof(POSVEL_T), Np, inFile);
  fread(&vxArr[0], sizeof(POSVEL_T), Np, inFile);

  fread(&yArr[0], sizeof(POSVEL_T), Np, inFile);
  fread(&vyArr[0], sizeof(POSVEL_T), Np, inFile);

  fread(&zArr[0], sizeof(POSVEL_T), Np, inFile);
  fread(&vzArr[0], sizeof(POSVEL_T), Np, inFile);

  fread(&phiArr[0], sizeof(POSVEL_T), Np, inFile);

  fread(&idArr[0], sizeof(ID_T), Np, inFile);

  fread(&maskArr[0], sizeof(MASK_T), Np, inFile);

  fclose(inFile);

  return 0;
}

int main(int argc, char *argv[]) {
	int i;
	const char *fn;
	if (argc < 2) {
		fprintf(stderr, "Usage: %s [restart file name]\n", argv[0]);
		return 1;
        }

	fn = argv[1];
	if (readRestart(fn) != 0) {
          return 1;
        }

	printf("# %d particles from %s\n", Np, fn);
	printf("#%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n",
		"x", "y", "z", "vx", "vy", "vz",
		"phi", "id", "mask");
	for (i = 0; i < Np; ++i) {
		printf("%g\t%g\t%g\t%g\t%g\t%g\t%g\t%d\t%hx\n",
			xArr[i], yArr[i], zArr[i],
			vxArr[i], vyArr[i], vzArr[i],
			phiArr[i], idArr[i], maskArr[i]);
	}

	return 0;
}

