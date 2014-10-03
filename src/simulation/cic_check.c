#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>

#define MAX_STR_LEN 1024

int main(int argc, char *argv[]) {
  char *refName, *inName, tmpStr[MAX_STR_LEN];
  FILE *refFile, *inFile;
  int ngRef, ngIn, NgRef, NgIn;
  int i;
  float absMin=1.0e10, absMax=0.0;
  float *refArr;
  int inIndx, absMaxIndx;
  float inVal, absDiff, absMaxVal, absMaxRef;

  if(argc < 5) {
    fprintf(stderr,"USAGE: cic_check <refName> <NgRef> <inName> <NgIn>\n");
    exit(-1);
  }

  refName = argv[1];
  NgRef = atoi(argv[2]);
  inName = argv[3];
  NgIn = atoi(argv[4]);

  //NgRef = ngRef*ngRef*ngRef;
  //NgIn = ngIn*ngIn*ngIn;

  refArr = (float *)malloc(NgRef*sizeof(float));

  refFile = fopen(refName,"r");
  for(i=0; i<NgRef; i++) {
    fscanf(refFile,"%s",tmpStr);
    fscanf(refFile,"%s",tmpStr);
    fscanf(refFile,"%s",tmpStr);
    fscanf(refFile,"%s",tmpStr);

    fscanf(refFile,"%s",tmpStr);
    refArr[i] = atof(tmpStr);
  }
  fclose(refFile);

  inFile = fopen(inName,"r");
  for(i=0; i<NgIn; i++) {
    fscanf(inFile,"%s",tmpStr);
    inIndx = atoi(tmpStr);

    fscanf(inFile,"%s",tmpStr);
    fscanf(inFile,"%s",tmpStr);
    fscanf(inFile,"%s",tmpStr);

    fscanf(inFile,"%s",tmpStr);
    inVal = atof(tmpStr);

    absDiff = fabs(inVal - refArr[inIndx]);
    if(absDiff > absMax) {
      absMax = absDiff;
      absMaxIndx = inIndx;
      absMaxVal = inVal;
      absMaxRef = refArr[inIndx];
    }
    if(absDiff < absMin) {
      absMin = absDiff;
    }
  }
  fclose(inFile);

  printf("absMin = %f\n",absMin);
  printf("absMax = %f (%d, %f, %f)\n",absMax,absMaxIndx,absMaxRef,absMaxVal);

  free(refArr);

  return 0;
}
