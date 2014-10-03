#include "Domain.h"

#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#define FUDGE 0.0

using namespace std;

float Domain::m_grid2phys_pos = -1.0;
float Domain::m_phys2grid_pos = -1.0;
float Domain::m_grid2phys_vel = -1.0;
float Domain::m_phys2grid_vel = -1.0;

float Domain::m_rL = -1.0;
int Domain::m_ng = -1;
float Domain::m_oL = -1.0;
int Domain::m_ng_overload = -1;
int Domain::m_Ng_local_alive = -1;
int Domain::m_Ng_local_total = -1;
int Domain::m_maxNpLocal = -1;

float Domain::m_rL_local_alive[DIMENSION];
float Domain::m_rL_local_total[DIMENSION];
int Domain::m_ng_local_alive[DIMENSION];
int Domain::m_ng_local_total[DIMENSION];

int Domain::m_corner_grid_alive[DIMENSION];
int Domain::m_corner_grid_total[DIMENSION];
float Domain::m_corner_phys_alive[DIMENSION];
float Domain::m_corner_phys_total[DIMENSION];

void Domain::initialize(const Basedata &indat) {
  int i;

  int decompSize[DIMENSION];
  int myPosition[DIMENSION];
  Partition::getDecompSize(decompSize);
  Partition::getMyPosition(myPosition);

  m_rL = indat.rL();
  m_ng = indat.ng();

  m_phys2grid_pos = m_ng/m_rL;
  m_grid2phys_pos = m_rL/m_ng;

  //extra factor of H0 in Mpc/h to cancel km/s
  m_phys2grid_vel = m_ng/(100.0*m_rL);
  m_grid2phys_vel = 100.0*m_rL/m_ng;

  m_oL = indat.oL();
  m_ng_overload = static_cast< int >(ceilf(m_oL*m_phys2grid_pos));
  m_oL = m_ng_overload*m_grid2phys_pos;

  m_Ng_local_alive = 1;
  m_Ng_local_total = 1;
  for(i=0; i<DIMENSION; i++) {

    if(m_ng%decompSize[i] != 0) {
      fprintf(stderr,"ERROR: Domain decompSize does not divide ng\n");
      exit(-1);
    }

    m_rL_local_alive[i] = m_rL/decompSize[i];
    m_rL_local_total[i] = m_rL_local_alive[i] + 2.0*m_oL;
    
    m_ng_local_alive[i] = m_ng/decompSize[i];
    m_ng_local_total[i] = m_ng_local_alive[i] + 2*m_ng_overload + 1;

    m_Ng_local_alive *= m_ng_local_alive[i];
    m_Ng_local_total *= m_ng_local_total[i];

    m_corner_grid_alive[i] = myPosition[i]*m_ng_local_alive[i];
    m_corner_grid_total[i] = m_corner_grid_alive[i] - m_ng_overload;

    m_corner_phys_alive[i] = m_corner_grid_alive[i]*m_grid2phys_pos;
    m_corner_phys_total[i] = m_corner_grid_total[i]*m_grid2phys_pos;
  }

  int np = indat.np();
  float ppc = 1.0*np/m_ng;
  ppc = ppc*ppc*ppc;
  m_maxNpLocal = static_cast< int >(ceilf((1.0+FUDGE)*m_Ng_local_total*ppc));
}

void Domain::rL_local_alive(float rla[]) {
  int i;
  for(i=0; i<DIMENSION; i++)
    rla[i] = m_rL_local_alive[i];
  return;
}

void Domain::rL_local_total(float rlt[]) {
  int i;
  for(i=0; i<DIMENSION; i++)
    rlt[i] = m_rL_local_total[i];
  return;
}

void Domain::ng_local_alive(int ngla[]) {
  int i;
  for(i=0; i<DIMENSION; i++)
    ngla[i] = m_ng_local_alive[i];
  return;
}

void Domain::ng_local_total(int nglt[]) {
  int i;
  for(i=0; i<DIMENSION; i++)
    nglt[i] = m_ng_local_total[i];
  return;
}

void Domain::corner_grid_alive(int corner[]) {
  int i;
  for(i=0; i<DIMENSION; i++)
    corner[i] = m_corner_grid_alive[i];
  return;
}

void Domain::corner_grid_total(int corner[]) {
  int i;
  for(i=0; i<DIMENSION; i++)
    corner[i] = m_corner_grid_total[i];
  return;
}

void Domain::corner_phys_alive(float corner[]) {
  int i;
  for(i=0; i<DIMENSION; i++)
    corner[i] = m_corner_phys_alive[i];
  return;
}

void Domain::corner_phys_total(float corner[]) {
  int i;
  for(i=0; i<DIMENSION; i++)
    corner[i] = m_corner_phys_total[i];
  return;
}
