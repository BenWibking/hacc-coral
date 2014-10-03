#include "Particles.h"
#include "RCBForceTree.h"
#include "RCOForceTree.h"


int my_posix_memalign(void **memptr, size_t alignment, size_t size) {
  int ret;

  ret = posix_memalign(memptr, alignment, size);
  assert(ret==0);

  return ret;
}
//int my_posix_memalign(void **memptr, size_t alignment, size_t size);



/*
FGrid::FGrid() :
  m_b(0.72),
  m_c(0.01),
  m_d(0.27),
  m_e(0.0001),
  m_f(360.0),
  m_g(100.0),
  m_h(0.67),
  m_l(17.0),
  m_rmax(3.116326355)
{};



float FGrid::fgor(float r) {
  float f0 = m_c + 2.0/3.0*m_b*m_b*m_b;
  float r2 = r*r;
  float r4 = r2*r2;
  float r6 = r4*r2;
  float coshbr = coshf(m_b*r);
  float r3fgor = tanhf(m_b*r) - m_b*r/coshbr/coshbr
    + m_c*r*r2*(1.0 + m_d*r2)*expf(-1.0*m_d*r2)
    + m_e*r2*(m_f*r2 + m_g*r4 + m_l*r6)*expf(-1.0*m_h*r2);
  float rdiv = r + 1.0*(r<=0.0);
  float rdiv3 = rdiv*rdiv*rdiv;
  return (r3fgor/rdiv3 + (r<=0.0)*f0)*(r<m_rmax)*(r>=0.0);
}



void FGrid::fgor_r2_interp(int nInterp, float **r2, float **f) {
  my_posix_memalign((void **)r2, MEM_ALIGN, nInterp*sizeof(float) );
  my_posix_memalign((void **)f, MEM_ALIGN, nInterp*sizeof(float) );

  double dr2 = (m_rmax*m_rmax)/(nInterp-1.0);
  for(int i=0; i<nInterp; i++) {
    (*r2)[i] = i*dr2;
    (*f)[i] = fgor(sqrt(i*dr2));
  }

  return;
}



FGridEvalFit::FGridEvalFit(FGrid *fg) {
  m_fg = fg;
}



float FGridEvalFit::eval(float r2) {
  return m_fg->fgor(r2);
}



FGridEvalInterp::FGridEvalInterp(FGrid *fg, int nInterp) {
  m_nInterp = nInterp;
  fg->fgor_r2_interp(m_nInterp, &m_r2, &m_f);
  m_r2min = m_r2[0];
  m_r2max = m_r2[m_nInterp-1];
  m_dr2 = (m_r2max - m_r2min)/(m_nInterp - 1.0);
}



FGridEvalInterp::~FGridEvalInterp() {
  free(m_r2);
  free(m_f);
}


float FGridEvalInterp::eval(float r2) {
  int inRange, indx;
  float inRangef;
  inRange = (r2 > m_r2min)*(r2 < m_r2max);
  inRangef = 1.0*inRange;
  indx = int((r2 - m_r2min)/m_dr2)*inRange;
  return inRangef*(m_f[indx]+(r2-m_r2[indx])/m_dr2*(m_f[indx+1]-m_f[indx]));
}
*/



void Particles::forceInterp(uint32_t nInterp) {
  /*
  delete m_fl;
  delete m_fgore;
  */

  m_fgore = new FGridEvalPoly(m_fg);
  m_fl = new ForceLawSR(m_fgore, m_rsm);

  return;
}

void Particles::allocField() {
  //m_field = new GRID_T[Domain::Ng_local_total() + 1];
  my_posix_memalign( (void **)&m_field, MEM_ALIGN, 
		     (Domain::Ng_local_total()+1)*sizeof(GRID_T) );
}

Particles::Particles(const Basedata & bdata, MC3Options & options) :
  m_xArr(NULL),
  m_yArr(NULL),
  m_zArr(NULL),
  m_vxArr(NULL),
  m_vyArr(NULL),
  m_vzArr(NULL),
  m_phiArr(NULL),
  m_idArr(NULL),
  m_maskArr(NULL),
  m_field(NULL),
  m_Np_local_total( 0 ),
  m_alpha(-1.0),
  m_gpscal(-1.0),
  m_nsub(-1),
  m_edge(-1.0),
  m_rsm(-1.0),
  m_fsrrmax(-1.0),
  m_cmsize(-1.0),
  m_openAngle(-1.0),
  m_skipStreamQ(0),
  m_skipKickSRQ(0),
  m_useFastTreeEval(0),
  m_useRCBTree(0),
  m_rcbTreeExtraLevels(0),
  m_rcbTreePPN(0),
  m_rcbTreeTaskPartMin(0),
  m_mpiio(0)
{
  m_alpha = bdata.alpha();
  m_gpscal = ( static_cast< float >( bdata.ng() ) ) /
    ( static_cast< float >( bdata.np() ) );

  m_nsub = bdata.nsub();
  m_edge = bdata.edge();
  m_rsm = bdata.rsm();
  m_cmsize = bdata.cmsize();
  m_openAngle = bdata.openAngle();

  m_fg = new FGrid();
  m_fsrrmax = m_fg->rmax();

  for(int i=0; i<N_POSVEL_T; i++)
    m_pvData[i] = NULL;

  int nglt[DIMENSION];
  Domain::ng_local_total(nglt);
  m_cm = new CMLite( (int)ceilf(nglt[0]/m_cmsize),
		     (int)ceilf(nglt[1]/m_cmsize),
		     (int)ceilf(nglt[2]/m_cmsize) );

  m_fgore = new FGridEvalFit(m_fg);
  m_fl = new ForceLawSR(m_fgore, m_rsm);

  if(options.interpQ()) {
    forceInterp(options.nInterp());
  }

  if(options.polyQ()) {
    m_fgore = new FGridEvalPoly(m_fg);
    m_fl = new ForceLawSR(m_fgore, m_rsm);
  }

  m_skipStreamQ = options.skipStreamQ();
  m_skipKickSRQ = options.skipKickSRQ();

  m_useFastTreeEval = options.useFastTreeEval();
  m_useRCBTree = options.useRCBTree();
  m_rcbTreeExtraLevels = options.rcbTreeExtraLevels();
  if (!m_rcbTreeExtraLevels) {
    m_rcbTreeExtraLevels = 2;
  }

  m_rcbTreePPN = options.rcbTreePPN();
  if (!m_rcbTreePPN) {
    if (m_useRCBTree == 2 || m_useRCBTree == 5) {
      // default for the monopole mode
      m_rcbTreePPN = 2;
    } else {
      // default for the quadrupole mode
      m_rcbTreePPN = 12;
    }
  }

  m_rcbTreeTaskPartMin = options.rcbTreeTaskPartMin();
  if (!m_rcbTreeTaskPartMin) {
    m_rcbTreeTaskPartMin = 128;
  }

  m_mpiio = options.mpiio();
}



Particles::~Particles() {
  dropParticles();
  //if(m_field) delete m_field;
  if(m_field) free(m_field);
  delete m_cm;
}



void Particles::dropParticles() {

  for(int i=0; i<N_POSVEL_T; i++) {
    if(m_pvData[i]) {
      free(m_pvData[i]);
      m_pvData[i] = NULL;
    }
  }
  updatePointers();

  if(m_idArr) {
    free(m_idArr);
    m_idArr = NULL;
  }

  if(m_maskArr) {
    free(m_maskArr);
    m_maskArr = NULL;
  }

  m_Np_local_total = 0;

  m_coords_localQ = 0;

  return;
}



void Particles::updatePointers() {
  m_xArr = m_pvData[0];
  m_yArr = m_pvData[1];
  m_zArr = m_pvData[2];
  m_vxArr = m_pvData[3];
  m_vyArr = m_pvData[4];
  m_vzArr = m_pvData[5];
  m_phiArr = m_pvData[6];
  m_massArr = m_pvData[7];

  return;
}

void Particles::updatePointers2() {
  m_pvData[0] = m_xArr;
  m_pvData[1] = m_yArr;
  m_pvData[2] = m_zArr;
  m_pvData[3] = m_vxArr;
  m_pvData[4] = m_vyArr;
  m_pvData[5] = m_vzArr;
  m_pvData[6] = m_phiArr;
  m_pvData[7] = m_massArr;

  return;
}



void Particles::allocParticles(int Np) {
  dropParticles();
  m_Np_local_total = Np;

  for(int i=0; i<N_POSVEL_T; i++)
    my_posix_memalign((void **)&(m_pvData[i]), MEM_ALIGN, Np*sizeof(POSVEL_T));
  updatePointers();

  my_posix_memalign( (void **)&m_idArr, MEM_ALIGN, Np*sizeof(ID_T) );
  my_posix_memalign( (void **)&m_maskArr, MEM_ALIGN, Np*sizeof(MASK_T) );

  return;
}



void Particles::copyAllFromVectors(vector<POSVEL_T> **xx,
                                   vector<POSVEL_T> **yy,
                                   vector<POSVEL_T> **zz,
                                   vector<POSVEL_T> **vx,
                                   vector<POSVEL_T> **vy,
                                   vector<POSVEL_T> **vz,
                                   vector<POSVEL_T> **phi,
                                   vector<ID_T> **id,
                                   vector<MASK_T> **mask,
                                   float anow,
                                   int deleteVectorsQ) {
  int Np = (*xx)->size();

  dropParticles();

  m_Np_local_total = Np;

  vector<POSVEL_T> *mass = new vector<POSVEL_T>;
  mass->reserve(Np);
  for(int i=0; i<Np; i++)
    mass->push_back(1.0);

  vector<POSVEL_T> *pv[N_POSVEL_T];
  pv[0] = *xx;
  pv[1] = *yy;
  pv[2] = *zz;
  pv[3] = *vx;
  pv[4] = *vy;
  pv[5] = *vz;
  pv[6] = *phi;
  pv[7] = mass;

  for(int i=0; i<N_POSVEL_T; i++) {
    my_posix_memalign((void **)&(m_pvData[i]), MEM_ALIGN, Np*sizeof(POSVEL_T));
    memcpy(&(m_pvData[i][0]), &(*pv[i])[0], Np*sizeof(POSVEL_T));
    if(deleteVectorsQ)
      delete pv[i];
  }
  updatePointers();

  my_posix_memalign((void **)&m_idArr, MEM_ALIGN, Np*sizeof(ID_T));
  memcpy(&m_idArr[0], &(**id)[0], Np*sizeof(ID_T));
  if(deleteVectorsQ)
    delete *id;

  my_posix_memalign((void **)&m_maskArr, MEM_ALIGN, Np*sizeof(MASK_T));
  memcpy(&m_maskArr[0], &(**mask)[0], Np*sizeof(MASK_T));
  if(deleteVectorsQ)
    delete *mask;

  if(deleteVectorsQ) {
    //*xx = *yy = *zz = *vx = *vy = *vz = *phi = *id = *mask = NULL;
    *xx = NULL;
    *yy = NULL;
    *zz = NULL;
    *vx = NULL;
    *vy = NULL;
    *vz = NULL;
    *phi = NULL;
    *id = NULL;
    *mask = NULL;
  } else {
    delete mass;
  }

  m_coords_localQ = 0;
  coords_global2local(anow);

  return;
}



//MAYBE SHOULD CHECK BOUNDS FOR EACH DIMENSION
inline int 
Particles::array_index(int xx, 
		       int yy, 
		       int zz, 
		       int ng[DIMENSION], 
		       int lo[DIMENSION],
		       int hi[DIMENSION],
		       int safe) {
/*                       
  int indx, inbounds=1;
  indx = _l_array_index(xx,yy,zz,ng[1],ng[2]);
  inbounds *= (xx >= lo[0]);
  inbounds *= (xx < hi[0]);
  inbounds *= (yy >= lo[1]);
  inbounds *= (yy < hi[1]);
  inbounds *= (zz >= lo[2]);
  inbounds *= (zz < hi[2]);
  return indx*inbounds + safe*(1-inbounds);
*/
  if ( ( xx >= lo[0] ) && ( xx < hi[0] ) && ( yy >= lo[1] ) && ( yy < hi[1] ) && ( zz >= lo[2] ) && ( zz < hi[2] ) )
      return ( xx * ng[1] + yy ) * ng[2] + zz;
  else
      return safe;
}



#define MAX(a, b) (a > b ? a : b)
#define MIN(a, b) (a < b ? a : b)

void Particles::resortParticles() {
  assert(m_coords_localQ == 1);

  int indx;
  int xx, yy, zz, cnt;

  //allocated
  int *permArr;
  void *tmpArr;

  //cast
  int *cntArr, *indxArr;
  POSVEL_T *srcPV[N_POSVEL_T], *tmpPV;
  MASK_T *tmpMask;
  ID_T *tmpID;

  int cumulative, tmp, Np, ng[DIMENSION], Ng;
  int zero[DIMENSION] = {0,0,0};

  Np = m_Np_local_total;
  Domain::ng_local_total(ng);
  Ng = Domain::Ng_local_total();

  my_posix_memalign( (void **)&permArr, MEM_ALIGN, Np*sizeof(int) );
  my_posix_memalign(&tmpArr,MEM_ALIGN,Np*MAX(sizeof(ID_T),sizeof(POSVEL_T)));
  indxArr = (int *)tmpArr;
  cntArr = (int *)m_field;

  memset(cntArr, 0, (Ng+1)*sizeof(int));
  memset(indxArr, 0, Np*sizeof(int));
  memset(permArr, 0, Np*sizeof(int));

  for(int i=0; i<Np; i++) {
    xx = static_cast< int >(FLOOR(m_xArr[i]));
    yy = static_cast< int >(FLOOR(m_yArr[i]));
    zz = static_cast< int >(FLOOR(m_zArr[i]));
    indx = array_index(xx, yy, zz, ng, zero, ng, Ng);
    indxArr[i] = indx;
    cntArr[indx]++;
  }

  m_Np_last = m_Np_local_total - cntArr[Ng];

  cumulative = Np;
  for(int i=Ng; i>=0; i--) {
    tmp = cntArr[i];
    cntArr[i] = cumulative-tmp;
    cumulative -= tmp;
  }

  for(int i=0; i<Np; i++) {
    permArr[i] = cntArr[indxArr[i]];
    cntArr[indxArr[i]]++;
  }

  srcPV[0] = m_xArr;
  srcPV[1] = m_yArr;
  srcPV[2] = m_zArr;
  srcPV[3] = m_vxArr;
  srcPV[4] = m_vyArr;
  srcPV[5] = m_vzArr;
  srcPV[6] = m_phiArr;
  srcPV[7] = m_massArr;

  tmpPV = (POSVEL_T *)tmpArr;
  for(int j=0; j<N_POSVEL_T; j++) {
    for(int i=0; i<Np; i++)
      tmpPV[permArr[i]] = srcPV[j][i];
    memcpy(srcPV[j], tmpPV, Np*sizeof(POSVEL_T));
  }

  tmpID = (ID_T *)tmpArr;
  for(int i=0; i<Np; i++)
    tmpID[permArr[i]] = m_idArr[i];
  memcpy(m_idArr, tmpID, Np*sizeof(ID_T));

  tmpMask = (MASK_T *)tmpArr;
  for(int i=0; i<Np; i++)
    tmpMask[permArr[i]] = m_maskArr[i];
  memcpy(m_maskArr, tmpMask, Np*sizeof(MASK_T));

  memset(cntArr, 0, (Ng+1)*sizeof(int));
  free(permArr);
  free(tmpArr);

  return;
}



void Particles::buildChainingMesh() {
  assert(m_coords_localQ == 1);

  int indx;
  int xx, yy, zz, cnt;

  //allocated
  int *permArr;
  void *tmpArr;

  //cast
  int *cntArr, *indxArr;
  POSVEL_T *srcPV[N_POSVEL_T], *tmpPV;
  MASK_T *tmpMask;
  ID_T *tmpID;

  int cumulative, tmp, Np, ng[DIMENSION], Ng;
  int zero[DIMENSION] = {0,0,0};

  Np = m_Np_local_total;

  ng[0] = m_cm->ng[0];
  ng[1] = m_cm->ng[1];
  ng[2] = m_cm->ng[2];
  Ng = m_cm->Ng;

  my_posix_memalign( (void **)&permArr, MEM_ALIGN, Np*sizeof(int) );
  my_posix_memalign(&tmpArr,MEM_ALIGN,Np*MAX(sizeof(ID_T),sizeof(POSVEL_T)));
  indxArr = (int *)tmpArr;
  cntArr = (int *)m_field;

  memset(cntArr, 0, (Ng+1)*sizeof(int));
  memset(indxArr, 0, Np*sizeof(int));
  memset(permArr, 0, Np*sizeof(int));

  int nNg = 0;
  for(int i=0; i<Np; i++) {
    xx = static_cast< int >(FLOOR(m_xArr[i]/m_cmsize));
    yy = static_cast< int >(FLOOR(m_yArr[i]/m_cmsize));
    zz = static_cast< int >(FLOOR(m_zArr[i]/m_cmsize));
    indx = array_index(xx, yy, zz, ng, zero, ng, Ng);
    nNg += (indx==Ng);
    indxArr[i] = indx;
    cntArr[indx]++;
  }

  cumulative = Np - nNg;
  m_cm->indxhi[Ng] = Np;
  m_cm->indxlo[Ng] = cumulative;
  for(int i=Ng-1; i>=0; i--) {
    tmp = cntArr[i];
    cntArr[i] = cumulative-tmp;
    m_cm->indxhi[i] = cumulative;
    cumulative -= tmp;
    m_cm->indxlo[i] = cumulative;
  }

  for(int i=0; i<Np; i++) {
    permArr[i] = cntArr[indxArr[i]];
    cntArr[indxArr[i]]++;
  }

  srcPV[0] = m_xArr;
  srcPV[1] = m_yArr;
  srcPV[2] = m_zArr;
  srcPV[3] = m_vxArr;
  srcPV[4] = m_vyArr;
  srcPV[5] = m_vzArr;
  srcPV[6] = m_phiArr;
  srcPV[7] = m_massArr;

  tmpPV = (POSVEL_T *)tmpArr;
  for(int j=0; j<N_POSVEL_T; j++) {
    for(int i=0; i<Np; i++)
      tmpPV[permArr[i]] = srcPV[j][i];
    memcpy(srcPV[j], tmpPV, Np*sizeof(POSVEL_T));
  }

  tmpID = (ID_T *)tmpArr;
  for(int i=0; i<Np; i++)
    tmpID[permArr[i]] = m_idArr[i];
  memcpy(m_idArr, tmpID, Np*sizeof(ID_T));

  tmpMask = (MASK_T *)tmpArr;
  for(int i=0; i<Np; i++)
    tmpMask[permArr[i]] = m_maskArr[i];
  memcpy(m_maskArr, tmpMask, Np*sizeof(MASK_T));

  memset(cntArr, 0, Ng*sizeof(int));
  free(permArr);
  free(tmpArr);

  return;
}



void Particles::cic() {
  assert(m_coords_localQ == 1);

  int ix, iy, iz, ip, jp, kp;
  int nn, ng[DIMENSION], Ng, np;
  int zero[DIMENSION] = {0,0,0};
  POSVEL_T xx, yy, zz;
  POSVEL_T ab, de, gh;
  POSVEL_T c;
  POSVEL_T *x, *y, *z;
  GRID_T *rhoArr;
  int safe;

  Domain::ng_local_total(ng);
  Ng = Domain::Ng_local_total();
  np = m_Np_local_total;
  c = m_gpscal*m_gpscal*m_gpscal;

  x = &m_xArr[0];
  y = &m_yArr[0];
  z = &m_zArr[0];
  rhoArr = m_field;
  safe = Ng;

  memset(&rhoArr[0], 0, sizeof(GRID_T) * (Ng+1) );

  for(nn=0; nn < np; nn++) {
    xx = x[nn];
    yy = y[nn];
    zz = z[nn];

    ix = static_cast< int >(FLOOR(xx));
    iy = static_cast< int >(FLOOR(yy));
    iz = static_cast< int >(FLOOR(zz));

    ip = ix+1;
    jp = iy+1;
    kp = iz+1;

    ab = 1.0 + (ix-xx);
    de = 1.0 + (iy-yy);
    gh = 1.0 + (iz-zz);

    rhoArr[array_index(ix,iy,iz,ng,zero,ng,safe)]+=c*ab*de*gh;
    rhoArr[array_index(ix,jp,iz,ng,zero,ng,safe)]+=c*ab*(1.0-de)*gh;
    rhoArr[array_index(ix,jp,kp,ng,zero,ng,safe)]+=c*ab*(1.0-de)*(1.0-gh);
    rhoArr[array_index(ix,iy,kp,ng,zero,ng,safe)]+=c*ab*de*(1.0-gh);
    rhoArr[array_index(ip,iy,kp,ng,zero,ng,safe)]+=c*(1.0-ab)*de*(1.0-gh);
    rhoArr[array_index(ip,jp,kp,ng,zero,ng,safe)]+=c*(1.0-ab)*(1.0-de)*(1.0-gh);
    rhoArr[array_index(ip,jp,iz,ng,zero,ng,safe)]+=c*(1.0-ab)*(1.0-de)*gh;
    rhoArr[array_index(ip,iy,iz,ng,zero,ng,safe)]+=c*(1.0-ab)*de*gh;
  }

  return;
}



void Particles::inverse_cic(float tau, float fscal, int comp) {
  assert(m_coords_localQ == 1);

  int ix, iy, iz, ip, jp, kp;
  int nn, ng[DIMENSION], Ng, np;
  int zero[DIMENSION] = {0,0,0};
  POSVEL_T xx, yy, zz;
  POSVEL_T ab, de, gh;
  POSVEL_T f;
  POSVEL_T *x, *y, *z;
  POSVEL_T *vels[DIMENSION+1];
  GRID_T *grad_phi;
  int safe;

  Domain::ng_local_total(ng);
  Ng = Domain::Ng_local_total();
  np = m_Np_local_total;

  x = &m_xArr[0];
  y = &m_yArr[0];
  z = &m_zArr[0];

  vels[0] = &m_vxArr[0];
  vels[1] = &m_vyArr[0];
  vels[2] = &m_vzArr[0];

  vels[3] = &m_phiArr[0];

  grad_phi = m_field;
  safe = Ng;
  grad_phi[safe] = 0.0;

#ifdef _OPENMP
#pragma omp parallel for private(xx,yy,zz,ix,iy,iz,ip,jp,kp,ab,de,gh,f)
#endif
  for(nn=0; nn < np; nn++) {
    xx = x[nn];
    yy = y[nn];
    zz = z[nn];

    ix = static_cast< int >(FLOOR(xx));
    iy = static_cast< int >(FLOOR(yy));
    iz = static_cast< int >(FLOOR(zz));

    ip = ix+1;
    jp = iy+1;
    kp = iz+1;

    ab = 1.0 + (ix-xx);
    de = 1.0 + (iy-yy);
    gh = 1.0 + (iz-zz);

    f = 0;

    f += grad_phi[array_index(ix,iy,iz,ng,zero,ng,safe)]*ab*de*gh;
    f += grad_phi[array_index(ix,jp,iz,ng,zero,ng,safe)]*ab*(1.0-de)*gh;
    f += grad_phi[array_index(ix,jp,kp,ng,zero,ng,safe)]*ab*(1.0-de)*(1.0-gh);
    f += grad_phi[array_index(ix,iy,kp,ng,zero,ng,safe)]*ab*de*(1.0-gh);
    f += grad_phi[array_index(ip,iy,kp,ng,zero,ng,safe)]*(1.0-ab)*de*(1.0-gh);
    f += grad_phi[array_index(ip,jp,kp,ng,zero,ng,safe)]*(1.0-ab)*(1.0-de)*(1.0-gh);
    f += grad_phi[array_index(ip,jp,iz,ng,zero,ng,safe)]*(1.0-ab)*(1.0-de)*gh;
    f += grad_phi[array_index(ip,iy,iz,ng,zero,ng,safe)]*(1.0-ab)*de*gh;

    vels[comp][nn] += f*fscal*tau;
  }

  return;
}



void Particles::inverse_cic_potential() {
  assert(m_coords_localQ == 1);

  //set m_phiArr=0
  memset(m_phiArr, 0, m_Np_local_total*sizeof(GRID_T));

  //call inverse_cic
  inverse_cic(1.0, 1.0, 3);

  return;
}



void Particles::map1(float pp, float tau, float adot) {
  assert(m_coords_localQ == 1);

  POSVEL_T *x, *y, *z, *vx, *vy, *vz;
  int i;

  x = m_xArr;
  y = m_yArr;
  z = m_zArr;
  vx = m_vxArr;
  vy = m_vyArr;
  vz = m_vzArr;

  const float pf = powf(pp, (1.0 + 1.0 / m_alpha) );
  const float prefactor = 1.0 / (m_alpha * adot * pf);

#ifdef _OPENMP
#pragma omp parallel for
#endif
  for (i = 0; i < m_Np_local_total; i++) {
    x[i] = x[i] + prefactor * tau * vx[i];
    y[i] = y[i] + prefactor * tau * vy[i];
    z[i] = z[i] + prefactor * tau * vz[i];
  }

  return;
}



void Particles::coords_global2local(float anow) {
  int i;
  float phys2grid_pos, phys2grid_vel, pre_vel;
  int corner[DIMENSION];

  if(m_coords_localQ == 0) {

    phys2grid_pos = Domain::phys2grid_pos();
    phys2grid_vel = Domain::phys2grid_vel();
    Domain::corner_grid_total(corner);
    
    pre_vel = phys2grid_vel*anow*anow;
    
    for(i=0; i<m_Np_local_total; i++) {
      m_xArr[i] *= phys2grid_pos;
      m_xArr[i] -= corner[0];
      
      m_yArr[i] *= phys2grid_pos;
      m_yArr[i] -= corner[1];
      
      m_zArr[i] *= phys2grid_pos;
      m_zArr[i] -= corner[2];
      
      m_vxArr[i] *= pre_vel;
      m_vyArr[i] *= pre_vel;
      m_vzArr[i] *= pre_vel;
    }

    m_coords_localQ = 1;
  }

  return;
}



void Particles::coords_local2global(float anow) {
  int i;
  float grid2phys_pos, grid2phys_vel, pre_vel;
  int corner[DIMENSION];

  if(m_coords_localQ == 1) {

    grid2phys_pos = Domain::grid2phys_pos();
    grid2phys_vel = Domain::grid2phys_vel();
    Domain::corner_grid_total(corner);
    
    pre_vel = grid2phys_vel/anow/anow;
    
    for(i=0; i<m_Np_local_total; i++) {
      m_xArr[i] += corner[0];
      m_xArr[i] *= grid2phys_pos;
      
      m_yArr[i] += corner[1];
      m_yArr[i] *= grid2phys_pos;
      
      m_zArr[i] += corner[2];
      m_zArr[i] *= grid2phys_pos;
      
      m_vxArr[i] *= pre_vel;
      m_vyArr[i] *= pre_vel;
      m_vzArr[i] *= pre_vel;
    }

    m_coords_localQ = 0;
  }

  return;
}



int Particles::Np_local_alive() {
  assert(m_coords_localQ == 1);

  int npla=0, i, ngla[DIMENSION];
  POSVEL_T lo, hi[DIMENSION];

  lo = 1.0*Domain::ng_overload();
  Domain::ng_local_alive(ngla);
  for(i=0; i<DIMENSION; i++)
    hi[i] = lo + 1.0*ngla[i];

  for(i=0; i<m_Np_local_total; i++) {
    npla += (m_xArr[i]>=lo)*(m_xArr[i]<hi[0])*(m_yArr[i]>=lo)*(m_yArr[i]<hi[1])*(m_zArr[i]>=lo)*(m_zArr[i]<hi[2]);
  }

  return npla;
}



void Particles::writeAliveHCosmo( const char *outName, float anow ) {
  int i;
  int Npla = Np_local_alive();

  coords_local2global(anow);

  float ca[DIMENSION];
  Domain::corner_phys_alive(ca);
  float rla[DIMENSION];
  Domain::rL_local_alive(rla);
  float lo[DIMENSION], hi[DIMENSION];
  for(i=0; i<DIMENSION; i++) {
    lo[i] = ca[i];
    hi[i] = ca[i] + rla[i];
  }

  if(m_mpiio) {
    assert(sizeof(ID_T) == 8);

    vector<POSVEL_T> *xx = new vector<POSVEL_T>;
    vector<POSVEL_T> *vx = new vector<POSVEL_T>;
    vector<POSVEL_T> *yy = new vector<POSVEL_T>;
    vector<POSVEL_T> *vy = new vector<POSVEL_T>;
    vector<POSVEL_T> *zz = new vector<POSVEL_T>;
    vector<POSVEL_T> *vz = new vector<POSVEL_T>;
    vector<POSVEL_T> *phi = new vector<POSVEL_T>;
    vector<MASK_T> *mask = new vector<MASK_T>;
    vector<ID_T> *id = new vector<ID_T>;

    copyAliveIntoVectors(xx, yy, zz, vx, vy, vz,
			 phi, id, mask, anow, 
			 Npla, 0);

    RestartIO *wr = new RestartIO(IO_WRITE_RESTART, (char *)outName, MPI_COMM_WORLD);
    wr->WriteRestart(Npla, 
		     (POSVEL_T *)&(xx[0]),
		     (POSVEL_T *)&(yy[0]),
		     (POSVEL_T *)&(zz[0]),
		     (POSVEL_T *)&(vx[0]),
		     (POSVEL_T *)&(vy[0]),
		     (POSVEL_T *)&(vz[0]),
		     (POSVEL_T *)&(phi[0]),
		     (ID_T *)&(id[0]),
		     (MASK_T *)&(mask[0]));
    delete wr;

    delete xx;
    delete vx;
    delete yy;
    delete vy;
    delete zz;
    delete vz;
    delete phi;
    delete mask;
    delete id;
  } else {

#ifdef SINGLE_RANK_OUTPUT
    if(Partition::getMyProc() == 0) {
#endif
      
      FILE *outFile = fopen(outName, "wb");
      for (i = 0; i < m_Np_local_total; i++) {
	if( (m_xArr[i]>=lo[0])*(m_xArr[i]<hi[0])*
	    (m_yArr[i]>=lo[1])*(m_yArr[i]<hi[1])*
	    (m_zArr[i]>=lo[2])*(m_zArr[i]<hi[2]) ) {
	  fwrite(&m_xArr[i], sizeof(POSVEL_T), 1, outFile);
	  fwrite(&m_vxArr[i], sizeof(POSVEL_T), 1, outFile);
	  fwrite(&m_yArr[i], sizeof(POSVEL_T), 1, outFile);
	  fwrite(&m_vyArr[i], sizeof(POSVEL_T), 1, outFile);
	  fwrite(&m_zArr[i], sizeof(POSVEL_T), 1, outFile);
	  fwrite(&m_vzArr[i], sizeof(POSVEL_T), 1, outFile);
	  fwrite(&m_phiArr[i], sizeof(POSVEL_T), 1, outFile);
	  fwrite(&m_idArr[i], sizeof(ID_T), 1, outFile);
	}
      }
      fclose(outFile); 
      
#ifdef SINGLE_RANK_OUTPUT
    }
#endif

  }

  coords_global2local(anow);

  return;
}



void Particles::copyAliveIntoVectors(vector<POSVEL_T> *xx,
                                     vector<POSVEL_T> *yy,
                                     vector<POSVEL_T> *zz,
                                     vector<POSVEL_T> *vx,
                                     vector<POSVEL_T> *vy,
                                     vector<POSVEL_T> *vz,
                                     vector<POSVEL_T> *phi,
                                     vector<ID_T> *id,
                                     vector<MASK_T> *mask,
                                     float anow,
                                     int npr,
                                     int dropParticlesQ)
{
  coords_local2global(anow);

  float ca[DIMENSION];
  Domain::corner_phys_alive(ca);
  float rla[DIMENSION];
  Domain::rL_local_alive(rla);
  float lo[DIMENSION], hi[DIMENSION];
  for(int i=0; i<DIMENSION; i++) {
    lo[i] = ca[i];
    hi[i] = ca[i] + rla[i];
  }

  int Np = m_Np_local_total;
  uint8_t *aa;
  my_posix_memalign((void **)&aa, MEM_ALIGN, Np*sizeof(uint8_t));

  for(int i=0; i<Np; i++)
    aa[i] = (m_xArr[i]>=lo[0])*(m_xArr[i]<hi[0])*(m_yArr[i]>=lo[1])*(m_yArr[i]<hi[1])*(m_zArr[i]>=lo[2])*(m_zArr[i]<hi[2]);

  vector<POSVEL_T> *mass = new vector<POSVEL_T>;
  mass->reserve(Np);
  for(int i=0; i<Np; i++)
    mass->push_back(1.0);

  vector<POSVEL_T> *pv[N_POSVEL_T];
  pv[0] = xx;
  pv[1] = yy;
  pv[2] = zz;
  pv[3] = vx;
  pv[4] = vy;
  pv[5] = vz;
  pv[6] = phi;
  pv[7] = mass;

  for(int j=0; j<N_POSVEL_T; j++) {
    pv[j]->reserve(npr);
    for(int i=0; i<Np; i++)
      if(aa[i])
        pv[j]->push_back(m_pvData[j][i]);
    if(dropParticlesQ) {
      free(m_pvData[j]);
      m_pvData[j] = NULL;
    }
  }
  updatePointers();
  id->reserve(npr);
  for(int i=0; i<Np; i++)
    if(aa[i])
      id->push_back(m_idArr[i]);
  if(dropParticlesQ) {
    free(m_idArr);
    m_idArr = NULL;
  }

  mask->reserve(npr);
  for(int i=0; i<Np; i++)
    if(aa[i])
      mask->push_back(m_maskArr[i]);
  if(dropParticlesQ) {
    free(m_maskArr);
    m_maskArr = NULL;
  }

  if(dropParticlesQ)
    dropParticles();
  else
    coords_global2local(anow);

  free(aa);

  delete mass;

  return;
}



void Particles::writeRawAscii( const char *fname )
{
  FILE* out_fp;
  out_fp = fopen(fname, "w");

  int i;
  for (i = 0; i < m_Np_local_total; i++) {
    fprintf(out_fp, "%f %f %f %f %f %f %f %d\n",
	    m_xArr[i],
	    m_yArr[i],
	    m_zArr[i],
	    m_vxArr[i],
	    m_vyArr[i],
	    m_vzArr[i],
	    m_phiArr[i],
	    m_idArr[i]);
  }
  fclose(out_fp);

  return;
}



void Particles::writeRestart( const char *outName )
{
  assert(m_coords_localQ == 1);

  if(m_mpiio) {
    assert(sizeof(ID_T) == 8);
    RestartIO *wr = new RestartIO(IO_WRITE_RESTART, (char *)outName, MPI_COMM_WORLD);
    int Nplt = m_Np_local_total;
    wr->WriteRestart(Nplt, m_xArr, m_yArr, m_zArr, m_vxArr, m_vyArr, m_vzArr, m_phiArr, m_idArr, m_maskArr);
    delete wr;
  } else {
    FILE *outFile = fopen(outName, "wb");
    fwrite(&m_Np_local_total, sizeof(int), 1, outFile);
    fwrite(&m_xArr[0], sizeof(POSVEL_T), m_Np_local_total, outFile);
    fwrite(&m_vxArr[0], sizeof(POSVEL_T), m_Np_local_total, outFile);
    fwrite(&m_yArr[0], sizeof(POSVEL_T), m_Np_local_total, outFile);
    fwrite(&m_vyArr[0], sizeof(POSVEL_T), m_Np_local_total, outFile);
    fwrite(&m_zArr[0], sizeof(POSVEL_T), m_Np_local_total, outFile);
    fwrite(&m_vzArr[0], sizeof(POSVEL_T), m_Np_local_total, outFile);
    fwrite(&m_phiArr[0], sizeof(POSVEL_T), m_Np_local_total, outFile);
    fwrite(&m_idArr[0], sizeof(ID_T), m_Np_local_total, outFile);
    fwrite(&m_maskArr[0], sizeof(MASK_T), m_Np_local_total, outFile);
    fclose(outFile); 
  }

  return;
}



void Particles::readRestart( const char *inName ) {

  if(m_mpiio) {
    assert(sizeof(ID_T) == 8);

    RestartIO *rr = new RestartIO(IO_READ_RESTART, (char *)inName, MPI_COMM_WORLD);
    m_Np_local_total = rr->ReadRestart(m_xArr, m_yArr, m_zArr, m_vxArr, m_vyArr, m_vzArr, m_phiArr, m_idArr, m_maskArr);
    delete rr;

    my_posix_memalign( (void **)&m_massArr, MEM_ALIGN, m_Np_local_total*sizeof(POSVEL_T) );
    updatePointers2();
  } else {
    FILE *inFile = fopen(inName, "rb");
    int Np;
    fread(&Np, sizeof(int), 1, inFile);
    allocParticles(Np);
    fread(&m_xArr[0], sizeof(POSVEL_T), Np, inFile);
    fread(&m_vxArr[0], sizeof(POSVEL_T), Np, inFile);
    fread(&m_yArr[0], sizeof(POSVEL_T), Np, inFile);
    fread(&m_vyArr[0], sizeof(POSVEL_T), Np, inFile);
    fread(&m_zArr[0], sizeof(POSVEL_T), Np, inFile);
    fread(&m_vzArr[0], sizeof(POSVEL_T), Np, inFile);
    fread(&m_phiArr[0], sizeof(POSVEL_T), Np, inFile);
    fread(&m_idArr[0], sizeof(ID_T), Np, inFile);
    fread(&m_maskArr[0], sizeof(MASK_T), Np, inFile);  
    fclose(inFile);
  }

  for(int i=0; i<m_Np_local_total; i++)
    m_massArr[i] = 1.0;

  m_coords_localQ = 1;

  return;
}



vector<unsigned int>* Particles::aliveIndices() {
  assert(m_coords_localQ == 1);

  vector<unsigned int> *ai = new vector<unsigned int>;
  ai->reserve(Np_local_alive());

  int ngla[DIMENSION];
  POSVEL_T lo, hi[DIMENSION];

  lo = 1.0*Domain::ng_overload();
  Domain::ng_local_alive(ngla);
  for(int i=0; i<DIMENSION; i++)
    hi[i] = lo + 1.0*ngla[i];

  for(int i=0; i<m_Np_local_total; i++)
    if( (m_xArr[i]>=lo)*(m_xArr[i]<hi[0])*(m_yArr[i]>=lo)*(m_yArr[i]<hi[1])*(m_zArr[i]>=lo)*(m_zArr[i]<hi[2]) )
      ai->push_back(i);

  return ai;
}



/*
void Particles::writeRawBin( const char *outName )
{
  STATUS_T status;
  status=0;

  FILE *outFile = fopen(outName, "wb");
  int i;
  for (i = 0; i < m_Np_local_total; i++) {
    fwrite(&m_xArr[i], sizeof(FLOAT_T), 1, outFile);
    fwrite(&m_vxArr[i], sizeof(FLOAT_T), 1, outFile);
    fwrite(&m_yArr[i], sizeof(FLOAT_T), 1, outFile);
    fwrite(&m_vyArr[i], sizeof(FLOAT_T), 1, outFile);
    fwrite(&m_zArr[i], sizeof(FLOAT_T), 1, outFile);
    fwrite(&m_vzArr[i], sizeof(FLOAT_T), 1, outFile);
    fwrite(&status, sizeof(STATUS_T), 1, outFile);
    fwrite(&m_idArr[i], sizeof(INDEX_T), 1, outFile);
  }
  fclose(outFile); 

  return;
}
*/



void Particles::subCycle(TimeStepper *gts) {
  // SimpleTimings::TimerRef t_cm = SimpleTimings::getTimer("cm");
  // SimpleTimings::TimerRef t_map1 = SimpleTimings::getTimer("map1");

  double stepFraction = 1.0/m_nsub;

  for(int step=0; step < m_nsub; step++) {
    //half stream
    // SimpleTimings::startTimer(t_map1);
    if(!m_skipStreamQ)
      map1(gts->pp(), stepFraction*gts->tau2(), gts->adot());
    // SimpleTimings::stopTimerStats(t_map1);

    //kick
    if(!m_skipKickSRQ)
      map2(gts, stepFraction);

    //half stream
    // SimpleTimings::startTimer(t_map1);
    if(!m_skipStreamQ)
      map1(gts->pp(), stepFraction*gts->tau2(), gts->adot());
    // SimpleTimings::stopTimerStats(t_map1);
  }

  return;
}



void Particles::map2(TimeStepper *ts, TS_FLOAT stepFraction) {
  int Np = m_Np_local_total;

  //local subvolume dimensions in grid units
  int nglt[DIMENSION];
  Domain::ng_local_total(nglt);

  //limits for BH tree build
  float ngltree[DIMENSION];
  ngltree[0] = 1.0*MAX( MAX( nglt[0], nglt[1] ), nglt[2] );
  ngltree[2] = ngltree[1] = ngltree[0];
  float zero[DIMENSION] = {0.0, 0.0, 0.0};

  //limits for particle velocity updates
  POSVEL_T xlo, ylo, zlo, xhi, yhi, zhi;
  xlo = m_edge;
  ylo = m_edge;
  zlo = m_edge;
  xhi = 1.0*nglt[0] - m_edge;
  yhi = 1.0*nglt[1] - m_edge;
  zhi = 1.0*nglt[2] - m_edge;

  POSVEL_T lo[DIMENSION] = { xlo, ylo, zlo };
  POSVEL_T hi[DIMENSION] = { xhi, yhi, zhi };

  POSVEL_T divscal, pi, c;
  divscal = m_gpscal*m_gpscal*m_gpscal;
  pi = 4.0*atanf(1.0);
  c = divscal/4.0/pi*ts->fscal()*ts->tau()*stepFraction;



  // SimpleTimings::TimerRef t_map2_srt = SimpleTimings::getTimer("map2s");
  // SimpleTimings::TimerRef t_map2_bld = SimpleTimings::getTimer("map2b");
  // SimpleTimings::TimerRef t_map2_wlk = SimpleTimings::getTimer("map2w");

  //move out of bounds particles to end of arrays
  // SimpleTimings::startTimer(t_map2_srt);
  resortParticles();
  // SimpleTimings::stopTimerStats(t_map2_srt);


  //only use tree for in bounds particlesx
  Np = m_Np_last;

#if 0
  for (int i = 0; i < Np; ++i) {
    printf("%d: %f %f %f\n", i, m_xArr[i], m_yArr[i], m_zArr[i]);
  }
#endif

  for(int i=0; i<m_Np_local_total; i++)
    m_massArr[i] = 1.0;

  if (m_useRCBTree == 5) {
    // SimpleTimings::startTimer(t_map2_bld);
    RCOMonopoleForceTree *sft = new RCOMonopoleForceTree(zero,
  				      ngltree,
                                      lo, hi,
  				      Np,
  				      m_xArr,
  				      m_yArr,
  				      m_zArr,
  				      m_vxArr,
  				      m_vyArr,
  				      m_vzArr,
  				      m_massArr,
                                      m_phiArr,
                                      m_idArr,
                                      m_maskArr,
  				      1.0,
                                      m_fsrrmax,
                                      m_openAngle,
                                      m_rcbTreePPN,
                                      m_rcbTreeExtraLevels,
  				      m_fl,
  				      c);
    // SimpleTimings::stopTimerStats(t_map2_bld);
    delete sft;
    return;
  } else if (m_useRCBTree == 4) {
    // SimpleTimings::startTimer(t_map2_bld);
    RCOQuadrupoleForceTree *sft = new RCOQuadrupoleForceTree(zero,
  				      ngltree,
                                      lo, hi,
  				      Np,
  				      m_xArr,
  				      m_yArr,
  				      m_zArr,
  				      m_vxArr,
  				      m_vyArr,
  				      m_vzArr,
  				      m_massArr,
                                      m_phiArr,
                                      m_idArr,
                                      m_maskArr,
  				      1.0,
                                      m_fsrrmax,
                                      m_openAngle,
                                      m_rcbTreePPN,
                                      m_rcbTreeExtraLevels,
  				      m_fl,
  				      c);
    // SimpleTimings::stopTimerStats(t_map2_bld);
    delete sft;
    return;
  } else if (m_useRCBTree == 2) {
    // SimpleTimings::startTimer(t_map2_bld);
    RCBMonopoleForceTree *sft = new RCBMonopoleForceTree(zero,
  				      ngltree,
                                      lo, hi,
  				      Np,
  				      m_xArr,
  				      m_yArr,
  				      m_zArr,
  				      m_vxArr,
  				      m_vyArr,
  				      m_vzArr,
  				      m_massArr,
                                      m_phiArr,
                                      m_idArr,
                                      m_maskArr,
  				      1.0,
                                      m_fsrrmax,
                                      m_rsm,
                                      m_openAngle,
                                      m_rcbTreePPN,
                                      m_rcbTreeExtraLevels,
                                      m_rcbTreeTaskPartMin,
  				      m_fl,
  				      c);
    // SimpleTimings::stopTimerStats(t_map2_bld);
    delete sft;
    return;
  } else if (m_useRCBTree) {
    // SimpleTimings::startTimer(t_map2_bld);
    RCBQuadrupoleForceTree *sft = new RCBQuadrupoleForceTree(zero,
  				      ngltree,
                                      lo, hi,
  				      Np,
  				      m_xArr,
  				      m_yArr,
  				      m_zArr,
  				      m_vxArr,
  				      m_vyArr,
  				      m_vzArr,
  				      m_massArr,
                                      m_phiArr,
                                      m_idArr,
                                      m_maskArr,
  				      1.0,
                                      m_fsrrmax,
                                      m_rsm,
                                      m_openAngle,
                                      m_rcbTreePPN,
                                      m_rcbTreeExtraLevels,
                                      m_rcbTreeTaskPartMin,
  				      m_fl,
  				      c);
    // SimpleTimings::stopTimerStats(t_map2_bld);
    delete sft;
    return;
  }

  //build tree
  // SimpleTimings::startTimer(t_map2_bld);
  BHForceTree *bhft = new BHForceTree(zero,
				      ngltree,
				      Np,
				      m_xArr,
				      m_yArr,
				      m_zArr,
				      m_vxArr,
				      m_vyArr,
				      m_vzArr,
				      m_massArr,
				      1.0,
				      m_fl,
				      c);
  // SimpleTimings::stopTimerStats(t_map2_bld);

  //loop over particles and gadget tree walk

  int numthreads=1, threadid=0;

#ifdef _OPENMP
  numthreads = omp_get_max_threads();
#endif

  vector<POSVEL_T>** xInteractArr = new vector<POSVEL_T>*[numthreads];
  vector<POSVEL_T>** yInteractArr = new vector<POSVEL_T>*[numthreads];
  vector<POSVEL_T>** zInteractArr = new vector<POSVEL_T>*[numthreads];
  vector<POSVEL_T>** mInteractArr = new vector<POSVEL_T>*[numthreads];

  double* tblArr = new double[numthreads];
  double* telArr = new double[numthreads];
  int* threadused = new int[numthreads];

  for(int i=0; i<numthreads; i++) {
    xInteractArr[i] = new vector<POSVEL_T>;
    yInteractArr[i] = new vector<POSVEL_T>;
    zInteractArr[i] = new vector<POSVEL_T>;
    mInteractArr[i] = new vector<POSVEL_T>;

    tblArr[i] = 0.0;
    telArr[i] = 0.0;
    threadused[i] = 0;
  }

  // SimpleTimings::startTimer(t_map2_wlk);
#ifdef _OPENMP
#pragma omp parallel for private(threadid)
#endif
  for(int i=0; i<Np; i++) {
#ifdef _OPENMP
    threadid = omp_get_thread_num();
#endif
    threadused[threadid] = 1;

    POSVEL_T updateq = 1.0;
    updateq *= (m_xArr[i] > xlo)*(m_xArr[i] < xhi);
    updateq *= (m_yArr[i] > ylo)*(m_yArr[i] < yhi);
    updateq *= (m_zArr[i] > zlo)*(m_zArr[i] < zhi);

    double tbl = 0.0, tel = 0.0;

    if(updateq) {
      if(m_useFastTreeEval)
	bhft->treeForceGadgetTopDownFast2(i, m_openAngle, m_fsrrmax,
					  xInteractArr[threadid],
					  yInteractArr[threadid],
					  zInteractArr[threadid],
					  mInteractArr[threadid],
					  &tbl, &tel);
      else
	bhft->treeForceGadgetTopDown(i, m_openAngle, m_fsrrmax);
    }

    tblArr[threadid] += tbl;
    telArr[threadid] += tel;
  }

  // SimpleTimings::stopTimerStats(t_map2_wlk);

  double tblmax, tblmin, tblavg, tbldiff;
  double telmax, telmin, telavg, teldiff;
  int numthreadsused = 0;

  tblmax = tblmin = tblArr[0];
  tblavg = tbldiff = 0.0;

  telmax = telmin = telArr[0];
  telavg = teldiff = 0.0;

  for(int i=0; i<numthreads; i++) {
    if(threadused[i]!=0) {
      numthreadsused++;
      tblmax = tblArr[i] > tblmax ? tblArr[i] : tblmax;
      tblmin = tblArr[i] < tblmin ? tblArr[i] : tblmin;
      tblavg += tblArr[i];
      telmax = telArr[i] > telmax ? telArr[i] : telmax;
      telmin = telArr[i] < telmin ? telArr[i] : telmin;
      telavg += telArr[i];
    }
  }
  tblavg *= 1.0/numthreadsused;
  telavg *= 1.0/numthreadsused;
  tbldiff = tblmax-tblmin;
  teldiff = telmax-telmin;

  // SimpleTimings::TimerRef t_blmax = SimpleTimings::getTimer("blmax");
  // SimpleTimings::TimerRef t_blmin = SimpleTimings::getTimer("blmin");
  // SimpleTimings::TimerRef t_blavg = SimpleTimings::getTimer("blavg");
  // SimpleTimings::TimerRef t_bldif = SimpleTimings::getTimer("bldif");
  // SimpleTimings::TimerRef t_elmax = SimpleTimings::getTimer("elmax");
  // SimpleTimings::TimerRef t_elmin = SimpleTimings::getTimer("elmin");
  // SimpleTimings::TimerRef t_elavg = SimpleTimings::getTimer("elavg");
  // SimpleTimings::TimerRef t_eldif = SimpleTimings::getTimer("eldif");

  // SimpleTimings::fakeTimerStats(t_blmax, tblmax);
  // SimpleTimings::fakeTimerStats(t_blmin, tblmin);
  // SimpleTimings::fakeTimerStats(t_blavg, tblavg);
  // SimpleTimings::fakeTimerStats(t_bldif, tbldiff);
  // SimpleTimings::fakeTimerStats(t_elmax, telmax);
  // SimpleTimings::fakeTimerStats(t_elmin, telmin);
  // SimpleTimings::fakeTimerStats(t_elavg, telavg);
  // SimpleTimings::fakeTimerStats(t_eldif, teldiff);

  delete bhft;
  for(int i=0; i<numthreads; i++) {
    delete xInteractArr[i];
    delete yInteractArr[i];
    delete zInteractArr[i];
    delete mInteractArr[i];
  }
  delete xInteractArr;
  delete yInteractArr;
  delete zInteractArr;
  delete mInteractArr;

  delete tblArr;
  delete telArr;
  delete threadused;

  return;
}



void Particles::subCycleCM(TimeStepper *gts) {

  return;
}



void Particles::map2CM(TimeStepper *ts, TS_FLOAT stepFraction) {
  
  return;
}



CMLite::CMLite(int ngx, int ngy, int ngz) {
  ng[0] = ngx;
  ng[1] = ngy;
  ng[2] = ngz;
  Ng = ngx*ngy*ngz;
  indxlo = new int[Ng+1];
  indxhi = new int[Ng+1];
}


CMLite::~CMLite() {
  delete [] indxlo;
  delete [] indxhi;
}
