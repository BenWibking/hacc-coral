#ifndef DOMAIN_H
#define DOMAIN_H

#include "Definition.h"
#include "Partition.h"
#include "Basedata.h"

using namespace std;

class Domain {

 public:

  Domain() {};
  ~Domain() {};

  static void initialize(const Basedata &);

  static float grid2phys_pos() { return m_grid2phys_pos; }
  static float phys2grid_pos() { return m_phys2grid_pos; }
  static float grid2phys_vel() { return m_grid2phys_vel; }
  static float phys2grid_vel() { return m_phys2grid_vel; }

  //scalars
  static float rL() { return m_rL; }
  static int ng() { return m_ng; }
  static float oL() { return m_oL; }
  static int ng_overload() { return m_ng_overload; }
  static int Ng_local_alive() { return m_Ng_local_alive; }
  static int Ng_local_total() { return m_Ng_local_total; }
  static int maxNpLocal() { return m_maxNpLocal; }

  //direction dependent
  static void rL_local_alive(float[]);
  static void rL_local_total(float[]);
  static void ng_local_alive(int[]);
  static void ng_local_total(int[]);
  static void corner_grid_alive(int[]);
  static void corner_grid_total(int[]);
  static void corner_phys_alive(float[]);
  static void corner_phys_total(float[]);

 private:

  static float m_grid2phys_pos;
  static float m_phys2grid_pos;
  static float m_grid2phys_vel;
  static float m_phys2grid_vel;

  //scalar
  static float m_rL;
  static int m_ng;
  static float m_oL;
  static int m_ng_overload;
  static int m_Ng_local_alive;
  static int m_Ng_local_total;
  static int m_maxNpLocal;

  //direction dependent
  static float m_rL_local_alive[DIMENSION];
  static float m_rL_local_total[DIMENSION];
  static int m_ng_local_alive[DIMENSION];
  static int m_ng_local_total[DIMENSION];
  static int m_corner_grid_alive[DIMENSION];
  static int m_corner_grid_total[DIMENSION];
  static float m_corner_phys_alive[DIMENSION];
  static float m_corner_phys_total[DIMENSION];
};



#endif
