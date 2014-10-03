#ifndef TIMESTEPPER_H
#define TIMESTEPPER_H

#include <cmath>

#define TS_FLOAT double

class TimeStepper {
 public:

  TimeStepper(TS_FLOAT alpha_, TS_FLOAT ain_, TS_FLOAT afin_,
	      int nsteps_, TS_FLOAT omegatot_);
  ~TimeStepper();

  void advanceHalfStep();
  void advanceFullStep();
  void reverseHalfStep();
  void reverseFullStep();

  TS_FLOAT aa() { return m_aa; }
  TS_FLOAT pp() { return m_pp; }
  TS_FLOAT zz() { return m_zz; }
  TS_FLOAT alpha() { return m_alpha; }
  TS_FLOAT tau() { return m_tau; }
  TS_FLOAT tau2() { return m_tau2; }
  TS_FLOAT adot() { return m_adot; }
  TS_FLOAT omegatot() { return m_omegatot; }
  TS_FLOAT ain() { return m_ain; }
  TS_FLOAT afin() { return m_afin; }
  TS_FLOAT pin() { return m_pin; }
  TS_FLOAT pfin() { return m_pfin; }
  TS_FLOAT zin() { return m_zin; }
  TS_FLOAT zfin() { return m_zfin; }
  int nsteps() { return m_nsteps; }
  TS_FLOAT phiscal() { return m_phiscal; }
  TS_FLOAT fscal() { return m_fscal; }

 private:

  TimeStepper();
  TimeStepper( const TimeStepper& );
  TimeStepper& operator = (const TimeStepper& );

  void set_adot();
  void set_scal();

  TS_FLOAT m_aa;
  TS_FLOAT m_pp;
  TS_FLOAT m_zz;
  TS_FLOAT m_alpha;
  TS_FLOAT m_tau;
  TS_FLOAT m_tau2;
  TS_FLOAT m_adot;
  TS_FLOAT m_omegatot;
  TS_FLOAT m_ain;
  TS_FLOAT m_afin;
  TS_FLOAT m_pin;
  TS_FLOAT m_pfin;
  TS_FLOAT m_zin;
  TS_FLOAT m_zfin;
  int m_nsteps;
  TS_FLOAT m_phiscal;
  TS_FLOAT m_fscal;
};

#endif
