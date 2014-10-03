#include "TimeStepper.h"

TimeStepper::TimeStepper(TS_FLOAT alpha_, TS_FLOAT ain_, TS_FLOAT afin_,
			 int nsteps_, TS_FLOAT omegatot_) :
  m_aa(-1.0),
  m_pp(-1.0),
  m_alpha(-1.0),
  m_tau(-1.0),
  m_tau2(-1.0),
  m_adot(-1.0),
  m_omegatot(-1.0),
  m_ain(-1.0),
  m_afin(-1.0),
  m_pin(-1.0),
  m_pfin(-1.0),
  m_nsteps(-1),
  m_phiscal(-1.0),
  m_fscal(-1.0)
{
  m_alpha = alpha_;
  m_ain = ain_;
  m_afin = afin_;
  m_nsteps = nsteps_;
  m_omegatot = omegatot_;

  m_pin = pow(m_ain, m_alpha);
  m_pfin = pow(m_afin, m_alpha);
  m_zin = 1.0/m_ain - 1.0;
  m_zfin = 1.0/m_afin - 1.0;
  m_tau = (m_pfin - m_pin)/(1.0*m_nsteps);
  m_tau2 = 0.5*m_tau;

  m_pp = m_pin;
  m_aa = m_ain;
  m_zz = m_zin;

  set_adot();
  set_scal();
}

TimeStepper::~TimeStepper() {

}

void TimeStepper::set_adot() {
  TS_FLOAT pp1 = pow(m_pp, 3.0/m_alpha);
  TS_FLOAT tmp = m_omegatot + (1.0-m_omegatot)*pp1;
  tmp = tmp / pow(m_pp, 1.0/m_alpha);
  m_adot = sqrt(tmp);
  return;
}

void TimeStepper::set_scal() {
  set_adot();
  m_phiscal = 1.5*m_omegatot/pow(m_pp, 1.0/m_alpha);
  m_fscal = m_phiscal/(m_alpha*m_adot*pow(m_pp, 1.0 - (1.0/m_alpha)));
  return;
}

void TimeStepper::advanceHalfStep() {
  m_pp += m_tau2;
  m_aa = pow(m_pp, 1.0/m_alpha);
  m_zz = 1.0/m_aa - 1.0;
  set_adot();
  set_scal();
  return;
}

void TimeStepper::reverseHalfStep() {
  m_pp -= m_tau2;
  m_aa = pow(m_pp, 1.0/m_alpha);
  m_zz = 1.0/m_aa - 1.0;
  set_adot();
  set_scal();
  return;
}

void TimeStepper::advanceFullStep() {
  advanceHalfStep();
  advanceHalfStep();
  return;
}

void TimeStepper::reverseFullStep() {
  reverseHalfStep();
  reverseHalfStep();
  return;
}
