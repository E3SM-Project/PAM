#ifndef _THERMO_H_
#define _THERMO_H_

#include "common.h"
#include <math.h>
// this defines thermodynamics
// it is the internal energy as a function of predicted variables (alpha, entropic variable, concentrations)
// along with it's various derivatives!
// and other useful thermodynamic relationships

// also do other thermodynamic potentials (key is enthalpy) here for p variants

// for entropic variable look at entropy, potential temperature, virtual potential temperature and potential enthalpy
// with either constant kappa or "unapprox"
// and play around with/without ice phase 


// HOW DO WE HANDLE SETTING INITIAL CONDITIONS
// A GIVEN TEST CASE RETURNS v(usually 0),p,T,qn as a function of position
// Then Hamiltonian or thermo gives alpha, theta as functions of p,T,qn?
// And then Hamiltonian can get rho_n,Theta from these?
// Won't be exactly hydrostatically balanced, but that is ok
// Our scheme doesn't support a state of exact discrete hydrostatic balance anyways...

//FIX THESE!
struct thermo_constants {
real Rd = 231.1;
real Cpd = 231.1;
real Cvd = 231.1;
real Rv = 231.1;
real Cpv = 231.1;
real Cvv = 231.1;
real Cl = 231.1;
real Ci = 231.1;
real pr = 231.1;
real Lvr = 231.1;
real Lfr = 231.1;
real Tr = 231.1;
real gamma_d = Cpd/Cvd;
real kappa_d = Rd/Cpd;
real delta_d = Rd/Cvd;
};
thermo_constants cst;

class ThermoPotential {
  
public: 
  virtual real YAKL_INLINE compute_U(real alpha, real theta, real qd, real qv, real ql, real qi)
  {};

  virtual real YAKL_INLINE compute_dUdalpha(real alpha, real theta, real qd, real qv, real ql, real qi)
  {};

  virtual real YAKL_INLINE compute_dUdtheta(real alpha, real theta, real qd, real qv, real ql, real qi)
  {};

  virtual real YAKL_INLINE compute_dUdqd(real alpha, real theta, real qd, real qv, real ql, real qi)
  {};

  virtual real YAKL_INLINE compute_dUdqv(real alpha, real theta, real qd, real qv, real ql, real qi)
  {};

  virtual real YAKL_INLINE compute_dUdql(real alpha, real theta, real qd, real qv, real ql, real qi)
  {};
  
  virtual real YAKL_INLINE compute_dUdqi(real alpha, real theta, real qd, real qv, real ql, real qi)
  {};
  
  virtual real YAKL_INLINE compute_alpha(real p, real T, real qd, real qv, real ql, real qi)
  {};

  virtual real YAKL_INLINE compute_theta(real p, real T, real qd, real qv, real ql, real qi)
  {};
};

// This ignores any q arguments, as expected
class IdealGas_Pottemp : public ThermoPotential {
public: 

  real YAKL_INLINE compute_U(real alpha, real theta, real qd, real qv, real ql, real qi)
  {
    return cst.Cvd * pow(theta,cst.gamma_d) * pow(cst.Rd/(alpha*cst.pr),cst.delta_d);
  };

  real YAKL_INLINE compute_dUdalpha(real alpha, real theta, real qd, real qv, real ql, real qi)
  {
    return - cst.pr * pow(theta * cst.Rd/(alpha * cst.pr),cst.gamma_d);
  };

  real YAKL_INLINE compute_dUdtheta(real alpha, real theta, real qd, real qv, real ql, real qi)
  {
    return cst.Cpd * pow(theta * cst.Rd/(alpha * cst.pr),cst.delta_d);
  };

  real YAKL_INLINE compute_dUdqd(real alpha, real theta, real qd, real qv, real ql, real qi)
  {return 0;};

  real YAKL_INLINE compute_dUdqv(real alpha, real theta, real qd, real qv, real ql, real qi)
  {return 0;};

  real YAKL_INLINE compute_dUdql(real alpha, real theta, real qd, real qv, real ql, real qi)
  {return 0;};
  
  real YAKL_INLINE compute_dUdqi(real alpha, real theta, real qd, real qv, real ql, real qi)
  {return 0;};
  
  real YAKL_INLINE compute_alpha(real p, real T, real qd, real qv, real ql, real qi)
  {
    return cst.Rd * T / p;
  };

  real YAKL_INLINE compute_theta(real p, real T, real qd, real qv, real ql, real qi)
  {
    return T * pow(cst.pr/p, cst.kappa_d);
  };
};



// FINISH THESE
// This ignores any q arguments, as expected
class IdealGas_Entropy : public ThermoPotential {
public: 

  real YAKL_INLINE compute_U(real alpha, real theta, real qd, real qv, real ql, real qi)
  {};

  real YAKL_INLINE compute_dUdalpha(real alpha, real theta, real qd, real qv, real ql, real qi)
  {};

  real YAKL_INLINE compute_dUdtheta(real alpha, real theta, real qd, real qv, real ql, real qi)
  {};

  real YAKL_INLINE compute_dUdqd(real alpha, real theta, real qd, real qv, real ql, real qi)
  {};

  real YAKL_INLINE compute_dUdqv(real alpha, real theta, real qd, real qv, real ql, real qi)
  {};

  real YAKL_INLINE compute_dUdql(real alpha, real theta, real qd, real qv, real ql, real qi)
  {};
  
  real YAKL_INLINE compute_dUdqi(real alpha, real theta, real qd, real qv, real ql, real qi)
  {};
  
  real YAKL_INLINE compute_alpha(real p, real T, real qd, real qv, real ql, real qi)
  {};

  real YAKL_INLINE compute_theta(real p, real T, real qd, real qv, real ql, real qi)
  {};
  
};







class ConstantKappa_VirtualPottemp : public ThermoPotential {
public: 

  real YAKL_INLINE compute_U(real alpha, real theta, real qd, real qv, real ql, real qi)
  {
    return cst.Cvd * pow(theta,cst.gamma_d) * pow(cst.Rd/(alpha*cst.pr),cst.delta_d) -cst.Cvd*(qd * cst.Rd + qv * cst.Rv)/cst.Rd*cst.Tr + qv*cst.Lvr - qv*cst.Rv*cst.Tr - qi*cst.Lfr;
  };

  real YAKL_INLINE compute_dUdalpha(real alpha, real theta, real qd, real qv, real ql, real qi)
  {
    return - cst.pr * pow(theta * cst.Rd/(alpha * cst.pr),cst.gamma_d);
  };

  real YAKL_INLINE compute_dUdtheta(real alpha, real theta, real qd, real qv, real ql, real qi)
  {
    return cst.Cpd * pow(theta * cst.Rd/(alpha * cst.pr),cst.delta_d);
  };

  real YAKL_INLINE compute_dUdqd(real alpha, real theta, real qd, real qv, real ql, real qi)
  {return -cst.Cvd * cst.Tr;};

  real YAKL_INLINE compute_dUdqv(real alpha, real theta, real qd, real qv, real ql, real qi)
  {return -cst.Cvd*cst.Rv/cst.Rd*cst.Tr + cst.Lvr - cst.Rv*cst.Tr;};

  real YAKL_INLINE compute_dUdql(real alpha, real theta, real qd, real qv, real ql, real qi)
  {return 0;};
  
  real YAKL_INLINE compute_dUdqi(real alpha, real theta, real qd, real qv, real ql, real qi)
  {return -cst.Lfr;};
  
  real YAKL_INLINE compute_alpha(real p, real T, real qd, real qv, real ql, real qi)
  {
    return (qd * cst.Rd + qv * cst.Rv) * T / p;
  };

  real YAKL_INLINE compute_theta(real p, real T, real qd, real qv, real ql, real qi)
  {
    return (qd * cst.Rd + qv * cst.Rv) * T / cst.Rd * pow(cst.pr/p, cst.kappa_d);
  };
  
};

// CONSTANT KAPPA WITH ENTROPY IS AN UGLY MESS BECAUSE CHEMICAL POTENTIALS ARE HIGHLY NON-TRIVIAL
// UNAPPROX WITH POTTEMP OR ENTROPY IS AN UGLY MESS BECAUSE CHEMICAL POTENTIALS ARE HIGHLY NON-TRIVIAL 
// so skip them for now...
// maybe unapprox with pottemp is slightly more tractable? main issue is chemical potentials, otherwise it looks VERY similar to ideal gas pottemp with moist constants...
// eventually do them all for comparison purposes...can check correctness by comparing different entropic variables, should be the "SAME"

#endif
