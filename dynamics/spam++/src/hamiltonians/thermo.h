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
real constexpr Rd = 231.1;
real constexpr Cpd = 231.1;
real constexpr Cvd = 231.1;
real constexpr pr = 231.1;
real constexpr Lvr = 231.1;
real constexpr Lfr = 231.1;
real constexpr Tr = 231.1;

real constexpr gamma_d = Cpd/Cvd;
real constexpr kappa_d = Rd/Cpd;
real constexpr delta_d = Rd/Cvd;


// This ignores any q arguments, as expected
class idealGas_pottemp {
  
  real YAKL_INLINE compute_U(real alpha, real theta, real qd, real qv, real ql, real qi)
  {
    return Cvd * pow(theta,gamma_d) * pow(Rd/(alpha*pr),delta_d);
  };

  real YAKL_INLINE compute_dUdalpha(real alpha, real theta, real qd, real qv, real ql, real qi)
  {
    return - pr * pow(theta * Rd/(alpha * pr),gamma_d);
  };

  real YAKL_INLINE compute_dUdtheta(real alpha, real theta, real qd, real qv, real ql, real qi)
  {
    return Cpd * pow(theta * Rd/(alpha * pr),delta_d);
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
    return Rd * T / p;
  };

  real YAKL_INLINE compute_theta(real p, real T, real qd, real qv, real ql, real qi)
  {
    return T * pow(pr/p, kappa_d);
  };
};



// FINISH THESE
// This ignores any q arguments, as expected
class idealGas_entropy {
  
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







class constantkappa_virtualpottemp {
  real YAKL_INLINE compute_U(real alpha, real theta, real qd, real qv, real ql, real qi)
  {
    return Cvd * pow(theta,gamma_d) * pow(Rd/(alpha*pr),delta_d) -Cvd*(qd * Rd + qv * Rv)/Rd*Tr + qv*Lvr - qv*Rv*Tr - qi*Lfr;
  };

  real YAKL_INLINE compute_dUdalpha(real alpha, real theta, real qd, real qv, real ql, real qi)
  {
    return - pr * pow(theta * Rd/(alpha * pr),gamma_d);
  };

  real YAKL_INLINE compute_dUdtheta(real alpha, real theta, real qd, real qv, real ql, real qi)
  {
    return Cpd * pow(theta * Rd/(alpha * pr),delta_d);
  };

  real YAKL_INLINE compute_dUdqd(real alpha, real theta, real qd, real qv, real ql, real qi)
  {return -Cvd * Tr;};

  real YAKL_INLINE compute_dUdqv(real alpha, real theta, real qd, real qv, real ql, real qi)
  {return -Cvd*Rv/Rd*Tr +Lvr - Rv*Tr;};

  real YAKL_INLINE compute_dUdql(real alpha, real theta, real qd, real qv, real ql, real qi)
  {return 0;};
  
  real YAKL_INLINE compute_dUdqi(real alpha, real theta, real qd, real qv, real ql, real qi)
  {return -Lfr;};
  
  real YAKL_INLINE compute_alpha(real p, real T, real qd, real qv, real ql, real qi)
  {
    return (qd * Rd + qv * Rv) * T / p;
  };

  real YAKL_INLINE compute_theta(real p, real T, real qd, real qv, real ql, real qi)
  {
    return (qd * Rd + qv * Rv) * T / Rd * pow(pr/p, kappa_d);
  };
  
};

// CONSTANT KAPPA WITH ENTROPY IS AN UGLY MESS BECAUSE CHEMICAL POTENTIALS ARE HIGHLY NON-TRIVIAL
// UNAPPROX WITH POTTEMP OR ENTROPY IS AN UGLY MESS BECAUSE CHEMICAL POTENTIALS ARE HIGHLY NON-TRIVIAL 
// so skip them for now...
// maybe unapprox with pottemp is slightly more tractable? main issue is chemical potentials, otherwise it looks VERY similar to ideal gas pottemp with moist constants...
// eventually do them all for comparison purposes...can check correctness by comparing different entropic variables, should be the "SAME"

#endif
