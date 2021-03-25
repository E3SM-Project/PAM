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

struct thermo_constants {
real Rd = 287.0;
real Rv = 461.0;
real Cpd = 1004.0;
real Cpv = 1885.0;
real Cvd = 717.0; 
real Cvv = 1424.0;
real Cl = 4186.0;
real Ci = 2050.0;

real pr = 1000.0 * 100.;
real Tr = 273.15; //0 C
real Lv0 = 3.1285 * pow(10.,6);
real Lvr = Lv0 + (Cpv-Cl)*Tr; // latent heat at Tr,pr
real Lfr = 333.55 * pow(10.,6); // NOT SURE ABOUT THIS VALUE...
real gamma_d = Cpd/Cvd;
real kappa_d = Rd/Cpd;
real delta_d = Rd/Cvd;
};

class ThermoPotential {
  
public: 
  thermo_constants cst;

  virtual real YAKL_INLINE compute_U(real alpha, real entropic_var, real qd, real qv, real ql, real qi)
  {};

  virtual real YAKL_INLINE compute_dUdalpha(real alpha, real entropic_var, real qd, real qv, real ql, real qi)
  {};

  virtual real YAKL_INLINE compute_dUdentropic_var_density(real alpha, real entropic_var, real qd, real qv, real ql, real qi)
  {};

  virtual real YAKL_INLINE compute_dUdqd(real alpha, real entropic_var, real qd, real qv, real ql, real qi)
  {};

  virtual real YAKL_INLINE compute_dUdqv(real alpha, real entropic_var, real qd, real qv, real ql, real qi)
  {};

  virtual real YAKL_INLINE compute_dUdql(real alpha, real entropic_var, real qd, real qv, real ql, real qi)
  {};
  
  virtual real YAKL_INLINE compute_dUdqi(real alpha, real entropic_var, real qd, real qv, real ql, real qi)
  {};

  virtual real YAKL_INLINE compute_H(real p, real entropic_var, real qd, real qv, real ql, real qi)
  {};

  virtual real YAKL_INLINE compute_dHdp(real p, real entropic_var, real qd, real qv, real ql, real qi)
  {};

  virtual real YAKL_INLINE compute_dHdentropic_var_density(real p, real entropic_var, real qd, real qv, real ql, real qi)
  {};

  virtual real YAKL_INLINE compute_dHdqd(real p, real entropic_var, real qd, real qv, real ql, real qi)
  {};

  virtual real YAKL_INLINE compute_dHdqv(real p, real entropic_var, real qd, real qv, real ql, real qi)
  {};

  virtual real YAKL_INLINE compute_dHdql(real p, real entropic_var, real qd, real qv, real ql, real qi)
  {};
  
  virtual real YAKL_INLINE compute_dHdqi(real p, real entropic_var, real qd, real qv, real ql, real qi)
  {};
  
  virtual real YAKL_INLINE compute_alpha(real p, real T, real qd, real qv, real ql, real qi)
  {};

  virtual real YAKL_INLINE compute_entropic_var(real p, real T, real qd, real qv, real ql, real qi)
  {};
};

// This ignores any q arguments, as expected
class IdealGas_Pottemp : public ThermoPotential {
public: 

  real YAKL_INLINE compute_U(real alpha, real entropic_var, real qd, real qv, real ql, real qi)
  {
    return cst.Cvd * pow(entropic_var,cst.gamma_d) * pow(cst.Rd/(alpha*cst.pr),cst.delta_d);
  };

  real YAKL_INLINE compute_dUdalpha(real alpha, real entropic_var, real qd, real qv, real ql, real qi)
  {
    return - cst.pr * pow(entropic_var * cst.Rd/(alpha * cst.pr),cst.gamma_d);
  };

  real YAKL_INLINE compute_dUdentropic_var_density(real alpha, real entropic_var, real qd, real qv, real ql, real qi)
  {
    return cst.Cpd * pow(entropic_var * cst.Rd/(alpha * cst.pr),cst.delta_d);
  };

  real YAKL_INLINE compute_dUdqd(real alpha, real entropic_var, real qd, real qv, real ql, real qi)
  {return 0;};

  real YAKL_INLINE compute_dUdqv(real alpha, real entropic_var, real qd, real qv, real ql, real qi)
  {return 0;};

  real YAKL_INLINE compute_dUdql(real alpha, real entropic_var, real qd, real qv, real ql, real qi)
  {return 0;};
  
  real YAKL_INLINE compute_dUdqi(real alpha, real entropic_var, real qd, real qv, real ql, real qi)
  {return 0;};

//FIX
  real YAKL_INLINE compute_H(real p, real entropic_var, real qd, real qv, real ql, real qi)
  {};

//FIX
  real YAKL_INLINE compute_dHdp(real p, real entropic_var, real qd, real qv, real ql, real qi)
  {};

//FIX
  real YAKL_INLINE compute_dHdentropic_var_density(real p, real entropic_var, real qd, real qv, real ql, real qi)
  {};

  real YAKL_INLINE compute_dHdqd(real p, real entropic_var, real qd, real qv, real ql, real qi)
  {return 0;};

  real YAKL_INLINE compute_dHdqv(real p, real entropic_var, real qd, real qv, real ql, real qi)
  {return 0;};

  real YAKL_INLINE compute_dHdql(real p, real entropic_var, real qd, real qv, real ql, real qi)
  {return 0;};
  
  real YAKL_INLINE compute_dHdqi(real p, real entropic_var, real qd, real qv, real ql, real qi)
  {return 0;};
  
  real YAKL_INLINE compute_alpha(real p, real T, real qd, real qv, real ql, real qi)
  {
    return cst.Rd * T / p;
  };

  real YAKL_INLINE compute_entropic_var(real p, real T, real qd, real qv, real ql, real qi)
  {
    return T * pow(cst.pr/p, cst.kappa_d);
  };
};



// This ignores any q arguments, as expected
class IdealGas_Entropy : public ThermoPotential {
public: 

  real YAKL_INLINE compute_U(real alpha, real entropic_var, real qd, real qv, real ql, real qi)
  {};

  real YAKL_INLINE compute_dUdalpha(real alpha, real entropic_var, real qd, real qv, real ql, real qi)
  {};

  real YAKL_INLINE compute_dUdentropic_var_density(real alpha, real entropic_var, real qd, real qv, real ql, real qi)
  {};

  real YAKL_INLINE compute_dUdqd(real alpha, real entropic_var, real qd, real qv, real ql, real qi)
  {return 0;};

  real YAKL_INLINE compute_dUdqv(real alpha, real entropic_var, real qd, real qv, real ql, real qi)
  {return 0;};

  real YAKL_INLINE compute_dUdql(real alpha, real entropic_var, real qd, real qv, real ql, real qi)
  {return 0;};
  
  real YAKL_INLINE compute_dUdqi(real alpha, real entropic_var, real qd, real qv, real ql, real qi)
  {return 0;};

  real YAKL_INLINE compute_H(real p, real entropic_var, real qd, real qv, real ql, real qi)
  {};

  real YAKL_INLINE compute_dHdp(real p, real entropic_var, real qd, real qv, real ql, real qi)
  {};

  real YAKL_INLINE compute_dHdentropic_var_density(real p, real entropic_var, real qd, real qv, real ql, real qi)
  {};

  real YAKL_INLINE compute_dHdqd(real p, real entropic_var, real qd, real qv, real ql, real qi)
  {return 0;};

  real YAKL_INLINE compute_dHdqv(real p, real entropic_var, real qd, real qv, real ql, real qi)
  {return 0;};

  real YAKL_INLINE compute_dHdql(real p, real entropic_var, real qd, real qv, real ql, real qi)
  {return 0;};
  
  real YAKL_INLINE compute_dHdqi(real p, real entropic_var, real qd, real qv, real ql, real qi)
  {return 0;};

  real YAKL_INLINE compute_alpha(real p, real T, real qd, real qv, real ql, real qi)
  {};

  real YAKL_INLINE compute_entropic_var(real p, real T, real qd, real qv, real ql, real qi)
  {};

};







class ConstantKappa_VirtualPottemp : public ThermoPotential {
public: 

  real YAKL_INLINE compute_U(real alpha, real entropic_var, real qd, real qv, real ql, real qi)
  {
    return cst.Cvd * pow(entropic_var,cst.gamma_d) * pow(cst.Rd/(alpha*cst.pr),cst.delta_d) -cst.Cvd*(qd * cst.Rd + qv * cst.Rv)/cst.Rd*cst.Tr + qv*cst.Lvr - qv*cst.Rv*cst.Tr - qi*cst.Lfr;
  };

  real YAKL_INLINE compute_dUdalpha(real alpha, real entropic_var, real qd, real qv, real ql, real qi)
  {
    return - cst.pr * pow(entropic_var * cst.Rd/(alpha * cst.pr),cst.gamma_d);
  };

  real YAKL_INLINE compute_dUdentropic_var_density(real alpha, real entropic_var, real qd, real qv, real ql, real qi)
  {
    return cst.Cpd * pow(entropic_var * cst.Rd/(alpha * cst.pr),cst.delta_d);
  };

  real YAKL_INLINE compute_dUdqd(real alpha, real entropic_var, real qd, real qv, real ql, real qi)
  {return -cst.Cvd * cst.Tr;};

  real YAKL_INLINE compute_dUdqv(real alpha, real entropic_var, real qd, real qv, real ql, real qi)
  {return -cst.Cvd*cst.Rv/cst.Rd*cst.Tr + cst.Lvr - cst.Rv*cst.Tr;};

  real YAKL_INLINE compute_dUdql(real alpha, real entropic_var, real qd, real qv, real ql, real qi)
  {return 0;};
  
  real YAKL_INLINE compute_dUdqi(real alpha, real entropic_var, real qd, real qv, real ql, real qi)
  {return -cst.Lfr;};

//FIX
  real YAKL_INLINE compute_H(real p, real entropic_var, real qd, real qv, real ql, real qi)
  {return 0;};

//FIX
  real YAKL_INLINE compute_dHdp(real p, real entropic_var, real qd, real qv, real ql, real qi)
  {return 0;};

//FIX
  real YAKL_INLINE compute_dHdentropic_var_density(real p, real entropic_var, real qd, real qv, real ql, real qi)
  {return 0;};

  real YAKL_INLINE compute_dHdqd(real p, real entropic_var, real qd, real qv, real ql, real qi)
  {return -cst.Cvd * cst.Tr;};

  real YAKL_INLINE compute_dHdqv(real p, real entropic_var, real qd, real qv, real ql, real qi)
  {return -cst.Cvd*cst.Rv/cst.Rd*cst.Tr + cst.Lvr - cst.Rv*cst.Tr;};

  real YAKL_INLINE compute_dHdql(real p, real entropic_var, real qd, real qv, real ql, real qi)
  {return 0;};
  
  real YAKL_INLINE compute_dHdqi(real p, real entropic_var, real qd, real qv, real ql, real qi)
  {return -cst.Lfr;};
  
  real YAKL_INLINE compute_alpha(real p, real T, real qd, real qv, real ql, real qi)
  {
    return (qd * cst.Rd + qv * cst.Rv) * T / p;
  };

  real YAKL_INLINE compute_entropic_var(real p, real T, real qd, real qv, real ql, real qi)
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
