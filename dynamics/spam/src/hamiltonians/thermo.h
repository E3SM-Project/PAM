#ifndef _THERMO_H_
#define _THERMO_H_

#include "common.h"
#include <math.h>
// this defines thermodynamics
// it is the internal energy as a function of predicted variables (alpha, entropic variable, concentrations)
// along with it's various derivatives!
// and other useful thermodynamic relationships
// also does the same thing for enthalpy which is a function of (p, entropic variable, concentrations)
// it also includes a "solve" routine for p as a function of (rho, entropic variable, concentrations), which is what shows up in the p-variants of FC
// for all the thermodynamic potentials we consider, this solve is actually explicit since the expression can be inverted

// thermo is also responsible, in the CE/MCE equations, of computing alpha and entropic variable from p,T
// this is used in setting the initial conditions

// the currently implemented thermodynamics are
// 1) ideal gas with potential temperature as entropic variable
// 2) ideal gas with specific entropy as entropic variable
// 3) moist air with zero-volume condensates using constant kappa approximation with virtual potential temperature as entropic variable

// soon we will add
// 4) moist air with zero-volume condensates using constant kappa approximation with entropy as entropic variable
// 5) moist air with zero-volume condensates with specific entropy as entropic variable
// 6) moist air with zero-volume condensates with potential temperature as entropic variable

// In all cases of moist air we do not assume any equilibrium between the phases of water

// Some future possibilities are: non-ideal gas ie temperature dependent heat capaicty? condensates with volume? equilibrium between water phases? others?


// These are basic thermodynamic constants used in all of the thermodynamics below
// They can be overwritten when setting certain initial conditions, if a test case specifies different values
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
  
  virtual real YAKL_INLINE solve_p(real rho, real entropic_var, real qd, real qv, real ql, real qi)
  {};

  virtual real YAKL_INLINE compute_T(real alpha, real entropic_var, real qd, real qv, real ql, real qi)
  {};

  virtual real YAKL_INLINE compute_entropic_var_from_T(real alpha, real T, real qd, real qv, real ql, real qi)
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

  real YAKL_INLINE compute_H(real p, real entropic_var, real qd, real qv, real ql, real qi)
  { 
    return cst.Cpd * entropic_var * pow(p/cst.pr,cst.kappa_d);
  };

  real YAKL_INLINE compute_dHdp(real p, real entropic_var, real qd, real qv, real ql, real qi)
  { 
    return cst.Rd * entropic_var / p * pow(p/cst.pr,cst.kappa_d);
  };

  real YAKL_INLINE compute_dHdentropic_var_density(real p, real entropic_var, real qd, real qv, real ql, real qi)
  { 
    return cst.Cpd * pow(p/cst.pr,cst.kappa_d);
  };

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
  
  real YAKL_INLINE solve_p(real rho, real entropic_var, real qd, real qv, real ql, real qi)
  { 
    return cst.pr * pow(entropic_var * rho * cst.Rd/cst.pr,cst.gamma_d);
  };

  real YAKL_INLINE compute_T(real alpha, real entropic_var, real qd, real qv, real ql, real qi)
  {
    real p = cst.pr * pow(entropic_var * cst.Rd / (alpha * cst.pr),cst.gamma_d);
    return alpha * p / cst.Rd ;
  };

  real YAKL_INLINE compute_entropic_var_from_T(real alpha, real T, real qd, real qv, real ql, real qi)
  {
    real p = cst.Rd  * T / alpha;
    return cst.Rd  * T / cst.Rd * pow(cst.pr / p, cst.kappa_d);
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

  real YAKL_INLINE solve_p(real rho, real entropic_var, real qd, real qv, real ql, real qi)
  {};
};






// THE INTERNAL ENERGY HERE IS NOT BEING CONSERVED PROPERLY, IN THE SENSE THAT IT'S VALUE OSCILLATIONS STRONGLY...
// POSSIBLY DUE TO LARGE VARIATIONS INDUCED BY -cst.Cvd*(qd * cst.Rd + qv * cst.Rv)/cst.Rd*cst.Tr?
// not entirely sure here, need to done more checking!

// ALSO, FIX LATENT HEAT DEFS TO MATCH CHANGES IN PAPER

class ConstantKappa_VirtualPottemp : public ThermoPotential {
public: 

  real YAKL_INLINE compute_U(real alpha, real entropic_var, real qd, real qv, real ql, real qi)
  {
    return cst.Cvd * pow(entropic_var,cst.gamma_d) * pow(cst.Rd/(alpha*cst.pr),cst.delta_d);
    //-cst.Cvd*(qd * cst.Rd + qv * cst.Rv)/cst.Rd*cst.Tr + qv*cst.Lvr - qv*cst.Rv*cst.Tr - qi*cst.Lfr;
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
  {
    return 0;
    //return -cst.Cvd * cst.Tr;
  };

  real YAKL_INLINE compute_dUdqv(real alpha, real entropic_var, real qd, real qv, real ql, real qi)
  {
    return 0;
    //return -cst.Cvd*cst.Rv/cst.Rd*cst.Tr + cst.Lvr - cst.Rv*cst.Tr;
  };

  real YAKL_INLINE compute_dUdql(real alpha, real entropic_var, real qd, real qv, real ql, real qi)
  {return 0;};
  
  real YAKL_INLINE compute_dUdqi(real alpha, real entropic_var, real qd, real qv, real ql, real qi)
  {
    return 0;
    //return -cst.Lfr;
  };

  real YAKL_INLINE compute_H(real p, real entropic_var, real qd, real qv, real ql, real qi)
  { 
    return cst.Cpd * entropic_var * pow(p/cst.pr,cst.kappa_d) -cst.Cpd*(qd * cst.Rd + qv * cst.Rv)/cst.Rd*cst.Tr + qv*cst.Lvr + qd*cst.Rd*cst.Tr - qi*cst.Lfr;
  };

  real YAKL_INLINE compute_dHdp(real p, real entropic_var, real qd, real qv, real ql, real qi)
  { 
    return cst.Rd * entropic_var / p * pow(p/cst.pr,cst.kappa_d);
  };

  real YAKL_INLINE compute_dHdentropic_var_density(real p, real entropic_var, real qd, real qv, real ql, real qi)
  { 
    return cst.Cpd * pow(p/cst.pr,cst.kappa_d);
  };

  real YAKL_INLINE compute_dHdqd(real p, real entropic_var, real qd, real qv, real ql, real qi)
  {
    return -cst.Cvd * cst.Tr;
  };

  real YAKL_INLINE compute_dHdqv(real p, real entropic_var, real qd, real qv, real ql, real qi)
  {
    return -cst.Cvd*cst.Rv/cst.Rd*cst.Tr + cst.Lvr - cst.Rv*cst.Tr;
  };

  real YAKL_INLINE compute_dHdql(real p, real entropic_var, real qd, real qv, real ql, real qi)
  {return 0;};
  
  real YAKL_INLINE compute_dHdqi(real p, real entropic_var, real qd, real qv, real ql, real qi)
  {
    return -cst.Lfr;
  };
  
  real YAKL_INLINE compute_alpha(real p, real T, real qd, real qv, real ql, real qi)
  {
    return (qd * cst.Rd + qv * cst.Rv) * T / p;
  };

  real YAKL_INLINE compute_entropic_var(real p, real T, real qd, real qv, real ql, real qi)
  {
    return (qd * cst.Rd + qv * cst.Rv) * T / cst.Rd * pow(cst.pr/p, cst.kappa_d);
  };
  
  real YAKL_INLINE solve_p(real rho, real entropic_var, real qd, real qv, real ql, real qi)
  { 
    return cst.pr * pow(entropic_var * rho * cst.Rd/cst.pr,cst.gamma_d);
  };

  real YAKL_INLINE compute_T(real alpha, real entropic_var, real qd, real qv, real ql, real qi)
  {
    real Rstar = cst.Rd * qd + cst.Rv * qv;
    real p = cst.pr * pow(entropic_var * cst.Rd / (alpha * cst.pr),cst.gamma_d);
    return alpha * p / Rstar;
  };

  real YAKL_INLINE compute_entropic_var_from_T(real alpha, real T, real qd, real qv, real ql, real qi)
  {
    real Rstar = cst.Rd * qd + cst.Rv * qv;
    real p = Rstar * T / alpha;
    return Rstar * T / cst.Rd * pow(cst.pr / p, cst.kappa_d);
  };
};



class ConstantKappa_Entropy : public ThermoPotential {
public: 

  real YAKL_INLINE compute_U(real alpha, real entropic_var, real qd, real qv, real ql, real qi)
  {};

  real YAKL_INLINE compute_dUdalpha(real alpha, real entropic_var, real qd, real qv, real ql, real qi)
  {};

  real YAKL_INLINE compute_dUdentropic_var_density(real alpha, real entropic_var, real qd, real qv, real ql, real qi)
  {};

  real YAKL_INLINE compute_dUdqd(real alpha, real entropic_var, real qd, real qv, real ql, real qi)
  {};

  real YAKL_INLINE compute_dUdqv(real alpha, real entropic_var, real qd, real qv, real ql, real qi)
  {};

  real YAKL_INLINE compute_dUdql(real alpha, real entropic_var, real qd, real qv, real ql, real qi)
  {};
  
  real YAKL_INLINE compute_dUdqi(real alpha, real entropic_var, real qd, real qv, real ql, real qi)
  {};

  real YAKL_INLINE compute_H(real p, real entropic_var, real qd, real qv, real ql, real qi)
  {};

  real YAKL_INLINE compute_dHdp(real p, real entropic_var, real qd, real qv, real ql, real qi)
  {};

  real YAKL_INLINE compute_dHdentropic_var_density(real p, real entropic_var, real qd, real qv, real ql, real qi)
  {};

  real YAKL_INLINE compute_dHdqd(real p, real entropic_var, real qd, real qv, real ql, real qi)
  {};

  real YAKL_INLINE compute_dHdqv(real p, real entropic_var, real qd, real qv, real ql, real qi)
  {};

  real YAKL_INLINE compute_dHdql(real p, real entropic_var, real qd, real qv, real ql, real qi)
  {};
  
  real YAKL_INLINE compute_dHdqi(real p, real entropic_var, real qd, real qv, real ql, real qi)
  {};

  real YAKL_INLINE compute_alpha(real p, real T, real qd, real qv, real ql, real qi)
  {};

  real YAKL_INLINE compute_entropic_var(real p, real T, real qd, real qv, real ql, real qi)
  {};

  real YAKL_INLINE solve_p(real rho, real entropic_var, real qd, real qv, real ql, real qi)
  {};
};


class Unapprox_Pottemp : public ThermoPotential {
public: 

  real YAKL_INLINE compute_U(real alpha, real entropic_var, real qd, real qv, real ql, real qi)
  {};

  real YAKL_INLINE compute_dUdalpha(real alpha, real entropic_var, real qd, real qv, real ql, real qi)
  {};

  real YAKL_INLINE compute_dUdentropic_var_density(real alpha, real entropic_var, real qd, real qv, real ql, real qi)
  {};

  real YAKL_INLINE compute_dUdqd(real alpha, real entropic_var, real qd, real qv, real ql, real qi)
  {};

  real YAKL_INLINE compute_dUdqv(real alpha, real entropic_var, real qd, real qv, real ql, real qi)
  {};

  real YAKL_INLINE compute_dUdql(real alpha, real entropic_var, real qd, real qv, real ql, real qi)
  {};
  
  real YAKL_INLINE compute_dUdqi(real alpha, real entropic_var, real qd, real qv, real ql, real qi)
  {};

  real YAKL_INLINE compute_H(real p, real entropic_var, real qd, real qv, real ql, real qi)
  {};

  real YAKL_INLINE compute_dHdp(real p, real entropic_var, real qd, real qv, real ql, real qi)
  {};

  real YAKL_INLINE compute_dHdentropic_var_density(real p, real entropic_var, real qd, real qv, real ql, real qi)
  {};

  real YAKL_INLINE compute_dHdqd(real p, real entropic_var, real qd, real qv, real ql, real qi)
  {};

  real YAKL_INLINE compute_dHdqv(real p, real entropic_var, real qd, real qv, real ql, real qi)
  {};

  real YAKL_INLINE compute_dHdql(real p, real entropic_var, real qd, real qv, real ql, real qi)
  {};
  
  real YAKL_INLINE compute_dHdqi(real p, real entropic_var, real qd, real qv, real ql, real qi)
  {};

  real YAKL_INLINE compute_alpha(real p, real T, real qd, real qv, real ql, real qi)
  {};

  real YAKL_INLINE compute_entropic_var(real p, real T, real qd, real qv, real ql, real qi)
  {};

  real YAKL_INLINE solve_p(real rho, real entropic_var, real qd, real qv, real ql, real qi)
  {};
};

class Unapprox_Entropy : public ThermoPotential {
public: 

  real YAKL_INLINE compute_U(real alpha, real entropic_var, real qd, real qv, real ql, real qi)
  {};

  real YAKL_INLINE compute_dUdalpha(real alpha, real entropic_var, real qd, real qv, real ql, real qi)
  {};

  real YAKL_INLINE compute_dUdentropic_var_density(real alpha, real entropic_var, real qd, real qv, real ql, real qi)
  {};

  real YAKL_INLINE compute_dUdqd(real alpha, real entropic_var, real qd, real qv, real ql, real qi)
  {};

  real YAKL_INLINE compute_dUdqv(real alpha, real entropic_var, real qd, real qv, real ql, real qi)
  {};

  real YAKL_INLINE compute_dUdql(real alpha, real entropic_var, real qd, real qv, real ql, real qi)
  {};
  
  real YAKL_INLINE compute_dUdqi(real alpha, real entropic_var, real qd, real qv, real ql, real qi)
  {};

  real YAKL_INLINE compute_H(real p, real entropic_var, real qd, real qv, real ql, real qi)
  {};

  real YAKL_INLINE compute_dHdp(real p, real entropic_var, real qd, real qv, real ql, real qi)
  {};

  real YAKL_INLINE compute_dHdentropic_var_density(real p, real entropic_var, real qd, real qv, real ql, real qi)
  {};

  real YAKL_INLINE compute_dHdqd(real p, real entropic_var, real qd, real qv, real ql, real qi)
  {};

  real YAKL_INLINE compute_dHdqv(real p, real entropic_var, real qd, real qv, real ql, real qi)
  {};

  real YAKL_INLINE compute_dHdql(real p, real entropic_var, real qd, real qv, real ql, real qi)
  {};
  
  real YAKL_INLINE compute_dHdqi(real p, real entropic_var, real qd, real qv, real ql, real qi)
  {};

  real YAKL_INLINE compute_alpha(real p, real T, real qd, real qv, real ql, real qi)
  {};

  real YAKL_INLINE compute_entropic_var(real p, real T, real qd, real qv, real ql, real qi)
  {};

  real YAKL_INLINE solve_p(real rho, real entropic_var, real qd, real qv, real ql, real qi)
  {};
};

#endif
