#pragma once

#include "common.h"
#include <math.h>
// this defines thermodynamics
// it is the internal energy as a function of predicted variables (alpha,
// entropic variable, concentrations) along with it's various derivatives! and
// other useful thermodynamic relationships also does the same thing for
// enthalpy which is a function of (p, entropic variable, concentrations) it
// also includes a "solve" routine for p as a function of (rho, entropic
// variable, concentrations), which is what shows up in the p-variants of FC for
// all the thermodynamic potentials we consider, this solve is actually explicit
// since the expression can be inverted

// thermo is also responsible, in the CE/MCE equations, of computing alpha and
// entropic variable from p,T this is used in setting the initial conditions

// the currently implemented thermodynamics are
// 1) ideal gas with potential temperature as entropic variable
// 2) ideal gas with specific entropy as entropic variable
// 3) moist air with zero-volume condensates using constant kappa approximation
// with virtual potential temperature as entropic variable

// soon we will add
// 4) moist air with zero-volume condensates using constant kappa approximation
// with entropy as entropic variable 5) moist air with zero-volume condensates
// with specific entropy as entropic variable 6) moist air with zero-volume
// condensates with potential temperature as entropic variable

// In all cases of moist air we do not assume any equilibrium between the phases
// of water

// Some future possibilities are: non-ideal gas ie temperature dependent heat
// capaicty? condensates with volume? equilibrium between water phases? others?

// These are basic thermodynamic constants used in all of the thermodynamics
// below They can be overwritten when setting certain initial conditions, if a
// test case specifies different values

namespace pamc {

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
  real Tr = 273.15; // 0 C
  real Lv0 = 3.1285 * pow(10., 6.);
  real Lvr = Lv0 + (Cpv - Cl) * Tr; // latent heat at Tr,pr
  real Lfr = 333.55 * pow(10., 6.); // NOT SURE ABOUT THIS VALUE...
  real gamma_d = Cpd / Cvd;
  real kappa_d = Rd / Cpd;
  real delta_d = Rd / Cvd;
};

class ThermoNone {

public:
  static constexpr bool moist_species_decouple_from_dynamics = true;
  thermo_constants cst;
};

// This ignores any q arguments, as expected
class IdealGas_Pottemp {
public:
  static constexpr bool moist_species_decouple_from_dynamics = true;
  thermo_constants cst;

  real YAKL_INLINE compute_U(real alpha, real entropic_var, real qd, real qv,
                             real ql, real qi) const {
    return cst.Cvd * pow(entropic_var, cst.gamma_d) *
           pow(cst.Rd / (alpha * cst.pr), cst.delta_d);
  };

  real YAKL_INLINE compute_dUdalpha(real alpha, real entropic_var, real qd,
                                    real qv, real ql, real qi) const {
    return -cst.pr * pow(entropic_var * cst.Rd / (alpha * cst.pr), cst.gamma_d);
  };

  real YAKL_INLINE compute_dUdentropic_var(real alpha, real entropic_var,
                                           real qd, real qv, real ql,
                                           real qi) const {
    return cst.Cpd * pow(entropic_var * cst.Rd / (alpha * cst.pr), cst.delta_d);
  };

  real YAKL_INLINE compute_dUdqd(real alpha, real entropic_var, real qd,
                                 real qv, real ql, real qi) const {
    return 0;
  };

  real YAKL_INLINE compute_dUdqv(real alpha, real entropic_var, real qd,
                                 real qv, real ql, real qi) const {
    return 0;
  };

  real YAKL_INLINE compute_dUdql(real alpha, real entropic_var, real qd,
                                 real qv, real ql, real qi) const {
    return 0;
  };

  real YAKL_INLINE compute_dUdqi(real alpha, real entropic_var, real qd,
                                 real qv, real ql, real qi) const {
    return 0;
  };

  real YAKL_INLINE compute_H(real p, real entropic_var, real qd, real qv,
                             real ql, real qi) const {
    return cst.Cpd * entropic_var * pow(p / cst.pr, cst.kappa_d);
  };

  real YAKL_INLINE compute_dHdp(real p, real entropic_var, real qd, real qv,
                                real ql, real qi) const {
    return cst.Rd * entropic_var / p * pow(p / cst.pr, cst.kappa_d);
  };

  real YAKL_INLINE compute_dHdentropic_var(real p, real entropic_var, real qd,
                                           real qv, real ql, real qi) const {
    return cst.Cpd * pow(p / cst.pr, cst.kappa_d);
  };

  real YAKL_INLINE compute_dHdqd(real p, real entropic_var, real qd, real qv,
                                 real ql, real qi) const {
    return 0;
  };

  real YAKL_INLINE compute_dHdqv(real p, real entropic_var, real qd, real qv,
                                 real ql, real qi) const {
    return 0;
  };

  real YAKL_INLINE compute_dHdql(real p, real entropic_var, real qd, real qv,
                                 real ql, real qi) const {
    return 0;
  };

  real YAKL_INLINE compute_dHdqi(real p, real entropic_var, real qd, real qv,
                                 real ql, real qi) const {
    return 0;
  };

  real YAKL_INLINE compute_alpha(real p, real T, real qd, real qv, real ql,
                                 real qi) const {
    return cst.Rd * T / p;
  };

  real YAKL_INLINE compute_entropic_var_from_p_T(real p, real T, real qd,
                                                 real qv, real ql,
                                                 real qi) const {
    // real entropy =  cst.Cpd * log(T / cst.Tr) - cst.Rd * log(p / cst.pr);
    // return theta(entropy);
    return T * pow(cst.pr / p, cst.kappa_d);
  };

  real YAKL_INLINE solve_p(real rho, real entropic_var, real qd, real qv,
                           real ql, real qi) const {
    return cst.pr * pow(entropic_var * rho * cst.Rd / cst.pr, cst.gamma_d);
  };

  real YAKL_INLINE compute_T_from_alpha(real alpha, real entropic_var, real qd,
                                        real qv, real ql, real qi) const {
    real p =
        cst.pr * pow(entropic_var * cst.Rd / (alpha * cst.pr), cst.gamma_d);
    return alpha * p / cst.Rd;
  }

  real YAKL_INLINE compute_T_from_p(real p, real entropic_var, real qd, real qv,
                                    real ql, real qi) const {
    return pow(p / cst.pr, cst.kappa_d) * entropic_var;
  };

  real YAKL_INLINE compute_entropic_var_from_alpha_T(real alpha, real T,
                                                     real qd, real qv, real ql,
                                                     real qi) const {
    real p = cst.Rd * T / alpha;
    return cst.Rd * T / cst.Rd * pow(cst.pr / p, cst.kappa_d);
  }

  real YAKL_INLINE compute_dpdentropic_var(real alpha, real entropic_var,
                                           real qd, real qv, real ql,
                                           real qi) const {
    real rho = 1 / alpha;
    real p = solve_p(rho, entropic_var, qd, qv, ql, qi);
    return cst.gamma_d * p / entropic_var;
  }

  real YAKL_INLINE compute_soundspeed(real alpha, real entropic_var, real qd,
                                      real qv, real ql, real qi) const {

    real rho = 1 / alpha;
    real p = solve_p(rho, entropic_var, qd, qv, ql, qi);
    return sqrt(cst.gamma_d * p * alpha);
  }
};

// This ignores any q arguments, as expected
class IdealGas_Entropy {
public:
  static constexpr bool moist_species_decouple_from_dynamics = true;
  thermo_constants cst;

  real YAKL_INLINE compute_U(real alpha, real entropic_var, real qd, real qv,
                             real ql, real qi) const {
    return cst.Cvd * cst.Tr *
           pow(alpha * cst.pr / (cst.Rd * cst.Tr), -cst.delta_d) *
           exp(entropic_var / cst.Cvd);
  };

  real YAKL_INLINE compute_dUdalpha(real alpha, real entropic_var, real qd,
                                    real qv, real ql, real qi) const {

    real U = compute_U(alpha, entropic_var, qd, qv, ql, qi);
    return -cst.Rd / cst.Cvd * U / alpha;
  };

  real YAKL_INLINE compute_dUdentropic_var(real alpha, real entropic_var,
                                           real qd, real qv, real ql,
                                           real qi) const {
    real U = compute_U(alpha, entropic_var, qd, qv, ql, qi);
    return U / cst.Cvd;
  };

  real YAKL_INLINE compute_dUdqd(real alpha, real entropic_var, real qd,
                                 real qv, real ql, real qi) const {
    return 0;
  };

  real YAKL_INLINE compute_dUdqv(real alpha, real entropic_var, real qd,
                                 real qv, real ql, real qi) const {
    return 0;
  };

  real YAKL_INLINE compute_dUdql(real alpha, real entropic_var, real qd,
                                 real qv, real ql, real qi) const {
    return 0;
  };

  real YAKL_INLINE compute_dUdqi(real alpha, real entropic_var, real qd,
                                 real qv, real ql, real qi) const {
    return 0;
  };

  real YAKL_INLINE compute_H(real p, real entropic_var, real qd, real qv,
                             real ql, real qi) const {
    return cst.Cpd * cst.Tr * pow(p / cst.pr, cst.kappa_d) *
           exp(entropic_var / cst.Cpd);
  };

  real YAKL_INLINE compute_dHdp(real p, real entropic_var, real qd, real qv,
                                real ql, real qi) const {
    real H = compute_H(p, entropic_var, qd, qv, ql, qi);
    return cst.Rd * H / (p * cst.Cpd);
  }

  real YAKL_INLINE compute_dHdentropic_var(real p, real entropic_var, real qd,
                                           real qv, real ql, real qi) const {
    real H = compute_H(p, entropic_var, qd, qv, ql, qi);
    return H / cst.Cpd;
  };

  real YAKL_INLINE compute_dHdqd(real p, real entropic_var, real qd, real qv,
                                 real ql, real qi) const {
    return 0;
  };

  real YAKL_INLINE compute_dHdqv(real p, real entropic_var, real qd, real qv,
                                 real ql, real qi) const {
    return 0;
  };

  real YAKL_INLINE compute_dHdql(real p, real entropic_var, real qd, real qv,
                                 real ql, real qi) const {
    return 0;
  };

  real YAKL_INLINE compute_dHdqi(real p, real entropic_var, real qd, real qv,
                                 real ql, real qi) const {
    return 0;
  };

  real YAKL_INLINE compute_alpha(real p, real T, real qd, real qv, real ql,
                                 real qi) const {
    return cst.Rd * T / p;
  };

  real YAKL_INLINE compute_entropic_var_from_p_T(real p, real T, real qd,
                                                 real qv, real ql,
                                                 real qi) const {
    return cst.Cpd * log(T / cst.Tr) - cst.Rd * log(p / cst.pr);
  };

  real YAKL_INLINE solve_p(real rho, real entropic_var, real qd, real qv,
                           real ql, real qi) const {
    real alpha = 1 / rho;
    real U = compute_U(alpha, entropic_var, qd, qv, ql, qi);
    return cst.Rd / cst.Cvd * U / alpha;
  };

  real YAKL_INLINE compute_entropic_var_from_alpha_T(real alpha, real T,
                                                     real qd, real qv, real ql,
                                                     real qi) const {
    // real alpha_r = cst.Rd * cst.Tr / cst.pr;
    // cst.Cvd * log(T / cst.Tr) + cst.Rd * log(alpha / alpha_r);

    real p = cst.Rd * T / alpha;
    return compute_entropic_var_from_p_T(p, T, qd, qv, ql, qi);
  }

  real YAKL_INLINE compute_T_from_alpha(real alpha, real entropic_var, real qd,
                                        real qv, real ql, real qi) const {

    real U = compute_U(alpha, entropic_var, qd, qv, ql, qi);
    return U / cst.Cvd;
  }

  real YAKL_INLINE compute_T_from_p(real p, real entropic_var, real qd, real qv,
                                    real ql, real qi) const {
    return compute_dHdentropic_var(p, entropic_var, qd, qv, ql, qi);
  };

  real YAKL_INLINE compute_dpdentropic_var(real alpha, real entropic_var,
                                           real qd, real qv, real ql,
                                           real qi) const {
    real dUds = compute_dUdentropic_var(alpha, entropic_var, qd, qv, ql, qi);
    return cst.Rd / cst.Cvd * dUds / alpha;
  }

  real YAKL_INLINE compute_soundspeed(real alpha, real entropic_var, real qd,
                                      real qv, real ql, real qi) const {

    real rho = 1 / alpha;
    real p = solve_p(rho, entropic_var, qd, qv, ql, qi);
    return sqrt(cst.gamma_d * p * alpha);
  }
};

class ConstantKappa_VirtualPottemp {
public:
  static constexpr bool moist_species_decouple_from_dynamics = true;

  thermo_constants cst;
  real YAKL_INLINE compute_U(real alpha, real entropic_var, real qd, real qv,
                             real ql, real qi) const {
    const real Rstar = qd * cst.Rd + qv * cst.Rv;
    return cst.Cvd * pow(entropic_var, cst.gamma_d) *
               pow(cst.Rd / (alpha * cst.pr), cst.delta_d) -
           cst.Cvd * Rstar / cst.Rd * cst.Tr - qv * cst.Rv * cst.Tr +
           qv * (cst.Lvr + cst.Lfr) + ql * cst.Lfr;
  };

  real YAKL_INLINE compute_dUdalpha(real alpha, real entropic_var, real qd,
                                    real qv, real ql, real qi) const {
    return -cst.pr * pow(entropic_var * cst.Rd / (alpha * cst.pr), cst.gamma_d);
  };

  real YAKL_INLINE compute_dUdentropic_var(real alpha, real entropic_var,
                                           real qd, real qv, real ql,
                                           real qi) const {
    return cst.Cpd * pow(entropic_var * cst.Rd / (alpha * cst.pr), cst.delta_d);
  };

  real YAKL_INLINE compute_dUdqd(real alpha, real entropic_var, real qd,
                                 real qv, real ql, real qi) const {
    return -cst.Cvd * cst.Tr;
  };

  real YAKL_INLINE compute_dUdqv(real alpha, real entropic_var, real qd,
                                 real qv, real ql, real qi) const {
    return -cst.Cvd * cst.Rv / cst.Rd * cst.Tr + cst.Lvr + cst.Lfr -
           cst.Rv * cst.Tr;
  };

  real YAKL_INLINE compute_dUdql(real alpha, real entropic_var, real qd,
                                 real qv, real ql, real qi) const {
    return cst.Lfr;
  };

  real YAKL_INLINE compute_dUdqi(real alpha, real entropic_var, real qd,
                                 real qv, real ql, real qi) const {
    return 0;
  };

  real YAKL_INLINE compute_H(real p, real entropic_var, real qd, real qv,
                             real ql, real qi) const {
    const real Rstar = qd * cst.Rd + qv * cst.Rv;
    return cst.Cpd * entropic_var * pow(p / cst.pr, cst.kappa_d) -
           cst.Cpd * Rstar / cst.Rd * cst.Tr + qd * cst.Rd * cst.Tr +
           qv * (cst.Lvr + cst.Lfr) + ql * cst.Lfr;
  };

  real YAKL_INLINE compute_dHdp(real p, real entropic_var, real qd, real qv,
                                real ql, real qi) const {
    return cst.Rd * entropic_var / p * pow(p / cst.pr, cst.kappa_d);
  };

  real YAKL_INLINE compute_dHdentropic_var(real p, real entropic_var, real qd,
                                           real qv, real ql, real qi) const {
    return cst.Cpd * pow(p / cst.pr, cst.kappa_d);
  };

  real YAKL_INLINE compute_dHdqd(real p, real entropic_var, real qd, real qv,
                                 real ql, real qi) const {
    return -cst.Cpd * cst.Tr + cst.Rd * cst.Tr;
  };

  real YAKL_INLINE compute_dHdqv(real p, real entropic_var, real qd, real qv,
                                 real ql, real qi) const {
    return -cst.Cpd * cst.Rv / cst.Rd * cst.Tr + cst.Lvr + cst.Lfr;
  };

  real YAKL_INLINE compute_dHdql(real p, real entropic_var, real qd, real qv,
                                 real ql, real qi) const {
    return cst.Lfr;
  };

  real YAKL_INLINE compute_dHdqi(real p, real entropic_var, real qd, real qv,
                                 real ql, real qi) const {
    return 0;
  };

  real YAKL_INLINE compute_alpha(real p, real T, real qd, real qv, real ql,
                                 real qi) const {
    return (qd * cst.Rd + qv * cst.Rv) * T / p;
  };

  real YAKL_INLINE compute_entropic_var_from_p_T(real p, real T, real qd,
                                                 real qv, real ql,
                                                 real qi) const {
    return (qd * cst.Rd + qv * cst.Rv) * T / cst.Rd *
           pow(cst.pr / p, cst.kappa_d);
  };

  real YAKL_INLINE solve_p(real rho, real entropic_var, real qd, real qv,
                           real ql, real qi) const {
    return cst.pr * pow(entropic_var * rho * cst.Rd / cst.pr, cst.gamma_d);
  };

  real YAKL_INLINE compute_T_from_alpha(real alpha, real entropic_var, real qd,
                                        real qv, real ql, real qi) const {
    real Rstar = cst.Rd * qd + cst.Rv * qv;
    real p =
        cst.pr * pow(entropic_var * cst.Rd / (alpha * cst.pr), cst.gamma_d);
    return alpha * p / Rstar;
  };

  real YAKL_INLINE compute_T_from_p(real p, real entropic_var, real qd, real qv,
                                    real ql, real qi) const {
    real Rstar = cst.Rd * qd + cst.Rv * qv;
    return pow(p / cst.pr, cst.kappa_d) * entropic_var * cst.Rd / Rstar;
  };

  real YAKL_INLINE compute_entropic_var_from_alpha_T(real alpha, real T,
                                                     real qd, real qv, real ql,
                                                     real qi) const {
    real Rstar = cst.Rd * qd + cst.Rv * qv;
    real p = Rstar * T / alpha;
    return Rstar * T / cst.Rd * pow(cst.pr / p, cst.kappa_d);
  };

  real YAKL_INLINE compute_dpdentropic_var(real alpha, real entropic_var,
                                           real qd, real qv, real ql,
                                           real qi) const {
    real rho = 1 / alpha;
    real p = solve_p(rho, entropic_var, qd, qv, ql, qi);
    return cst.gamma_d * p / entropic_var;
  }

  real YAKL_INLINE compute_soundspeed(real alpha, real entropic_var, real qd,
                                      real qv, real ql, real qi) const {

    real rho = 1 / alpha;
    real p = solve_p(rho, entropic_var, qd, qv, ql, qi);
    return sqrt(cst.gamma_d * p * alpha);
  }
};

class ConstantKappa_Entropy {
public:
  static constexpr bool moist_species_decouple_from_dynamics = false;
  thermo_constants cst;
  // real YAKL_INLINE compute_U(real alpha, real entropic_var, real qd, real qv,
  //                            real ql, real qi) const {};

  // real YAKL_INLINE compute_dUdalpha(real alpha, real entropic_var, real qd,
  //                                   real qv, real ql, real qi) const {};

  // real YAKL_INLINE compute_dUdentropic_var(real alpha, real entropic_var,
  //                                          real qd, real qv, real ql,
  //                                          real qi) const {};

  // real YAKL_INLINE compute_dUdqd(real alpha, real entropic_var, real qd,
  //                                real qv, real ql, real qi) const {};

  // real YAKL_INLINE compute_dUdqv(real alpha, real entropic_var, real qd,
  //                                real qv, real ql, real qi) const {};

  // real YAKL_INLINE compute_dUdql(real alpha, real entropic_var, real qd,
  //                                real qv, real ql, real qi) const {};

  // real YAKL_INLINE compute_dUdqi(real alpha, real entropic_var, real qd,
  //                                real qv, real ql, real qi) const {};

  // real YAKL_INLINE compute_H(real p, real entropic_var, real qd, real qv,
  //                            real ql, real qi) const {};

  // real YAKL_INLINE compute_dHdp(real p, real entropic_var, real qd, real qv,
  //                               real ql, real qi) const {};

  // real YAKL_INLINE compute_dHdentropic_var(real p, real entropic_var, real
  // qd,
  //                                          real qv, real ql, real qi) const
  //                                          {};

  // real YAKL_INLINE compute_dHdqd(real p, real entropic_var, real qd, real qv,
  //                                real ql, real qi) const {};

  // real YAKL_INLINE compute_dHdqv(real p, real entropic_var, real qd, real qv,
  //                                real ql, real qi) const {};

  // real YAKL_INLINE compute_dHdql(real p, real entropic_var, real qd, real qv,
  //                                real ql, real qi) const {};

  // real YAKL_INLINE compute_dHdqi(real p, real entropic_var, real qd, real qv,
  //                                real ql, real qi) const {};

  // real YAKL_INLINE compute_alpha(real p, real T, real qd, real qv, real ql,
  //                                real qi) const {};

  // real YAKL_INLINE compute_entropic_var_from_p_T(real p, real T, real qd,
  //                                                real qv, real ql,
  //                                                real qi) const {};

  // real YAKL_INLINE solve_p(real rho, real entropic_var, real qd, real qv,
  //                          real ql, real qi) const {};
};

class Unapprox_Pottemp {
public:
  static constexpr bool moist_species_decouple_from_dynamics = false;
  thermo_constants cst;
  // real YAKL_INLINE compute_U(real alpha, real entropic_var, real qd, real qv,
  //                            real ql, real qi) const {};

  // real YAKL_INLINE compute_dUdalpha(real alpha, real entropic_var, real qd,
  //                                   real qv, real ql, real qi) const {};

  // real YAKL_INLINE compute_dUdentropic_var(real alpha, real entropic_var,
  //                                          real qd, real qv, real ql,
  //                                          real qi) const {};

  // real YAKL_INLINE compute_dUdqd(real alpha, real entropic_var, real qd,
  //                                real qv, real ql, real qi) const {};

  // real YAKL_INLINE compute_dUdqv(real alpha, real entropic_var, real qd,
  //                                real qv, real ql, real qi) const {};

  // real YAKL_INLINE compute_dUdql(real alpha, real entropic_var, real qd,
  //                                real qv, real ql, real qi) const {};

  // real YAKL_INLINE compute_dUdqi(real alpha, real entropic_var, real qd,
  //                                real qv, real ql, real qi) const {};

  // real YAKL_INLINE compute_H(real p, real entropic_var, real qd, real qv,
  //                            real ql, real qi) const {};

  // real YAKL_INLINE compute_dHdp(real p, real entropic_var, real qd, real qv,
  //                               real ql, real qi) const {};

  // real YAKL_INLINE compute_dHdentropic_var(real p, real entropic_var, real
  // qd,
  //                                          real qv, real ql, real qi) const
  //                                          {};

  // real YAKL_INLINE compute_dHdqd(real p, real entropic_var, real qd, real qv,
  //                                real ql, real qi) const {};

  // real YAKL_INLINE compute_dHdqv(real p, real entropic_var, real qd, real qv,
  //                                real ql, real qi) const {};

  // real YAKL_INLINE compute_dHdql(real p, real entropic_var, real qd, real qv,
  //                                real ql, real qi) const {};

  // real YAKL_INLINE compute_dHdqi(real p, real entropic_var, real qd, real qv,
  //                                real ql, real qi) const {};

  // real YAKL_INLINE compute_alpha(real p, real T, real qd, real qv, real ql,
  //                                real qi) const {};

  // real YAKL_INLINE compute_entropic_var_from_p_T(real p, real T, real qd,
  //                                                real qv, real ql,
  //                                                real qi) const {};

  // real YAKL_INLINE solve_p(real rho, real entropic_var, real qd, real qv,
  //                          real ql, real qi) const {};
};

class Unapprox_Entropy {
public:
  static constexpr bool moist_species_decouple_from_dynamics = false;
  thermo_constants cst;
  // real YAKL_INLINE compute_U(real alpha, real entropic_var, real qd, real qv,
  //                            real ql, real qi) const {};

  // real YAKL_INLINE compute_dUdalpha(real alpha, real entropic_var, real qd,
  //                                   real qv, real ql, real qi) const {};

  // real YAKL_INLINE compute_dUdentropic_var(real alpha, real entropic_var,
  //                                          real qd, real qv, real ql,
  //                                          real qi) const {};

  // real YAKL_INLINE compute_dUdqd(real alpha, real entropic_var, real qd,
  //                                real qv, real ql, real qi) const {};

  // real YAKL_INLINE compute_dUdqv(real alpha, real entropic_var, real qd,
  //                                real qv, real ql, real qi) const {};

  // real YAKL_INLINE compute_dUdql(real alpha, real entropic_var, real qd,
  //                                real qv, real ql, real qi) const {};

  // real YAKL_INLINE compute_dUdqi(real alpha, real entropic_var, real qd,
  //                                real qv, real ql, real qi) const {};

  // real YAKL_INLINE compute_H(real p, real entropic_var, real qd, real qv,
  //                            real ql, real qi) const {};

  // real YAKL_INLINE compute_dHdp(real p, real entropic_var, real qd, real qv,
  //                               real ql, real qi) const {};

  // real YAKL_INLINE compute_dHdentropic_var(real p, real entropic_var, real
  // qd,
  //                                          real qv, real ql, real qi) const
  //                                          {};

  // real YAKL_INLINE compute_dHdqd(real p, real entropic_var, real qd, real qv,
  //                                real ql, real qi) const {};

  // real YAKL_INLINE compute_dHdqv(real p, real entropic_var, real qd, real qv,
  //                                real ql, real qi) const {};

  // real YAKL_INLINE compute_dHdql(real p, real entropic_var, real qd, real qv,
  //                                real ql, real qi) const {};

  // real YAKL_INLINE compute_dHdqi(real p, real entropic_var, real qd, real qv,
  //                                real ql, real qi) const {};

  // real YAKL_INLINE compute_alpha(real p, real T, real qd, real qv, real ql,
  //                                real qi) const {};

  // real YAKL_INLINE compute_entropic_var_from_p_T(real p, real T, real qd,
  //                                                real qv, real ql,
  //                                                real qi) const {};

  // real YAKL_INLINE solve_p(real rho, real entropic_var, real qd, real qv,
  //                          real ql, real qi) const {};
};

#ifdef PAMC_THERMONONE
using ThermoPotential = ThermoNone;
#elif PAMC_IDEAL_GAS_POTTEMP
using ThermoPotential = IdealGas_Pottemp;
#elif PAMC_IDEAL_GAS_ENTROPY
using ThermoPotential = IdealGas_Entropy;
#elif PAMC_CONST_KAPPA_VIRPOTTEMP
using ThermoPotential = ConstantKappa_VirtualPottemp;
#elif PAMC_UNAPPROX_POTTEMP
using ThermoPotential = Unapprox_Pottemp;
#elif PAMC_UNAPPROX_ENTROPY
using ThermoPotential = Unapprox_Entropy;
#endif
} // namespace pamc
