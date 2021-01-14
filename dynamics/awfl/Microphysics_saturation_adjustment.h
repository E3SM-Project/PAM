
#pragma once

#include "const.h"
#include "DataManager.h"

class Microphysics {
public:
  int static constexpr num_tracers = 2;

  struct Constants {
    real R_d               ;
    real cp_d              ;
    real cv_d              ;
    real gamma_d           ;
    real kappa_d           ;
    real R_v               ;
    real cp_v              ;
    real cv_v              ;
    real p0                ;
    real C0_d              ;
    real p0_nkappa_d       ;
  };

  int  tracer_index_vapor;

  Constants constants;

  SArray<real,1,num_tracers> tracer_IDs; // tracer index for microphysics tracers

  int static constexpr ID_V = 0;  // Local index for water vapor
  int static constexpr ID_C = 1;  // Local index for cloud liquid



  Microphysics() {
    constants.R_d         = 287.;
    constants.cp_d        = 1004.;
    constants.cv_d        = constants.cp_d - constants.R_d;
    constants.gamma_d     = constants.cp_d / constants.cv_d;
    constants.kappa_d     = constants.R_d  / constants.cp_d;
    constants.R_v         = 461.;
    constants.cp_v        = 1859;
    constants.cv_v        = constants.R_v - constants.cp_v;
    constants.p0          = 1.e5;
    constants.p0_nkappa_d = pow( constants.p0 , -constants.kappa_d );
    constants.C0_d        = pow( constants.R_d * constants.p0_nkappa_d , constants.gamma_d );

    tracer_index_vapor = -1;
  }



  static constexpr int get_num_tracers() {
    return num_tracers;
  }



  template <class DC>
  void init(std::string inFile , DC &dycore , DataManager &dm) {
    // Register tracers in the dycore
    //                                          name             description      positive   adds mass
    tracer_IDs(ID_V) = dycore.add_tracer(dm , "water_vapor"  , "Water Vapor"  , true     , true);
    tracer_IDs(ID_C) = dycore.add_tracer(dm , "cloud_liquid" , "Cloud liquid" , true     , true);
    tracer_index_vapor = tracer_IDs(ID_V);
  }



  // Returns saturation vapor pressure
  YAKL_INLINE real saturation_vapor_pressure(real temp) const {
    real tc = temp - 273.15;
    return 610.94 * exp( 17.625*tc / (243.04+tc) );
  }



  YAKL_INLINE real latent_heat_condensation(real temp) const {
    real tc = temp - 273.15;
    return (2500.8 - 2.36*tc + 0.0016*tc*tc - 0.00006*tc*tc*tc)*1000;
  }



  YAKL_INLINE real R_moist(real rho, real rho_v, Constants const &cn) const {
    real rho_d = rho - rho_v;
    real q_d = rho_d / rho;
    real q_v = rho_v / rho;
    return cn.R_d * q_d + cn.R_v * q_v;
  }



  YAKL_INLINE real cp_moist(real rho, real rho_v, Constants const &cn) const {
    return R_moist(rho, rho_v, cn) / cn.R_d * cn.cp_d;
  }



  YAKL_INLINE real cv_moist(real rho, real rho_v, Constants const &cn) const {
    return R_moist(rho, rho_v, cn) / cn.R_d * cn.cv_d;
  }



  YAKL_INLINE real pressure_from_rho_theta(real rho_theta, Constants const &cn) const {
    return cn.C0_d * pow( rho_theta , cn.gamma_d );
  }



  YAKL_INLINE real pressure_from_temp(real rho , real rho_v , real temp, Constants const &cn) const {
    real rho_d = rho - rho_v;
    return rho_d*cn.R_d*temp + rho_v*cn.R_v*temp;
  }



  YAKL_INLINE real theta_from_temp(real rho , real rho_v , real temp, Constants const &cn) const {
    real p = pressure_from_temp(rho, rho_v, temp, cn);
    return pow( p/cn.C0_d , 1./cn.gamma_d ) / rho;
  }



  YAKL_INLINE real temp_from_rho_theta(real rho , real rho_v , real rho_theta, Constants const &cn) const {
    real p = pressure_from_rho_theta(rho_theta, cn);
    real R = R_moist(rho, rho_v, cn);
    return p / rho / R;
  }



  // Computes the state (vapor density, cloud liquid density, and temperature) that achieves
  // a factor of "ratio" of the saturated state
  YAKL_INLINE void compute_adjusted_state(real rho , real &rho_v , real &rho_c , real &rho_theta, Constants const &cn) const {
    // Define a tolerance for convergence
    real tol = 1.e-6;
    if (std::is_same<real,double>::value) tol = 1.e-13;

    // Temperature before adjustment (before latent heating or cooling)
    real temp = temp_from_rho_theta(rho, rho_v, rho_theta, cn);

    // Saturation vapor pressure at this temperature
    real svp = saturation_vapor_pressure( temp );

    // Vapor pressure at this temperature
    real pv = rho_v * cn.R_v * temp;
    
    // If we're super-saturated, we need to condense until saturation is reached
    if        (pv > svp) {
      ////////////////////////////////////////////////////////
      // Bisection method
      ////////////////////////////////////////////////////////
      // Set bounds on how much mass to condense
      real cond1  = 0;     // Minimum amount we can condense out
      real cond2 = rho_v;  // Maximum amount we can condense out

      bool keep_iterating = true;
      while (keep_iterating) {
        real rho_cond = (cond1 + cond2) / 2;                 // How much water vapor to condense for this iteration
        real rv_loc = rho_v - rho_cond;                      // New vapor density
        real rc_loc = rho_c + rho_cond;                      // New cloud liquid density
        real Lv = latent_heat_condensation(temp);            // Compute latent heat of condensation
        real R  = R_moist (rho, rv_loc, cn);                 // New moist gas constant
        real cp = cp_moist(rho, rv_loc, cn);                 // New moist specific heat at constant pressure
        real temp_loc = temp + rho_cond*Lv/(rho*cp);         // New temperature after condensation
        real svp_loc = saturation_vapor_pressure(temp_loc);  // New saturation vapor pressure after condensation
        real pv_loc = rv_loc * cn.R_v * temp_loc;            // New vapor pressure after condensation
        // If we're supersaturated still, we need to condense out more water vapor
        // otherwise, we need to condense out less water vapor
        if (pv_loc > svp_loc) {
          cond1 = rho_cond;
        } else {
          cond2 = rho_cond;
        }
        // If we've converged, then we can stop iterating
        if (abs(cond2-cond1) <= tol) {
          rho_v = rv_loc;
          rho_c = rc_loc;
          rho_theta = rho * theta_from_temp(rho, rho_v, temp_loc, cn);
          keep_iterating = false;
        }
      }

    // If we are unsaturated and have cloud liquid
    } else if (pv < svp && rho_c > 0) {
      // If there's cloud, evaporate enough to achieve saturation
      // or all of it if there isn't enough to reach saturation
      ////////////////////////////////////////////////////////
      // Bisection method
      ////////////////////////////////////////////////////////
      // Set bounds on how much mass to evaporate
      real evap1 = 0;     // minimum amount we can evaporate
      real evap2 = rho_c; // maximum amount we can evaporate

      bool keep_iterating = true;
      while (keep_iterating) {
        real rho_evap = (evap1 + evap2) / 2;                 // How much water vapor to evapense
        real rv_loc = rho_v + rho_evap;                      // New vapor density
        real rc_loc = rho_c - rho_evap;                      // New cloud liquid density
        real Lv = latent_heat_condensation(temp);            // Compute latent heat of condensation for water
        real R  = R_moist (rho, rv_loc, cn);                 // New moist gas constant
        real cp = cp_moist(rho, rv_loc, cn);                 // New moist specific heat
        real temp_loc = temp - rho_evap*Lv/(rho*cp);         // New temperature after evaporation
        real svp_loc = saturation_vapor_pressure(temp_loc);  // New saturation vapor pressure after evaporation
        real pv_loc = rv_loc * cn.R_v * temp_loc;            // New vapor pressure after evaporation
        // If we're unsaturated still, we need to evaporate out more water vapor
        // otherwise, we need to evaporate out less water vapor
        if (pv_loc < svp_loc) {
          evap1 = rho_evap;
        } else {
          evap2 = rho_evap;
        }
        // If we've converged, then we can stop iterating
        if (abs(evap2-evap1) <= tol) {
          rho_v = rv_loc;
          rho_c = rc_loc;
          rho_theta = rho * theta_from_temp(rho, rho_v, temp_loc, cn);
          keep_iterating = false;
        }
      }
    }
  }



  void timeStep( DataManager &dm , real dt ) const {
    auto rho       = dm.get_collapsed<real>("density");
    auto rho_theta = dm.get_collapsed<real>("density_theta");
    auto rho_v     = dm.get_collapsed<real>("water_vapor");
    auto rho_c     = dm.get_collapsed<real>("cloud_liquid");

    auto &constants = this->constants;

    int num_cells = rho.totElems();
    parallel_for( num_cells , YAKL_LAMBDA (int i) {
      compute_adjusted_state( rho(i) , rho_v(i) , rho_c(i) , rho_theta(i) , constants);
    });
  }


};


