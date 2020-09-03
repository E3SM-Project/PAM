
#pragma once

#include "const.h"
#include "Profiles.h"
#include "DataManager.h"

class Microphysics_saturation_adjustment {
public:
  int static constexpr num_tracers = 2;

  real R_d        ;
  real cp_d       ;
  real cv_d       ;
  real gamma_d    ;
  real kappa_d    ;
  real R_v        ;
  real cp_v       ;
  real cv_v       ;
  real p0         ;
  real C0_d       ;
  int  tracer_index_vapor;
  real p0_nkappa_d;



  Microphysics_saturation_adjustment() {
    R_d     = 287.;
    cp_d    = 1004.;
    cv_d    = cp_d-R_d;
    gamma_d = cp_d / cv_d;
    kappa_d = R_d / cp_d;
    R_v     = 461.;
    cp_v    = 1859;
    cv_v    = R_v-cp_v;
    p0      = 1.e5;
    p0_nkappa_d = pow( p0 , -kappa_d );
    C0_d        = pow( R_d * p0_nkappa_d , gamma_d );
    tracer_index_vapor = -1;
  }



  static constexpr int get_num_tracers() {
    return num_tracers;
  }



  template <class SP>
  void init(std::string inFile , SP &spaceOp , DataManager &dm) {
    // Register tracers in the dycore
    //                                           name             description      positive   adds mass
    tracer_index_vapor = spaceOp.add_tracer(dm , "water_vapor"  , "Water Vapor"  , true     , true);
    int dummy          = spaceOp.add_tracer(dm , "cloud_liquid" , "Cloud liquid" , true     , true);
  }



  template <class SP> void init_tracers(SP &spaceOp, DataManager &dm) const {
    auto init_vapor_mass = YAKL_LAMBDA (real x, real y, real z, real xlen, real ylen, real zlen,
                                        real rho, real rho_theta)->real {
      real pert = profiles::ellipsoid_linear(x,y,z  ,  xlen/2,ylen/2,2000  ,  2000,2000,2000  ,  0.8);
      real temp = temp_from_rho_theta(rho , 0 , rho_theta);
      real svp  = saturation_vapor_pressure(temp);
      real p_v  = pert*svp;
      real r_v  = p_v / (R_v*temp);
      return r_v;
    };
    auto init_cloud_mass = YAKL_LAMBDA (real x, real y, real z, real xlen, real ylen, real zlen,
                                        real rho, real rho_theta)->real {
      return 0;
    };

    spaceOp.init_tracer_by_location("water_vapor"  , init_vapor_mass , dm , *this);
    spaceOp.init_tracer_by_location("cloud_liquid" , init_cloud_mass , dm , *this);
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



  YAKL_INLINE real R_moist(real rho, real rho_v) const {
    real rho_d = rho - rho_v;
    real q_d = rho_d / rho;
    real q_v = rho_v / rho;
    return R_d * q_d + R_v * q_v;
  }



  YAKL_INLINE real cp_moist(real rho, real rho_v) const {
    return R_moist(rho, rho_v) / R_d * cp_d;
  }



  YAKL_INLINE real cv_moist(real rho, real rho_v) const {
    return R_moist(rho, rho_v) / R_d * cv_d;
  }



  YAKL_INLINE real pressure_C0_moist(real rho, real rho_v) const {
    return pow( R_moist(rho, rho_v) * p0_nkappa_d , gamma_d );
  }



  YAKL_INLINE real pressure_from_rho_theta(real rho, real rho_v, real rho_theta) const {
    real C0 = pressure_C0_moist(rho, rho_v);
    return C0 * pow( rho_theta , gamma_d );
  }



  YAKL_INLINE real pressure_from_temp(real rho , real rho_v , real temp) const {
    real rho_d = rho - rho_v;
    return rho_d*R_d*temp + rho_v*R_v*temp;
  }



  YAKL_INLINE real theta_from_temp(real rho , real rho_v , real temp) const {
    real p = pressure_from_temp(rho, rho_v, temp);
    return temp * pow( p0/p , kappa_d );
  }



  YAKL_INLINE real temp_from_rho_theta(real rho , real rho_v , real rho_theta) const {
    real p = pressure_from_rho_theta(rho, rho_v, rho_theta);
    return (rho_theta/rho) * pow( p/p0 , kappa_d );
  }



  // Computes the state (vapor density, cloud liquid density, and temperature) that achieves
  // a factor of "ratio" of the saturated state
  YAKL_INLINE void compute_adjusted_state(real rho , real &rho_v , real &rho_c , real &rho_theta) const {
    // Define a tolerance for convergence
    real tol = 1.e-6;
    if (std::is_same<real,double>::value) tol = 1.e-13;

    // Temperature before adjustment (before latent heating or cooling)
    real temp = temp_from_rho_theta(rho, rho_v, rho_theta);

    // Saturation vapor pressure at this temperature
    real svp = saturation_vapor_pressure( temp );

    // Vapor pressure at this temperature
    real pv = rho_v * R_v * temp;
    
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
        real R  = R_moist (rho, rv_loc);                     // New moist gas constant
        real cp = cp_moist(rho, rv_loc);                     // New moist specific heat at constant pressure
        real temp_loc = temp + rho_cond*Lv/(rho*cp);         // New temperature after condensation
        real svp_loc = saturation_vapor_pressure(temp_loc);  // New saturation vapor pressure after condensation
        real pv_loc = rv_loc * R_v * temp_loc;               // New vapor pressure after condensation
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
          rho_theta = rho * theta_from_temp(rho, rho_v, temp_loc);
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
        real R  = R_moist (rho, rv_loc);                     // New moist gas constant
        real cp = cp_moist(rho, rv_loc);                     // New moist specific heat
        real temp_loc = temp - rho_evap*Lv/(rho*cp);         // New temperature after evaporation
        real svp_loc = saturation_vapor_pressure(temp_loc);  // New saturation vapor pressure after evaporation
        real pv_loc = rv_loc * R_v * temp_loc;               // New vapor pressure after evaporation
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
          rho_theta = rho * theta_from_temp(rho, rho_v, temp_loc);
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

    int num_cells = rho.totElems();
    parallel_for( num_cells , YAKL_LAMBDA (int i) {
      compute_adjusted_state( rho(i) , rho_v(i) , rho_c(i) , rho_theta(i) );
    });
  }


};


