
#pragma once

#include "pam_const.h"
#include "pam_coupler.h"


YAKL_INLINE static real saturation_vapor_pressure(real temp) {
  real tc = temp - 273.15;
  return 610.94 * exp( 17.625*tc / (243.04+tc) );
}


YAKL_INLINE static real latent_heat_condensation(real temp) {
  real tc = temp - 273.15;
  return (2500.8 - 2.36*tc + 0.0016*tc*tc - 0.00006*tc*tc*tc)*1000;
}


YAKL_INLINE static real cp_moist(real rho_d, real rho_v, real rho_c, real cp_d, real cp_v, real cp_l) {
  // For the moist specific heat, ignore other species than water vapor and cloud droplets
  real rho = rho_d + rho_v + rho_c;
  return rho_d / rho * cp_d  +  rho_v / rho * cp_v  +  rho_c / rho * cp_l;
}

// Compute an instantaneous adjustment of sub or super saturation
YAKL_INLINE static void compute_adjusted_state(real rho, real rho_d , real &rho_v , real &rho_c , real &temp,
                                               real R_v , real cp_d , real cp_v , real cp_l) {
  // Define a tolerance for convergence
  real tol = 1.e-6;

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
      real rho_cond = (cond1 + cond2) / 2;                    // How much water vapor to condense for this iteration
      real rv_loc = max( 0._fp , rho_v - rho_cond );          // New vapor density
      real rc_loc = max( 0._fp , rho_c + rho_cond );          // New cloud liquid density
      real Lv = latent_heat_condensation(temp);               // Compute latent heat of condensation
      real cp = cp_moist(rho_d,rv_loc,rc_loc,cp_d,cp_v,cp_l); // New moist specific heat at constant pressure
      real temp_loc = temp + rho_cond*Lv/(rho*cp);            // New temperature after condensation
      real svp_loc = saturation_vapor_pressure(temp_loc);     // New saturation vapor pressure after condensation
      real pv_loc = rv_loc * R_v * temp_loc;                  // New vapor pressure after condensation
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
        temp  = temp_loc;
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
      real rho_evap = (evap1 + evap2) / 2;                    // How much water vapor to evapense
      real rv_loc = max( 0._fp , rho_v + rho_evap );          // New vapor density
      real rc_loc = max( 0._fp , rho_c - rho_evap );          // New cloud liquid density
      real Lv = latent_heat_condensation(temp);               // Compute latent heat of condensation for water
      real cp = cp_moist(rho_d,rv_loc,rc_loc,cp_d,cp_v,cp_l); // New moist specific heat
      real temp_loc = temp - rho_evap*Lv/(rho*cp);            // New temperature after evaporation
      real svp_loc = saturation_vapor_pressure(temp_loc);     // New saturation vapor pressure after evaporation
      real pv_loc = rv_loc * R_v * temp_loc;                  // New vapor pressure after evaporation
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
        temp  = temp_loc;
        keep_iterating = false;
      }
    }
  }
}


inline void saturation_adjustment( PamCoupler &coupler , real dt ) {
  using yakl::intrinsics::size;
  real1d rho_d = coupler.dm.get_collapsed<real>("density_dry");
  real1d temp  = coupler.dm.get_collapsed<real>("temp");
  real1d rho_v = coupler.dm.get_collapsed<real>("water_vapor");
  real1d rho_c;
  std::string micro_scheme = coupler.get_option<std::string>("micro");
  if      (micro_scheme == "kessler") { rho_c = coupler.dm.get_collapsed<real>( "cloud_liquid" ); }
  else if (micro_scheme == "p3"     ) { rho_c = coupler.dm.get_collapsed<real>( "cloud_water"  ); }
  else { endrun("ERROR: saturation_adjustment.h only currently supports kessler and p3 microphysics"); }
  auto tracer_names = coupler.get_tracer_names();
  pam::MultiField<real,1> massy_tracers;
  for (int tr=0; tr < tracer_names.size(); tr++) {
    std::string desc;
    bool        found, positive, adds_mass;
    coupler.get_tracer_info(tracer_names[tr],desc,found,positive,adds_mass);
    if (adds_mass) massy_tracers.add_field( coupler.dm.get_collapsed<real>(tracer_names[tr]) );
  }
  real R_v  = coupler.R_v ;
  real cp_d = coupler.cp_d;
  real cp_v = coupler.cp_v;
  real cp_l = 4188.0;
  parallel_for( SimpleBounds<1>(size(rho_d)) , YAKL_LAMBDA (int i) {
    real rho = rho_d(i);
    for (int tr=0; tr < massy_tracers.get_num_fields(); tr++) { rho += massy_tracers(tr,i); }
    compute_adjusted_state(rho, rho_d(i), rho_v(i), rho_c(i), temp(i), R_v, cp_d, cp_v, cp_l);
  });
}
