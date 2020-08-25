
#include "const.h"

class PhysicsSaturationAdjustment {
  real constexpr R_d  = 287.;
  real constexpr cp_d = 1004.;
  real constexpr cv_d = cp_d-R_d;
  real constexpr R_v  = 461.;
  real constexpr cp_v = 1859;
  real constexpr cv_v = R_v-cp_v;
  real constexpr p0   = 1.e5;


  // Returns saturation vapor pressure
  YAKL_INLINE real computeSatVapPress(real temp) {
    real tc = temp - 273.15;
    return 610.94 * exp( 17.625*tc / (243.04+tc) );
  }


  // Returns the latent heat of condensation
  YAKL_INLINE real computeLatentCondense(real temp) {
    real tc = temp - 273.15;
    return (2500.8 - 2.36*tc + 0.0016*tc*tc - 0.00006*tc*tc*tc)*1000;
  }


  // Computes moist ideal gas constant and moist specific heat at constant pressure
  YAKL_INLINE void computeR_cp(real rho, real rho_vapor, real &R, real &cp) {
    real q_d = (rho - rho_vapor) / rho; // Wet mixing ratio of dry air
    real q_v = rho_vapor / rho;         // Wet mixing ratio of water vapor
    R  = q_d*R_d  + q_v*R_v ;
    cp = q_d*cp_d + q_v*cp_v;
  }


  // Returns total pressure
  YAKL_INLINE real computePressFromTheta(real rho , real rho_vapor , real rho_theta) {
    real R, cp;
    computeR_cp(rho, rho_vapor, R, cp);
    real gamma = cp / cv;
    return p = pow( R*pow( p0 , -R/cp ) * rho_theta , gamma );
  }


  // Returns total pressure
  YAKL_INLINE real computePressFromTemp(real rho , real rho_vapor , real temp) {
    real rho_d = rho - rho_vapor;
    return rho_d*R_d*temp + rho_vapor*R_v*temp;
  }


  // Returns total pressure
  YAKL_INLINE real computeThetaFromTemp(real rho , real rho_vapor , real temp) {
    real p = computePressFromTemp(rho, rho_vapor, temp);
    real R, cp;
    computeR_cp(rho, rho_vapor, R, cp);
    return temp * pow( p0/p , R/cp );
  }


  // Computes temperature
  YAKL_INLINE real computeTempFromTheta(real rho, real rho_vapor, real rho_theta) {
    real p = computePressFromTheta(rho, rho_vapor, rho_theta);
    real R, cp;
    computeR_cp(rho, rho_vapor, R, cp);
    return (rho_theta/rho) * pow( p/p0 , R/cp );
  }


  // Computes the state (vapor density, cloud liquid density, and temperature) that achieves
  // a factor of "ratio" of the saturated state
  YAKL_INLINE void computeAdjustedState(real rho , real &rho_vapor , real &rho_cloud , real &rho_theta) {
    real temp = computeTempFromTheta(rho, rho_vapor, rho_theta);
    real svp = computeSatVapPress( temp );
    real pv = rho_vapor * R_v * temp;
    
    if        (pv > svp) {  // If we are super-saturated
      // Condense enough water vapor to achieve saturation
      ////////////////////////////////////////////////////////
      // Bisection method
      ////////////////////////////////////////////////////////
      // Set bounds on how much mass to condense
      real cond1  = 0;
      real cond2 = rho_vapor;

      bool keepIterating = true;
      while (keepIterating) {
        real rho_cond = (cond1 + cond2) / 2; // How much water vapor to condense
        real rvLoc = rho_vapor - rho_cond;
        real rcLoc = rho_cloud + rho_cond;
        real Lv = computeLatentCondense(temp);    // Compute latent heat of condensation for water
        real tempLoc = temp + rho_cond*Lv/(rho*cp);
        real svpLoc = computeSatVapPress(tempLoc);
        // If we're supersaturated still, we need to condense out more water vapor
        // otherwise, we need to condense out less water vapor
        if (rvLoc > svpLoc) {
          cond1 = rho_cond;
        } else {
          cond2 = rho_cond;
        }
        // If we've converged, then we can stop iterating
        if (abs(cond2-cond1) <= 1.e-6) {
          rho_vapor = rvLoc;
          rho_cloud = rcLoc;
          rho_theta = rho * computeThetaFromTemp(rho, rho_vapor, tempLoc);
          keepIterating = false;
        }
      }
    } else if (pv < svp) {  // If we are unsaturated
      // If there's cloud, evaporate enough to achieve saturation
      // or all of it if there isn't enough to reach saturation
      ////////////////////////////////////////////////////////
      // Bisection method
      ////////////////////////////////////////////////////////
      // Set bounds on how much mass to evaporate
      real evap1  = 0;
      real evap2 = rho_cloud;

      bool keepIterating = true;
      while (keepIterating) {
        real rho_evap = (evap1 + evap2) / 2; // How much water vapor to evapense
        real rvLoc = rho_vapor + rho_evap;
        real rcLoc = rho_cloud - rho_evap;
        real Lv = computeLatentCondense(temp);    // Compute latent heat of condensation for water
        real tempLoc = temp - rho_evap*Lv/(rho*cp);
        real svpLoc = computeSatVapPress(tempLoc);
        // If we're unsaturated still, we need to evaporate out more water vapor
        // otherwise, we need to evaporate out less water vapor
        if (rvLoc < svpLoc) {
          evap1 = rho_evap;
        } else {
          evap2 = rho_evap;
        }
        // If we've converged, then we can stop iterating
        if (abs(evap2-evap1) <= 1.e-6) {
          rho_vapor = rvLoc;
          rho_cloud = rcLoc;
          rho_theta = rho * computeThetaFromTemp(rho, rho_vapor, tempLoc);
          keepIterating = false;
        }
      }
    }
  }


};


