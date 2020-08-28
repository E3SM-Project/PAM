
#pragma once

#include "const.h"
#include "Profiles.h"

class PhysicsSaturationAdjustment {
public:
  int static constexpr numTracers = 2;

  real static constexpr R_d  = 287.;
  real static constexpr cp_d = 1004.;
  real static constexpr cv_d = cp_d-R_d;
  real static constexpr gamma_d = cp_d / cv_d;
  real static constexpr kappa_d = R_d / cp_d;
  real static constexpr R_v  = 461.;
  real static constexpr cp_v = 1859;
  real static constexpr cv_v = R_v-cp_v;
  real static constexpr p0   = 1.e5;
  real C0_d;

  int dataSpec;
  int static constexpr DATA_SPEC_THERMAL_MOIST = 2;

  int static constexpr ID_V = 0;
  int static constexpr ID_C = 1;

  struct DynState {
    real rho;
    real w;
    real rho_theta;
  };


  typedef SArray<real,1,numTracers> MicroTracers;


  PhysicsSaturationAdjustment() {
    C0_d = pow( R_d*pow( p0 , -kappa_d ) , gamma_d );
  }


  void init(std::string inFile) {
  }


  template <class SP> void addTracers(SP &spaceOp) const {
    //                name             description      positive   adds mass
    spaceOp.addTracer("water_vapor"  , "Water Vapor"  , true     , true);
    spaceOp.addTracer("cloud_liquid" , "Cloud liquid" , true     , true);
  }


  template <class SP> void initTracers(SP &spaceOp, typename SP::TracerArr &tracers) const {
    auto initVaporMass = YAKL_LAMBDA (real x, real y, real z, real xlen, real ylen, real zlen, DynState const &state)->real {
      real pert = profiles::ellipsoid_linear(x,y,z  ,  xlen/2,ylen/2,2000  ,  2000,2000,2000  ,  0.8);
      real temp = tempFromRhoTheta(state.rho , state.rho_theta);
      real svp  = saturationVaporPressure(temp);
      real p_v  = pert*svp;
      real r_v  = p_v / (R_v*temp);
      return r_v;
    };
    auto initCloudMass = YAKL_LAMBDA (real x, real y, real z, real xlen, real ylen, real zlen, DynState const &state)->real {
      return 0;
    };

    spaceOp.initTracerByLocation("water_vapor"  , initVaporMass , tracers , *this);
    spaceOp.initTracerByLocation("cloud_liquid" , initCloudMass , tracers , *this);
  }


  // Returns saturation vapor pressure
  YAKL_INLINE real saturationVaporPressure(real temp) const {
    real tc = temp - 273.15;
    return 610.94 * exp( 17.625*tc / (243.04+tc) );
  }


  YAKL_INLINE real gasConstant(real rho, real rho_v, real rho_c) const {
    real rho_d = rho - rho_v - rho_c;
    real q_d = rho_d / rho;
    real q_v = rho_v / rho;
    return R_d*q_d + R_v*q_v;
  }


  YAKL_INLINE real specificHeatPressure(real rho, real rho_v, real rho_c) const {
    real rho_d = rho - rho_v - rho_c;
    real q_d = rho_d / rho;
    real q_v = rho_v / rho;
    return cp_d*q_d + cp_v*q_v;
  }


  // Returns the latent heat of condensation
  YAKL_INLINE real latentHeatCondensation(real temp) const {
    real tc = temp - 273.15;
    return (2500.8 - 2.36*tc + 0.0016*tc*tc - 0.00006*tc*tc*tc)*1000;
  }


  YAKL_INLINE real pressureFromRhoTheta(real rho_theta) const {
    return C0_d * pow( rho_theta , gamma_d );
  }


  YAKL_INLINE real pressureFromTemp(real rho , real rho_v , real rho_c , real temp) const {
    real rho_d = rho - rho_v - rho_c;
    return rho_d*R_d*temp + rho_v*R_v*temp;
  }


  YAKL_INLINE real thetaFromTemp(real rho , real rho_v , real rho_c , real temp) const {
    real p = pressureFromTemp(rho, rho_v, rho_c, temp);
    return temp * pow( p0/p , kappa_d );
  }


  YAKL_INLINE real thetaFromTemp(real rho , MicroTracers const &tracers , real temp) const {
    real p = pressureFromTemp(rho, tracers(ID_V), tracers(ID_C), temp);
    return temp * pow( p0/p , kappa_d );
  }


  YAKL_INLINE real tempFromRhoTheta(real rho , real rho_theta) const {
    real p = pressureFromRhoTheta(rho_theta);
    return (rho_theta/rho) * pow( p/p0 , kappa_d );
  }


  // Computes the state (vapor density, cloud liquid density, and temperature) that achieves
  // a factor of "ratio" of the saturated state
  YAKL_INLINE void computeAdjustedState(real rho_d , real &rho_v , real &rho_c , real &rho_theta) const {
    real tol = 1.e-6;
    if (std::is_same<real,double>::value) tol = 1.e-13;

    real rho = rho_d + rho_v + rho_c;

    real temp = tempFromRhoTheta(rho, rho_theta);
    real svp = saturationVaporPressure( temp );
    real pv = rho_v * R_v * temp;   // water vapor pressure
    
    if        (pv > svp) {  // If we are super-saturated
      // Condense enough water vapor to achieve saturation
      ////////////////////////////////////////////////////////
      // Bisection method
      ////////////////////////////////////////////////////////
      // Set bounds on how much mass to condense
      real cond1  = 0;
      real cond2 = rho_v;

      bool keepIterating = true;
      while (keepIterating) {
        real rho_cond = (cond1 + cond2) / 2;                // How much water vapor to condense
        real rvLoc = rho_v - rho_cond;                      // New vapor density
        real rcLoc = rho_c + rho_cond;                      // New cloud liquid density
        real Lv = latentHeatCondensation(temp);             // Compute latent heat of condensation for water
        real R = gasConstant          (rho, rvLoc, rcLoc);  // New moist gas constant
        real cp = specificHeatPressure(rho, rvLoc, rcLoc);  // New moist specific heat at constant pressure
        real tempLoc = temp + rho_cond*Lv/(rho*cp);         // New temperature after condensation
        real svpLoc = saturationVaporPressure(tempLoc);     // New saturation vapor pressure after condensation
        real pvLoc = rvLoc * R_v * tempLoc;                 // New vapor pressure after condensation
        // If we're supersaturated still, we need to condense out more water vapor
        // otherwise, we need to condense out less water vapor
        if (pvLoc > svpLoc) {
          cond1 = rho_cond;
        } else {
          cond2 = rho_cond;
        }
        // If we've converged, then we can stop iterating
        if (abs(cond2-cond1) <= tol) {
          rho_v = rvLoc;
          rho_c = rcLoc;
          rho = rho_d + rho_v + rho_c;
          rho_theta = rho * thetaFromTemp(rho, rho_v, rho_c, tempLoc);
          keepIterating = false;
        }
      }
    } else if (pv < svp && rho_c > 0) {  // If we are unsaturated and have cloud liquid
      // If there's cloud, evaporate enough to achieve saturation
      // or all of it if there isn't enough to reach saturation
      ////////////////////////////////////////////////////////
      // Bisection method
      ////////////////////////////////////////////////////////
      // Set bounds on how much mass to evaporate
      real evap1  = 0;
      real evap2 = rho_c;

      bool keepIterating = true;
      while (keepIterating) {
        real rho_evap = (evap1 + evap2) / 2;                 // How much water vapor to evapense
        real rvLoc = rho_v + rho_evap;                       // New vapor density
        real rcLoc = rho_c - rho_evap;                       // New cloud liquid density
        real Lv = latentHeatCondensation(temp);              // Compute latent heat of condensation for water
        real R  = gasConstant         (rho, rvLoc, rcLoc);   // New moist gas constant
        real cp = specificHeatPressure(rho, rvLoc, rcLoc);   // New moist specific heat
        real tempLoc = temp - rho_evap*Lv/(rho*cp);          // New temperature after evaporation
        real svpLoc = saturationVaporPressure(tempLoc);      // New saturation vapor pressure after evaporation
        real pvLoc = rvLoc * R_v * tempLoc;                  // New vapor pressure after evaporation
        // If we're unsaturated still, we need to evaporate out more water vapor
        // otherwise, we need to evaporate out less water vapor
        if (pvLoc < svpLoc) {
          evap1 = rho_evap;
        } else {
          evap2 = rho_evap;
        }
        // If we've converged, then we can stop iterating
        if (abs(evap2-evap1) <= tol) {
          rho_v = rvLoc;
          rho_c = rcLoc;
          rho = rho_d + rho_v + rho_c;
          rho_theta = rho * thetaFromTemp(rho, rho_v, rho_c, tempLoc);
          keepIterating = false;
        }
      }
    }
  }


};


