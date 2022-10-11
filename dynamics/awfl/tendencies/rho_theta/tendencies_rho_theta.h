
#pragma once

#include "pam_coupler.h"
#include "MultipleFields.h"
#include "awfl_const.h"
#include "TransformMatrices.h"
#include "idealized_profiles.h"

namespace awfl {
namespace tendencies_rho_theta {

  int static constexpr num_state   = 5;
  int static constexpr max_tracers = 50;
  // For indexing into the state and state tendency arrays
  int static constexpr idR = 0;  // density
  int static constexpr idU = 1;  // rho*u
  int static constexpr idV = 2;  // rho*v
  int static constexpr idW = 3;  // rho*w
  int static constexpr idT = 4;  // rho*potential temperature



  // Hydrostasis data
  struct Hydrostasis {
    real hydrostasis_parameters_sum;  // Sum of the current parameters (to see if it's changed)
    real3d hyDensSten;                // A stencil around each cell of hydrostatic density
    real3d hyDensThetaSten;           // A stencil around each cell of hydrostatic density * potential temperature
    real3d hyDensGLL;                 // GLL point values of hydrostatic background density in each cell
    real3d hyDensThetaGLL;            // GLL point values of hydrostatic background density*potential temperature
    void allocate_and_initialize(pam::PamCoupler const &coupler) {
      auto nz   = coupler.get_nz  ();
      auto nens = coupler.get_nens();
      hyDensSten      = real3d("hyDensSten       ",nz,ord,nens);
      hyDensThetaSten = real3d("hyDensThetaSten  ",nz,ord,nens);
      hyDensGLL       = real3d("hyDensGLL        ",nz,ngll,nens);
      hyDensThetaGLL  = real3d("hyDensThetaGLL   ",nz,ngll,nens);
      hydrostasis_parameters_sum = 0;
    }
  };



  YAKL_INLINE static real hydrostatic_dens_theta( realConst3d hy_params , real z , real z0 , real dz ,
                                                  int k, int iens , real C0 , real gamma ) {
    real p = pam::hydrostatic_pressure( hy_params , z , z0 , dz , k , iens );
    // p = C0*(rho*theta)^gamma
    return pow(p/C0,1._fp/gamma);
  }



  void convert_dynamics_to_coupler_state( pam::PamCoupler &coupler ,
                                          realConst5d state ,
                                          realConst5d tracers );



  void convert_coupler_state_to_dynamics( pam::PamCoupler const &coupler ,
                                          real5d const &state ,
                                          real5d const &tracers ,
                                          Hydrostasis &hydrostasis );



  void init_idealized_state_and_tracers( pam::PamCoupler &coupler );



  real5d createStateArr(pam::PamCoupler const &coupler);



  real5d createTracerArr(pam::PamCoupler const &coupler);



  real5d createStateTendArr(pam::PamCoupler const &coupler);



  real5d createTracerTendArr(pam::PamCoupler const &coupler);



  std::vector<real> compute_mass( pam::PamCoupler const &coupler , realConst5d state , realConst5d tracers );

}
}


