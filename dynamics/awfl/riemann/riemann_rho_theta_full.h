
#pragma once

#include "pam_coupler.h"
#include "awfl_const.h"

namespace awfl {

  void riemann_rho_theta_full_x(pam::PamCoupler const &coupler, real6d const & state_limits  , real5d const &state_flux  ,
                                                                real6d const & tracer_limits , real5d const &tracer_flux );

  void riemann_rho_theta_full_y(pam::PamCoupler const &coupler, real6d const & state_limits  , real5d const &state_flux  ,
                                                                real6d const & tracer_limits , real5d const &tracer_flux );

  void riemann_rho_theta_full_z(pam::PamCoupler const &coupler, real6d const & state_limits  , real5d const &state_flux  ,
                                                                real6d const & tracer_limits , real5d const &tracer_flux );

}


