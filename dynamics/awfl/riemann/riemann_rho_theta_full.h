
#pragma once

#include "pam_coupler.h"
#include "awfl_const.h"

namespace awfl {

  void riemann_rho_theta_full_x( pam::PamCoupler const &coupler ,
                                 real6d const & state_limits    ,
                                 real5d const &state_flux       ,
                                 real6d const & tracer_limits   ,
                                 real5d const &tracer_flux      );

  void riemann_rho_theta_full_y( pam::PamCoupler const &coupler ,
                                 real6d const & state_limits    ,
                                 real5d const &state_flux       ,
                                 real6d const & tracer_limits   ,
                                 real5d const &tracer_flux      );

  void riemann_rho_theta_full_z( pam::PamCoupler const &coupler ,
                                 real6d const & state_limits    ,
                                 real5d const &state_flux       ,
                                 real6d const & tracer_limits   ,
                                 real5d const &tracer_flux      );

  void riemann_rho_theta_full_xyz( pam::PamCoupler const &coupler  ,
                                   real6d const & state_limits_x   ,
                                   real5d const &state_flux_x      ,
                                   real6d const & tracer_limits_x  ,
                                   real5d const &tracer_flux_x     ,
                                   real6d const & state_limits_y   ,
                                   real5d const &state_flux_y      ,
                                   real6d const & tracer_limits_y  ,
                                   real5d const &tracer_flux_y     ,
                                   real6d const & state_limits_z   ,
                                   real5d const &state_flux_z      ,
                                   real6d const & tracer_limits_z  ,
                                   real5d const &tracer_flux_z     );

}


