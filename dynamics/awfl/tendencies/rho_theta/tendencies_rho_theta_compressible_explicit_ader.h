
#pragma once

#include "tendencies_rho_theta.h"
#include "reconstruct.h"
#include "riemann_rho_theta_full.h"
#include "diff_trans_tracer.h"
#include "diff_trans_rho_theta_full.h"
#include "compute_time_average.h"
#include "fct_positivity.h"


namespace awfl {
namespace tendencies_rho_theta {
namespace compressible_explicit_ader {



  void compute_tendencies_x( pam::PamCoupler const &coupler,
                             real5d const &state   , real5d const &state_tend  ,
                             real5d const &tracers , real5d const &tracer_tend ,
                             awfl::Recon const &recon , real &dt );



  void compute_tendencies_y( pam::PamCoupler const &coupler ,
                             real5d const &state   , real5d const &state_tend  ,
                             real5d const &tracers , real5d const &tracer_tend ,
                             awfl::Recon const &recon , real &dt );



  void compute_tendencies_z( pam::PamCoupler const &coupler ,
                             real5d const &state   , real5d const &state_tend       ,
                             real5d const &tracers , real5d const &tracer_tend      ,
                             awfl::Recon const &recon ,
                             awfl::tendencies_rho_theta::Hydrostasis const &hydrostasis ,
                             real &dt );



}
}
}


