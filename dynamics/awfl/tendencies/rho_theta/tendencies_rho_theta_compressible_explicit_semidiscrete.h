
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
namespace compressible_explicit_semidiscrete {


  void compute_tendencies_xyz( pam::PamCoupler const &coupler ,
                               realConst5d state   , real5d const &state_tend  ,
                               realConst5d tracers , real5d const &tracer_tend ,
                               realConst5d tracers_start ,
                               awfl::Recon const &recon ,
                               awfl::tendencies_rho_theta::Hydrostasis const &hydrostasis ,
                               real dt );


}
}
}


