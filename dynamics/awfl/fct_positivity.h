
#pragma once

#include "pam_coupler.h"
#include "awfl_const.h"

namespace awfl {
  void fct_positivity_x( pam::PamCoupler const &coupler, realConst5d tracers, real5d const &tracer_flux ,
                         boolConst1d tracer_pos, real dt);

  void fct_positivity_y( pam::PamCoupler const &coupler, realConst5d tracers, real5d const &tracer_flux ,
                         boolConst1d tracer_pos, real dt);

  void fct_positivity_z( pam::PamCoupler const &coupler, realConst5d tracers, real5d const &tracer_flux ,
                         boolConst1d tracer_pos, real dt);

  void fct_positivity_xyz( pam::PamCoupler const &coupler ,
                           realConst5d tracers            ,
                           real5d const &tracer_flux_x    ,
                           real5d const &tracer_flux_y    ,
                           real5d const &tracer_flux_z    ,
                           boolConst1d tracer_pos, real dt);
}


