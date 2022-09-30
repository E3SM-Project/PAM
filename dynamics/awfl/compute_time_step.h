
#pragma once

#include "pam_coupler.h"

namespace awfl {
  real compute_time_step(pam::PamCoupler const &coupler, real cfl);
}

