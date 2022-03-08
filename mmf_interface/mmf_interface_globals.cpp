
#include "pam_coupler.h"
#include <vector>

namespace mmf_interface {
  // This is a vector because CPU threaded regions will require a different PamCoupler instance for each thread
  std::vector<pam::PamCoupler> couplers;

  std::function<void()> gcm_initialize = [] () { yakl::yakl_throw("ERROR: user has not set gcm_initialize()"); };
  std::function<void()> gcm_tend       = [] () { yakl::yakl_throw("ERROR: user has not set gcm_tend      ()"); };
  std::function<void()> gcm_finalize   = [] () { yakl::yakl_throw("ERROR: user has not set gcm_finalize  ()"); };
}


