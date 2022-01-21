
#include "pam_coupler.h"
#include <vector>

// This is a vector because CPU threaded regions will require a different PamCoupler instance for each thread
std::vector<pam::PamCoupler> pam_interface_couplers;


