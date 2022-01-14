
#include "PamCoupler.h"
#include <vector>

// This is a vector because CPU threaded regions will require a different PamCoupler instance for each thread
std::vector<PamCoupler> pam_interface_couplers;


