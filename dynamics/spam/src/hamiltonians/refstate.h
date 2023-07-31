#pragma once

#include "profiles.h"
#include "topology.h"

namespace pamc {

struct ReferenceState_SWE {
  real ref_height;

#ifdef PAMC_EXTRUDED
  Profile dens;
  Profile geop;
  Profile q_di;
  Profile q_pi;
  Profile rho_di;
  Profile rho_pi;
  Profile Nsq_pi;
  Profile B;
#endif
  bool is_initialized = false;

  template <class VS>
  void initialize(const Topology &primal_topology,
                  const Topology &dual_topology) {
#ifdef PAMC_EXTRUDED
    this->dens.initialize(dual_topology, "ref dens", 1, 1, VS::ndensity);
    this->geop.initialize(dual_topology, "ref geop", 1, 1, 1);
    this->rho_pi.initialize(primal_topology, "refrho_pi", 0, 0, 1);
    this->q_pi.initialize(primal_topology, "refq_pi", 0, 0, VS::ndensity);
    this->rho_di.initialize(dual_topology, "refrho_di", 0, 0, 1);
    this->q_di.initialize(dual_topology, "refq_di", 0, 0, VS::ndensity);
    this->Nsq_pi.initialize(primal_topology, "refNsq_pi", 0, 0, 1);
#endif

    this->is_initialized = true;
  }
};

struct ReferenceState_Euler {
  Profile dens;
  Profile geop;
  Profile q_di;
  Profile q_pi;
  Profile rho_di;
  Profile rho_pi;
  Profile Nsq_pi;
  Profile B;
  bool is_initialized = false;

  template <class VS>
  void initialize(const Topology &primal_topology,
                  const Topology &dual_topology) {

    this->dens.initialize(dual_topology, "ref dens", 1, 1, VS::ndensity);
    this->geop.initialize(dual_topology, "ref geop", 1, 1, 1);
    this->rho_pi.initialize(primal_topology, "refrho_pi", 0, 0, 1);
    this->q_pi.initialize(primal_topology, "refq_pi", 0, 0, VS::ndensity);
    this->rho_di.initialize(dual_topology, "refrho_di", 0, 0, 1);
    this->q_di.initialize(dual_topology, "refq_di", 0, 0, VS::ndensity);
    this->Nsq_pi.initialize(primal_topology, "refNsq_pi", 0, 0, 1);
    this->B.initialize(dual_topology, "ref B", 1, 1, VS::ndensity_active);

    this->is_initialized = true;
  }
};

#if defined PAMC_SWE || defined PAMC_TSWE
using ReferenceState = ReferenceState_SWE;
#else
using ReferenceState = ReferenceState_Euler;
#endif
} // namespace pamc
