#pragma once

#include "profiles.h"
#include "topology.h"

namespace pamc {

struct ReferenceState_SWE {
  real ref_height;

#ifdef PAMC_EXTRUDED
  Profile dens;
  Profile geop;
  Profile q_pi;
  Profile q_di;
  Profile rho_pi;
  Profile rho_di;
  Profile pres_pi;
  Profile pres_di;
  Profile Nsq_pi;
  Profile v;
  Profile B;
#endif
  bool is_initialized = false;

  template <class VS>
  void initialize(const Topology &primal_topology,
                  const Topology &dual_topology) {
#ifdef PAMC_EXTRUDED
    this->dens.initialize(dual_topology, "ref dens", 1, 1, VS::ndensity);
    this->geop.initialize(dual_topology, "ref geop", 1, 1, 1);
    this->q_pi.initialize(primal_topology, "refq_pi", 0, 0, VS::ndensity);
    this->q_di.initialize(dual_topology, "refq_di", 0, 0, VS::ndensity);
    this->v.initialize(primal_topology, "ref v", 1, 0, 1);
    this->rho_pi.initialize(primal_topology, "refrho_pi", 0, 0, 1);
    this->rho_di.initialize(dual_topology, "refrho_di", 0, 0, 1);
    this->pres_pi.initialize(primal_topology, "refp_pi", 0, 0, 1);
    this->pres_di.initialize(dual_topology, "refp_di", 0, 0, 1);
    this->Nsq_pi.initialize(primal_topology, "refNsq_pi", 0, 0, 1);
    this->B.initialize(dual_topology, "ref B", 1, 1, VS::ndensity_active);
#endif

    this->is_initialized = true;
  }
};

struct ReferenceState_Euler {
  Profile dens;
  Profile geop;
  Profile q_pi;
  Profile q_di;
  Profile rho_pi;
  Profile rho_di;
  Profile pres_pi;
  Profile pres_di;
  Profile Nsq_pi;
  Profile v;
  Profile B;
  bool is_initialized = false;

  template <class VS>
  void initialize(const Topology &primal_topology,
                  const Topology &dual_topology) {

    this->dens.initialize(dual_topology, "ref dens", 1, 1, VS::ndensity);
    this->geop.initialize(dual_topology, "ref geop", 1, 1, 1);
    this->q_pi.initialize(primal_topology, "refq_pi", 0, 0, VS::ndensity);
    this->q_di.initialize(dual_topology, "refq_di", 0, 0, VS::ndensity);
    this->rho_pi.initialize(primal_topology, "refrho_pi", 0, 0, 1);
    this->rho_di.initialize(dual_topology, "refrho_di", 0, 0, 1);
    this->pres_pi.initialize(primal_topology, "refp_pi", 0, 0, 1);
    this->pres_di.initialize(dual_topology, "refp_di", 0, 0, 1);
    this->Nsq_pi.initialize(primal_topology, "refNsq_pi", 0, 0, 1);
    this->v.initialize(primal_topology, "ref v", 1, 0, 1);
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
