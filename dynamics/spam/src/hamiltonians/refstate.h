#pragma once

#include "profiles.h"
#include "topology.h"

struct ReferenceState_SWE {
  real ref_height;

#ifdef _EXTRUDED
  Profile dens;
  Profile geop;
  Profile q_di;
  Profile q_pi;
  Profile rho_di;
  Profile rho_pi;
  Profile Nsq_pi;
  Profile v;
#endif
  bool is_initialized = false;

  template <class VS>
  void initialize(const Topology &primal_topology,
                  const Topology &dual_topology) {
#ifdef _EXTRUDED
    this->dens.initialize(dual_topology, "ref dens", 1, 1, VS::ndensity);
    this->geop.initialize(dual_topology, "ref geop", 1, 1, 1);
    this->rho_pi.initialize(primal_topology, "refrho_pi", 0, 0, 1);
    this->q_pi.initialize(primal_topology, "refq_pi", 0, 0, VS::ndensity);
    this->rho_di.initialize(dual_topology, "refrho_di", 0, 0, 1);
    this->q_di.initialize(dual_topology, "refq_di", 0, 0, VS::ndensity);
    this->Nsq_pi.initialize(primal_topology, "refNsq_pi", 0, 0, 1);
    this->v.initialize(primal_topology, "ref v", 1, 0, 1);
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
  Profile v;
#ifdef FORCE_REFSTATE_HYDROSTATIC_BALANCE
  Profile B;
#endif
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
    this->v.initialize(primal_topology, "ref v", 1, 0, 1);

#ifdef FORCE_REFSTATE_HYDROSTATIC_BALANCE
    this->B.initialize(dual_topology, "ref B", 1, 1, VS::ndensity_active);
#endif

    this->is_initialized = true;
  }
};

#if defined _SWE || defined _TSWE
using ReferenceState = ReferenceState_SWE;
#else
using ReferenceState = ReferenceState_Euler;
#endif
