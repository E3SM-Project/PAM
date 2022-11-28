#pragma once

#include "Microphysics.h"
#include "SGS.h"

uint constexpr ntracers_dycore = 0;
uint constexpr ntracers_active =
    0; // applies only for swe/tswe, determines how many of the tracers are
       // dynamically active

//////////////////////////////////////////////////////////////////////////////

// forces reference state to be in perfect hydrostatic balance by subtracting
// the hydrostatic balance equation evaluated at the reference state in
// the velocity tendency
constexpr bool force_refstate_hydrostatic_balance = true;

// Number of Dimensions
uint constexpr ndims = 1;

// Set dens sizes
// uint constexpr ntracers = ntracers_dycore + ntracers_physics;
// Tracers are stored as dycore_tracers, physics_tracers
// Dycore tracers never add mass

#ifdef _SWE
uint constexpr ndensity_dycore = 1;
uint constexpr ntracers_physics = 0;
#elif _TSWE
uint constexpr ndensity_dycore = 2;
uint constexpr ntracers_physics = 0;
#elif defined _CE || defined _CEp
uint constexpr ndensity_dycore = 2;
uint constexpr ntracers_physics = 0;
// Here we have assumed that the micro has at least defined the 3 key tracers:
// mass of (cloud) vapor/liquid/ice (the last might be zero, depending on the
// microphysics)
#elif defined _MCErho || defined _MCErhop || defined _MCErhod ||               \
    defined _MCErhodp
uint constexpr ndensity_dycore = 2;
uint constexpr ntracers_physics =
    Microphysics::get_num_tracers() + SGS::get_num_tracers();
// the number of moist variables that PAM thermodynamic formulae assume
// (qd, qv, qc, qi)
uint constexpr nmoist = 4;
#endif
// ADD ANELASTIC + MOIST ANELASTIC

uint constexpr ndensity_nophysics = ndensity_dycore + ntracers_dycore;
uint constexpr ndensity = ndensity_dycore + ntracers_dycore + ntracers_physics;

#if defined _MCErho && defined _CONST_KAPPA_VIRPOTTEMP
bool constexpr tracers_decouple_from_dynamics = true;
uint constexpr ndensity_B = ndensity_dycore;
#else
bool constexpr tracers_decouple_from_dynamics = false;
uint constexpr ndensity_B = ndensity;
#endif

// Number of variables
// v, w, dens, densfct
uint constexpr nprognostic = 3;
#define VVAR 0
#define WVAR 1
#define DENSVAR 2

// hs, coriolis
uint constexpr nconstant = 2;
#define HSVAR 0
#define CORIOLISXZVAR 1

// functional derivatives = F, FW, B, K, he, hew
// primal grid reconstruction stuff- U, W, dens0, edgerecon, recon,
// vertedgerecon, vertrecon fct stuff- Phi, Mf, edgeflux Q/W STUFF?

uint constexpr nauxiliary = 26;

#define FVAR 0
#define BVAR 1
#define KVAR 2
#define UVAR 3
#define FWVAR 4
#define UWVAR 5

#define DENS0VAR 6
#define DENSRECONVAR 7
#define DENSEDGERECONVAR 8
#define DENSVERTRECONVAR 9
#define DENSVERTEDGERECONVAR 10

#define EDGEFLUXVAR 11
#define VERTEDGEFLUXVAR 12
#define MFVAR 13

#define QXZ0VAR 14
#define QXZRECONVAR 15
#define QXZVERTRECONVAR 16
#define QXZEDGERECONVAR 17
#define QXZVERTEDGERECONVAR 18
#define QXZFLUXVAR 19
#define QXZVERTFLUXVAR 20

#define FXZ0VAR 21
#define CORIOLISXZRECONVAR 22
#define CORIOLISXZEDGERECONVAR 23
#define CORIOLISXZVERTRECONVAR 24
#define CORIOLISXZVERTEDGERECONVAR 25

// track total densities, dens min/max, energy (total, K, P, I), PV, PE
uint constexpr nstats = 6;

#define DENSSTAT 0
#define DENSMINSTAT 1
#define DENSMAXSTAT 2
#define ESTAT 3
#define PVSTAT 4
#define PESTAT 5

// ADD ANELASTIC + MOIST ANELASTIC
//  #if defined _CEp || defined _MCErhop || defined _MCErhodp
//  #define PVAR 18
//  #endif

class ModelParameters : public Parameters {
public:
  std::string initdataStr;
  std::string tracerdataStr[ntracers_dycore];
  bool dycore_tracerpos[ntracers_dycore];
  bool acoustic_balance;
  bool uniform_vertical;
  real entropicvar_diffusion_coeff;
  real velocity_diffusion_coeff;

  realConst2d zint;
};
