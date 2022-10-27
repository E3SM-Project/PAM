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
#endif
// ADD ANELASTIC + MOIST ANELASTIC

uint constexpr ndensity_nophysics = ndensity_dycore + ntracers_dycore;
uint constexpr ndensity = ndensity_dycore + ntracers_dycore + ntracers_physics;

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

uint constexpr nauxiliary = 28;

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

#define PHIVAR 11
#define PHIVERTVAR 12
#define EDGEFLUXVAR 13
#define VERTEDGEFLUXVAR 14
#define MFVAR 15

#define QXZ0VAR 16
#define QXZRECONVAR 17
#define QXZVERTRECONVAR 18
#define QXZEDGERECONVAR 19
#define QXZVERTEDGERECONVAR 20
#define QXZFLUXVAR 21
#define QXZVERTFLUXVAR 22

#define FXZ0VAR 23
#define CORIOLISXZRECONVAR 24
#define CORIOLISXZEDGERECONVAR 25
#define CORIOLISXZVERTRECONVAR 26
#define CORIOLISXZVERTEDGERECONVAR 27

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
  real entropicvar_diffusion_coeff;
  real velocity_diffusion_coeff;
};
