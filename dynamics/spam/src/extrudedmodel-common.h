#pragma once

#include "Microphysics.h"
#include "SGS.h"

uint constexpr ntracers_dycore = 0;
uint constexpr ntracers_active =
    0; // applies only for swe/tswe, determines how many of the tracers are
       // dynamically active

//////////////////////////////////////////////////////////////////////////////

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
uint constexpr nconstant = 12;
#define HSVAR 0
#define CORIOLISXZVAR 1
#define REFDENSVAR 2
#define REFDENS0VAR 3
#define REFNSQ0VAR 4
#define BDENS0VAR 5
#define REFDENSRECONVAR 6
#define REFDENSEDGERECONVAR 7
#define REFDENSVERTRECONVAR 8
#define REFDENSVERTEDGERECONVAR 9
#define REFHEVAR 10
#define REFHEWVAR 11

// functional derivatives = F, FW, B, K, he, hew
// primal grid reconstruction stuff- U, W, dens0, edgerecon, recon,
// vertedgerecon, vertrecon fct stuff- Phi, Mf, edgeflux Q/W STUFF?

uint constexpr nauxiliary = 34;

#define FVAR 0
#define BVAR 1
#define KVAR 2
#define HEVAR 3
#define UVAR 4
#define FWVAR 5
#define HEWVAR 6
#define UWVAR 7

#define DENS0VAR 8
#define DENSRECONVAR 9
#define DENSEDGERECONVAR 10
#define DENSVERTRECONVAR 11
#define DENSVERTEDGERECONVAR 12

#define PHIVAR 13
#define PHIVERTVAR 14
#define EDGEFLUXVAR 15
#define VERTEDGEFLUXVAR 16
#define MFVAR 17

#define QXZ0VAR 18
#define QXZRECONVAR 19
#define QXZVERTRECONVAR 20
#define QXZEDGERECONVAR 21
#define QXZVERTEDGERECONVAR 22
#define QXZFLUXVAR 23
#define QXZVERTFLUXVAR 24

#define FTVAR 25
#define FTWVAR 26

#define FXZ0VAR 27
#define CORIOLISXZRECONVAR 28
#define CORIOLISXZEDGERECONVAR 29
#define CORIOLISXZVERTRECONVAR 30
#define CORIOLISXZVERTEDGERECONVAR 31

#define FVAR2 32
#define FWVAR2 33

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
};
