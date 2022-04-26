#pragma once

uint constexpr ntracers_dycore = 6;
uint constexpr ntracers_active = 6; //applies only for swe/tswe, determines how many of the tracers are dynamically active

//////////////////////////////////////////////////////////////////////////////

// Number of Dimensions
uint constexpr ndims = 2;

// Set dens sizes
//uint constexpr ntracers = ntracers_dycore + ntracers_physics;
//Tracers are stored as dycore_tracers, physics_tracers
//Dycore tracers never add mass

uint constexpr ntracers_physics = 0;

#ifdef _SWE
uint constexpr ndensity_dycore = 1;
#elif _TSWE
uint constexpr ndensity_dycore = 2;
#elif defined _CE || defined _CEp
uint constexpr ndensity_dycore = 2;
// Here we have assumed that the micro has at least defined the 3 key tracers: mass of (cloud) vapor/liquid/ice (the last might be zero, depending on the microphysics)
#elif defined _MCErho || defined _MCErhop || defined _MCErhod || defined _MCErhodp
uint constexpr ndensity_dycore = 2;
#endif
//ADD ANELASTIC + MOIST ANELASTIC

uint constexpr ndensity_nophysics = ndensity_dycore + ntracers_dycore;
uint constexpr ndensity = ndensity_dycore + ntracers_dycore;

// Number of variables
// v, dens
uint constexpr nprognostic = 2;
#define VVAR 0
#define DENSVAR 1

//This is needed in order for varset to compile for layermodel- it shouldn't ever be used so the -1 is okay (and will trigger an easy to catch bug if it is!)
#define WVAR -1

// hs, coriolis
uint constexpr nconstant = 2;
#define HSVAR 0
#define CORIOLISVAR 1

//functional derivatives = F, B, K, he
//dual grid reconstruction stuff- q0, f0, FT, qedgerecon, qrecon, coriolisedgercon, coriolisrecon
//primal grid reconstruction stuff- U, dens0, edgerecon, recon
//fct stuff- Phi, Mf, edgeflux

#if defined _CEp || defined _MCErhop || defined _MCErhodp
uint constexpr nauxiliary = 19;
#else
uint constexpr nauxiliary = 18;
#endif

#define FVAR 0
#define BVAR 1
#define KVAR 2
#define HEVAR 3
#define UVAR 4

#define DENS0VAR 5
#define DENSRECONVAR 6
#define DENSEDGERECONVAR 7

#define FTVAR 8
#define Q0VAR 9
#define F0VAR 10
#define QRECONVAR 11
#define QEDGERECONVAR 12
#define CORIOLISRECONVAR 13
#define CORIOLISEDGERECONVAR 14

#define PHIVAR 15
#define EDGEFLUXVAR 16
#define MFVAR 17

//ADD ANELASTIC + MOIST ANELASTIC
#if defined _CEp || defined _MCErhop || defined _MCErhodp
#define PVAR 18
#endif

// q, associated concentration 0-forms for den

uint constexpr ndiagnostic = 2;
#define QDIAGVAR 0
#define DENSLDIAGVAR 1

//track total densities, dens min/max, densfct min/max, energy (total, K, P, I), PV, PE,
uint constexpr nstats = 6;

#define DENSSTAT 0
#define DENSMINSTAT 1
#define DENSMAXSTAT 2
#define ESTAT 3
#define PVSTAT 4
#define PESTAT 5

class ModelParameters : public Parameters
{
public: 
  std::string initdataStr;
  std::string tracerdataStr[ntracers_dycore];
  bool dycore_tracerpos[ntracers_dycore];
  bool acoustic_balance;
};
