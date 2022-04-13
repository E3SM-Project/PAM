#pragma once

uint constexpr ntracers_nofct = 0;
uint constexpr ntracers_fct = 0;

//////////////////////////////////////////////////////////////////////////////

// Number of Dimensions
uint constexpr ndims = 1;

// Set dens and ntracers/ntracersfct sizes
uint constexpr ntracers = ntracers_nofct + ntracers_fct;

#ifdef _SWE
uint constexpr ndensity_nofct = 1 + ntracers_nofct;
uint constexpr ndensity_fct = ntracers_fct;
#elif _TSWE
uint constexpr ndensity_nofct = 2 + ntracers_nofct;
uint constexpr ndensity_fct = ntracers_fct;
#elif defined _CE || defined _CEp
uint constexpr ndensity_nofct = 2 + ntracers_nofct;
uint constexpr ndensity_fct = ntracers_fct;
#elif defined _MCErho || defined _MCErhop || defined _MCErhod || defined _MCErhodp
uint constexpr ndensity_nofct = 2 + ntracers_nofct;
uint constexpr ndensity_fct = 3 + ntracers_fct;
#endif
//ADD ANELASTIC + MOIST ANELASTIC

uint constexpr ndensity = ndensity_nofct + ndensity_fct;

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

//functional derivatives = F, FW, B, K, he, hew
//primal grid reconstruction stuff- U, W, dens0, edgerecon, recon, vertedgerecon, vertrecon
//fct stuff- Phi, Mf, edgeflux
//Q/W STUFF?

uint constexpr nauxiliary = 32;

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

// q, associated concentration 0-forms for den

uint constexpr ndiagnostic = 2;
#define DENSLDIAGVAR 0
#define QXZDIAGVAR 1

//track total densities, dens min/max, energy (total, K, P, I), PV, PE
uint constexpr nstats = 6;

#define DENSSTAT 0
#define DENSMINSTAT 1
#define DENSMAXSTAT 2
#define ESTAT 3
#define PVSTAT 4
#define PESTAT 5

//ADD ANELASTIC + MOIST ANELASTIC
// #if defined _CEp || defined _MCErhop || defined _MCErhodp
// #define PVAR 18
// #endif

class ModelParameters : public Parameters
{
public: 
  std::string initdataStr;
  std::string tracerdataStr[ntracers];
  bool acoustic_balance;
};
