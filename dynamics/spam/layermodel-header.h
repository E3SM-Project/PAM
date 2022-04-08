#pragma once

// Number of Dimensions
uint constexpr ndims = 2;

uint constexpr ntracers_nofct = 3;
uint constexpr ntracers_fct = 3;
//uint constexpr MAXTRACERS = 6;

//enum class DATA_INIT { DOUBLEVORTEX, RB, MRB, LRB, MLRB, DENSITYWAVE };
//enum class TRACER_INIT { GAUSSIAN, SQUARE, DOUBLESQUARE };

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
// v, dens
uint constexpr nprognostic = 2;
#define VVAR 0
#define DENSVAR 1

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