#pragma once

#include "Microphysics.h"
#include "SGS.h"

uint constexpr ntracers_dycore = 0;

//////////////////////////////////////////////////////////////////////////////

// forces reference state to be in perfect hydrostatic balance by subtracting
// the hydrostatic balance equation evaluated at the reference state in
// the velocity tendency
#if !defined(_SWE) && !defined(_TSWE)
#define FORCE_REFSTATE_HYDROSTATIC_BALANCE
#endif

// for debugging anelastic
#if defined _AN || defined _MAN
#define CHECK_ANELASTIC_CONSTRAINT
#endif

// Number of Dimensions
uint constexpr ndims = 1;

#include "params.h"

// Number of variables
// v, w, dens, densfct
uint constexpr nprognostic = 3;
#define VVAR 0
#define WVAR 1
#define DENSVAR 2

// hs, coriolis
uint constexpr nconstant = ndims > 1 ? 3 : 2;
#define HSVAR 0
#define CORIOLISHZVAR 1
#define CORIOLISXYVAR 2

// functional derivatives = F, FW, B, K, he, hew
// primal grid reconstruction stuff- U, W, dens0, edgerecon, recon,
// vertedgerecon, vertrecon fct stuff- Phi, Mf, edgeflux Q/W STUFF?

uint constexpr nauxiliary = ndims > 1 ? 37 : 30;

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

#define QHZVAR 14
#define QHZRECONVAR 15
#define QHZVERTRECONVAR 16
#define QHZEDGERECONVAR 17
#define QHZVERTEDGERECONVAR 18
#define QHZFLUXVAR 19
#define QHZVERTFLUXVAR 20

#define FHZVAR 21
#define CORIOLISHZRECONVAR 22
#define CORIOLISHZEDGERECONVAR 23
#define CORIOLISHZVERTRECONVAR 24
#define CORIOLISHZVERTEDGERECONVAR 25

#define F2VAR 26
#define FW2VAR 27
#define FTVAR 28
#define FTWVAR 29

// 3d auxiliary variables
#define QXYVAR 30
#define QXYRECONVAR 31
#define QXYEDGERECONVAR 32
#define FXYVAR 33
#define CORIOLISXYRECONVAR 34
#define CORIOLISXYEDGERECONVAR 35
#define FTXYVAR 36

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
  std::string tracerdataStr[ntracers_dycore + GPU_PAD];
  bool dycore_tracerpos[ntracers_dycore + GPU_PAD];
  bool acoustic_balance;
  bool uniform_vertical;
  real entropicvar_diffusion_coeff;
  real velocity_diffusion_coeff;

  realConst2d zint;
};
