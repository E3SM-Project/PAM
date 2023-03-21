#pragma once

#include "Microphysics.h"
#include "SGS.h"
#include "params.h"

uint constexpr ntracers_dycore = 0;

//////////////////////////////////////////////////////////////////////////////

// forces reference state to be in perfect hydrostatic balance by subtracting
// the hydrostatic balance equation evaluated at the reference state in
// the velocity tendency
#define FORCE_REFSTATE_HYDROSTATIC_BALANCE

// for debugging anelastic
#if defined _AN || defined _MAN
#define CHECK_ANELASTIC_CONSTRAINT
#endif

// Number of Dimensions
uint constexpr ndims = 1;

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
  std::string tracerdataStr[ntracers_dycore + GPU_PAD];
  bool dycore_tracerpos[ntracers_dycore + GPU_PAD];
  bool acoustic_balance;
  bool uniform_vertical;
  real entropicvar_diffusion_coeff;
  real velocity_diffusion_coeff;

  realConst2d zint;
};
