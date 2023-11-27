#pragma once

#include "Microphysics.h"
#include "SGS.h"

namespace pamc {

uint constexpr ntracers_dycore = 0;

//////////////////////////////////////////////////////////////////////////////

// Number of Dimensions
uint constexpr ndims = PAMC_NDIMS;
} // namespace pamc

#include "params.h"

namespace pamc {

// Number of variables
// v, w, dens, densfct
uint constexpr nprognostic = 3;
int constexpr VVAR = 0;
int constexpr WVAR = 1;
int constexpr DENSVAR = 2;

// hs, coriolis
uint constexpr nconstant = ndims > 1 ? 3 : 2;
int constexpr HSVAR = 0;
int constexpr CORIOLISHZVAR = 1;
int constexpr CORIOLISXYVAR = 2;

// functional derivatives = F, FW, B, K, he, hew
// primal grid reconstruction stuff- U, W, dens0, edgerecon, recon,
// vertedgerecon, vertrecon fct stuff- Phi, Mf, edgeflux Q/W STUFF?

uint constexpr nauxiliary = ndims > 1 ? 39 : 32;

int constexpr FVAR = 0;
int constexpr BVAR = 1;
int constexpr KVAR = 2;
int constexpr UVAR = 3;
int constexpr FWVAR = 4;
int constexpr UWVAR = 5;

int constexpr DENS0VAR = 6;
int constexpr DENSRECONVAR = 7;
int constexpr DENSEDGERECONVAR = 8;
int constexpr DENSVERTRECONVAR = 9;
int constexpr DENSVERTEDGERECONVAR = 10;

int constexpr EDGEFLUXVAR = 11;
int constexpr VERTEDGEFLUXVAR = 12;
int constexpr MFVAR = 13;

int constexpr QHZVAR = 14;
int constexpr QHZRECONVAR = 15;
int constexpr QHZVERTRECONVAR = 16;
int constexpr QHZEDGERECONVAR = 17;
int constexpr QHZVERTEDGERECONVAR = 18;
int constexpr QHZFLUXVAR = 19;
int constexpr QHZVERTFLUXVAR = 20;

int constexpr FHZVAR = 21;
int constexpr CORIOLISHZRECONVAR = 22;
int constexpr CORIOLISHZEDGERECONVAR = 23;
int constexpr CORIOLISHZVERTRECONVAR = 24;
int constexpr CORIOLISHZVERTEDGERECONVAR = 25;

int constexpr F2VAR = 26;
int constexpr FW2VAR = 27;
int constexpr FTVAR = 28;
int constexpr FTWVAR = 29;

int constexpr FDIFFVAR = 30;
int constexpr FWDIFFVAR = 31;

// 3d auxiliary variables
int constexpr QXYVAR = 32;
int constexpr QXYRECONVAR = 33;
int constexpr QXYEDGERECONVAR = 34;
int constexpr FXYVAR = 35;
int constexpr CORIOLISXYRECONVAR = 36;
int constexpr CORIOLISXYEDGERECONVAR = 37;
int constexpr FTXYVAR = 38;

// track total densities, dens min/max, energy (total, K, P, I), PV, PE
uint constexpr nstats = 6;

int constexpr DENSSTAT = 0;
int constexpr DENSMINSTAT = 1;
int constexpr DENSMAXSTAT = 2;
int constexpr ESTAT = 3;
int constexpr PVSTAT = 4;
int constexpr PESTAT = 5;

class ModelParameters : public Parameters {
public:
  std::string init_data;
  std::string init_dycore_tracer[ntracers_dycore + GPU_PAD];
  bool dycore_tracer_pos[ntracers_dycore + GPU_PAD];
  bool acoustic_balance;
  bool uniform_vertical;
  real scalar_horiz_diffusion_coeff;
  real scalar_vert_diffusion_coeff;
  real velocity_vort_horiz_diffusion_coeff;
  real velocity_vort_vert_diffusion_coeff;
  real velocity_div_horiz_diffusion_coeff;
  real velocity_div_vert_diffusion_coeff;
  bool scalar_diffusion_subtract_refstate;
  bool velocity_diffusion_subtract_refstate;

  // forces reference state to be in perfect hydrostatic balance by subtracting
  // the hydrostatic balance equation evaluated at the reference state in
  // the velocity tendency
  bool force_refstate_hydrostatic_balance;
  bool check_anelastic_constraint;

  realConst2d zint;
};
} // namespace pamc
