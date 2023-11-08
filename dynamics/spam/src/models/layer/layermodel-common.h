#pragma once

namespace pamc {

uint constexpr ntracers_dycore = 6;

//////////////////////////////////////////////////////////////////////////////

// Number of Dimensions
uint constexpr ndims = 2;
} // namespace pamc

#include "params.h"

namespace pamc {
// Number of variables
// v, dens
uint constexpr nprognostic = 2;
int constexpr VVAR = 0;
int constexpr DENSVAR = 1;

// This is needed in order for varset to compile for layermodel- it shouldn't
// ever be used so the -1 is okay (and will trigger an easy to catch bug if it
// is!)
int constexpr WVAR = -1;

// hs, coriolis
uint constexpr nconstant = 2;
int constexpr HSVAR = 0;
int constexpr CORIOLISVAR = 1;

int constexpr MASSDENSINDX = 0;

// functional derivatives = F, B, K, he
// dual grid reconstruction stuff- q0, f0, FT, qedgerecon, qrecon,
// coriolisedgercon, coriolisrecon primal grid reconstruction stuff- U, dens0,
// edgerecon, recon fct stuff- Phi, Mf, edgeflux

#if defined PAMC_CEp || defined PAMC_MCErhop || defined PAMC_MCErhodp
uint constexpr nauxiliary = 19;
#else
uint constexpr nauxiliary = 18;
#endif

int constexpr FVAR = 0;
int constexpr BVAR = 1;
int constexpr KVAR = 2;
int constexpr HEVAR = 3;
int constexpr UVAR = 4;

int constexpr DENS0VAR = 5;
int constexpr DENSRECONVAR = 6;
int constexpr DENSEDGERECONVAR = 7;

int constexpr FTVAR = 8;
int constexpr Q0VAR = 9;
int constexpr F0VAR = 10;
int constexpr QRECONVAR = 11;
int constexpr QEDGERECONVAR = 12;
int constexpr CORIOLISRECONVAR = 13;
int constexpr CORIOLISEDGERECONVAR = 14;

int constexpr EDGEFLUXVAR = 15;
int constexpr MFVAR = 16;
int constexpr F2VAR = 17;

// track total densities, dens min/max, densfct min/max, energy (total, K, P,
// I), PV, PE,
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
};
} // namespace pamc
