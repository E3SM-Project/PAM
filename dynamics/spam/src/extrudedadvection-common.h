#pragma once

uint constexpr ntnofctdofs = 3;
uint constexpr ntfctdofs = 3;
uint constexpr nQdofs = 3;


////////////////////////////////////////


// Number of Dimensions
uint constexpr ndims = 1;

uint constexpr ntdofs = ntnofctdofs + ntfctdofs;

// Number of variables
// dens
uint constexpr nprognostic = 2;
#define DENSVAR 0
#define QXZVAR 1

uint constexpr nconstant = 7;
#define VVAR 0
#define WVAR 1
#define UVAR 2
#define UWVAR 3
#define VTVAR 4
#define WTVAR 5
#define D2VAR 6

//functional derivatives = F, FW, B, K, he, hew
//primal grid reconstruction stuff- U, W, dens0, edgerecon, recon, vertedgerecon, vertrecon
//fct stuff- Phi, Mf, edgeflux
//Q/W STUFF?

uint constexpr nauxiliary = 17;

#define DENS0VAR 0
#define DENSRECONVAR 1
#define DENSEDGERECONVAR 2
#define DENSVERTRECONVAR 3
#define DENSVERTEDGERECONVAR 4

#define PHIVAR 5
#define PHIVERTVAR 6
#define EDGEFLUXVAR 7
#define VERTEDGEFLUXVAR 8
#define MFVAR 9

#define QXZ0VAR 10
#define QXZRECONVAR 11
#define QXZVERTRECONVAR 12
#define QXZEDGERECONVAR 13
#define QXZVERTEDGERECONVAR 14
#define QXZFLUXVAR 15
#define QXZVERTFLUXVAR 16

// q, associated concentration 0-forms for den

uint constexpr ndiagnostic = 9;
#define DENSLDIAGVAR 0
#define QXZDIAGVAR 1

#define QXZRECONDIAGVAR 2
#define QXZVERTRECONDIAGVAR 3
#define QXZEDGERECONDIAGVAR 4
#define QXZVERTEDGERECONDIAGVAR 5
#define QXZFLUXDIAGVAR 6
#define QXZVERTFLUXDIAGVAR 7
#define D2DIAGVAR 8

//track total densities, dens min/max
uint constexpr nstats = 4;

#define DENSSTAT 0
#define DENSMINSTAT 1
#define DENSMAXSTAT 2
#define QXZSTAT 3

class ModelParameters : public Parameters
{
public: 
  std::string windInitStr;
  std::string TInitStr[ntdofs];
  std::string QInitStr[nQdofs];
};