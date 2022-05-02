#pragma once

uint constexpr ntnofctdofs = 3;
uint constexpr ntfctdofs = 3;
uint constexpr nQdofs = 3;


////////////////////////////////////////


// Number of Dimensions
uint constexpr ndims = 2;

uint constexpr ntdofs = ntnofctdofs + ntfctdofs;

// Number of variables
uint constexpr nprognostic = 2;
uint constexpr nconstant = 3;
uint constexpr nauxiliary = 10;
uint constexpr ndiagnostic = 2;
uint constexpr nstats = 4;

#define TVAR 0
#define QVAR 1

#define VVAR 0
#define UVAR 1
#define UTVAR 2

#define T0VAR 0
#define TRECONVAR 1
#define TEDGERECONVAR 2
#define PHIVAR 3
#define EDGEFLUXVAR 4
#define MFVAR 5
#define Q0VAR 6
#define QRECONVAR 7
#define QEDGERECONVAR 8
#define QFLUXVAR 9

#define TDIAGVAR 0
#define QDIAGVAR 1

#define TMASSSTAT 0
#define TMINSTAT 1
#define TMAXSTAT 2
#define QMASSSTAT 3

class ModelParameters : public Parameters
{
public: 
  std::string windInitStr;
  std::string TInitStr[ntdofs];
  std::string QInitStr[nQdofs];
};