#pragma once

// Number of Dimensions
uint constexpr ndims = 2;

uint constexpr ntnofctdofs = 3;
uint constexpr ntfctdofs = 3;
uint constexpr nQdofs = 3;
uint constexpr ntdofs = ntnofctdofs + ntfctdofs;

// Number of variables
uint constexpr nprognostic = 3;
uint constexpr nconstant = 3;
uint constexpr nauxiliary = 13;
uint constexpr ndiagnostic = 3;
uint constexpr nstats = 7;

#define TVAR 0
#define TFCTVAR 1
#define QVAR 2

#define VVAR 0
#define UVAR 1
#define UTVAR 2

#define T0VAR 0
#define TRECONVAR 1
#define TEDGERECONVAR 2
#define TFCT0VAR 3
#define TFCTRECONVAR 4
#define TFCTEDGERECONVAR 5
#define PHIVAR 6
#define EDGEFLUXVAR 7
#define MFVAR 8
#define Q0VAR 9
#define QRECONVAR 10
#define QEDGERECONVAR 11
#define QFLUXVAR 12

#define TDIAGVAR 0
#define TFCTDIAGVAR 1
#define QDIAGVAR 2

#define TMASSSTAT 0
#define TMINSTAT 1
#define TMAXSTAT 2
#define TFCTMASSSTAT 3
#define TFCTMINSTAT 4
#define TFCTMAXSTAT 5
#define QMASSSTAT 6