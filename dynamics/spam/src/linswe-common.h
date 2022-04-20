#pragma once

uint constexpr ntracers_nofct = 0;
uint constexpr ntracers_fct = 0;

//////////////////////////////////////////////////////////////////////////////

// Number of Dimensions
uint constexpr ndims = 2;

// Set dens and ntracers/ntracersfct sizes
uint constexpr ntracers = ntracers_nofct + ntracers_fct;

uint constexpr ndensity_nofct = 1 + ntracers_nofct;
uint constexpr ndensity_fct = ntracers_fct;

uint constexpr ndensity = ndensity_nofct + ndensity_fct;

// Number of variables
// v, h
uint constexpr nprognostic = 2;
#define VVAR 0
#define HVAR 1

// G, HREF
uint constexpr nconstant = 4;
#define GRAVVAR 0
#define HREFVAR 1
#define HREFCVAR 2
#define HSVAR 3

// u = H * v, h0 = I * h
uint constexpr nauxiliary = 2;
#define UVAR 0
#define H0VAR 1

uint constexpr ndiagnostic = 4;
#define H0DIAGVAR 0
#define HE0DIAGVAR 1
#define HECDIAGVAR 2
#define VEDIAGVAR 3

// track total densities, dens min/max, densfct min/max, energy (total, K, P,
// I), PV, PE,
uint constexpr nstats = 1;

#define ESTAT 0

class ModelParameters : public Parameters {
public:
  std::string initdataStr;
};
