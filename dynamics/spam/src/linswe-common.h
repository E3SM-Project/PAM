#pragma once

uint constexpr ntracers_dycore = 0;
uint constexpr ntracers_active = 0;
uint constexpr ntracers_physics = 0;
uint constexpr ndensity_dycore = 1;
uint constexpr ndensity_nophysics = ndensity_dycore + ntracers_dycore;
uint constexpr ndensity = ndensity_dycore + ntracers_dycore;


//////////////////////////////////////////////////////////////////////////////

// Number of Dimensions
uint constexpr ndims = 2;


// Number of variables
// v, h
uint constexpr nprognostic = 2;
#define VVAR 0
#define DENSVAR 1

//This is needed in order for varset to compile for layermodel- it shouldn't ever be used so the -1 is okay (and will trigger an easy to catch bug if it is!)
#define WVAR -1

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
  std::string tracerdataStr[ntracers_dycore];
  bool dycore_tracerpos[ntracers_dycore];
};
