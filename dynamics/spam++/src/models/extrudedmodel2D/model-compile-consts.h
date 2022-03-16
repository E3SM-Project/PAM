
#ifndef _MODEL_COMPILE_CONSTS_H_
#define _MODEL_COMPILE_CONSTS_H_

// Number of Dimensions
uint constexpr ndims = 1;

// These are all set in CMakeLists.txt by reading model_consts.build
uint constexpr ntracers_nofct = _NTRACERS_NOFCT;
uint constexpr ntracers_fct = _NTRACERS_FCT;
uint constexpr ntracers = ntracers_nofct + ntracers_fct;

// Set dens and ntracers/ntracersfct sizes
#ifdef _SWE
uint constexpr ndensity_nofct = 1 + ntracers_nofct;
uint constexpr ndensity_fct = ntracers_fct;
#endif
#ifdef _TSWE
uint constexpr ndensity_nofct = 2 + ntracers_nofct;
uint constexpr ndensity_fct = ntracers_fct;
#endif
#ifdef _CE
uint constexpr ndensity_nofct = 2 + ntracers_nofct;
uint constexpr ndensity_fct = ntracers_fct;
#endif
#ifdef _MCE
uint constexpr ndensity_nofct = 2 + ntracers_nofct;
uint constexpr ndensity_fct = 3 + ntracers_fct;
#endif
#ifdef _AN
uint constexpr ndensity_nofct = 1 + ntracers_nofct;
uint constexpr ndensity_fct = ntracers_fct;
#endif
#ifdef _MAN
uint constexpr ndensity_nofct = 1 + ntracers_nofct;
uint constexpr ndensity_fct = 3 + ntracers_fct;
#endif

uint constexpr ndensity = ndensity_nofct + ndensity_fct;

uint constexpr MAXTRACERS = 6;

enum class DATA_INIT { RB, MRB, LRB, MLRB, DOUBLEVORTEX }; //DOUBLEVORTEX
enum class TRACER_INIT { GAUSSIAN, SQUARE, DOUBLESQUARE };

class ModelParameters : public Parameters {
public:
  DATA_INIT data_init_cond;
  TRACER_INIT tracer_init_cond[MAXTRACERS];
  
  real g;
  bool acoustic_balance;
};

#endif
