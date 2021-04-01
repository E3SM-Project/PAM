
#ifndef _MODEL_COMPILE_CONSTS_H_
#define _MODEL_COMPILE_CONSTS_H_

// Number of Dimensions
uint constexpr ndims = 1;

// These are all set in CMakeLists.txt by reading model_consts.build
uint constexpr ntracers = _NTRACERS;
uint constexpr ntracers_fct = _NTRACERS_FCT;

// Set dens, densfct and ntracers/ntracersfct sizes
#ifdef _SWE
uint constexpr ndensity = 1 + ntracers;
uint constexpr ndensityfct = ntracers_fct;
#endif
#ifdef _TSWE
uint constexpr ndensity = 2 + ntracers;
uint constexpr ndensityfct = ntracers_fct;
#endif
#ifdef _CE
uint constexpr ndensity = 2 + ntracers;
uint constexpr ndensityfct = ntracers_fct;
#endif
#ifdef _MCE
uint constexpr ndensity = 2 + ntracers;
uint constexpr ndensityfct = 3 + ntracers_fct;
#endif
#ifdef _AN
uint constexpr ndensity = 1 + ntracers;
uint constexpr ndensityfct = ntracers_fct;
#endif
#ifdef _MAN
uint constexpr ndensity = 1 + ntracers;
uint constexpr ndensityfct = 3 + ntracers_fct;
#endif

uint constexpr MAXTRACERS = 3;

enum class DATA_INIT { DOUBLEVORTEX, RB, MRB, LRB, MLRB };
enum class TRACER_INIT { GAUSSIAN, SQUARE, DOUBLESQUARE };

class ModelParameters : public Parameters {
public:
  DATA_INIT data_init_cond;
  TRACER_INIT tracer_init_cond[MAXTRACERS];
  TRACER_INIT tracerFCT_init_cond[MAXTRACERS];
  
  real g;
  bool acoustic_balance;
};

#endif
