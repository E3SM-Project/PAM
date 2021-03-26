
#ifndef _MODEL_COMPILE_CONSTS_H_
#define _MODEL_COMPILE_CONSTS_H_

//EVENTUALLY THIS NEEDS TO BE MORE CLEVERLY SET AT COMPILE TIME
#define _IDEAL_GAS_POTTEMP


// Number of Dimensions
uint constexpr ndims = 2;

uint constexpr ntracers = NTRACERS;
uint constexpr ntracers_fct = NTRACERS_FCT;

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

uint constexpr MAXTRACERS = 3;

enum class DATA_INIT { DOUBLEVORTEX, RB, LRB, MLRB };
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
