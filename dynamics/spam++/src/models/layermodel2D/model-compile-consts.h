
#ifndef _MODEL_COMPILE_CONSTS_H_
#define _MODEL_COMPILE_CONSTS_H_

//EVENTUALLY THIS NEEDS TO BE MORE CLEVERLY SET AT COMPILE TIME
//ALONG WITH A HOST OF OTHER THINGS IN common.h ie Q type, recon orders, etc.

#define _IDEAL_GAS_POTTEMP

#ifdef _SWE
uint constexpr ntracers = NTRACERS;
uint constexpr ntracers_fct = NTRACERS_FCT;
uint constexpr ndensity = 1 + ntracers;
uint constexpr ndensityfct = ntracers_fct;
#endif
#ifdef _TSWE
uint constexpr ntracers = NTRACERS;
uint constexpr ntracers_fct = NTRACERS_FCT;
uint constexpr ndensity = 2 + ntracers;
uint constexpr ndensityfct = ntracers_fct;
#endif
#ifdef _CE
uint constexpr ndensity = 2;
uint constexpr ndensityfct = 0;
uint constexpr ntracers = 0;
uint constexpr ntracers_fct = 0;
#endif
// MAYBE FURTHER SPECIALIZE TO NO ICE CASE?
#ifdef _MCERHO
uint constexpr ndensity = 2;
uint constexpr ndensityfct = 3;
uint constexpr ntracers = 0;
uint constexpr ntracers_fct = 0;
#endif
#ifdef _MCERHOD
uint constexpr ndensity = 2;
uint constexpr ndensityfct = 3;
uint constexpr ntracers = 0;
uint constexpr ntracers_fct = 0;
#endif

uint constexpr MAXTRACERS = 3;

// GOLO needs a better name here
enum class DATA_INIT { DOUBLEVORTEX, RB, MRB };
enum class TRACER_INIT { GAUSSIAN, SQUARE, DOUBLESQUARE };

class ModelParameters : public Parameters {
public:
  DATA_INIT data_init_cond;
  TRACER_INIT tracer_init_cond[MAXTRACERS];
  TRACER_INIT tracerFCT_init_cond[MAXTRACERS];
  
  real g;
};

#endif
