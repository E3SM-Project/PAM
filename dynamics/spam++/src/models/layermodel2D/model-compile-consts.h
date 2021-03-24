
#ifndef _MODEL_COMPILE_CONSTS_H_
#define _MODEL_COMPILE_CONSTS_H_

// GOLO needs a better name here
enum class DATA_INIT { DOUBLEVORTEX, GOLO, TC2, TC5 };
enum class TRACER_INIT { GAUSSIAN, SQUARE, DOUBLESQUARE };

class ModelParameters : public Parameters {
public:
  DATA_INIT data_init_cond;
  TRACER_INIT tracer_init_cond[NTRACERS];
  TRACER_INIT tracerFCT_init_cond[NTRACERS_FCT];
};

#endif
