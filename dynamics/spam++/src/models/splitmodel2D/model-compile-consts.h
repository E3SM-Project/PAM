
#ifndef _MODEL_COMPILE_CONSTS_H_
#define _MODEL_COMPILE_CONSTS_H_

// GOLO needs a better name here
enum class DATA_INIT { RB, HGW, NHGW };

class ModelParameters : public Parameters {
public:
  DATA_INIT data_init_cond;
};

#endif
