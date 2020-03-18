
#ifndef _MODEL_COMPILE_CONSTS_H_
#define _MODEL_COMPILE_CONSTS_H_

enum class DATA_INIT { GAUSSIAN, VORTICES, SQUARE };
enum class WIND_INIT { UNIFORM_X, UNIFORM_Y, UNIFORM_Z, UNIFORM_XY, DEFORMATIONAL };
enum class DIVERGENCE_TYPE { DF, CQ };

uint constexpr nqdofs = 1;
DIVERGENCE_TYPE constexpr div_type  = DIVERGENCE_TYPE::CQ;

class ModelParameters : public Parameters {
public:
  DATA_INIT data_init_cond;
  WIND_INIT wind_init_cond;
};

#endif
