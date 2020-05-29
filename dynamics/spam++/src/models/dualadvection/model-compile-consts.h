
#ifndef _MODEL_COMPILE_CONSTS_H_
#define _MODEL_COMPILE_CONSTS_H_

enum class DATA_INIT { GAUSSIAN, VORTICES, SQUARE };
enum class WIND_INIT { UNIFORM_X, UNIFORM_Y, UNIFORM_Z, UNIFORM_XY, DEFORMATIONAL };

uint constexpr nqdofs = NTRACERS;

class ModelParameters : public Parameters {
public:
  DATA_INIT data_init_cond[nqdofs];
  WIND_INIT wind_init_cond;
};

#endif
