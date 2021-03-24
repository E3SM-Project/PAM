
#ifndef _MODEL_COMPILE_CONSTS_H_
#define _MODEL_COMPILE_CONSTS_H_


uint constexpr ntdofs = NTRACERS;
uint constexpr ntfctdofs = NTRACERS_FCT;
uint constexpr nQdofs = NTRACERS_Q;

enum class DATA_INIT { GAUSSIAN, VORTICES, SQUARE };
enum class WIND_INIT { UNIFORM_X, UNIFORM_Y, UNIFORM_Z, UNIFORM_XY, DEFORMATIONAL };

class ModelParameters : public Parameters {
public:
  DATA_INIT data_init_cond[ntdofs];
  DATA_INIT dataFCT_init_cond[ntfctdofs];
  DATA_INIT dataQ_init_cond[nQdofs];
  WIND_INIT wind_init_cond;
};

#endif
