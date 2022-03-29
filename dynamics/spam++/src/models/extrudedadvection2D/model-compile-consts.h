
#ifndef _MODEL_COMPILE_CONSTS_H_
#define _MODEL_COMPILE_CONSTS_H_

// Number of Dimensions
uint constexpr ndims = 1;

uint constexpr ntnofctdofs = _NTRACERS_NOFCT;
uint constexpr ntfctdofs = _NTRACERS_FCT;
uint constexpr nQdofs = _NTRACERS_Q;
uint constexpr ntdofs = ntnofctdofs + ntfctdofs;

enum class DATA_INIT { GAUSSIAN, VORTICES, SQUARE, DOUBLESQUARE };
enum class WIND_INIT { UNIFORM_X, UNIFORM_Y, UNIFORM_XY, DEFORMATIONAL, DOUBLEVORTEX };

#define MAXTTRACERS 6
#define MAXQTRACERS 3

class ModelParameters : public Parameters {
public:
  DATA_INIT data_init_cond[MAXTTRACERS];
  DATA_INIT dataQ_init_cond[MAXQTRACERS];
  WIND_INIT wind_init_cond;
};

#endif
