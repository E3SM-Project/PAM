
#ifndef _MODEL_COMPILE_CONSTS_H_
#define _MODEL_COMPILE_CONSTS_H_


enum class DATA_INIT { GAUSSIAN, VORTICES, SQUARE };
enum class WIND_INIT { UNIFORM_X, UNIFORM_Y, UNIFORM_Z, DEFORMATIONAL };

DATA_INIT data_init_cond = DATA_INIT::GAUSSIAN;
WIND_INIT wind_init_cond = WIND_INIT::UNIFORM_X;

uint constexpr nqdofs = 1;

#endif
