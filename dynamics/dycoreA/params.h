
#ifndef _PARAMS_H_
#define _PARAMS_H_

#include "const.h"

int constexpr TIME_ADER     = 1;
int constexpr TIME_ADER_IMP = 2;
int constexpr TIME_SSPRK3   = 3;

int constexpr EQN_THETA_CONS = 1;
int constexpr EQN_THETA_PRIM = 2;

int constexpr DATA_INIT_COLLISION = 1;
int constexpr DATA_INIT_THERMAL   = 2;
int constexpr DATA_INIT_STRAKA    = 3;

extern real outFreq;
extern int  timeMethod;
extern int  strakaVisc;

#endif
