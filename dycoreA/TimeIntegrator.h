
#ifndef _TIMEINTEGRATOR_H_
#define _TIMEINTEGRATOR_H_

#include "const.h"
#include "Parallel.h"
#include "Domain.h"
#include "Tendencies.h"

class TimeIntegrator {

  realArr stateTmp;
  realArr tend;
  realArr tendTmp;
  Tendencies tendencies;
  int dsSwitch;

public :


  void initialize(Domain &dom);


  void stepForward(realArr &state, Domain &dom, Exchange &exch, Parallel const &par);


  void stepForwardADER(realArr &state, Domain &dom, Exchange &exch, Parallel const &par);


  void stepForwardSSPRK3(realArr &state, Domain const &dom, Exchange &exch, Parallel const &par);


  void applyTendencies(realArr &state2, real const c0, realArr const &state0,
                                              real const c1, realArr const &state1,
                                              real const ct, realArr const &tend, Domain const &dom);


  void appendTendencies(realArr &tend, realArr const &tendTmp, Domain const &dom);


  void applyHeatingCooling(realArr &state, Parallel const &par, Domain const &dom);

};

#endif
