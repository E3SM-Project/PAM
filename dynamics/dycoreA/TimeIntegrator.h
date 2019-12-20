
#ifndef _TIMEINTEGRATOR_H_
#define _TIMEINTEGRATOR_H_

#include "const.h"
#include "Parallel.h"
#include "Domain.h"
#include "Tendencies.h"

class TimeIntegrator {

  realArr stateTmp;
  realArr tend;
  Tendencies tendencies;
  int dsSwitch;

public :


  void initialize(Domain &dom);


  void stepForward(realArr &state, Domain &dom, Exchange &exch, Parallel const &par);


  void stepForwardADER(realArr &state, Domain &dom, Exchange &exch, Parallel const &par);


  void applyTendencies(realArr &state2, realArr const &tend, Domain const &dom);


  void applyTendencies(realArr &stateFinal, realArr const &state0, real dt, realArr const &tend, Domain const &dom);


};

#endif
