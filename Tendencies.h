
#ifndef _TENDENCIES_H_
#define _TENDENCIES_H_

#include "const.h"
#include "Parallel.h"
#include "SArray.h"
#include "Riemann.h"
#include "Domain.h"
#include "Exchange.h"
#include "TendenciesThetaConsADER.h"
#include "TendenciesThetaConsSD.h"
#include "TendenciesThetaPrimADER.h"
#include "TendenciesThetaPrimSD.h"

class Tendencies {
  TendenciesThetaConsADER  tendenciesThetaConsADER;
  TendenciesThetaConsSD    tendenciesThetaConsSD  ;
  TendenciesThetaPrimADER  tendenciesThetaPrimADER;
  TendenciesThetaPrimSD    tendenciesThetaPrimSD  ;

  int const useTendThetaConsADER = 1;
  int const useTendThetaConsSD   = 2;
  int const useTendThetaPrimADER = 3;
  int const useTendThetaPrimSD   = 4;

  int useTend;

public :


  void initialize(Domain const &dom);


  void compEulerTend_X(realArr &state, Domain const &dom, Exchange &exch, Parallel const &par, realArr &tend);


  void compEulerTend_Y(realArr &state, Domain const &dom, Exchange &exch, Parallel const &par, realArr &tend);


  void compEulerTend_Z(realArr &state, Domain const &dom, Exchange &exch, Parallel const &par, realArr &tend);


  void compEulerTend_S(realArr &state, Domain const &dom, Exchange &exch, Parallel const &par, realArr &tend);


  void compStrakaTend(realArr &state, Domain const &dom, Exchange &exch, Parallel const &par, realArr &tend);

};

#endif

