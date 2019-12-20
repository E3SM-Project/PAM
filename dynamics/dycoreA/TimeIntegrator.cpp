
#include "TimeIntegrator.h"


void TimeIntegrator::initialize(Domain &dom) {
  tend = realArr("tend",numState,dom.nz,dom.ny,dom.nx);
  tendencies.initialize(dom);
  dsSwitch = 1;
}


void TimeIntegrator::stepForward(realArr &state, Domain &dom, Exchange &exch, Parallel const &par) {
  stepForwardADER(state, dom, exch, par);
}


void TimeIntegrator::stepForwardADER(realArr &state, Domain &dom, Exchange &exch, Parallel const &par) {
  if (dsSwitch) {
    dsSwitch = 0;
    tendencies.compEulerTend_X(state, dom, exch, par, tend);
    applyTendencies( state , tend, dom);
    if (!dom.run2d) {
      tendencies.compEulerTend_Y(state, dom, exch, par, tend);
      applyTendencies( state , tend, dom);
    }
    // dom.dt /= 2;
    // tendencies.compEulerTend_Z(state, dom, tend);
    // applyTendencies( state , tend, dom);
    // tendencies.compEulerTend_Z(state, dom, tend);
    // applyTendencies( state , tend, dom);
    // dom.dt *= 2;
    // tendencies.compEulerTend_S(state, dom, tend);
    // applyTendencies( state , tend, dom);
  } else {
    dsSwitch = 1;
    // tendencies.compEulerTend_S(state, dom, tend);
    // applyTendencies( state , tend, dom);
    // dom.dt /= 2;
    // tendencies.compEulerTend_Z(state, dom, tend);
    // applyTendencies( state , tend, dom);
    // tendencies.compEulerTend_Z(state, dom, tend);
    // applyTendencies( state , tend, dom);
    // dom.dt *= 2;
    if (!dom.run2d) {
      tendencies.compEulerTend_Y(state, dom, exch, par, tend);
      applyTendencies( state , tend, dom);
    }
    tendencies.compEulerTend_X(state, dom, exch, par, tend);
    applyTendencies( state , tend, dom);
  }
}



void TimeIntegrator::applyTendencies(realArr &state2, realArr const &tend, Domain const &dom) {
  // for (int l=0; l<numState; l++) {
  //   for (int k=0; k<dom.nz; k++) {
  //     for (int j=0; j<dom.ny; j++) {
  //       for (int i=0; i<dom.nx; i++) {
  yakl::parallel_for( numState,dom.nz,dom.ny,dom.nx , YAKL_LAMBDA (int l, int k, int j, int i) {
    state2(l,hs+k,hs+j,hs+i) += dom.dt * tend(l,k,j,i);
  });
}


