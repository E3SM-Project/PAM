
#include "TimeIntegrator.h"


//////////////////////////////////////////////////////////////////////////////////////////////////
// Do allocations, initialize Tendencies object, and initialize Strang splitting direction switch
// 
// INPUTS
//   dom: The Domain object
//////////////////////////////////////////////////////////////////////////////////////////////////
void TimeIntegrator::initialize(Domain &dom) {
  stateTmp = realArr("tend",numState,dom.nz+2*hs,dom.ny+2*hs,dom.nx+2*hs);
  tend     = realArr("tend",numState,dom.nz,dom.ny,dom.nx);
  tendencies.initialize(dom);
  dsSwitch = 1;
}



//////////////////////////////////////////////////////////////////////////////////////////////////
// Perform a time step
// 
// INPUTS 
//   state: The state vector: rho,u,v,w,theta. dims  (numState,nz+2*hs,ny+2*hs,nx+2*hs)
//   dom: The Domain class object   
//   par: The Parallel class object
// 
// OUTPUTS
//   state: The state vector unchanged except for halos 
//   dom: The Domain class object   
//   exch: The Exchange class object for halo and edge exchanges
//////////////////////////////////////////////////////////////////////////////////////////////////
void TimeIntegrator::stepForward(realArr &state, Domain &dom, Exchange &exch, Parallel const &par) {
  stepForwardADER(state, dom, exch, par);
}



//////////////////////////////////////////////////////////////////////////////////////////////////
// Perform a time step using Differential Transforms in time for single-stage, single-step
// high-order-accurate time integration. Uses a second-order-accurate alternating Strang splitting
// 
// INPUTS 
//   state: The state vector: rho,u,v,w,theta. dims  (numState,nz+2*hs,ny+2*hs,nx+2*hs)
//   dom: The Domain class object   
//   par: The Parallel class object
// 
// OUTPUTS
//   state: The state vector unchanged except for halos 
//   dom: The Domain class object   
//   exch: The Exchange class object for halo and edge exchanges
//////////////////////////////////////////////////////////////////////////////////////////////////
void TimeIntegrator::stepForwardADER(realArr &state, Domain &dom, Exchange &exch, Parallel const &par) {
  if (dsSwitch) {
    dsSwitch = 0;
    tendencies.compEulerTend_X(state, dom, exch, par, tend);
    applyTendencies( state , tend, dom);
    if (!dom.run2d) {
      tendencies.compEulerTend_Y(state, dom, exch, par, tend);
      applyTendencies( state , tend, dom);
    }
    tendencies.compEulerTend_Z(state, dom, tend);
    applyTendencies( state , tend, dom);
    tendencies.compEulerTend_S(state, dom, tend);
    applyTendencies( state , tend, dom);
  } else {
    dsSwitch = 1;
    tendencies.compEulerTend_S(state, dom, tend);
    applyTendencies( state , tend, dom);
    tendencies.compEulerTend_Z(state, dom, tend);
    applyTendencies( state , tend, dom);
    if (!dom.run2d) {
      tendencies.compEulerTend_Y(state, dom, exch, par, tend);
      applyTendencies( state , tend, dom);
    }
    tendencies.compEulerTend_X(state, dom, exch, par, tend);
    applyTendencies( state , tend, dom);
  }
}



//////////////////////////////////////////////////////////////////////////////////////////////////
// Add tendencies to the state
// 
// INPUTS 
//   tend: Time tendencies of the state vector, averaged over the time step
//   dom: The Domain class object   
// 
// OUTPUTS
//   state: The state vector
//////////////////////////////////////////////////////////////////////////////////////////////////
void TimeIntegrator::applyTendencies(realArr &state, realArr const &tend, Domain const &dom) {
  // for (int l=0; l<numState; l++) {
  //   for (int k=0; k<dom.nz; k++) {
  //     for (int j=0; j<dom.ny; j++) {
  //       for (int i=0; i<dom.nx; i++) {
  yakl::parallel_for( numState*dom.nz*dom.ny*dom.nx , YAKL_LAMBDA (int iGlob) {
    int l, k, j, i;
    yakl::unpackIndices( iGlob , numState,dom.nz,dom.ny,dom.nx , l,k,j,i );

    state(l,hs+k,hs+j,hs+i) += dom.dt * tend(l,k,j,i);
  });
}



//////////////////////////////////////////////////////////////////////////////////////////////////
// Perform a semi-discrete update for a low-storage three-stage Runge-Kutta method
// 
// INPUTS 
//   state0: The state at the beginning of the multi-stage time step
//   tend: Time tendencies of the state vector
//   dom: The Domain class object   
//   dt: The time step to apply to the tendencies for this stage
// 
// OUTPUTS
//   stateFinal: The state vector at the next stage
//////////////////////////////////////////////////////////////////////////////////////////////////
void TimeIntegrator::applyTendencies(realArr &stateFinal, realArr const &state0, real dt,
                                     realArr const &tend, Domain const &dom) {
  // for (int l=0; l<numState; l++) {
  //   for (int k=0; k<dom.nz; k++) {
  //     for (int j=0; j<dom.ny; j++) {
  //       for (int i=0; i<dom.nx; i++) {
  yakl::parallel_for( numState*dom.nz*dom.ny*dom.nx , YAKL_LAMBDA (int iGlob) {
    int l, k, j, i;
    yakl::unpackIndices( iGlob , numState,dom.nz,dom.ny,dom.nx , l,k,j,i );

    stateFinal(l,hs+k,hs+j,hs+i) = state0(l,hs+k,hs+j,hs+i) + dt * tend(l,k,j,i);
  });
}



