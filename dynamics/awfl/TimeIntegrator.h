
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



  //////////////////////////////////////////////////////////////////////////////////////////////////
  // Do allocations, initialize Tendencies object, and initialize Strang splitting direction switch
  // 
  // INPUTS
  //   dom: The Domain object
  //////////////////////////////////////////////////////////////////////////////////////////////////
  void initialize(Domain &dom);



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
  void stepForward(realArr &state, Domain &dom, Exchange &exch, Parallel const &par);



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
  void stepForwardADER(realArr &state, Domain &dom, Exchange &exch, Parallel const &par);



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
  void applyTendencies(realArr &state2, realArr const &tend, Domain const &dom);



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
  void applyTendencies(realArr &stateFinal, realArr const &state0, real dt, realArr const &tend, Domain const &dom);

};

#endif
