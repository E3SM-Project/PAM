
#include "Tendencies.h"



void Tendencies::initialize(Domain const &dom) {
  TransformMatrices<real> trans;

  stateLimits = realArr("stateLimits",numState,2,dom.nz+1,dom.ny+1,dom.nx+1);
  flux_r      = realArr("flux_r"                ,dom.nz+1,dom.ny+1,dom.nx+1);
  flux_re     = realArr("flux_re"               ,dom.nz+1,dom.ny+1,dom.nx+1);

  // Setup the matrix to transform a stenicl (or coefs) into tord derivative GLL points
  SArray<real,ord,ord> s2c_ho;
  SArray<real,ord,ord> c2d_ho;
  trans.sten_to_coefs (s2c_ho);
  trans.coefs_to_deriv(c2d_ho);
  trans.coefs_to_gll_lower( to_gll );
  if (dom.doWeno) {
    to_derivX_gll = (to_gll * c2d_ho) / dom.dx;
    to_derivY_gll = (to_gll * c2d_ho) / dom.dy;
    to_derivZ_gll = (to_gll * c2d_ho) / dom.dz;
  } else {
    to_derivX_gll = (to_gll * c2d_ho * s2c_ho) / dom.dx;
    to_derivY_gll = (to_gll * c2d_ho * s2c_ho) / dom.dy;
    to_derivZ_gll = (to_gll * c2d_ho * s2c_ho) / dom.dz;
  }

  SArray<real,tord,tord> g2c, c2d, c2g;
  trans.gll_to_coefs  (g2c);
  trans.coefs_to_deriv(c2d);
  trans.coefs_to_gll  (c2g);
  aderDerivX = (c2g * c2d * g2c) / dom.dx;
  aderDerivY = (c2g * c2d * g2c) / dom.dy;
  aderDerivZ = (c2g * c2d * g2c) / dom.dz;

  // Setup the matrix to transform a stencil of ord cell averages into tord GLL points
  if (dom.doWeno) {
    trans.coefs_to_gll_lower( to_gll );
  } else {
    trans.sten_to_gll_lower( to_gll );
  }

  trans.weno_sten_to_coefs(wenoRecon);

  trans.get_gll_weights(gllWts);

  wenoSetIdealSigma(wenoIdl,wenoSigma);
}



///////////////////////////////////////////////////////////////////////////////////////////////////////
// Compute time tendencies for rho, u, v, w, and rho*e using ADER Differential Transform time
// stepping in the x-direction
// 
// INPUTS
//   state: The state vector: rho,u,v,w,rho*e. dims  (numState,nz+2*hs,ny+2*hs,nx+2*hs)
//   dom: The Domain class object   
//   par: The Parallel class object
// 
// OUTPUTS
//   state: The state vector unchanged except for halos 
//   exch: The Exchange class object for halo and edge exchanges
//   tend: The time tendencies of the state vector   (numstate,nz,ny,nx)
///////////////////////////////////////////////////////////////////////////////////////////////////////
void Tendencies::compEulerTend_X(realArr &state, Domain const &dom, Exchange &exch, Parallel const &par, 
                                 realArr &tend) {
  // Create local references of class arrays so that C++ Lambdas will properly capture them
  auto &stateLimits   = this->stateLimits  ;
  auto &flux_r        = this->flux_r       ;
  auto &flux_re       = this->flux_re      ;
  auto &gllWts        = this->gllWts       ;
  auto &to_gll        = this->to_gll       ;
  auto &to_derivX_gll = this->to_derivX_gll;
  auto &aderDerivX    = this->aderDerivX   ;
  auto &wenoRecon     = this->wenoRecon    ;
  auto &wenoIdl       = this->wenoIdl      ;
  auto &wenoSigma     = this->wenoSigma    ;

  // Exchange halos in the x-direction
  exch.haloInit      ();
  exch.haloPackN_x   (dom, state, numState);
  exch.haloExchange_x(dom, par);
  exch.haloUnpackN_x (dom, state, numState);

  ///////////////////////////////////////////////////////////////////////////////////////////////////////
  // For each cell:
  // (1) Reconstruct state values & spatial derivatives for rho, u, v, w, and p at tord GLL points
  //     across the cell
  // (2) Compute time derivatives for rho, u, v, w, and p at each GLL point using Differential
  //     Transforms in time. Save tendencies for u, v, and w as a by product of this operation
  // (3) Compute time averages of rho, u, v, w, and p, and utend, vtend, and wtend
  // (4) Compute local flux difference splitting tendencies for u, v, and w using spatial
  //     quadrature over tord GLL points
  // (5) Save the time-averaged state value estimates at the left and right of the cell for use in the
  //     upwind Riemann solver in the next step
  ///////////////////////////////////////////////////////////////////////////////////////////////////////
  // for (int k=0; k<dom.nz; k++) {
  //   for (int j=0; j<dom.ny; j++) {
  //     for (int i=0; i<dom.nx; i++) {
  yakl::parallel_for( dom.nz*dom.ny*dom.nx , YAKL_LAMBDA ( int const iGlob ) {
    int k, j, i;
    yakl::unpackIndices(iGlob,dom.nz,dom.ny,dom.nx,k,j,i);
    SArray<real,numState,tord,tord> stateDTs;  // GLL state DTs (var,time,space)
    SArray<real,numState,tord,tord> derivDTs;  // GLL deriv DTs (var,time,space)
    SArray<real         ,tord,tord> utend   ;  // DTs of u RHS      (time,space)
    SArray<real         ,tord,tord> vtend   ;  // DTs of v RHS      (time,space)
    SArray<real         ,tord,tord> wtend   ;  // DTs of w RHS      (time,space)

    // Compute tord GLL points of the fluid state and spatial derivative
    for (int l=0; l<numState; l++) {
      SArray<real,ord> stencil;
      SArray<real,tord> gllPts;
      // Store the stencil values
      // REMINDER: idT stands for idThermodynamic. It's meant to be general
      if (l != idT) {  // rho perturbation, u, v, and w
        for (int ii=0; ii<ord; ii++) {
          stencil(ii) = state(l,hs+k,hs+j,i+ii);
        }
      } else {         // pressure perturbation
        for (int ii=0; ii<ord; ii++) {
          real r  = state(idR,hs+k,hs+j,i+ii) + dom.hyDensCells  (hs+k);
          real u  = state(idU,hs+k,hs+j,i+ii);
          real v  = state(idV,hs+k,hs+j,i+ii);
          real w  = state(idW,hs+k,hs+j,i+ii);
          real re = state(idT,hs+k,hs+j,i+ii) + dom.hyEnergyCells(hs+k);
          real ke = 0.5_fp*r*(u*u+v*v+w*w);
          real p = RD/CV*(re-ke);
          stencil(ii) = p - dom.hyPressureCells(hs+k);
        }
      }

      // Reconstruct and store GLL points of the state values
      reconStencil(stencil, gllPts, dom.doWeno, wenoRecon, to_gll       , wenoIdl, wenoSigma);
      for (int ii=0; ii<tord; ii++) { stateDTs(l,0,ii) = gllPts(ii); }

      // Reconstruct and store GLL points of the state derivatives
      reconStencil(stencil, gllPts, dom.doWeno, wenoRecon, to_derivX_gll, wenoIdl, wenoSigma);
      for (int ii=0; ii<tord; ii++) { derivDTs(l,0,ii) = gllPts(ii); }
    }
    // Add hydrostasis to density and pressure to make them the full quantities
    for (int ii=0; ii<tord; ii++) {
      stateDTs(idR,0,ii) += dom.hyDensCells    (hs+k);
      stateDTs(idT,0,ii) += dom.hyPressureCells(hs+k);
    }

    // Compute tord-1 time derivatives of the state, state spatial derivatives, 
    // u RHS, v RHS, and w RHS using temporal Differential Transforms
    diffTransformEulerX( stateDTs, derivDTs, utend, vtend, wtend, aderDerivX );

    // Compute the time-average and store into the zeroth time index
    timeAvg( stateDTs , dom );
    timeAvg( utend    , dom );
    timeAvg( vtend    , dom );
    timeAvg( wtend    , dom );

    // Compute the local tendency contribution for high-order flux difference
    // splitting for wind via quadrature
    // tend_local = int( RHS , x , x_(i-1/2) , x_(i+1/2) )
    tend(idU,k,j,i) = 0;
    tend(idV,k,j,i) = 0;
    tend(idW,k,j,i) = 0;
    for (int ii=0; ii<tord; ii++) {
      tend(idU,k,j,i) += gllWts(ii) * utend(0,ii);
      tend(idV,k,j,i) += gllWts(ii) * vtend(0,ii);
      tend(idW,k,j,i) += gllWts(ii) * wtend(0,ii);
    }

    // Store the state vector in stateLimits to compute upwind forcing
    for (int l=0; l<numState; l++) {
      stateLimits(l,1,k,j,i  ) = stateDTs(l,0,0     ); // Store the left  cell edge state estimates
      stateLimits(l,0,k,j,i+1) = stateDTs(l,0,tord-1); // Store the Right cell edge state estimates
    }
  });

  //Reconcile the edge fluxes via MPI exchange
  exch.haloInit      ();
  exch.edgePackN_x   (dom, stateLimits, numState);
  exch.edgeExchange_x(dom, par);
  exch.edgeUnpackN_x (dom, stateLimits, numState);

  ///////////////////////////////////////////////////////////////////////////////////////////////////////
  // For each x-direction cell interface:
  // (1) Split the interface flux difference into leftward and rightward propagating waves (f-waves)
  //     (A+)*(q_R - q_L)   and   (A-)*(q_R - q_L)
  // (2) Compute the upwind flux vector for rho and rho*e using the upwind characteristic state
  ///////////////////////////////////////////////////////////////////////////////////////////////////////
  // for (int k=0; k<dom.nz; k++) {
  //   for (int j=0; j<dom.ny; j++) {
  //     for (int i=0; i<dom.nx+1; i++) {
  yakl::parallel_for( dom.nz*dom.ny*(dom.nx+1) , YAKL_LAMBDA (int const iGlob) {
    int k, j, i;
    yakl::unpackIndices(iGlob,dom.nz,dom.ny,dom.nx+1,k,j,i);

    ////////////////////////////////////////////////////////////////////////////
    // Split the flux difference for the wind updates
    ////////////////////////////////////////////////////////////////////////////
    // Compute the average state at the interface
    real r = 0.5_fp * ( stateLimits(idR,1,k,j,i) + stateLimits(idR,0,k,j,i)); // rho
    real u = 0.5_fp * ( stateLimits(idU,1,k,j,i) + stateLimits(idU,0,k,j,i)); // u
    real p = 0.5_fp * ( stateLimits(idT,1,k,j,i) + stateLimits(idT,0,k,j,i)); // p
    real cs2 = GAMMA*p/r;    // speed of sound squared
    real cs = sqrt(cs2);     // speed of sound
    // Compute the state jump across the interface
    real dr = stateLimits(idR,1,k,j,i) - stateLimits(idR,0,k,j,i);
    real du = stateLimits(idU,1,k,j,i) - stateLimits(idU,0,k,j,i);
    real dv = stateLimits(idV,1,k,j,i) - stateLimits(idV,0,k,j,i);
    real dw = stateLimits(idW,1,k,j,i) - stateLimits(idW,0,k,j,i);
    real dp = stateLimits(idT,1,k,j,i) - stateLimits(idT,0,k,j,i);
    // state at the left side of the interface
    real r1 = stateLimits(idR,0,k,j,i);
    real u1 = stateLimits(idU,0,k,j,i);
    real v1 = stateLimits(idV,0,k,j,i);
    real w1 = stateLimits(idW,0,k,j,i);
    real p1 = stateLimits(idT,0,k,j,i);
    // state at the right side of the interface
    real r2 = stateLimits(idR,1,k,j,i);
    real u2 = stateLimits(idU,1,k,j,i);
    real v2 = stateLimits(idV,1,k,j,i);
    real w2 = stateLimits(idW,1,k,j,i);
    real p2 = stateLimits(idT,1,k,j,i);
    // Block to force compiler to release df from the stack after the block
    {
      // Compute the product of the flux Jacobian and the state jump across the interface (A*dq)
      SArray<real,numState> df;
      df(0) = u*dr + r*du;
      df(1) = u*du + dp/r;
      df(2) = u*dv;
      df(3) = u*dw;
      df(4) = u*dp + GAMMA*p*du;
      // Zero out the stateLimits space for this spatial index for idU, idV, and idW
      stateLimits(idU,0,k,j,i) = 0;
      stateLimits(idV,0,k,j,i) = 0;
      stateLimits(idW,0,k,j,i) = 0;
      stateLimits(idU,1,k,j,i) = 0;
      stateLimits(idV,1,k,j,i) = 0;
      stateLimits(idW,1,k,j,i) = 0;
      // Wave 1 (u-cs): presumed always leftward  propagating (no shocks)
      stateLimits(idU,0,k,j,i) += (-cs/r) * ( -r/(2*cs)*df(1) + df(4)/(2*cs2) );
      // Wave 2 (u+cs): presumed always rightward propagating (no shocks)
      stateLimits(idU,1,k,j,i) += ( cs/r) * (  r/(2*cs)*df(1) + df(4)/(2*cs2) );
      // Wave 3 does only affects density, so it's ignored
      // Waves 4 and 5 (u): 
      // If u > zero, it's rightward propagating, otherwise leftward
      // No need to worry about zero wind speed becaue then the wave is zero anyway
      if (u > 0) {
        stateLimits(idV,1,k,j,i) += df(2); // wave 4 (u)
        stateLimits(idW,1,k,j,i) += df(3); // wave 5 (u)
      } else {
        stateLimits(idV,0,k,j,i) += df(2); // wave 4 (u)
        stateLimits(idW,0,k,j,i) += df(3); // wave 5 (u)
      }
    }

    ////////////////////////////////////////////////////////////////////////////
    // Compute the upwind Flux vector for mass and energy
    ////////////////////////////////////////////////////////////////////////////
    // We can re-use the r, u, p, cs2, and cs calculated earlier
    // Store upwind state based on wind velocity
    real ru, uu, vu, wu, pu;
    if (u > 0) {
      ru = r1;  uu = u1;  vu = v1;  wu = w1;  pu = p1;
    } else {
      ru = r2;  uu = u2;  vu = v2;  wu = w2;  pu = p2;
    }
    // The next two sections compute the upwind state vector defined by rho, u, v, w, p
    // First, compute upwind characteristic variables
    SArray<real,numState> chu;  // upwind characteristic variables
    chu(0) = -r/(2*cs)*u2 + p2/(2*cs2); // u-cs wave: assuming no shocks
    chu(1) =  r/(2*cs)*u1 + p1/(2*cs2); // u+cs wave: assuming no shocks
    chu(2) = ru - pu/cs2;
    chu(3) = vu;
    chu(4) = wu;
    // Next, compute the upwind state based on upwind characteristic variables
    ru =     chu(0)   +     chu(1)   + chu(2);
    uu = -cs*chu(0)/r +  cs*chu(1)/r;
    vu =                                       chu(3);
    wu =                                               chu(4);
    pu = cs2*chu(0)   + cs2*chu(1);
    // Finally, compute the upwind flux based on the upwind state
    real keu = 0.5_fp*ru*(uu*uu+vu*vu+wu*wu); // upwind kinetic energy
    real reu = pu*CV/RD + keu;                // upwind rho*e
    flux_r (k,j,i) = ru*uu;                   // upwind mass flux
    flux_re(k,j,i) = uu*reu + uu*pu;          // upwind energy flux
  });

  ///////////////////////////////////////////////////////////////////////////////////////////////////////
  // For each cell:
  // (1) Append the u, v, and w tendencies with the flux difference splitting waves entering the cell
  //     domain
  // (2) Compute the rho and rho*e tendencies using the upwind flux vectors
  ///////////////////////////////////////////////////////////////////////////////////////////////////////
  //   for (int k=0; k<dom.nz; k++) {
  //     for (int j=0; j<dom.ny; j++) {
  //       for (int i=0; i<dom.nx; i++) {
  yakl::parallel_for( dom.nz*dom.ny*dom.nx , YAKL_LAMBDA (int const iGlob) {
    int l, k, j, i;
    yakl::unpackIndices(iGlob,dom.nz,dom.ny,dom.nx,k,j,i);
    // Flux vector form update for mass and energy
    tend(idR,k,j,i)  = - ( flux_r (k,j,i+1) - flux_r (k,j,i) ) / dom.dx;  // mass tendency
    tend(idT,k,j,i)  = - ( flux_re(k,j,i+1) - flux_re(k,j,i) ) / dom.dx;  // energy tendency
    // Flux difference splitting form update for velocities
    tend(idU,k,j,i) += - ( stateLimits(idU,1,k,j,i) + stateLimits(idU,0,k,j,i+1) ) / dom.dx;  // u tendency
    tend(idV,k,j,i) += - ( stateLimits(idV,1,k,j,i) + stateLimits(idV,0,k,j,i+1) ) / dom.dx;  // v tendency
    tend(idW,k,j,i) += - ( stateLimits(idW,1,k,j,i) + stateLimits(idW,0,k,j,i+1) ) / dom.dx;  // w tendency
  });
}



///////////////////////////////////////////////////////////////////////////////////////////////////////
// Compute time tendencies for rho, u, v, w, and rho*e using ADER Differential Transform time
// stepping in the y-direction
// 
// INPUTS
//   state: The state vector: rho,u,v,w,rho*e. dims  (numState,nz+2*hs,ny+2*hs,nx+2*hs)
//   dom: The Domain class object   
//   par: The Parallel class object
// 
// OUTPUTS
//   state: The state vector unchanged except for halos 
//   exch: The Exchange class object for halo and edge exchanges
//   tend: The time tendencies of the state vector   (numstate,nz,ny,nx)
///////////////////////////////////////////////////////////////////////////////////////////////////////
void Tendencies::compEulerTend_Y(realArr &state, Domain const &dom, Exchange &exch, Parallel const &par, 
                                 realArr &tend) {
  // Create local references of class arrays so that C++ Lambdas will properly capture them
  auto &stateLimits   = this->stateLimits  ;
  auto &flux_r        = this->flux_r       ;
  auto &flux_re       = this->flux_re      ;
  auto &gllWts        = this->gllWts       ;
  auto &to_gll        = this->to_gll       ;
  auto &to_derivY_gll = this->to_derivY_gll;
  auto &aderDerivY    = this->aderDerivY   ;
  auto &wenoRecon     = this->wenoRecon    ;
  auto &wenoIdl       = this->wenoIdl      ;
  auto &wenoSigma     = this->wenoSigma    ;

  // Exchange halos in the y-direction
  exch.haloInit      ();
  exch.haloPackN_y   (dom, state, numState);
  exch.haloExchange_y(dom, par);
  exch.haloUnpackN_y (dom, state, numState);

  ///////////////////////////////////////////////////////////////////////////////////////////////////////
  // For each cell:
  // (1) Reconstruct state values & spatial derivatives for rho, u, v, w, and p at tord GLL points
  //     across the cell
  // (2) Compute time derivatives for rho, u, v, w, and p at each GLL point using Differential
  //     Transforms in time. Save tendencies for u, v, and w as a by product of this operation
  // (3) Compute time averages of rho, u, v, w, and p, and utend, vtend, and wtend
  // (4) Compute local flux difference splitting tendencies for u, v, and w using spatial
  //     quadrature over tord GLL points
  // (5) Save the time-averaged state value estimates at the left and right of the cell for use in the
  //     upwind Riemann solver in the next step
  ///////////////////////////////////////////////////////////////////////////////////////////////////////
  // for (int k=0; k<dom.nz; k++) {
  //   for (int j=0; j<dom.ny; j++) {
  //     for (int i=0; i<dom.nx; i++) {
  yakl::parallel_for( dom.nz*dom.ny*dom.nx , YAKL_LAMBDA ( int const iGlob ) {
    int k, j, i;
    yakl::unpackIndices(iGlob,dom.nz,dom.ny,dom.nx,k,j,i);
    SArray<real,numState,tord,tord> stateDTs;  // GLL state DTs (var,time,space)
    SArray<real,numState,tord,tord> derivDTs;  // GLL deriv DTs (var,time,space)
    SArray<real         ,tord,tord> utend   ;  // DTs of u RHS      (time,space)
    SArray<real         ,tord,tord> vtend   ;  // DTs of v RHS      (time,space)
    SArray<real         ,tord,tord> wtend   ;  // DTs of w RHS      (time,space)

    // Compute tord GLL points of the fluid state and spatial derivative
    for (int l=0; l<numState; l++) {
      SArray<real,ord> stencil;
      SArray<real,tord> gllPts;
      // Store the stencil values
      // We're reconstructing perturbation rho, but it doesn't affect the horizontal derivative
      if (l != idT) {  // rho perturbation, u, v, and w
        for (int ii=0; ii<ord; ii++) {
          stencil(ii) = state(l,hs+k,j+ii,hs+i);
        }
      } else {         // pressure perturbation
        for (int ii=0; ii<ord; ii++) {
          real r  = state(idR,hs+k,j+ii,hs+i) + dom.hyDensCells  (hs+k);
          real u  = state(idU,hs+k,j+ii,hs+i);
          real v  = state(idV,hs+k,j+ii,hs+i);
          real w  = state(idW,hs+k,j+ii,hs+i);
          real re = state(idT,hs+k,j+ii,hs+i) + dom.hyEnergyCells(hs+k);
          real ke = 0.5_fp*r*(u*u+v*v+w*w);
          real p = RD/CV*(re-ke);
          stencil(ii) = p - dom.hyPressureCells(hs+k);
        }
      }

      // Reconstruct and store GLL points of the state values
      reconStencil(stencil, gllPts, dom.doWeno, wenoRecon, to_gll       , wenoIdl, wenoSigma);
      for (int ii=0; ii<tord; ii++) { stateDTs(l,0,ii) = gllPts(ii); }

      // Reconstruct and store GLL points of the state derivatives
      reconStencil(stencil, gllPts, dom.doWeno, wenoRecon, to_derivX_gll, wenoIdl, wenoSigma);
      for (int ii=0; ii<tord; ii++) { derivDTs(l,0,ii) = gllPts(ii); }
    }
    // Add hydrostasis to density and pressure to make them the full quantities
    for (int ii=0; ii<tord; ii++) {
      stateDTs(idR,0,ii) += dom.hyDensCells    (hs+k);
      stateDTs(idT,0,ii) += dom.hyPressureCells(hs+k);
    }

    // Compute tord-1 time derivatives of the state, state spatial derivatives, 
    // u RHS, v RHS, and w RHS using temporal Differential Transforms
    diffTransformEulerY( stateDTs, derivDTs, utend, vtend, wtend, aderDerivY );

    // Compute the time-average and store into the zeroth time index
    timeAvg( stateDTs , dom );
    timeAvg( utend    , dom );
    timeAvg( vtend    , dom );
    timeAvg( wtend    , dom );

    // Compute the local tendency contribution for high-order flux difference
    // splitting for wind via quadrature
    // int( RHS , y , y_(j-1/2) , y_(j+1/2) )
    tend(idU,k,j,i) = 0;
    tend(idV,k,j,i) = 0;
    tend(idW,k,j,i) = 0;
    for (int ii=0; ii<tord; ii++) {
      tend(idU,k,j,i) += gllWts(ii) * utend(0,ii);
      tend(idV,k,j,i) += gllWts(ii) * vtend(0,ii);
      tend(idW,k,j,i) += gllWts(ii) * wtend(0,ii);
    }

    // Store the state vector in stateLimits to compute upwind forcing
    for (int l=0; l<numState; l++) {
      stateLimits(l,1,k,j  ,i) = stateDTs(l,0,0     ); // Store the left  cell edge state estimates
      stateLimits(l,0,k,j+1,i) = stateDTs(l,0,tord-1); // Store the Right cell edge state estimates
    }
  });

  //Reconcile the edge fluxes via MPI exchange.
  exch.haloInit      ();
  exch.edgePackN_y   (dom, stateLimits, numState);
  exch.edgeExchange_y(dom, par);
  exch.edgeUnpackN_y (dom, stateLimits, numState);

  ///////////////////////////////////////////////////////////////////////////////////////////////////////
  // For each y-direction cell interface:
  // (1) Split the interface flux difference into leftward and rightward propagating waves (f-waves)
  //     (A+)*(q_R - q_L)   and   (A-)*(q_R - q_L)
  // (2) Compute the upwind flux vector for rho and rho*e using the upwind characteristic state
  ///////////////////////////////////////////////////////////////////////////////////////////////////////
  // for (int k=0; k<dom.nz; k++) {
  //   for (int j=0; j<dom.ny+1; j++) {
  //     for (int i=0; i<dom.nx; i++) {
  yakl::parallel_for( dom.nz*(dom.ny+1)*dom.nx , YAKL_LAMBDA (int const iGlob) {
    int k, j, i;
    yakl::unpackIndices(iGlob,dom.nz,dom.ny+1,dom.nx,k,j,i);

    ////////////////////////////////////////////////////////////////////////////
    // Split the flux difference for the wind updates
    ////////////////////////////////////////////////////////////////////////////
    // Compute the average state at the interface
    real r = 0.5_fp * ( stateLimits(idR,1,k,j,i) + stateLimits(idR,0,k,j,i)); // rho
    real v = 0.5_fp * ( stateLimits(idV,1,k,j,i) + stateLimits(idV,0,k,j,i)); // v
    real p = 0.5_fp * ( stateLimits(idT,1,k,j,i) + stateLimits(idT,0,k,j,i)); // p
    real cs2 = GAMMA*p/r;    // speed of sound squared
    real cs = sqrt(cs2);     // speed of sound
    // Compute the state jump across the interface
    real dr = stateLimits(idR,1,k,j,i) - stateLimits(idR,0,k,j,i);
    real du = stateLimits(idU,1,k,j,i) - stateLimits(idU,0,k,j,i);
    real dv = stateLimits(idV,1,k,j,i) - stateLimits(idV,0,k,j,i);
    real dw = stateLimits(idW,1,k,j,i) - stateLimits(idW,0,k,j,i);
    real dp = stateLimits(idT,1,k,j,i) - stateLimits(idT,0,k,j,i);
    // state at the left side of the interface
    real r1 = stateLimits(idR,0,k,j,i);
    real u1 = stateLimits(idU,0,k,j,i);
    real v1 = stateLimits(idV,0,k,j,i);
    real w1 = stateLimits(idW,0,k,j,i);
    real p1 = stateLimits(idT,0,k,j,i);
    // state at the right side of the interface
    real r2 = stateLimits(idR,1,k,j,i);
    real u2 = stateLimits(idU,1,k,j,i);
    real v2 = stateLimits(idV,1,k,j,i);
    real w2 = stateLimits(idW,1,k,j,i);
    real p2 = stateLimits(idT,1,k,j,i);
    // Block to force compiler to release df from the stack after the block
    {
      // Compute the product of the flux Jacobian and the state jump across the interface (A*dq)
      SArray<real,numState> df;
      df(0) = v*dr + r*dv;
      df(1) = v*du;
      df(2) = v*dv + dp/r;
      df(3) = v*dw;
      df(4) = v*dp + GAMMA*p*dv;
      // Zero out the stateLimits space for this spatial index for idU, idV, and idW
      stateLimits(idU,0,k,j,i) = 0;
      stateLimits(idV,0,k,j,i) = 0;
      stateLimits(idW,0,k,j,i) = 0;
      stateLimits(idU,1,k,j,i) = 0;
      stateLimits(idV,1,k,j,i) = 0;
      stateLimits(idW,1,k,j,i) = 0;
      // Wave 1 (v-cs): presumed always leftward  propagating (no shocks)
      stateLimits(idV,0,k,j,i) += (-cs/r) * ( -r/(2*cs)*df(2) + df(4)/(2*cs2) );
      // Wave 2 (v+cs): presumed always rightward propagating (no shocks)
      stateLimits(idV,1,k,j,i) += ( cs/r) * (  r/(2*cs)*df(2) + df(4)/(2*cs2) );
      // Wave 3 does only affects density, so it's ignored
      // Waves 4 and 5 (v): 
      // If v > zero, it's rightward propagating, otherwise leftward
      // No need to worry about zero wind speed becaue then the wave is zero anyway
      if (v > 0) {
        stateLimits(idU,1,k,j,i) += df(1);  // wave 4 (v)
        stateLimits(idW,1,k,j,i) += df(3);  // wave 5 (v)
      } else {
        stateLimits(idU,0,k,j,i) += df(1);  // wave 4 (v)
        stateLimits(idW,0,k,j,i) += df(3);  // wave 5 (v)
      }
    }

    ////////////////////////////////////////////////////////////////////////////
    // Compute the upwind Flux vector for mass and energy
    ////////////////////////////////////////////////////////////////////////////
    // We can re-use the r, v, p, cs2, and cs calculated earlier
    // Store upwind state based on wind velocity
    real ru, uu, vu, wu, pu;
    if (v > 0) {
      ru = r1;  uu = u1;  vu = v1;  wu = w1;  pu = p1;
    } else {
      ru = r2;  uu = u2;  vu = v2;  wu = w2;  pu = p2;
    }
    // The next two sections compute the upwind state vector defined by rho, u, v, w, p
    // First, compute upwind characteristic variables
    SArray<real,numState> chu;  // upwind characteristic variables
    chu(0) = -r/(2*cs)*v2 + p2/(2*cs2); // v-cs wave: assuming no shocks
    chu(1) =  r/(2*cs)*v1 + p1/(2*cs2); // v+cs wave: assuming no shocks
    chu(2) = ru - pu/cs2;
    chu(3) = uu;
    chu(4) = wu;
    // Next, compute the upwind state based on upwind characteristic variables
    ru =     chu(0)   +     chu(1)   + chu(2);
    uu =                                       chu(3);
    vu = -cs*chu(0)/r +  cs*chu(1)/r;
    wu =                                               chu(4);
    pu = cs2*chu(0)   + cs2*chu(1);
    // Finally, compute the upwind flux based on the upwind state
    real keu = 0.5_fp*ru*(uu*uu+vu*vu+wu*wu); // upwind kinetic energy
    real reu = pu*CV/RD + keu;                // upwind rho*e
    flux_r (k,j,i) = ru*vu;                   // upwind mass flux
    flux_re(k,j,i) = vu*reu + vu*pu;          // upwind energy flux
  });

  ///////////////////////////////////////////////////////////////////////////////////////////////////////
  // For each cell:
  // (1) Append the u, v, and w tendencies with the flux difference splitting waves entering the cell
  //     domain
  // (2) Compute the rho and rho*e tendencies using the upwind flux vectors
  ///////////////////////////////////////////////////////////////////////////////////////////////////////
  //   for (int k=0; k<dom.nz; k++) {
  //     for (int j=0; j<dom.ny; j++) {
  //       for (int i=0; i<dom.nx; i++) {
  yakl::parallel_for( dom.nz*dom.ny*dom.nx , YAKL_LAMBDA (int const iGlob) {
    int l, k, j, i;
    yakl::unpackIndices(iGlob,dom.nz,dom.ny,dom.nx,k,j,i);
    // Flux vector form update for mass and energy
    tend(idR,k,j,i)  = - ( flux_r (k,j+1,i) - flux_r (k,j,i) ) / dom.dy;  // mass tendency
    tend(idT,k,j,i)  = - ( flux_re(k,j+1,i) - flux_re(k,j,i) ) / dom.dy;  // energy tendency
    // Flux difference splitting form update for velocities
    tend(idU,k,j,i) += - ( stateLimits(idU,1,k,j,i) + stateLimits(idU,0,k,j+1,i) ) / dom.dy;  // u tendency
    tend(idV,k,j,i) += - ( stateLimits(idV,1,k,j,i) + stateLimits(idV,0,k,j+1,i) ) / dom.dy;  // v tendency
    tend(idW,k,j,i) += - ( stateLimits(idW,1,k,j,i) + stateLimits(idW,0,k,j+1,i) ) / dom.dy;  // w tendency
  });
}



///////////////////////////////////////////////////////////////////////////////////////////////////////
// Compute time tendencies for rho, u, v, w, and rho*e using ADER Differential Transform time
// stepping in the z-direction
// 
// INPUTS
//   state: The state vector: rho,u,v,w,rho*e. dims  (numState,nz+2*hs,ny+2*hs,nx+2*hs)
//   dom: The Domain class object   
// 
// OUTPUTS
//   state: The state vector unchanged except for halos 
//   tend: The time tendencies of the state vector   (numstate,nz,ny,nx)
///////////////////////////////////////////////////////////////////////////////////////////////////////
void Tendencies::compEulerTend_Z(realArr &state, Domain const &dom, realArr &tend) {
  // Create local references of class arrays so that C++ Lambdas will properly capture them
  auto &stateLimits   = this->stateLimits  ;
  auto &flux_r        = this->flux_r       ;
  auto &flux_re       = this->flux_re      ;
  auto &gllWts        = this->gllWts       ;
  auto &to_gll        = this->to_gll       ;
  auto &to_derivZ_gll = this->to_derivZ_gll;
  auto &aderDerivZ    = this->aderDerivZ   ;
  auto &wenoRecon     = this->wenoRecon    ;
  auto &wenoIdl       = this->wenoIdl      ;
  auto &wenoSigma     = this->wenoSigma    ;

  // Apply boundaries in the z-direction
  stateBoundariesZ( state , dom );

  ///////////////////////////////////////////////////////////////////////////////////////////////////////
  // For each cell:
  // (1) Reconstruct state values & spatial derivatives for rho, u, v, w, and p at tord GLL points
  //     across the cell
  // (2) Compute time derivatives for rho, u, v, w, and p at each GLL point using Differential
  //     Transforms in time. Save tendencies for u, v, and w as a by product of this operation
  // (3) Compute time averages of rho, u, v, w, and p, and utend, vtend, and wtend
  // (4) Compute local flux difference splitting tendencies for u, v, and w using spatial
  //     quadrature over tord GLL points
  // (5) Save the time-averaged state value estimates at the left and right of the cell for use in the
  //     upwind Riemann solver in the next step
  ///////////////////////////////////////////////////////////////////////////////////////////////////////
  // for (int k=0; k<dom.nz; k++) {
  //   for (int j=0; j<dom.ny; j++) {
  //     for (int i=0; i<dom.nx; i++) {
  yakl::parallel_for( dom.nz*dom.ny*dom.nx , YAKL_LAMBDA ( int const iGlob ) {
    int k, j, i;
    yakl::unpackIndices(iGlob,dom.nz,dom.ny,dom.nx,k,j,i);
    SArray<real,numState,tord,tord> stateDTs;  // GLL state DTs (var,time,space)
    SArray<real,numState,tord,tord> derivDTs;  // GLL deriv DTs (var,time,space)
    SArray<real         ,tord,tord> utend   ;  // DTs of u RHS      (time,space)
    SArray<real         ,tord,tord> vtend   ;  // DTs of v RHS      (time,space)
    SArray<real         ,tord,tord> wtend   ;  // DTs of w RHS      (time,space)

    // Compute tord GLL points of the fluid state and spatial derivative
    for (int l=0; l<numState; l++) {
      SArray<real,ord> stencil;
      SArray<real,tord> gllPts;
      // Store the stencil values
      if (l != idT) {  // rho perturbation, u, v, and w
        for (int ii=0; ii<ord; ii++) {
          stencil(ii) = state(l,k+ii,hs+j,hs+i);
        }
      } else {         // pressure perturbation
        for (int ii=0; ii<ord; ii++) {
          real r  = state(idR,k+ii,hs+j,hs+i) + dom.hyDensCells  (k+ii);
          real u  = state(idU,k+ii,hs+j,hs+i);
          real v  = state(idV,k+ii,hs+j,hs+i);
          real w  = state(idW,k+ii,hs+j,hs+i);
          real re = state(idT,k+ii,hs+j,hs+i) + dom.hyEnergyCells(k+ii);
          real ke = 0.5_fp*r*(u*u+v*v+w*w);
          real p = RD/CV*(re-ke);
          stencil(ii) = p - dom.hyPressureCells(k+ii);
        }
      }

      // Reconstruct and store GLL points of the state values
      reconStencil(stencil, gllPts, dom.doWeno, wenoRecon, to_gll       , wenoIdl, wenoSigma);
      for (int ii=0; ii<tord; ii++) { stateDTs(l,0,ii) = gllPts(ii); }

      // Reconstruct and store GLL points of the state derivatives
      reconStencil(stencil, gllPts, dom.doWeno, wenoRecon, to_derivX_gll, wenoIdl, wenoSigma);
      for (int ii=0; ii<tord; ii++) { derivDTs(l,0,ii) = gllPts(ii); }
    }
    // Add hydrostasis to density and pressure to make them the full quantities
    for (int ii=0; ii<tord; ii++) {
      stateDTs(idR,0,ii) += dom.hyDensGLL    (k,ii);
      stateDTs(idT,0,ii) += dom.hyPressureGLL(k,ii);
    }

    // Compute tord-1 time derivatives of the state, state spatial derivatives, 
    // u RHS, v RHS, and w RHS using temporal Differential Transforms
    SArray<real,tord> dph;
    for (int ii=0; ii<tord; ii++) { dph(ii) = dom.hyPressureDerivGLL(k,ii); }
    diffTransformEulerZ( stateDTs, derivDTs, utend, vtend, wtend, dph, aderDerivZ );

    // Compute the time-average and store into the zeroth time index
    timeAvg( stateDTs , dom );
    timeAvg( utend    , dom );
    timeAvg( vtend    , dom );
    timeAvg( wtend    , dom );

    // Compute the local tendency contribution for high-order flux difference
    // splitting for wind via quadrature
    // int( RHS , z , z_(k-1/2) , z_(k+1/2) )
    tend(idU,k,j,i) = 0;
    tend(idV,k,j,i) = 0;
    tend(idW,k,j,i) = 0;
    for (int ii=0; ii<tord; ii++) {
      tend(idU,k,j,i) += gllWts(ii) * utend(0,ii);
      tend(idV,k,j,i) += gllWts(ii) * vtend(0,ii);
      tend(idW,k,j,i) += gllWts(ii) * wtend(0,ii);
    }

    // Store the state vector in fwaves to compute fwaves from cell-interface state jumps
    for (int l=0; l<numState; l++) {
      stateLimits(l,1,k  ,j,i) = stateDTs(l,0,0     ); // Store the left  cell edge state estimates
      stateLimits(l,0,k+1,j,i) = stateDTs(l,0,tord-1); // Store the Right cell edge state estimates
    }
  });

  // Apply boundaries for fluxes in the z-direction
  edgeBoundariesZ( stateLimits , dom );

  ///////////////////////////////////////////////////////////////////////////////////////////////////////
  // For each z-direction cell interface:
  // (1) Split the interface flux difference into leftward and rightward propagating waves (f-waves)
  //     (A+)*(q_R - q_L)   and   (A-)*(q_R - q_L)
  // (2) Compute the upwind flux vector for rho and rho*e using the upwind characteristic state
  ///////////////////////////////////////////////////////////////////////////////////////////////////////
  // for (int k=0; k<dom.nz+1; k++) {
  //   for (int j=0; j<dom.ny; j++) {
  //     for (int i=0; i<dom.nx; i++) {
  yakl::parallel_for( (dom.nz+1)*dom.ny*dom.nx , YAKL_LAMBDA (int const iGlob) {
    int k, j, i;
    yakl::unpackIndices(iGlob,dom.nz+1,dom.ny,dom.nx,k,j,i);

    ////////////////////////////////////////////////////////////////////////////
    // Split the flux difference for the wind updates
    ////////////////////////////////////////////////////////////////////////////
    // Compute the average state at the interface
    real r = 0.5_fp * ( stateLimits(idR,1,k,j,i) + stateLimits(idR,0,k,j,i)); // rho
    real w = 0.5_fp * ( stateLimits(idW,1,k,j,i) + stateLimits(idW,0,k,j,i)); // w
    real p = 0.5_fp * ( stateLimits(idT,1,k,j,i) + stateLimits(idT,0,k,j,i)); // p
    real cs2 = GAMMA*p/r;    // speed of sound squared
    real cs = sqrt(cs2);     // speed of sound
    // Compute the state jump across the interface
    real dr = stateLimits(idR,1,k,j,i) - stateLimits(idR,0,k,j,i);
    real du = stateLimits(idU,1,k,j,i) - stateLimits(idU,0,k,j,i);
    real dv = stateLimits(idV,1,k,j,i) - stateLimits(idV,0,k,j,i);
    real dw = stateLimits(idW,1,k,j,i) - stateLimits(idW,0,k,j,i);
    real dp = stateLimits(idT,1,k,j,i) - stateLimits(idT,0,k,j,i);
    // state at the left side of the interface
    real r1 = stateLimits(idR,0,k,j,i);
    real u1 = stateLimits(idU,0,k,j,i);
    real v1 = stateLimits(idV,0,k,j,i);
    real w1 = stateLimits(idW,0,k,j,i);
    real p1 = stateLimits(idT,0,k,j,i);
    // state at the right side of the interface
    real r2 = stateLimits(idR,1,k,j,i);
    real u2 = stateLimits(idU,1,k,j,i);
    real v2 = stateLimits(idV,1,k,j,i);
    real w2 = stateLimits(idW,1,k,j,i);
    real p2 = stateLimits(idT,1,k,j,i);
    // Block to force compiler to release df from the stack after the block
    {
      // Compute the product of the flux Jacobian and the state jump across the interface (A*dq)
      SArray<real,numState> df;
      df(0) = w*dr + r*dw;
      df(1) = w*du;
      df(2) = w*dv;
      df(3) = w*dw + dp/r;
      df(4) = w*dp + GAMMA*p*dw;
      // Zero out the stateLimits space for this spatial index for idU, idV, and idW
      stateLimits(idU,0,k,j,i) = 0;
      stateLimits(idV,0,k,j,i) = 0;
      stateLimits(idW,0,k,j,i) = 0;
      stateLimits(idU,1,k,j,i) = 0;
      stateLimits(idV,1,k,j,i) = 0;
      stateLimits(idW,1,k,j,i) = 0;
      // Wave 1 (w-cs): presumed always leftward  propagating (no shocks)
      stateLimits(idW,0,k,j,i) += (-cs/r) * ( -r/(2*cs)*df(3) + df(4)/(2*cs2) );
      // Wave 2 (w+cs): presumed always rightward propagating (no shocks)
      stateLimits(idW,1,k,j,i) += ( cs/r) * (  r/(2*cs)*df(3) + df(4)/(2*cs2) );
      // Wave 3 does only affects density, so it's ignored
      // Waves 4 and 5 (w): 
      // If w > zero, it's rightward propagating, otherwise leftward
      // No need to worry about zero wind speed becaue then the wave is zero anyway
      if (w > 0) {
        stateLimits(idU,1,k,j,i) += df(1);  // wave 4 (w)
        stateLimits(idV,1,k,j,i) += df(2);  // wave 5 (w)
      } else {
        stateLimits(idU,0,k,j,i) += df(1);  // wave 4 (w)
        stateLimits(idV,0,k,j,i) += df(2);  // wave 5 (w)
      }
    }

    ////////////////////////////////////////////////////////////////////////////
    // Compute the upwind Flux vector for mass and energy
    ////////////////////////////////////////////////////////////////////////////
    // We can re-use the r, u, p, cs2, and cs calculated earlier
    // Store upwind state based on wind velocity
    real ru, uu, vu, wu, pu;
    if (w > 0) {
      ru = r1;  uu = u1;  vu = v1;  wu = w1;  pu = p1;
    } else {
      ru = r2;  uu = u2;  vu = v2;  wu = w2;  pu = p2;
    }
    // The next two sections compute the upwind state vector defined by rho, u, v, w, p
    // First, compute upwind characteristic variables
    SArray<real,numState> chu;  // upwind characteristic variables
    chu(0) = -r/(2*cs)*w2 + p2/(2*cs2); // w-cs wave: assuming no shocks
    chu(1) =  r/(2*cs)*w1 + p1/(2*cs2); // w+cs wave: assuming no shocks
    chu(2) = ru - pu/cs2;
    chu(3) = uu;
    chu(4) = vu;
    // Next, compute the upwind state based on upwind characteristic variables
    ru =     chu(0)   +     chu(1)   + chu(2);
    uu =                                       chu(3);
    vu =                                               chu(4);
    wu = -cs*chu(0)/r +  cs*chu(1)/r;
    pu = cs2*chu(0)   + cs2*chu(1);
    // Finally, compute the upwind flux based on the upwind state
    real keu = 0.5_fp*ru*(uu*uu+vu*vu+wu*wu); // upwind kinetic energy
    real reu = pu*CV/RD + keu;                // upwind rho*e
    flux_r (k,j,i) = ru*wu;                   // upwind mass flux
    flux_re(k,j,i) = wu*reu + wu*pu;          // upwind energy flux
  });

  ///////////////////////////////////////////////////////////////////////////////////////////////////////
  // For each cell:
  // (1) Append the u, v, and w tendencies with the flux difference splitting waves entering the cell
  //     domain
  // (2) Compute the rho and rho*e tendencies using the upwind flux vectors
  ///////////////////////////////////////////////////////////////////////////////////////////////////////
  //   for (int k=0; k<dom.nz; k++) {
  //     for (int j=0; j<dom.ny; j++) {
  //       for (int i=0; i<dom.nx; i++) {
  yakl::parallel_for( dom.nz*dom.ny*dom.nx , YAKL_LAMBDA (int const iGlob) {
    int l, k, j, i;
    yakl::unpackIndices(iGlob,dom.nz,dom.ny,dom.nx,k,j,i);
    // Flux vector form update for mass and energy
    tend(idR,k,j,i)  = - ( flux_r (k+1,j,i) - flux_r (k,j,i) ) / dom.dz;  // mass tendency
    tend(idT,k,j,i)  = - ( flux_re(k+1,j,i) - flux_re(k,j,i) ) / dom.dz;  // energy tendency
    // Flux difference splitting form update for velocities
    tend(idU,k,j,i) += - ( stateLimits(idU,1,k,j,i) + stateLimits(idU,0,k+1,j,i) ) / dom.dz;  // u tendency
    tend(idV,k,j,i) += - ( stateLimits(idV,1,k,j,i) + stateLimits(idV,0,k+1,j,i) ) / dom.dz;  // v tendency
    tend(idW,k,j,i) += - ( stateLimits(idW,1,k,j,i) + stateLimits(idW,0,k+1,j,i) ) / dom.dz;  // w tendency
  });
}



////////////////////////////////////////////////////////////////////////////////////////////
// Compute the source term tendencies
// 
// INPUTS
//   state: The fluid state vector, rho,u,v,w,rho*e; dims(numState,nz+2*hs,ny+2*hs,nx+2*hs)
//   dom: The Domain class object
// 
// OUTPUTS
//   tend: State tendency; dims(numStat,nz,ny,nx)
////////////////////////////////////////////////////////////////////////////////////////////
void Tendencies::compEulerTend_S(realArr const &state, Domain const &dom, realArr &tend) {
  // for (int k=0; k<dom.nz; k++) {
  //   for (int j=0; j<dom.ny; j++) {
  //     for (int i=0; i<dom.nx; i++) {
  yakl::parallel_for( dom.nz*dom.ny*dom.nx , YAKL_LAMBDA (int const iGlob) {
    int k, j, i;
    yakl::unpackIndices(iGlob,dom.nz,dom.ny,dom.nx,k,j,i);
    real rp = state(idR,hs+k,hs+j,hs+i); // density perturbation
    real rh = dom.hyDensCells(hs+k);     // hydrostatic density
    tend(idR,k,j,i) = 0;
    tend(idU,k,j,i) = 0;
    tend(idV,k,j,i) = 0;
    tend(idW,k,j,i) = -rp/(rp+rh)*GRAV;
    tend(idT,k,j,i) = 0;
  });
}



void Tendencies::compStrakaTend(realArr &state, Domain const &dom, Exchange &exch, Parallel const &par, 
                                realArr &tend) {
}


////////////////////////////////////////////////////////////////////////////////////////////
// Fill in vertical halo cells with boundary forcing data. vertical velocity is zero, and
// all other variables replicate the last domain value to mimick zero gradient free-slip
// 
// INPUTS
//   state: The fluid state vector, rho,u,v,w,rho*e; dims(numState,nz+2*hs,ny+2*hs,nx+2*hs)
//   dom: The Domain class object
// 
// OUTPUTS
//   state: fluid state with halos filled with boundary forcing data
////////////////////////////////////////////////////////////////////////////////////////////
void Tendencies::stateBoundariesZ(realArr &state, Domain const &dom) {
  // for (int j=0; j<dom.ny; j++) {
  //   for (int i=0; i<dom.nx; i++) {
  //     for (int ii=0; ii<hs; ii++) {
  yakl::parallel_for( dom.ny*dom.nx*hs , YAKL_LAMBDA (int const iGlob) {
    int j, i, ii;
    yakl::unpackIndices(iGlob,dom.ny,dom.nx,hs,j,i,ii);
    state(idR,ii,hs+j,hs+i) = state(idR,hs,hs+j,hs+i);
    state(idU,ii,hs+j,hs+i) = state(idU,hs,hs+j,hs+i);
    state(idV,ii,hs+j,hs+i) = state(idV,hs,hs+j,hs+i);
    state(idW,ii,hs+j,hs+i) = 0;
    state(idT,ii,hs+j,hs+i) = state(idT,hs,hs+j,hs+i);

    state(idR,dom.nz+hs+ii,hs+j,hs+i) = state(idR,dom.nz+hs-1,hs+j,hs+i);
    state(idU,dom.nz+hs+ii,hs+j,hs+i) = state(idU,dom.nz+hs-1,hs+j,hs+i);
    state(idV,dom.nz+hs+ii,hs+j,hs+i) = state(idV,dom.nz+hs-1,hs+j,hs+i);
    state(idW,dom.nz+hs+ii,hs+j,hs+i) = 0;
    state(idT,dom.nz+hs+ii,hs+j,hs+i) = state(idT,dom.nz+hs-1,hs+j,hs+i);
  });
}


////////////////////////////////////////////////////////////////////////////////////////////
// Fill in stateLimits vertical domain edge data with boundary forcing
// Vertical velocity is zero, and all other variables replicate the last domain value to 
// mimick zero gradient free-slip
// 
// INPUTS
//   state: stateLimits for rho,u,v,w,p ;  dims(numState,2,nz+1,ny+1,nx+1)
//   dom: The Domain class object
// 
// OUTPUTS
//   state: stateLimits with vertical domain edge data filled in
////////////////////////////////////////////////////////////////////////////////////////////
void Tendencies::edgeBoundariesZ(realArr &stateLimits, Domain const &dom) {
  // for (int j=0; j<dom.ny; j++) {
  //   for (int i=0; i<dom.nx; i++) {
  yakl::parallel_for( dom.ny*dom.nx , YAKL_LAMBDA (int const iGlob) {
    int j, i;
    yakl::unpackIndices(iGlob,dom.ny,dom.nx,j,i);
    stateLimits(idR,0,0     ,j,i) = stateLimits(idR,1,0     ,j,i);
    stateLimits(idU,0,0     ,j,i) = stateLimits(idU,1,0     ,j,i);
    stateLimits(idV,0,0     ,j,i) = stateLimits(idV,1,0     ,j,i);
    stateLimits(idW,0,0     ,j,i) = 0;
    stateLimits(idW,1,0     ,j,i) = 0;
    stateLimits(idT,0,0     ,j,i) = stateLimits(idT,1,0     ,j,i);

    stateLimits(idR,1,dom.nz,j,i) = stateLimits(idR,0,dom.nz,j,i);
    stateLimits(idU,1,dom.nz,j,i) = stateLimits(idU,0,dom.nz,j,i);
    stateLimits(idV,1,dom.nz,j,i) = stateLimits(idV,0,dom.nz,j,i);
    stateLimits(idW,0,dom.nz,j,i) = 0;
    stateLimits(idW,1,dom.nz,j,i) = 0;
    stateLimits(idT,1,dom.nz,j,i) = stateLimits(idT,0,dom.nz,j,i);
  });
}



