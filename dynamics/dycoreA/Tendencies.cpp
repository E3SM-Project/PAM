
#include "Tendencies.h"



////////////////////////////////////////////////////////////////////////////////////////////////////
// Allocate arrays, pre-compute transformation arrays, WENO stuff, and GLL weights
// 
// INPUTS
//   dom: The Domain object
////////////////////////////////////////////////////////////////////////////////////////////////////
void Tendencies::initialize(Domain const &dom) {
  TransformMatrices<real> trans;

  stateLimits = realArr("stateLimits",numState+1,2,dom.nz+1,dom.ny+1,dom.nx+1);

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
  s2d2gX = ( to_gll * c2d_ho * s2c_ho ) / dom.dx;
  s2d2gY = ( to_gll * c2d_ho * s2c_ho ) / dom.dy;
  s2d2gZ = ( to_gll * c2d_ho * s2c_ho ) / dom.dz;

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

  trans.sten_to_gll_lower(s2g);

  trans.weno_sten_to_coefs(wenoRecon);

  trans.get_gll_weights(gllWts);

  wenoSetIdealSigma(wenoIdl,wenoSigma);
}



///////////////////////////////////////////////////////////////////////////////////////////////////////
// Compute time tendencies for rho, u, v, w, and theta using ADER Differential Transform time
// stepping in the x-direction
// 
// INPUTS
//   state: The state vector: rho,u,v,w,theta. dims  (numState,nz+2*hs,ny+2*hs,nx+2*hs)
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
  // (1) Reconstruct state values & spatial derivatives for rho, u, v, w, and theta at tord GLL points
  //     across the cell
  // (2) Compute time derivatives for rho, u, v, w, and theta at each GLL point using Differential
  //     Transforms in time. Save tendencies as a by product of this operation
  // (3) Compute time averages of the state and the state time tendencies
  // (4) Compute local flux difference splitting tendencies using spatial quadrature over tord GLL points
  // (5) Save the time-averaged state value estimates at the left and right of the cell for use in the
  //     upwind Riemann solver in the next step
  ///////////////////////////////////////////////////////////////////////////////////////////////////////
  // for (int k=0; k<dom.nz; k++) {
  //   for (int j=0; j<dom.ny; j++) {
  //     for (int i=0; i<dom.nx; i++) {
  yakl::parallel_for( dom.nz*dom.ny*dom.nx , YAKL_LAMBDA ( int const iGlob ) {
    int k, j, i;
    yakl::unpackIndices(iGlob,dom.nz,dom.ny,dom.nx,k,j,i);

    SArray<real,numState,tord,tord> stateDTs;  // GLL state DTs    (var,time,space)
    SArray<real,numState,tord,tord> derivDTs;  // GLL deriv DTs    (var,time,space)
    SArray<real,numState,tord,tord> tendDTs;   // GLL tendency DTs (var,time,space)

    // Compute tord GLL points of the fluid state and spatial derivative
    for (int l=0; l<numState; l++) {
      SArray<real,ord> stencil;
      SArray<real,tord> gllPts;

      // Store the stencil values
      for (int ii=0; ii<ord; ii++) {
        stencil(ii) = state(l,hs+k,hs+j,i+ii);
      }

      if (l == idT) {
        // Reconstruct and store GLL points of the state values
        reconStencil(stencil, gllPts, dom.doWeno, wenoRecon, to_gll       , wenoIdl, wenoSigma);
        for (int ii=0; ii<tord; ii++) { stateDTs(l,0,ii) = gllPts(ii); }

        // Reconstruct and store GLL points of the state derivatives
        reconStencil(stencil, gllPts, dom.doWeno, wenoRecon, to_derivX_gll, wenoIdl, wenoSigma);
        for (int ii=0; ii<tord; ii++) { derivDTs(l,0,ii) = gllPts(ii); }
      } else {
        // Reconstruct and store GLL points of the state values
        reconStencil(stencil, gllPts, 0         , wenoRecon, s2g          , wenoIdl, wenoSigma);
        for (int ii=0; ii<tord; ii++) { stateDTs(l,0,ii) = gllPts(ii); }

        // Reconstruct and store GLL points of the state derivatives
        reconStencil(stencil, gllPts, 0         , wenoRecon, s2d2gX       , wenoIdl, wenoSigma);
        for (int ii=0; ii<tord; ii++) { derivDTs(l,0,ii) = gllPts(ii); }
      }
    }
    // Add hydrostasis to density and potential temperature to make them the full quantities
    for (int ii=0; ii<tord; ii++) {
      stateDTs(idR,0,ii) += dom.hyDensCells (hs+k);
      stateDTs(idT,0,ii) += dom.hyThetaCells(hs+k);
    }

    // Compute tord-1 time derivatives of the state, state spatial derivatives, 
    // and state tendencies using temporal Differential Transforms
    diffTransformEulerX( stateDTs, derivDTs, tendDTs, aderDerivX );

    // Compute the time-average and store into the zeroth time index
    timeAvg( stateDTs , dom );
    timeAvg( tendDTs  , dom );

    // Compute the local tendency contribution for high-order flux difference
    // splitting for wind via quadrature
    // tend_local = int( RHS , x , x_(i-1/2) , x_(i+1/2) )
    for (int l=0; l<numState; l++) {
      tend(l,k,j,i) = 0;
      for (int ii=0; ii<tord; ii++) {
        tend(l,k,j,i) += gllWts(ii) * tendDTs(l,0,ii);
      }
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
    real t = 0.5_fp * ( stateLimits(idT,1,k,j,i) + stateLimits(idT,0,k,j,i)); // theta
    real p = C0*pow(r*t,GAMMA);
    real cs2 = GAMMA*p/r;    // speed of sound squared
    real cs = sqrt(cs2);     // speed of sound
    // Compute the state jump across the interface
    real dr = stateLimits(idR,1,k,j,i) - stateLimits(idR,0,k,j,i);
    real du = stateLimits(idU,1,k,j,i) - stateLimits(idU,0,k,j,i);
    real dv = stateLimits(idV,1,k,j,i) - stateLimits(idV,0,k,j,i);
    real dw = stateLimits(idW,1,k,j,i) - stateLimits(idW,0,k,j,i);
    real dt = stateLimits(idT,1,k,j,i) - stateLimits(idT,0,k,j,i);
    // Compute the product of the flux Jacobian and the state jump across the interface (A*dq)
    SArray<real,numState> df;
    df(0) = u*dr + r*du;
    df(1) = u*du + cs2/r*dr + cs2/t*dt;
    df(2) = u*dv;
    df(3) = u*dw;
    df(4) = u*dt;
    // Zero out the stateLimits space for this spatial index
    for (int l=0; l<numState; l++) {
      stateLimits(l,0,k,j,i) = 0;
      stateLimits(l,1,k,j,i) = 0;
    }
    real ch;
    // Wave 1 (u-cs): presumed negatively propagating (no shocks)
    ch = 0.5_fp*df(0) - r/(2*cs)*df(1) + r/(2*t)*df(4);
    stateLimits(idR,0,k,j,i) += ch;
    stateLimits(idU,0,k,j,i) += -cs/r * ch;
    // Wave 2 (u+cs): presumed positively propagating (no shocks)
    ch = 0.5_fp*df(0) + r/(2*cs)*df(1) + r/(2*t)*df(4);
    stateLimits(idR,1,k,j,i) += ch;
    stateLimits(idU,1,k,j,i) +=  cs/r * ch;
    if (u > 0) {
      stateLimits(idR,1,k,j,i) += -r/t*df(4); // Wave 3 (u)
      stateLimits(idT,1,k,j,i) += df(4);      // Wave 3 (u)
      stateLimits(idV,1,k,j,i) += df(2);      // Wave 4 (u)
      stateLimits(idW,1,k,j,i) += df(3);      // Wave 5 (u)
    } else {
      stateLimits(idR,0,k,j,i) += -r/t*df(4); // Wave 3 (u)
      stateLimits(idT,0,k,j,i) += df(4);      // Wave 3 (u)
      stateLimits(idV,0,k,j,i) += df(2);      // Wave 4 (u)
      stateLimits(idW,0,k,j,i) += df(3);      // Wave 5 (u)
    }
  });

  ////////////////////////////////////////////////////////////////////////////////////////////
  // Add interface fluctuation contributions to each cell, combining the postiive-propatating
  // waves from the left interfact and the necative-propagating waves from the right interface
  ////////////////////////////////////////////////////////////////////////////////////////////
  // for (int l=0; l<numState; l++) {
  //   for (int k=0; k<dom.nz; k++) {
  //     for (int j=0; j<dom.ny; j++) {
  //       for (int i=0; i<dom.nx; i++) {
  yakl::parallel_for( numState*dom.nz*dom.ny*dom.nx , YAKL_LAMBDA (int const iGlob) {
    int l, k, j, i;
    yakl::unpackIndices(iGlob,numState,dom.nz,dom.ny,dom.nx,l,k,j,i);

    tend(l,k,j,i) += - ( stateLimits(l,1,k,j,i) + stateLimits(l,0,k,j,i+1) ) / dom.dx;
  });
}



///////////////////////////////////////////////////////////////////////////////////////////////////////
// Compute time tendencies for rho, u, v, w, and theta using ADER Differential Transform time
// stepping in the x-direction
// 
// INPUTS
//   state: The state vector: rho,u,v,w,theta. dims  (numState,nz+2*hs,ny+2*hs,nx+2*hs)
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
  // (1) Reconstruct state values & spatial derivatives for rho, u, v, w, and theta at tord GLL points
  //     across the cell
  // (2) Compute time derivatives for rho, u, v, w, and theta at each GLL point using Differential
  //     Transforms in time. Save tendencies as a by product of this operation
  // (3) Compute time averages of the state and the state time tendencies
  // (4) Compute local flux difference splitting tendencies using spatial quadrature over tord GLL points
  // (5) Save the time-averaged state value estimates at the left and right of the cell for use in the
  //     upwind Riemann solver in the next step
  ///////////////////////////////////////////////////////////////////////////////////////////////////////
  // for (int k=0; k<dom.nz; k++) {
  //   for (int j=0; j<dom.ny; j++) {
  //     for (int i=0; i<dom.nx; i++) {
  yakl::parallel_for( dom.nz*dom.ny*dom.nx , YAKL_LAMBDA ( int const iGlob ) {
    int k, j, i;
    yakl::unpackIndices(iGlob,dom.nz,dom.ny,dom.nx,k,j,i);

    SArray<real,numState,tord,tord> stateDTs;  // GLL state DTs    (var,time,space)
    SArray<real,numState,tord,tord> derivDTs;  // GLL deriv DTs    (var,time,space)
    SArray<real,numState,tord,tord> tendDTs;   // GLL tendency DTs (var,time,space)

    // Compute tord GLL points of the fluid state and spatial derivative
    for (int l=0; l<numState; l++) {
      SArray<real,ord> stencil;
      SArray<real,tord> gllPts;

      // Store the stencil values
      for (int ii=0; ii<ord; ii++) {
        stencil(ii) = state(l,hs+k,j+ii,hs+i);
      }

      if (l == idT) {
        // Reconstruct and store GLL points of the state values
        reconStencil(stencil, gllPts, dom.doWeno, wenoRecon, to_gll       , wenoIdl, wenoSigma);
        for (int ii=0; ii<tord; ii++) { stateDTs(l,0,ii) = gllPts(ii); }

        // Reconstruct and store GLL points of the state derivatives
        reconStencil(stencil, gllPts, dom.doWeno, wenoRecon, to_derivY_gll, wenoIdl, wenoSigma);
        for (int ii=0; ii<tord; ii++) { derivDTs(l,0,ii) = gllPts(ii); }
      } else {
        // Reconstruct and store GLL points of the state values
        reconStencil(stencil, gllPts, 0         , wenoRecon, s2g          , wenoIdl, wenoSigma);
        for (int ii=0; ii<tord; ii++) { stateDTs(l,0,ii) = gllPts(ii); }

        // Reconstruct and store GLL points of the state derivatives
        reconStencil(stencil, gllPts, 0         , wenoRecon, s2d2gY       , wenoIdl, wenoSigma);
        for (int ii=0; ii<tord; ii++) { derivDTs(l,0,ii) = gllPts(ii); }
      }
    }
    // Add hydrostasis to density and potential temperature to make them the full quantities
    for (int ii=0; ii<tord; ii++) {
      stateDTs(idR,0,ii) += dom.hyDensCells (hs+k);
      stateDTs(idT,0,ii) += dom.hyThetaCells(hs+k);
    }

    // Compute tord-1 time derivatives of the state, state spatial derivatives, 
    // and state tendencies using temporal Differential Transforms
    diffTransformEulerY( stateDTs, derivDTs, tendDTs, aderDerivY );

    // Compute the time-average and store into the zeroth time index
    timeAvg( stateDTs , dom );
    timeAvg( tendDTs  , dom );

    // Compute the local tendency contribution for high-order flux difference
    // splitting for wind via quadrature
    // tend_local = int( RHS , y , y_(j-1/2) , y_(j+1/2) )
    for (int l=0; l<numState; l++) {
      tend(l,k,j,i) = 0;
      for (int ii=0; ii<tord; ii++) {
        tend(l,k,j,i) += gllWts(ii) * tendDTs(l,0,ii);
      }
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
    real t = 0.5_fp * ( stateLimits(idT,1,k,j,i) + stateLimits(idT,0,k,j,i)); // theta
    real p = C0*pow(r*t,GAMMA);
    real cs2 = GAMMA*p/r;    // speed of sound squared
    real cs = sqrt(cs2);     // speed of sound
    // Compute the state jump across the interface
    real dr = stateLimits(idR,1,k,j,i) - stateLimits(idR,0,k,j,i);
    real du = stateLimits(idU,1,k,j,i) - stateLimits(idU,0,k,j,i);
    real dv = stateLimits(idV,1,k,j,i) - stateLimits(idV,0,k,j,i);
    real dw = stateLimits(idW,1,k,j,i) - stateLimits(idW,0,k,j,i);
    real dt = stateLimits(idT,1,k,j,i) - stateLimits(idT,0,k,j,i);
    // Compute the product of the flux Jacobian and the state jump across the interface (A*dq)
    SArray<real,numState> df;
    df(0) = v*dr + r*dv;
    df(1) = v*du;
    df(2) = v*dv + cs2/r*dr + cs2/t*dt;
    df(3) = v*dw;
    df(4) = v*dt;
    // Zero out the stateLimits space for this spatial index
    for (int l=0; l<numState; l++) {
      stateLimits(l,0,k,j,i) = 0;
      stateLimits(l,1,k,j,i) = 0;
    }
    real ch;
    // Wave 1 (v-cs): presumed negatively propagating (no shocks)
    ch = 0.5_fp*df(0) - r/(2*cs)*df(2) + r/(2*t)*df(4);
    stateLimits(idR,0,k,j,i) += ch;
    stateLimits(idV,0,k,j,i) += -cs/r * ch;
    // Wave 2 (v+cs): presumed positively propagating (no shocks)
    ch = 0.5_fp*df(0) + r/(2*cs)*df(2) + r/(2*t)*df(4);
    stateLimits(idR,1,k,j,i) += ch;
    stateLimits(idV,1,k,j,i) +=  cs/r * ch;
    if (v > 0) {
      stateLimits(idR,1,k,j,i) += -r/t*df(4); // Wave 3 (v)
      stateLimits(idT,1,k,j,i) += df(4);      // Wave 3 (v)
      stateLimits(idU,1,k,j,i) += df(1);      // Wave 4 (v)
      stateLimits(idW,1,k,j,i) += df(3);      // Wave 5 (v)
    } else {
      stateLimits(idR,0,k,j,i) += -r/t*df(4); // Wave 3 (v)
      stateLimits(idT,0,k,j,i) += df(4);      // Wave 3 (v)
      stateLimits(idU,0,k,j,i) += df(1);      // Wave 4 (v)
      stateLimits(idW,0,k,j,i) += df(3);      // Wave 5 (v)
    }
  });

  ////////////////////////////////////////////////////////////////////////////////////////////
  // Add interface fluctuation contributions to each cell, combining the postiive-propatating
  // waves from the left interfact and the necative-propagating waves from the right interface
  ////////////////////////////////////////////////////////////////////////////////////////////
  // for (int l=0; l<numState; l++) {
  //   for (int k=0; k<dom.nz; k++) {
  //     for (int j=0; j<dom.ny; j++) {
  //       for (int i=0; i<dom.nx; i++) {
  yakl::parallel_for( numState*dom.nz*dom.ny*dom.nx , YAKL_LAMBDA (int const iGlob) {
    int l, k, j, i;
    yakl::unpackIndices(iGlob,numState,dom.nz,dom.ny,dom.nx,l,k,j,i);

    tend(l,k,j,i) += - ( stateLimits(l,1,k,j,i) + stateLimits(l,0,k,j+1,i) ) / dom.dy;
  });
}



///////////////////////////////////////////////////////////////////////////////////////////////////////
// Compute time tendencies for rho, u, v, w, and theta using ADER Differential Transform time
// stepping in the x-direction
// 
// INPUTS
//   state: The state vector: rho,u,v,w,theta. dims  (numState,nz+2*hs,ny+2*hs,nx+2*hs)
//   dom: The Domain class object   
// 
// OUTPUTS
//   state: The state vector unchanged except for halos 
//   tend: The time tendencies of the state vector   (numstate,nz,ny,nx)
///////////////////////////////////////////////////////////////////////////////////////////////////////
void Tendencies::compEulerTend_Z(realArr &state, Domain const &dom, realArr &tend) {
  // Create local references of class arrays so that C++ Lambdas will properly capture them
  auto &stateLimits   = this->stateLimits  ;
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
  // (1) Reconstruct state values & spatial derivatives for rho, u, v, w, theta, and pressure at tord GLL
  //     points across the cell
  // (2) Compute time derivatives for rho, u, v, w, theta, and pressure at each GLL point using
  //     Differential Transforms in time. Save tendencies for rho, u, v, w, and theta as a by product of
  //     this operation
  // (3) Compute time averages of the state and the state time tendencies
  // (4) Compute local flux difference splitting tendencies using spatial quadrature over tord GLL points
  // (5) Save the time-averaged state value estimates at the left and right of the cell for use in the
  //     upwind Riemann solver in the next step
  ///////////////////////////////////////////////////////////////////////////////////////////////////////
  // for (int k=0; k<dom.nz; k++) {
  //   for (int j=0; j<dom.ny; j++) {
  //     for (int i=0; i<dom.nx; i++) {
  yakl::parallel_for( dom.nz*dom.ny*dom.nx , YAKL_LAMBDA ( int const iGlob ) {
    int k, j, i;
    yakl::unpackIndices(iGlob,dom.nz,dom.ny,dom.nx,k,j,i);

    SArray<real,numState+1,tord,tord> stateDTs;  // GLL state DTs    (var,time,space)
    SArray<real,numState+1,tord,tord> derivDTs;  // GLL deriv DTs    (var,time,space)
    SArray<real,numState  ,tord,tord> tendDTs;   // GLL tendency DTs (var,time,space)

    // Compute tord GLL points of the fluid state and spatial derivative
    for (int l=0; l<numState+1; l++) {
      SArray<real,ord> stencil;
      SArray<real,tord> gllPts;
      // Store the stencil values
      for (int ii=0; ii<ord; ii++) {
        if (l < numState) {
          stencil(ii) = state(l,k+ii,hs+j,hs+i);
        } else {
          real r = state(idR,k+ii,hs+j,hs+i) + dom.hyDensCells (k+ii);
          real t = state(idT,k+ii,hs+j,hs+i) + dom.hyThetaCells(k+ii);
          stencil(ii) = C0*pow(r*t,GAMMA) - dom.hyPressureCells(k+ii);
        }
      }

      if (l == idT) {
        // Reconstruct and store GLL points of the state values
        reconStencil(stencil, gllPts, dom.doWeno, wenoRecon, to_gll       , wenoIdl, wenoSigma);
        for (int ii=0; ii<tord; ii++) { stateDTs(l,0,ii) = gllPts(ii); }

        // Reconstruct and store GLL points of the state derivatives
        reconStencil(stencil, gllPts, dom.doWeno, wenoRecon, to_derivZ_gll, wenoIdl, wenoSigma);
        for (int ii=0; ii<tord; ii++) { derivDTs(l,0,ii) = gllPts(ii); }
      } else {
        // Reconstruct and store GLL points of the state values
        reconStencil(stencil, gllPts, 0         , wenoRecon, s2g          , wenoIdl, wenoSigma);
        for (int ii=0; ii<tord; ii++) { stateDTs(l,0,ii) = gllPts(ii); }

        // Reconstruct and store GLL points of the state derivatives
        reconStencil(stencil, gllPts, 0         , wenoRecon, s2d2gZ       , wenoIdl, wenoSigma);
        for (int ii=0; ii<tord; ii++) { derivDTs(l,0,ii) = gllPts(ii); }
      }
    }
    // Add hydrostasis to density, theta, and pressure to make them the full quantities
    for (int ii=0; ii<tord; ii++) {
      stateDTs(idR,0,ii) += dom.hyDensGLL    (k,ii);
      stateDTs(idT,0,ii) += dom.hyThetaGLL   (k,ii);
      stateDTs(idP,0,ii) += dom.hyPressureGLL(k,ii);
    }
    // Enforce boundary conditions on the reconstructed vertical wind
    if (k == dom.nz-1) { stateDTs(idW,0,tord-1) = 0; }
    if (k == 0       ) { stateDTs(idW,0,0     ) = 0; }

    // Compute tord-1 time derivatives of the state, state spatial derivatives, 
    // and state tendencies using temporal Differential Transforms
    diffTransformEulerZ( stateDTs, derivDTs, tendDTs, aderDerivZ );

    // Compute the time-average and store into the zeroth time index
    timeAvg( stateDTs , dom );
    timeAvg( tendDTs  , dom );

    // Enforce boundary conditions on the time-averaged vertical wind
    if (k == dom.nz-1) { stateDTs(idW,0,tord-1) = 0; }
    if (k == 0       ) { stateDTs(idW,0,0     ) = 0; }

    // Compute the local tendency contribution for high-order flux difference
    // splitting for wind via quadrature
    // tend_local = int( RHS , z , z_(k-1/2) , z_(k+1/2) )
    for (int l=0; l<numState; l++) {
      tend(l,k,j,i) = 0;
      for (int ii=0; ii<tord; ii++) {
        tend(l,k,j,i) += gllWts(ii) * tendDTs(l,0,ii);
      }
    }

    // Store the state vector in fwaves to compute fwaves from cell-interface state jumps
    for (int l=0; l<numState+1; l++) {
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
    real r = 0.5_fp * ( stateLimits(idR,1,k,j,i) + stateLimits(idR,0,k,j,i) ); // rho
    real w = 0.5_fp * ( stateLimits(idW,1,k,j,i) + stateLimits(idW,0,k,j,i) ); // w
    real t = 0.5_fp * ( stateLimits(idT,1,k,j,i) + stateLimits(idT,0,k,j,i) ); // theta
    real p = 0.5_fp * ( stateLimits(idP,1,k,j,i) + stateLimits(idP,0,k,j,i) ); // pressure
    real cs2 = GAMMA*p/r;    // speed of sound squared
    real cs = sqrt(cs2);     // speed of sound
    // Compute the state jump across the interface
    real dr = stateLimits(idR,1,k,j,i) - stateLimits(idR,0,k,j,i);
    real du = stateLimits(idU,1,k,j,i) - stateLimits(idU,0,k,j,i);
    real dv = stateLimits(idV,1,k,j,i) - stateLimits(idV,0,k,j,i);
    real dw = stateLimits(idW,1,k,j,i) - stateLimits(idW,0,k,j,i);
    real dt = stateLimits(idT,1,k,j,i) - stateLimits(idT,0,k,j,i);
    real dp = stateLimits(idP,1,k,j,i) - stateLimits(idP,0,k,j,i);
    // Compute the product of the flux Jacobian and the state jump across the interface (A*dq)
    SArray<real,numState> df;
    df(0) = w*dr + r*dw;
    df(1) = w*du;
    df(2) = w*dv;
    df(3) = w*dw + dp/r;
    df(4) = w*dt;
    // Zero out the stateLimits space for this spatial index
    for (int l=0; l<numState; l++) {
      stateLimits(l,0,k,j,i) = 0;
      stateLimits(l,1,k,j,i) = 0;
    }
    real ch;
    // Wave 1 (w-cs): presumed negatively propagating (no shocks)
    ch = 0.5_fp*df(0) - r/(2*cs)*df(3) + r/(2*t)*df(4);
    stateLimits(idR,0,k,j,i) += ch;
    stateLimits(idW,0,k,j,i) += -cs/r * ch;
    // Wave 2 (w+cs): presumed positively propagating (no shocks)
    ch = 0.5_fp*df(0) + r/(2*cs)*df(3) + r/(2*t)*df(4);
    stateLimits(idR,1,k,j,i) += ch;
    stateLimits(idW,1,k,j,i) +=  cs/r * ch;
    if (w > 0) {
      stateLimits(idR,1,k,j,i) += -r/t*df(4); // Wave 3 (w)
      stateLimits(idT,1,k,j,i) += df(4);      // Wave 3 (w)
      stateLimits(idU,1,k,j,i) += df(1);      // Wave 4 (w)
      stateLimits(idV,1,k,j,i) += df(2);      // Wave 5 (w)
    } else {
      stateLimits(idR,0,k,j,i) += -r/t*df(4); // Wave 3 (w)
      stateLimits(idT,0,k,j,i) += df(4);      // Wave 3 (w)
      stateLimits(idU,0,k,j,i) += df(1);      // Wave 4 (w)
      stateLimits(idV,0,k,j,i) += df(2);      // Wave 5 (w)
    }
  });

  ////////////////////////////////////////////////////////////////////////////////////////////
  // Add interface fluctuation contributions to each cell, combining the postiive-propatating
  // waves from the left interfact and the necative-propagating waves from the right interface
  ////////////////////////////////////////////////////////////////////////////////////////////
  // for (int l=0; l<numState; l++) {
  //   for (int k=0; k<dom.nz; k++) {
  //     for (int j=0; j<dom.ny; j++) {
  //       for (int i=0; i<dom.nx; i++) {
  yakl::parallel_for( numState*dom.nz*dom.ny*dom.nx , YAKL_LAMBDA (int const iGlob) {
    int l, k, j, i;
    yakl::unpackIndices(iGlob,numState,dom.nz,dom.ny,dom.nx,l,k,j,i);

    tend(l,k,j,i) += - ( stateLimits(l,1,k,j,i) + stateLimits(l,0,k+1,j,i) ) / dom.dz;
  });
}



////////////////////////////////////////////////////////////////////////////////////////////
// Compute the source term tendencies
// 
// INPUTS
//   state: The fluid state vector, rho,u,v,w,theta; dims(numState,nz+2*hs,ny+2*hs,nx+2*hs)
//   dom: The Domain class object
// 
// OUTPUTS
//   tend: State tendency; dims(numState,nz,ny,nx)
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



// void Tendencies::compStrakaTend(realArr &state, Domain const &dom, Exchange &exch, Parallel const &par, 
//                                 realArr &tend) {
// }



////////////////////////////////////////////////////////////////////////////////////////////
// Fill in vertical halo cells with boundary forcing data. vertical velocity is zero, and
// all other variables replicate the last domain value to mimick zero gradient free-slip
// 
// INPUTS
//   state: The fluid state vector, rho,u,v,w,theta; dims(numState,nz+2*hs,ny+2*hs,nx+2*hs)
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
//   state: stateLimits for rho,u,v,w,theta,pressure ;  dims(numState+1,2,nz+1,ny+1,nx+1)
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
    stateLimits(idP,0,0     ,j,i) = stateLimits(idP,1,0     ,j,i);

    stateLimits(idR,1,dom.nz,j,i) = stateLimits(idR,0,dom.nz,j,i);
    stateLimits(idU,1,dom.nz,j,i) = stateLimits(idU,0,dom.nz,j,i);
    stateLimits(idV,1,dom.nz,j,i) = stateLimits(idV,0,dom.nz,j,i);
    stateLimits(idW,0,dom.nz,j,i) = 0;
    stateLimits(idW,1,dom.nz,j,i) = 0;
    stateLimits(idT,1,dom.nz,j,i) = stateLimits(idT,0,dom.nz,j,i);
    stateLimits(idP,1,dom.nz,j,i) = stateLimits(idP,0,dom.nz,j,i);
  });
}



