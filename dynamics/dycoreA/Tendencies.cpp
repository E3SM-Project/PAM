
#include "Tendencies.h"



void Tendencies::initialize(Domain const &dom) {
  TransformMatrices<real> trans;

  stateLimits = realArr("stateLimits",numState,2,dom.nz+1,dom.ny+1,dom.nx+1);
  flux_r      = realArr("flux_r"                ,dom.nz+1,dom.ny+1,dom.nx+1);
  flux_re     = realArr("flux_re"               ,dom.nz+1,dom.ny+1,dom.nx+1);
  dph_glob    = realArr("dph_glob"              ,dom.nz,tord);

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

  // Compute vertical derivative of hydrostatic pressure
  for (int k=0; k<dom.nz; k++) {
    SArray<real,ord> stencil;
    SArray<real,tord> gllPts;
    for (int ii=0; ii<ord; ii++) {
      stencil(ii) = dom.hyPressureCells(k+ii);
    }
    // Reconstruct and store GLL points of the pressure perturbation spatial derivatives
    reconStencil(stencil, gllPts, dom.doWeno, wenoRecon, to_derivX_gll, wenoIdl, wenoSigma);
    for (int ii=0; ii<tord; ii++) { dph_glob(k,ii) = gllPts(ii); }
  }
}



void Tendencies::compEulerTend_X(realArr &state, Domain const &dom, Exchange &exch, Parallel const &par, realArr &tend) {
  auto &fwaves      = this->fwaves     ;
  auto &src         = this->src        ;
  auto &gllWts      = this->gllWts     ;
  auto &to_gll      = this->to_gll     ;
  auto &wenoRecon   = this->wenoRecon  ;
  auto &wenoIdl     = this->wenoIdl    ;
  auto &wenoSigma   = this->wenoSigma  ;
  auto &aderDerivX  = this->aderDerivX ;

  // Exchange halos in the x-direction
  exch.haloInit      ();
  exch.haloPackN_x   (dom, state, numState);
  exch.haloExchange_x(dom, par);
  exch.haloUnpackN_x (dom, state, numState);

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
      if (l != idT) {  // rho, u, v, and w
        for (int ii=0; ii<ord; ii++) {
          stencil(ii) = state(l,hs+k,hs+j,i+ii);
        }
      } else {         // pressure
        for (int ii=0; ii<ord; ii++) {
          real r  = state(idR,hs+k,hs+j,i+ii) + dom.hyDensCells  (hs+k);
          real u  = state(idU,hs+k,hs+j,i+ii);
          real v  = state(idV,hs+k,hs+j,i+ii);
          real w  = state(idW,hs+k,hs+j,i+ii);
          real re = state(idT,hs+k,hs+j,i+ii) + dom.hyEnergyCells(hs+k);
          real ke = 0.5_fp*r*(u*u+v*v+w*w);
          real p = RD/CV*(re-ke);
          stencil(ii) = p;
        }
      }

      // Reconstruct and store GLL points of the state values
      reconStencil(stencil, gllPts, dom.doWeno, wenoRecon, to_gll       , wenoIdl, wenoSigma);
      for (int ii=0; ii<tord; ii++) { stateDTs(l,0,ii) = gllPts(ii); }

      // Reconstruct and store GLL points of the state derivatives
      reconStencil(stencil, gllPts, dom.doWeno, wenoRecon, to_derivX_gll, wenoIdl, wenoSigma);
      for (int ii=0; ii<tord; ii++) { derivDTs(l,0,ii) = gllPts(ii); }
    }
    for (int ii=0; ii<tord; ii++) {
      stateDTs(idR,0,ii) += dom.hyDensCells  (hs+k);
    }

    // Compute tord-1 time derivatives of the state, state spatial derivatives, mass flux,
    // energy flux, u RHS, v RHS, and w RHS using temporal Differential Transforms
    diffTransformEulerPrimX( stateDTs, derivDTs, utend, vtend, wtend, aderDerivX );

    // Compute the time-average and store into the zeroth time index
    timeAvg( stateDTs , dom );
    timeAvg( utend    , dom );
    timeAvg( vtend    , dom );
    timeAvg( wtend    , dom );

    // Compute the cell-centered tendencies for wind via quadrature because it
    // is using a high-order flux difference splitting
    // int( RHS , x , x_(i-1/2) , x_(i+1/2) )
    tend(idU,hs+k,hs+j,hs+i) = 0;
    tend(idV,hs+k,hs+j,hs+i) = 0;
    tend(idW,hs+k,hs+j,hs+i) = 0;
    for (int ii=0; ii<tord; ii++) {
      tend(idU,hs+k,hs+j,hs+i) += gllWts(ii) * utend(0,ii);
      tend(idV,hs+k,hs+j,hs+i) += gllWts(ii) * utend(0,ii);
      tend(idW,hs+k,hs+j,hs+i) += gllWts(ii) * utend(0,ii);
    }

    // Store the state vector in fwaves to compute fwaves from cell-interface state jumps
    for (int l=0; l<numState; l++) {
      stateLimits(l,1,k,j,i  ) = stateDTs(l,0,0     ); // Store the left cell edge state estimates
      stateLimits(l,0,k,j,i+1) = stateDTs(l,0,tord-1); // Store the Right cell edge state estimates
    }
  });

  //Reconcile the edge fluxes via MPI exchange.
  exch.haloInit      ();
  exch.edgePackN_x   (dom, stateLimits, numState);
  exch.edgeExchange_x(dom, par);
  exch.edgeUnpackN_x (dom, stateLimits, numState);

  // for (int k=0; k<dom.nz; k++) {
  //   for (int j=0; j<dom.ny; j++) {
  //     for (int i=0; i<dom.nx+1; i++) {
  yakl::parallel_for( dom.nz*dom.ny*(dom.nx+1) , YAKL_LAMBDA (int const iGlob) {
    int k, j, i;
    yakl::unpackIndices(iGlob,dom.nz,dom.ny,dom.nx+1,k,j,i);

    ////////////////////////////////////////////////////////////////////////////
    // Flux difference splitting for the wind updates
    ////////////////////////////////////////////////////////////////////////////
    // Compute the average state at the interface
    real r = 0.5_fp * ( stateLimits(idR,1,k,j,i) + stateLimits(idR,0,k,j,i)); // rho
    real u = 0.5_fp * ( stateLimits(idU,1,k,j,i) + stateLimits(idU,0,k,j,i)); // u
    real p = 0.5_fp * ( stateLimits(idT,1,k,j,i) + stateLimits(idT,0,k,j,i)); // p
    real cs2 = GAMMA*p/r;    // speed of sound squared
    real cs = sqrt(cs2);     // speed of sound
    // Compute the state jump across the interface
    real dr = 0.5_fp * ( stateLimits(idR,1,k,j,i) - stateLimits(idR,0,k,j,i));
    real du = 0.5_fp * ( stateLimits(idU,1,k,j,i) - stateLimits(idU,0,k,j,i));
    real dv = 0.5_fp * ( stateLimits(idV,1,k,j,i) - stateLimits(idV,0,k,j,i));
    real dw = 0.5_fp * ( stateLimits(idW,1,k,j,i) - stateLimits(idW,0,k,j,i));
    real dp = 0.5_fp * ( stateLimits(idT,1,k,j,i) - stateLimits(idT,0,k,j,i));
    // state at the left side of the interface
    real r1 = stateLimits(idR,0,k,j,i);
    real u1 = stateLimits(idU,0,k,j,i);
    real v1 = stateLimits(idV,0,k,j,i);
    real w1 = stateLimits(idW,0,k,j,i);
    real p1 = stateLimits(idP,0,k,j,i);
    // state at the right side of the interface
    real r2 = stateLimits(idR,1,k,j,i);
    real u2 = stateLimits(idU,1,k,j,i);
    real v2 = stateLimits(idV,1,k,j,i);
    real w2 = stateLimits(idW,1,k,j,i);
    real p2 = stateLimits(idP,1,k,j,i);
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
      // Wave 1 (u-cs): presumed always leftward propagating (no shocks)
      stateLimits(idU,0,k,j,i) += (-cs/r) * ( -r/(2*cs)*df(1) + df(4)/(2*cs2) );
      // Wave 2 (u+cs): presumed always rightward propagating (no shocks)
      stateLimits(idU,1,k,j,i) += ( cs/r) * (  r/(2*cs)*df(1) + df(4)/(2*cs2) );
      // Wave 3 does only affects density, so it's ignored
      // Waves 4 and 5 (u): 
      real w4 = df(2);
      real w5 = df(3);
      // If u > zero, it's rightward propagating, otherwise leftward
      // No need to worry about zero wind speed becaue then the wave is zero anyway
      if (u > 0) {
        stateLimits(idV,1,k,j,i) += w4;
        stateLimits(idW,1,k,j,i) += w5;
      } else {
        stateLimits(idV,0,k,j,i) += w4;
        stateLimits(idW,0,k,j,i) += w5;
      }
    }

    ////////////////////////////////////////////////////////////////////////////
    // Flux vector splitting for mass and energy
    ////////////////////////////////////////////////////////////////////////////
    // We can re-use the r, u, p, cs2, and cs calculated earlier
    // Store upwind state based on wind velocity
    real ru, uu, vu, wu, pu;
    if (u > 0) {
      ru = r1;  uu = u1;  vu = v1;  wu = w1;  pu = p1;
    } else {
      ru = r2;  uu = u2;  vu = v2;  wu = w2;  pu = p2;
    }
    // Now compute upwind characteristic variabes
    SArray<real,numState> chu;  // upwind characteristic variables
    chu(0) = -r/(2*cs)*u2 + p2/(2*cs2); // u-cs wave: assuming no shocks
    chu(1) =  r/(2*cs)*u1 + p1/(2*cs2); // u+cs wave: assuming no shocks
    chu(2) = ru - pu/cs2;
    chu(3) = vu;
    chu(4) = wu;
    // Now compute the upwind state based on upwind characteristic variables
    ru =     chu(0)   +     chu(1)   + chu(2);
    uu = -cs*chu(0)/r +  cs*chu(1)/r;
    vu =                                       chu(3);
    wu =                                               chu(4);
    pu = cs2*chu(0)   + cs2*chu(1)
    // Now compute the upwind flux based on the upwind state
    real keu = 0.5_fp*ru*(uu*uu+vu*vu+wu*wu); // upwind kinetic energy
    real reu = pu*CV/RD + keu;                // upwind rho*e
    flux_r (k,j,i) = ru*uu;                   // upwind mass flux
    flux_re(k,j,i) = uu*reu + uu*pu;          // upwind energy flux
  });

  // Apply the fwaves to the tendencies
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
    tend(idU,k,j,i) += - ( fwaves(idU,1,k,j,i) + fwaves(idU,0,k,j,i+1) ) / dom.dx;  // u tendency
    tend(idV,k,j,i) += - ( fwaves(idV,1,k,j,i) + fwaves(idV,0,k,j,i+1) ) / dom.dx;  // v tendency
    tend(idW,k,j,i) += - ( fwaves(idW,1,k,j,i) + fwaves(idW,0,k,j,i+1) ) / dom.dx;  // w tendency
  });
}


void Tendencies::compEulerTend_Y(realArr &state, Domain const &dom, Exchange &exch, Parallel const &par, realArr &tend) {
  auto &fwaves      = this->fwaves     ;
  auto &src         = this->src        ;
  auto &gllWts      = this->gllWts     ;
  auto &to_gll      = this->to_gll     ;
  auto &wenoRecon   = this->wenoRecon  ;
  auto &wenoIdl     = this->wenoIdl    ;
  auto &wenoSigma   = this->wenoSigma  ;
  auto &aderDerivY  = this->aderDerivY ;

  // Exchange halos in the y-direction
  exch.haloInit      ();
  exch.haloPackN_y   (dom, state, numState);
  exch.haloExchange_y(dom, par);
  exch.haloUnpackN_y (dom, state, numState);

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
      if (l != idT) {  // rho, u, v, and w
        for (int ii=0; ii<ord; ii++) {
          stencil(ii) = state(l,hs+k,j+ii,hs+i);
        }
      } else {               // pressure
        for (int ii=0; ii<ord; ii++) {
          real r  = state(idR,hs+k,j+ii,hs+i) + dom.hyDensCells  (hs+k);
          real u  = state(idU,hs+k,j+ii,hs+i);
          real v  = state(idV,hs+k,j+ii,hs+i);
          real w  = state(idW,hs+k,j+ii,hs+i);
          real re = state(idT,hs+k,j+ii,hs+i) + dom.hyEnergyCells(hs+k);
          real ke = 0.5_fp*r*(u*u+v*v+w*w);
          real p = RD/CV*(re-ke);
          stencil(ii) = p;
        }
      }

      // Reconstruct and store GLL points of the state values
      reconStencil(stencil, gllPts, dom.doWeno, wenoRecon, to_gll       , wenoIdl, wenoSigma);
      for (int ii=0; ii<tord; ii++) { stateDTs(l,0,ii) = gllPts(ii); }

      // Reconstruct and store GLL points of the state derivatives
      reconStencil(stencil, gllPts, dom.doWeno, wenoRecon, to_derivX_gll, wenoIdl, wenoSigma);
      for (int ii=0; ii<tord; ii++) { derivDTs(l,0,ii) = gllPts(ii); }
    }
    for (int ii=0; ii<tord; ii++) {
      stateDTs(idR,0,ii) += dom.hyDensCells  (hs+k);
    }

    // Compute tord-1 time derivatives of the state, state spatial derivatives, mass flux,
    // energy flux, u RHS, v RHS, and w RHS using temporal Differential Transforms
    diffTransformEulerPrimY( stateDTs, derivDTs, utend, vtend, wtend, aderDerivY );

    // Compute the time-average and store into the zeroth time index
    timeAvg( stateDTs , dom );
    timeAvg( utend    , dom );
    timeAvg( vtend    , dom );
    timeAvg( wtend    , dom );

    // Compute the cell-centered tendencies for wind via quadrature because it
    // is using a high-order flux difference splitting
    // int( RHS , x , x_(i-1/2) , x_(i+1/2) )
    tend(idU,hs+k,hs+j,hs+i) = 0;
    tend(idV,hs+k,hs+j,hs+i) = 0;
    tend(idW,hs+k,hs+j,hs+i) = 0;
    for (int ii=0; ii<tord; ii++) {
      tend(idU,hs+k,hs+j,hs+i) += gllWts(ii) * utend(0,ii);
      tend(idV,hs+k,hs+j,hs+i) += gllWts(ii) * utend(0,ii);
      tend(idW,hs+k,hs+j,hs+i) += gllWts(ii) * utend(0,ii);
    }

    // Store the state vector in fwaves to compute fwaves from cell-interface state jumps
    for (int l=0; l<numState; l++) {
      stateLimits(l,1,k,j,i  ) = stateDTs(l,0,0     ); // Store the left cell edge state estimates
      stateLimits(l,0,k,j,i+1) = stateDTs(l,0,tord-1); // Store the Right cell edge state estimates
    }
  });

  //Reconcile the edge fluxes via MPI exchange.
  exch.haloInit      ();
  exch.edgePackN_y   (dom, stateLimits, numState);
  exch.edgeExchange_y(dom, par);
  exch.edgeUnpackN_y (dom, stateLimits, numState);

  // for (int k=0; k<dom.nz; k++) {
  //   for (int j=0; j<dom.ny+1; j++) {
  //     for (int i=0; i<dom.nx; i++) {
  yakl::parallel_for( dom.nz*(dom.ny+1)*dom.nx , YAKL_LAMBDA (int const iGlob) {
    int k, j, i;
    yakl::unpackIndices(iGlob,dom.nz,dom.ny+1,dom.nx,k,j,i);

    ////////////////////////////////////////////////////////////////////////////
    // Flux difference splitting for the wind updates
    ////////////////////////////////////////////////////////////////////////////
    // Compute the average state at the interface
    real r = 0.5_fp * ( stateLimits(idR,1,k,j,i) + stateLimits(idR,0,k,j,i)); // rho
    real v = 0.5_fp * ( stateLimits(idV,1,k,j,i) + stateLimits(idV,0,k,j,i)); // v
    real p = 0.5_fp * ( stateLimits(idT,1,k,j,i) + stateLimits(idT,0,k,j,i)); // p
    real cs2 = GAMMA*p/r;    // speed of sound squared
    real cs = sqrt(cs2);     // speed of sound
    // Compute the state jump across the interface
    real dr = 0.5_fp * ( stateLimits(idR,1,k,j,i) - stateLimits(idR,0,k,j,i));
    real du = 0.5_fp * ( stateLimits(idU,1,k,j,i) - stateLimits(idU,0,k,j,i));
    real dv = 0.5_fp * ( stateLimits(idV,1,k,j,i) - stateLimits(idV,0,k,j,i));
    real dw = 0.5_fp * ( stateLimits(idW,1,k,j,i) - stateLimits(idW,0,k,j,i));
    real dp = 0.5_fp * ( stateLimits(idT,1,k,j,i) - stateLimits(idT,0,k,j,i));
    // state at the left side of the interface
    real r1 = stateLimits(idR,0,k,j,i);
    real u1 = stateLimits(idU,0,k,j,i);
    real v1 = stateLimits(idV,0,k,j,i);
    real w1 = stateLimits(idW,0,k,j,i);
    real p1 = stateLimits(idP,0,k,j,i);
    // state at the right side of the interface
    real r2 = stateLimits(idR,1,k,j,i);
    real u2 = stateLimits(idU,1,k,j,i);
    real v2 = stateLimits(idV,1,k,j,i);
    real w2 = stateLimits(idW,1,k,j,i);
    real p2 = stateLimits(idP,1,k,j,i);
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
      // Wave 1 (u-cs): presumed always leftward propagating (no shocks)
      stateLimits(idV,0,k,j,i) += (-cs/r) * ( -r/(2*cs)*df(2) + df(4)/(2*cs2) );
      // Wave 2 (u+cs): presumed always rightward propagating (no shocks)
      stateLimits(idV,1,k,j,i) += ( cs/r) * (  r/(2*cs)*df(2) + df(4)/(2*cs2) );
      // Wave 3 does only affects density, so it's ignored
      // Waves 4 and 5 (u): 
      real w4 = df(1);
      real w5 = df(3);
      // If u > zero, it's rightward propagating, otherwise leftward
      // No need to worry about zero wind speed becaue then the wave is zero anyway
      if (u > 0) {
        stateLimits(idU,1,k,j,i) += w4;
        stateLimits(idW,1,k,j,i) += w5;
      } else {
        stateLimits(idU,0,k,j,i) += w4;
        stateLimits(idW,0,k,j,i) += w5;
      }
    }

    ////////////////////////////////////////////////////////////////////////////
    // Flux vector splitting for mass and energy
    ////////////////////////////////////////////////////////////////////////////
    // We can re-use the r, u, p, cs2, and cs calculated earlier
    // Store upwind state based on wind velocity
    real ru, uu, vu, wu, pu;
    if (u > 0) {
      ru = r1;  uu = u1;  vu = v1;  wu = w1;  pu = p1;
    } else {
      ru = r2;  uu = u2;  vu = v2;  wu = w2;  pu = p2;
    }
    // Now compute upwind characteristic variabes
    SArray<real,numState> chu;  // upwind characteristic variables
    chu(0) = -r/(2*cs)*v2 + p2/(2*cs2); // u-cs wave: assuming no shocks
    chu(1) =  r/(2*cs)*v1 + p1/(2*cs2); // u+cs wave: assuming no shocks
    chu(2) = ru - pu/cs2;
    chu(3) = uu;
    chu(4) = wu;
    // Now compute the upwind state based on upwind characteristic variables
    ru =     chu(0)   +     chu(1)   + chu(2);
    uu =                                       chu(3);
    vu = -cs*chu(0)/r +  cs*chu(1)/r;
    wu =                                               chu(4);
    pu = cs2*chu(0)   + cs2*chu(1)
    // Now compute the upwind flux based on the upwind state
    real keu = 0.5_fp*ru*(uu*uu+vu*vu+wu*wu); // upwind kinetic energy
    real reu = pu*CV/RD + keu;                // upwind rho*e
    flux_r (k,j,i) = ru*vu;                   // upwind mass flux
    flux_re(k,j,i) = vu*reu + vu*pu;          // upwind energy flux
  });

  // Apply the fwaves to the tendencies
  //   for (int k=0; k<dom.nz; k++) {
  //     for (int j=0; j<dom.ny; j++) {
  //       for (int i=0; i<dom.nx; i++) {
  yakl::parallel_for( dom.nz*dom.ny*dom.nx , YAKL_LAMBDA (int const iGlob) {
    int l, k, j, i;
    yakl::unpackIndices(iGlob,dom.nz,dom.ny,dom.nx,k,j,i);
    // Flux vector form update for mass and energy
    tend(idR,k,j,i)  = - ( flux_r (k,j,i+1) - flux_r (k,j,i) ) / dom.dy;  // mass tendency
    tend(idT,k,j,i)  = - ( flux_re(k,j,i+1) - flux_re(k,j,i) ) / dom.dy;  // energy tendency
    // Flux difference splitting form update for velocities
    tend(idU,k,j,i) += - ( fwaves(idU,1,k,j,i) + fwaves(idU,0,k,j,i+1) ) / dom.dy;  // u tendency
    tend(idV,k,j,i) += - ( fwaves(idV,1,k,j,i) + fwaves(idV,0,k,j,i+1) ) / dom.dy;  // v tendency
    tend(idW,k,j,i) += - ( fwaves(idW,1,k,j,i) + fwaves(idW,0,k,j,i+1) ) / dom.dy;  // w tendency
  });
}



void Tendencies::compEulerTend_Z(realArr &state, Domain const &dom, Exchange &exch, Parallel const &par, realArr &tend) {
  auto &fwaves      = this->fwaves     ;
  auto &src         = this->src        ;
  auto &gllWts      = this->gllWts     ;
  auto &to_gll      = this->to_gll     ;
  auto &wenoRecon   = this->wenoRecon  ;
  auto &wenoIdl     = this->wenoIdl    ;
  auto &wenoSigma   = this->wenoSigma  ;
  auto &aderDerivZ  = this->aderDerivZ ;

  // Apply boundaries in the z-direction
  stateBoundariesZ( state , dom );

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
      if (l == idR) {         // rho
        for (int ii=0; ii<ord; ii++) {
          stencil(ii) = state(idR,k+ii,hs+j,hs+i) + dom.hyDensCells(k+ii);
        }
      } else if (l != idT) {  // u, v, and w
        for (int ii=0; ii<ord; ii++) {
          stencil(ii) = state(l,k+ii,hs+j,hs+i);
        }
      } else {                // pressure
        for (int ii=0; ii<ord; ii++) {
          real r  = state(idR,k+ii,hs+j,hs+i) + dom.hyDensCells  (k+ii);
          real u  = state(idU,k+ii,hs+j,hs+i);
          real v  = state(idV,k+ii,hs+j,hs+i);
          real w  = state(idW,k+ii,hs+j,hs+i);
          real re = state(idT,k+ii,hs+j,hs+i) + dom.hyEnergyCells(k+ii);
          real ke = 0.5_fp*r*(u*u+v*v+w*w);
          real p = RD/CV*(re-ke);
          stencil(ii) = p;
        }
      }

      // Reconstruct and store GLL points of the state values
      reconStencil(stencil, gllPts, dom.doWeno, wenoRecon, to_gll       , wenoIdl, wenoSigma);
      for (int ii=0; ii<tord; ii++) { stateDTs(l,0,ii) = gllPts(ii); }

      // Reconstruct and store GLL points of the state derivatives
      reconStencil(stencil, gllPts, dom.doWeno, wenoRecon, to_derivX_gll, wenoIdl, wenoSigma);
      for (int ii=0; ii<tord; ii++) { derivDTs(l,0,ii) = gllPts(ii); }
    }

    // Compute tord-1 time derivatives of the state, state spatial derivatives, mass flux,
    // energy flux, u RHS, v RHS, and w RHS using temporal Differential Transforms
    diffTransformEulerPrimY( stateDTs, derivDTs, utend, vtend, wtend, dph, aderDerivZ );

    // Compute the time-average and store into the zeroth time index
    timeAvg( stateDTs , dom );
    timeAvg( utend    , dom );
    timeAvg( vtend    , dom );
    timeAvg( wtend    , dom );

    // Compute the cell-centered tendencies for wind via quadrature because it
    // is using a high-order flux difference splitting
    // int( RHS , x , x_(i-1/2) , x_(i+1/2) )
    tend(idU,hs+k,hs+j,hs+i) = 0;
    tend(idV,hs+k,hs+j,hs+i) = 0;
    tend(idW,hs+k,hs+j,hs+i) = 0;
    for (int ii=0; ii<tord; ii++) {
      tend(idU,hs+k,hs+j,hs+i) += gllWts(ii) * utend(0,ii);
      tend(idV,hs+k,hs+j,hs+i) += gllWts(ii) * utend(0,ii);
      tend(idW,hs+k,hs+j,hs+i) += gllWts(ii) * utend(0,ii);
    }

    // Store the state vector in fwaves to compute fwaves from cell-interface state jumps
    for (int l=0; l<numState; l++) {
      stateLimits(l,1,k,j,i  ) = stateDTs(l,0,0     ); // Store the left cell edge state estimates
      stateLimits(l,0,k,j,i+1) = stateDTs(l,0,tord-1); // Store the Right cell edge state estimates
    }
  });

  // Apply boundaries for fluxes in the z-direction
  edgeBoundariesZ( stateLimits , dom );

  // for (int k=0; k<dom.nz; k++) {
  //   for (int j=0; j<dom.ny; j++) {
  //     for (int i=0; i<dom.nx+1; i++) {
  yakl::parallel_for( dom.nz*dom.ny*(dom.nx+1) , YAKL_LAMBDA (int const iGlob) {
    int k, j, i;
    yakl::unpackIndices(iGlob,dom.nz,dom.ny,dom.nx+1,k,j,i);

    ////////////////////////////////////////////////////////////////////////////
    // Flux difference splitting for the wind updates
    ////////////////////////////////////////////////////////////////////////////
    // Compute the average state at the interface
    real r = 0.5_fp * ( stateLimits(idR,1,k,j,i) + stateLimits(idR,0,k,j,i)); // rho
    real w = 0.5_fp * ( stateLimits(idW,1,k,j,i) + stateLimits(idW,0,k,j,i)); // w
    real p = 0.5_fp * ( stateLimits(idT,1,k,j,i) + stateLimits(idT,0,k,j,i)); // p
    real cs2 = GAMMA*p/r;    // speed of sound squared
    real cs = sqrt(cs2);     // speed of sound
    // Compute the state jump across the interface
    real dr = 0.5_fp * ( stateLimits(idR,1,k,j,i) - stateLimits(idR,0,k,j,i));
    real du = 0.5_fp * ( stateLimits(idU,1,k,j,i) - stateLimits(idU,0,k,j,i));
    real dv = 0.5_fp * ( stateLimits(idV,1,k,j,i) - stateLimits(idV,0,k,j,i));
    real dw = 0.5_fp * ( stateLimits(idW,1,k,j,i) - stateLimits(idW,0,k,j,i));
    real dp = 0.5_fp * ( stateLimits(idT,1,k,j,i) - stateLimits(idT,0,k,j,i));
    // state at the left side of the interface
    real r1 = stateLimits(idR,0,k,j,i);
    real u1 = stateLimits(idU,0,k,j,i);
    real v1 = stateLimits(idV,0,k,j,i);
    real w1 = stateLimits(idW,0,k,j,i);
    real p1 = stateLimits(idP,0,k,j,i);
    // state at the right side of the interface
    real r2 = stateLimits(idR,1,k,j,i);
    real u2 = stateLimits(idU,1,k,j,i);
    real v2 = stateLimits(idV,1,k,j,i);
    real w2 = stateLimits(idW,1,k,j,i);
    real p2 = stateLimits(idP,1,k,j,i);
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
      // Wave 1 (u-cs): presumed always leftward propagating (no shocks)
      stateLimits(idV,0,k,j,i) += (-cs/r) * ( -r/(2*cs)*df(2) + df(4)/(2*cs2) );
      // Wave 2 (u+cs): presumed always rightward propagating (no shocks)
      stateLimits(idV,1,k,j,i) += ( cs/r) * (  r/(2*cs)*df(2) + df(4)/(2*cs2) );
      // Wave 3 does only affects density, so it's ignored
      // Waves 4 and 5 (u): 
      real w4 = df(1);
      real w5 = df(3);
      // If u > zero, it's rightward propagating, otherwise leftward
      // No need to worry about zero wind speed becaue then the wave is zero anyway
      if (u > 0) {
        stateLimits(idU,1,k,j,i) += w4;
        stateLimits(idW,1,k,j,i) += w5;
      } else {
        stateLimits(idU,0,k,j,i) += w4;
        stateLimits(idW,0,k,j,i) += w5;
      }
    }

    ////////////////////////////////////////////////////////////////////////////
    // Flux vector splitting for mass and energy
    ////////////////////////////////////////////////////////////////////////////
    // We can re-use the r, u, p, cs2, and cs calculated earlier
    // Store upwind state based on wind velocity
    real ru, uu, vu, wu, pu;
    if (u > 0) {
      ru = r1;  uu = u1;  vu = v1;  wu = w1;  pu = p1;
    } else {
      ru = r2;  uu = u2;  vu = v2;  wu = w2;  pu = p2;
    }
    // Now compute upwind characteristic variabes
    SArray<real,numState> chu;  // upwind characteristic variables
    chu(0) = -r/(2*cs)*v2 + p2/(2*cs2); // u-cs wave: assuming no shocks
    chu(1) =  r/(2*cs)*v1 + p1/(2*cs2); // u+cs wave: assuming no shocks
    chu(2) = ru - pu/cs2;
    chu(3) = uu;
    chu(4) = wu;
    // Now compute the upwind state based on upwind characteristic variables
    ru =     chu(0)   +     chu(1)   + chu(2);
    uu =                                       chu(3);
    vu = -cs*chu(0)/r +  cs*chu(1)/r;
    wu =                                               chu(4);
    pu = cs2*chu(0)   + cs2*chu(1)
    // Now compute the upwind flux based on the upwind state
    real keu = 0.5_fp*ru*(uu*uu+vu*vu+wu*wu); // upwind kinetic energy
    real reu = pu*CV/RD + keu;                // upwind rho*e
    flux_r (k,j,i) = ru*vu;                   // upwind mass flux
    flux_re(k,j,i) = vu*reu + vu*pu;          // upwind energy flux
  });

  // Apply the fwaves to the tendencies
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
    tend(idU,k,j,i) += - ( fwaves(idU,1,k,j,i) + fwaves(idU,0,k,j,i+1) ) / dom.dx;  // u tendency
    tend(idV,k,j,i) += - ( fwaves(idV,1,k,j,i) + fwaves(idV,0,k,j,i+1) ) / dom.dx;  // v tendency
    tend(idW,k,j,i) += - ( fwaves(idW,1,k,j,i) + fwaves(idW,0,k,j,i+1) ) / dom.dx;  // w tendency
  });
}



void Tendencies::compEulerTend_S(realArr &state, Domain const &dom, Exchange &exch, Parallel const &par, realArr &tend) {
  // Form the tendencies
  // for (int k=0; k<dom.nz; k++) {
  //   for (int j=0; j<dom.ny; j++) {
  //     for (int i=0; i<dom.nx; i++) {
  yakl::parallel_for( dom.nz*dom.ny*dom.nx , YAKL_LAMBDA (int const iGlob) {
    int k, j, i;
    yakl::unpackIndices(iGlob,dom.nz,dom.ny,dom.nx,k,j,i);
    tend(idR ,k,j,i) = 0;
    tend(idU ,k,j,i) = 0;
    tend(idV ,k,j,i) = 0;
    tend(idW ,k,j,i) = -state(idR,hs+k,hs+j,hs+i) / (state(idR,hs+k,hs+j,hs+i) + dom.hyDensCells(hs+k)) * GRAV;
    tend(idRT,k,j,i) = 0;
  });
}



void Tendencies::compStrakaTend(realArr &state, Domain const &dom, Exchange &exch, Parallel const &par, realArr &tend) {
  // //Exchange halos in the x-direction
  // exch.haloInit      ();
  // exch.haloPackN_x   (dom, state, numState);
  // exch.haloExchange_x(dom, par);
  // exch.haloUnpackN_x (dom, state, numState);

  // //Exchange halos in the y-direction
  // exch.haloInit      ();
  // exch.haloPackN_x   (dom, state, numState);
  // exch.haloExchange_x(dom, par);
  // exch.haloUnpackN_x (dom, state, numState);
  // 
  // // Boundaries for the fluid state in the z-direction
  // stateBoundariesZ(state, dom);
}


void TendenciesThetaPrimADER::stateBoundariesZ(realArr &state, Domain const &dom) {
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


void TendenciesThetaPrimADER::edgeBoundariesZ(realArr &stateLimits, Domain const &dom) {
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



