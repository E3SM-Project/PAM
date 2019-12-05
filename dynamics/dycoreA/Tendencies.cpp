
#include "Tendencies.h"



void Tendencies::initialize(Domain const &dom) {
  TransformMatrices<real> trans;

  fluxLimits  = realArr("fluxLimits",numState,2,dom.nz+1,dom.ny+1,dom.nx+1);
  stateLimits = realArr("srcLimits" ,numState,2,dom.nz+1,dom.ny+1,dom.nx+1);
  flux        = realArr("flux"      ,numState  ,dom.nz+1,dom.ny+1,dom.nx+1);
  src         = realArr("src"       ,dom.nz,dom.ny,dom.nx);

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
    SArray<real,numState,tord,tord> fluxDTs;   // GLL flux DTs (var,time,space)

    // Compute tord GLL points of the fluid state and spatial derivative
    for (int l=0; l<numState; l++) {
      SArray<real,ord> stencil;
      SArray<real,tord> gllPts;
      // Store the stencil values
      for (int ii=0; ii<ord; ii++) { stencil(ii) = state(l,hs+k,hs+j,i+ii); }
      //Reconstruct and store GLL points of the state values
      reconStencil(stencil, gllPts, dom.doWeno, wenoRecon, to_gll, wenoIdl, wenoSigma);
      for (int ii=0; ii<tord; ii++) { stateDTs(l,0,ii) = gllPts(ii); }
    }
    for (int ii=0; ii<tord; ii++) {
      stateDTs(idR,0,ii) += dom.hyDensCells    (hs+k);
      stateDTs(idP,0,ii) += dom.hyPressureCells(hs+k);
    }

    // Compute tord-1 time derivatives of the state and spatial derivs
    // using temporal Differential Transforms
    diffTransformEulerPrimX( stateDTs , fluxDTs , aderDerivX );

    // Compute the temporal time-average and store into the zeroth time index
    timeAvg( stateDTs , dom );
    timeAvg( fluxDTs  , dom );

    // Store the state vector in fwaves to compute fwaves from cell-interface state jumps
    for (int l=0; l<numState; l++) {
      // Store the left cell edge state and flux estimates
      stateLimits(l,1,k,j,i  ) = stateDTs(l,0,0);
      fluxLimits (l,1,k,j,i  ) = fluxDTs (l,0,0);

      // Store the Right cell edge state and flux estimates
      stateLimits(l,0,k,j,i+1) = stateDTs(l,0,tord-1);
      fluxLimits (l,0,k,j,i+1) = fluxDTs (l,0,tord-1);
    }
  });

  //Reconcile the edge fluxes via MPI exchange.
  exch.haloInit      ();
  exch.edgePackN_x   (dom, stateLimits, numState);
  exch.edgePackN_x   (dom, fluxLimits , numState);
  exch.edgeExchange_x(dom, par);
  exch.edgeUnpackN_x (dom, stateLimits, numState);
  exch.edgeUnpackN_x (dom, fluxLimits , numState);

  // Compute the fwaves from the cell interface jumps
  // for (int k=0; k<dom.nz; k++) {
  //   for (int j=0; j<dom.ny; j++) {
  //     for (int i=0; i<dom.nx+1; i++) {
  yakl::parallel_for( dom.nz*dom.ny*(dom.nx+1) , YAKL_LAMBDA (int const iGlob) {
    SArray<real,numState> s1, s2, f1, f2, upw;
    for (int l=0; l<numState; l++) {
      s1(l) = stateLimits(l,0,k,j,i);
      s2(l) = stateLimits(l,1,k,j,i);
      f1(l) = fluxLimits (l,0,k,j,i);
      f2(l) = fluxLimits (l,1,k,j,i);
    }
    riemannX(s1, s2, f1, f2, upw);
    for (int l=0; l<numState; l++) {
      flux(l,k,j,i) = upw(l);
    }
  });

  // Apply the fwaves to the tendencies
  // for (int l=0; l<numState; l++) {
  //   for (int k=0; k<dom.nz; k++) {
  //     for (int j=0; j<dom.ny; j++) {
  //       for (int i=0; i<dom.nx; i++) {
  yakl::parallel_for( numState*dom.nz*dom.ny*dom.nx , YAKL_LAMBDA (int const iGlob) {
    int l, k, j, i;
    yakl::unpackIndices(iGlob,numState,dom.nz,dom.ny,dom.nx,l,k,j,i);
    tend(l,k,j,i) += - ( fwaves(l,1,k,j,i) + fwaves(l,0,k,j,i+1) ) / dom.dx;
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
    SArray<real,numState,tord,tord> derivDTs;  // GLL state spatial derivative DTs (var,time,space)

    // Compute tord GLL points of the fluid state and spatial derivative
    for (int l=0; l<numState; l++) {
      SArray<real,ord> stencil;
      SArray<real,tord> gllPts;
      // Store the stencil values
      for (int ii=0; ii<ord; ii++) { stencil(ii) = state(l,hs+k,j+ii,hs+i); }

      //Reconstruct and store GLL points of the state values
      reconStencil(stencil, gllPts, dom.doWeno, wenoRecon, to_gll, wenoIdl, wenoSigma);
      for (int ii=0; ii<tord; ii++) { stateDTs(l,0,ii) = gllPts(ii); }

      //Reconstruct and store GLL points of the state spatial derivative values
      reconStencil(stencil, gllPts, dom.doWeno, wenoRecon, to_derivY_gll, wenoIdl, wenoSigma);
      for (int ii=0; ii<tord; ii++) { derivDTs(l,0,ii) = gllPts(ii); }
    }
    for (int ii=0; ii<tord; ii++) {
      stateDTs(idR ,0,ii) += dom.hyDensCells     (hs+k);
      stateDTs(idRT,0,ii) += dom.hyDensThetaCells(hs+k);
    }

    // Compute tord-1 time derivatives of the state and spatial derivs
    // using temporal Differential Transforms
    diffTransformEulerPrimY( stateDTs , derivDTs , aderDerivY );

    // Compute the temporal time-average and store into the zeroth time index
    timeAvg( stateDTs );
    timeAvg( derivDTs );

    // Compute tend = -int( (dg/dq)*(dq/dy) , y , -dy/2 , dy/2 ) / dy 
    for (int l=0; l<numState; l++) { tend(l,k,j,i) = 0; }
    for (int ii=0; ii<tord; ii++) {
      real qw  = gllWts(ii);
      real drdy  = derivDTs(idR ,0,ii);
      real dudy  = derivDTs(idU ,0,ii);
      real dvdy  = derivDTs(idV ,0,ii);
      real dwdy  = derivDTs(idW ,0,ii);
      real drtdy = derivDTs(idRT,0,ii);
      real r   = stateDTs(idR ,0,ii);
      real v   = stateDTs(idV ,0,ii);
      real t   = stateDTs(idRT,0,ii) / r;
      real p   = C0*pow(r*t,GAMMA);
      real cs2 = GAMMA*p/r;
      tend(idR ,k,j,i) += -qw*( v*drdy          + r  *dvdy                          ) / dom.dy;
      tend(idU ,k,j,i) += -qw*(          v*dudy                                     ) / dom.dy;
      tend(idV ,k,j,i) += -qw*(                 + v  *dvdy          + cs2/r/t*drtdy ) / dom.dy;
      tend(idW ,k,j,i) += -qw*(                            + v*dwdy                 ) / dom.dy;
      tend(idRT,k,j,i) += -qw*(                   r*t*dvdy          + v      *drtdy ) / dom.dy;
    }

    // Store the state vector in fwaves to compute fwaves from cell-interface state jumps
    for (int l=0; l<numState; l++) {
      fwaves(l,1,k,j  ,i) = stateDTs(l,0,0     );
      fwaves(l,0,k,j+1,i) = stateDTs(l,0,tord-1);
    }
  });

  // Reconcile the edge state via MPI exchange.
  exch.haloInit      ();
  exch.edgePackN_y   (dom, fwaves, numState);
  exch.edgeExchange_y(dom, par);
  exch.edgeUnpackN_y (dom, fwaves, numState);

  // Compute the fwaves from the cell interface jumps
  // for (int k=0; k<dom.nz; k++) {
  //   for (int j=0; j<dom.ny; j++) {
  //     for (int i=0; i<dom.nx+1; i++) {
  yakl::parallel_for( dom.nz*(dom.ny+1)*dom.nx) , YAKL_LAMBDA (int const iGlob) {
    int k, j, i;
    yakl::unpackIndices(iGlob,dom.nz,dom.ny+1,dom.nx,k,j,i);
    // Compute averaged values for the flux Jacobian diagonalization
    real r   = 0.5_fp * ( fwaves(idR,0,k,j,i) + fwaves(idR,1,k,j,i) );
    real v   = 0.5_fp * ( fwaves(idV,0,k,j,i) + fwaves(idV,1,k,j,i) );
    real t   = 0.5_fp * ( fwaves(idT,0,k,j,i) + fwaves(idT,1,k,j,i) );
    real p   = C0*pow(r*t,GAMMA);
    real cs2 = GAMMA*p/r;
    real cs  = sqrt(cs2);

    // Compute the state jump over the interface
    SArray<real,numState> dq;
    for (int l=0; l<numState; l++) {
      dq(l) = fwaves(l,1,k,j,i) - fwaves(l,0,k,j,i);
    }

    // Compute df = A*dq
    SArray<real,numState> df;
    df(idR ) = v*dq(idR)             + r  *dq(idV)                       dq(idRT);
    df(idU ) =             v*dq(idU)                                     dq(idRT);
    df(idV ) =                       + v  *dq(idV)             + cs2/r/t*dq(idRT);
    df(idW ) =                                     + v*dq(idW)           dq(idRT);
    df(idRT) =                         r*t*dq(idV)             + v      *dq(idRT);

    // Compute characteristic variables (L*df)
    SArray<real,numState> ch;
    ch(0) = -r/(2*cs)*df(idV) + df(idRT)/(2*t);
    ch(1) =  r/(2*cs)*df(idV) + df(idRT)/(2*t);
    ch(2) = df(idR)           - df(idRT)/t;
    ch(3) = df(idU);
    ch(4) = df(idW);

    // Compute fwaves
    for (int l=0; l<numState; l++) {
      fwaves(l,0,k,j,i) = 0;
      fwaves(l,1,k,j,i) = 0;
    }

    // First wave (u-cs); always negative wave speed
    fwaves(idR ,0,k,j,i) +=       ch(0);
    fwaves(idV ,0,k,j,i) += -cs/r*ch(0);
    fwaves(idRT,0,k,j,i) += t    *ch(0);

    // Second wave (u+cs); always positive wave speed
    fwaves(idR ,1,k,j,i) +=       ch(1);
    fwaves(idV ,1,k,j,i) +=  cs/r*ch(1);
    fwaves(idRT,1,k,j,i) += t    *ch(1);

    if (u > 0) {
      // Third wave
      fwaves(idR,1,k,j,i) += ch(2);
      // Fourth wave
      fwaves(idU,1,k,j,i) += ch(3);
      // Fifth Wave
      fwaves(idW,1,k,j,i) += ch(4);
    } else {
      // Third wave
      fwaves(idR,0,k,j,i) += ch(2);
      // Fourth wave
      fwaves(idU,0,k,j,i) += ch(3);
      // Fifth Wave
      fwaves(idW,0,k,j,i) += ch(4);
    }
  });

  // Apply the fwaves to the tendencies
  // for (int l=0; l<numState; l++) {
  //   for (int k=0; k<dom.nz; k++) {
  //     for (int j=0; j<dom.ny; j++) {
  //       for (int i=0; i<dom.nx; i++) {
  yakl::parallel_for( numState*dom.nz*dom.ny*dom.nx , YAKL_LAMBDA (int const iGlob) {
    int l, k, j, i;
    yakl::unpackIndices(iGlob,numState,dom.nz,dom.ny,dom.nx,l,k,j,i);
    tend(l,k,j,i) += - ( fwaves(l,1,k,j,i) + fwaves(l,0,k,j+1,i) ) / dom.dy;
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

  // Pre-compute pressure perturbation
  // for (int k=0; k<dom.nz; k++) {
  //   for (int j=0; j<dom.ny; j++) {
  //     for (int i=0; i<dom.nx; i++) {
  yakl::parallel_for( dom.nz*dom.ny*dom.nx , YAKL_LAMBDA ( int const iGlob ) {
    int k, j, i;
    yakl::unpackIndices(iGlob,dom.nz,dom.ny,dom.nx,k,j,i);
    real rt = state(idRT,hs+k,hs+j,hs+i) + dom.hyDensThetaCells(hs+k);
    pressure(hs+k,j,i) = C0*pow(rt,GAMMA) - dom.hyPressureCells(hs+k);
  });

  // Exchange halos in the y-direction
  stateBoundariesZ( state , dom );

  // for (int k=0; k<dom.nz; k++) {
  //   for (int j=0; j<dom.ny; j++) {
  //     for (int i=0; i<dom.nx; i++) {
  yakl::parallel_for( dom.nz*dom.ny*dom.nx , YAKL_LAMBDA ( int const iGlob ) {
    int k, j, i;
    yakl::unpackIndices(iGlob,dom.nz,dom.ny,dom.nx,k,j,i);
    SArray<real,numState+1,tord,tord> stateDTs;  // GLL state DTs (var,time,space)
    SArray<real,numState+1,tord,tord> derivDTs;  // GLL state spatial derivative DTs (var,time,space)

    // Compute tord GLL points of the fluid state and spatial derivative
    for (int l=0; l<numState+1; l++) {
      SArray<real,ord> stencil;
      SArray<real,tord> gllPts;
      // Store the stencil values
      if (l < numState) {
        for (int ii=0; ii<ord; ii++) { stencil(ii) = state(l,k+ii,hs+j,hs+i); }
      } else {
        for (int ii=0; ii<ord; ii++) { stencil(ii) = pressure(k+ii); }
      }

      //Reconstruct and store GLL points of the state values
      reconStencil(stencil, gllPts, dom.doWeno, wenoRecon, to_gll, wenoIdl, wenoSigma);
      for (int ii=0; ii<tord; ii++) { stateDTs(l,0,ii) = gllPts(ii); }

      //Reconstruct and store GLL points of the state spatial derivative values
      reconStencil(stencil, gllPts, dom.doWeno, wenoRecon, to_derivZ_gll, wenoIdl, wenoSigma);
      for (int ii=0; ii<tord; ii++) { derivDTs(l,0,ii) = gllPts(ii); }
    }
    for (int ii=0; ii<tord; ii++) {
      stateDTs(idR ,0,ii) += dom.hyDensGLL     (k,ii);
      stateDTs(idRT,0,ii) += dom.hyDensThetaGLL(k,ii);
    }

    // Compute tord-1 time derivatives of the state and spatial derivs
    // using temporal Differential Transforms
    diffTransformEulerPrimZ( stateDTs , derivDTs , aderDerivZ );

    // Compute the temporal time-average and store into the zeroth time index
    timeAvg( stateDTs );
    timeAvg( derivDTs );

    // Compute tend = -int( (dh/dq)*(dq/dz) , z , -dz/2 , dz/2 ) / dz 
    for (int l=0; l<numState; l++) { tend(l,k,j,i) = 0; }
    for (int ii=0; ii<tord; ii++) {
      real qw  = gllWts(ii);
      real drdz  = derivDTs(idR ,0,ii);
      real dudz  = derivDTs(idU ,0,ii);
      real dvdz  = derivDTs(idV ,0,ii);
      real dwdz  = derivDTs(idW ,0,ii);
      real drtdz = derivDTs(idRT,0,ii);
      real dpdz  = derivDTs(idP ,0,ii);
      real r   = stateDTs(idR ,0,ii);
      real w   = stateDTs(idW ,0,ii);
      real rt  = stateDTs(idRT,0,ii);
      tend(idR ,k,j,i) += -qw*( w*drdz                 + r*dwdz                    ) / dom.dz;
      tend(idU ,k,j,i) += -qw*(        w*dudz                                      ) / dom.dz;
      tend(idV ,k,j,i) += -qw*(               + w*dvdz                             ) / dom.dz;
      tend(idW ,k,j,i) += -qw*(                          w*dwdz           + dpdz/r ) / dom.dz;
      tend(idRT,k,j,i) += -qw*(                         rt*dwdz + w*drtdz          ) / dom.dz;
    }

    // Store the state vector in fwaves to compute fwaves from cell-interface state jumps
    for (int l=0; l<numState+1; l++) {
      fwaves(l,1,k  ,j,i) = stateDTs(l,0,0     );
      fwaves(l,0,k+1,j,i) = stateDTs(l,0,tord-1);
    }
  });

  edgeBoundariesZ( fwaves , dom );

  // Compute the fwaves from the cell interface jumps
  // for (int k=0; k<dom.nz; k++) {
  //   for (int j=0; j<dom.ny; j++) {
  //     for (int i=0; i<dom.nx+1; i++) {
  yakl::parallel_for( (dom.nz+1)*dom.ny*dom.nx) , YAKL_LAMBDA (int const iGlob) {
    int k, j, i;
    yakl::unpackIndices(iGlob,dom.nz+1,dom.ny,dom.nx,k,j,i);
    // Compute averaged values for the flux Jacobian diagonalization
    real r   = 0.5_fp * ( fwaves(idR,0,k,j,i) + fwaves(idR,1,k,j,i) );
    real w   = 0.5_fp * ( fwaves(idW,0,k,j,i) + fwaves(idW,1,k,j,i) );
    real t   = 0.5_fp * ( fwaves(idT,0,k,j,i) + fwaves(idT,1,k,j,i) );
    real p   = C0*pow(r*t,GAMMA);
    real cs2 = GAMMA*p/r;
    real cs  = sqrt(cs2);

    // Compute the state jump over the interface
    SArray<real,numState+1> dq;
    for (int l=0; l<numState+1; l++) {
      dq(l) = fwaves(l,1,k,j,i) - fwaves(l,0,k,j,i);
    }

    // Compute df = A*dq
    SArray<real,numState+1> df;
    df(idR ) = w*dq(idR)                     + r    *dq(idW)                               ;
    df(idU ) =           w*dq(idU)                                                         ;
    df(idV ) =                     w*dq(idV)                                               ;
    df(idW ) =                                 w    *dq(idW)                    + dq(idP)/r;
    df(idRT) =                                 r*t  *dq(idW) + w      *dq(idRT)            ;
    df(idP ) =                                 cs2*r*dq(idW) + cs2*w/t*dq(idRT)            ;

    // Compute characteristic variables (L*df)
    SArray<real,numState> ch;
    ch(0) = -r/(2*cs)*df(idV) - w*df(idRT)/(2*(cs*t-t*w)) + df(idP)/(2*(cs2-cs*w));
    ch(1) =  r/(2*cs)*df(idV) + w*df(idRT)/(2*(cs*t+t*w)) + df(idP)/(2*(cs2+cs*w));
    ch(2) = df(idR)           -   df(idRT)/t;
    ch(3) = df(idU);
    ch(4) = df(idV);

    // Compute fwaves
    for (int l=0; l<numState; l++) {
      fwaves(l,0,k,j,i) = 0;
      fwaves(l,1,k,j,i) = 0;
    }

    // First wave (w-cs); always negative wave speed
    fwaves(idR ,0,k,j,i) +=       ch(0);
    fwaves(idW ,0,k,j,i) += -cs/r*ch(0);
    fwaves(idRT,0,k,j,i) += t    *ch(0);

    // Second wave (w+cs); always positive wave speed
    fwaves(idR ,1,k,j,i) +=       ch(1);
    fwaves(idW ,1,k,j,i) +=  cs/r*ch(1);
    fwaves(idRT,1,k,j,i) += t    *ch(1);

    if (u > 0) {
      // Third wave
      fwaves(idR,1,k,j,i) += ch(2);
      // Fourth wave
      fwaves(idU,1,k,j,i) += ch(3);
      // Fifth Wave
      fwaves(idV,1,k,j,i) += ch(4);
    } else {
      // Third wave
      fwaves(idR,0,k,j,i) += ch(2);
      // Fourth wave
      fwaves(idU,0,k,j,i) += ch(3);
      // Fifth Wave
      fwaves(idV,0,k,j,i) += ch(4);
    }
  });

  // Apply the fwaves to the tendencies
  // for (int l=0; l<numState; l++) {
  //   for (int k=0; k<dom.nz; k++) {
  //     for (int j=0; j<dom.ny; j++) {
  //       for (int i=0; i<dom.nx; i++) {
  yakl::parallel_for( numState*dom.nz*dom.ny*dom.nx , YAKL_LAMBDA (int const iGlob) {
    int l, k, j, i;
    yakl::unpackIndices(iGlob,numState,dom.nz,dom.ny,dom.nx,l,k,j,i);
    tend(l,k,j,i) += - ( fwaves(l,1,k,j,i) + fwaves(l,0,k+1,j,i) ) / dom.dz;
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
  //Exchange halos in the x-direction
  exch.haloInit      ();
  exch.haloPackN_x   (dom, state, numState);
  exch.haloExchange_x(dom, par);
  exch.haloUnpackN_x (dom, state, numState);

  //Exchange halos in the y-direction
  exch.haloInit      ();
  exch.haloPackN_x   (dom, state, numState);
  exch.haloExchange_x(dom, par);
  exch.haloUnpackN_x (dom, state, numState);
  
  // Boundaries for the fluid state in the z-direction
  stateBoundariesZ(state, dom);
}


void TendenciesThetaPrimADER::stateBoundariesZ(realArr &state, Domain const &dom) {
  // for (int j=0; j<dom.ny; j++) {
  //   for (int i=0; i<dom.nx; i++) {
  //     for (int ii=0; ii<hs; ii++) {
  yakl::parallel_for( dom.ny*dom.nx*hs , YAKL_LAMBDA (int const iGlob) {
    int j, i, ii;
    yakl::unpackIndices(iGlob,dom.ny,dom.nx,hs,j,i,ii);
    state(idR ,ii,hs+j,hs+i) = state(idR ,hs,hs+j,hs+i);
    state(idU ,ii,hs+j,hs+i) = state(idU ,hs,hs+j,hs+i);
    state(idV ,ii,hs+j,hs+i) = state(idV ,hs,hs+j,hs+i);
    state(idW ,ii,hs+j,hs+i) = 0;
    state(idRT,ii,hs+j,hs+i) = state(idRT,hs,hs+j,hs+i);
    pressure(  ii,hs+j,hs+i) = pressure(  hs,hs+j,hs+i);

    state(idR ,dom.nz+hs+ii,hs+j,hs+i) = state(idR ,dom.nz+hs-1,hs+j,hs+i);
    state(idU ,dom.nz+hs+ii,hs+j,hs+i) = state(idU ,dom.nz+hs-1,hs+j,hs+i);
    state(idV ,dom.nz+hs+ii,hs+j,hs+i) = state(idV ,dom.nz+hs-1,hs+j,hs+i);
    state(idW ,dom.nz+hs+ii,hs+j,hs+i) = 0;
    state(idRT,dom.nz+hs+ii,hs+j,hs+i) = state(idRT,dom.nz+hs-1,hs+j,hs+i);
    pressure(  dom.nz+hs+ii,hs+j,hs+i) = pressure(  dom.nz+hs-1,hs+j,hs+i);
  });
}


void TendenciesThetaPrimADER::edgeBoundariesZ(realArr &stateLimits, Domain const &dom) {
  // for (int j=0; j<dom.ny; j++) {
  //   for (int i=0; i<dom.nx; i++) {
  yakl::parallel_for( dom.ny*dom.nx , YAKL_LAMBDA (int const iGlob) {
    int j, i;
    yakl::unpackIndices(iGlob,dom.ny,dom.nx,j,i);
    stateLimits(idR ,0,0     ,j,i) = stateLimits(idR ,1,0     ,j,i);
    stateLimits(idU ,0,0     ,j,i) = stateLimits(idU ,1,0     ,j,i);
    stateLimits(idV ,0,0     ,j,i) = stateLimits(idV ,1,0     ,j,i);
    stateLimits(idW ,0,0     ,j,i) = 0;
    stateLimits(idW ,1,0     ,j,i) = 0;
    stateLimits(idRT,0,0     ,j,i) = stateLimits(idRT,1,0     ,j,i);
    stateLimits(idP ,0,0     ,j,i) = stateLimits(idP ,1,0     ,j,i);

    stateLimits(idR ,1,dom.nz,j,i) = stateLimits(idR ,0,dom.nz,j,i);
    stateLimits(idU ,1,dom.nz,j,i) = stateLimits(idU ,0,dom.nz,j,i);
    stateLimits(idV ,1,dom.nz,j,i) = stateLimits(idV ,0,dom.nz,j,i);
    stateLimits(idW ,0,dom.nz,j,i) = 0;
    stateLimits(idW ,1,dom.nz,j,i) = 0;
    stateLimits(idRT,1,dom.nz,j,i) = stateLimits(idRT,0,dom.nz,j,i);
    stateLimits(idP ,1,dom.nz,j,i) = stateLimits(idP ,0,dom.nz,j,i);
  });
}



