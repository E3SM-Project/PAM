
#include "TendenciesThetaConsADER.h"


void TendenciesThetaConsADER::initialize(Domain const &dom) {
  TransformMatrices<real> trans;

  fluxLimits  = realArr("fluxLimits",numState,2,dom.nz+1,dom.ny+1,dom.nx+1);
  stateLimits = realArr("srcLimits" ,numState,2,dom.nz+1,dom.ny+1,dom.nx+1);
  flux        = realArr("flux"      ,numState  ,dom.nz+1,dom.ny+1,dom.nx+1);
  src         = realArr("src"       ,dom.nz,dom.ny,dom.nx);
  stateGLL    = realArr("stateGLL"  ,numState,dom.nz,dom.ny,dom.nx,tord);

  // Setup the matrix to transform a stencil of ord cell averages into tord GLL points
  if (dom.doWeno) {
    trans.coefs_to_gll_lower( to_gll );
  } else {
    trans.sten_to_gll_lower( to_gll );
  }

  trans.weno_sten_to_coefs(wenoRecon);

  SArray<real,tord,tord> g2c, c2d, c2g;
  trans.gll_to_coefs  (g2c);
  trans.coefs_to_deriv(c2d);
  trans.coefs_to_gll  (c2g);
  aderDerivX = (c2g * c2d * g2c) / dom.dx;
  aderDerivY = (c2g * c2d * g2c) / dom.dy;
  aderDerivZ = (c2g * c2d * g2c) / dom.dz;

  trans.get_gll_weights(gllWts);

  wenoSetIdealSigma(wenoIdl,wenoSigma);
}


void TendenciesThetaConsADER::compEulerTend_X(realArr &state, Domain const &dom, Exchange &exch, Parallel const &par, realArr &tend) {
  auto &wenoRecon   = this->wenoRecon  ;
  auto &to_gll      = this->to_gll     ;
  auto &stateLimits = this->stateLimits;
  auto &fluxLimits  = this->fluxLimits ;
  auto &wenoIdl     = this->wenoIdl    ;
  auto &wenoSigma   = this->wenoSigma  ;
  auto &aderDerivX  = this->aderDerivX ;
  auto &stateGLL    = this->stateGLL   ;
  auto &flux        = this->flux       ;

  //Exchange halos in the x-direction
  exch.haloInit      ();
  exch.haloPackN_x   (dom, state, numState);
  exch.haloExchange_x(dom, par);
  exch.haloUnpackN_x (dom, state, numState);

  // Reconstruct to tord GLL points in the x-direction
  // for (int k=0; k<dom.nz; k++) {
  //   for (int j=0; j<dom.ny; j++) {
  //     for (int i=0; i<dom.nx; i++) {
  yakl::parallel_for( dom.nz,dom.ny,dom.nx , YAKL_LAMBDA (int k, int j, int i) {
    // Compute tord GLL points of the state vector
    for (int l=0; l<numState; l++) {
      SArray<real,ord> stencil;
      SArray<real,tord> gllPts;
      for (int ii=0; ii<ord; ii++) { stencil(ii) = state(l,hs+k,hs+j,i+ii); }
      reconStencil(stencil, gllPts, dom.doWeno, wenoRecon, to_gll, wenoIdl, wenoSigma);
      for (int ii=0; ii<tord; ii++) { stateGLL(l,k,j,i,ii) = gllPts(ii); }
    }
  });

  // for (int k=0; k<dom.nz; k++) {
  //   for (int j=0; j<dom.ny; j++) {
  //     for (int i=0; i<dom.nx; i++) {
  yakl::parallel_for( dom.nz,dom.ny,dom.nx , YAKL_LAMBDA (int k, int j, int i) {
    SArray<real,numState,tord,tord> stateDTs;  // GLL state values
    SArray<real,numState,tord,tord> fluxDTs;   // GLL flux values

    for (int l=0; l<numState; l++) {
      for (int ii=0; ii<tord; ii++) { stateDTs(l,0,ii) = stateGLL(l,k,j,i,ii); }
    }

    for (int ii=0; ii<tord; ii++) {
      stateDTs(idR ,0,ii) += dom.hyDensCells     (hs+k);
      stateDTs(idRT,0,ii) += dom.hyDensThetaCells(hs+k);
    }

    // Compute DTs of the state and flux, and collapse down into a time average
    diffTransformEulerConsX( stateDTs , fluxDTs , aderDerivX );
    timeAvg( stateDTs , dom );
    timeAvg( fluxDTs  , dom );

    // Store state and flux limits into a globally indexed array
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

  // Riemann solver
  // for (int k=0; k<dom.nz; k++) {
  //   for (int j=0; j<dom.ny; j++) {
  //     for (int i=0; i<dom.nx+1; i++) {
  yakl::parallel_for( dom.nz,dom.ny,(dom.nx+1) , YAKL_LAMBDA (int k, int j, int i) {
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
      // flux(l,k,j,i) = 0.5_fp * ( f2(l) + f1(l) - dom.cfl*dom.dx/dom.dt * (s2(l) - s1(l)) );
    }
  });

  // Form the tendencies
  // for (int l=0; l<numState; l++) {
  //   for (int k=0; k<dom.nz; k++) {
  //     for (int j=0; j<dom.ny; j++) {
  //       for (int i=0; i<dom.nx; i++) {
  yakl::parallel_for( numState,dom.nz,dom.ny,dom.nx , YAKL_LAMBDA (int l, int k, int j , int i) {
    tend(l,k,j,i) = - ( flux(l,k,j,i+1) - flux(l,k,j,i) ) / dom.dx;
  });
}


void TendenciesThetaConsADER::compEulerTend_Y(realArr &state, Domain const &dom, Exchange &exch, Parallel const &par, realArr &tend) {
  auto &wenoRecon   = this->wenoRecon  ;
  auto &to_gll      = this->to_gll     ;
  auto &stateLimits = this->stateLimits;
  auto &fluxLimits  = this->fluxLimits ;
  auto &wenoIdl     = this->wenoIdl    ;
  auto &wenoSigma   = this->wenoSigma  ;
  auto &aderDerivY  = this->aderDerivY ;
  auto &stateGLL    = this->stateGLL   ;
  auto &flux        = this->flux       ;

  //Exchange halos in the y-direction
  exch.haloInit      ();
  exch.haloPackN_y   (dom, state, numState);
  exch.haloExchange_y(dom, par);
  exch.haloUnpackN_y (dom, state, numState);

  // Reconstruct to tord GLL points in the y-direction
  // for (int k=0; k<dom.nz; k++) {
  //   for (int j=0; j<dom.ny; j++) {
  //     for (int i=0; i<dom.nx; i++) {
  yakl::parallel_for( dom.nz,dom.ny,dom.nx , YAKL_LAMBDA (int k, int j, int i) {
    // Compute GLL points from cell averages
    for (int l=0; l<numState; l++) {
      SArray<real,ord> stencil;
      SArray<real,tord> gllPts;
      for (int ii=0; ii<ord; ii++) { stencil(ii) = state(l,hs+k,j+ii,hs+i); }
      reconStencil(stencil, gllPts, dom.doWeno, wenoRecon, to_gll, wenoIdl, wenoSigma);
      for (int ii=0; ii<tord; ii++) { stateGLL(l,k,j,i,ii) = gllPts(ii); }
    }
  });

  // for (int k=0; k<dom.nz; k++) {
  //   for (int j=0; j<dom.ny; j++) {
  //     for (int i=0; i<dom.nx; i++) {
  yakl::parallel_for( dom.nz,dom.ny,dom.nx , YAKL_LAMBDA (int k, int j, int i) {
    SArray<real,numState,tord,tord> stateDTs;  // GLL state values
    SArray<real,numState,tord,tord> fluxDTs;   // GLL flux values

    // Compute GLL points from cell averages
    for (int l=0; l<numState; l++) {
      for (int ii=0; ii<tord; ii++) { stateDTs(l,0,ii) = stateGLL(l,k,j,i,ii); }
    }
    for (int ii=0; ii<tord; ii++) {
      stateDTs(idR ,0,ii) += dom.hyDensCells     (hs+k);
      stateDTs(idRT,0,ii) += dom.hyDensThetaCells(hs+k);
    }

    // Compute DTs of the state and flux, and collapse down into a time average
    diffTransformEulerConsY( stateDTs , fluxDTs , aderDerivY );
    timeAvg( stateDTs , dom );
    timeAvg( fluxDTs  , dom );

    // Store state and flux limits into a globally indexed array
    for (int l=0; l<numState; l++) {
      // Store the left cell edge state and flux estimates
      stateLimits(l,1,k,j  ,i) = stateDTs(l,0,0);
      fluxLimits (l,1,k,j  ,i) = fluxDTs (l,0,0);

      // Store the Right cell edge state and flux estimates
      stateLimits(l,0,k,j+1,i) = stateDTs(l,0,tord-1);
      fluxLimits (l,0,k,j+1,i) = fluxDTs (l,0,tord-1);
    }
  });

  //Reconcile the edge fluxes via MPI exchange.
  exch.haloInit      ();
  exch.edgePackN_y   (dom, stateLimits, numState);
  exch.edgePackN_y   (dom, fluxLimits , numState);
  exch.edgeExchange_y(dom, par);
  exch.edgeUnpackN_y (dom, stateLimits, numState);
  exch.edgeUnpackN_y (dom, fluxLimits , numState);

  // Riemann solver
  // for (int k=0; k<dom.nz; k++) {
  //   for (int j=0; j<dom.ny+1; j++) {
  //     for (int i=0; i<dom.nx; i++) {
  yakl::parallel_for( dom.nz,(dom.ny+1),dom.nx , YAKL_LAMBDA (int k, int j, int i) {
    SArray<real,numState> s1, s2, f1, f2, upw;
    for (int l=0; l<numState; l++) {
      s1(l) = stateLimits(l,0,k,j,i);
      s2(l) = stateLimits(l,1,k,j,i);
      f1(l) = fluxLimits (l,0,k,j,i);
      f2(l) = fluxLimits (l,1,k,j,i);
    }
    riemannY(s1, s2, f1, f2, upw);
    for (int l=0; l<numState; l++) {
      flux(l,k,j,i) = upw(l);
      // flux(l,k,j,i) = 0.5_fp * ( f2(l) + f1(l) - dom.cfl*dom.dy/dom.dt * (s2(l) - s1(l)) );
    }
  });

  // Form the tendencies
  // for (int l=0; l<numState; l++) {
  //   for (int k=0; k<dom.nz; k++) {
  //     for (int j=0; j<dom.ny; j++) {
  //       for (int i=0; i<dom.nx; i++) {
  yakl::parallel_for( numState,dom.nz,dom.ny,dom.nx , YAKL_LAMBDA (int l, int k, int j, int i) {
    tend(l,k,j,i) = - ( flux(l,k,j+1,i) - flux(l,k,j,i) ) / dom.dy;
  });
}


void TendenciesThetaConsADER::compEulerTend_Z(realArr &state, Domain const &dom, Exchange &exch, Parallel const &par, realArr &tend) {
  auto &wenoRecon   = this->wenoRecon  ;
  auto &to_gll      = this->to_gll     ;
  auto &stateLimits = this->stateLimits;
  auto &fluxLimits  = this->fluxLimits ;
  auto &wenoIdl     = this->wenoIdl    ;
  auto &wenoSigma   = this->wenoSigma  ;
  auto &aderDerivZ  = this->aderDerivZ ;
  auto &stateGLL    = this->stateGLL   ;
  auto &src         = this->src        ;
  auto &gllWts      = this->gllWts     ;
  auto &flux        = this->flux       ;

  // Boundaries for the fluid state in the z-direction
  stateBoundariesZ(state, dom);

  // Reconstruct tord GLL points in the z-direction
  // for (int k=0; k<dom.nz; k++) {
  //   for (int j=0; j<dom.ny; j++) {
  //     for (int i=0; i<dom.nx; i++) {
  yakl::parallel_for( dom.nz,dom.ny,dom.nx , YAKL_LAMBDA (int k, int j, int i) {
    // Compute GLL points from cell averages
    for (int l=0; l<numState; l++) {
      SArray<real,ord> stencil;
      SArray<real,tord> gllPts;
      for (int ii=0; ii<ord; ii++) { stencil(ii) = state(l,k+ii,hs+j,hs+i); }
      reconStencil(stencil, gllPts, dom.doWeno, wenoRecon, to_gll, wenoIdl, wenoSigma);
      for (int ii=0; ii<tord; ii++) { stateGLL(l,k,j,i,ii) = gllPts(ii); }
    }
  });

  // for (int k=0; k<dom.nz; k++) {
  //   for (int j=0; j<dom.ny; j++) {
  //     for (int i=0; i<dom.nx; i++) {
  yakl::parallel_for( dom.nz,dom.ny,dom.nx , YAKL_LAMBDA (int k, int j, int i) {
    SArray<real,numState,tord,tord> stateDTs;  // GLL state values
    SArray<real,numState,tord,tord> fluxDTs;   // GLL flux values
    SArray<real,tord,tord> sourceDTs;   // GLL source values
    SArray<real,tord> hyRHOT;
    SArray<real,tord> hyRHO;

    // Compute GLL points from cell averages
    for (int l=0; l<numState; l++) {
      for (int ii=0; ii<tord; ii++) { stateDTs(l,0,ii) = stateGLL(l,k,j,i,ii); }
    }
    for (int ii=0; ii<tord; ii++) {
      stateDTs(idR ,0,ii) += dom.hyDensGLL     (k,ii);
      stateDTs(idRT,0,ii) += dom.hyDensThetaGLL(k,ii);
      hyRHOT(ii) = dom.hyDensThetaGLL(k,ii);
      hyRHO (ii) = dom.hyDensGLL     (k,ii);
    }

    // Boundary conditions
    if (k == 0       ) { stateDTs(idRW,0,0     ) = 0; }
    if (k == dom.nz-1) { stateDTs(idRW,0,tord-1) = 0; }

    // Compute DTs of the state and flux, and collapse down into a time average
    diffTransformEulerConsZ( stateDTs , fluxDTs , sourceDTs , aderDerivZ , hyRHOT, hyRHO );
    timeAvg( stateDTs  , dom );
    timeAvg( fluxDTs   , dom );
    timeAvg( sourceDTs , dom );

    // Boundary conditions
    if (k == 0       ) { stateDTs(idRW,0,0     ) = 0; }
    if (k == dom.nz-1) { stateDTs(idRW,0,tord-1) = 0; }

    // Store state and flux limits into a globally indexed array
    for (int l=0; l<numState; l++) {
      // Store the left cell edge state and flux estimates
      stateLimits(l,1,k  ,j,i) = stateDTs(l,0,0);
      fluxLimits (l,1,k  ,j,i) = fluxDTs (l,0,0);

      // Store the Right cell edge state and flux estimates
      stateLimits(l,0,k+1,j,i) = stateDTs(l,0,tord-1);
      fluxLimits (l,0,k+1,j,i) = fluxDTs (l,0,tord-1);
    }
    src(k,j,i) = 0;
    for (int ii=0; ii<tord; ii++) {
      src(k,j,i) += sourceDTs(0,ii) * gllWts(ii);
    }
  });

  // Apply boundary conditions to fluxes and state values
  edgeBoundariesZ(stateLimits, fluxLimits, dom);

  // Riemann solver
  // for (int k=0; k<dom.nz+1; k++) {
  //   for (int j=0; j<dom.ny; j++) {
  //     for (int i=0; i<dom.nx; i++) {
  yakl::parallel_for( (dom.nz+1),dom.ny,dom.nx , YAKL_LAMBDA (int k, int j, int i) {
    SArray<real,numState> s1, s2, f1, f2, upw;
    for (int l=0; l<numState; l++) {
      s1(l) = stateLimits(l,0,k,j,i);
      s2(l) = stateLimits(l,1,k,j,i);
      f1(l) = fluxLimits (l,0,k,j,i);
      f2(l) = fluxLimits (l,1,k,j,i);
    }
    riemannZ(s1, s2, f1, f2, upw);
    for (int l=0; l<numState; l++) {
      flux(l,k,j,i) = upw(l);
      // flux(l,k,j,i) = 0.5_fp * ( f2(l) + f1(l) - dom.cfl*dom.dz/dom.dt * (s2(l) - s1(l)) );
    }
  });

  // Form the tendencies
  // for (int l=0; l<numState; l++) {
  //   for (int k=0; k<dom.nz; k++) {
  //     for (int j=0; j<dom.ny; j++) {
  //       for (int i=0; i<dom.nx; i++) {
  yakl::parallel_for( numState,dom.nz,dom.ny,dom.nx , YAKL_LAMBDA (int l, int k, int j, int i) {
    tend(l,k,j,i) = - ( flux(l,k+1,j,i) - flux(l,k,j,i) ) / dom.dz;
    if (l==idRW) {
      tend(l,k,j,i) += src(k,j,i);
    }
  });
}


void TendenciesThetaConsADER::compEulerTend_S(realArr &state, Domain const &dom, Exchange &exch, Parallel const &par, realArr &tend) {
  // for (int l=0; l<numState; l++) {
  //   for (int k=0; k<dom.nz; k++) {
  //     for (int j=0; j<dom.ny; j++) {
  //       for (int i=0; i<dom.nx; i++) {
  yakl::parallel_for( numState,dom.nz,dom.ny,dom.nx , YAKL_LAMBDA (int l, int k, int j, int i) {
    tend(l,k,j,i) = 0;
  });
}


void TendenciesThetaConsADER::stateBoundariesZ(realArr &state, Domain const &dom) {
  // for (int j=0; j<dom.ny; j++) {
  //   for (int i=0; i<dom.nx; i++) {
  //     for (int ii=0; ii<hs; ii++) {
  yakl::parallel_for( dom.ny,dom.nx,hs , YAKL_LAMBDA (int j, int i, int ii) {
    state(idR ,ii,hs+j,hs+i) = state(idR ,hs,hs+j,hs+i);
    state(idRU,ii,hs+j,hs+i) = state(idRU,hs,hs+j,hs+i);
    state(idRV,ii,hs+j,hs+i) = state(idRV,hs,hs+j,hs+i);
    state(idRW,ii,hs+j,hs+i) = 0;
    state(idRT,ii,hs+j,hs+i) = state(idRT,hs,hs+j,hs+i);

    state(idR ,dom.nz+hs+ii,hs+j,hs+i) = state(idR ,dom.nz+hs-1,hs+j,hs+i);
    state(idRU,dom.nz+hs+ii,hs+j,hs+i) = state(idRU,dom.nz+hs-1,hs+j,hs+i);
    state(idRV,dom.nz+hs+ii,hs+j,hs+i) = state(idRV,dom.nz+hs-1,hs+j,hs+i);
    state(idRW,dom.nz+hs+ii,hs+j,hs+i) = 0;
    state(idRT,dom.nz+hs+ii,hs+j,hs+i) = state(idRT,dom.nz+hs-1,hs+j,hs+i);
  });
}


void TendenciesThetaConsADER::edgeBoundariesZ(realArr &stateLimits, realArr &fluxLimits, Domain const &dom) {
  // for (int j=0; j<dom.ny; j++) {
  //   for (int i=0; i<dom.nx; i++) {
  yakl::parallel_for( dom.ny,dom.nx , YAKL_LAMBDA (int j, int i) {
    stateLimits(idR ,0,0     ,j,i) = stateLimits(idR ,1,0     ,j,i);
    stateLimits(idRU,0,0     ,j,i) = stateLimits(idRU,1,0     ,j,i);
    stateLimits(idRV,0,0     ,j,i) = stateLimits(idRV,1,0     ,j,i);
    stateLimits(idRW,0,0     ,j,i) = 0;
    stateLimits(idRW,1,0     ,j,i) = 0;
    stateLimits(idRT,0,0     ,j,i) = stateLimits(idRT,1,0     ,j,i);

    stateLimits(idR ,1,dom.nz,j,i) = stateLimits(idR ,0,dom.nz,j,i);
    stateLimits(idRU,1,dom.nz,j,i) = stateLimits(idRU,0,dom.nz,j,i);
    stateLimits(idRV,1,dom.nz,j,i) = stateLimits(idRV,0,dom.nz,j,i);
    stateLimits(idRW,0,dom.nz,j,i) = 0;
    stateLimits(idRW,1,dom.nz,j,i) = 0;
    stateLimits(idRT,1,dom.nz,j,i) = stateLimits(idRT,0,dom.nz,j,i);

    fluxLimits(idR ,0,0     ,j,i) = fluxLimits(idR ,1,0     ,j,i);
    fluxLimits(idRU,0,0     ,j,i) = fluxLimits(idRU,1,0     ,j,i);
    fluxLimits(idRV,0,0     ,j,i) = fluxLimits(idRV,1,0     ,j,i);
    fluxLimits(idRW,0,0     ,j,i) = fluxLimits(idRW,1,0     ,j,i);
    fluxLimits(idRT,0,0     ,j,i) = fluxLimits(idRT,1,0     ,j,i);

    fluxLimits(idR ,1,dom.nz,j,i) = fluxLimits(idR ,0,dom.nz,j,i);
    fluxLimits(idRU,1,dom.nz,j,i) = fluxLimits(idRU,0,dom.nz,j,i);
    fluxLimits(idRV,1,dom.nz,j,i) = fluxLimits(idRV,0,dom.nz,j,i);
    fluxLimits(idRW,1,dom.nz,j,i) = fluxLimits(idRW,0,dom.nz,j,i);
    fluxLimits(idRT,1,dom.nz,j,i) = fluxLimits(idRT,0,dom.nz,j,i);
  });
}


void TendenciesThetaConsADER::compStrakaTend(realArr &state, Domain const &dom, Exchange &exch, Parallel const &par, realArr &tend) {
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

  // for (int k=0; k<dom.nz; k++) {
  //   for (int j=0; j<dom.ny; j++) {
  //     for (int i=0; i<dom.nx; i++) {
  yakl::parallel_for( dom.nz,dom.ny,dom.nx , YAKL_LAMBDA (int k, int j, int i) {
    real r = ( state(idR,hs+k,hs+j,hs+i) + dom.hyDensCells(hs+k) );
    SArray<real,numState,3> sten;

    tend(idR,k,j,i) = 0.;

    real dx2 = dom.dx*dom.dx;
    for (int l=1; l<numState; l++) {
      for (int ii=-1; ii<=1; ii++) {
        if (l == idRT) {
          sten(l,ii+1) = ( state(l,hs+k,hs+j,hs+i+ii) + dom.hyDensThetaCells(hs+k) ) / ( state(idR,hs+k,hs+j,hs+i+ii) + dom.hyDensCells(hs+k) );
        } else {
          sten(l,ii+1) = state(l,hs+k,hs+j,hs+i+ii) / ( state(idR,hs+k,hs+j,hs+i+ii) + dom.hyDensCells(hs+k) );
        }
      }
      tend(l,k,j,i)  = 75 * r * ( sten(l,2) - 2*sten(l,1) + sten(l,0) ) / dx2;
    }

    real dy2 = dom.dy*dom.dy;
    for (int l=1; l<numState; l++) {
      for (int ii=-1; ii<=1; ii++) {
        if (l == idRT) {
          sten(l,ii+1) = ( state(l,hs+k,hs+j+ii,hs+i) + dom.hyDensThetaCells(hs+k) ) / ( state(idR,hs+k,hs+j+ii,hs+i) + dom.hyDensCells(hs+k) );
        } else {
          sten(l,ii+1) = state(l,hs+k,hs+j+ii,hs+i) / ( state(idR,hs+k,hs+j+ii,hs+i) + dom.hyDensCells(hs+k) );
        }
      }
      tend(l,k,j,i) += 75 * r * ( sten(l,2) - 2*sten(l,1) + sten(l,0) ) / dy2;
    }

    real dz2 = dom.dz*dom.dz;
    for (int l=1; l<numState; l++) {
      for (int ii=-1; ii<=1; ii++) {
        if (l == idRT) {
          sten(l,ii+1) = ( state(l,hs+k+ii,hs+j,hs+i) + dom.hyDensThetaCells(hs+k+ii) ) / ( state(idR,hs+k+ii,hs+j,hs+i) + dom.hyDensCells(hs+k+ii) );
        } else {
          sten(l,ii+1) = state(l,hs+k+ii,hs+j,hs+i) / ( state(idR,hs+k+ii,hs+j,hs+i) + dom.hyDensCells(hs+k+ii) );
        }
      }
      tend(l,k,j,i) += 75 * r * ( sten(l,2) - 2*sten(l,1) + sten(l,0) ) / dz2;
    }
  });
}



