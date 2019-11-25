
#include "TendenciesThetaConsSD.h"



void TendenciesThetaConsSD::initialize(Domain const &dom) {
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

  trans.get_gll_weights(gllWts);

  wenoSetIdealSigma(wenoIdl,wenoSigma);

}


void TendenciesThetaConsSD::compEulerTend_X(realArr &state, Domain const &dom, Exchange &exch, Parallel const &par, realArr &tend) {
  auto &stateLimits = this->stateLimits;
  auto &fluxLimits  = this->fluxLimits ;
  auto &flux        = this->flux       ;
  auto &src         = this->src        ;
  auto &stateGLL    = this->stateGLL   ;
  auto &gllWts      = this->gllWts     ;
  auto &to_gll      = this->to_gll     ;
  auto &wenoRecon   = this->wenoRecon  ;
  auto &wenoIdl     = this->wenoIdl    ;
  auto &wenoSigma   = this->wenoSigma  ;

  //Exchange halos in the x-direction
  exch.haloInit      ();
  exch.haloPackN_x   (dom, state, numState);
  exch.haloExchange_x(dom, par);
  exch.haloUnpackN_x (dom, state, numState);

  // Reconstruct to tord GLL points in the x-direction
  // for (int k=0; k<dom.nz; k++) {
  //   for (int j=0; j<dom.ny; j++) {
  //     for (int i=0; i<dom.nx; i++) {
  yakl::parallel_for( dom.nz*dom.ny*dom.nx , YAKL_LAMBDA (int const iGlob) {
    int k, j, i;
    yakl::unpackIndices(iGlob,dom.nz,dom.ny,dom.nx,k,j,i);
    SArray<real,numState,tord> gllState;  // GLL state values
    SArray<real,numState,tord> gllFlux;   // GLL flux values

    // Compute tord GLL points of the state vector
    for (int l=0; l<numState; l++) {
      SArray<real,ord> stencil;
      SArray<real,tord> gllPts;
      for (int ii=0; ii<ord; ii++) { stencil(ii) = state(l,hs+k,hs+j,i+ii); }
      reconStencil(stencil, gllPts, dom.doWeno, wenoRecon, to_gll, wenoIdl, wenoSigma);
      for (int ii=0; ii<tord; ii++) { gllState(l,ii) = gllPts(ii); }
    }
    for (int ii=0; ii<tord; ii++) {
      gllState(idR ,ii) += dom.hyDensCells     (hs+k);
      gllState(idRT,ii) += dom.hyDensThetaCells(hs+k);
    }

    // Compute fluxes and at the GLL points
    for (int ii=0; ii<tord; ii++) {
      real r = gllState(idR ,ii);
      real u = gllState(idRU,ii) / r;
      real v = gllState(idRV,ii) / r;
      real w = gllState(idRW,ii) / r;
      real t = gllState(idRT,ii) / r;
      real p = C0 * pow( r*t , GAMMA );

      gllFlux(idR ,ii) = r*u;
      gllFlux(idRU,ii) = r*u*u + p;
      gllFlux(idRV,ii) = r*u*v;
      gllFlux(idRW,ii) = r*u*w;
      gllFlux(idRT,ii) = r*u*t;
    }

    // Store state and flux limits into a globally indexed array
    for (int l=0; l<numState; l++) {
      // Store the left cell edge state and flux estimates
      stateLimits(l,1,k,j,i  ) = gllState(l,0);
      fluxLimits (l,1,k,j,i  ) = gllFlux (l,0);

      // Store the Right cell edge state and flux estimates
      stateLimits(l,0,k,j,i+1) = gllState(l,tord-1);
      fluxLimits (l,0,k,j,i+1) = gllFlux (l,tord-1);
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
  yakl::parallel_for( dom.nz*dom.ny*(dom.nx+1) , YAKL_LAMBDA (int const iGlob) {
    int k, j, i;
    yakl::unpackIndices(iGlob,dom.nz,dom.ny,dom.nx+1,k,j,i);
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
  yakl::parallel_for( numState*dom.nz*dom.ny*dom.nx , YAKL_LAMBDA (int const iGlob) {
    int l, k, j, i;
    yakl::unpackIndices(iGlob,numState,dom.nz,dom.ny,dom.nx,l,k,j,i);
    tend(l,k,j,i) = - ( flux(l,k,j,i+1) - flux(l,k,j,i) ) / dom.dx;
  });
}


void TendenciesThetaConsSD::compEulerTend_Y(realArr &state, Domain const &dom, Exchange &exch, Parallel const &par, realArr &tend) {
  auto &stateLimits = this->stateLimits;
  auto &fluxLimits  = this->fluxLimits ;
  auto &flux        = this->flux       ;
  auto &src         = this->src        ;
  auto &stateGLL    = this->stateGLL   ;
  auto &gllWts      = this->gllWts     ;
  auto &to_gll      = this->to_gll     ;
  auto &wenoRecon   = this->wenoRecon  ;
  auto &wenoIdl     = this->wenoIdl    ;
  auto &wenoSigma   = this->wenoSigma  ;

  //Exchange halos in the y-direction
  exch.haloInit      ();
  exch.haloPackN_y   (dom, state, numState);
  exch.haloExchange_y(dom, par);
  exch.haloUnpackN_y (dom, state, numState);

  // Reconstruct to tord GLL points in the y-direction
  // for (int k=0; k<dom.nz; k++) {
  //   for (int j=0; j<dom.ny; j++) {
  //     for (int i=0; i<dom.nx; i++) {
  yakl::parallel_for( dom.nz*dom.ny*dom.nx , YAKL_LAMBDA (int const iGlob) {
    int k, j, i;
    yakl::unpackIndices(iGlob,dom.nz,dom.ny,dom.nx,k,j,i);
    SArray<real,numState,tord> gllState;
    SArray<real,numState,tord> gllFlux;

    // Compute GLL points from cell averages
    for (int l=0; l<numState; l++) {
      SArray<real,ord> stencil;
      SArray<real,tord> gllPts;
      for (int ii=0; ii<ord; ii++) { stencil(ii) = state(l,hs+k,j+ii,hs+i); }
      reconStencil(stencil, gllPts, dom.doWeno, wenoRecon, to_gll, wenoIdl, wenoSigma);
      for (int ii=0; ii<tord; ii++) { gllState(l,ii) = gllPts(ii); }
    }
    for (int ii=0; ii<tord; ii++) {
      gllState(idR ,ii) += dom.hyDensCells     (hs+k);
      gllState(idRT,ii) += dom.hyDensThetaCells(hs+k);
    }

    // Compute fluxes and at the GLL points
    for (int ii=0; ii<tord; ii++) {
      real r = gllState(idR ,ii);
      real u = gllState(idRU,ii) / r;
      real v = gllState(idRV,ii) / r;
      real w = gllState(idRW,ii) / r;
      real t = gllState(idRT,ii) / r;
      real p = C0 * pow( r*t , GAMMA );

      gllFlux(idR ,ii) = r*v;
      gllFlux(idRU,ii) = r*v*u;
      gllFlux(idRV,ii) = r*v*v + p;
      gllFlux(idRW,ii) = r*v*w;
      gllFlux(idRT,ii) = r*v*t;
    }

    // Store state and flux limits into a globally indexed array
    for (int l=0; l<numState; l++) {
      // Store the left cell edge state and flux estimates
      stateLimits(l,1,k,j  ,i) = gllState(l,0);
      fluxLimits (l,1,k,j  ,i) = gllFlux (l,0);

      // Store the Right cell edge state and flux estimates
      stateLimits(l,0,k,j+1,i) = gllState(l,tord-1);
      fluxLimits (l,0,k,j+1,i) = gllFlux (l,tord-1);
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
  yakl::parallel_for( dom.nz*(dom.ny+1)*dom.nx , YAKL_LAMBDA (int const iGlob) {
    int k, j, i;
    yakl::unpackIndices(iGlob,dom.nz,dom.ny+1,dom.nx,k,j,i);
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
  yakl::parallel_for( numState*dom.nz*dom.ny*dom.nx , YAKL_LAMBDA (int const iGlob) {
    int l, k, j, i;
    yakl::unpackIndices(iGlob,numState,dom.nz,dom.ny,dom.nx,l,k,j,i);
    tend(l,k,j,i) = - ( flux(l,k,j+1,i) - flux(l,k,j,i) ) / dom.dy;
  });
}


void TendenciesThetaConsSD::compEulerTend_Z(realArr &state, Domain const &dom, Exchange &exch, Parallel const &par, realArr &tend) {
  auto &stateLimits = this->stateLimits;
  auto &fluxLimits  = this->fluxLimits ;
  auto &flux        = this->flux       ;
  auto &src         = this->src        ;
  auto &stateGLL    = this->stateGLL   ;
  auto &gllWts      = this->gllWts     ;
  auto &to_gll      = this->to_gll     ;
  auto &wenoRecon   = this->wenoRecon  ;
  auto &wenoIdl     = this->wenoIdl    ;
  auto &wenoSigma   = this->wenoSigma  ;

  // Boundaries for the fluid state in the z-direction
  stateBoundariesZ(state, dom);

  // Reconstruct to tord GLL points in the x-direction
  // for (int k=0; k<dom.nz; k++) {
  //   for (int j=0; j<dom.ny; j++) {
  //     for (int i=0; i<dom.nx; i++) {
  yakl::parallel_for( dom.nz*dom.ny*dom.nx , YAKL_LAMBDA (int const iGlob) {
    int k, j, i;
    yakl::unpackIndices(iGlob,dom.nz,dom.ny,dom.nx,k,j,i);
    SArray<real,numState,tord> gllState;
    SArray<real,numState,tord> gllFlux;

    // Compute GLL points from cell averages
    for (int l=0; l<numState; l++) {
      SArray<real,ord> stencil;
      SArray<real,tord> gllPts;
      for (int ii=0; ii<ord; ii++) { stencil(ii) = state(l,k+ii,hs+j,hs+i); }
      reconStencil(stencil, gllPts, dom.doWeno, wenoRecon, to_gll, wenoIdl, wenoSigma);
      for (int ii=0; ii<tord; ii++) { gllState(l,ii) = gllPts(ii); }
    }
    for (int ii=0; ii<tord; ii++) {
      gllState(idR ,ii) += dom.hyDensGLL     (k,ii);
      gllState(idRT,ii) += dom.hyDensThetaGLL(k,ii);
    }

    // Boundary conditions
    if (k == 0       ) { gllState(idRW,0     ) = 0; }
    if (k == dom.nz-1) { gllState(idRW,tord-1) = 0; }

    // Compute fluxes and at the GLL points
    for (int ii=0; ii<tord; ii++) {
      real r = gllState(idR ,ii);
      real u = gllState(idRU,ii) / r;
      real v = gllState(idRV,ii) / r;
      real w = gllState(idRW,ii) / r;
      real t = gllState(idRT,ii) / r;
      real p = C0 * pow( r*t , GAMMA );

      gllFlux(idR ,ii) = r*w;
      gllFlux(idRU,ii) = r*w*u;
      gllFlux(idRV,ii) = r*w*v;
      gllFlux(idRW,ii) = r*w*w + p - C0*pow(dom.hyDensThetaGLL(k,ii),GAMMA);
      gllFlux(idRT,ii) = r*w*t;
    }

    // Store state and flux limits into a globally indexed array
    for (int l=0; l<numState; l++) {
      // Store the left cell edge state and flux estimates
      stateLimits(l,1,k  ,j,i) = gllState(l,0);
      fluxLimits (l,1,k  ,j,i) = gllFlux (l,0);

      // Store the Right cell edge state and flux estimates
      stateLimits(l,0,k+1,j,i) = gllState(l,tord-1);
      fluxLimits (l,0,k+1,j,i) = gllFlux (l,tord-1);
    }
  });

  // Apply boundary conditions to fluxes and state values
  edgeBoundariesZ(stateLimits, fluxLimits, dom);

  // Riemann solver
  // for (int k=0; k<dom.nz+1; k++) {
  //   for (int j=0; j<dom.ny; j++) {
  //     for (int i=0; i<dom.nx; i++) {
  yakl::parallel_for( (dom.nz+1)*dom.ny*dom.nx , YAKL_LAMBDA (int const iGlob) {
    int k, j, i;
    yakl::unpackIndices(iGlob,dom.nz+1,dom.ny,dom.nx,k,j,i);
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
  yakl::parallel_for( numState*dom.nz*dom.ny*dom.nx , YAKL_LAMBDA (int const iGlob) {
    int l, k, j, i;
    yakl::unpackIndices(iGlob,numState,dom.nz,dom.ny,dom.nx,l,k,j,i);
    tend(l,k,j,i) = - ( flux(l,k+1,j,i) - flux(l,k,j,i) ) / dom.dz;
    if (l==idRW) {
      tend(l,k,j,i) += src(k,j,i);
    }
  });
}


void TendenciesThetaConsSD::compEulerTend_S(realArr &state, Domain const &dom, Exchange &exch, Parallel const &par, realArr &tend) {
  // Form the tendencies
  // for (int k=0; k<dom.nz; k++) {
  //   for (int j=0; j<dom.ny; j++) {
  //     for (int i=0; i<dom.nx; i++) {
  yakl::parallel_for( dom.nz*dom.ny*dom.nx , YAKL_LAMBDA (int const iGlob) {
    int k, j, i;
    yakl::unpackIndices(iGlob,dom.nz,dom.ny,dom.nx,k,j,i);
    tend(idR ,k,j,i) = 0;
    tend(idRU,k,j,i) = 0;
    tend(idRV,k,j,i) = 0;
    tend(idRW,k,j,i) = -state(idR,hs+k,hs+j,hs+i) * GRAV;
    tend(idRT,k,j,i) = 0;
  });
}


void TendenciesThetaConsSD::stateBoundariesZ(realArr &state, Domain const &dom) {
  // for (int j=0; j<dom.ny; j++) {
  //   for (int i=0; i<dom.nx; i++) {
  //     for (int ii=0; ii<hs; ii++) {
  yakl::parallel_for( dom.ny*dom.nx*hs , YAKL_LAMBDA (int const iGlob) {
    int j, i, ii;
    yakl::unpackIndices(iGlob,dom.ny,dom.nx,hs,j,i,ii);
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


void TendenciesThetaConsSD::edgeBoundariesZ(realArr &stateLimits, realArr &fluxLimits, Domain const &dom) {
  // for (int j=0; j<dom.ny; j++) {
  //   for (int i=0; i<dom.nx; i++) {
  yakl::parallel_for( dom.ny*dom.nx , YAKL_LAMBDA (int const iGlob) {
    int j, i;
    yakl::unpackIndices(iGlob,dom.ny,dom.nx,j,i);
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


void TendenciesThetaConsSD::compStrakaTend(realArr &state, Domain const &dom, Exchange &exch, Parallel const &par, realArr &tend) {
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
  yakl::parallel_for( dom.nz*dom.ny*dom.nx , YAKL_LAMBDA (int const iGlob) {
    int k, j, i;
    yakl::unpackIndices(iGlob,dom.nz,dom.ny,dom.nx,k,j,i);
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


