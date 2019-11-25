
#include "TendenciesThetaPrimSD.h"


void TendenciesThetaPrimSD::initialize(Domain const &dom) {
  TransformMatrices<real> trans;

  stateLimits = realArr("srcLimits" ,numState,2,dom.nz+1,dom.ny+1,dom.nx+1);
  fwaves      = realArr("fwaves"    ,numState,2,dom.nz+1,dom.ny+1,dom.nx+1);
  src         = realArr("src"       ,dom.nz,dom.ny,dom.nx);
  stateGLL    = realArr("stateGLL"  ,numState,dom.nz,dom.ny,dom.nx,tord);

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


void TendenciesThetaPrimSD::compEulerTend_X(realArr &state, Domain const &dom, Exchange &exch, Parallel const &par, realArr &tend) {
  auto &stateLimits = this->stateLimits;
  auto &fwaves      = this->fwaves     ;
  auto &src         = this->src        ;
  auto &stateGLL    = this->stateGLL   ;
  auto &gllWts      = this->gllWts     ;
  auto &to_gll      = this->to_gll     ;
  auto &wenoRecon   = this->wenoRecon  ;
  auto &wenoIdl     = this->wenoIdl    ;
  auto &wenoSigma   = this->wenoSigma  ;

  // Exchange halos in the x-direction
  exch.haloInit      ();
  exch.haloPackN_x   (dom, state, numState);
  exch.haloExchange_x(dom, par);
  exch.haloUnpackN_x (dom, state, numState);

  // Compute tend = -A*(qR - qL)/dx, and store cell-edge state vectors
  // for (int k=0; k<dom.nz; k++) {
  //   for (int j=0; j<dom.ny; j++) {
  //     for (int i=0; i<dom.nx; i++) {
  yakl::parallel_for( dom.nz*dom.ny*dom.nx , YAKL_LAMBDA ( int const iGlob ) {
    int k, j, i;
    yakl::unpackIndices(iGlob,dom.nz,dom.ny,dom.nx,k,j,i);
    SArray<real,numState,tord> gllState;  // GLL state values

    // Compute tord GLL points of the state vector
    for (int l=0; l<numState; l++) {
      SArray<real,ord> stencil;
      SArray<real,tord> gllPts;
      for (int ii=0; ii<ord; ii++) { stencil(ii) = state(l,hs+k,hs+j,i+ii); }
      reconStencil(stencil, gllPts, dom.doWeno, wenoRecon, to_gll, wenoIdl, wenoSigma);
      for (int ii=0; ii<tord; ii++) { gllState(l,ii) = gllPts(ii); }
    }
    for (int ii=0; ii<tord; ii++) {
      gllState(idR,ii) += dom.hyDensCells (hs+k);
      gllState(idT,ii) += dom.hyThetaCells(hs+k);
    }

    // Compute dq   (qR - qL)
    SArray<real,numState> dq;
    for (int l=0; l<numState; l++) {
      dq(l) = gllState(l,tord-1) - gllState(l,0);
    }

    // Compute cell-average-based values for flux Jacobian, A
    real r = state(idR,hs+k,hs+j,hs+i) + dom.hyDensCells (hs+k);
    real u = state(idU,hs+k,hs+j,hs+i);
    real t = state(idT,hs+k,hs+j,hs+i) + dom.hyThetaCells(hs+k);
    real p = C0*pow(r*t,GAMMA);
    real cs2 = GAMMA*p/r;

    // Compute tend = -A*dq/dx (A is sparse, so this is more efficient to do by hand)
    tend(0,k,j,i) = - ( u    *dq(0) + r*dq(1)                                   ) / dom.dx;
    tend(1,k,j,i) = - ( cs2/r*dq(0) + u*dq(1)                     + cs2/t*dq(4) ) / dom.dx;
    tend(2,k,j,i) = - (                       + u*dq(2)                         ) / dom.dx;
    tend(3,k,j,i) = - (                                 + u*dq(3)               ) / dom.dx;
    tend(4,k,j,i) = - (                                           + u    *dq(4) ) / dom.dx;

    // Store the state vector in stateLimits to compute fwaves from cell-interface state jumps
    for (int l=0; l<numState; l++) {
      stateLimits(l,1,k,j,i  ) = gllState(l,0     );
      stateLimits(l,0,k,j,i+1) = gllState(l,tord-1);
    }
  });

  // Reconcile the edge state via MPI exchange.
  exch.haloInit      ();
  exch.edgePackN_x   (dom, stateLimits, numState);
  exch.edgeExchange_x(dom, par);
  exch.edgeUnpackN_x (dom, stateLimits, numState);

  // Compute the fwaves from the cell interface jumps
  // for (int k=0; k<dom.nz; k++) {
  //   for (int j=0; j<dom.ny; j++) {
  //     for (int i=0; i<dom.nx+1; i++) {
  yakl::parallel_for( dom.nz*dom.ny*(dom.nx+1) , YAKL_LAMBDA (int const iGlob) {
    int k, j, i;
    yakl::unpackIndices(iGlob,dom.nz,dom.ny,dom.nx+1,k,j,i);
    // Compute averaged values for the flux Jacobian diagonalization
    real r = 0.5_fp * ( stateLimits(idR,0,k,j,i) + stateLimits(idR,1,k,j,i) );
    real u = 0.5_fp * ( stateLimits(idU,0,k,j,i) + stateLimits(idU,1,k,j,i) );
    real t = 0.5_fp * ( stateLimits(idT,0,k,j,i) + stateLimits(idT,1,k,j,i) );
    real p = C0*pow(r*t,GAMMA);
    real cs = sqrt(GAMMA*p/r);
    real cs2 = cs*cs;

    // Compute the state jump over the interface
    SArray<real,numState> dq;
    for (int l=0; l<numState; l++) {
      dq(l) = stateLimits(l,1,k,j,i) - stateLimits(l,0,k,j,i);
    }

    // Compute df = A*dq
    SArray<real,numState> df;
    df(0) = u    *dq(0) + r*dq(1)                                  ;
    df(1) = cs2/r*dq(0) + u*dq(1)                     + cs2/t*dq(4);
    df(2) =                       + u*dq(2)                        ;
    df(3) =                                 + u*dq(3)              ;
    df(4) =                                           + u    *dq(4);

    // Compute characteristic variables (L*dq)
    SArray<real,numState> ch;
    ch(0) = 0.5_fp*df(0) - r/(2*cs)*df(1) + r/(2*t)*df(4);
    ch(1) = 0.5_fp*df(0) + r/(2*cs)*df(1) + r/(2*t)*df(4);
    ch(2) =                                 -r/t   *df(4);
    ch(3) = df(2);
    ch(4) = df(3);

    // Compute fwaves
    for (int l=0; l<numState; l++) {
      fwaves(l,0,k,j,i) = 0;
      fwaves(l,1,k,j,i) = 0;
    }

    // First wave (u-cs); always negative wave speed
    fwaves(0,0,k,j,i) += ch(0);
    fwaves(1,0,k,j,i) += -cs/r*ch(0);

    // Second wave (u+cs); always positive wave speed
    fwaves(0,1,k,j,i) += ch(1);
    fwaves(1,1,k,j,i) += cs/r*ch(1);

    if (u > 0) {
      // Third wave
      fwaves(0,1,k,j,i) += ch(2);
      fwaves(4,1,k,j,i) += -t/r*ch(2);
      // Fourth wave
      fwaves(2,1,k,j,i) += ch(3);
      // Fifth Wave
      fwaves(3,1,k,j,i) += ch(4);
    } else {
      // Third wave
      fwaves(0,0,k,j,i) += ch(2);
      fwaves(4,0,k,j,i) += -t/r*ch(2);
      // Fourth wave
      fwaves(2,0,k,j,i) += ch(3);
      // Fifth Wave
      fwaves(3,0,k,j,i) += ch(4);
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


void TendenciesThetaPrimSD::compEulerTend_Y(realArr &state, Domain const &dom, Exchange &exch, Parallel const &par, realArr &tend) {
  auto &stateLimits = this->stateLimits;
  auto &fwaves      = this->fwaves     ;
  auto &src         = this->src        ;
  auto &stateGLL    = this->stateGLL   ;
  auto &gllWts      = this->gllWts     ;
  auto &to_gll      = this->to_gll     ;
  auto &wenoRecon   = this->wenoRecon  ;
  auto &wenoIdl     = this->wenoIdl    ;
  auto &wenoSigma   = this->wenoSigma  ;

  // Exchange halos in the y-direction
  exch.haloInit      ();
  exch.haloPackN_y   (dom, state, numState);
  exch.haloExchange_y(dom, par);
  exch.haloUnpackN_y (dom, state, numState);

  // Compute tend = -A*(qR - qL)/dy, and store cell-edge state vectors
  // for (int k=0; k<dom.nz; k++) {
  //   for (int j=0; j<dom.ny; j++) {
  //     for (int i=0; i<dom.nx; i++) {
  yakl::parallel_for( dom.nz*dom.ny*dom.nx , YAKL_LAMBDA (int const iGlob) {
    int k, j, i;
    yakl::unpackIndices(iGlob,dom.nz,dom.ny,dom.nx,k,j,i);
    SArray<real,numState,tord> gllState;  // GLL state values

    // Compute tord GLL points of the state vector
    for (int l=0; l<numState; l++) {
      SArray<real,ord> stencil;
      SArray<real,tord> gllPts;
      for (int ii=0; ii<ord; ii++) { stencil(ii) = state(l,hs+k,j+ii,hs+i); }
      reconStencil(stencil, gllPts, dom.doWeno, wenoRecon, to_gll, wenoIdl, wenoSigma);
      for (int ii=0; ii<tord; ii++) { gllState(l,ii) = gllPts(ii); }
    }
    for (int ii=0; ii<tord; ii++) {
      gllState(idR,ii) += dom.hyDensCells (hs+k);
      gllState(idT,ii) += dom.hyThetaCells(hs+k);
    }

    // Compute dq   (qR - qL)
    SArray<real,numState> dq;
    for (int l=0; l<numState; l++) {
      dq(l) = gllState(l,tord-1) - gllState(l,0);
    }

    // Compute cell-average-based values for flux Jacobian, A
    real r = state(idR,hs+k,hs+j,hs+i) + dom.hyDensCells (hs+k);
    real v = state(idV,hs+k,hs+j,hs+i);
    real t = state(idT,hs+k,hs+j,hs+i) + dom.hyThetaCells(hs+k);
    real p = C0*pow(r*t,GAMMA);
    real cs2 = GAMMA*p/r;

    // Compute tend = -A*dq/dx (A is sparse, so this is more efficient to do by hand)
    tend(0,k,j,i) = - ( v    *dq(0)           + r*dq(2)                         ) / dom.dy;
    tend(1,k,j,i) = - (               v*dq(1)                                   ) / dom.dy;
    tend(2,k,j,i) = - ( cs2/r*dq(0)           + v*dq(2)           + cs2/t*dq(4) ) / dom.dy;
    tend(3,k,j,i) = - (                                 + v*dq(3)               ) / dom.dy;
    tend(4,k,j,i) = - (                                           + v    *dq(4) ) / dom.dy;

    // Store the state vector in stateLimits to compute fwaves from cell-interface state jumps
    for (int l=0; l<numState; l++) {
      stateLimits(l,1,k,j  ,i) = gllState(l,0     );
      stateLimits(l,0,k,j+1,i) = gllState(l,tord-1);
    }
  });

  // Reconcile the edge state via MPI exchange.
  exch.haloInit      ();
  exch.edgePackN_y   (dom, stateLimits, numState);
  exch.edgeExchange_y(dom, par);
  exch.edgeUnpackN_y (dom, stateLimits, numState);

  // Compute the fwaves from the cell interface jumps
  // for (int k=0; k<dom.nz; k++) {
  //   for (int j=0; j<dom.ny+1; j++) {
  //     for (int i=0; i<dom.nx; i++) {
  yakl::parallel_for( dom.nz*(dom.ny+1)*dom.nx , YAKL_LAMBDA (int const iGlob) {
    int k, j, i;
    yakl::unpackIndices(iGlob,dom.nz,dom.ny+1,dom.nx,k,j,i);
    // Compute averaged values for the flux Jacobian diagonalization
    real r = 0.5_fp * ( stateLimits(idR,0,k,j,i) + stateLimits(idR,1,k,j,i) );
    real v = 0.5_fp * ( stateLimits(idV,0,k,j,i) + stateLimits(idV,1,k,j,i) );
    real t = 0.5_fp * ( stateLimits(idT,0,k,j,i) + stateLimits(idT,1,k,j,i) );
    real p = C0*pow(r*t,GAMMA);
    real cs = sqrt(GAMMA*p/r);
    real cs2 = cs*cs;

    // Compute the state jump over the interface
    SArray<real,numState> dq;
    for (int l=0; l<numState; l++) {
      dq(l) = stateLimits(l,1,k,j,i) - stateLimits(l,0,k,j,i);
    }

    // Compute df = A*dq
    SArray<real,numState> df;
    df(0) = v    *dq(0)           + r*dq(2)                        ;
    df(1) =             + v*dq(1)                                  ;
    df(2) = cs2/r*dq(0)           + v*dq(2)           + cs2/t*dq(4);
    df(3) =                                 + v*dq(3)              ;
    df(4) =                                           + v    *dq(4);

    // Compute characteristic variables (L*dq)
    SArray<real,numState> ch;
    ch(0) = 0.5_fp*df(0) - r/(2*cs)*df(2) + r/(2*t)*df(4);
    ch(1) = 0.5_fp*df(0) + r/(2*cs)*df(2) + r/(2*t)*df(4);
    ch(2) =                                 -r/t   *df(4);
    ch(3) = df(1);
    ch(4) = df(3);

    // Compute fwaves
    for (int l=0; l<numState; l++) {
      fwaves(l,0,k,j,i) = 0;
      fwaves(l,1,k,j,i) = 0;
    }

    // First wave (v-cs); always negative wave speed
    fwaves(0,0,k,j,i) += ch(0);
    fwaves(2,0,k,j,i) += -cs/r*ch(0);

    // Second wave (v+cs); always positive wave speed
    fwaves(0,1,k,j,i) += ch(1);
    fwaves(2,1,k,j,i) += cs/r*ch(1);

    if (v > 0) {
      // Third wave
      fwaves(0,1,k,j,i) += ch(2);
      fwaves(4,1,k,j,i) += -t/r*ch(2);
      // Fourth wave
      fwaves(1,1,k,j,i) += ch(3);
      // Fifth Wave
      fwaves(3,1,k,j,i) += ch(4);
    } else {
      // Third wave
      fwaves(0,0,k,j,i) += ch(2);
      fwaves(4,0,k,j,i) += -t/r*ch(2);
      // Fourth wave
      fwaves(1,0,k,j,i) += ch(3);
      // Fifth Wave
      fwaves(3,0,k,j,i) += ch(4);
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


void TendenciesThetaPrimSD::compEulerTend_Z(realArr &state, Domain const &dom, Exchange &exch, Parallel const &par, realArr &tend) {
  auto &stateLimits = this->stateLimits;
  auto &fwaves      = this->fwaves     ;
  auto &src         = this->src        ;
  auto &stateGLL    = this->stateGLL   ;
  auto &gllWts      = this->gllWts     ;
  auto &to_gll      = this->to_gll     ;
  auto &wenoRecon   = this->wenoRecon  ;
  auto &wenoIdl     = this->wenoIdl    ;
  auto &wenoSigma   = this->wenoSigma  ;

  // for (int l=0; l<numState; l++) {
  //   for (int k=0; k<dom.nz; k++) {
  //     for (int j=0; j<dom.ny; j++) {
  //       for (int i=0; i<dom.nx; i++) {
  yakl::parallel_for( numState*dom.nz*dom.ny*dom.nx , YAKL_LAMBDA (int const iGlob) {
    int l, k, j, i;
    yakl::unpackIndices(iGlob,numState,dom.nz,dom.ny,dom.nx,l,k,j,i);
    tend(l,k,j,i)  = 0;
  });

  // Apply BCs to state boundaries
  stateBoundariesZ(state, dom);

  //////////////////////////////////////////////////////////////////////////
  // COMPUTE FAST TENDENCIES
  //////////////////////////////////////////////////////////////////////////
  
  // Compute tend = -A*(qR - qL)/dz, and store cell-edge state vectors for fast solve
  // for (int k=0; k<dom.nz; k++) {
  //   for (int j=0; j<dom.ny; j++) {
  //     for (int i=0; i<dom.nx; i++) {
  yakl::parallel_for( dom.nz*dom.ny*dom.nx , YAKL_LAMBDA (int const iGlob) {
    int k, j, i;
    yakl::unpackIndices(iGlob,dom.nz,dom.ny,dom.nx,k,j,i);
    SArray<real,3,tord> gllState;  // GLL state values
    SArray<real,ord> stencil;
    SArray<real,tord> gllPts;
    // Reconstruct density
    for (int ii=0; ii<ord; ii++) { stencil(ii) = state(idR,k+ii,hs+j,hs+i); }
    reconStencil(stencil, gllPts, dom.doWeno, wenoRecon, to_gll, wenoIdl, wenoSigma);
    for (int ii=0; ii<tord; ii++) { gllState(0,ii) = gllPts(ii) + dom.hyDensGLL(k,ii); }

    // Reconstruct w
    for (int ii=0; ii<ord; ii++) { stencil(ii) = state(idW,k+ii,hs+j,hs+i); }
    reconStencil(stencil, gllPts, dom.doWeno, wenoRecon, to_gll, wenoIdl, wenoSigma);
    for (int ii=0; ii<tord; ii++) { gllState(1,ii) = gllPts(ii); }

    // Reconstruct theta
    for (int ii=0; ii<ord; ii++) { stencil(ii) = state(idT,k+ii,hs+j,hs+i); }
    reconStencil(stencil, gllPts, dom.doWeno, wenoRecon, to_gll, wenoIdl, wenoSigma);
    for (int ii=0; ii<tord; ii++) { gllState(2,ii) = gllPts(ii) + dom.hyThetaGLL(k,ii); }

    // Replace theta with perturbation pressure
    for (int ii=0; ii<tord; ii++) {
      real r = gllState(0,ii);
      real t = gllState(2,ii);
      gllState(2,ii) = C0*pow(r*t,GAMMA) - dom.hyPressureGLL(k,ii);
    }

    // Compute dq   (qR - qL)
    SArray<real,3> dq;
    for (int l=0; l<3; l++) {
      dq(l) = gllState(l,tord-1) - gllState(l,0);
    }

    // Compute cell-average-based values for flux Jacobian, A
    real r = state(idR,hs+k,hs+j,hs+i) + dom.hyDensCells(hs+k);

    // Compute tend = -A*dq/dx (A is sparse, so this is more efficient to do by hand)
    tend(idR,k,j,i) = - ( r*dq(1) ) / dom.dz;
    tend(idW,k,j,i) = - ( dq(2)/r ) / dom.dz;

    // Store the state vector in stateLimits to compute fwaves from cell-interface state jumps
    for (int l=0; l<3; l++) {
      stateLimits(l,1,k  ,j,i) = gllState(l,0     );
      stateLimits(l,0,k+1,j,i) = gllState(l,tord-1);
    }
  });

  // Enforce boundary conditions on the state limits
  // for (int j=0; j<dom.ny; j++) {
  //   for (int i=0; i<dom.nx; i++) {
  yakl::parallel_for( dom.ny*dom.nx , YAKL_LAMBDA (int const iGlob) {
    int j, i;
    yakl::unpackIndices(iGlob,dom.ny,dom.nx,j,i);
    stateLimits(0,0,0     ,j,i) = stateLimits(0,1,0     ,j,i);
    stateLimits(1,0,0     ,j,i) = 0;
    stateLimits(1,1,0     ,j,i) = 0;
    stateLimits(2,0,0     ,j,i) = stateLimits(2,1,0     ,j,i);

    stateLimits(0,1,dom.nz,j,i) = stateLimits(0,0,dom.nz,j,i);
    stateLimits(1,0,dom.nz,j,i) = 0;
    stateLimits(1,1,dom.nz,j,i) = 0;
    stateLimits(2,1,dom.nz,j,i) = stateLimits(2,0,dom.nz,j,i);
  });

  // Compute the fwaves from the cell interface jumps
  // for (int k=0; k<dom.nz+1; k++) {
  //   for (int j=0; j<dom.ny; j++) {
  //     for (int i=0; i<dom.nx; i++) {
  yakl::parallel_for( (dom.nz+1)*dom.ny*dom.nx , YAKL_LAMBDA (int const iGlob) {
    int k, j, i;
    yakl::unpackIndices(iGlob,dom.nz+1,dom.ny,dom.nx,k,j,i);
    // Compute averaged values for the flux Jacobian diagonalization
    real r = 0.5_fp * ( stateLimits(0,0,k,j,i) + stateLimits(0,1,k,j,i) );
    real p = 0.5_fp * ( stateLimits(2,0,k,j,i) + stateLimits(2,1,k,j,i) );
    if (k < dom.nz) {
      p += dom.hyPressureGLL(k,0);
    } else {
      p += dom.hyPressureGLL(dom.nz-1,tord-1);
    }
    real cs = sqrt(GAMMA*p/r);
    real cs2 = cs*cs;

    // Compute the state jump over the interface
    SArray<real,3> dq;
    for (int l=0; l<3; l++) {
      dq(l) = stateLimits(l,1,k,j,i) - stateLimits(l,0,k,j,i);
    }

    // Compute df = A*dq
    SArray<real,3> df;
    df(0) = r*dq(1);
    df(1) = dq(2)/r;
    df(2) = r*cs2*dq(1);

    // Compute characteristic variables (L*df)
    SArray<real,2> ch;
    ch(0) = -r*df(1)/(2*cs) + df(2)/(2*cs2);
    ch(1) =  r*df(1)/(2*cs) + df(2)/(2*cs2);

    // Compute fwaves
    // First wave (-cs); always negative wave speed
    fwaves(0,0,k,j,i) = (ch(0));
    fwaves(1,0,k,j,i) = (-cs*ch(0)/r);

    // Second wave (+cs); always positive wave speed
    fwaves(0,1,k,j,i) = (ch(1));
    fwaves(1,1,k,j,i) = ( cs*ch(1)/r);
  });

  // Apply the fwaves to the tendencies
  // for (int k=0; k<dom.nz; k++) {
  //   for (int j=0; j<dom.ny; j++) {
  //     for (int i=0; i<dom.nx; i++) {
  yakl::parallel_for( dom.nz*dom.ny*dom.nx , YAKL_LAMBDA (int const iGlob) {
    int k, j, i;
    yakl::unpackIndices(iGlob,dom.nz,dom.ny,dom.nx,k,j,i);
    tend(idR,k,j,i) += - ( fwaves(0,1,k,j,i) + fwaves(0,0,k+1,j,i) ) / dom.dz;
    tend(idW,k,j,i) += - ( fwaves(1,1,k,j,i) + fwaves(1,0,k+1,j,i) ) / dom.dz;
  });

  //////////////////////////////////////////////////////////////////////////
  // COMPUTE SLOOOOWWWW TENDENCIES
  //////////////////////////////////////////////////////////////////////////
  
  // Compute tend = -A*(qR - qL)/dz, and store cell-edge state vectors for fast solve
  // for (int k=0; k<dom.nz; k++) {
  //   for (int j=0; j<dom.ny; j++) {
  //     for (int i=0; i<dom.nx; i++) {
  yakl::parallel_for( dom.nz*dom.ny*dom.nx , YAKL_LAMBDA (int const iGlob) {
    int k, j, i;
    yakl::unpackIndices(iGlob,dom.nz,dom.ny,dom.nx,k,j,i);
    SArray<real,numState,tord> gllState;  // GLL state values

    // Compute tord GLL points of the state vector
    for (int l=0; l<numState; l++) {
      SArray<real,ord> stencil;
      SArray<real,tord> gllPts;
      for (int ii=0; ii<ord; ii++) { stencil(ii) = state(l,k+ii,hs+j,hs+i); }
      reconStencil(stencil, gllPts, dom.doWeno, wenoRecon, to_gll, wenoIdl, wenoSigma);
      for (int ii=0; ii<tord; ii++) { gllState(l,ii) = gllPts(ii); }
    }
    for (int ii=0; ii<tord; ii++) {
      gllState(idR,ii) += dom.hyDensGLL (k,ii);
      gllState(idT,ii) += dom.hyThetaGLL(k,ii);
    }

    // Compute dq   (qR - qL)
    SArray<real,numState> dq;
    for (int l=0; l<numState; l++) {
      dq(l) = gllState(l,tord-1) - gllState(l,0);
    }

    // Compute cell-average-based values for flux Jacobian, A
    real w = state(idW,hs+k,hs+j,hs+i);

    // Compute tend = -A*dq/dx (A is sparse, so this is more efficient to do by hand)
    tend(0,k,j,i) += - ( w*dq(0) ) / dom.dz;
    tend(1,k,j,i) += - ( w*dq(1) ) / dom.dz;
    tend(2,k,j,i) += - ( w*dq(2) ) / dom.dz;
    tend(3,k,j,i) += - ( w*dq(3) ) / dom.dz;
    tend(4,k,j,i) += - ( w*dq(4) ) / dom.dz;

    // Store the state vector in stateLimits to compute fwaves from cell-interface state jumps
    for (int l=0; l<numState; l++) {
      stateLimits(l,1,k  ,j,i) = gllState(l,0     );
      stateLimits(l,0,k+1,j,i) = gllState(l,tord-1);
    }
  });

  // Enforce boundary conditions
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

  // Compute the fwaves from the cell interface jumps
  // for (int k=0; k<dom.nz+1; k++) {
  //   for (int j=0; j<dom.ny; j++) {
  //     for (int i=0; i<dom.nx; i++) {
  yakl::parallel_for( (dom.nz+1)*dom.ny*dom.nx , YAKL_LAMBDA (int const iGlob) {
    int k, j, i;
    yakl::unpackIndices(iGlob,dom.nz+1,dom.ny,dom.nx,k,j,i);
    // Compute averaged values for the flux Jacobian diagonalization
    real w = 0.5_fp * ( stateLimits(idW,0,k,j,i) + stateLimits(idW,1,k,j,i) );

    // Compute the state jump over the interface
    SArray<real,numState> dq;
    for (int l=0; l<numState; l++) {
      dq(l) = stateLimits(l,1,k,j,i) - stateLimits(l,0,k,j,i);
    }

    // Compute df = A*dq
    SArray<real,numState> df;
    df(0) = w*dq(0);
    df(1) = w*dq(1);
    df(2) = w*dq(2);
    df(3) = w*dq(3);
    df(4) = w*dq(4);

    // Compute fwaves
    for (int l=0; l<numState; l++) {
      fwaves(l,0,k,j,i) = 0;
      fwaves(l,1,k,j,i) = 0;
    }

    if (w > 0) {
      fwaves(0,1,k,j,i) += df(0);
      fwaves(1,1,k,j,i) += df(1);
      fwaves(2,1,k,j,i) += df(2);
      fwaves(3,1,k,j,i) += df(3);
      fwaves(4,1,k,j,i) += df(4);
    } else {
      fwaves(0,0,k,j,i) += df(0);
      fwaves(1,0,k,j,i) += df(1);
      fwaves(2,0,k,j,i) += df(2);
      fwaves(3,0,k,j,i) += df(3);
      fwaves(4,0,k,j,i) += df(4);
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


void TendenciesThetaPrimSD::compEulerTend_S(realArr &state, Domain const &dom, Exchange &exch, Parallel const &par, realArr &tend) {
  // Form the tendencies
  // for (int k=0; k<dom.nz; k++) {
  //   for (int j=0; j<dom.ny; j++) {
  //     for (int i=0; i<dom.nx; i++) {
  yakl::parallel_for( dom.nz*dom.ny*dom.nx , YAKL_LAMBDA (int const iGlob) {
    int k, j, i;
    yakl::unpackIndices(iGlob,dom.nz,dom.ny,dom.nx,k,j,i);
    tend(idR ,k,j,i) = 0;
    tend(idU,k,j,i) = 0;
    tend(idV,k,j,i) = 0;
    tend(idW,k,j,i) = -state(idR,hs+k,hs+j,hs+i) / (state(idR,hs+k,hs+j,hs+i) + dom.hyDensCells(hs+k)) * GRAV;
    tend(idT,k,j,i) = 0;
  });
}


void TendenciesThetaPrimSD::stateBoundariesZ(realArr &state, Domain const &dom) {
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


void TendenciesThetaPrimSD::edgeBoundariesZ(realArr &stateLimits, realArr &fluxLimits, Domain const &dom) {
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


void TendenciesThetaPrimSD::compStrakaTend(realArr &state, Domain const &dom, Exchange &exch, Parallel const &par, realArr &tend) {
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


