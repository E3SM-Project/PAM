#pragma once

#include "common.h"
#include "TransformMatrices.h"

template<uint ord, uint tord=2, uint hs=(ord-1)/2> YAKL_INLINE void map_weights( SArray<real,1,hs+2> const &idl , SArray<real,1,hs+2> &wts ) {
  // Map the weights for quicker convergence. WARNING: Ideal weights must be (0,1) before mapping
  for (int i=0; i<hs+2; i++) {
    wts(i) = wts(i) * ( idl(i) + idl(i)*idl(i) - 3._fp*idl(i)*wts(i) + wts(i)*wts(i) ) / ( idl(i)*idl(i) + wts(i) * ( 1._fp - 2._fp * idl(i) ) );
  }
}


template<uint ord, uint tord=2, uint hs=(ord-1)/2> YAKL_INLINE void convexify( SArray<real,1,hs+2> &wts ) {
  real sum = 0._fp;
  real const eps = 1.0e-20;
  for (int i=0; i<hs+2; i++) { sum += wts(i); }
  for (int i=0; i<hs+2; i++) { wts(i) /= (sum + eps); }
}


template<uint ord, uint tord=2, uint hs=(ord-1)/2> YAKL_INLINE void wenoSetIdealSigma(SArray<real,1,hs+2> &idl, real &sigma) {
//   if        (ord == 3) {
//   sigma = 0.1_fp;
//   idl(0) = 1._fp;
//   idl(1) = 1._fp;
//   idl(2) = 100._fp;
// } else if (ord == 5) {
//   sigma = 0.1_fp;
//   idl(0) = 1._fp;
//   idl(1) = 100._fp;
//   idl(2) = 1._fp;
//   idl(3) = 1000._fp;
// } else if (ord == 7) {
//   sigma = 0.01_fp;
//   idl(0) = 1._fp;
//   idl(1) = 20._fp;
//   idl(2) = 20._fp;
//   idl(3) = 1._fp;
//   idl(4) = 400._fp;
// } else if (ord == 9) {
//   sigma = 0.1_fp;
//   idl(0) = 1._fp;
//   idl(1) = 18._fp;
//   idl(2) = 76._fp;
//   idl(3) = 18._fp;
//   idl(4) = 1._fp;
//   idl(5) = 5832._fp;

// These are new tunings from Matt
  if        (ord == 3) {
        sigma = 0.0343557947899881_fp;
        idl(0) = 1._fp;
        idl(1) = 1._fp;
        idl(2) = 1224.61619926508_fp;
      } else if (ord == 5) {
        sigma = 0.73564225445964_fp;
        idl(0) = 1._fp;
        idl(1) = 73.564225445964_fp;
        idl(2) = 1._fp;
        idl(3) = 1584.89319246111_fp;
      } else if (ord == 7) {
        sigma = 0.125594321575479_fp;
        idl(0) = 1._fp;
        idl(1) = 7.35642254459641_fp;
        idl(2) = 7.35642254459641_fp;
        idl(3) = 1._fp;
        idl(4) = 794.328234724281_fp;
      } else if (ord == 9) {
        sigma = 0.0288539981181442_fp;
        idl(0) = 1._fp;
        idl(1) = 2.15766927997459_fp;
        idl(2) = 2.40224886796286_fp;
        idl(3) = 2.15766927997459_fp;
        idl(4) = 1._fp;
        idl(5) = 1136.12697719888_fp;
  } else if (ord == 11) {
    // These aren't tuned!!!
    sigma = 0.1_fp;
    idl(0) = 1._fp;
    idl(1) = 1._fp;
    idl(2) = 1._fp;
    idl(3) = 1._fp;
    idl(4) = 1._fp;
    idl(5) = 1._fp;
    idl(6) = 1._fp;
  } else if (ord == 13) {
    // These aren't tuned!!!
    sigma = 0.1_fp;
    idl(0) = 1._fp;
    idl(1) = 1._fp;
    idl(2) = 1._fp;
    idl(3) = 1._fp;
    idl(4) = 1._fp;
    idl(5) = 1._fp;
    idl(6) = 1._fp;
    idl(7) = 1._fp;
  } else if (ord == 15) {
    // These aren't tuned!!!
    sigma = 0.1_fp;
    idl(0) = 1._fp;
    idl(1) = 1._fp;
    idl(2) = 1._fp;
    idl(3) = 1._fp;
    idl(4) = 1._fp;
    idl(5) = 1._fp;
    idl(6) = 1._fp;
    idl(7) = 1._fp;
    idl(8) = 1._fp;
  }
  convexify<ord>( idl );
}



template<uint ord, uint tord=2, uint hs=(ord-1)/2> YAKL_INLINE void perform_weno_recon( SArray<real,3,ord,ord,ord> const &recon , SArray<real,1,ord> const &u , SArray<real,1,hs+2> const &idl , SArray<real,2,hs+2,ord> &a ) {
  // Init to zero
  for (int j=0; j<hs+2; j++) {
    for (int i=0; i<ord; i++) {
      a(j,i) = 0._fp;
    }
  }

  // Compute three quadratic polynomials (left, center, and right) and the high-order polynomial
  for(int i=0; i<hs+1; i++) {
    for (int ii=0; ii<hs+1; ii++) {
      for (int s=0; s<hs+1; s++) {
        a(i,ii) += recon(i,s,ii) * u(i+s);
      }
    }
  }
  for (int ii=0; ii<ord; ii++) {
    for (int s=0; s<ord; s++) {
      a(hs+1,ii) += recon(hs+1,s,ii) * u(s);
    }
  }

  // Compute "bridge" polynomial
  for (int i=0; i<hs+1; i++) {
    for (int ii=0; ii<hs+1; ii++) {
      a(hs+1,ii) -= idl(i)*a(i,ii);
    }
  }
  for (int ii=0; ii<ord; ii++) {
    a(hs+1,ii) /= idl(hs+1);
  }
}



template<uint ord, uint tord=2, uint hs=(ord-1)/2> YAKL_INLINE void compute_weno_weights( SArray<real,2,hs+2,ord> const &a , SArray<real,1,hs+2> const &idl , real const sigma , SArray<real,1,hs+2> &wts ) {
  SArray<real,1,hs+2> tv;
  SArray<real,1,hs+1> lotmp;
  SArray<real,1,ord > hitmp;
  real lo_avg;
  real const eps = 1.0e-20;

  //TransformMatrices<real> transform;

  // Compute total variation of all candidate polynomials
  for (int i=0; i<hs+1; i++) {
    for (int ii=0; ii<hs+1; ii++) {
      lotmp(ii) = a(i,ii);
    }
    tv(i) = TransformMatrices::coefs_to_tv(lotmp);
  }
  for (int ii=0; ii<ord; ii++) {
    hitmp(ii) = a(hs+1,ii);
  }
  tv(hs+1) = TransformMatrices::coefs_to_tv(hitmp);

  // Reduce the bridge polynomial TV to something closer to the other TV values
  lo_avg = 0._fp;
  for (int i=0; i<hs+1; i++) {
    lo_avg += tv(i);
  }
  lo_avg /= hs+1;
  tv(hs+1) = lo_avg + ( tv(hs+1) - lo_avg ) * sigma;

  // WENO weights are proportional to the inverse of TV**2 and then re-confexified
  for (int i=0; i<hs+2; i++) {
    wts(i) = idl(i) / ( tv(i)*tv(i) + eps );
  }
  convexify<ord>(wts);

  // Map WENO weights for sharper fronts and less sensitivity to "eps"
  map_weights<ord>(idl,wts);
  convexify<ord>(wts);
}



template<uint ord, uint tord=2, uint hs=(ord-1)/2> YAKL_INLINE void apply_weno_weights( SArray<real,2,hs+2,ord> const &a , SArray<real,1,hs+2> const &wts , SArray<real,1,ord> &aw ) {
  // WENO polynomial is the weighted sum of candidate polynomials using WENO weights instead of ideal weights
  for (int i=0; i<ord; i++) {
    aw(i) = 0._fp;
  }
  for (int i=0; i<hs+2; i++) {
    for (int ii=0; ii<ord; ii++) {
      aw(ii) += wts(i) * a(i,ii);
    }
  }
}



template<uint ord, uint tord=2, uint hs=(ord-1)/2> YAKL_INLINE void compute_weno_coefs( SArray<real,3,ord,ord,ord> const &recon , SArray<real,1,ord> const &u , SArray<real,1,ord> &aw , SArray<real,1,hs+2> const &idl , real const sigma ) {
  SArray<real,2,hs+2,ord> a;
  SArray<real,1,hs+2> wts;
  perform_weno_recon( recon , u , idl , a );
  compute_weno_weights( a , idl , sigma , wts );
  apply_weno_weights( a , wts , aw );
}



template<uint ord, uint tord=2, uint hs=(ord-1)/2> YAKL_INLINE void weno_recon_and_weights( SArray<real,3,ord,ord,ord> const &recon , SArray<real,1,ord> const &u , SArray<real,1,hs+2> const &idl , real const sigma , SArray<real,1,hs+2> &wts ) {
  SArray<real,2,hs+2,ord> a;
  perform_weno_recon( recon , u , idl , a );
  compute_weno_weights( a , idl , sigma , wts );
}



template<uint ord, uint tord=2, uint hs=(ord-1)/2> YAKL_INLINE void weno_recon_and_apply( SArray<real,3,ord,ord,ord> const &recon , SArray<real,1,ord> const &u , SArray<real,1,hs+2> const &idl , SArray<real,1,hs+2> const &wts , SArray<real,1,ord> &aw ) {
  SArray<real,2,hs+2,ord> a;
  perform_weno_recon( recon , u , idl , a );
  apply_weno_weights( a , wts , aw );
}



// Transform ord stencil cell averages into tord GLL point values
template<uint ord, uint tord=2, uint hs=(ord-1)/2> YAKL_INLINE void reconStencil(SArray<real,1,ord> const &stencil, SArray<real,1,tord> &gll,
                              SArray<real,3,ord,ord,ord> const &wenoRecon, SArray<real,2,ord,tord> const &to_gll,
                              SArray<real,1,hs+2> const &wenoIdl, real wenoSigma) {
  SArray<real,1,ord> coefs;
    compute_weno_coefs<ord>(wenoRecon,stencil,coefs,wenoIdl,wenoSigma);

  for (int ii=0; ii<tord; ii++) {
    gll(ii) = 0.;
    for (int s=0; s<ord; s++) {
      gll(ii) += to_gll(s,ii) * coefs(s);
    }
  }
}


template<uint ndofs, uint nd, uint ord, uint tord=2, uint hs=(ord-1)/2> void YAKL_INLINE weno_func(SArray<real,3,ndofs,nd,2> &edgerecon, SArray<real,3,ndofs,nd,ord> const &dens,
  SArray<real,3,ord,ord,ord> const &wenoRecon, SArray<real,2,ord,2> const &to_gll, SArray<real,1,hs+2> const &wenoIdl, real wenoSigma) {

  SArray<real,1,ord> stencil;
  SArray<real,1,2> gllPts;

    for (int l=0; l<ndofs; l++) {
      for (int d=0; d<nd; d++) {

      for (int ii=0; ii<ord; ii++) { stencil(ii) = dens(l,d,ii); }
      reconStencil<ord>(stencil, gllPts, wenoRecon, to_gll, wenoIdl, wenoSigma);
       edgerecon(l,d,0) = gllPts(0);
       edgerecon(l,d,1) = gllPts(1);
     }
      }
}