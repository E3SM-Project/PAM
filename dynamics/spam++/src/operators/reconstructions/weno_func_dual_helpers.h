
#pragma once

#include "common.h"
#include "TransformMatrices.h"

// This catches the case of dual_reconstruction_order being even (ie using CFV), which otherwise breaks
uint constexpr dual_ord = ((dual_reconstruction_order % 2) == 0) ? 3 : dual_reconstruction_order;
uint constexpr dual_hs       = (dual_ord-1)/2;
uint constexpr dual_tord = 2;

YAKL_INLINE void map_weights_dual( SArray<real,dual_hs+2> const &idl , SArray<real,dual_hs+2> &wts ) {
  // Map the weights for quicker convergence. WARNING: Ideal weights must be (0,1) before mapping
  for (int i=0; i<dual_hs+2; i++) {
    wts(i) = wts(i) * ( idl(i) + idl(i)*idl(i) - 3._fp*idl(i)*wts(i) + wts(i)*wts(i) ) / ( idl(i)*idl(i) + wts(i) * ( 1._fp - 2._fp * idl(i) ) );
  }
}



YAKL_INLINE void convexify_dual( SArray<real,dual_hs+2> &wts ) {
  real sum = 0._fp;
  real const eps = 1.0e-20;
  for (int i=0; i<dual_hs+2; i++) { sum += wts(i); }
  for (int i=0; i<dual_hs+2; i++) { wts(i) /= (sum + eps); }
}



YAKL_INLINE void wenoSetIdealSigma_dual(SArray<real,dual_hs+2> &idl, real &sigma) {
  if        (dual_ord == 3) {
        sigma = 0.0343557947899881_fp;
        idl(0) = 1._fp;
        idl(1) = 1._fp;
        idl(2) = 1224.61619926508_fp;
      } else if (dual_ord == 5) {
        sigma = 0.73564225445964_fp;
        idl(0) = 1._fp;
        idl(1) = 73.564225445964_fp;
        idl(2) = 1._fp;
        idl(3) = 1584.89319246111_fp;
      } else if (dual_ord == 7) {
        sigma = 0.125594321575479_fp;
        idl(0) = 1._fp;
        idl(1) = 7.35642254459641_fp;
        idl(2) = 7.35642254459641_fp;
        idl(3) = 1._fp;
        idl(4) = 794.328234724281_fp;
      } else if (dual_ord == 9) {
        sigma = 0.0288539981181442_fp;
        idl(0) = 1._fp;
        idl(1) = 2.15766927997459_fp;
        idl(2) = 2.40224886796286_fp;
        idl(3) = 2.15766927997459_fp;
        idl(4) = 1._fp;
        idl(5) = 1136.12697719888_fp;
  } else if (dual_ord == 11) {
    // These aren't tuned!!!
    sigma = 0.1_fp;
    idl(0) = 1._fp;
    idl(1) = 1._fp;
    idl(2) = 1._fp;
    idl(3) = 1._fp;
    idl(4) = 1._fp;
    idl(5) = 1._fp;
    idl(6) = 1._fp;
  } else if (dual_ord == 13) {
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
  } else if (dual_ord == 15) {
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
  convexify_dual( idl );
}



YAKL_INLINE void perform_weno_recon_dual( SArray<real,dual_ord,dual_ord,dual_ord> const &recon , SArray<real,dual_ord> const &u , SArray<real,dual_hs+2> const &idl , SArray<real,dual_hs+2,dual_ord> &a ) {
  // Init to zero
  for (int j=0; j<dual_hs+2; j++) {
    for (int i=0; i<dual_ord; i++) {
      a(j,i) = 0._fp;
    }
  }

  // Compute three quadratic polynomials (left, center, and right) and the high-order polynomial
  for(int i=0; i<dual_hs+1; i++) {
    for (int ii=0; ii<dual_hs+1; ii++) {
      for (int s=0; s<dual_hs+1; s++) {
        a(i,ii) += recon(i,s,ii) * u(i+s);
      }
    }
  }
  for (int ii=0; ii<dual_ord; ii++) {
    for (int s=0; s<dual_ord; s++) {
      a(dual_hs+1,ii) += recon(dual_hs+1,s,ii) * u(s);
    }
  }

  // Compute "bridge" polynomial
  for (int i=0; i<dual_hs+1; i++) {
    for (int ii=0; ii<dual_hs+1; ii++) {
      a(dual_hs+1,ii) -= idl(i)*a(i,ii);
    }
  }
  for (int ii=0; ii<dual_ord; ii++) {
    a(dual_hs+1,ii) /= idl(dual_hs+1);
  }
}



YAKL_INLINE void compute_weno_weights_dual( SArray<real,dual_hs+2,dual_ord> const &a , SArray<real,dual_hs+2> const &idl , real const sigma , SArray<real,dual_hs+2> &wts ) {
  SArray<real,dual_hs+2> tv;
  SArray<real,dual_hs+1> lotmp;
  SArray<real,dual_ord > hitmp;
  real lo_avg;
  real const eps = 1.0e-20;

  TransformMatrices<real> transform;

  // Compute total variation of all candidate polynomials
  for (int i=0; i<dual_hs+1; i++) {
    for (int ii=0; ii<dual_hs+1; ii++) {
      lotmp(ii) = a(i,ii);
    }
    tv(i) = transform.coefs_to_tv(lotmp);
  }
  for (int ii=0; ii<dual_ord; ii++) {
    hitmp(ii) = a(dual_hs+1,ii);
  }
  tv(dual_hs+1) = transform.coefs_to_tv(hitmp);

  // Reduce the bridge polynomial TV to something closer to the other TV values
  lo_avg = 0._fp;
  for (int i=0; i<dual_hs+1; i++) {
    lo_avg += tv(i);
  }
  lo_avg /= dual_hs+1;
  tv(dual_hs+1) = lo_avg + ( tv(dual_hs+1) - lo_avg ) * sigma;

  // WENO weights are proportional to the inverse of TV**2 and then re-confexified
  for (int i=0; i<dual_hs+2; i++) {
    wts(i) = idl(i) / ( tv(i)*tv(i) + eps );
  }
  convexify_dual(wts);

  // Map WENO weights for sharper fronts and less sensitivity to "eps"
  map_weights_dual(idl,wts);
  convexify_dual(wts);
}



YAKL_INLINE void apply_weno_weights_dual( SArray<real,dual_hs+2,dual_ord> const &a , SArray<real,dual_hs+2> const &wts , SArray<real,dual_ord> &aw ) {
  // WENO polynomial is the weighted sum of candidate polynomials using WENO weights instead of ideal weights
  for (int i=0; i<dual_ord; i++) {
    aw(i) = 0._fp;
  }
  for (int i=0; i<dual_hs+2; i++) {
    for (int ii=0; ii<dual_ord; ii++) {
      aw(ii) += wts(i) * a(i,ii);
    }
  }
}



YAKL_INLINE void compute_weno_coefs_dual( SArray<real,dual_ord,dual_ord,dual_ord> const &recon , SArray<real,dual_ord> const &u , SArray<real,dual_ord> &aw , SArray<real,dual_hs+2> const &idl , real const sigma ) {
  SArray<real,dual_hs+2,dual_ord> a;
  SArray<real,dual_hs+2> wts;
  perform_weno_recon_dual( recon , u , idl , a );
  compute_weno_weights_dual( a , idl , sigma , wts );
  apply_weno_weights_dual( a , wts , aw );
}



YAKL_INLINE void weno_recon_and_weights_dual( SArray<real,dual_ord,dual_ord,dual_ord> const &recon , SArray<real,dual_ord> const &u , SArray<real,dual_hs+2> const &idl , real const sigma , SArray<real,dual_hs+2> &wts ) {
  SArray<real,dual_hs+2,dual_ord> a;
  perform_weno_recon_dual( recon , u , idl , a );
  compute_weno_weights_dual( a , idl , sigma , wts );
}



YAKL_INLINE void weno_recon_and_apply_dual( SArray<real,dual_ord,dual_ord,dual_ord> const &recon , SArray<real,dual_ord> const &u , SArray<real,dual_hs+2> const &idl , SArray<real,dual_hs+2> const &wts , SArray<real,dual_ord> &aw ) {
  SArray<real,dual_hs+2,dual_ord> a;
  perform_weno_recon_dual( recon , u , idl , a );
  apply_weno_weights_dual( a , wts , aw );
}



// Transform ord stencil cell averages into tord GLL point values
YAKL_INLINE void reconStencil_dual(SArray<real,dual_ord> const &stencil, SArray<real,dual_tord> &gll,
                              SArray<real,dual_ord,dual_ord,dual_ord> const &wenoRecon, SArray<real,dual_ord,dual_tord> const &to_gll,
                              SArray<real,dual_hs+2> const &wenoIdl, real wenoSigma) {
  SArray<real,dual_ord> coefs;
    compute_weno_coefs_dual(wenoRecon,stencil,coefs,wenoIdl,wenoSigma);

  for (int ii=0; ii<dual_tord; ii++) {
    gll(ii) = 0.;
    for (int s=0; s<dual_ord; s++) {
      gll(ii) += to_gll(s,ii) * coefs(s);
    }
  }
}
