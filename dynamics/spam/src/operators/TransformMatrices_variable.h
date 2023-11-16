
#pragma once

#include <math.h>

namespace pamc {
namespace TransformMatrices_variable {

template <index_t ord>
YAKL_INLINE void coefs_to_sten_variable(SArray<real, 1, ord + 1> const &locs,
                                        SArray<real, 2, ord, ord> &rslt) {
  // Create the Vandermonde matrix
  SArray<real, 1, ord + 1> locs_pwr;
  // Initialize power of locations
  for (int i = 0; i < ord + 1; i++) {
    locs_pwr(i) = locs(i);
  }
  // Store first column of the matrix
  for (int i = 0; i < ord; i++) {
    rslt(0, i) = 1;
  }
  for (int i = 1; i < ord; i++) {
    for (int j = 0; j < ord + 1; j++) {
      locs_pwr(j) *= locs(j);
    }
    for (int j = 0; j < ord; j++) {
      rslt(i, j) = 1. / (i + 1.) * (locs_pwr(j) - locs_pwr(j + 1)) /
                   (locs(j) - locs(j + 1));
    }
  }
}

template <index_t ord>
YAKL_INLINE void sten_to_coefs_variable(SArray<real, 1, ord + 1> const &locs,
                                        SArray<real, 2, ord, ord> &rslt) {
  using yakl::intrinsics::matinv_ge;

  // Get coefs to stencil matrix
  SArray<real, 2, ord, ord> c2s;
  coefs_to_sten_variable(locs, c2s);

  // Invert to get sten_to_coefs
  rslt = yakl::intrinsics::matinv_ge(c2s);
}

template <index_t ord>
YAKL_INLINE void
weno_lower_sten_to_coefs(SArray<real, 1, ord + 1> const &locs,
                         SArray<real, 3, (ord - 1) / 2 + 1, (ord - 1) / 2 + 1,
                                (ord - 1) / 2 + 1> &weno_recon) {
  int constexpr hs = (ord - 1) / 2;

  SArray<real, 2, hs + 1, hs + 1> recon_lo;
  SArray<real, 1, hs + 2> locs_lo;

  // Create low-order matrices
  for (int i = 0; i < hs + 1; i++) {
    for (int ii = 0; ii < hs + 2; ii++) {
      locs_lo(ii) = locs(i + ii);
    }
    sten_to_coefs_variable(locs_lo, recon_lo);
    for (int jj = 0; jj < hs + 1; jj++) {
      for (int ii = 0; ii < hs + 1; ii++) {
        weno_recon(i, jj, ii) = recon_lo(jj, ii);
      }
    }
  }
}

} // namespace TransformMatrices_variable
} // namespace pamc
