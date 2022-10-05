
#pragma once

#include "awfl_const.h"

namespace awfl {

  YAKL_INLINE void reconstruct_gll_values( SArray<real,1,ord> const &stencil                     ,
                                           SArray<real,1,ngll> &gll                              ,
                                           SArray<real,2,ord,ngll> const &coefs_to_gll           ,
                                           SArray<real,2,ord,ngll> const &sten_to_gll            ,
                                           SArray<real,2,ord,ord>  const &sten_to_coefs          ,
                                           SArray<real,3,hs+1,hs+1,hs+1> const &weno_recon_lower ,
                                           SArray<real,1,hs+2> const &idl                        ,
                                           real sigma, bool doweno ) {
    if (doweno) {

      // Reconstruct values
      SArray<real,1,ord> wenoCoefs;
      weno::compute_weno_coefs<ord>( weno_recon_lower , sten_to_coefs , stencil , wenoCoefs , idl , sigma );
      // Transform ord weno coefficients into ngll GLL points
      for (int ii=0; ii<ngll; ii++) {
        real tmp = 0;
        for (int s=0; s < ord; s++) {
          tmp += coefs_to_gll(s,ii) * wenoCoefs(s);
        }
        gll(ii) = tmp;
      }

    } else {

      // Transform ord stencil cell averages into ngll GLL points
      for (int ii=0; ii<ngll; ii++) {
        real tmp = 0;
        for (int s=0; s < ord; s++) {
          tmp += sten_to_gll(s,ii) * stencil(s);
        }
        gll(ii) = tmp;
      }

    } // if doweno
  }

}


