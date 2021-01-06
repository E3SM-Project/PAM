
#pragma once

#include "const.h"
#include <math.h>
#include <Eigen/Dense>
using yakl::SArray;

namespace TransformMatrices_variable {


  template <unsigned int ord>
  inline void coefs_to_sten_variable(SArray<double,1,ord+1> const &locs ,
                                     SArray<double,2,ord,ord> &rslt) {
    // Create the Vandermonde matrix
    SArray<double,1,ord+1> locs_pwr;
    // Initialize power of locations
    for (int i=0; i < ord+1; i++) {
      locs_pwr(i) = locs(i);
    }
    // Store first column of the matrix
    for (int i=0; i < ord; i++) {
      rslt(0,i) = 1;
    }
    for (int i=1; i < ord; i++) {
      for (int j=0; j < ord+1; j++) {
        locs_pwr(j) *= locs(j);
      }
      for (int j=0; j < ord; j++) {
        rslt(i,j) = 1./(i+1.) * (locs_pwr(j) - locs_pwr(j+1)) / (locs(j)-locs(j+1));
      }
    }
  }


  template <unsigned int ord>
  inline void sten_to_coefs_variable(SArray<double,1,ord+1> const &locs ,
                                     SArray<double,2,ord,ord> &rslt) {
    // Get coefs to stencil matrix
    SArray<double,2,ord,ord> c2s;
    coefs_to_sten_variable(locs , c2s);

    // Invert to get sten_to_coefs
    Eigen::Matrix<double,ord,ord,Eigen::RowMajor> c2s_mat(c2s.data());
    auto c2s_inv_mat = c2s_mat.fullPivLu().inverse();

    // Store the Eigen matrix into an SArray matrix
    for (int j=0; j < ord; j++) {
      for (int i=0; i < ord; i++) {
        rslt(j,i) = c2s_inv_mat(j,i);
      }
    }
  }


}


