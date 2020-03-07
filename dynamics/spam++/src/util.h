
#ifndef _UTIL_H_
#define _UTIL_H_

#include "fields.h"
#include "geometry.h"

// COULD DO F (*init_cond_func)(real, real, real) for some templated arbitrary class F
// might be useful for 1/2 forms in 2D/3D where we need to return vectors...

template<int ndims, int nquadx, int nquady, int nquadz> void set_n_form(real (*init_cond_func)(real, real, real), Field<ndims> &x, Geometry<ndims, nquadx, nquady, nquadz> &geometry, int ndof) {

  if (ndims == 1) { geometry.set_edge_values(init_cond_func, x, ndof); }
  if (ndims == 2) { geometry.set_surface_values(init_cond_func, x, ndof); }
  if (ndims == 3) { geometry.set_volume_values(init_cond_func, x, ndof); }

}

template<int ndims, int nquadx, int nquady, int nquadz> void set_n_minus_1_form(real (*init_cond_func)(real, real, real), Field<ndims> &x, Geometry<ndims, nquadx, nquady, nquadz> &geometry, int ndof) {

  if (ndims == 1) { geometry.set_point_values(init_cond_func, x, ndof) ;}
  if (ndims == 2) { geometry.set_edge_values(init_cond_func, x, ndof); }
  if (ndims == 3) { geometry.set_surface_values(init_cond_func, x, ndof); }

}


#endif
