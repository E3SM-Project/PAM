
#ifndef _UNIFORM_GEOMETRY_H_
#define _UNIFORM_GEOMETRY_H_


#include "common.h"
#include "field.h"
#include "topology.h"


template<int ndims, int ic_quad_pts_x, int ic_quad_pts_y, int ic_quad_pts_z> class Geometry {

public:

  SArray<real,ic_quad_pts_x> x_quad_pts;
  SArray<real,ic_quad_pts_x> x_quad_wts;

  SArray<real,ic_quad_pts_y> y_quad_pts;
  SArray<real,ic_quad_pts_y> y_quad_wts;

  SArray<real,ic_quad_pts_z> z_quad_pts;
  SArray<real,ic_quad_pts_z> z_quad_wts;

  void initialize(Topology &topology, real xlen, real xcent, real ylen = 0.0, real ycent = 0.0, real zlen = 0.0, real ycent = 0.0);
  template <class F> YAKL_INLINE void set_zero_form_values(F const &initial_value_function, Field &field, int ndof);
  template <class F> YAKL_INLINE void set_one_form_values(F const &initial_value_function, Field &field, int ndof);
  template <class F> YAKL_INLINE void set_two_form_values(F const &initial_value_function, Field &field, int ndof);
  template <class F> YAKL_INLINE void set_three_form_values(F const &initial_value_function, Field &field, int ndof);

};


class UniformRectangularGeometry : Geometry {
public:
  real dx, dy, dz;
  real Lx, Ly, Lz;
  real xc, yc, zc;
};


void set_quad_pts_wts(int npts, SArray &pts, SArray &wts);

#endif
