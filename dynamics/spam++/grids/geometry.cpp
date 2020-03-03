




#include "uniform_geometry.h"

  void set_quad_pts_wts(int npts, SArray &pts, SArray &wts) {

    if (npts == 1)
    {
      pts(0) = 0.5;
      wts(0) = 1.0;
    }

    if (npts == 3)
    {
      pts(0) = 0.112701665379258311482073460022;
      pts(1) = 0.500000000000000000000000000000;
      pts(2) = 0.887298334620741688517926539980;

      wts(0) = 0.277777777777777777777777777779;
      wts(1) = 0.444444444444444444444444444444;
      wts(2) = 0.277777777777777777777777777779;
    }


  }


  void UniformRectangularGeometry::initialize(Topology &topology, real xlen, real xcent, real ylen = 0.0, real ycent = 0.0, real zlen = 0.0, real ycent = 0.0) {

  Lx = xlen;
  Ly = ylen;
  Lz = zlen;

  xc = xcent;
  yc = ycent;
  zc = zcent;

  dx = Lx/topology.n_cells_x;
  dy = Ly/topology.n_cells_y;
  dz = Lz/topology.n_cells_z;

  if (ndims == 1)
  {
    Ly = 0.0;
    Lz = 0.0;
    yc = 0.0;
    zc = 0.0;
    dy = 0.0;
    dz = 0.0;
  }

  if (ndims == 2)
  {
    Lz = 0.0;
    zc = 0.0;
    dz = 0.0;
  }

  set_quad_pts_wts(ic_quad_pts_x, x_quad_pts, x_quad_wts);
  set_quad_pts_wts(ic_quad_pts_y, y_quad_pts, y_quad_wts);
  set_quad_pts_wts(ic_quad_pts_z, z_quad_pts, z_quad_wts);
}


// We assume here that each degree of freedom is a discrete differential form
// 0 forms are sampled at vertices
// 1 forms are integrated along lines
// 2 forms are integrated over faces
// 3 forms are integrated over volumes

// 0 forms are sampled at vertices
// 1 forms are integrated along lines or over faces (primal vs. dual)
// 2 forms are integrated over volumes

// 0 forms are sampled at vertices
// 1 forms are integrated over volumes

template <class F> YAKL_INLINE void UniformRectangularGeometry::set_zero_form_values(F const &initial_value_function, Field &field, int ndof) {
  int is = topology.is;
  int js = topology.js;
  int ks = topology.ks;

  yakl::parallel_for("SetZeroForm", topology.n_cells, YAKL_LAMBDA (int iGlob) {
    int k, j, i;
    yakl::unpackIndices(iGlob, topology.n_cells_z, topology.n_cells_y, topology.n_cells_x ,k, j, i);
      field.data(ndof, k+ks, j+js, i+is) = initial_value_function(i*dx + xc - Lx/2., j*dy + yc - Ly/2., k*dz + zc - Lz/2.);
  });
}

template <class F> YAKL_INLINE void UniformRectangularGeometry::set_three_form_values(F const &initial_value_function, Field &field, int ndof) {
  int is = topology.is;
  int js = topology.js;
  int ks = topology.ks;
  real ll_corner_x, ll_corner_y, ll_corner_z;
  real tempval;
  int offset = field.get_offset(3);

  yakl::parallel_for("SetThreeForm", topology.n_cells, YAKL_LAMBDA (int iGlob) {
    int k, j, i;
    yakl::unpackIndices(iGlob, topology.n_cells_z, topology.n_cells_y, topology.n_cells_x ,k, j, i);
    tempval = 0.0;

    ll_corner_x = i*dx + xc - Lx/2.;
    ll_corner_y = j*dy + yc - Ly/2.;
    ll_corner_z = k*dz + zc - Lz/2.;
    for (int nqx=0; nqx<ic_quad_pts_x; nqx++) {
      for (int nqy=0; nqy<ic_quad_pts_y; nqy++) {
        for (int nqz=0; nqz<ic_quad_pts_z; nqz++) {
          tempval = tempval + initial_value_function(ll_corner_x + x_quad_pts(nqx)*dx, ll_corner_y + y_quad_pts(nqy)*dy, ll_corner_z + z_quad_pts(nqz)*dz) * x_quad_wts(nqx) * y_quad_wts(nqy) * z_quad_wts(nqz);
          } } }
      field.data(ndof+offset, k+ks, j+js, i+is) = tempval;
  });
}



// NOT SURE EDGE INDEXING IS CORRECT HERE!
template <class F> YAKL_INLINE void UniformRectangularGeometry::set_one_form_values(F const &initial_value_function, Field &field, int ndof) {
  int is = topology.is;
  int js = topology.js;
  int ks = topology.ks;
  real ll_corner_x, ll_corner_y, ll_corner_z;
  real tempval;
  int offset = field.get_offset(1);

  yakl::parallel_for("SetOneForm", topology.n_cells , YAKL_LAMBDA (int iGlob) {
    int k, j, i;
    yakl::unpackIndices(iGlob, topology.n_cells_z, topology.n_cells_y, topology.n_cells_x ,k, j, i);
    tempval = 0.0;

    ll_corner_x = i*dx + xc - Lx/2.;
    ll_corner_y = j*dy + yc - Ly/2.;
    ll_corner_z = k*dz + zc - Lz/2.;

    if (ndim == 1) {
      ll_corner_y = 0.0;
      ll_corner_z = 0.0;
    }

    if (ndim == 2) {
        ll_corner_z = 0.0;
    }

        // x-edge
        for (int nqx=0; nqx<ic_quad_pts_x; nqx++) {
          tempval = tempval + initial_value_function(ll_corner_x + x_quad_pts(nqx)*dx, ll_corner_y, ll_corner_z) * x_quad_wts(nqx);
        }
        field.data(ndof+offset, k+ks, j+js, i+is) = tempval;

        // y-edge
        if (ndim >=2) {
        for (int nqy=0; nqy<ic_quad_pts_y; nqy++) {
          tempval = tempval + initial_value_function(ll_corner_x, ll_corner_y + y_quad_pts(nqy)*dy, ll_corner_z) * y_quad_wts(nqy);
        }
        field.data(ndof+field.ndofs1+offset, k+ks, j+js, i+is) = tempval;
        }

        if (ndim ==3) {
        // z-edge
        for (int nqz=0; nqy<ic_quad_pts_z; nqz++) {
          tempval = tempval + initial_value_function(ll_corner_x, ll_corner_y, ll_corner_z + z_quad_pts(nqz)*dz) * z_quad_wts(nqz);
        }
        field.data(ndof+2*field.ndofs1+offset, k+ks, j+js, i+is) = tempval;
      }

  });
}


// NOT SURE SURFACE INDEXING IS CORRECT HERE!
template <class F> YAKL_INLINE void set_two_form_values(F const &initial_value_function, Field &field, int ndof) {
  int is = topology.is;
  int js = topology.js;
  int ks = topology.ks;
  real ll_corner_x, ll_corner_y, ll_corner_z;
  real tempval;
  int offset = field.get_offset(2);

  yakl::parallel_for("SetTwoForm", topology.n_cells , YAKL_LAMBDA (int iGlob) {
    int k, j, i;
    yakl::unpackIndices(iGlob, topology.n_cells_z, topology.n_cells_y, topology.n_cells_x ,k, j, i);
    tempval = 0.0;

    ll_corner_x = i*dx + xc - Lx/2.;
    ll_corner_y = j*dy + yc - Ly/2.;
    ll_corner_z = k*dz + zc - Lz/2.;

    if (ndim == 2) {
        ll_corner_z = 0.0;
    }

        // xy surface
        for (int nqx=0; nqx<ic_quad_pts_x; nqx++) {
          for (int nqy=0; nqy<ic_quad_pts_y; nqy++) {
            tempval = tempval + initial_value_function(ll_corner_x + x_quad_pts(nqx)*dx, ll_corner_y + y_quad_pts(nqy)*dy, ll_corner_z) * x_quad_wts(nqx) * y_quad_wts(nqy);
        }}
        field.data(ndof+offset, k+ks, j+js, i+is) = tempval;

        if (ndim ==3) {
        // xz surface
        for (int nqx=0; nqx<ic_quad_pts_x; nqx++) {
          for (int nqy=0; nqy<ic_quad_pts_y; nqy++) {
            tempval = tempval + initial_value_function(ll_corner_x + x_quad_pts(nqx)*dx, ll_corner_y, ll_corner_z + z_quad_pts(nqz)*dz) * x_quad_wts(nqx) * z_quad_wts(nqz);
        }}
        field.data(ndof+field.ndofs2+offset, k+ks, j+js, i+is) = tempval;

        // yz surface
        for (int nqy=0; nqy<ic_quad_pts_y; nqy++) {
          for (int nqz=0; nqz<ic_quad_pts_z; nqz++) {
          tempval = tempval + initial_value_function(ll_corner_x, ll_corner_y + y_quad_pts(nqy)*dy, ll_corner_z + z_quad_pts(nqz)*dz) * z_quad_wts(nqz) * y_quad_wts(nqy);
        }}
        field.data(ndof+2*field.ndofs2+offset, k+ks, j+js, i+is) = tempval;
      }

  });
}
