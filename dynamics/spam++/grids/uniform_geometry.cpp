




#include "uniform_geometry.h"

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

  if (ic_quad_pts_x == 3)
  {
    x_quad_pts(0) = 0.112701665379258311482073460022;
    x_quad_pts(1) = 0.500000000000000000000000000000;
    x_quad_pts(2) = 0.887298334620741688517926539980;

    x_quad_wts(0) = 0.277777777777777777777777777779;
    x_quad_wts(1) = 0.444444444444444444444444444444;
    x_quad_wts(2) = 0.277777777777777777777777777779;
  }

  if (ic_quad_pts_y == 3)
  {
    y_quad_pts(0) = 0.112701665379258311482073460022;
    y_quad_pts(1) = 0.500000000000000000000000000000;
    y_quad_pts(2) = 0.887298334620741688517926539980;

    y_quad_wts(0) = 0.277777777777777777777777777779;
    y_quad_wts(1) = 0.444444444444444444444444444444;
    y_quad_wts(2) = 0.277777777777777777777777777779;
  }

  if (ic_quad_pts_z == 3)
  {
    z_quad_pts(0) = 0.112701665379258311482073460022;
    z_quad_pts(1) = 0.500000000000000000000000000000;
    z_quad_pts(2) = 0.887298334620741688517926539980;

    z_quad_wts(0) = 0.277777777777777777777777777779;
    z_quad_wts(1) = 0.444444444444444444444444444444;
    z_quad_wts(2) = 0.277777777777777777777777777779;
  }

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

  yakl::parallel_for( "SetZeroForm" , topology.n_cells , YAKL_LAMBDA (int iGlob) {
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

  yakl::parallel_for( "SetThreeForm" , topology.n_cells , YAKL_LAMBDA (int iGlob) {
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





}




// ALL BROKEN IN MULTIPLE DIMENSIONS- THERE ARE MULTIPLE EDGES TO CONSIDER HERE
template <class F> YAKL_INLINE void UniformRectangularGeometry::set_one_form_values(F const &initial_value_function, Field &field, int ndof) {
  int is = topology.cell_is;
  int js = topology.cell_js;
  int ks = topology.cell_ks;
  real x, y, z, wq;
  real tempval;
  int offset = field.get_offset(1);
  // FIX THIS TO REFER TO PRIMAL CELLS, ETC.

  yakl::parallel_for( "SetOneForm" , topology.ncells , YAKL_LAMBDA (int iGlob) {
    int k, j, i;
    yakl::unpackIndices( iGlob , topology.ncells_z,topology.ncells_y,topology.ncells_x ,k,j,i);
    tempval = 0.0;
    // BROKEN- SHOULD DO EACH DIRECTION INDEP I THINK...
    for (int nq=0; nq<ic_quad_pts_x; nq++) {
      geom.get_line_coords_weights(nq, i, j, k, &x, &y, &z, &wq);
      // SOMEWHAT BROKEN- DOESNT TAKE INTO ACCOUNT VECTOR QUANTITY AND NORMAL/TANGENTIAL DIRECTION FOR 2D and 3D MESHES..
      // ALSO IN 2D WE CAN INTEGRATE ALONG LINE OR TANGENTIAL TO IT! IE PRIMAL VS. DUAL...
      if (ndim >= 2) {
        tempval = tempval + initial_value_function(x,y,z)*wq;
      }
      else {
        tempval = tempval + initial_value_function(x,y,z)*wq;
      }
    }
    field.data(ndof+offset, k+ks, j+js, i+is) = tempval;
  });
}

// ALL BROKEN IN MULTIPLE DIMENSIONS- THERE ARE MULTIPLE FACES TO CONSIDER HERE
template <class F> YAKL_INLINE void set_two_form_values(F const &initial_value_function, Field &field, int ndof) {
  int is = topology.cell_is;
  int js = topology.cell_js;
  int ks = topology.cell_ks;
  real x, y, z, wq;
  real tempval;
  int offset = field.get_offset(2);

  // FIX THIS TO REFER TO PRIMAL CELLS, ETC.
  yakl::parallel_for( "SetTwoForm" , topology.ncells , YAKL_LAMBDA (int iGlob) {
    int k, j, i;
    yakl::unpackIndices( iGlob , topology.ncells_z,topology.ncells_y,topology.ncells_x ,k,j,i);
    tempval = 0.0;
    // BROKEN- SHOULD DO EACH DIRECTION INDEP I THINK...
    for (int nq=0; nq<ic_quad_pts_x*ic_quad_pts_y; nq++) {
        geom.get_surface_coords_weights(nq, i, j, k, &x, &y, &z, &wq);

        // SOMEWHAT BROKEN- DOESNT TAKE INTO ACCOUNT VECTOR QUANTITY AND NORMAL/TANGENTIAL DIRECTION FOR 3D MESHES..
        if (ndims == 3) {
          tempval = tempval + initial_value_function(x,y,z)*wq;
        }
        else {
          tempval = tempval + initial_value_function(x,y,z)*wq;
        }
      }
      field.data(ndof+offset, k+ks, j+js, i+is) = tempval;
  });
}
