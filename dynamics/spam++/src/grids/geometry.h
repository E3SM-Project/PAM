
#ifndef _UNIFORM_GEOMETRY_H_
#define _UNIFORM_GEOMETRY_H_


#include "common.h"
#include "fields.h"
#include "topology.h"


// REALLY HERE WE SHOULD PROBABLY DO SPECIALIZATION BY DIMENSION FOR EVERYTHING....

template<int ndims> struct coords {};

template<> struct coords<1>
{
  real x;
  real y;
  real z;
};

template<> struct coords<2>
{
  real x;
  real y;
  real z;
};

template<> struct coords<3>
{
  real x;
  real y;
  real z;
};


template<uint ndims, uint nquadx, uint nquady, uint nquadz> class Geometry {

public:

  const Topology<ndims> *topology;

  bool is_initialized;
  Geometry();
  Geometry( const Geometry<ndims,nquadx,nquady,nquadz> &geom) = delete;
  Geometry& operator=( const Geometry<ndims,nquadx,nquady,nquadz> &geom) = delete;

  SArray<real,nquadx> x_quad_pts_ref;
  SArray<real,nquadx> x_quad_wts_ref;
  SArray<real,nquady> y_quad_pts_ref;
  SArray<real,nquady> y_quad_wts_ref;
  SArray<real,nquadz> z_quad_pts_ref;
  SArray<real,nquadz> z_quad_wts_ref;

  void YAKL_INLINE get_point_quad_pts_wts(int i, int j, int k, SArray<coords<ndims>,1> &quad_pts_phys, SArray<real,1> &quad_wts_phys) {};

  void YAKL_INLINE get_volume_quad_pts_wts(int i, int j, int k, SArray<coords<ndims>,nquadx,nquady,nquadz> &quad_pts_phys, SArray<real,nquadx,nquady,nquadz> &quad_wts_phys) {};

  void YAKL_INLINE get_edge_quad_pts_wts(int i, int j, int k,
    SArray<coords<ndims>,nquadx> &x_quad_pts_phys, SArray<real,nquadx> &x_quad_wts_phys,
    SArray<coords<ndims>,nquady> &y_quad_pts_phys, SArray<real,nquady> &y_quad_wts_phys,
    SArray<coords<ndims>,nquadz> &z_quad_pts_phys, SArray<real,nquadz> &z_quad_wts_phys) {};

  void YAKL_INLINE get_edge_normals(int i, int j, int k, SArray<coords<ndims>,nquadx> &x_normals, SArray<coords<ndims>,nquady> &y_normals, SArray<coords<ndims>,nquadz> &z_normals) {};
  void YAKL_INLINE get_edge_tangents(int i, int j, int k, SArray<coords<ndims>,nquadx> &x_tangents, SArray<coords<ndims>,nquady> &y_tangents, SArray<coords<ndims>,nquadz> &z_tangents) {};
  void YAKL_INLINE get_surface_normals(int i, int j, int k, SArray<coords<ndims>,nquadx,nquady> &xy_normals, SArray<coords<ndims>,nquady,nquadz> &yz_normals, SArray<coords<ndims>,nquadx,nquadz> &xz_normals) {};


  void YAKL_INLINE get_surface_quad_pts_wts(int i, int j, int k,
      SArray<coords<ndims>,nquadx,nquady> &xy_quad_pts_phys, SArray<real,nquadx,nquady> &xy_quad_wts_phys,
      SArray<coords<ndims>,nquady,nquadz> &yz_quad_pts_phys, SArray<real,nquady,nquadz> &yz_quad_wts_phys,
      SArray<coords<ndims>,nquadx,nquadz> &xz_quad_pts_phys, SArray<real,nquadx,nquadz> &xz_quad_wts_phys) {};

  void initialize(const Topology<ndims> &topo);

  YAKL_INLINE void set_point_values(real (*initial_value_function)(real, real, real), Field<ndims> &field, int ndof);
  YAKL_INLINE void set_volume_values(real (*initial_value_function)(real, real, real), Field<ndims> &field, int ndof);


  // THESE MUST CHANGE FOR EDGE/SURFACES...MAYBE RETURN A COORD STRUCT?
  YAKL_INLINE void set_edge_values(real (*initial_value_function)(real, real, real), Field<ndims> &field, int ndof);
  YAKL_INLINE void set_surface_values(real (*initial_value_function)(real, real, real), Field<ndims> &field, int ndof);


};




template<uint ndims, uint nquadx, uint nquady, uint nquadz> class UniformRectangularGeometry: public Geometry<ndims,nquadx,nquady,nquadz> {

public:

  real dx, dy, dz;
  real Lx, Ly, Lz;
  real xc, yc, zc;
  void initialize(const Topology<ndims> &topo, real dxx, real dyy, real dzz, real xlen, real ylen, real zlen, real xcent, real ycent, real zcent);
  void printinfo();
  void YAKL_INLINE get_point_quad_pts_wts(int i, int j, int k, SArray<coords<ndims>,1> &quad_pts_phys, SArray<real,1> &quad_wts_phys);
  void YAKL_INLINE get_volume_quad_pts_wts(int i, int j, int k, SArray<coords<ndims>,nquadx,nquady,nquadz> &quad_pts_phys, SArray<real,nquadx,nquady,nquadz> &quad_wts_phys);
  void  YAKL_INLINE get_edge_quad_pts_wts(int i, int j, int k,
    SArray<coords<ndims>,nquadx> &x_quad_pts_phys, SArray<real,nquadx> &x_quad_wts_phys,
    SArray<coords<ndims>,nquady> &y_quad_pts_phys, SArray<real,nquady> &y_quad_wts_phys,
    SArray<coords<ndims>,nquadz> &z_quad_pts_phys, SArray<real,nquadz> &z_quad_wts_phys);
  void YAKL_INLINE get_surface_quad_pts_wts(int i, int j, int k,
    SArray<coords<ndims>,nquadx,nquady> &xy_quad_pts_phys, SArray<real,nquadx,nquady> &xy_quad_wts_phys,
    SArray<coords<ndims>,nquady,nquadz> &yz_quad_pts_phys, SArray<real,nquady,nquadz> &yz_quad_wts_phys,
    SArray<coords<ndims>,nquadx,nquadz> &xz_quad_pts_phys, SArray<real,nquadx,nquadz> &xz_quad_wts_phys);
};



// EVENTUALLY DO TEMPLATE PARTIAL SPECIALIZATION HERE
// WITH A NOT IMPLEMENTED ERROR FOR THE DEFAULT
template <class T, uint npts> void set_ref_quad_pts_wts(SArray<T,npts> &pts, SArray<T,npts> &wts) {

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

    template<uint ndims, uint nquadx, uint nquady, uint nquadz> Geometry<ndims,nquadx,nquady,nquadz>::Geometry()
    {
      this->is_initialized = false;
      std::cout << "CREATED GEOMETRY\n";
    }

  template<uint ndims, uint nquadx, uint nquady, uint nquadz> void Geometry<ndims,nquadx,nquady,nquadz>::initialize(const Topology<ndims> &topo)
  {

  this->topology = &topo;

  set_ref_quad_pts_wts<real,nquadx>(this->x_quad_pts_ref, this->x_quad_wts_ref);
  set_ref_quad_pts_wts<real,nquady>(this->y_quad_pts_ref, this->y_quad_wts_ref);
  set_ref_quad_pts_wts<real,nquadz>(this->z_quad_pts_ref, this->z_quad_wts_ref);

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

// x edge is the edge where x varies and y,z are constant
// xy surface is the surface where x and y vary and z is constant
// etc.

// in 2D edges are stored y,x  = (U,V)
// in 3D edges are stored z,y,x
// in 3D surfaces are stored yz, xz, xy = (U,V,W)


template<uint ndims, uint nquadx, uint nquady, uint nquadz> YAKL_INLINE void Geometry<ndims,nquadx,nquady,nquadz>::set_point_values(real (*initial_value_function)(real, real, real), Field<ndims> &field, int ndof)
{

  int is = this->topology->is;
  int js = this->topology->js;
  int ks = this->topology->ks;

  SArray<coords<ndims>,1> quad_pts_phys;
  SArray<real,1> quad_wts_phys;

  yakl::parallel_for("SetPointValues", this->topology->n_cells, YAKL_LAMBDA (int iGlob) {
    int k, j, i;
    yakl::unpackIndices(iGlob, this->topology->n_cells_z, this->topology->n_cells_y, this->topology->n_cells_x ,k, j, i);

      get_point_quad_pts_wts(i, j, k, quad_pts_phys, quad_wts_phys);
      field.data(ndof, k+ks, j+js, i+is) = initial_value_function(quad_pts_phys(0).x, quad_pts_phys(0).y, quad_pts_phys(0).z) * quad_wts_phys(0);
  });
}


template<uint ndims, uint nquadx, uint nquady, uint nquadz> YAKL_INLINE void Geometry<ndims,nquadx,nquady,nquadz>::set_volume_values(real (*initial_value_function)(real, real, real), Field<ndims> &field, int ndof)
{

  SArray<coords<ndims>,nquadx,nquady,nquadz> quad_pts_phys;
  SArray<real,nquadx,nquady,nquadz> quad_wts_phys;

  int is = this->topology->is;
  int js = this->topology->js;
  int ks = this->topology->ks;
  real tempval;
  int offset = field.get_offset(3);

  yakl::parallel_for("SetVolumeValues", this->topology->n_cells, YAKL_LAMBDA (int iGlob) {
    int k, j, i;
    yakl::unpackIndices(iGlob, this->topology->n_cells_z, this->topology->n_cells_y, this->topology->n_cells_x ,k, j, i);

    get_volume_quad_pts_wts(i, j, k, quad_pts_phys, quad_wts_phys);
    tempval = 0.0;
    for (int nqx=0; nqx<nquadx; nqx++) {
      for (int nqy=0; nqy<nquady; nqy++) {
        for (int nqz=0; nqz<nquadz; nqz++) {
          tempval = tempval + initial_value_function(quad_pts_phys(nqx, nqy, nqz).x, quad_pts_phys(nqx, nqy, nqz).y, quad_pts_phys(nqx, nqy, nqz).z) * quad_wts_phys(nqx, nqy, nqz);
          } } }
      field.data(ndof+offset, k+ks, j+js, i+is) = tempval;
  });
}


// THESE ARE ACTUALLY BROKEN
// REALLY SHOULD BE SPECIALIZED FOR DIFFERENT DIMENSIONS ie 1-forms in 1D, 2D, 3D all different
// 2-forms in 2D/3D are also different

template<uint ndims, uint nquadx, uint nquady, uint nquadz> YAKL_INLINE void Geometry<ndims,nquadx,nquady,nquadz>::set_edge_values(real (*initial_value_function)(real, real, real), Field<ndims> &field, int ndof)
{

  SArray<coords<ndims>,nquadx> x_quad_pts_phys;
  SArray<real,nquadx> x_quad_wts_phys;
  SArray<coords<ndims>,nquady> y_quad_pts_phys;
  SArray<real,nquady> y_quad_wts_phys;
  SArray<coords<ndims>,nquadz> z_quad_pts_phys;
  SArray<real,nquadz> z_quad_wts_phys;

  int is = this->topology->is;
  int js = this->topology->js;
  int ks = this->topology->ks;
  real tempval;
  int offset = field.get_offset(1);

  yakl::parallel_for("SetEdgeValues", this->topology->n_cells , YAKL_LAMBDA (int iGlob) {
    int k, j, i;
    yakl::unpackIndices(iGlob, this->topology->n_cells_z, this->topology->n_cells_y, this->topology->n_cells_x ,k, j, i);


    get_edge_quad_pts_wts(i, j, k, x_quad_pts_phys, x_quad_wts_phys, y_quad_pts_phys, y_quad_wts_phys, z_quad_pts_phys, z_quad_wts_phys);

    if (ndims ==3) {
    // z-edge
    tempval = 0.0;
    for (int nqz=0; nqz<nquadz; nqz++) {
      tempval = tempval + initial_value_function(z_quad_pts_phys(nqz).x, z_quad_pts_phys(nqz).y, z_quad_pts_phys(nqz).z) * z_quad_wts_phys(nqz);
    }
    field.data(ndof+offset, k+ks, j+js, i+is) = tempval;
  }

  // y-edge
  if (ndims >=2) {
    tempval = 0.0;
  for (int nqy=0; nqy<nquady; nqy++) {
    tempval = tempval + initial_value_function(y_quad_pts_phys(nqy).x, y_quad_pts_phys(nqy).y, y_quad_pts_phys(nqy).z) * y_quad_wts_phys(nqy);
  }
  field.data(ndof+(ndims-2)*field.ndof1+offset, k+ks, j+js, i+is) = tempval;
  }

        // x-edge
        tempval = 0.0;
        for (int nqx=0; nqx<nquadx; nqx++) {
          tempval = tempval + initial_value_function(x_quad_pts_phys(nqx).x, x_quad_pts_phys(nqx).y, x_quad_pts_phys(nqx).z) * x_quad_wts_phys(nqx);
          //std::cout << i << " " << j << " " << k << " " << nqx << " " << x_quad_pts_phys(0,nqx) << " " << x_quad_pts_phys(1,nqx) << " " << x_quad_pts_phys(2,nqx) << " " << tempval << "\n";
        }
        field.data(ndof+(ndims-1)*field.ndof1+offset, k+ks, j+js, i+is) = tempval;




  });
}



template<uint ndims, uint nquadx, uint nquady, uint nquadz> YAKL_INLINE void Geometry<ndims,nquadx,nquady,nquadz>::set_surface_values(real (*initial_value_function)(real, real, real), Field<ndims> &field, int ndof)
{

  SArray<coords<ndims>,nquadx,nquady> xy_quad_pts_phys;
  SArray<real,nquadx,nquady> xy_quad_wts_phys;
  SArray<coords<ndims>,nquady,nquadz> yz_quad_pts_phys;
  SArray<real,nquady,nquadz> yz_quad_wts_phys;
  SArray<coords<ndims>,nquadx,nquadz> xz_quad_pts_phys;
  SArray<real,nquadx,nquadz> xz_quad_wts_phys;

  int is = this->topology->is;
  int js = this->topology->js;
  int ks = this->topology->ks;
  real tempval;
  int offset = field.get_offset(2);

  yakl::parallel_for("SetSurfaceValues", this->topology->n_cells , YAKL_LAMBDA (int iGlob) {
    int k, j, i;
    yakl::unpackIndices(iGlob, this->topology->n_cells_z, this->topology->n_cells_y, this->topology->n_cells_x ,k, j, i);

    get_surface_quad_pts_wts(i, j, k, xy_quad_pts_phys, xy_quad_wts_phys, yz_quad_pts_phys, yz_quad_wts_phys, xz_quad_pts_phys, xz_quad_wts_phys);

    if (ndims ==3) {
    // yz surface
    tempval = 0.0;
    for (int nqy=0; nqy<nquady; nqy++) {
      for (int nqz=0; nqz<nquadz; nqz++) {
      tempval = tempval + initial_value_function(yz_quad_pts_phys(nqy,nqz).x, yz_quad_pts_phys(nqy,nqz).y, yz_quad_pts_phys(nqy,nqz).z) * yz_quad_wts_phys(nqy,nqz);
    }}
    field.data(ndof+offset, k+ks, j+js, i+is) = tempval;

        // xz surface
        tempval = 0.0;
        for (int nqx=0; nqx<nquadx; nqx++) {
          for (int nqz=0; nqz<nquadz; nqz++) {
            tempval = tempval + initial_value_function(xz_quad_pts_phys(nqx,nqz).x, xz_quad_pts_phys(nqx,nqz).y, xz_quad_pts_phys(nqx,nqz).z) * xz_quad_wts_phys(nqx,nqz);
        }}
        field.data(ndof+field.ndof2+offset, k+ks, j+js, i+is) = tempval;
      }

      // xy surface
      tempval = 0.0;
      for (int nqx=0; nqx<nquadx; nqx++) {
        for (int nqy=0; nqy<nquady; nqy++) {
          tempval = tempval + initial_value_function(xy_quad_pts_phys(nqx,nqy).x, xy_quad_pts_phys(nqx,nqy).y, xy_quad_pts_phys(nqx,nqy).z) * xy_quad_wts_phys(nqx,nqy);
      }}
      if (ndims ==3) { field.data(ndof+2*field.ndof2+offset, k+ks, j+js, i+is) = tempval; }
      if (ndims ==2) { field.data(ndof+offset, k+ks, j+js, i+is) = tempval; }

  });
}












  template<uint ndims, uint nquadx, uint nquady, uint nquadz> void UniformRectangularGeometry<ndims,nquadx,nquady,nquadz>::initialize(const Topology<ndims> &topo, real dxx, real dyy, real dzz, real xlen, real ylen, real zlen, real xcent, real ycent, real zcent)
  {

  Geometry<ndims,nquadx,nquady,nquadz>::initialize(topo);

  this->Lx = xlen;
  this->Ly = ylen;
  this->Lz = zlen;

  this->xc = xcent;
  this->yc = ycent;
  this->zc = zcent;

  this->dx = dxx;
  this->dy = dyy;
  this->dz = dzz;

  if (ndims == 1)
  {
    this->Ly = 0.0;
    this->Lz = 0.0;
    this->yc = 0.0;
    this->zc = 0.0;
    this->dy = 0.0;
    this->dz = 0.0;
  }

  if (ndims == 2)
  {
    this->Lz = 0.0;
    this->zc = 0.0;
    this->dz = 0.0;
  }

  this->is_initialized = true;


}

template<uint ndims, uint nquadx, uint nquady, uint nquadz> void UniformRectangularGeometry<ndims,nquadx,nquady,nquadz>::printinfo()
{
  std::cout << "uniform rectangular geometry info\n" << std::flush;
  std::cout << "Lx " << this->Lx << " Ly " << this->Ly << " Lz " << this->Lz << "\n" << std::flush;
  std::cout << "xc " << this->xc << " yc " << this->yc << " zc " << this->zc << "\n" << std::flush;
  std::cout << "dx " << this->dx << " dy " << this->dy << " dz " << this->dz << "\n" << std::flush;
}


template<uint ndims, uint nquadx, uint nquady, uint nquadz> void YAKL_INLINE UniformRectangularGeometry<ndims,nquadx,nquady,nquadz>::get_point_quad_pts_wts(int i, int j, int k, SArray<coords<ndims>,1> &quad_pts_phys, SArray<real,1> &quad_wts_phys)
{
  quad_pts_phys(0).x = i*this->dx + this->xc - this->Lx/2.;
  quad_pts_phys(0).y = j*this->dy + this->yc - this->Ly/2.;
  quad_pts_phys(0).z = k*this->dz + this->zc - this->Lz/2.;
  quad_wts_phys(0) = 1.;
}


template<uint ndims, uint nquadx, uint nquady, uint nquadz> void YAKL_INLINE UniformRectangularGeometry<ndims,nquadx,nquady,nquadz>::get_volume_quad_pts_wts(int i, int j, int k, SArray<coords<ndims>,nquadx,nquady,nquadz> &quad_pts_phys, SArray<real,nquadx,nquady,nquadz> &quad_wts_phys)
{
  real ll_corner_x = i*this->dx + this->xc - this->Lx/2.;
  real ll_corner_y = j*this->dy + this->yc - this->Ly/2.;
  real ll_corner_z = k*this->dz + this->zc - this->Lz/2.;
  for (int nqx=0; nqx<nquadx; nqx++) {
    for (int nqy=0; nqy<nquady; nqy++) {
      for (int nqz=0; nqz<nquadz; nqz++) {
        quad_pts_phys(nqx, nqy, nqz).x = ll_corner_x + this->x_quad_pts_ref(nqx) * this->dx;
        quad_pts_phys(nqx, nqy, nqz).y = ll_corner_y + this->y_quad_pts_ref(nqy) * this->dy;
        quad_pts_phys(nqx, nqy, nqz).z = ll_corner_z + this->z_quad_pts_ref(nqz) * this->dz;
        quad_wts_phys(nqx, nqy, nqz) = this->x_quad_wts_ref(nqx) * this->y_quad_wts_ref(nqy) * this->z_quad_wts_ref(nqz) * this->dx * this->dy * this->dz;
      }}}

}

template<uint ndims, uint nquadx, uint nquady, uint nquadz> void YAKL_INLINE UniformRectangularGeometry<ndims,nquadx,nquady,nquadz>::get_edge_quad_pts_wts(int i, int j, int k,
  SArray<coords<ndims>,nquadx> &x_quad_pts_phys, SArray<real,nquadx> &x_quad_wts_phys,
  SArray<coords<ndims>,nquady> &y_quad_pts_phys, SArray<real,nquady> &y_quad_wts_phys,
  SArray<coords<ndims>,nquadz> &z_quad_pts_phys, SArray<real,nquadz> &z_quad_wts_phys)
{
  real ll_corner_x = i*this->dx + this->xc - this->Lx/2.;
  real ll_corner_y = j*this->dy + this->yc - this->Ly/2.;
  real ll_corner_z = k*this->dz + this->zc - this->Lz/2.;

  for (int nqx=0; nqx<nquadx; nqx++) {
      x_quad_pts_phys(nqx).x = ll_corner_x + this->x_quad_pts_ref(nqx) * this->dx;
      x_quad_pts_phys(nqx).y = ll_corner_y;
      x_quad_pts_phys(nqx).z = ll_corner_z;
      x_quad_wts_phys(nqx) = this->x_quad_wts_ref(nqx) * this->dx;
  }

  for (int nqy=0; nqy<nquady; nqy++) {
      y_quad_pts_phys(nqy).x = ll_corner_x;
      y_quad_pts_phys(nqy).y = ll_corner_y + this->y_quad_pts_ref(nqy) * this->dy;
      y_quad_pts_phys(nqy).z = ll_corner_z;
      y_quad_wts_phys(nqy) = this->y_quad_wts_ref(nqy) * this->dy;
  }

  for (int nqz=0; nqz<nquadz; nqz++) {
      z_quad_pts_phys(nqz).x = ll_corner_x;
      z_quad_pts_phys(nqz).y = ll_corner_y;
      z_quad_pts_phys(nqz).z = ll_corner_z + this->z_quad_pts_ref(nqz) * this->dz;
      z_quad_wts_phys(nqz) = this->z_quad_wts_ref(nqz) * this->dz;
  }

}

template<uint ndims, uint nquadx, uint nquady, uint nquadz> void YAKL_INLINE UniformRectangularGeometry<ndims,nquadx,nquady,nquadz>::get_surface_quad_pts_wts(int i, int j, int k,
  SArray<coords<ndims>,nquadx,nquady> &xy_quad_pts_phys, SArray<real,nquadx,nquady> &xy_quad_wts_phys,
  SArray<coords<ndims>,nquady,nquadz> &yz_quad_pts_phys, SArray<real,nquady,nquadz> &yz_quad_wts_phys,
  SArray<coords<ndims>,nquadx,nquadz> &xz_quad_pts_phys, SArray<real,nquadx,nquadz> &xz_quad_wts_phys)
{
  real ll_corner_x = i*this->dx + this->xc - this->Lx/2.;
  real ll_corner_y = j*this->dy + this->yc - this->Ly/2.;
  real ll_corner_z = k*this->dz + this->zc - this->Lz/2.;

  for (int nqx=0; nqx<nquadx; nqx++) {
    for (int nqy=0; nqy<nquady; nqy++) {
        xy_quad_pts_phys(nqx, nqy).x = ll_corner_x + this->x_quad_pts_ref(nqx) * this->dx;
        xy_quad_pts_phys(nqx, nqy).y = ll_corner_y + this->y_quad_pts_ref(nqy) * this->dy;
        xy_quad_pts_phys(nqx, nqy).z = ll_corner_z;
        xy_quad_wts_phys(nqx, nqy) = this->x_quad_wts_ref(nqx) * this->y_quad_wts_ref(nqy) * this->dx * this->dy;
  }}

    for (int nqy=0; nqy<nquady; nqy++) {
      for (int nqz=0; nqz<nquadz; nqz++) {
        yz_quad_pts_phys(nqy, nqz).x = ll_corner_x;
        yz_quad_pts_phys(nqy, nqz).y = ll_corner_y + this->y_quad_pts_ref(nqy) * this->dy;
        yz_quad_pts_phys(nqy, nqz).z = ll_corner_z + this->z_quad_pts_ref(nqz) * this->dz;
        yz_quad_wts_phys(nqy, nqz) = this->y_quad_wts_ref(nqy) * this->z_quad_wts_ref(nqz) * this->dy * this->dz;
  }}

  for (int nqx=0; nqx<nquadx; nqx++) {
      for (int nqz=0; nqz<nquadz; nqz++) {
        xz_quad_pts_phys(nqx, nqz).x = ll_corner_x + this->x_quad_pts_ref(nqx) * this->dx;
        xz_quad_pts_phys(nqx, nqz).y = ll_corner_y;
        xz_quad_pts_phys(nqx, nqz).z = ll_corner_z + this->z_quad_pts_ref(nqz) * this->dz;
        xz_quad_wts_phys(nqx, nqz) = this->x_quad_wts_ref(nqx) * this->z_quad_wts_ref(nqz) * this->dx * this->dz;
  }}

}









#endif
