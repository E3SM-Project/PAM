
#ifndef _UNIFORM_GEOMETRY_H_
#define _UNIFORM_GEOMETRY_H_


#include "common.h"
#include "fields.h"
#include "topology.h"
#include "model.h"


template<int number_of_dims> struct coords {};
template<int number_of_dims> struct vec {};
template<uint number_of_dims, uint nquadx, uint nquady, uint nquadz> class Geometry {};
template<uint number_of_dims, uint nquadx, uint nquady, uint nquadz> class UniformRectangularTwistedGeometry: public Geometry<ndims,nquadx,nquady,nquadz> {};
template<uint number_of_dims, uint nquadx, uint nquady, uint nquadz> class UniformRectangularStraightGeometry: public Geometry<ndims,nquadx,nquady,nquadz> {};



template <class T, uint npts> void set_ref_quad_pts_wts(SArray<T,npts> &pts, SArray<T,npts> &wts)
{
if (npts == 1)
{
    pts(0) = 0.5;
    wts(0) = 1.0;
}

if (npts == 2)
{
    pts(0) = -1./(2.*sqrt(3.)) + 0.5;
    pts(1) = 1./(2.*sqrt(3.)) + 0.5;

    wts(0) = 0.5;
    wts(1) = 0.5;
}

if (npts == 3)
{
    pts(0) = -1./2.*sqrt(3./5.) + 0.5;
    pts(1) = 0.5;
    pts(2) = 1./2.*sqrt(3./5.) + 0.5;

    wts(0) = 5./18.;
    wts(1) = 4./9.;
    wts(2) = 5./18.;
}
//   {
//     pts(0) = 0.112701665379258311482073460022;
//     pts(1) = 0.500000000000000000000000000000;
//     pts(2) = 0.887298334620741688517926539980;
//
//     wts(0) = 0.277777777777777777777777777779;
//     wts(1) = 0.444444444444444444444444444444;
//     wts(2) = 0.277777777777777777777777777779;
//   }

if (npts == 4)
{
    pts(0) = -sqrt(3./7. + 2./7. * sqrt(6./5.)) + 0.5;
    pts(1) = -sqrt(3./7. - 2./7. * sqrt(6./5.)) + 0.5;
    pts(2) = sqrt(3./7. - 2./7. * sqrt(6./5.)) + 0.5;
    pts(3) = sqrt(3./7. + 2./7. * sqrt(6./5.)) + 0.5;

    wts(0) = (18. - sqrt(30.))/72.;
    wts(1) = (18. + sqrt(30.))/72.;
    wts(2) = (18. + sqrt(30.))/72.;
    wts(3) = (18. - sqrt(30.))/72.;
}

if (npts == 5)
{
    pts(0) = -1./3.*sqrt(5. + 2. * sqrt(10./7.)) + 0.5;
    pts(1) = -1./3.*sqrt(5. - 2. * sqrt(10./7.)) + 0.5;
    pts(2) = 0.5;
    pts(3) = 1./3.*sqrt(5. - 2. * sqrt(10./7.)) + 0.5;
    pts(4) = 1./3.*sqrt(5. + 2. * sqrt(10./7.)) + 0.5;

    wts(0) = (322. - 13.*sqrt(70.))/1800.;
    wts(1) = (322. + 13.*sqrt(70.))/1800.;
    wts(2) = 64./225.;
    wts(3) = (322. + 13.*sqrt(70.))/1800.;
    wts(4) = (322. - 13.*sqrt(70.))/1800.;
}

}






// *********************   1D  ******************************/

template<> struct coords<1>
{
  real x=0.;
};

template<uint nquadx, uint nquady, uint nquadz> class Geometry<1,nquadx,nquady,nquadz> {

public:

  const Topology *topology;
  bool is_initialized;
  Geometry();
  Geometry( const Geometry<1,nquadx,nquady,nquadz> &geom) = delete;
  Geometry& operator=( const Geometry<1,nquadx,nquady,nquadz> &geom) = delete;

  SArray<real,nquadx> x_quad_pts_ref;
  SArray<real,nquadx> x_quad_wts_ref;

  void initialize(const Topology &topo);

  virtual void YAKL_INLINE get_0form_quad_pts_wts(int i, SArray<coords<1>,1> &quad_pts_phys, SArray<real,1> &quad_wts_phys) {};
  virtual void YAKL_INLINE get_1form_quad_pts_wts(int i, SArray<coords<1>,nquadx> &quad_pts_phys, SArray<real,nquadx> &quad_wts_phys) {};
  //virtual real YAKL_INLINE get_area_0form(int k, int j, int i) {};
  //virtual real YAKL_INLINE get_area_1form(int k, int j, int i) {};
  virtual real YAKL_INLINE get_area_lform(int l, int d, int k, int j, int i) {};

  YAKL_INLINE void set_0form_values(real (*initial_value_function)(real), Field &field, int ndof);
  YAKL_INLINE void set_1form_values(real (*initial_value_function)(real), Field &field, int ndof);

};

template<uint nquadx, uint nquady, uint nquadz> class UniformRectangularStraightGeometry<1,nquadx,nquady,nquadz>: public Geometry<1,nquadx,nquady,nquadz> {
public:

  real dx, dy, dz;
  real Lx, Ly, Lz;
  real xc, yc, zc;
  void initialize(const Topology &topo, const ModelParameters &params);

  void printinfo();

  void YAKL_INLINE get_0form_quad_pts_wts(int i, SArray<coords<1>,1> &quad_pts_phys, SArray<real,1> &quad_wts_phys);
  void YAKL_INLINE get_1form_quad_pts_wts(int i, SArray<coords<1>,nquadx> &quad_pts_phys, SArray<real,nquadx> &quad_wts_phys);
  //real YAKL_INLINE get_area_0form(int k, int j, int i);
  //real YAKL_INLINE get_area_1form(int k, int j, int i);
  real YAKL_INLINE get_area_lform(int l, int d, int k, int j, int i);

};

template<uint nquadx, uint nquady, uint nquadz> class UniformRectangularTwistedGeometry<1,nquadx,nquady,nquadz>: public Geometry<1,nquadx,nquady,nquadz> {
public:

  real dx, dy, dz;
  real Lx, Ly, Lz;
  real xc, yc, zc;
  void initialize(const Topology &topo, const ModelParameters &params);

  void printinfo();

  void YAKL_INLINE get_0form_quad_pts_wts(int i, SArray<coords<1>,1> &quad_pts_phys, SArray<real,1> &quad_wts_phys);
  void YAKL_INLINE get_1form_quad_pts_wts(int i, SArray<coords<1>,nquadx> &quad_pts_phys, SArray<real,nquadx> &quad_wts_phys);
  //real YAKL_INLINE get_area_0form(int k, int j, int i);
  //real YAKL_INLINE get_area_1form(int k, int j, int i);
  real YAKL_INLINE get_area_lform(int l, int d, int k, int j, int i);

};












template<uint nquadx, uint nquady, uint nquadz> Geometry<1,nquadx,nquady,nquadz>::Geometry()
    {
      this->is_initialized = false;
      std::cout << "CREATED GEOMETRY\n";
    }

  template<uint nquadx, uint nquady, uint nquadz> void Geometry<1,nquadx,nquady,nquadz>::initialize(const Topology &topo)
  {

  this->topology = &topo;
  set_ref_quad_pts_wts<real,nquadx>(this->x_quad_pts_ref, this->x_quad_wts_ref);
}


template<uint nquadx, uint nquady, uint nquadz> YAKL_INLINE void Geometry<1,nquadx,nquady,nquadz>::set_0form_values(real (*initial_value_function)(real), Field &field, int ndof)
{

int is = this->topology->is;
int js = this->topology->js;
int ks = this->topology->ks;

SArray<coords<1>,1> quad_pts_phys;
SArray<real,1> quad_wts_phys;
int offset = field.get_offset(1);

yakl::parallel_for("Set0FormValues", this->topology->n_cells, YAKL_LAMBDA (int iGlob) {
  int k, j, i;
  yakl::unpackIndices(iGlob, this->topology->n_cells_z, this->topology->n_cells_y, this->topology->n_cells_x ,k, j, i);
    get_0form_quad_pts_wts(i, quad_pts_phys, quad_wts_phys);
    field.data(ndof+offset, k+ks, j+js, i+is) = initial_value_function(quad_pts_phys(0).x) * quad_wts_phys(0);
});
}

template<uint nquadx, uint nquady, uint nquadz> YAKL_INLINE void Geometry<1,nquadx,nquady,nquadz>::set_1form_values(real (*initial_value_function)(real), Field &field, int ndof)
{

SArray<coords<1>,nquadx> quad_pts_phys;
SArray<real,nquadx> quad_wts_phys;

int is = this->topology->is;
int js = this->topology->js;
int ks = this->topology->ks;
real tempval;

yakl::parallel_for("Set1FormValues", this->topology->n_cells, YAKL_LAMBDA (int iGlob) {
  int k, j, i;
  yakl::unpackIndices(iGlob, this->topology->n_cells_z, this->topology->n_cells_y, this->topology->n_cells_x ,k, j, i);

  get_1form_quad_pts_wts(i, quad_pts_phys, quad_wts_phys);
  tempval = 0.0;
  for (int nqx=0; nqx<nquadx; nqx++) {
        tempval = tempval + initial_value_function(quad_pts_phys(nqx).x) * quad_wts_phys(nqx);
        }
    field.data(ndof, k+ks, j+js, i+is) = tempval;
});
}



template<uint nquadx, uint nquady, uint nquadz> void UniformRectangularStraightGeometry<1,nquadx,nquady,nquadz>::initialize(const Topology &topo, const ModelParameters &params)
{

Geometry<1,nquadx,nquady,nquadz>::initialize(topo);

this->Lx = params.xlen;
this->Ly = 0.;
this->Lz = 0.;

this->xc = params.xc;
this->yc = 0.;
this->zc = 0.;

this->dx = params.xlen/params.nx_glob;
this->dy = 0.;
this->dz = 0.;

this->is_initialized = true;
}

template<uint nquadx, uint nquady, uint nquadz> void UniformRectangularTwistedGeometry<1,nquadx,nquady,nquadz>::initialize(const Topology &topo, const ModelParameters &params)
{

Geometry<1,nquadx,nquady,nquadz>::initialize(topo);

this->Lx = params.xlen;
this->Ly = 0.;
this->Lz = 0.;

this->xc = params.xc;
this->yc = 0.;
this->zc = 0.;

this->dx = params.xlen/params.nx_glob;
this->dy = 0.;
this->dz = 0.;

this->is_initialized = true;
}

template<uint nquadx, uint nquady, uint nquadz> void UniformRectangularStraightGeometry<1,nquadx,nquady,nquadz>::printinfo()
{
  std::cout << "uniform rectangular geometry 1D info: straight\n" << std::flush;
  std::cout << "Lx " << this->Lx << " Ly " << this->Ly << " Lz " << this->Lz << "\n" << std::flush;
  std::cout << "xc " << this->xc << " yc " << this->yc << " zc " << this->zc << "\n" << std::flush;
  std::cout << "dx " << this->dx << " dy " << this->dy << " dz " << this->dz << "\n" << std::flush;
}

template<uint nquadx, uint nquady, uint nquadz> void UniformRectangularTwistedGeometry<1,nquadx,nquady,nquadz>::printinfo()
{
  std::cout << "uniform rectangular geometry 1D info: twisted\n" << std::flush;
  std::cout << "Lx " << this->Lx << " Ly " << this->Ly << " Lz " << this->Lz << "\n" << std::flush;
  std::cout << "xc " << this->xc << " yc " << this->yc << " zc " << this->zc << "\n" << std::flush;
  std::cout << "dx " << this->dx << " dy " << this->dy << " dz " << this->dz << "\n" << std::flush;
}



template<uint nquadx, uint nquady, uint nquadz> real YAKL_INLINE UniformRectangularStraightGeometry<1,nquadx,nquady,nquadz>::get_area_lform(int l, int d, int k, int j, int i)
{
  if (l == 0) {return 1.;}
  if (l == 1) {return this->dx;}
}

template<uint nquadx, uint nquady, uint nquadz> real YAKL_INLINE UniformRectangularTwistedGeometry<1,nquadx,nquady,nquadz>::get_area_lform(int l, int d, int k, int j, int i)
{
  if (l == 0) {return 1.;}
  if (l == 1) {return this->dx;}
}

// template<uint nquadx, uint nquady, uint nquadz> real YAKL_INLINE UniformRectangularStraightGeometry<1,nquadx,nquady,nquadz>::get_area_1form(int k, int j, int i)
// {
//   return this->dx;
// }
// 
// template<uint nquadx, uint nquady, uint nquadz> real YAKL_INLINE UniformRectangularTwistedGeometry<1,nquadx,nquady,nquadz>::get_area_1form(int k, int j, int i)
// {
//   return this->dx;
// }
// 
// template<uint nquadx, uint nquady, uint nquadz> real YAKL_INLINE UniformRectangularStraightGeometry<1,nquadx,nquady,nquadz>::get_area_0form(int k, int j, int i)
// {
//   return 1.0;
// }
// 
// template<uint nquadx, uint nquady, uint nquadz> real YAKL_INLINE UniformRectangularTwistedGeometry<1,nquadx,nquady,nquadz>::get_area_0form(int k, int j, int i)
// {
//   return 1.0;
// }

template<uint nquadx, uint nquady, uint nquadz> void YAKL_INLINE UniformRectangularTwistedGeometry<1,nquadx,nquady,nquadz>::get_0form_quad_pts_wts(int i, SArray<coords<1>,1> &quad_pts_phys, SArray<real,1> &quad_wts_phys)
{
  quad_pts_phys(0).x = (i+this->topology->i_beg)*this->dx + this->xc - this->Lx/2.;
  quad_wts_phys(0) = 1.;
}

template<uint nquadx, uint nquady, uint nquadz> void YAKL_INLINE UniformRectangularTwistedGeometry<1,nquadx,nquady,nquadz>::get_1form_quad_pts_wts(int i, SArray<coords<1>,nquadx> &quad_pts_phys, SArray<real,nquadx> &quad_wts_phys)
{
  real ll_corner_x = (i+this->topology->i_beg)*this->dx + this->xc - this->Lx/2.;
  for (int nqx=0; nqx<nquadx; nqx++)
  {
    quad_pts_phys(nqx).x = ll_corner_x + this->x_quad_pts_ref(nqx) * this->dx;
    quad_wts_phys(nqx) = this->x_quad_wts_ref(nqx) * this->dx;
  }
}

template<uint nquadx, uint nquady, uint nquadz> void YAKL_INLINE UniformRectangularStraightGeometry<1,nquadx,nquady,nquadz>::get_0form_quad_pts_wts(int i, SArray<coords<1>,1> &quad_pts_phys, SArray<real,1> &quad_wts_phys)
{
  quad_pts_phys(0).x = (i+this->topology->i_beg)*this->dx + this->xc - this->Lx/2. + this->dx/2.;
  quad_wts_phys(0) = 1.;
}

template<uint nquadx, uint nquady, uint nquadz> void YAKL_INLINE UniformRectangularStraightGeometry<1,nquadx,nquady,nquadz>::get_1form_quad_pts_wts(int i, SArray<coords<1>,nquadx> &quad_pts_phys, SArray<real,nquadx> &quad_wts_phys)
{
  real ll_corner_x = (i+this->topology->i_beg)*this->dx + this->xc - this->Lx/2. - this->dx/2.;
  for (int nqx=0; nqx<nquadx; nqx++)
  {
    quad_pts_phys(nqx).x = ll_corner_x + this->x_quad_pts_ref(nqx) * this->dx;
    quad_wts_phys(nqx) = this->x_quad_wts_ref(nqx) * this->dx;
  }
}






// *********************   2D  ******************************/

enum class LINE_INTEGRAL_TYPE { TANGENT, NORMAL };


template<> struct coords<2>
{
  real x=0., y=0.;
};

template<> struct vec<2>
{
  real u=0., v=0.;
};

template<uint nquadx, uint nquady, uint nquadz> class Geometry<2,nquadx,nquady,nquadz> {

public:

  const Topology *topology;
  bool is_initialized;
  bool straight;
  
  Geometry();
  Geometry( const Geometry<2,nquadx,nquady,nquadz> &geom) = delete;
  Geometry& operator=( const Geometry<2,nquadx,nquady,nquadz> &geom) = delete;

  SArray<real,nquadx> x_quad_pts_ref;
  SArray<real,nquadx> x_quad_wts_ref;
  SArray<real,nquady> y_quad_pts_ref;
  SArray<real,nquady> y_quad_wts_ref;

  void initialize(const Topology &topo);

  virtual void YAKL_INLINE get_0form_quad_pts_wts(int i, int j, SArray<coords<2>,1> &quad_pts_phys, SArray<real,1> &quad_wts_phys) {};
  virtual void YAKL_INLINE get_1form_quad_pts_wts(int i, int j, SArray<coords<2>,nquadx> &x_quad_pts_phys, SArray<real,nquadx> &x_quad_wts_phys, SArray<coords<2>,nquady> &y_quad_pts_phys, SArray<real,nquady> &y_quad_wts_phys) {};
  virtual void YAKL_INLINE get_2form_quad_pts_wts(int i, int j, SArray<coords<2>,nquadx,nquady> &quad_pts_phys, SArray<real,nquadx,nquady> &quad_wts_phys) {};
  virtual void YAKL_INLINE get_edge_tangents(int i, int j, SArray<vec<2>,nquadx> &x_tangents, SArray<vec<2>,nquady> &y_tangents) {};
  virtual void YAKL_INLINE get_edge_normals(int i, int j, SArray<vec<2>,nquadx> &x_normals, SArray<vec<2>,nquady> &y_normals) {};

  //virtual real YAKL_INLINE get_area_0form(int k, int j, int i) {};
  //virtual real YAKL_INLINE get_area_1form(int l, int k, int j, int i) {};
  //virtual real YAKL_INLINE get_area_2form(int k, int j, int i) {};
  virtual real YAKL_INLINE get_area_lform(int l, int d, int k, int j, int i) {};

  YAKL_INLINE void set_0form_values(real (*initial_value_function)(real, real), Field &field, int ndof);
  YAKL_INLINE void set_1form_values(vec<2> (*initial_value_function)(real, real), Field &field, int ndof, LINE_INTEGRAL_TYPE line_type);
  YAKL_INLINE void set_2form_values(real (*initial_value_function)(real, real), Field &field, int ndof);

};

template<uint nquadx, uint nquady, uint nquadz> class UniformRectangularTwistedGeometry<2,nquadx,nquady,nquadz>: public Geometry<2,nquadx,nquady,nquadz> {


public:

  real dx, dy, dz;
  real Lx, Ly, Lz;
  real xc, yc, zc;
  void initialize(const Topology &topo, const ModelParameters &params);

  void printinfo();

  void YAKL_INLINE get_0form_quad_pts_wts(int i, int j, SArray<coords<2>,1> &quad_pts_phys, SArray<real,1> &quad_wts_phys);
  void YAKL_INLINE get_1form_quad_pts_wts(int i, int j, SArray<coords<2>,nquadx> &x_quad_pts_phys, SArray<real,nquadx> &x_quad_wts_phys, SArray<coords<2>,nquady> &y_quad_pts_phys, SArray<real,nquady> &y_quad_wts_phys);
  void YAKL_INLINE get_2form_quad_pts_wts(int i, int j, SArray<coords<2>,nquadx,nquady> &quad_pts_phys, SArray<real,nquadx,nquady> &quad_wts_phys);
  void YAKL_INLINE get_edge_tangents(int i, int j, SArray<vec<2>,nquadx> &x_tangents, SArray<vec<2>,nquady> &y_tangents);
  void YAKL_INLINE get_edge_normals(int i, int j, SArray<vec<2>,nquadx> &x_normals, SArray<vec<2>,nquady> &y_normals);

  //real YAKL_INLINE get_area_0form(int k, int j, int i);
  //real YAKL_INLINE get_area_1form(int l, int k, int j, int i);
  //real YAKL_INLINE get_area_2form(int k, int j, int i);
  real YAKL_INLINE get_area_lform(int l, int d, int k, int j, int i);

};

template<uint nquadx, uint nquady, uint nquadz> class UniformRectangularStraightGeometry<2,nquadx,nquady,nquadz>: public Geometry<2,nquadx,nquady,nquadz> {


public:

  real dx, dy, dz;
  real Lx, Ly, Lz;
  real xc, yc, zc;
  void initialize(const Topology &topo, const ModelParameters &params);

  void printinfo();

  void YAKL_INLINE get_0form_quad_pts_wts(int i, int j, SArray<coords<2>,1> &quad_pts_phys, SArray<real,1> &quad_wts_phys);
  void YAKL_INLINE get_1form_quad_pts_wts(int i, int j, SArray<coords<2>,nquadx> &x_quad_pts_phys, SArray<real,nquadx> &x_quad_wts_phys, SArray<coords<2>,nquady> &y_quad_pts_phys, SArray<real,nquady> &y_quad_wts_phys);
  void YAKL_INLINE get_2form_quad_pts_wts(int i, int j, SArray<coords<2>,nquadx,nquady> &quad_pts_phys, SArray<real,nquadx,nquady> &quad_wts_phys);
  void YAKL_INLINE get_edge_tangents(int i, int j, SArray<vec<2>,nquadx> &x_tangents, SArray<vec<2>,nquady> &y_tangents);
  void YAKL_INLINE get_edge_normals(int i, int j, SArray<vec<2>,nquadx> &x_normals, SArray<vec<2>,nquady> &y_normals);

  //real YAKL_INLINE get_area_0form(int k, int j, int i);
  //real YAKL_INLINE get_area_1form(int l, int k, int j, int i);
  //real YAKL_INLINE get_area_2form(int k, int j, int i);
  real YAKL_INLINE get_area_lform(int l, int d, int k, int j, int i);
  
};





template<uint nquadx, uint nquady, uint nquadz> Geometry<2,nquadx,nquady,nquadz>::Geometry()
    {
      this->is_initialized = false;
      std::cout << "CREATED GEOMETRY\n";
    }

  template<uint nquadx, uint nquady, uint nquadz> void Geometry<2,nquadx,nquady,nquadz>::initialize(const Topology &topo)
  {

  this->topology = &topo;
  set_ref_quad_pts_wts<real,nquadx>(this->x_quad_pts_ref, this->x_quad_wts_ref);
  set_ref_quad_pts_wts<real,nquady>(this->y_quad_pts_ref, this->y_quad_wts_ref);
}


  template<uint nquadx, uint nquady, uint nquadz> YAKL_INLINE void Geometry<2,nquadx,nquady,nquadz>::set_0form_values(real (*initial_value_function)(real, real), Field &field, int ndof)
{

  int is = this->topology->is;
  int js = this->topology->js;
  int ks = this->topology->ks;

  SArray<coords<2>,1> quad_pts_phys;
  SArray<real,1> quad_wts_phys;

  yakl::parallel_for("Set0FormValues", this->topology->n_cells, YAKL_LAMBDA (int iGlob) {
    int k, j, i;
    yakl::unpackIndices(iGlob, this->topology->n_cells_z, this->topology->n_cells_y, this->topology->n_cells_x ,k, j, i);
      get_0form_quad_pts_wts(i, j, quad_pts_phys, quad_wts_phys);
      field.data(ndof, k+ks, j+js, i+is) = initial_value_function(quad_pts_phys(0).x, quad_pts_phys(0).y) * quad_wts_phys(0);
  });
}

  template<uint nquadx, uint nquady, uint nquadz> YAKL_INLINE void Geometry<2,nquadx,nquady,nquadz>::set_2form_values(real (*initial_value_function)(real, real), Field &field, int ndof)
{

  SArray<coords<2>,nquadx,nquady> quad_pts_phys;
  SArray<real,nquadx,nquady> quad_wts_phys;

  int is = this->topology->is;
  int js = this->topology->js;
  int ks = this->topology->ks;
  real tempval;
  int offset = field.get_offset(2);

  yakl::parallel_for("Set2FormValues", this->topology->n_cells, YAKL_LAMBDA (int iGlob) {
    int k, j, i;
    yakl::unpackIndices(iGlob, this->topology->n_cells_z, this->topology->n_cells_y, this->topology->n_cells_x ,k, j, i);

    get_2form_quad_pts_wts(i, j, quad_pts_phys, quad_wts_phys);
    tempval = 0.0;
    for (int nqx=0; nqx<nquadx; nqx++) {
      for (int nqy=0; nqy<nquady; nqy++) {
          tempval = tempval + initial_value_function(quad_pts_phys(nqx, nqy).x, quad_pts_phys(nqx, nqy).y) * quad_wts_phys(nqx, nqy);
        }}
      field.data(ndof+offset, k+ks, j+js, i+is) = tempval;
  });
}



  template<uint nquadx, uint nquady, uint nquadz>  YAKL_INLINE void Geometry<2,nquadx,nquady,nquadz>::set_1form_values(vec<2> (*initial_value_function)(real, real), Field &field, int ndof, LINE_INTEGRAL_TYPE line_type)
  {
      SArray<coords<2>,nquadx> x_quad_pts_phys;
      SArray<real,nquadx> x_quad_wts_phys;
      SArray<coords<2>,nquady> y_quad_pts_phys;
      SArray<real,nquady> y_quad_wts_phys;

      SArray<vec<2>,nquadx> x_line_vec;
      SArray<vec<2>,nquady> y_line_vec;

      vec<2> initval;

      int is = this->topology->is;
      int js = this->topology->js;
      int ks = this->topology->ks;
      real tempval;
      int offset = field.get_offset(1);

      //twisted edges
      int yedge_offset = 0;
      int xedge_offset = field.ndof1;
      
      //compared to twisted edges in 2D, straight edges are stored x/y (V, -U) instead of y/x
      if (this->straight) {
        yedge_offset = field.ndof1;
        xedge_offset = 0;
      }

      yakl::parallel_for("Set1FormValues", this->topology->n_cells , YAKL_LAMBDA (int iGlob) {
        int k, j, i;
        yakl::unpackIndices(iGlob, this->topology->n_cells_z, this->topology->n_cells_y, this->topology->n_cells_x ,k, j, i);

        get_1form_quad_pts_wts(i, j, x_quad_pts_phys, x_quad_wts_phys, y_quad_pts_phys, y_quad_wts_phys);

        if (line_type == LINE_INTEGRAL_TYPE::TANGENT) {get_edge_tangents(i, j, x_line_vec, y_line_vec);}
        if (line_type == LINE_INTEGRAL_TYPE::NORMAL) {get_edge_normals(i, j, x_line_vec, y_line_vec);}

      // y-edge
        tempval = 0.0;
      for (int nqy=0; nqy<nquady; nqy++) {
        initval = initial_value_function(y_quad_pts_phys(nqy).x, y_quad_pts_phys(nqy).y);
        tempval = tempval + (initval.u * y_line_vec(nqy).u + initval.v * y_line_vec(nqy).v) * y_quad_wts_phys(nqy);
      }
      field.data(ndof+offset+yedge_offset, k+ks, j+js, i+is) = tempval;

            // x-edge
            tempval = 0.0;
            for (int nqx=0; nqx<nquadx; nqx++) {
              initval = initial_value_function(x_quad_pts_phys(nqx).x, x_quad_pts_phys(nqx).y);
              tempval = tempval + (initval.u * x_line_vec(nqx).u + initval.v * x_line_vec(nqx).v) * x_quad_wts_phys(nqx);
            }
            field.data(ndof+xedge_offset+offset, k+ks, j+js, i+is) = tempval;
      });
  }

template<uint nquadx, uint nquady, uint nquadz> void UniformRectangularTwistedGeometry<2,nquadx,nquady,nquadz>::initialize(const Topology &topo, const ModelParameters &params)
{

Geometry<2,nquadx,nquady,nquadz>::initialize(topo);

this->Lx = params.xlen;
this->Ly = params.ylen;
this->Lz = 0.;

this->xc = params.xc;
this->yc = params.yc;
this->zc = 0.;

this->dx = params.xlen/params.nx_glob;
this->dy = params.ylen/params.ny_glob;
this->dz = 0.;

this->straight = false;

this->is_initialized = true;
}

template<uint nquadx, uint nquady, uint nquadz> void UniformRectangularTwistedGeometry<2,nquadx,nquady,nquadz>::printinfo()
{
  std::cout << "uniform rectangular geometry 2D info: twisted \n" << std::flush;
  std::cout << "Lx " << this->Lx << " Ly " << this->Ly << " Lz " << this->Lz << "\n" << std::flush;
  std::cout << "xc " << this->xc << " yc " << this->yc << " zc " << this->zc << "\n" << std::flush;
  std::cout << "dx " << this->dx << " dy " << this->dy << " dz " << this->dz << "\n" << std::flush;
}


template<uint nquadx, uint nquady, uint nquadz> real YAKL_INLINE UniformRectangularTwistedGeometry<2,nquadx,nquady,nquadz>::get_area_lform(int l, int d, int k, int j, int i)
{
  if (l == 0) {return 1.;}
  if (l == 1) {
    if (d==0) {return this->dy;}
    if (d==1) {return this->dx;}
    }
  if (l == 2) {return this->dx * this->dy;}
}

// template<uint nquadx, uint nquady, uint nquadz> real YAKL_INLINE UniformRectangularTwistedGeometry<2,nquadx,nquady,nquadz>::get_area_2form(int k, int j, int i)
// {
//   return this->dx * this->dy;
// }
// 
// template<uint nquadx, uint nquady, uint nquadz> real YAKL_INLINE UniformRectangularTwistedGeometry<2,nquadx,nquady,nquadz>::get_area_0form(int k, int j, int i)
// {
//   return 1.;
// }
// 
// template<uint nquadx, uint nquady, uint nquadz> real YAKL_INLINE UniformRectangularTwistedGeometry<2,nquadx,nquady,nquadz>::get_area_1form(int l, int k, int j, int i)
// {
//   if (l==0) {return this->dy;};
//   if (l==1) {return this->dx;};
// }

template<uint nquadx, uint nquady, uint nquadz> void YAKL_INLINE UniformRectangularTwistedGeometry<2,nquadx,nquady,nquadz>::get_0form_quad_pts_wts(int i, int j, SArray<coords<2>,1> &quad_pts_phys, SArray<real,1> &quad_wts_phys)
{
  quad_pts_phys(0).x = (i+this->topology->i_beg)*this->dx + this->xc - this->Lx/2.;
  quad_pts_phys(0).y = (j+this->topology->j_beg)*this->dy + this->yc - this->Ly/2.;
  quad_wts_phys(0) = 1.;
}

template<uint nquadx, uint nquady, uint nquadz> void YAKL_INLINE UniformRectangularTwistedGeometry<2,nquadx,nquady,nquadz>::get_2form_quad_pts_wts(int i, int j, SArray<coords<2>,nquadx,nquady> &quad_pts_phys, SArray<real,nquadx,nquady> &quad_wts_phys)
{
  real ll_corner_x = (i+this->topology->i_beg)*this->dx + this->xc - this->Lx/2.;
  real ll_corner_y = (j+this->topology->j_beg)*this->dy + this->yc - this->Ly/2.;
  for (int nqx=0; nqx<nquadx; nqx++)
  for (int nqy=0; nqy<nquady; nqy++)
  {{
    quad_pts_phys(nqx, nqy).x = ll_corner_x + this->x_quad_pts_ref(nqx) * this->dx;
    quad_pts_phys(nqx, nqy).y = ll_corner_y + this->y_quad_pts_ref(nqy) * this->dy;
    quad_wts_phys(nqx, nqy) = this->x_quad_wts_ref(nqx) * this->dx * this->y_quad_wts_ref(nqy) * this->dy;
  }}
}

template<uint nquadx, uint nquady, uint nquadz> void YAKL_INLINE UniformRectangularTwistedGeometry<2,nquadx,nquady,nquadz>::get_1form_quad_pts_wts(int i, int j, SArray<coords<2>,nquadx> &x_quad_pts_phys, SArray<real,nquadx> &x_quad_wts_phys, SArray<coords<2>,nquady> &y_quad_pts_phys, SArray<real,nquady> &y_quad_wts_phys)
{
    real ll_corner_x = (i+this->topology->i_beg)*this->dx + this->xc - this->Lx/2.;
    real ll_corner_y = (j+this->topology->j_beg)*this->dy + this->yc - this->Ly/2.;
    for (int nqx=0; nqx<nquadx; nqx++) {
        x_quad_pts_phys(nqx).x = ll_corner_x + this->x_quad_pts_ref(nqx) * this->dx;
        x_quad_pts_phys(nqx).y = ll_corner_y;
        x_quad_wts_phys(nqx) = this->x_quad_wts_ref(nqx) * this->dx;
    }
    for (int nqy=0; nqy<nquady; nqy++) {
        y_quad_pts_phys(nqy).x = ll_corner_x;
        y_quad_pts_phys(nqy).y = ll_corner_y + this->y_quad_pts_ref(nqy) * this->dy;
        y_quad_wts_phys(nqy) = this->y_quad_wts_ref(nqy) * this->dy;
    }
}


template<uint nquadx, uint nquady, uint nquadz> void UniformRectangularStraightGeometry<2,nquadx,nquady,nquadz>::initialize(const Topology &topo, const ModelParameters &params)
{

Geometry<2,nquadx,nquady,nquadz>::initialize(topo);

this->Lx = params.xlen;
this->Ly = params.ylen;
this->Lz = 0.;

this->xc = params.xc;
this->yc = params.yc;
this->zc = 0.;

// CHECK THESE...
this->dx = params.xlen/params.nx_glob;
this->dy = params.ylen/params.ny_glob;
this->dz = 0.;

this->straight = true;

this->is_initialized = true;
}

template<uint nquadx, uint nquady, uint nquadz> void UniformRectangularStraightGeometry<2,nquadx,nquady,nquadz>::printinfo()
{
  std::cout << "uniform rectangular geometry 2D info: straight \n" << std::flush;
  std::cout << "Lx " << this->Lx << " Ly " << this->Ly << " Lz " << this->Lz << "\n" << std::flush;
  std::cout << "xc " << this->xc << " yc " << this->yc << " zc " << this->zc << "\n" << std::flush;
  std::cout << "dx " << this->dx << " dy " << this->dy << " dz " << this->dz << "\n" << std::flush;
}

template<uint nquadx, uint nquady, uint nquadz> real YAKL_INLINE UniformRectangularStraightGeometry<2,nquadx,nquady,nquadz>::get_area_lform(int l, int d, int k, int j, int i)
{
  if (l == 0) {return 1.;}
  if (l == 1) {
    if (d==0) {return this->dx;}
    if (d==1) {return this->dy;}
    }
  if (l == 2) {return this->dx * this->dy;}
}


// template<uint nquadx, uint nquady, uint nquadz> real YAKL_INLINE UniformRectangularStraightGeometry<2,nquadx,nquady,nquadz>::get_area_2form(int k, int j, int i)
// {
//   return this->dx * this->dy;
// }
// 
// template<uint nquadx, uint nquady, uint nquadz> real YAKL_INLINE UniformRectangularStraightGeometry<2,nquadx,nquady,nquadz>::get_area_0form(int k, int j, int i)
// {
//   return 1.;
// }
// 
// template<uint nquadx, uint nquady, uint nquadz> real YAKL_INLINE UniformRectangularStraightGeometry<2,nquadx,nquady,nquadz>::get_area_1form(int l, int k, int j, int i)
// {
//   if (l==0) {return this->dx;};
//   if (l==1) {return this->dy;};
// }

template<uint nquadx, uint nquady, uint nquadz> void YAKL_INLINE UniformRectangularStraightGeometry<2,nquadx,nquady,nquadz>::get_0form_quad_pts_wts(int i, int j, SArray<coords<2>,1> &quad_pts_phys, SArray<real,1> &quad_wts_phys)
{
  quad_pts_phys(0).x = (i+this->topology->i_beg)*this->dx + this->xc - this->Lx/2. + this->dx/2.;
  quad_pts_phys(0).y = (j+this->topology->j_beg)*this->dy + this->yc - this->Ly/2. + this->dy/2.;
  quad_wts_phys(0) = 1.;
}

template<uint nquadx, uint nquady, uint nquadz> void YAKL_INLINE UniformRectangularStraightGeometry<2,nquadx,nquady,nquadz>::get_2form_quad_pts_wts(int i, int j, SArray<coords<2>,nquadx,nquady> &quad_pts_phys, SArray<real,nquadx,nquady> &quad_wts_phys)
{
    real ll_corner_x = (i+this->topology->i_beg)*this->dx + this->xc - this->Lx/2. - this->dx/2.;
    real ll_corner_y = (j+this->topology->j_beg)*this->dy + this->yc - this->Ly/2. - this->dy/2.;
    for (int nqx=0; nqx<nquadx; nqx++)
    for (int nqy=0; nqy<nquady; nqy++)
    {{
      quad_pts_phys(nqx, nqy).x = ll_corner_x + this->x_quad_pts_ref(nqx) * this->dx;
      quad_pts_phys(nqx, nqy).y = ll_corner_y + this->y_quad_pts_ref(nqy) * this->dy;
      quad_wts_phys(nqx, nqy) = this->x_quad_wts_ref(nqx) * this->dx * this->y_quad_wts_ref(nqy) * this->dy;
    }}
}

template<uint nquadx, uint nquady, uint nquadz> void YAKL_INLINE UniformRectangularStraightGeometry<2,nquadx,nquady,nquadz>::get_1form_quad_pts_wts(int i, int j, SArray<coords<2>,nquadx> &x_quad_pts_phys, SArray<real,nquadx> &x_quad_wts_phys, SArray<coords<2>,nquady> &y_quad_pts_phys, SArray<real,nquady> &y_quad_wts_phys)
{
      real ll_corner_x = (i+this->topology->i_beg)*this->dx + this->xc - this->Lx/2. - this->dx/2.;
      real ll_corner_y = (j+this->topology->j_beg)*this->dy + this->yc - this->Ly/2. - this->dy/2.;
      for (int nqx=0; nqx<nquadx; nqx++) {
          x_quad_pts_phys(nqx).x = ll_corner_x + this->x_quad_pts_ref(nqx) * this->dx;
          x_quad_pts_phys(nqx).y = ll_corner_y + this->dy;
          x_quad_wts_phys(nqx) = this->x_quad_wts_ref(nqx) * this->dx;
      }
      for (int nqy=0; nqy<nquady; nqy++) {
          y_quad_pts_phys(nqy).x = ll_corner_x + this->dx;
          y_quad_pts_phys(nqy).y = ll_corner_y + this->y_quad_pts_ref(nqy) * this->dy;
          y_quad_wts_phys(nqy) = this->y_quad_wts_ref(nqy) * this->dy;
      }
}


template<uint nquadx, uint nquady, uint nquadz> void YAKL_INLINE UniformRectangularTwistedGeometry<2,nquadx,nquady,nquadz>::get_edge_tangents(int i, int j, SArray<vec<2>,nquadx> &x_tangents, SArray<vec<2>,nquady> &y_tangents)
{
    for (int nqx=0; nqx<nquadx; nqx++)
    {
        x_tangents(nqx).u = -1.;
        x_tangents(nqx).v = 0.;
    }
    for (int nqy=0; nqy<nquady; nqy++)
    {
      y_tangents(nqy).u = 0.;
      y_tangents(nqy).v = 1.;
    }
}


template<uint nquadx, uint nquady, uint nquadz> void YAKL_INLINE UniformRectangularTwistedGeometry<2,nquadx,nquady,nquadz>::get_edge_normals(int i, int j, SArray<vec<2>,nquadx> &x_normals, SArray<vec<2>,nquady> &y_normals)
{
    for (int nqx=0; nqx<nquadx; nqx++)
    {
        x_normals(nqx).u = 0.;
        x_normals(nqx).v = 1.;
    }
    for (int nqy=0; nqy<nquady; nqy++)
    {
      y_normals(nqy).u = 1.;
      y_normals(nqy).v = 0.;
    }
}

// For straight edges in 2D, the twisted tangent is the straight normal, and the straight normal is the twisted tangent

template<uint nquadx, uint nquady, uint nquadz> void YAKL_INLINE UniformRectangularStraightGeometry<2,nquadx,nquady,nquadz>::get_edge_tangents(int i, int j, SArray<vec<2>,nquadx> &x_tangents, SArray<vec<2>,nquady> &y_tangents)
{
    for (int nqx=0; nqx<nquadx; nqx++)
    {
        x_tangents(nqx).u = 1.;
        x_tangents(nqx).v = 0.;
    }
    for (int nqy=0; nqy<nquady; nqy++)
    {
      y_tangents(nqy).u = 0.;
      y_tangents(nqy).v = 1.;
    }
}


template<uint nquadx, uint nquady, uint nquadz> void YAKL_INLINE UniformRectangularStraightGeometry<2,nquadx,nquady,nquadz>::get_edge_normals(int i, int j, SArray<vec<2>,nquadx> &x_normals, SArray<vec<2>,nquady> &y_normals)
{
    for (int nqx=0; nqx<nquadx; nqx++)
    {
        x_normals(nqx).u = 0.;
        x_normals(nqx).v = 1.;
    }
    for (int nqy=0; nqy<nquady; nqy++)
    {
      y_normals(nqy).u = -1.;
      y_normals(nqy).v = 0.;
    }
}


  // *********************   3D  ******************************/
// Unclear if this is really useful
// All of our 3D models used extruded...

//   template<> struct coords<3>
//   {
//     real x=0., y=0., z=0.;
//   };
// 
//   template<> struct vec<3>
//   {
//     real u=0., v=0., w=0.;
//   };
// 
// 
//   template<uint nquadx, uint nquady, uint nquadz> class Geometry<3,nquadx,nquady,nquadz> {
// 
//   public:
// 
//     const Topology *topology;
//     bool is_initialized;
//     Geometry();
//     Geometry( const Geometry<3,nquadx,nquady,nquadz> &geom) = delete;
//     Geometry& operator=( const Geometry<3,nquadx,nquady,nquadz> &geom) = delete;
// 
//     SArray<real,nquadx> x_quad_pts_ref;
//     SArray<real,nquadx> x_quad_wts_ref;
//     SArray<real,nquady> y_quad_pts_ref;
//     SArray<real,nquady> y_quad_wts_ref;
//     SArray<real,nquadz> z_quad_pts_ref;
//     SArray<real,nquadz> z_quad_wts_ref;
// 
//     void initialize(const Topology &topo);
// 
//     virtual void YAKL_INLINE get_primal_0form_quad_pts_wts(int i, int j, int k, SArray<coords<3>,1> &quad_pts_phys, SArray<real,1> &quad_wts_phys) {};
//     virtual void YAKL_INLINE get_primal_1form_quad_pts_wts(int i, int j, int k,
//       SArray<coords<3>,nquadx> &x_quad_pts_phys, SArray<real,nquadx> &x_quad_wts_phys,
//       SArray<coords<3>,nquady> &y_quad_pts_phys, SArray<real,nquady> &y_quad_wts_phys,
//       SArray<coords<3>,nquadz> &z_quad_pts_phys, SArray<real,nquadz> &z_quad_wts_phys) {};
//     virtual void YAKL_INLINE get_primal_2form_quad_pts_wts(int i, int j, int k,
//       SArray<coords<3>,nquadx,nquady> &xy_quad_pts_phys, SArray<real,nquadx,nquady> &xy_quad_wts_phys,
//       SArray<coords<3>,nquady,nquadz> &yz_quad_pts_phys, SArray<real,nquady,nquadz> &yz_quad_wts_phys,
//       SArray<coords<3>,nquadx,nquadz> &xz_quad_pts_phys, SArray<real,nquadz,nquadz> &xz_quad_wts_phys) {};
//     virtual void YAKL_INLINE get_primal_3form_quad_pts_wts(int i, int j, int k, SArray<coords<3>,nquadx,nquady,nquadz> &quad_pts_phys, SArray<real,nquadx,nquady,nquadz> &quad_wts_phys) {};
//     virtual void YAKL_INLINE get_primal_edge_tangents(int i, int j, int k, SArray<vec<3>,nquadx> &x_tangents, SArray<vec<3>,nquady> &y_tangents, SArray<vec<3>,nquadz> &z_tangents) {};
//     virtual void YAKL_INLINE get_primal_surface_normals(int i, int j, int k, SArray<vec<3>,nquadx,nquady> &xy_normals, SArray<vec<3>,nquady,nquadz> &yz_normals, SArray<vec<3>,nquadx,nquadz> &xz_normals) {};
// 
//     virtual void YAKL_INLINE get_dual_0form_quad_pts_wts(int i, int j, int k, SArray<coords<3>,1> &quad_pts_phys, SArray<real,1> &quad_wts_phys) {};
//     virtual void YAKL_INLINE get_dual_1form_quad_pts_wts(int i, int j, int k,
//       SArray<coords<3>,nquadx> &x_quad_pts_phys, SArray<real,nquadx> &x_quad_wts_phys,
//       SArray<coords<3>,nquady> &y_quad_pts_phys, SArray<real,nquady> &y_quad_wts_phys,
//       SArray<coords<3>,nquadz> &z_quad_pts_phys, SArray<real,nquadz> &z_quad_wts_phys) {};
//     virtual void YAKL_INLINE get_dual_2form_quad_pts_wts(int i, int j, int k,
//       SArray<coords<3>,nquadx,nquady> &xy_quad_pts_phys, SArray<real,nquadx,nquady> &xy_quad_wts_phys,
//       SArray<coords<3>,nquady,nquadz> &yz_quad_pts_phys, SArray<real,nquady,nquadz> &yz_quad_wts_phys,
//       SArray<coords<3>,nquadx,nquadz> &xz_quad_pts_phys, SArray<real,nquadz,nquadz> &xz_quad_wts_phys) {};
//     virtual void YAKL_INLINE get_dual_3form_quad_pts_wts(int i, int j, int k, SArray<coords<3>,nquadx,nquady,nquadz> &quad_pts_phys, SArray<real,nquadx,nquady,nquadz> &quad_wts_phys) {};
//     virtual void YAKL_INLINE get_dual_edge_tangents(int i, int j, int k, SArray<vec<3>,nquadx> &x_tangents, SArray<vec<3>,nquady> &y_tangents, SArray<vec<3>,nquadz> &z_tangents) {};
//     virtual void YAKL_INLINE get_dual_surface_normals(int i, int j, int k, SArray<vec<3>,nquadx,nquady> &xy_normals, SArray<vec<3>,nquady,nquadz> &yz_normals, SArray<vec<3>,nquadx,nquadz> &xz_normals) {};
// 
// 
//     virtual real YAKL_INLINE get_J_cell(int k, int j, int i) {};
//     virtual real YAKL_INLINE get_J_dual_cell(int k, int j, int i) {};
//     virtual real YAKL_INLINE get_H_edge(int l, int k, int j, int i) {};
// 
//     YAKL_INLINE void set_primal_0form_values(real (*initial_value_function)(real, real, real), Field &field, int ndof);
//     YAKL_INLINE void set_primal_1form_values(vec<3> (*initial_value_function)(real, real, real), Field &field, int ndof);
//     YAKL_INLINE void set_primal_2form_values(vec<3> (*initial_value_function)(real, real, real), Field &field, int ndof);
//     YAKL_INLINE void set_primal_3form_values(real (*initial_value_function)(real, real, real), Field &field, int ndof);
// 
//     YAKL_INLINE void set_dual_0form_values(real (*initial_value_function)(real, real, real), Field &field, int ndof);
//     YAKL_INLINE void set_dual_1form_values(vec<3> (*initial_value_function)(real, real, real), Field &field, int ndof);
//     YAKL_INLINE void set_dual_2form_values(vec<3> (*initial_value_function)(real, real, real), Field &field, int ndof);
//     YAKL_INLINE void set_dual_3form_values(real (*initial_value_function)(real, real, real), Field &field, int ndof);
// 
//   };
// 
//   template<uint nquadx, uint nquady, uint nquadz> class UniformRectangularGeometry<3,nquadx,nquady,nquadz>: public Geometry<3,nquadx,nquady,nquadz> {
// 
// 
//   public:
// 
//     real dx, dy, dz;
//     real Lx, Ly, Lz;
//     real xc, yc, zc;
//     void initialize(const Topology &topo, const ModelParameters &params);
// 
//     void printinfo();
// 
//     void YAKL_INLINE get_primal_0form_quad_pts_wts(int i, int j, int k, SArray<coords<3>,1> &quad_pts_phys, SArray<real,1> &quad_wts_phys);
//     void YAKL_INLINE get_primal_1form_quad_pts_wts(int i, int j, int k,
//       SArray<coords<3>,nquadx> &x_quad_pts_phys, SArray<real,nquadx> &x_quad_wts_phys,
//       SArray<coords<3>,nquady> &y_quad_pts_phys, SArray<real,nquady> &y_quad_wts_phys,
//       SArray<coords<3>,nquadz> &z_quad_pts_phys, SArray<real,nquadz> &z_quad_wts_phys);
//     void YAKL_INLINE get_primal_2form_quad_pts_wts(int i, int j, int k,
//       SArray<coords<3>,nquadx,nquady> &xy_quad_pts_phys, SArray<real,nquadx,nquady> &xy_quad_wts_phys,
//       SArray<coords<3>,nquady,nquadz> &yz_quad_pts_phys, SArray<real,nquady,nquadz> &yz_quad_wts_phys,
//       SArray<coords<3>,nquadx,nquadz> &xz_quad_pts_phys, SArray<real,nquadz,nquadz> &xz_quad_wts_phys);
//     void YAKL_INLINE get_primal_3form_quad_pts_wts(int i, int j, int k, SArray<coords<3>,nquadx,nquady,nquadz> &quad_pts_phys, SArray<real,nquadx,nquady,nquadz> &quad_wts_phys);
//     void YAKL_INLINE get_primal_edge_tangents(int i, int j, int k, SArray<vec<3>,nquadx> &x_tangents, SArray<vec<3>,nquady> &y_tangents, SArray<vec<3>,nquadz> &z_tangents);
//     void YAKL_INLINE get_primal_surface_normals(int i, int j, int k, SArray<vec<3>,nquadx,nquady> &xy_normals, SArray<vec<3>,nquady,nquadz> &yz_normals, SArray<vec<3>,nquadx,nquadz> &xz_normals);
// 
//     void YAKL_INLINE get_dual_0form_quad_pts_wts(int i, int j, int k, SArray<coords<3>,1> &quad_pts_phys, SArray<real,1> &quad_wts_phys);
//     void YAKL_INLINE get_dual_1form_quad_pts_wts(int i, int j, int k,
//       SArray<coords<3>,nquadx> &x_quad_pts_phys, SArray<real,nquadx> &x_quad_wts_phys,
//       SArray<coords<3>,nquady> &y_quad_pts_phys, SArray<real,nquady> &y_quad_wts_phys,
//       SArray<coords<3>,nquadz> &z_quad_pts_phys, SArray<real,nquadz> &z_quad_wts_phys);
//     void YAKL_INLINE get_dual_2form_quad_pts_wts(int i, int j, int k,
//       SArray<coords<3>,nquadx,nquady> &xy_quad_pts_phys, SArray<real,nquadx,nquady> &xy_quad_wts_phys,
//       SArray<coords<3>,nquady,nquadz> &yz_quad_pts_phys, SArray<real,nquady,nquadz> &yz_quad_wts_phys,
//       SArray<coords<3>,nquadx,nquadz> &xz_quad_pts_phys, SArray<real,nquadz,nquadz> &xz_quad_wts_phys);
//     void YAKL_INLINE get_dual_3form_quad_pts_wts(int i, int j, int k, SArray<coords<3>,nquadx,nquady,nquadz> &quad_pts_phys, SArray<real,nquadx,nquady,nquadz> &quad_wts_phys);
//     void YAKL_INLINE get_dual_edge_tangents(int i, int j, int k, SArray<vec<3>,nquadx> &x_tangents, SArray<vec<3>,nquady> &y_tangents, SArray<vec<3>,nquadz> &z_tangents);
//     void YAKL_INLINE get_dual_surface_normals(int i, int j, int k, SArray<vec<3>,nquadx,nquady> &xy_normals, SArray<vec<3>,nquady,nquadz> &yz_normals, SArray<vec<3>,nquadx,nquadz> &xz_normals);
// 
//     real YAKL_INLINE get_J_cell(int k, int j, int i);
//     real YAKL_INLINE get_J_dual_cell(int k, int j, int i);
//     real YAKL_INLINE get_H_edge(int l, int k, int j, int i);
// 
//   };
// 
// 
// 
// 
// 
// 
// 
//   template<uint nquadx, uint nquady, uint nquadz> Geometry<3,nquadx,nquady,nquadz>::Geometry()
//       {
//         this->is_initialized = false;
//         std::cout << "CREATED GEOMETRY\n";
//       }
// 
//     template<uint nquadx, uint nquady, uint nquadz> void Geometry<3,nquadx,nquady,nquadz>::initialize(const Topology &topo)
//     {
// 
//     this->topology = &topo;
//     set_ref_quad_pts_wts<real,nquadx>(this->x_quad_pts_ref, this->x_quad_wts_ref);
//     set_ref_quad_pts_wts<real,nquady>(this->y_quad_pts_ref, this->y_quad_wts_ref);
//     set_ref_quad_pts_wts<real,nquadz>(this->z_quad_pts_ref, this->z_quad_wts_ref);
//   }
// 
// 
//     template<uint nquadx, uint nquady, uint nquadz> YAKL_INLINE void Geometry<3,nquadx,nquady,nquadz>::set_primal_0form_values(real (*initial_value_function)(real, real, real), Field &field, int ndof)
//   {
// 
//     int is = this->topology->is;
//     int js = this->topology->js;
//     int ks = this->topology->ks;
// 
//     SArray<coords<3>,1> quad_pts_phys;
//     SArray<real,1> quad_wts_phys;
// 
//     yakl::parallel_for("Set0FormValues", this->topology->n_cells, YAKL_LAMBDA (int iGlob) {
//       int k, j, i;
//       yakl::unpackIndices(iGlob, this->topology->n_cells_z, this->topology->n_cells_y, this->topology->n_cells_x ,k, j, i);
//         get_primal_0form_quad_pts_wts(i, j, quad_pts_phys, quad_wts_phys);
//         field.data(ndof, k+ks, j+js, i+is) = initial_value_function(quad_pts_phys(0).x, quad_pts_phys(0).y, quad_pts_phys(0).z) * quad_wts_phys(0);
//     });
//   }
// 
//     template<uint nquadx, uint nquady, uint nquadz> YAKL_INLINE void Geometry<3,nquadx,nquady,nquadz>::set_primal_3form_values(real (*initial_value_function)(real, real, real), Field &field, int ndof)
//   {
// 
//     SArray<coords<3>,nquadx,nquady,nquadz> quad_pts_phys;
//     SArray<real,nquadx,nquady,nquadz> quad_wts_phys;
// 
//     int is = this->topology->is;
//     int js = this->topology->js;
//     int ks = this->topology->ks;
//     real tempval;
//     int offset = field.get_offset(3);
// 
//     yakl::parallel_for("Set3FormValues", this->topology->n_cells, YAKL_LAMBDA (int iGlob) {
//       int k, j, i;
//       yakl::unpackIndices(iGlob, this->topology->n_cells_z, this->topology->n_cells_y, this->topology->n_cells_x ,k, j, i);
// 
//       get_primal_3form_quad_pts_wts(i, j, k, quad_pts_phys, quad_wts_phys);
//       tempval = 0.0;
//       for (int nqx=0; nqx<nquadx; nqx++) {
//         for (int nqy=0; nqy<nquady; nqy++) {
//           for (int nqz=0; nqz<nquadz; nqz++) {
//             tempval = tempval + initial_value_function(quad_pts_phys(nqx, nqy, nqz).x, quad_pts_phys(nqx, nqy, nqz).y, quad_pts_phys(nqx, nqy, nqz).z) * quad_wts_phys(nqx, nqy, nqz);
//           }}}
//         field.data(ndof+offset, k+ks, j+js, i+is) = tempval;
//     });
//   }
// 
// 
// 
//     template<uint nquadx, uint nquady, uint nquadz>  YAKL_INLINE void Geometry<3,nquadx,nquady,nquadz>::set_primal_1form_values(vec<3> (*initial_value_function)(real, real, real), Field &field, int ndof)
//     {
//         SArray<coords<3>,nquadx> x_quad_pts_phys;
//         SArray<real,nquadx> x_quad_wts_phys;
//         SArray<coords<3>,nquady> y_quad_pts_phys;
//         SArray<real,nquady> y_quad_wts_phys;
//         SArray<coords<3>,nquadz> z_quad_pts_phys;
//         SArray<real,nquadz> z_quad_wts_phys;
// 
//         SArray<vec<3>,nquadx> x_tangents;
//         SArray<vec<3>,nquady> y_tangents;
//         SArray<vec<3>,nquadz> z_tangents;
// 
//         vec<3> initval;
// 
//         int is = this->topology->is;
//         int js = this->topology->js;
//         int ks = this->topology->ks;
//         real tempval;
//         int offset = field.get_offset(1);
// 
//         yakl::parallel_for("Set1FormValues", this->topology->n_cells , YAKL_LAMBDA (int iGlob) {
//           int k, j, i;
//           yakl::unpackIndices(iGlob, this->topology->n_cells_z, this->topology->n_cells_y, this->topology->n_cells_x ,k, j, i);
// 
// 
//           get_primal_1form_quad_pts_wts(i, j, k, x_quad_pts_phys, x_quad_wts_phys, y_quad_pts_phys, y_quad_wts_phys, z_quad_pts_phys, z_quad_wts_phys);
//           get_primal_edge_tangents(i, j, k, x_tangents, y_tangents, z_tangents);
// 
//           // z-edge
//             tempval = 0.0;
//           for (int nqz=0; nqz<nquadz; nqz++) {
//             initval = initial_value_function(z_quad_pts_phys(nqz).x, z_quad_pts_phys(nqz).y, z_quad_pts_phys(nqz).z);
//             tempval = tempval + (initval.u * z_tangents(nqz).u + initval.v * z_tangents(nqz).v + initval.w * z_tangents(nqz).w) * z_quad_wts_phys(nqz);
//           }
//           field.data(ndof+offset, k+ks, j+js, i+is) = tempval;
// 
//         // y-edge
//           tempval = 0.0;
//         for (int nqy=0; nqy<nquady; nqy++) {
//           initval = initial_value_function(y_quad_pts_phys(nqy).x, y_quad_pts_phys(nqy).y, y_quad_pts_phys(nqy).z);
//           tempval = tempval + (initval.u * y_tangents(nqy).u + initval.v * y_tangents(nqy).v + initval.w * y_tangents(nqy).w) * y_quad_wts_phys(nqy);
//         }
//         field.data(ndof+field.ndof1+offset, k+ks, j+js, i+is) = tempval;
// 
//               // x-edge
//               tempval = 0.0;
//               for (int nqx=0; nqx<nquadx; nqx++) {
//                 initval = initial_value_function(x_quad_pts_phys(nqx).x, x_quad_pts_phys(nqx).y, x_quad_pts_phys(nqx).z);
//                 tempval = tempval + (initval.u * x_tangents(nqx).u + initval.v * x_tangents(nqx).v + initval.w * x_tangents(nqx).w) * x_quad_wts_phys(nqx);
//               }
//               field.data(ndof+2*field.ndof1+offset, k+ks, j+js, i+is) = tempval;
//         });
//     }
// 
// 
//   template<uint nquadx, uint nquady, uint nquadz> YAKL_INLINE void Geometry<3,nquadx,nquady,nquadz>::set_primal_2form_values(vec<3> (*initial_value_function)(real, real, real), Field &field, int ndof)
//   {
// 
//   SArray<coords<3>,nquadx,nquady> xy_quad_pts_phys;
//   SArray<real,nquadx,nquady> xy_quad_wts_phys;
//   SArray<coords<3>,nquady,nquadz> yz_quad_pts_phys;
//   SArray<real,nquady,nquadz> yz_quad_wts_phys;
//   SArray<coords<3>,nquadx,nquadz> xz_quad_pts_phys;
//   SArray<real,nquadx,nquadz> xz_quad_wts_phys;
// 
//   SArray<vec<3>,nquadx,nquady> xy_normals;
//   SArray<vec<3>,nquady,nquadz> yz_normals;
//   SArray<vec<3>,nquadx,nquadz> xz_normals;
// 
//   vec<3> initval;
// 
//   int is = this->topology->is;
//   int js = this->topology->js;
//   int ks = this->topology->ks;
//   real tempval;
//   int offset = field.get_offset(2);
// 
//   yakl::parallel_for("Set2FormValues", this->topology->n_cells , YAKL_LAMBDA (int iGlob) {
//     int k, j, i;
//     yakl::unpackIndices(iGlob, this->topology->n_cells_z, this->topology->n_cells_y, this->topology->n_cells_x ,k, j, i);
// 
//     get_primal_2form_quad_pts_wts(i, j, k, xy_quad_pts_phys, xy_quad_wts_phys, yz_quad_pts_phys, yz_quad_wts_phys, xz_quad_pts_phys, xz_quad_wts_phys);
//     get_primal_surface_normals(i, j, k, xy_normals, yz_normals, xz_normals);
// 
//     // yz surface
//     tempval = 0.0;
//     for (int nqy=0; nqy<nquady; nqy++) {
//       for (int nqz=0; nqz<nquadz; nqz++) {
//         initval = initial_value_function(yz_quad_pts_phys(nqy,nqz).x, yz_quad_pts_phys(nqy,nqz).y, yz_quad_pts_phys(nqy,nqz).z);
//         tempval = tempval + (initval.u * yz_normals(nqy, nqz).u + initval.v * yz_normals(nqy, nqz).v + initval.w * yz_normals(nqy, nqz).w) * yz_quad_wts_phys(nqy, nqz);
//     }}
//     field.data(ndof+offset, k+ks, j+js, i+is) = tempval;
// 
//         // xz surface
//         tempval = 0.0;
//         for (int nqx=0; nqx<nquadx; nqx++) {
//           for (int nqz=0; nqz<nquadz; nqz++) {
//             initval = initial_value_function(xz_quad_pts_phys(nqx,nqz).x, xz_quad_pts_phys(nqx,nqz).y, xz_quad_pts_phys(nqx,nqz).z);
//             tempval = tempval + (initval.u * xz_normals(nqx, nqz).u + initval.v * xz_normals(nqx, nqz).v + initval.w * xz_normals(nqx, nqz).w)  * xz_quad_wts_phys(nqx,nqz);
//         }}
//         field.data(ndof+field.ndof2+offset, k+ks, j+js, i+is) = tempval;
// 
// 
//       // xy surface
//       tempval = 0.0;
//       for (int nqx=0; nqx<nquadx; nqx++) {
//         for (int nqy=0; nqy<nquady; nqy++) {
//           initval = initial_value_function(xy_quad_pts_phys(nqx,nqy).x, xy_quad_pts_phys(nqx,nqy).y, xy_quad_pts_phys(nqx,nqy).z);
//           tempval = tempval + (initval.u * xy_normals(nqx, nqy).u + initval.v * xy_normals(nqx, nqy).v + initval.w * xy_normals(nqx, nqy).w)  * xy_quad_wts_phys(nqx,nqy);
//       }}
//       field.data(ndof+2*field.ndof2+offset, k+ks, j+js, i+is) = tempval;
// 
//   });
// }
// 
// 
// template<uint nquadx, uint nquady, uint nquadz> YAKL_INLINE void Geometry<3,nquadx,nquady,nquadz>::set_dual_0form_values(real (*initial_value_function)(real, real, real), Field &field, int ndof)
// {
// 
// int is = this->topology->is;
// int js = this->topology->js;
// int ks = this->topology->ks;
// 
// SArray<coords<3>,1> quad_pts_phys;
// SArray<real,1> quad_wts_phys;
// int offset = field.get_offset(3);
// 
// yakl::parallel_for("Set0FormValues", this->topology->n_cells, YAKL_LAMBDA (int iGlob) {
//   int k, j, i;
//   yakl::unpackIndices(iGlob, this->topology->n_cells_z, this->topology->n_cells_y, this->topology->n_cells_x ,k, j, i);
//     get_dual_0form_quad_pts_wts(i, j, quad_pts_phys, quad_wts_phys);
//     field.data(ndof+offset, k+ks, j+js, i+is) = initial_value_function(quad_pts_phys(0).x, quad_pts_phys(0).y, quad_pts_phys(0).z) * quad_wts_phys(0);
// });
// }
// 
// template<uint nquadx, uint nquady, uint nquadz> YAKL_INLINE void Geometry<3,nquadx,nquady,nquadz>::set_dual_3form_values(real (*initial_value_function)(real, real, real), Field &field, int ndof)
// {
// 
// SArray<coords<3>,nquadx,nquady,nquadz> quad_pts_phys;
// SArray<real,nquadx,nquady,nquadz> quad_wts_phys;
// 
// int is = this->topology->is;
// int js = this->topology->js;
// int ks = this->topology->ks;
// real tempval;
// 
// yakl::parallel_for("Set3FormValues", this->topology->n_cells, YAKL_LAMBDA (int iGlob) {
//   int k, j, i;
//   yakl::unpackIndices(iGlob, this->topology->n_cells_z, this->topology->n_cells_y, this->topology->n_cells_x ,k, j, i);
// 
//   get_dual_3form_quad_pts_wts(i, j, k, quad_pts_phys, quad_wts_phys);
//   tempval = 0.0;
//   for (int nqx=0; nqx<nquadx; nqx++) {
//     for (int nqy=0; nqy<nquady; nqy++) {
//       for (int nqz=0; nqz<nquadz; nqz++) {
//         tempval = tempval + initial_value_function(quad_pts_phys(nqx, nqy, nqz).x, quad_pts_phys(nqx, nqy, nqz).y, quad_pts_phys(nqx, nqy, nqz).z) * quad_wts_phys(nqx, nqy, nqz);
//       }}}
//     field.data(ndof, k+ks, j+js, i+is) = tempval;
// });
// }
// 
// 
// 
// template<uint nquadx, uint nquady, uint nquadz>  YAKL_INLINE void Geometry<3,nquadx,nquady,nquadz>::set_dual_1form_values(vec<3> (*initial_value_function)(real, real, real), Field &field, int ndof)
// {
//     SArray<coords<3>,nquadx> x_quad_pts_phys;
//     SArray<real,nquadx> x_quad_wts_phys;
//     SArray<coords<3>,nquady> y_quad_pts_phys;
//     SArray<real,nquady> y_quad_wts_phys;
//     SArray<coords<3>,nquadz> z_quad_pts_phys;
//     SArray<real,nquadz> z_quad_wts_phys;
// 
//     SArray<vec<3>,nquadx> x_tangents;
//     SArray<vec<3>,nquady> y_tangents;
//     SArray<vec<3>,nquadz> z_tangents;
// 
//     vec<3> initval;
// 
//     int is = this->topology->is;
//     int js = this->topology->js;
//     int ks = this->topology->ks;
//     real tempval;
//     int offset = field.get_offset(2);
// 
//     yakl::parallel_for("Set1FormValues", this->topology->n_cells , YAKL_LAMBDA (int iGlob) {
//       int k, j, i;
//       yakl::unpackIndices(iGlob, this->topology->n_cells_z, this->topology->n_cells_y, this->topology->n_cells_x ,k, j, i);
// 
// 
//       get_dual_1form_quad_pts_wts(i, j, k, x_quad_pts_phys, x_quad_wts_phys, y_quad_pts_phys, y_quad_wts_phys, z_quad_pts_phys, z_quad_wts_phys);
//       get_dual_edge_tangents(i, j, k, x_tangents, y_tangents, z_tangents);
// 
//       // z-edge
//         tempval = 0.0;
//       for (int nqz=0; nqz<nquadz; nqz++) {
//         initval = initial_value_function(z_quad_pts_phys(nqz).x, z_quad_pts_phys(nqz).y, z_quad_pts_phys(nqz).z);
//         tempval = tempval + (initval.u * z_tangents(nqz).u + initval.v * z_tangents(nqz).v + initval.w * z_tangents(nqz).w) * z_quad_wts_phys(nqz);
//       }
//       field.data(ndof+2*field.ndof2+offset, k+ks, j+js, i+is) = tempval;
// 
//     // y-edge
//       tempval = 0.0;
//     for (int nqy=0; nqy<nquady; nqy++) {
//       initval = initial_value_function(y_quad_pts_phys(nqy).x, y_quad_pts_phys(nqy).y, y_quad_pts_phys(nqy).z);
//       tempval = tempval + (initval.u * y_tangents(nqy).u + initval.v * y_tangents(nqy).v + initval.w * y_tangents(nqy).w) * y_quad_wts_phys(nqy);
//     }
//     field.data(ndof+1*field.ndof2+offset, k+ks, j+js, i+is) = tempval;
// 
//           // x-edge
//           tempval = 0.0;
//           for (int nqx=0; nqx<nquadx; nqx++) {
//             initval = initial_value_function(x_quad_pts_phys(nqx).x, x_quad_pts_phys(nqx).y, x_quad_pts_phys(nqx).z);
//             tempval = tempval + (initval.u * x_tangents(nqx).u + initval.v * x_tangents(nqx).v + initval.w * x_tangents(nqx).w) * x_quad_wts_phys(nqx);
//           }
//           field.data(ndof+0*field.ndof2+offset, k+ks, j+js, i+is) = tempval;
//     });
// }
// 
// 
// template<uint nquadx, uint nquady, uint nquadz> YAKL_INLINE void Geometry<3,nquadx,nquady,nquadz>::set_dual_2form_values(vec<3> (*initial_value_function)(real, real, real), Field &field, int ndof)
// {
// 
// SArray<coords<3>,nquadx,nquady> xy_quad_pts_phys;
// SArray<real,nquadx,nquady> xy_quad_wts_phys;
// SArray<coords<3>,nquady,nquadz> yz_quad_pts_phys;
// SArray<real,nquady,nquadz> yz_quad_wts_phys;
// SArray<coords<3>,nquadx,nquadz> xz_quad_pts_phys;
// SArray<real,nquadx,nquadz> xz_quad_wts_phys;
// 
// SArray<vec<3>,nquadx,nquady> xy_normals;
// SArray<vec<3>,nquady,nquadz> yz_normals;
// SArray<vec<3>,nquadx,nquadz> xz_normals;
// 
// vec<3> initval;
// 
// int is = this->topology->is;
// int js = this->topology->js;
// int ks = this->topology->ks;
// real tempval;
// int offset = field.get_offset(1);
// 
// yakl::parallel_for("Set2FormValues", this->topology->n_cells , YAKL_LAMBDA (int iGlob) {
// int k, j, i;
// yakl::unpackIndices(iGlob, this->topology->n_cells_z, this->topology->n_cells_y, this->topology->n_cells_x ,k, j, i);
// 
// get_dual_2form_quad_pts_wts(i, j, k, xy_quad_pts_phys, xy_quad_wts_phys, yz_quad_pts_phys, yz_quad_wts_phys, xz_quad_pts_phys, xz_quad_wts_phys);
// get_dual_surface_normals(i, j, k, xy_normals, yz_normals, xz_normals);
// 
// // yz surface
// tempval = 0.0;
// for (int nqy=0; nqy<nquady; nqy++) {
//   for (int nqz=0; nqz<nquadz; nqz++) {
//     initval = initial_value_function(yz_quad_pts_phys(nqy,nqz).x, yz_quad_pts_phys(nqy,nqz).y, yz_quad_pts_phys(nqy,nqz).z);
//     tempval = tempval + (initval.u * yz_normals(nqy, nqz).u + initval.v * yz_normals(nqy, nqz).v + initval.w * yz_normals(nqy, nqz).w) * yz_quad_wts_phys(nqy, nqz);
// }}
// field.data(ndof+2*field.ndof1+offset, k+ks, j+js, i+is) = tempval;
// 
//     // xz surface
//     tempval = 0.0;
//     for (int nqx=0; nqx<nquadx; nqx++) {
//       for (int nqz=0; nqz<nquadz; nqz++) {
//         initval = initial_value_function(xz_quad_pts_phys(nqx,nqz).x, xz_quad_pts_phys(nqx,nqz).y, xz_quad_pts_phys(nqx,nqz).z);
//         tempval = tempval + (initval.u * xz_normals(nqx, nqz).u + initval.v * xz_normals(nqx, nqz).v + initval.w * xz_normals(nqx, nqz).w)  * xz_quad_wts_phys(nqx,nqz);
//     }}
//     field.data(ndof+1*field.ndof1+offset, k+ks, j+js, i+is) = tempval;
// 
// 
//   // xy surface
//   tempval = 0.0;
//   for (int nqx=0; nqx<nquadx; nqx++) {
//     for (int nqy=0; nqy<nquady; nqy++) {
//       initval = initial_value_function(xy_quad_pts_phys(nqx,nqy).x, xy_quad_pts_phys(nqx,nqy).y, xy_quad_pts_phys(nqx,nqy).z);
//       tempval = tempval + (initval.u * xy_normals(nqx, nqy).u + initval.v * xy_normals(nqx, nqy).v + initval.w * xy_normals(nqx, nqy).w)  * xy_quad_wts_phys(nqx,nqy);
//   }}
//   field.data(ndof+0*field.ndof1+offset, k+ks, j+js, i+is) = tempval;
// 
// });
// }
// 
// 
//   template<uint nquadx, uint nquady, uint nquadz> void UniformRectangularGeometry<3,nquadx,nquady,nquadz>::initialize(const Topology &topo, const ModelParameters &params)
//   {
// 
//   Geometry<3,nquadx,nquady,nquadz>::initialize(topo);
// 
//   this->Lx = params.xlen;
//   this->Ly = params.ylen;
//   this->Lz = params.zlen;
// 
//   this->xc = params.xc;
//   this->yc = params.yc;
//   this->zc = params.zc;
// 
//   this->dx = params.xlen/params.nx_glob;
//   this->dy = params.ylen/params.ny_glob;
//   this->dz = params.zlen/params.nz_glob;
// 
//   this->is_initialized = true;
//   }
// 
//   template<uint nquadx, uint nquady, uint nquadz> void UniformRectangularGeometry<3,nquadx,nquady,nquadz>::printinfo()
//   {
//     std::cout << "uniform rectangular geometry 3D info\n" << std::flush;
//     std::cout << "Lx " << this->Lx << " Ly " << this->Ly << " Lz " << this->Lz << "\n" << std::flush;
//     std::cout << "xc " << this->xc << " yc " << this->yc << " zc " << this->zc << "\n" << std::flush;
//     std::cout << "dx " << this->dx << " dy " << this->dy << " dz " << this->dz << "\n" << std::flush;
//   }
// 
//   template<uint nquadx, uint nquady, uint nquadz> real YAKL_INLINE UniformRectangularGeometry<3,nquadx,nquady,nquadz>::get_J_cell(int k, int j, int i)
//   {
//     return this->dx * this->dy * this->dz;
//   }
// 
// 
//   template<uint nquadx, uint nquady, uint nquadz> real YAKL_INLINE UniformRectangularGeometry<3,nquadx,nquady,nquadz>::get_J_dual_cell(int k, int j, int i)
//   {
//     return this->dx * this->dy * this->dz;
//   }
// 
//   template<uint nquadx, uint nquady, uint nquadz> real YAKL_INLINE UniformRectangularGeometry<3,nquadx,nquady,nquadz>::get_H_edge(int l, int k, int j, int i)
//   {
//     if (l==0) {return this->dy*this->dz/this->dx;};
//     if (l==1) {return this->dx*this->dz/this->dy;};
//     if (l==2) {return this->dx*this->dy/this->dz;};
//   }
// 
//   template<uint nquadx, uint nquady, uint nquadz> void YAKL_INLINE UniformRectangularGeometry<3,nquadx,nquady,nquadz>::get_primal_0form_quad_pts_wts(int i, int j, int k, SArray<coords<3>,1> &quad_pts_phys, SArray<real,1> &quad_wts_phys)
//   {
//     quad_pts_phys(0).x = (i+this->topology->i_beg)*this->dx + this->xc - this->Lx/2.;
//     quad_pts_phys(0).y = (j+this->topology->j_beg)*this->dy + this->yc - this->Ly/2.;
//     quad_pts_phys(0).z = (k+this->topology->k_beg)*this->dz + this->zc - this->Lz/2.;
//     quad_wts_phys(0) = 1.;
//   }
// 
//   template<uint nquadx, uint nquady, uint nquadz> void YAKL_INLINE UniformRectangularGeometry<3,nquadx,nquady,nquadz>::get_primal_3form_quad_pts_wts(int i, int j, int k, SArray<coords<3>,nquadx,nquady,nquadz> &quad_pts_phys, SArray<real,nquadx,nquady,nquadz> &quad_wts_phys)
//   {
//     real ll_corner_x = (i+this->topology->i_beg)*this->dx + this->xc - this->Lx/2.;
//     real ll_corner_y = (j+this->topology->j_beg)*this->dy + this->yc - this->Ly/2.;
//     real ll_corner_z = (k+this->topology->k_beg)*this->dz + this->zc - this->Lz/2.;
//     for (int nqx=0; nqx<nquadx; nqx++)
//     for (int nqy=0; nqy<nquady; nqy++)
//     for (int nqz=0; nqz<nquadz; nqz++)
//     {{{
//       quad_pts_phys(nqx, nqy, nqz).x = ll_corner_x + this->x_quad_pts_ref(nqx) * this->dx;
//       quad_pts_phys(nqx, nqy, nqz).y = ll_corner_y + this->y_quad_pts_ref(nqy) * this->dy;
//       quad_pts_phys(nqx, nqy, nqz).z = ll_corner_z + this->z_quad_pts_ref(nqz) * this->dz;
//       quad_wts_phys(nqx, nqy, nqz) = this->x_quad_wts_ref(nqx) * this->dx * this->y_quad_wts_ref(nqy) * this->dy * this->z_quad_wts_ref(nqz) * this->dz;
//     }}}
//   }
// 
//   template<uint nquadx, uint nquady, uint nquadz> void YAKL_INLINE UniformRectangularGeometry<3,nquadx,nquady,nquadz>::get_primal_1form_quad_pts_wts(int i, int j, int k,
//     SArray<coords<3>,nquadx> &x_quad_pts_phys, SArray<real,nquadx> &x_quad_wts_phys,
//     SArray<coords<3>,nquady> &y_quad_pts_phys, SArray<real,nquady> &y_quad_wts_phys,
//     SArray<coords<3>,nquadz> &z_quad_pts_phys, SArray<real,nquadz> &z_quad_wts_phys)
//   {
//       real ll_corner_x = (i+this->topology->i_beg)*this->dx + this->xc - this->Lx/2.;
//       real ll_corner_y = (j+this->topology->j_beg)*this->dy + this->yc - this->Ly/2.;
//       real ll_corner_z = (k+this->topology->k_beg)*this->dz + this->zc - this->Lz/2.;
//       for (int nqx=0; nqx<nquadx; nqx++) {
//           x_quad_pts_phys(nqx).x = ll_corner_x + this->x_quad_pts_ref(nqx) * this->dx;
//           x_quad_pts_phys(nqx).y = ll_corner_y;
//           x_quad_pts_phys(nqx).z = ll_corner_z;
//           x_quad_wts_phys(nqx) = this->x_quad_wts_ref(nqx) * this->dx;
//       }
//       for (int nqy=0; nqy<nquady; nqy++) {
//           y_quad_pts_phys(nqy).x = ll_corner_x;
//           y_quad_pts_phys(nqy).y = ll_corner_y + this->y_quad_pts_ref(nqy) * this->dy;
//           y_quad_pts_phys(nqy).z = ll_corner_z;
//           y_quad_wts_phys(nqy) = this->y_quad_wts_ref(nqy) * this->dy;
//       }
//       for (int nqz=0; nqz<nquadz; nqz++) {
//           z_quad_pts_phys(nqz).x = ll_corner_x;
//           z_quad_pts_phys(nqz).y = ll_corner_y;
//           z_quad_pts_phys(nqz).z = ll_corner_z + this->z_quad_pts_ref(nqz) * this->dz;
//           z_quad_wts_phys(nqz) = this->z_quad_wts_ref(nqz) * this->dz;
//       }
//   }
// 
// 
//   template<uint nquadx, uint nquady, uint nquadz> void YAKL_INLINE UniformRectangularGeometry<3,nquadx,nquady,nquadz>::get_primal_2form_quad_pts_wts(int i, int j, int k,
//     SArray<coords<3>,nquadx,nquady> &xy_quad_pts_phys, SArray<real,nquadx,nquady> &xy_quad_wts_phys,
//     SArray<coords<3>,nquady,nquadz> &yz_quad_pts_phys, SArray<real,nquady,nquadz> &yz_quad_wts_phys,
//     SArray<coords<3>,nquadx,nquadz> &xz_quad_pts_phys, SArray<real,nquadz,nquadz> &xz_quad_wts_phys)
//     {
//         real ll_corner_x = (i+this->topology->i_beg)*this->dx + this->xc - this->Lx/2.;
//         real ll_corner_y = (j+this->topology->j_beg)*this->dy + this->yc - this->Ly/2.;
//         real ll_corner_z = (k+this->topology->k_beg)*this->dz + this->zc - this->Lz/2.;
// 
//         for (int nqx=0; nqx<nquadx; nqx++) {
//           for (int nqy=0; nqy<nquady; nqy++) {
//               xy_quad_pts_phys(nqx, nqy).x = ll_corner_x + this->x_quad_pts_ref(nqx) * this->dx;
//               xy_quad_pts_phys(nqx, nqy).y = ll_corner_y + this->y_quad_pts_ref(nqy) * this->dy;
//               xy_quad_pts_phys(nqx, nqy).z = ll_corner_z;
//               xy_quad_wts_phys(nqx, nqy) = this->x_quad_wts_ref(nqx) * this->y_quad_wts_ref(nqy) * this->dx * this->dy;
//         }}
// 
//           for (int nqy=0; nqy<nquady; nqy++) {
//             for (int nqz=0; nqz<nquadz; nqz++) {
//               yz_quad_pts_phys(nqy, nqz).x = ll_corner_x;
//               yz_quad_pts_phys(nqy, nqz).y = ll_corner_y + this->y_quad_pts_ref(nqy) * this->dy;
//               yz_quad_pts_phys(nqy, nqz).z = ll_corner_z + this->z_quad_pts_ref(nqz) * this->dz;
//               yz_quad_wts_phys(nqy, nqz) = this->y_quad_wts_ref(nqy) * this->z_quad_wts_ref(nqz) * this->dy * this->dz;
//         }}
// 
//         for (int nqx=0; nqx<nquadx; nqx++) {
//             for (int nqz=0; nqz<nquadz; nqz++) {
//               xz_quad_pts_phys(nqx, nqz).x = ll_corner_x + this->x_quad_pts_ref(nqx) * this->dx;
//               xz_quad_pts_phys(nqx, nqz).y = ll_corner_y;
//               xz_quad_pts_phys(nqx, nqz).z = ll_corner_z + this->z_quad_pts_ref(nqz) * this->dz;
//               xz_quad_wts_phys(nqx, nqz) = this->x_quad_wts_ref(nqx) * this->z_quad_wts_ref(nqz) * this->dx * this->dz;
//         }}
//     }
// 
//   template<uint nquadx, uint nquady, uint nquadz> void YAKL_INLINE UniformRectangularGeometry<3,nquadx,nquady,nquadz>::get_primal_edge_tangents(int i, int j, int k, SArray<vec<3>,nquadx> &x_tangents, SArray<vec<3>,nquady> &y_tangents, SArray<vec<3>,nquadz> &z_tangents)
//   {
//       for (int nqx=0; nqx<nquadx; nqx++)
//       {
//           x_tangents(nqx).u = 1.;
//           x_tangents(nqx).v = 0.;
//           x_tangents(nqx).w = 0.;
//       }
//       for (int nqy=0; nqy<nquady; nqy++)
//       {
//         y_tangents(nqy).u = 0.;
//         y_tangents(nqy).v = 1.;
//         y_tangents(nqy).w = 0.;
//       }
//       for (int nqz=0; nqz<nquadz; nqz++)
//       {
//         z_tangents(nqz).u = 0.;
//         z_tangents(nqz).v = 0.;
//         z_tangents(nqz).w = 1.;
//       }
//   }
// 
// 
//   template<uint nquadx, uint nquady, uint nquadz> void YAKL_INLINE UniformRectangularGeometry<3,nquadx,nquady,nquadz>::get_primal_surface_normals(int i, int j, int k, SArray<vec<3>,nquadx,nquady> &xy_normals, SArray<vec<3>,nquady,nquadz> &yz_normals, SArray<vec<3>,nquadx,nquadz> &xz_normals)
//   {
//       for (int nqx=0; nqx<nquadx; nqx++) {
//         for (int nqy=0; nqy<nquady; nqy++) {
//             xy_normals(nqx, nqy).u = 0.;
//             xy_normals(nqx, nqy).v = 0.;
//             xy_normals(nqx, nqy).w = 1.;
//       }}
//         for (int nqy=0; nqy<nquady; nqy++) {
//           for (int nqz=0; nqz<nquadz; nqz++) {
//             yz_normals(nqy, nqz).u = 1;
//             yz_normals(nqy, nqz).v = 0;
//             yz_normals(nqy, nqz).w = 0;
//       }}
//       for (int nqx=0; nqx<nquadx; nqx++) {
//           for (int nqz=0; nqz<nquadz; nqz++) {
//             xz_normals(nqx, nqz).u = 0;
//             xz_normals(nqx, nqz).v = 1;
//             xz_normals(nqx, nqz).w = 0;
//       }}
//   }
// 
// 
// 
//   template<uint nquadx, uint nquady, uint nquadz> void YAKL_INLINE UniformRectangularGeometry<3,nquadx,nquady,nquadz>::get_dual_0form_quad_pts_wts(int i, int j, int k, SArray<coords<3>,1> &quad_pts_phys, SArray<real,1> &quad_wts_phys)
//   {
//     quad_pts_phys(0).x = (i+this->topology->i_beg)*this->dx + this->xc - this->Lx/2. + this->dx/2.;
//     quad_pts_phys(0).y = (j+this->topology->j_beg)*this->dy + this->yc - this->Ly/2. + this->dy/2.;
//     quad_pts_phys(0).z = (k+this->topology->k_beg)*this->dz + this->zc - this->Lz/2. + this->dz/2.;
//     quad_wts_phys(0) = 1.;
//   }
// 
//   template<uint nquadx, uint nquady, uint nquadz> void YAKL_INLINE UniformRectangularGeometry<3,nquadx,nquady,nquadz>::get_dual_3form_quad_pts_wts(int i, int j, int k, SArray<coords<3>,nquadx,nquady,nquadz> &quad_pts_phys, SArray<real,nquadx,nquady,nquadz> &quad_wts_phys)
//   {
//     real ll_corner_x = (i+this->topology->i_beg)*this->dx + this->xc - this->Lx/2. - this->dx/2.;
//     real ll_corner_y = (j+this->topology->j_beg)*this->dy + this->yc - this->Ly/2. - this->dy/2.;
//     real ll_corner_z = (k+this->topology->k_beg)*this->dz + this->zc - this->Lz/2. - this->dz/2.;
//     for (int nqx=0; nqx<nquadx; nqx++)
//     for (int nqy=0; nqy<nquady; nqy++)
//     for (int nqz=0; nqz<nquadz; nqz++)
//     {{{
//       quad_pts_phys(nqx, nqy, nqz).x = ll_corner_x + this->x_quad_pts_ref(nqx) * this->dx;
//       quad_pts_phys(nqx, nqy, nqz).y = ll_corner_y + this->y_quad_pts_ref(nqy) * this->dy;
//       quad_pts_phys(nqx, nqy, nqz).z = ll_corner_z + this->z_quad_pts_ref(nqz) * this->dz;
//       quad_wts_phys(nqx, nqy, nqz) = this->x_quad_wts_ref(nqx) * this->dx * this->y_quad_wts_ref(nqy) * this->dy * this->z_quad_wts_ref(nqz) * this->dz;
//     }}}
//   }
// 
//   template<uint nquadx, uint nquady, uint nquadz> void YAKL_INLINE UniformRectangularGeometry<3,nquadx,nquady,nquadz>::get_dual_1form_quad_pts_wts(int i, int j, int k,
//     SArray<coords<3>,nquadx> &x_quad_pts_phys, SArray<real,nquadx> &x_quad_wts_phys,
//     SArray<coords<3>,nquady> &y_quad_pts_phys, SArray<real,nquady> &y_quad_wts_phys,
//     SArray<coords<3>,nquadz> &z_quad_pts_phys, SArray<real,nquadz> &z_quad_wts_phys)
//   {
//       real ll_corner_x = (i+this->topology->i_beg)*this->dx + this->xc - this->Lx/2. - this->dx/2.;
//       real ll_corner_y = (j+this->topology->j_beg)*this->dy + this->yc - this->Ly/2. - this->dy/2.;
//       real ll_corner_z = (k+this->topology->k_beg)*this->dz + this->zc - this->Lz/2. - this->dz/2.;
//       for (int nqx=0; nqx<nquadx; nqx++) {
//           x_quad_pts_phys(nqx).x = ll_corner_x + this->x_quad_pts_ref(nqx) * this->dx;
//           x_quad_pts_phys(nqx).y = ll_corner_y + this->dy;
//           x_quad_pts_phys(nqx).z = ll_corner_z + this->dz;
//           x_quad_wts_phys(nqx) = this->x_quad_wts_ref(nqx) * this->dx;
//       }
//       for (int nqy=0; nqy<nquady; nqy++) {
//           y_quad_pts_phys(nqy).x = ll_corner_x + this->dx;
//           y_quad_pts_phys(nqy).y = ll_corner_y + this->y_quad_pts_ref(nqy) * this->dy;
//           y_quad_pts_phys(nqy).z = ll_corner_z + this->dz;
//           y_quad_wts_phys(nqy) = this->y_quad_wts_ref(nqy) * this->dy;
//       }
//       for (int nqz=0; nqz<nquadz; nqz++) {
//           z_quad_pts_phys(nqz).x = ll_corner_x + this->dx;
//           z_quad_pts_phys(nqz).y = ll_corner_y + this->dy;
//           z_quad_pts_phys(nqz).z = ll_corner_z + this->z_quad_pts_ref(nqz) * this->dz;
//           z_quad_wts_phys(nqz) = this->z_quad_wts_ref(nqz) * this->dz;
//       }
//   }
// 
// 
//   template<uint nquadx, uint nquady, uint nquadz> void YAKL_INLINE UniformRectangularGeometry<3,nquadx,nquady,nquadz>::get_dual_2form_quad_pts_wts(int i, int j, int k,
//     SArray<coords<3>,nquadx,nquady> &xy_quad_pts_phys, SArray<real,nquadx,nquady> &xy_quad_wts_phys,
//     SArray<coords<3>,nquady,nquadz> &yz_quad_pts_phys, SArray<real,nquady,nquadz> &yz_quad_wts_phys,
//     SArray<coords<3>,nquadx,nquadz> &xz_quad_pts_phys, SArray<real,nquadz,nquadz> &xz_quad_wts_phys)
//     {
//         real ll_corner_x = (i+this->topology->i_beg)*this->dx + this->xc - this->Lx/2. - this->dx/2.;
//         real ll_corner_y = (j+this->topology->j_beg)*this->dy + this->yc - this->Ly/2. - this->dy/2.;
//         real ll_corner_z = (k+this->topology->k_beg)*this->dz + this->zc - this->Lz/2. - this->dz/2.;
// 
//         for (int nqx=0; nqx<nquadx; nqx++) {
//           for (int nqy=0; nqy<nquady; nqy++) {
//               xy_quad_pts_phys(nqx, nqy).x = ll_corner_x + this->x_quad_pts_ref(nqx) * this->dx;
//               xy_quad_pts_phys(nqx, nqy).y = ll_corner_y + this->y_quad_pts_ref(nqy) * this->dy;
//               xy_quad_pts_phys(nqx, nqy).z = ll_corner_z + this->dz;
//               xy_quad_wts_phys(nqx, nqy) = this->x_quad_wts_ref(nqx) * this->y_quad_wts_ref(nqy) * this->dx * this->dy;
//         }}
// 
//           for (int nqy=0; nqy<nquady; nqy++) {
//             for (int nqz=0; nqz<nquadz; nqz++) {
//               yz_quad_pts_phys(nqy, nqz).x = ll_corner_x + this->dx;
//               yz_quad_pts_phys(nqy, nqz).y = ll_corner_y + this->y_quad_pts_ref(nqy) * this->dy;
//               yz_quad_pts_phys(nqy, nqz).z = ll_corner_z + this->z_quad_pts_ref(nqz) * this->dz;
//               yz_quad_wts_phys(nqy, nqz) = this->y_quad_wts_ref(nqy) * this->z_quad_wts_ref(nqz) * this->dy * this->dz;
//         }}
// 
//         for (int nqx=0; nqx<nquadx; nqx++) {
//             for (int nqz=0; nqz<nquadz; nqz++) {
//               xz_quad_pts_phys(nqx, nqz).x = ll_corner_x + this->x_quad_pts_ref(nqx) * this->dx;
//               xz_quad_pts_phys(nqx, nqz).y = ll_corner_y + this->dy;
//               xz_quad_pts_phys(nqx, nqz).z = ll_corner_z + this->z_quad_pts_ref(nqz) * this->dz;
//               xz_quad_wts_phys(nqx, nqz) = this->x_quad_wts_ref(nqx) * this->z_quad_wts_ref(nqz) * this->dx * this->dz;
//         }}
//     }
// 
//   template<uint nquadx, uint nquady, uint nquadz> void YAKL_INLINE UniformRectangularGeometry<3,nquadx,nquady,nquadz>::get_dual_edge_tangents(int i, int j, int k, SArray<vec<3>,nquadx> &x_tangents, SArray<vec<3>,nquady> &y_tangents, SArray<vec<3>,nquadz> &z_tangents)
//   {
//       for (int nqx=0; nqx<nquadx; nqx++)
//       {
//           x_tangents(nqx).u = 1.;
//           x_tangents(nqx).v = 0.;
//           x_tangents(nqx).w = 0.;
//       }
//       for (int nqy=0; nqy<nquady; nqy++)
//       {
//         y_tangents(nqy).u = 0.;
//         y_tangents(nqy).v = 1.;
//         y_tangents(nqy).w = 0.;
//       }
//       for (int nqz=0; nqz<nquadz; nqz++)
//       {
//         z_tangents(nqz).u = 0.;
//         z_tangents(nqz).v = 0.;
//         z_tangents(nqz).w = 1.;
//       }
//   }
// 
// 
//   template<uint nquadx, uint nquady, uint nquadz> void YAKL_INLINE UniformRectangularGeometry<3,nquadx,nquady,nquadz>::get_dual_surface_normals(int i, int j, int k, SArray<vec<3>,nquadx,nquady> &xy_normals, SArray<vec<3>,nquady,nquadz> &yz_normals, SArray<vec<3>,nquadx,nquadz> &xz_normals)
//   {
//       for (int nqx=0; nqx<nquadx; nqx++) {
//         for (int nqy=0; nqy<nquady; nqy++) {
//             xy_normals(nqx, nqy).u = 0.;
//             xy_normals(nqx, nqy).v = 0.;
//             xy_normals(nqx, nqy).w = 1.;
//       }}
//         for (int nqy=0; nqy<nquady; nqy++) {
//           for (int nqz=0; nqz<nquadz; nqz++) {
//             yz_normals(nqy, nqz).u = 1;
//             yz_normals(nqy, nqz).v = 0;
//             yz_normals(nqy, nqz).w = 0;
//       }}
//       for (int nqx=0; nqx<nquadx; nqx++) {
//           for (int nqz=0; nqz<nquadz; nqz++) {
//             xz_normals(nqx, nqz).u = 0;
//             xz_normals(nqx, nqz).v = 1;
//             xz_normals(nqx, nqz).w = 0;
//       }}
//   }


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

#endif
