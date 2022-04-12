#pragma once

#include "common.h"
#include "fields.h"
#include "topology.h"
#include "params.h"

template<int dim> struct coords {};
template<int dim> struct vec {};
template<int dim> struct vecext {};
template<int dim> struct coordsext {};

template<> struct coords<1>
{
  real x=0.0_fp;
};

template<> struct coords<2>
{
  real x=0.0_fp, y=0.0_fp;
};

template<> struct vec<2>
{
  real u=0.0_fp, v=0.0_fp;
};

template<> struct coordsext<2>
{
  real x=0.0_fp, z=0.0_fp;
};

template<> struct vecext<2>
{
  real u=0.0_fp, w=0.0_fp;
};

enum class LINE_INTEGRAL_TYPE { TANGENT, NORMAL };

template <class T, uint npts> void set_ref_quad_pts_wts(SArray<T,1,npts> &pts, SArray<T,1,npts> &wts)
{
if (npts == 1)
{
    pts(0) = 0.5_fp;
    wts(0) = 1.0_fp;
}

if (npts == 2)
{
    pts(0) = -1.0_fp/(2.0_fp*sqrt(3.)) + 0.5_fp;
    pts(1) = 1.0_fp/(2.0_fp*sqrt(3.)) + 0.5_fp;

    wts(0) = 0.5_fp;
    wts(1) = 0.5_fp;
}

if (npts == 3)
{
    pts(0) = -1.0_fp/2.0_fp*sqrt(3.0_fp/5.0_fp) + 0.5_fp;
    pts(1) = 0.5_fp;
    pts(2) = 0.5_fp*sqrt(3.0_fp/5.0_fp) + 0.5_fp;

    wts(0) = 5.0_fp/18.0_fp;
    wts(1) = 4.0_fp/9.0_fp;
    wts(2) = 5.0_fp/18.0_fp;
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
    pts(0) = -sqrt(3.0_fp/7.0_fp + 2.0_fp/7.0_fp * sqrt(6.0_fp/5.0_fp)) + 0.5_fp;
    pts(1) = -sqrt(3.0_fp/7.0_fp - 2.0_fp/7.0_fp * sqrt(6.0_fp/5.0_fp)) + 0.5_fp;
    pts(2) = sqrt(3.0_fp/7.0_fp - 2.0_fp/7.0_fp * sqrt(6.0_fp/5.0_fp)) + 0.5_fp;
    pts(3) = sqrt(3.0_fp/7.0_fp + 2.0_fp/7.0_fp * sqrt(6.0_fp/5.0_fp)) + 0.5_fp;

    wts(0) = (18.0_fp - sqrt(30.0_fp))/72.0_fp;
    wts(1) = (18.0_fp + sqrt(30.0_fp))/72.0_fp;
    wts(2) = (18.0_fp + sqrt(30.0_fp))/72.0_fp;
    wts(3) = (18.0_fp - sqrt(30.0_fp))/72.0_fp;
}

if (npts == 5)
{
    pts(0) = -1.0_fp/3.0_fp*sqrt(5.0_fp + 2.0_fp * sqrt(10.0_fp/7.0_fp)) + 0.5_fp;
    pts(1) = -1.0_fp/3.0_fp*sqrt(5.0_fp - 2.0_fp * sqrt(10.0_fp/7.0_fp)) + 0.5_fp;
    pts(2) = 0.5_fp;
    pts(3) = 1.0_fp/3.0_fp*sqrt(5.0_fp - 2.0_fp * sqrt(10.0_fp/7.0_fp)) + 0.5_fp;
    pts(4) = 1.0_fp/3.0_fp*sqrt(5.0_fp + 2.0_fp * sqrt(10.0_fp/7.0_fp)) + 0.5_fp;

    wts(0) = (322.0_fp - 13.0_fp*sqrt(70.0_fp))/1800.0_fp;
    wts(1) = (322.0_fp + 13.0_fp*sqrt(70.0_fp))/1800.0_fp;
    wts(2) = 64.0_fp/225.0_fp;
    wts(3) = (322.0_fp + 13.0_fp*sqrt(70.0_fp))/1800.0_fp;
    wts(4) = (322.0_fp - 13.0_fp*sqrt(70.0_fp))/1800.0_fp;
}

}

// EVENTUALLY SHOULD BE ABLE TO UNIFY THESE INTO A SINGLE GEOMETRY CLASS I THINK...

// *********************   2D Planar ******************************/

#ifndef _EXTRUDED   

class Geometry {

public:

  const Topology *topology;
  bool is_initialized;
  bool straight;
  
  Geometry();
  Geometry( const Geometry &geom) = delete;
  Geometry& operator=( const Geometry &geom) = delete;

  SArray<real,1,ic_quad_pts_x> x_quad_pts_ref;
  SArray<real,1,ic_quad_pts_x> x_quad_wts_ref;
  SArray<real,1,ic_quad_pts_y> y_quad_pts_ref;
  SArray<real,1,ic_quad_pts_y> y_quad_wts_ref;

  void initialize(const Topology &topo);

  virtual void YAKL_INLINE get_0form_quad_pts_wts(int i, int j, SArray<coords<2>,1,1> &quad_pts_phys, SArray<real,1,1> &quad_wts_phys) {};
  virtual void YAKL_INLINE get_1form_quad_pts_wts(int i, int j, SArray<coords<2>,1,ic_quad_pts_x> &x_quad_pts_phys, SArray<real,1,ic_quad_pts_x> &x_quad_wts_phys, SArray<coords<2>,1,ic_quad_pts_y> &y_quad_pts_phys, SArray<real,1,ic_quad_pts_y> &y_quad_wts_phys) {};
  virtual void YAKL_INLINE get_2form_quad_pts_wts(int i, int j, SArray<coords<2>,2,ic_quad_pts_x,ic_quad_pts_y> &quad_pts_phys, SArray<real,2,ic_quad_pts_x,ic_quad_pts_y> &quad_wts_phys) {};
  virtual void YAKL_INLINE get_edge_tangents(int i, int j, SArray<vec<2>,1,ic_quad_pts_x> &x_tangents, SArray<vec<2>,1,ic_quad_pts_y> &y_tangents) {};
  virtual void YAKL_INLINE get_edge_normals(int i, int j, SArray<vec<2>,1,ic_quad_pts_x> &x_normals, SArray<vec<2>,1,ic_quad_pts_y> &y_normals) {};

  virtual real YAKL_INLINE get_area_lform(int l, int d, int k, int j, int i) {};
  virtual real YAKL_INLINE get_area_00entity(int k, int j, int i) {};
  virtual real YAKL_INLINE get_area_11entity(int k, int j, int i) {};
  virtual real YAKL_INLINE get_area_10entity(int k, int j, int i) {};
  virtual real YAKL_INLINE get_area_01entity(int k, int j, int i) {};
  
  YAKL_INLINE void set_0form_values(real (*initial_value_function)(real, real), Field &field, int ndof);
  YAKL_INLINE void set_1form_values(vec<2> (*initial_value_function)(real, real), Field &field, int ndof, LINE_INTEGRAL_TYPE line_type);
  YAKL_INLINE void set_2form_values(real (*initial_value_function)(real, real), Field &field, int ndof);

  real dx, dy;
  real Lx, Ly;
  real xc, yc;
};

class UniformRectangularTwistedGeometry: public Geometry {


public:


  void initialize(const Topology &topo, const ModelParameters &params);

  void printinfo();

  void YAKL_INLINE get_0form_quad_pts_wts(int i, int j, SArray<coords<2>,1,1> &quad_pts_phys, SArray<real,1,1> &quad_wts_phys);
  void YAKL_INLINE get_1form_quad_pts_wts(int i, int j, SArray<coords<2>,1,ic_quad_pts_x> &x_quad_pts_phys, SArray<real,1,ic_quad_pts_x> &x_quad_wts_phys, SArray<coords<2>,1,ic_quad_pts_y> &y_quad_pts_phys, SArray<real,1,ic_quad_pts_y> &y_quad_wts_phys);
  void YAKL_INLINE get_2form_quad_pts_wts(int i, int j, SArray<coords<2>,2,ic_quad_pts_x,ic_quad_pts_y> &quad_pts_phys, SArray<real,2,ic_quad_pts_x,ic_quad_pts_y> &quad_wts_phys);
  void YAKL_INLINE get_edge_tangents(int i, int j, SArray<vec<2>,1,ic_quad_pts_x> &x_tangents, SArray<vec<2>,1,ic_quad_pts_y> &y_tangents);
  void YAKL_INLINE get_edge_normals(int i, int j, SArray<vec<2>,1,ic_quad_pts_x> &x_normals, SArray<vec<2>,1,ic_quad_pts_y> &y_normals);

  real YAKL_INLINE get_area_lform(int l, int d, int k, int j, int i);
  real YAKL_INLINE get_area_00entity(int k, int j, int i) {};
  real YAKL_INLINE get_area_11entity(int k, int j, int i) {};
  real YAKL_INLINE get_area_10entity(int k, int j, int i) {};
  real YAKL_INLINE get_area_01entity(int k, int j, int i) {};
};

class UniformRectangularStraightGeometry: public Geometry {


public:

  void initialize(const Topology &topo, const ModelParameters &params);

  void printinfo();

  void YAKL_INLINE get_0form_quad_pts_wts(int i, int j, SArray<coords<2>,1,1> &quad_pts_phys, SArray<real,1,1> &quad_wts_phys);
  void YAKL_INLINE get_1form_quad_pts_wts(int i, int j, SArray<coords<2>,1,ic_quad_pts_x> &x_quad_pts_phys, SArray<real,1,ic_quad_pts_x> &x_quad_wts_phys, SArray<coords<2>,1,ic_quad_pts_y> &y_quad_pts_phys, SArray<real,1,ic_quad_pts_y> &y_quad_wts_phys);
  void YAKL_INLINE get_2form_quad_pts_wts(int i, int j, SArray<coords<2>,2,ic_quad_pts_x,ic_quad_pts_y> &quad_pts_phys, SArray<real,2,ic_quad_pts_x,ic_quad_pts_y> &quad_wts_phys);
  void YAKL_INLINE get_edge_tangents(int i, int j, SArray<vec<2>,1,ic_quad_pts_x> &x_tangents, SArray<vec<2>,1,ic_quad_pts_y> &y_tangents);
  void YAKL_INLINE get_edge_normals(int i, int j, SArray<vec<2>,1,ic_quad_pts_x> &x_normals, SArray<vec<2>,1,ic_quad_pts_y> &y_normals);
  real YAKL_INLINE get_area_lform(int l, int d, int k, int j, int i);
  real YAKL_INLINE get_area_00entity(int k, int j, int i) {};
  real YAKL_INLINE get_area_11entity(int k, int j, int i) {};
  real YAKL_INLINE get_area_10entity(int k, int j, int i) {};
  real YAKL_INLINE get_area_01entity(int k, int j, int i) {};  
};





Geometry::Geometry()
    {
      this->is_initialized = false;
      std::cout << "CREATED GEOMETRY\n";
    }

void Geometry::initialize(const Topology &topo)
  {

  this->topology = &topo;
  set_ref_quad_pts_wts<real,ic_quad_pts_x>(this->x_quad_pts_ref, this->x_quad_wts_ref);
  set_ref_quad_pts_wts<real,ic_quad_pts_y>(this->y_quad_pts_ref, this->y_quad_wts_ref);
}


YAKL_INLINE void Geometry::set_0form_values(real (*initial_value_function)(real, real), Field &field, int ndof)
{

  int is = this->topology->is;
  int js = this->topology->js;
  int ks = this->topology->ks;


    parallel_for( Bounds<4>(this->topology->nl, this->topology->n_cells_y, this->topology->n_cells_x, this->topology->nens) , YAKL_LAMBDA(int k, int j, int i, int n) { 
      SArray<coords<2>,1,1> quad_pts_phys;
      SArray<real,1,1> quad_wts_phys;
      get_0form_quad_pts_wts(i, j, quad_pts_phys, quad_wts_phys);
      field.data(ndof, k+ks, j+js, i+is, n) = initial_value_function(quad_pts_phys(0).x, quad_pts_phys(0).y) * quad_wts_phys(0);
  });
}

  YAKL_INLINE void Geometry::set_2form_values(real (*initial_value_function)(real, real), Field &field, int ndof)
{



  int is = this->topology->is;
  int js = this->topology->js;
  int ks = this->topology->ks;

parallel_for( Bounds<4>(this->topology->nl, this->topology->n_cells_y, this->topology->n_cells_x, this->topology->nens) , YAKL_LAMBDA(int k, int j, int i, int n) { 
  SArray<coords<2>,2,ic_quad_pts_x,ic_quad_pts_y> quad_pts_phys;
  SArray<real,2,ic_quad_pts_x,ic_quad_pts_y> quad_wts_phys;
    get_2form_quad_pts_wts(i, j, quad_pts_phys, quad_wts_phys);
    real tempval = 0.0_fp;
    for (int nqx=0; nqx<ic_quad_pts_x; nqx++) {
      for (int nqy=0; nqy<ic_quad_pts_y; nqy++) {
          tempval = tempval + initial_value_function(quad_pts_phys(nqx, nqy).x, quad_pts_phys(nqx, nqy).y) * quad_wts_phys(nqx, nqy);
        }}
      field.data(ndof, k+ks, j+js, i+is, n) = tempval;
  });
}



  YAKL_INLINE void Geometry::set_1form_values(vec<2> (*initial_value_function)(real, real), Field &field, int ndof, LINE_INTEGRAL_TYPE line_type)
  {

      int is = this->topology->is;
      int js = this->topology->js;
      int ks = this->topology->ks;


      //twisted edges
      int yedge_offset = 0;
      int xedge_offset = field.ndofs;
      
      //compared to twisted edges in 2D, straight edges are stored x/y (V, -U) instead of y/x
      if (this->straight) {
        yedge_offset = field.ndofs;
        xedge_offset = 0;
      }

        parallel_for( Bounds<4>(this->topology->nl, this->topology->n_cells_y, this->topology->n_cells_x, this->topology->nens) , YAKL_LAMBDA(int k, int j, int i, int n) { 
          
          SArray<coords<2>,1,ic_quad_pts_x> x_quad_pts_phys;
          SArray<real,1,ic_quad_pts_x> x_quad_wts_phys;
          SArray<coords<2>,1,ic_quad_pts_y> y_quad_pts_phys;
          SArray<real,1,ic_quad_pts_y> y_quad_wts_phys;

          SArray<vec<2>,1,ic_quad_pts_x> x_line_vec;
          SArray<vec<2>,1,ic_quad_pts_y> y_line_vec;

          vec<2> initval;

        get_1form_quad_pts_wts(i, j, x_quad_pts_phys, x_quad_wts_phys, y_quad_pts_phys, y_quad_wts_phys);

        if (line_type == LINE_INTEGRAL_TYPE::TANGENT) {get_edge_tangents(i, j, x_line_vec, y_line_vec);}
        if (line_type == LINE_INTEGRAL_TYPE::NORMAL) {get_edge_normals(i, j, x_line_vec, y_line_vec);}

      // y-edge
        real tempvaly = 0.0_fp;
      for (int nqy=0; nqy<ic_quad_pts_y; nqy++) {
        initval = initial_value_function(y_quad_pts_phys(nqy).x, y_quad_pts_phys(nqy).y);
        tempvaly = tempvaly + (initval.u * y_line_vec(nqy).u + initval.v * y_line_vec(nqy).v) * y_quad_wts_phys(nqy);
      }
      field.data(ndof+yedge_offset, k+ks, j+js, i+is, n) = tempvaly;

            // x-edge
            real tempvalx = 0.0_fp;
            for (int nqx=0; nqx<ic_quad_pts_x; nqx++) {
              initval = initial_value_function(x_quad_pts_phys(nqx).x, x_quad_pts_phys(nqx).y);
              tempvalx = tempvalx + (initval.u * x_line_vec(nqx).u + initval.v * x_line_vec(nqx).v) * x_quad_wts_phys(nqx);
            }
            field.data(ndof+xedge_offset, k+ks, j+js, i+is, n) = tempvalx;
      });
  }

void UniformRectangularTwistedGeometry::initialize(const Topology &topo, const ModelParameters &params)
{

Geometry::initialize(topo);

this->Lx = params.xlen;
this->Ly = params.ylen;

this->xc = params.xc;
this->yc = params.yc;

this->dx = params.xlen/params.nx_glob;
this->dy = params.ylen/params.ny_glob;

this->straight = false;

this->is_initialized = true;
}

void UniformRectangularTwistedGeometry::printinfo()
{
  std::cout << "uniform rectangular geometry 2D info: twisted \n" << std::flush;
  std::cout << "Lx " << this->Lx << " Ly " << this->Ly << "\n" << std::flush;
  std::cout << "xc " << this->xc << " yc " << this->yc << "\n" << std::flush;
  std::cout << "dx " << this->dx << " dy " << this->dy << "\n" << std::flush;
}


real YAKL_INLINE UniformRectangularTwistedGeometry::get_area_lform(int l, int d, int k, int j, int i)
{
  if (l == 0) {return 1._fp;}
  if (l == 1) {
    if (d==0) {return this->dy;}
    if (d==1) {return this->dx;}
    }
  if (l == 2) {return this->dx * this->dy;}
}

void YAKL_INLINE UniformRectangularTwistedGeometry::get_0form_quad_pts_wts(int i, int j, SArray<coords<2>,1,1> &quad_pts_phys, SArray<real,1,1> &quad_wts_phys)
{
  quad_pts_phys(0).x = (i+this->topology->i_beg)*this->dx + this->xc - this->Lx * 0.5_fp;
  quad_pts_phys(0).y = (j+this->topology->j_beg)*this->dy + this->yc - this->Ly * 0.5_fp;
  quad_wts_phys(0) = 1._fp;
}

void YAKL_INLINE UniformRectangularTwistedGeometry::get_2form_quad_pts_wts(int i, int j, SArray<coords<2>,2,ic_quad_pts_x,ic_quad_pts_y> &quad_pts_phys, SArray<real,2,ic_quad_pts_x,ic_quad_pts_y> &quad_wts_phys)
{
  real ll_corner_x = (i+this->topology->i_beg)*this->dx + this->xc - this->Lx * 0.5_fp;
  real ll_corner_y = (j+this->topology->j_beg)*this->dy + this->yc - this->Ly * 0.5_fp;
  for (int nqx=0; nqx<ic_quad_pts_x; nqx++)
  for (int nqy=0; nqy<ic_quad_pts_y; nqy++)
  {{
    quad_pts_phys(nqx, nqy).x = ll_corner_x + this->x_quad_pts_ref(nqx) * this->dx;
    quad_pts_phys(nqx, nqy).y = ll_corner_y + this->y_quad_pts_ref(nqy) * this->dy;
    quad_wts_phys(nqx, nqy) = this->x_quad_wts_ref(nqx) * this->dx * this->y_quad_wts_ref(nqy) * this->dy;
  }}
}

void YAKL_INLINE UniformRectangularTwistedGeometry::get_1form_quad_pts_wts(int i, int j, SArray<coords<2>,1,ic_quad_pts_x> &x_quad_pts_phys, SArray<real,1,ic_quad_pts_x> &x_quad_wts_phys, SArray<coords<2>,1,ic_quad_pts_y> &y_quad_pts_phys, SArray<real,1,ic_quad_pts_y> &y_quad_wts_phys)
{
    real ll_corner_x = (i+this->topology->i_beg)*this->dx + this->xc - this->Lx * 0.5_fp;
    real ll_corner_y = (j+this->topology->j_beg)*this->dy + this->yc - this->Ly * 0.5_fp;
    for (int nqx=0; nqx<ic_quad_pts_x; nqx++) {
        x_quad_pts_phys(nqx).x = ll_corner_x + this->x_quad_pts_ref(nqx) * this->dx;
        x_quad_pts_phys(nqx).y = ll_corner_y;
        x_quad_wts_phys(nqx) = this->x_quad_wts_ref(nqx) * this->dx;
    }
    for (int nqy=0; nqy<ic_quad_pts_y; nqy++) {
        y_quad_pts_phys(nqy).x = ll_corner_x;
        y_quad_pts_phys(nqy).y = ll_corner_y + this->y_quad_pts_ref(nqy) * this->dy;
        y_quad_wts_phys(nqy) = this->y_quad_wts_ref(nqy) * this->dy;
    }
}


void UniformRectangularStraightGeometry::initialize(const Topology &topo, const ModelParameters &params)
{

Geometry::initialize(topo);

this->Lx = params.xlen;
this->Ly = params.ylen;

this->xc = params.xc;
this->yc = params.yc;

this->dx = params.xlen/params.nx_glob;
this->dy = params.ylen/params.ny_glob;

this->straight = true;

this->is_initialized = true;
}

 void UniformRectangularStraightGeometry::printinfo()
{
  std::cout << "uniform rectangular geometry 2D info: straight \n" << std::flush;
  std::cout << "Lx " << this->Lx << " Ly " << this->Ly << "\n" << std::flush;
  std::cout << "xc " << this->xc << " yc " << this->yc << "\n" << std::flush;
  std::cout << "dx " << this->dx << " dy " << this->dy << "\n" << std::flush;
}

real YAKL_INLINE UniformRectangularStraightGeometry::get_area_lform(int l, int d, int k, int j, int i)
{
  if (l == 0) {return 1._fp;}
  if (l == 1) {
    if (d==0) {return this->dx;}
    if (d==1) {return this->dy;}
    }
  if (l == 2) {return this->dx * this->dy;}
}

void YAKL_INLINE UniformRectangularStraightGeometry::get_0form_quad_pts_wts(int i, int j, SArray<coords<2>,1,1> &quad_pts_phys, SArray<real,1,1> &quad_wts_phys)
{
  quad_pts_phys(0).x = (i+this->topology->i_beg)*this->dx + this->xc - this->Lx* 0.5_fp + this->dx* 0.5_fp;
  quad_pts_phys(0).y = (j+this->topology->j_beg)*this->dy + this->yc - this->Ly* 0.5_fp + this->dy* 0.5_fp;
  quad_wts_phys(0) = 1.;
}

void YAKL_INLINE UniformRectangularStraightGeometry::get_2form_quad_pts_wts(int i, int j, SArray<coords<2>,2,ic_quad_pts_x,ic_quad_pts_y> &quad_pts_phys, SArray<real,2,ic_quad_pts_x,ic_quad_pts_y> &quad_wts_phys)
{
    real ll_corner_x = (i+this->topology->i_beg)*this->dx + this->xc - this->Lx* 0.5_fp - this->dx* 0.5_fp;
    real ll_corner_y = (j+this->topology->j_beg)*this->dy + this->yc - this->Ly* 0.5_fp - this->dy* 0.5_fp;
    for (int nqx=0; nqx<ic_quad_pts_x; nqx++)
    for (int nqy=0; nqy<ic_quad_pts_y; nqy++)
    {{
      quad_pts_phys(nqx, nqy).x = ll_corner_x + this->x_quad_pts_ref(nqx) * this->dx;
      quad_pts_phys(nqx, nqy).y = ll_corner_y + this->y_quad_pts_ref(nqy) * this->dy;
      quad_wts_phys(nqx, nqy) = this->x_quad_wts_ref(nqx) * this->dx * this->y_quad_wts_ref(nqy) * this->dy;
    }}
}

void YAKL_INLINE UniformRectangularStraightGeometry::get_1form_quad_pts_wts(int i, int j, SArray<coords<2>,1,ic_quad_pts_x> &x_quad_pts_phys, SArray<real,1,ic_quad_pts_x> &x_quad_wts_phys, SArray<coords<2>,1,ic_quad_pts_y> &y_quad_pts_phys, SArray<real,1,ic_quad_pts_y> &y_quad_wts_phys)
{
  real ll_corner_x = (i+this->topology->i_beg)*this->dx + this->xc - this->Lx* 0.5_fp - this->dx* 0.5_fp;
  real ll_corner_y = (j+this->topology->j_beg)*this->dy + this->yc - this->Ly* 0.5_fp - this->dy* 0.5_fp;
      for (int nqx=0; nqx<ic_quad_pts_x; nqx++) {
          x_quad_pts_phys(nqx).x = ll_corner_x + this->x_quad_pts_ref(nqx) * this->dx;
          x_quad_pts_phys(nqx).y = ll_corner_y + this->dy;
          x_quad_wts_phys(nqx) = this->x_quad_wts_ref(nqx) * this->dx;
      }
      for (int nqy=0; nqy<ic_quad_pts_y; nqy++) {
          y_quad_pts_phys(nqy).x = ll_corner_x + this->dx;
          y_quad_pts_phys(nqy).y = ll_corner_y + this->y_quad_pts_ref(nqy) * this->dy;
          y_quad_wts_phys(nqy) = this->y_quad_wts_ref(nqy) * this->dy;
      }
}


void YAKL_INLINE UniformRectangularTwistedGeometry::get_edge_tangents(int i, int j, SArray<vec<2>,1,ic_quad_pts_x> &x_tangents, SArray<vec<2>,1,ic_quad_pts_y> &y_tangents)
{
    for (int nqx=0; nqx<ic_quad_pts_x; nqx++)
    {
        x_tangents(nqx).u = -1._fp;
        x_tangents(nqx).v = 0._fp;
    }
    for (int nqy=0; nqy<ic_quad_pts_y; nqy++)
    {
      y_tangents(nqy).u = 0._fp;
      y_tangents(nqy).v = 1._fp;
    }
}


void YAKL_INLINE UniformRectangularTwistedGeometry::get_edge_normals(int i, int j, SArray<vec<2>,1,ic_quad_pts_x> &x_normals, SArray<vec<2>,1,ic_quad_pts_y> &y_normals)
{
    for (int nqx=0; nqx<ic_quad_pts_x; nqx++)
    {
        x_normals(nqx).u = 0._fp;
        x_normals(nqx).v = 1._fp;
    }
    for (int nqy=0; nqy<ic_quad_pts_y; nqy++)
    {
      y_normals(nqy).u = 1._fp;
      y_normals(nqy).v = 0._fp;
    }
}

// For straight edges in 2D, the twisted tangent is the straight normal, and the straight normal is the twisted tangent

void YAKL_INLINE UniformRectangularStraightGeometry::get_edge_tangents(int i, int j, SArray<vec<2>,1,ic_quad_pts_x> &x_tangents, SArray<vec<2>,1,ic_quad_pts_y> &y_tangents)
{
    for (int nqx=0; nqx<ic_quad_pts_x; nqx++)
    {
        x_tangents(nqx).u = 1._fp;
        x_tangents(nqx).v = 0._fp;
    }
    for (int nqy=0; nqy<ic_quad_pts_y; nqy++)
    {
      y_tangents(nqy).u = 0._fp;
      y_tangents(nqy).v = 1._fp;
    }
}


void YAKL_INLINE UniformRectangularStraightGeometry::get_edge_normals(int i, int j, SArray<vec<2>,1,ic_quad_pts_x> &x_normals, SArray<vec<2>,1,ic_quad_pts_y> &y_normals)
{
    for (int nqx=0; nqx<ic_quad_pts_x; nqx++)
    {
        x_normals(nqx).u = 0._fp;
        x_normals(nqx).v = 1._fp;
    }
    for (int nqy=0; nqy<ic_quad_pts_y; nqy++)
    {
      y_normals(nqy).u = -1._fp;
      y_normals(nqy).v = 0._fp;
    }
}

#endif







#ifdef _EXTRUDED   


// *********************   2D EXTRUDED  ******************************/
//HOW DO WE GENERALIZE THIS TO 2D+1D? ie 3D EXTRUDED?

//EVENTUALLY SWAP TO USE EXTRUDED COORDS ie x,"z"
//PART OF A BIGGER VARIABLE VERTICAL COORDINATE FIX IE REALLY SWAP EOMS TO A GENERIC "TERRAIN-FOLLOWING" EULERIAN VERTICAL COORDINATE...

class Geometry {

public:

  const Topology *topology;
  bool is_initialized;
  bool straight;
  
  Geometry();
  Geometry( const Geometry &geom) = delete;
  Geometry& operator=( const Geometry &geom) = delete;

  SArray<real,1,ic_quad_pts_x> x_quad_pts_ref;
  SArray<real,1,ic_quad_pts_x> x_quad_wts_ref;
  SArray<real,1,ic_quad_pts_z> v_quad_pts_ref;
  SArray<real,1,ic_quad_pts_z> v_quad_wts_ref;

  void initialize(const Topology &topo);

  virtual void YAKL_INLINE get_00form_quad_pts_wts(int i, int k, SArray<coordsext<2>,1,1> &quad_pts_phys, SArray<real,1,1> &quad_wts_phys) {};
  virtual void YAKL_INLINE get_11form_quad_pts_wts(int i, int k, SArray<coordsext<2>,2,ic_quad_pts_x,ic_quad_pts_z> &quad_pts_phys, SArray<real,2,ic_quad_pts_x,ic_quad_pts_z> &quad_wts_phys) {};

  virtual void YAKL_INLINE get_01form_quad_pts_wts(int i, int k, SArray<coordsext<2>,1,ic_quad_pts_z> &v_quad_pts_phys, SArray<real,1,ic_quad_pts_z> &v_quad_wts_phys) {};
  virtual void YAKL_INLINE get_01edge_tangents(int i, int k,  SArray<vecext<2>,1,ic_quad_pts_z> &v_tangents) {};
  virtual void YAKL_INLINE get_01edge_normals(int i, int k, SArray<vecext<2>,1,ic_quad_pts_z> &v_normals) {};
  
  virtual void YAKL_INLINE get_10form_quad_pts_wts(int i, int k, SArray<coordsext<2>,1,ic_quad_pts_x> &x_quad_pts_phys, SArray<real,1,ic_quad_pts_x> &x_quad_wts_phys) {};
  virtual void YAKL_INLINE get_10edge_tangents(int i, int k, SArray<vecext<2>,1,ic_quad_pts_x> &x_tangents) {};
  virtual void YAKL_INLINE get_10edge_normals(int i, int k, SArray<vecext<2>,1,ic_quad_pts_x> &x_normals) {};

  virtual real YAKL_INLINE get_area_00entity(int k, int j, int i) {};
  virtual real YAKL_INLINE get_area_11entity(int k, int j, int i) {};
  virtual real YAKL_INLINE get_area_10entity(int k, int j, int i) {};
  virtual real YAKL_INLINE get_area_01entity(int k, int j, int i) {};
  
  virtual real YAKL_INLINE get_area_lform(int l, int d, int k, int j, int i) {};

  YAKL_INLINE void set_00form_values(real (*initial_value_function)(real, real), Field &field, int ndof);
  YAKL_INLINE void set_11form_values(real (*initial_value_function)(real, real), Field &field, int ndof);
  YAKL_INLINE void set_10form_values(vecext<2> (*initial_value_function)(real, real), Field &field, int ndof, LINE_INTEGRAL_TYPE line_type);
  YAKL_INLINE void set_01form_values(vecext<2> (*initial_value_function)(real, real), Field &field, int ndof, LINE_INTEGRAL_TYPE line_type);

  real dx, dy, dz;
  real Lx, Ly, Lz;
  real xc, yc, zc;
  
};

class UniformRectangularTwistedExtrudedGeometry: public Geometry {


public:

  void initialize(const Topology &topo, const ModelParameters &params);

  void printinfo();

  void YAKL_INLINE get_00form_quad_pts_wts(int i, int k, SArray<coordsext<2>,1,1> &quad_pts_phys, SArray<real,1,1> &quad_wts_phys);
  void YAKL_INLINE get_11form_quad_pts_wts(int i, int k, SArray<coordsext<2>,2,ic_quad_pts_x,ic_quad_pts_z> &quad_pts_phys, SArray<real,2,ic_quad_pts_x,ic_quad_pts_z> &quad_wts_phys);

  void YAKL_INLINE get_01form_quad_pts_wts(int i, int k, SArray<coordsext<2>,1,ic_quad_pts_z> &v_quad_pts_phys, SArray<real,1,ic_quad_pts_z> &v_quad_wts_phys);
  void YAKL_INLINE get_01edge_tangents(int i, int k,  SArray<vecext<2>,1,ic_quad_pts_z> &v_tangents);
  void YAKL_INLINE get_01edge_normals(int i, int k, SArray<vecext<2>,1,ic_quad_pts_z> &v_normals);
  
  void YAKL_INLINE get_10form_quad_pts_wts(int i, int k, SArray<coordsext<2>,1,ic_quad_pts_x> &x_quad_pts_phys, SArray<real,1,ic_quad_pts_x> &x_quad_wts_phys);
  void YAKL_INLINE get_10edge_tangents(int i, int k, SArray<vecext<2>,1,ic_quad_pts_x> &x_tangents);
  void YAKL_INLINE get_10edge_normals(int i, int k, SArray<vecext<2>,1,ic_quad_pts_x> &x_normals);

  real YAKL_INLINE get_area_00entity(int k, int j, int i);
  real YAKL_INLINE get_area_11entity(int k, int j, int i);
  real YAKL_INLINE get_area_10entity(int k, int j, int i);
  real YAKL_INLINE get_area_01entity(int k, int j, int i);
  
  real YAKL_INLINE get_area_lform(int l, int d, int k, int j, int i) {};

};

class UniformRectangularStraightExtrudedGeometry: public Geometry {


public:

  void initialize(const Topology &topo, const ModelParameters &params);

  void printinfo();

  void YAKL_INLINE get_00form_quad_pts_wts(int i, int k, SArray<coordsext<2>,1,1> &quad_pts_phys, SArray<real,1,1> &quad_wts_phys);
  void YAKL_INLINE get_11form_quad_pts_wts(int i, int k, SArray<coordsext<2>,2,ic_quad_pts_x,ic_quad_pts_z> &quad_pts_phys, SArray<real,2,ic_quad_pts_x,ic_quad_pts_z> &quad_wts_phys);

  void YAKL_INLINE get_01form_quad_pts_wts(int i, int k, SArray<coordsext<2>,1,ic_quad_pts_z> &v_quad_pts_phys, SArray<real,1,ic_quad_pts_z> &v_quad_wts_phys);
  void YAKL_INLINE get_01edge_tangents(int i, int k,  SArray<vecext<2>,1,ic_quad_pts_z> &v_tangents);
  void YAKL_INLINE get_01edge_normals(int i, int k, SArray<vecext<2>,1,ic_quad_pts_z> &v_normals);
  
  void YAKL_INLINE get_10form_quad_pts_wts(int i, int k, SArray<coordsext<2>,1,ic_quad_pts_x> &x_quad_pts_phys, SArray<real,1,ic_quad_pts_x> &x_quad_wts_phys);
  void YAKL_INLINE get_10edge_tangents(int i, int k, SArray<vecext<2>,1,ic_quad_pts_x> &x_tangents);
  void YAKL_INLINE get_10edge_normals(int i, int k, SArray<vecext<2>,1,ic_quad_pts_x> &x_normals);

  real YAKL_INLINE get_area_00entity(int k, int j, int i);
  real YAKL_INLINE get_area_11entity(int k, int j, int i);
  real YAKL_INLINE get_area_10entity(int k, int j, int i);
  real YAKL_INLINE get_area_01entity(int k, int j, int i);
  real YAKL_INLINE get_area_lform(int l, int d, int k, int j, int i) {};
  
};




 Geometry::Geometry()
    {
      this->is_initialized = false;
      std::cout << "CREATED GEOMETRY\n";
    }

 void Geometry::initialize(const Topology &topo)
  {

  this->topology = &topo;
  set_ref_quad_pts_wts<real,ic_quad_pts_x>(this->x_quad_pts_ref, this->x_quad_wts_ref);
  set_ref_quad_pts_wts<real,ic_quad_pts_z>(this->v_quad_pts_ref, this->v_quad_wts_ref);
}

YAKL_INLINE void Geometry::set_00form_values(real (*initial_value_function)(real, real), Field &field, int ndof)
{
    int is = this->topology->is;
    int js = this->topology->js;
    int ks = this->topology->ks;
  
parallel_for( Bounds<4>(this->topology->ni, this->topology->n_cells_y, this->topology->n_cells_x, this->topology->nens) , YAKL_LAMBDA(int k, int j, int i, int n) { 
  SArray<coordsext<2>,1,1> quad_pts_phys;
  SArray<real,1,1> quad_wts_phys;
        get_00form_quad_pts_wts(i, k, quad_pts_phys, quad_wts_phys);
        field.data(ndof, k+ks, j+js, i+is, n) = initial_value_function(quad_pts_phys(0).x, quad_pts_phys(0).z) * quad_wts_phys(0);
    });
  
}
YAKL_INLINE void Geometry::set_11form_values(real (*initial_value_function)(real, real), Field &field, int ndof)
{

  
    int is = this->topology->is;
    int js = this->topology->js;
    int ks = this->topology->ks;

parallel_for( Bounds<4>(this->topology->nl, this->topology->n_cells_y, this->topology->n_cells_x, this->topology->nens) , YAKL_LAMBDA(int k, int j, int i, int n) { 
  
  SArray<coordsext<2>,2,ic_quad_pts_x,ic_quad_pts_z> quad_pts_phys;
  SArray<real,2,ic_quad_pts_x,ic_quad_pts_z> quad_wts_phys;
      real tempval = 0.0_fp;
      
      get_11form_quad_pts_wts(i, k, quad_pts_phys, quad_wts_phys);

      for (int nqx=0; nqx<ic_quad_pts_x; nqx++) {
        for (int nqz=0; nqz<ic_quad_pts_z; nqz++) {
            tempval = tempval + initial_value_function(quad_pts_phys(nqx, nqz).x, quad_pts_phys(nqx, nqz).z) * quad_wts_phys(nqx, nqz);
          }}

        field.data(ndof, k+ks, j+js, i+is, n) = tempval;
    });
}
  
YAKL_INLINE void Geometry::set_01form_values(vecext<2> (*initial_value_function)(real, real), Field &field, int ndof, LINE_INTEGRAL_TYPE line_type)
{

  
      int is = this->topology->is;
      int js = this->topology->js;
      int ks = this->topology->ks;
  

        parallel_for( Bounds<4>(this->topology->nl, this->topology->n_cells_y, this->topology->n_cells_x, this->topology->nens) , YAKL_LAMBDA(int k, int j, int i, int n) { 
          real tempval = 0.0_fp;
          SArray<coordsext<2>,1,ic_quad_pts_z> edge_quad_pts_phys;
          SArray<real,1,ic_quad_pts_z> edge_quad_wts_phys;
          SArray<vecext<2>,1,ic_quad_pts_z> edge_line_vec;
          vecext<2> initval;
          
        get_01form_quad_pts_wts(i, k, edge_quad_pts_phys, edge_quad_wts_phys);
  
        if (line_type == LINE_INTEGRAL_TYPE::TANGENT) {get_01edge_tangents(i, k, edge_line_vec);}
        if (line_type == LINE_INTEGRAL_TYPE::NORMAL) {get_01edge_normals(i, k, edge_line_vec);}
  
        tempval = 0.0_fp;
      for (int nqz=0; nqz<ic_quad_pts_z; nqz++) {
        initval = initial_value_function(edge_quad_pts_phys(nqz).x, edge_quad_pts_phys(nqz).z);
        tempval = tempval + (initval.u * edge_line_vec(nqz).u + initval.w * edge_line_vec(nqz).w) * edge_quad_wts_phys(nqz);
      }
      field.data(ndof, k+ks, j+js, i+is) = tempval;  
    });
}

  YAKL_INLINE void Geometry::set_10form_values(vecext<2> (*initial_value_function)(real, real), Field &field, int ndof, LINE_INTEGRAL_TYPE line_type)
{

  int is = this->topology->is;
  int js = this->topology->js;
  int ks = this->topology->ks;

    parallel_for( Bounds<4>(this->topology->ni, this->topology->n_cells_y, this->topology->n_cells_x, this->topology->nens) , YAKL_LAMBDA(int k, int j, int i, int n) { 
      SArray<coordsext<2>,1,ic_quad_pts_x> edge_quad_pts_phys;
      SArray<real,1,ic_quad_pts_x> edge_quad_wts_phys;
      SArray<vecext<2>,1,ic_quad_pts_x> edge_line_vec;
      vecext<2> initval;
    real tempval = 0.0_fp;
    
    get_10form_quad_pts_wts(i, k, edge_quad_pts_phys, edge_quad_wts_phys);

    if (line_type == LINE_INTEGRAL_TYPE::TANGENT) {get_10edge_tangents(i, k, edge_line_vec);}
    if (line_type == LINE_INTEGRAL_TYPE::NORMAL) {get_10edge_normals(i, k, edge_line_vec);}

    tempval = 0.0_fp;
  for (int nqx=0; nqx<ic_quad_pts_x; nqx++) {
    initval = initial_value_function(edge_quad_pts_phys(nqx).x, edge_quad_pts_phys(nqx).z);
    tempval = tempval + (initval.u * edge_line_vec(nqx).u + initval.w * edge_line_vec(nqx).w) * edge_quad_wts_phys(nqx);
  }
  field.data(ndof, k+ks, j+js, i+is, n) = tempval;
});
}

void UniformRectangularTwistedExtrudedGeometry::initialize(const Topology &topo, const ModelParameters &params)
{

Geometry::initialize(topo);

this->Lx = params.xlen;
this->Ly = 1._fp;
this->Lz = params.zlen;

this->xc = params.xc;
this->yc = 0.5_fp;
this->zc = params.zlen * 0.5_fp;

this->dx = params.xlen/params.nx_glob;
this->dy = 1._fp;
this->dz = params.zlen/params.nz;

this->straight = false;

this->is_initialized = true;
}

void UniformRectangularTwistedExtrudedGeometry::printinfo()
{
  std::cout << "uniform rectangular geometry 2D info: twisted \n" << std::flush;
  std::cout << "Lx " << this->Lx << " Ly " << this->Ly << " Lz " << this->Lz << "\n" << std::flush;
  std::cout << "xc " << this->xc << " yc " << this->yc << " zc " << this->zc << "\n" << std::flush;
  std::cout << "dx " << this->dx << " dy " << this->dy << " dz " << this->dz << "\n" << std::flush;
}

 real YAKL_INLINE UniformRectangularTwistedExtrudedGeometry::get_area_00entity(int k, int j, int i)
{ return 1.; }

real YAKL_INLINE UniformRectangularTwistedExtrudedGeometry::get_area_01entity(int k, int j, int i)
{ 
  if (k==this->topology->ks || k == this->topology->nl-1+this->topology->ks) //corrects for half-cells at top and bottom
  { return this->dz * 0.5_fp;}
  return this->dz;  
}

real YAKL_INLINE UniformRectangularTwistedExtrudedGeometry::get_area_10entity(int k, int j, int i)
{ return this->dx; }

real YAKL_INLINE UniformRectangularTwistedExtrudedGeometry::get_area_11entity(int k, int j, int i)
{ 
  if (k==this->topology->ks || k == this->topology->nl-1+this->topology->ks) //corrects for half-cells at top and bottom
  { return this->dx * this->dz * 0.5_fp;}
  return this->dx * this->dz; 
}


void YAKL_INLINE UniformRectangularTwistedExtrudedGeometry::get_00form_quad_pts_wts(int i, int k, SArray<coordsext<2>,1,1> &quad_pts_phys, SArray<real,1,1> &quad_wts_phys)
{
  quad_pts_phys(0).x = (i+this->topology->i_beg)*this->dx + this->xc - this->Lx * 0.5_fp;
  quad_pts_phys(0).z = k*this->dz - this->dz * 0.5_fp;
  
  //correct for half-cells
  if (k==0) {quad_pts_phys(0).z = 0.0_fp;}
  if (k==this->topology->ni-1) {quad_pts_phys(0).z = this->Lz;}
  
  quad_wts_phys(0) = 1.;
}

void YAKL_INLINE UniformRectangularTwistedExtrudedGeometry::get_11form_quad_pts_wts(int i, int k, SArray<coordsext<2>,2,ic_quad_pts_x,ic_quad_pts_z> &quad_pts_phys, SArray<real,2,ic_quad_pts_x,ic_quad_pts_z> &quad_wts_phys)
{
  real ll_corner_x = (i+this->topology->i_beg)*this->dx + this->xc - this->Lx * 0.5_fp;
  real ll_corner_z = k*this->dz - this->dz * 0.5_fp;
  real adj_dz = this->dz;
  
  //correct for half-cells
  if (k==0) {
    ll_corner_z = 0.0_fp;
    adj_dz = adj_dz * 0.5_fp;
  }
  if (k==this->topology->nl-1) {
    adj_dz = adj_dz * 0.5_fp;
  }

  for (int nqx=0; nqx<ic_quad_pts_x; nqx++)
  for (int nqz=0; nqz<ic_quad_pts_z; nqz++)
  {{
    quad_pts_phys(nqx, nqz).x = ll_corner_x + this->x_quad_pts_ref(nqx) * this->dx;
    quad_pts_phys(nqx, nqz).z = ll_corner_z + this->v_quad_pts_ref(nqz) * adj_dz;
    quad_wts_phys(nqx, nqz) = this->x_quad_wts_ref(nqx) * this->dx * this->v_quad_wts_ref(nqz) * adj_dz;
  }}
}

void YAKL_INLINE UniformRectangularTwistedExtrudedGeometry::get_01form_quad_pts_wts(int i, int k, SArray<coordsext<2>,1,ic_quad_pts_z> &v_quad_pts_phys, SArray<real,1,ic_quad_pts_z> &v_quad_wts_phys)
{
  real ll_corner_x = (i+this->topology->i_beg)*this->dx + this->xc - this->Lx * 0.5_fp;
  real ll_corner_z = k*this->dz - this->dz * 0.5_fp;
  real adj_dz = this->dz;
  
  //correct for half-cells
  if (k==0) {
    ll_corner_z = 0.0_fp;
    adj_dz = adj_dz * 0.5_fp;
  }
  if (k==this->topology->nl-1) {
    adj_dz = adj_dz * 0.5_fp;
  }
  
    for (int nqz=0; nqz<ic_quad_pts_z; nqz++) {
        v_quad_pts_phys(nqz).x = ll_corner_x;
        v_quad_pts_phys(nqz).z = ll_corner_z + this->v_quad_pts_ref(nqz) * adj_dz;
        v_quad_wts_phys(nqz) = this->v_quad_wts_ref(nqz) * adj_dz;
    }
}

void YAKL_INLINE UniformRectangularTwistedExtrudedGeometry::get_10form_quad_pts_wts(int i, int k, SArray<coordsext<2>,1,ic_quad_pts_x> &x_quad_pts_phys, SArray<real,1,ic_quad_pts_x> &x_quad_wts_phys)
{
  real ll_corner_x = (i+this->topology->i_beg)*this->dx + this->xc - this->Lx * 0.5_fp;
  real ll_corner_z = k*this->dz - this->dz * 0.5_fp;
  
  //correct for half-cells
  if (k==0) {ll_corner_z = 0.0_fp;}
  if (k==this->topology->ni-1) {ll_corner_z = this->Lz;}
  
    for (int nqx=0; nqx<ic_quad_pts_x; nqx++) {
        x_quad_pts_phys(nqx).x = ll_corner_x + this->x_quad_pts_ref(nqx) * this->dx;
        x_quad_pts_phys(nqx).z = ll_corner_z;
        x_quad_wts_phys(nqx) = this->x_quad_wts_ref(nqx) * this->dx;
    }
}


void UniformRectangularStraightExtrudedGeometry::initialize(const Topology &topo, const ModelParameters &params)
{

Geometry::initialize(topo);

this->Lx = params.xlen;
this->Ly = 0._fp;
this->Lz = params.zlen;

this->xc = params.xc;
this->yc = 0._fp;
this->zc = params.zlen/2.;

this->dx = params.xlen/params.nx_glob;
this->dy = 0._fp;
this->dz = params.zlen/params.nz;

this->straight = true;

this->is_initialized = true;
}

void UniformRectangularStraightExtrudedGeometry::printinfo()
{
  std::cout << "uniform rectangular geometry 2D info: straight \n" << std::flush;
  std::cout << "Lx " << this->Lx << " Ly " << this->Ly << " Lz " << this->Lz << "\n" << std::flush;
  std::cout << "xc " << this->xc << " yc " << this->yc << " zc " << this->zc << "\n" << std::flush;
  std::cout << "dx " << this->dx << " dy " << this->dy << " dz " << this->dz << "\n" << std::flush;
}



real YAKL_INLINE UniformRectangularStraightExtrudedGeometry::get_area_00entity(int k, int j, int i)
{return 1._fp;}

real YAKL_INLINE UniformRectangularStraightExtrudedGeometry::get_area_01entity(int k, int j, int i)
{return this->dz;}

real YAKL_INLINE UniformRectangularStraightExtrudedGeometry::get_area_10entity(int k, int j, int i)
{return this->dx;}

real YAKL_INLINE UniformRectangularStraightExtrudedGeometry::get_area_11entity(int k, int j, int i)
{return this->dx * this->dz;}

void YAKL_INLINE UniformRectangularStraightExtrudedGeometry::get_00form_quad_pts_wts(int i, int k, SArray<coordsext<2>,1,1> &quad_pts_phys, SArray<real,1,1> &quad_wts_phys)
{
  quad_pts_phys(0).x = (i+this->topology->i_beg)*this->dx + this->xc - this->Lx * 0.5_fp + this->dx * 0.5_fp;
  quad_pts_phys(0).z = (k)*this->dz;
  quad_wts_phys(0) = 1._fp;
}

void YAKL_INLINE UniformRectangularStraightExtrudedGeometry::get_11form_quad_pts_wts(int i, int k, SArray<coordsext<2>,2,ic_quad_pts_x,ic_quad_pts_z> &quad_pts_phys, SArray<real,2,ic_quad_pts_x,ic_quad_pts_z> &quad_wts_phys)
{
    real ll_corner_x = (i+this->topology->i_beg)*this->dx + this->xc - this->Lx * 0.5_fp - this->dx * 0.5_fp;
    real ll_corner_z = k*this->dz;
    for (int nqx=0; nqx<ic_quad_pts_x; nqx++)
    for (int nqz=0; nqz<ic_quad_pts_z; nqz++)
    {{
      quad_pts_phys(nqx, nqz).x = ll_corner_x + this->x_quad_pts_ref(nqx) * this->dx;
      quad_pts_phys(nqx, nqz).z = ll_corner_z + this->v_quad_pts_ref(nqz) * this->dz;
      quad_wts_phys(nqx, nqz) = this->x_quad_wts_ref(nqx) * this->dx * this->v_quad_wts_ref(nqz) * this->dz;
    }}
}

void YAKL_INLINE UniformRectangularStraightExtrudedGeometry::get_01form_quad_pts_wts(int i, int k, SArray<coordsext<2>,1,ic_quad_pts_z> &v_quad_pts_phys, SArray<real,1,ic_quad_pts_z> &v_quad_wts_phys)
{
  //Careful because 01 edge actually starts at same point as 00 vertex
  //This is slightly different than 2D layer situation...
    real ll_corner_x = (i+this->topology->i_beg)*this->dx + this->xc - this->Lx * 0.5_fp + this->dx * 0.5_fp;
    real ll_corner_z = k*this->dz;
    for (int nqz=0; nqz<ic_quad_pts_z; nqz++) {
        v_quad_pts_phys(nqz).x = ll_corner_x;
        v_quad_pts_phys(nqz).z = ll_corner_z + this->v_quad_pts_ref(nqz) * this->dz;
        v_quad_wts_phys(nqz) = this->v_quad_wts_ref(nqz) * this->dz;
    }
}

void YAKL_INLINE UniformRectangularStraightExtrudedGeometry::get_10form_quad_pts_wts(int i, int k, SArray<coordsext<2>,1,ic_quad_pts_x> &x_quad_pts_phys, SArray<real,1,ic_quad_pts_x> &x_quad_wts_phys)
{
  real ll_corner_x = (i+this->topology->i_beg)*this->dx + this->xc - this->Lx * 0.5_fp - this->dx * 0.5_fp;
  real ll_corner_z = k*this->dz;
  for (int nqx=0; nqx<ic_quad_pts_x; nqx++) {
      x_quad_pts_phys(nqx).x = ll_corner_x + this->x_quad_pts_ref(nqx) * this->dx;
      x_quad_pts_phys(nqx).z = ll_corner_z;
      x_quad_wts_phys(nqx) = this->x_quad_wts_ref(nqx) * this->dx;
  }
}
  
  
  
  
void YAKL_INLINE UniformRectangularTwistedExtrudedGeometry::get_10edge_tangents(int i, int k, SArray<vecext<2>,1,ic_quad_pts_x> &x_tangents)
{
    for (int nqx=0; nqx<ic_quad_pts_x; nqx++)
    {
        x_tangents(nqx).u = 1._fp;
        x_tangents(nqx).w = 0._fp;
    }
}

void YAKL_INLINE UniformRectangularTwistedExtrudedGeometry::get_01edge_tangents(int i, int k,SArray<vecext<2>,1,ic_quad_pts_z> &v_tangents)
{
  for (int nqz=0; nqz<ic_quad_pts_z; nqz++)
  {
    v_tangents(nqz).u = 0._fp;
    v_tangents(nqz).w = -1._fp;
  }
}

void YAKL_INLINE UniformRectangularTwistedExtrudedGeometry::get_10edge_normals(int i, int k, SArray<vecext<2>,1,ic_quad_pts_x> &x_normals)
{
    for (int nqx=0; nqx<ic_quad_pts_x; nqx++)
    {
        x_normals(nqx).u = 0._fp;
        x_normals(nqx).w = 1._fp;
    }
}

void YAKL_INLINE UniformRectangularTwistedExtrudedGeometry::get_01edge_normals(int i, int k, SArray<vecext<2>,1,ic_quad_pts_z> &v_normals)
{
  for (int nqz=0; nqz<ic_quad_pts_z; nqz++)
  {
    v_normals(nqz).u = 1._fp;
    v_normals(nqz).w = 0._fp;
  }
}
  
  
  
  
  
  
  
// For straight edges in 2D, the twisted tangent is the straight normal, and the straight normal is the twisted tangent

void YAKL_INLINE UniformRectangularStraightExtrudedGeometry::get_10edge_tangents(int i, int k, SArray<vecext<2>,1,ic_quad_pts_x> &x_tangents)
{
  for (int nqx=0; nqx<ic_quad_pts_x; nqx++)
  {
      x_tangents(nqx).u = 1._fp;
      x_tangents(nqx).w = 0._fp;
  }
}

void YAKL_INLINE UniformRectangularStraightExtrudedGeometry::get_01edge_tangents(int i, int k, SArray<vecext<2>,1,ic_quad_pts_z> &v_tangents)
{
  for (int nqz=0; nqz<ic_quad_pts_z; nqz++)
  {
    v_tangents(nqz).u = 0._fp;
    v_tangents(nqz).w = 1._fp;
  }
}

 void YAKL_INLINE UniformRectangularStraightExtrudedGeometry::get_10edge_normals(int i, int k, SArray<vecext<2>,1,ic_quad_pts_x> &x_normals)
{
  for (int nqx=0; nqx<ic_quad_pts_x; nqx++)
  {
      x_normals(nqx).u = 0._fp;
      x_normals(nqx).w = -1._fp;
  }
}

void YAKL_INLINE UniformRectangularStraightExtrudedGeometry::get_01edge_normals(int i, int k, SArray<vecext<2>,1,ic_quad_pts_z> &v_normals)
{
  for (int nqz=0; nqz<ic_quad_pts_z; nqz++)
  {
    v_normals(nqz).u = 1._fp;
    v_normals(nqz).w = 0._fp;
  }
}

#endif

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

