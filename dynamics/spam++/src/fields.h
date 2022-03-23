#ifndef _FIELDS_H_
#define _FIELDS_H_


#include "common.h"
#include "topology.h"
#include "exchange.h"
#include <string>


class Field {

public :

  real4d data;
  const Topology *topology;
  std::string name;
  int total_dofs; //total number of dofs in a "horiz" slice ie at each layer or interface
// basedof = 0,1,2 in 2D or 0,1 in 1D ie vertices,edges,cells or vertices,cells
// extdof = 0,1 (0 = interfaces, 1=layers)
  int basedof, extdof, ndofs;
  int _nz, _nloop, _nloop_halo;

  bool is_initialized;
  Field();
  //Field( const Field &f) = delete;
  //Field& operator=( const Field &f) = delete;
  void printinfo();
  void initialize(const Field &f, const std::string fieldName);
  void initialize(const Topology &topo, const std::string fieldName, int bdof, int edof, int nd);
  void copy(const Field & f);
  void waxpy(real alpha, const Field &x, const Field &y);
  void waxpbypcz(real alpha, real beta, real gamma, const Field &x, const Field &y, const Field &z);
  void zero();
  void zero(int ndof);
  void set(real val);
  void set(int ndof, real val);
  void set_top_bnd(real val);
  void set_top_bnd(int ndof, real val);
  void set_bottom_bnd(real val);
  void set_bottom_bnd(int ndof, real val);
  void set_bnd(real val);
  void set_bnd(int ndof, real val);
  // real sum();
  // real min();
  // real max();
  // real sum(int ndof);
  // real min(int ndof);
  // real max(int ndof);
};

    Field::Field()
    {
      this->is_initialized = false;
      std::cout << "CREATED FIELD\n";
    }

  void Field::printinfo()
  {
    std::cout << "field info " << this->name << "\n" << std::flush;
    std::cout <<  this->data << std::flush;
    std::cout << "total_dofs " << this->total_dofs << " basedof " << this->basedof << " extdof " << this->extdof << " ndofs " << this->ndofs  << "\n" << std::flush;
  }


  // creates a new Field f with same parameters as self, without copying data over
  void Field::initialize(const Field &f, const std::string fieldName)
  {
    this->initialize(*f.topology, fieldName, f.basedof, f.extdof, f.ndofs);
  }

  void Field::initialize(const Topology &topo, const std::string fieldName, int bdof, int edof, int nd)
  {

    this->topology = &topo;
    this->name = fieldName;
    this->basedof = bdof;
    this->extdof = edof;
    this->ndofs = nd;
    
    if (this->extdof == 0) {
      this->_nz = this->topology->ni;
      this->_nloop = this->topology->n_cells_interfaces;
      this->_nloop_halo = this->topology->n_cells_interfaces_with_halo;
    }
    if (this->extdof == 1) {
      this->_nz = this->topology->nl;
      this->_nloop = this->topology->n_cells_layers;
      this->_nloop_halo = this->topology->n_cells_layers_with_halo;
    }
    
    this->total_dofs = this->ndofs ;
    if (ndims == 2 && this->basedof==1) { this->total_dofs = 2*this->ndofs; } //2 edges per cell in 2D

      this->data = real4d(this->name.c_str(), this->total_dofs , this->_nz+2*this->topology->mirror_halo, this->topology->n_cells_y+2*this->topology->halosize_y, this->topology->n_cells_x+2*this->topology->halosize_x);


    this->zero();
    this->is_initialized = true;

  }

  void Field::set(real val)
  {
      parallel_for( Bounds<4>(this->total_dofs, this->_nz +2*this->topology->mirror_halo,this->topology->n_cells_y+2*this->topology->halosize_y, this->topology->n_cells_x+2*this->topology->halosize_x) , YAKL_LAMBDA(int ndof, int k, int j, int i) { 
        this->data(ndof, k, j, i) = val;
    });    
  }

  void Field::set(int ndof, real val)
  {
parallel_for( Bounds<3>(this->_nz +2*this->topology->mirror_halo,this->topology->n_cells_y+2*this->topology->halosize_y, this->topology->n_cells_x+2*this->topology->halosize_x) , YAKL_LAMBDA(int k, int j, int i) { 
        this->data(ndof, k, j, i) = val;
    });    
  }

  void Field::set_top_bnd(real val)
  {
      parallel_for( Bounds<3>(this->total_dofs,this->topology->n_cells_y+2*this->topology->halosize_y, this->topology->n_cells_x+2*this->topology->halosize_x) , YAKL_LAMBDA(int ndof, int j, int i) { 
        this->data(ndof, this->_nz + this->topology->mirror_halo, j, i) = val;
    });    
  }

  void Field::set_top_bnd(int ndof, real val)
  {
parallel_for( Bounds<2>(this->topology->n_cells_y+2*this->topology->halosize_y, this->topology->n_cells_x+2*this->topology->halosize_x) , YAKL_LAMBDA(int j, int i) { 
  this->data(ndof, this->_nz + this->topology->mirror_halo, j, i) = val;
    });    
  }

  void Field::set_bottom_bnd(real val)
  {
      parallel_for( Bounds<3>(this->total_dofs,this->topology->n_cells_y+2*this->topology->halosize_y, this->topology->n_cells_x+2*this->topology->halosize_x) , YAKL_LAMBDA(int ndof, int j, int i) { 
        this->data(ndof, this->topology->mirror_halo, j, i) = val;
    });    
  }

  void Field::set_bottom_bnd(int ndof, real val)
  {
parallel_for( Bounds<2>(this->topology->n_cells_y+2*this->topology->halosize_y, this->topology->n_cells_x+2*this->topology->halosize_x) , YAKL_LAMBDA(int j, int i) { 
  this->data(ndof, this->topology->mirror_halo, j, i) = val;
    });    
  }

  void Field::set_bnd(real val)
  {
    set_top_bnd(val);
    set_bottom_bnd(val);
  }
  
  void Field::set_bnd(int ndof, real val)
  {
    set_top_bnd(ndof, val);
    set_bottom_bnd(ndof, val);
  }
  
  void Field::zero()
  {
//    yakl::parallel_for("ZeroField", this->_nloop_halo, YAKL_LAMBDA (int iGlob) {
//      int k, j, i;
//      yakl::unpackIndices(iGlob, this->_nz +2*this->topology->mirror_halo, this->topology->n_cells_y+2*this->topology->halosize_y, this->topology->n_cells_x+2*this->topology->halosize_x, k, j, i);
parallel_for( Bounds<4>(this->total_dofs, this->_nz +2*this->topology->mirror_halo,this->topology->n_cells_y+2*this->topology->halosize_y, this->topology->n_cells_x+2*this->topology->halosize_x) , YAKL_LAMBDA(int ndof, int k, int j, int i) { 
    //  for (int ndof=0; ndof<this->total_dofs; ndof++) {
        this->data(ndof, k, j, i) = 0.0;
    //  }
    });
  }
  
  void Field::zero(int ndof)
  {
//    yakl::parallel_for("ZeroField", this->_nloop_halo, YAKL_LAMBDA (int iGlob) {
//      int k, j, i;
//      yakl::unpackIndices(iGlob, this->_nz +2*this->topology->mirror_halo, this->topology->n_cells_y+2*this->topology->halosize_y, this->topology->n_cells_x+2*this->topology->halosize_x, k, j, i);
parallel_for( Bounds<3>(this->_nz +2*this->topology->mirror_halo,this->topology->n_cells_y+2*this->topology->halosize_y, this->topology->n_cells_x+2*this->topology->halosize_x) , YAKL_LAMBDA(int k, int j, int i) { 
        this->data(ndof, k, j, i) = 0.0;
    });
  }
  
  // copies data from f into self
  void Field::copy(const Field & f)
  {

    int is = this->topology->is;
    int js = this->topology->js;
    int ks = this->topology->ks;
  //  yakl::parallel_for("CopyField", this->_nloop, YAKL_LAMBDA (int iGlob) {
  //    int k, j, i;
  //    yakl::unpackIndices(iGlob, this->_nz, this->topology->n_cells_y, this->topology->n_cells_x, k, j, i);
  parallel_for( Bounds<4>(this->total_dofs, this->_nz,this->topology->n_cells_y, this->topology->n_cells_x) , YAKL_LAMBDA(int ndof, int k, int j, int i) { 
  //    for (int ndof=0; ndof<this->total_dofs; ndof++) {
        this->data(ndof, k+ks, j+js, i+is) = f.data(ndof, k+ks, j+js, i+is);
    //  }
    });
  }

  // Computes w (self) = alpha x + y
  void Field::waxpy(real alpha, const Field &x, const Field &y)
  {

    int is = this->topology->is;
    int js = this->topology->js;
    int ks = this->topology->ks;
    //yakl::parallel_for("WAXPYField", this->_nloop, YAKL_LAMBDA (int iGlob) {
    //  int k, j, i;
    //  yakl::unpackIndices(iGlob, this->_nz, this->topology->n_cells_y, this->topology->n_cells_x, k, j, i);
    parallel_for( Bounds<4>(this->total_dofs, this->_nz,this->topology->n_cells_y, this->topology->n_cells_x) , YAKL_LAMBDA(int ndof, int k, int j, int i) { 
    //  for (int ndof=0; ndof<this->total_dofs; ndof++) {
        this->data(ndof, k+ks, j+js, i+is) = alpha * x.data(ndof, k+ks, j+js, i+is) + y.data(ndof, k+ks, j+js, i+is);
    //  }
    });
  }

  // Computes w (self) = alpha x + beta * y + gamma * z
  void Field::waxpbypcz(real alpha, real beta, real gamma, const Field &x, const Field &y, const Field &z)
  {

    int is = this->topology->is;
    int js = this->topology->js;
    int ks = this->topology->ks;
    //yakl::parallel_for("WAXPBYPCZField", this->_nloop, YAKL_LAMBDA (int iGlob) {
    //  int k, j, i;
    //  yakl::unpackIndices(iGlob, this->_nz, this->topology->n_cells_y, this->topology->n_cells_x, k, j, i);
    parallel_for( Bounds<4>(this->total_dofs, this->_nz,this->topology->n_cells_y, this->topology->n_cells_x) , YAKL_LAMBDA(int ndof, int k, int j, int i) { 
    //  for (int ndof=0; ndof<this->total_dofs; ndof++) {
        this->data(ndof, k+ks, j+js, i+is) = alpha * x.data(ndof, k+ks, j+js, i+is) + beta * y.data(ndof, k+ks, j+js, i+is) + gamma * z.data(ndof, k+ks, j+js, i+is);
    //  }
    });
  }
 // 
 //  // computes sum of field
 //  real Field::sum()
 //  {
 // 
 //    //int is = this->topology->is;
 //    //int js = this->topology->js;
 //    //int ks = this->topology->ks;
 //    //real sum = 0.0;
 // 
 //    //THIS SUMS OVER THE HALO ELEMENTS :(
 //    real sum = yakl::intrinsics::sum(this->data);
 // 
 //    //THIS CODE WAS ACTUALLY A REDUCTION, AND IT IS A MIRACLE IT EVER WORKED!
 //    //yakl::parallel_for("SumField", this->_nloop, YAKL_LAMBDA (int iGlob) {
 //    //  int k, j, i;
 //    //  yakl::unpackIndices(iGlob, this->_nz, this->topology->n_cells_y, this->topology->n_cells_x, k, j, i);
 //    //  for (int ndof=0; ndof<this->total_dofs; ndof++) {
 //    //    sum += this->data(ndof, k+ks, j+js, i+is);
 //    //  }
 //    //});
 //    return sum;
 //  }
 // 
 //  // computes sum of field
 //  real Field::sum(int ndof)
 //  {
 // 
 //    //THIS SUMS OVER THE HALO ELEMENTS, AND IS NOT SPECIFIC TO NDOF :(
 //    real sum = yakl::intrinsics::sum(this->data);
 // 
 //    //THIS CODE WAS ACTUALLY A REDUCTION, AND IT IS A MIRACLE IT EVER WORKED!
 //    // int is = this->topology->is;
 //    // int js = this->topology->js;
 //    // int ks = this->topology->ks;
 //    // real sum = 0.0;
 //    // yakl::parallel_for("SumField", this->_nloop, YAKL_LAMBDA (int iGlob) {
 //    //   int k, j, i;
 //    //   yakl::unpackIndices(iGlob, this->_nz, this->topology->n_cells_y, this->topology->n_cells_x, k, j, i);
 //    //     sum += this->data(ndof, k+ks, j+js, i+is);
 //    // });
 //    return sum;
 //  }
 // 
 //  // computes min of field
 // real Field::min()
 //  {
 //    //THIS SUMS OVER THE HALO ELEMENTS :(
 //    real min = yakl::intrinsics::minval(this->data);
 // 
 //    //THIS CODE WAS ACTUALLY A REDUCTION, AND IT IS A MIRACLE IT EVER WORKED!
 //    //int is = this->topology->is;
 //    //int js = this->topology->js;
 //    //int ks = this->topology->ks;
 //    // real min = this->data(0, 0, js, is);
 //    // yakl::parallel_for("MinField", this->_nloop, YAKL_LAMBDA (int iGlob) {
 //    //   int k, j, i;
 //    //   yakl::unpackIndices(iGlob, this->_nz, this->topology->n_cells_y, this->topology->n_cells_x, k, j, i);
 //    //   for (int ndof=0; ndof<this->total_dofs; ndof++) {
 //    //     min = mymin(this->data(ndof, k+ks, j+js, i+is), min);
 //    //   }
 //    // });
 //    return min;
 //  }
 // 
 //  // computes min of field
 // real Field::min(int ndof)
 //  {
 // 
 //    //THIS SUMS OVER THE HALO ELEMENTS, AND IS NOT SPECIFIC TO NDOF :(
 //    real min = yakl::intrinsics::minval(this->data);
 // 
 //    //THIS CODE WAS ACTUALLY A REDUCTION, AND IT IS A MIRACLE IT EVER WORKED!
 //    // int is = this->topology->is;
 //    // int js = this->topology->js;
 //    // int ks = this->topology->ks;
 //    // real min = this->data(ndof, 0, js, is);
 //    // yakl::parallel_for("MinField", this->_nloop, YAKL_LAMBDA (int iGlob) {
 //    //   int k, j, i;
 //    //   yakl::unpackIndices(iGlob, this->_nz, this->topology->n_cells_y, this->topology->n_cells_x, k, j, i);
 //    //     min = mymin(this->data(ndof, k+ks, j+js, i+is), min);
 //    // });
 //    return min;
 //  }
 // 
 //  // computes max of field
 //  real Field::max()
 //  {
 // 
 //    //THIS SUMS OVER THE HALO ELEMENTS :(
 //    real max = yakl::intrinsics::maxval(this->data);
 // 
 //    //THIS CODE WAS ACTUALLY A REDUCTION, AND IT IS A MIRACLE IT EVER WORKED!
 //    // int is = this->topology->is;
 //    // int js = this->topology->js;
 //    // int ks = this->topology->ks;
 //    // real max = this->data(0, 0, js, is);
 //    // yakl::parallel_for("MaxField", this->_nloop, YAKL_LAMBDA (int iGlob) {
 //    //   int k, j, i;
 //    //   yakl::unpackIndices(iGlob, this->_nz, this->topology->n_cells_y, this->topology->n_cells_x, k, j, i);
 //    //   for (int ndof=0; ndof<this->total_dofs; ndof++) {
 //    //     max = mymax(this->data(ndof, k+ks, j+js, i+is), max);
 //    //   }
 //    // });
 //    return max;
 //  }
 // 
 //  // computes max of field
 //  real Field::max(int ndof)
 //  {
 // 
 //    //THIS SUMS OVER THE HALO ELEMENTS, AND IS NOT SPECIFIC TO NDOF :(
 //    real max = yakl::intrinsics::maxval(this->data);
 // 
 //    //THIS CODE WAS ACTUALLY A REDUCTION, AND IT IS A MIRACLE IT EVER WORKED!
 //    // int is = this->topology->is;
 //    // int js = this->topology->js;
 //    // int ks = this->topology->ks;
 //    // real max = this->data(ndof, 0, js, is);
 //    // yakl::parallel_for("MaxField", this->_nloop, YAKL_LAMBDA (int iGlob) {
 //    //   int k, j, i;
 //    //   yakl::unpackIndices(iGlob, this->_nz, this->topology->n_cells_y, this->topology->n_cells_x, k, j, i);
 //    //     max = mymax(this->data(ndof, k+ks, j+js, i+is), max);
 //    // });
 //    return max;
 //  }

#endif
