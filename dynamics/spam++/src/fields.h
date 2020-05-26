#ifndef _FIELDS_H_
#define _FIELDS_H_


#include "common.h"
#include "topology.h"
#include "exchange.h"
#include <string>


class Field {

public :

  realArr data;
  const Topology *topology;
  std::string name;
  int total_dofs;
  int ndof0, ndof1, ndof2, ndof3;

  bool is_initialized;
  Field();
  Field( const Field &f) = delete;
  Field& operator=( const Field &f) = delete;
  int YAKL_INLINE get_offset(int ndof);
  void printinfo();
  void initialize(const Field &f, const std::string fieldName);
  void initialize(const Topology &topo, const std::string fieldName, int nd0, int nd1, int nd2, int nd3);
  void copy(const Field & f);
  void waxpy(real alpha, const Field &x, const Field &y);
  void waxpbypcz(real alpha, real beta, real gamma, const Field &x, const Field &y, const Field &z);
  real sum();
  real min();
  real max();

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
    std::cout << "total_dofs " << this->total_dofs << " ndof0 " << this->ndof0 << " ndof1 " << this->ndof1 << " ndof2 " << this->ndof2 << " ndof3 " << this->ndof3 << "\n" << std::flush;
  }



  // creates a new Field f with same parameters as self, without copying data over
  void Field::initialize(const Field &f, const std::string fieldName)
  {
    this->initialize(*f.topology, fieldName, f.ndof0, f.ndof1, f.ndof2, f.ndof3);
  }

  void Field::initialize(const Topology &topo, const std::string fieldName, int nd0, int nd1, int nd2, int nd3)
  {

    this->topology = &topo;
    this->name = fieldName;
    this->ndof0 = nd0;
    this->ndof1 = nd1;
    this->ndof2 = nd2;
    this->ndof3 = nd3;

    if (ndims == 1) { this->total_dofs = this->ndof0 + this->ndof1; }
    if (ndims == 2) { this->total_dofs = this->ndof0 + 2*this->ndof1 + this->ndof2; }
    if (ndims == 3) { this->total_dofs = this->ndof0 + 3*this->ndof1 + 3*this->ndof2 + this->ndof3; }

    this->data = realArr(this->name.c_str(), this->total_dofs , this->topology->n_cells_z+2*this->topology->halosize_z, this->topology->n_cells_y+2*this->topology->halosize_y, this->topology->n_cells_x+2*this->topology->halosize_x);


    yakl::parallel_for("ZeroField", this->topology->n_cells_with_halo, YAKL_LAMBDA (int iGlob) {
      int k, j, i;
      yakl::unpackIndices(iGlob, this->topology->n_cells_z+2*this->topology->halosize_z, this->topology->n_cells_y+2*this->topology->halosize_y, this->topology->n_cells_x+2*this->topology->halosize_x, k, j, i);
      for (int ndof=0; ndof<this->total_dofs; ndof++) {
        this->data(ndof, k, j, i) = 0.0;
      }
    });

    this->is_initialized = true;

  }

  int YAKL_INLINE Field::get_offset(int ndof)
  {
      if (ndof == 0)                  { return 0; }
      if (ndof == 1)                  { return this->ndof0; }
      if ((ndof == 2) && (ndims == 2)) { return this->ndof0 + 2*this->ndof1; }
      if ((ndof == 2) && (ndims == 3)) { return this->ndof0 + 3*this->ndof1; }
      if ((ndof == 3) && (ndims == 3)) { return this->ndof0 + 3*this->ndof1 + 3*this->ndof2; }
    }


  // copies data from f into self
  void Field::copy(const Field & f)
  {

    int is = this->topology->is;
    int js = this->topology->js;
    int ks = this->topology->ks;
    yakl::parallel_for("CopyField", this->topology->n_cells, YAKL_LAMBDA (int iGlob) {
      int k, j, i;
      yakl::unpackIndices(iGlob, this->topology->n_cells_z, this->topology->n_cells_y, this->topology->n_cells_x, k, j, i);
      for (int ndof=0; ndof<this->total_dofs; ndof++) {
        this->data(ndof, k+ks, j+js, i+is) = f.data(ndof, k+ks, j+js, i+is);
      }
    });
  }

  // Computes w (self) = alpha x + y
  void Field::waxpy(real alpha, const Field &x, const Field &y)
  {

    int is = this->topology->is;
    int js = this->topology->js;
    int ks = this->topology->ks;
    yakl::parallel_for("WAXPYField", this->topology->n_cells, YAKL_LAMBDA (int iGlob) {
      int k, j, i;
      yakl::unpackIndices(iGlob, this->topology->n_cells_z, this->topology->n_cells_y, this->topology->n_cells_x, k, j, i);
      for (int ndof=0; ndof<this->total_dofs; ndof++) {
        this->data(ndof, k+ks, j+js, i+is) = alpha * x.data(ndof, k+ks, j+js, i+is) + y.data(ndof, k+ks, j+js, i+is);
      }
    });
  }

  // Computes w (self) = alpha x + beta * y + gamma * z
  void Field::waxpbypcz(real alpha, real beta, real gamma, const Field &x, const Field &y, const Field &z)
  {

    int is = this->topology->is;
    int js = this->topology->js;
    int ks = this->topology->ks;
    yakl::parallel_for("WAXPBYPCZField", this->topology->n_cells, YAKL_LAMBDA (int iGlob) {
      int k, j, i;
      yakl::unpackIndices(iGlob, this->topology->n_cells_z, this->topology->n_cells_y, this->topology->n_cells_x, k, j, i);
      for (int ndof=0; ndof<this->total_dofs; ndof++) {
        this->data(ndof, k+ks, j+js, i+is) = alpha * x.data(ndof, k+ks, j+js, i+is) + beta * y.data(ndof, k+ks, j+js, i+is) + gamma * z.data(ndof, k+ks, j+js, i+is);
      }
    });
  }

  // computes sum of field
  real Field::sum()
  {

    int is = this->topology->is;
    int js = this->topology->js;
    int ks = this->topology->ks;
    real sum = 0.0;
    yakl::parallel_for("SumField", this->topology->n_cells, YAKL_LAMBDA (int iGlob) {
      int k, j, i;
      yakl::unpackIndices(iGlob, this->topology->n_cells_z, this->topology->n_cells_y, this->topology->n_cells_x, k, j, i);
      for (int ndof=0; ndof<this->total_dofs; ndof++) {
        sum += this->data(ndof, k+ks, j+js, i+is);
      }
    });
    return sum;
  }

  // computes min of field
 real Field::min()
  {

    int is = this->topology->is;
    int js = this->topology->js;
    int ks = this->topology->ks;
    real min = this->data(0, ks, js, is);
    yakl::parallel_for("MinField", this->topology->n_cells, YAKL_LAMBDA (int iGlob) {
      int k, j, i;
      yakl::unpackIndices(iGlob, this->topology->n_cells_z, this->topology->n_cells_y, this->topology->n_cells_x, k, j, i);
      for (int ndof=0; ndof<this->total_dofs; ndof++) {
        min = mymin(this->data(ndof, k+ks, j+js, i+is), min);
      }
    });
    return min;
  }

  // computes max of field
  real Field::max()
  {

    int is = this->topology->is;
    int js = this->topology->js;
    int ks = this->topology->ks;
    real max = this->data(0, ks, js, is);
    yakl::parallel_for("MaxField", this->topology->n_cells, YAKL_LAMBDA (int iGlob) {
      int k, j, i;
      yakl::unpackIndices(iGlob, this->topology->n_cells_z, this->topology->n_cells_y, this->topology->n_cells_x, k, j, i);
      for (int ndof=0; ndof<this->total_dofs; ndof++) {
        max = mymax(this->data(ndof, k+ks, j+js, i+is), max);
      }
    });
    return max;
  }

#endif
