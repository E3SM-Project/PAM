#pragma once

#include "common.h"
#include "topology.h"

struct FieldDescription {
  std::string name;
  Topology topology;
  int basedof, extdof, ndofs;
};

class Field {

public:
  real5d data;
  Topology topology;
  Exchange *m_exchange;
  std::string name;
  int total_dofs; // total number of dofs in a "horiz" slice ie at each layer or
                  // interface
  // basedof = 0,1,2 in 2D or 0,1 in 1D ie vertices,edges,cells or
  // vertices,cells extdof = 0,1 (0 = interfaces, 1=layers)
  int basedof, extdof, ndofs;
  int _nz; //, _nloop, _nloop_halo;

  bool is_initialized;
  Field();
  // Field( const Field &f) = delete;
  // Field& operator=( const Field &f) = delete;
  void printinfo();
  void initialize(const Field &f, const std::string fieldName);
  void initialize(const Topology &topo, Exchange *exchng,
                  const std::string fieldName, int bdof, int edof, int nd);
  void copy(const Field &f, const FIELDOP_EXTENT extent);
  void waxpy(real alpha, const Field &x, const Field &y,
             const FIELDOP_EXTENT extent);
  void waxpby(real alpha, real beta, const Field &x, const Field &y,
              const FIELDOP_EXTENT extent);
  void waxpbypcz(real alpha, real beta, real gamma, const Field &x,
                 const Field &y, const Field &z, const FIELDOP_EXTENT extent);
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
  void exchange();

  auto
  bounds(const FIELDOP_EXTENT extent = FIELDOP_EXTENT::WITHOUT_HALOS) const;
  auto
  offsets(const FIELDOP_EXTENT extent = FIELDOP_EXTENT::WITHOUT_HALOS) const;
  // real sum();
  // real min();
  // real max();
  // real sum(int ndof);
  // real min(int ndof);
  // real max(int ndof);
};

Field::Field() { this->is_initialized = false; }

void Field::printinfo() {
  std::cout << "field info " << this->name << "\n" << std::flush;
  std::cout << this->data << std::flush;
  std::cout << "total_dofs " << this->total_dofs << " basedof " << this->basedof
            << " extdof " << this->extdof << " ndofs " << this->ndofs << "\n"
            << std::flush;
}

// creates a new Field f with same parameters as self, without copying data over
void Field::initialize(const Field &f, const std::string fieldName) {
  this->initialize(f.topology, f.m_exchange, fieldName, f.basedof, f.extdof,
                   f.ndofs);
}

void Field::initialize(const Topology &topo, Exchange *exchng,
                       const std::string fieldName, int bdof, int edof,
                       int nd) {

  this->topology = topo;
  this->m_exchange = exchng;
  this->name = fieldName;
  this->basedof = bdof;
  this->extdof = edof;
  this->ndofs = nd;

  if (this->extdof == 0) {
    this->_nz = this->topology.ni;
    // this->_nloop = this->topology.n_cells_interfaces;
    // this->_nloop_halo = this->topology.n_cells_interfaces_with_halo;
  }
  if (this->extdof == 1) {
    this->_nz = this->topology.nl;
    // this->_nloop = this->topology.n_cells_layers;
    // this->_nloop_halo = this->topology.n_cells_layers_with_halo;
  }

  this->total_dofs = this->ndofs;
  if (ndims == 2 && this->basedof == 1) {
    this->total_dofs = 2 * this->ndofs;
  } // 2 edges per cell in 2D

  this->data = real5d(this->name.c_str(), this->total_dofs,
                      this->_nz + 2 * this->topology.mirror_halo,
                      this->topology.n_cells_y + 2 * this->topology.halosize_y,
                      this->topology.n_cells_x + 2 * this->topology.halosize_x,
                      this->topology.nens);

  this->zero();
  this->is_initialized = true;
}

auto Field::bounds(const FIELDOP_EXTENT extent) const {
  if (extent == FIELDOP_EXTENT::WITHOUT_HALOS) {
    return SimpleBounds<5>(this->total_dofs, this->_nz,
                           this->topology.n_cells_y, this->topology.n_cells_x,
                           this->topology.nens);
  } else if (extent == FIELDOP_EXTENT::WITH_HALOS) {
    return SimpleBounds<5>(
        this->total_dofs, this->_nz + 2 * this->topology.mirror_halo,
        this->topology.n_cells_y + 2 * this->topology.halosize_y,
        this->topology.n_cells_x + 2 * this->topology.halosize_x,
        this->topology.nens);
  }
}

auto Field::offsets(const FIELDOP_EXTENT extent) const {
  if (extent == FIELDOP_EXTENT::WITHOUT_HALOS) {
    return std::make_tuple(this->topology.is, this->topology.js,
                           this->topology.ks);
  } else if (extent == FIELDOP_EXTENT::WITH_HALOS) {
    return std::make_tuple(0, 0, 0);
  }
}

void Field::set(real val) {
  parallel_for(
      "Field set val",
      SimpleBounds<5>(this->total_dofs,
                      this->_nz + 2 * this->topology.mirror_halo,
                      this->topology.n_cells_y + 2 * this->topology.halosize_y,
                      this->topology.n_cells_x + 2 * this->topology.halosize_x,
                      this->topology.nens),
      YAKL_CLASS_LAMBDA(int ndof, int k, int j, int i, int n) {
        this->data(ndof, k, j, i, n) = val;
      });
}

void Field::set(int ndof, real val) {
  parallel_for(
      "Field set val ndof",
      SimpleBounds<4>(this->_nz + 2 * this->topology.mirror_halo,
                      this->topology.n_cells_y + 2 * this->topology.halosize_y,
                      this->topology.n_cells_x + 2 * this->topology.halosize_x,
                      this->topology.nens),
      YAKL_CLASS_LAMBDA(int k, int j, int i, int n) {
        this->data(ndof, k, j, i, n) = val;
      });
}

void Field::set_top_bnd(real val) {
  parallel_for(
      "Field set top bnd val",
      SimpleBounds<4>(this->total_dofs,
                      this->topology.n_cells_y + 2 * this->topology.halosize_y,
                      this->topology.n_cells_x + 2 * this->topology.halosize_x,
                      this->topology.nens),
      YAKL_CLASS_LAMBDA(int ndof, int j, int i, int n) {
        this->data(ndof, this->_nz + this->topology.mirror_halo, j, i, n) = val;
      });
}

void Field::set_top_bnd(int ndof, real val) {
  parallel_for(
      "Field set top bnd val ndof",
      SimpleBounds<3>(this->topology.n_cells_y + 2 * this->topology.halosize_y,
                      this->topology.n_cells_x + 2 * this->topology.halosize_x,
                      this->topology.nens),
      YAKL_CLASS_LAMBDA(int j, int i, int n) {
        this->data(ndof, this->_nz + this->topology.mirror_halo, j, i, n) = val;
      });
}

void Field::set_bottom_bnd(real val) {
  parallel_for(
      "Field set bottom val",
      SimpleBounds<4>(this->total_dofs,
                      this->topology.n_cells_y + 2 * this->topology.halosize_y,
                      this->topology.n_cells_x + 2 * this->topology.halosize_x,
                      this->topology.nens),
      YAKL_CLASS_LAMBDA(int ndof, int j, int i, int n) {
        this->data(ndof, this->topology.mirror_halo, j, i, n) = val;
      });
}

void Field::set_bottom_bnd(int ndof, real val) {
  parallel_for(
      "Field set bottom bnd val ndof",
      SimpleBounds<3>(this->topology.n_cells_y + 2 * this->topology.halosize_y,
                      this->topology.n_cells_x + 2 * this->topology.halosize_x,
                      this->topology.nens),
      YAKL_CLASS_LAMBDA(int j, int i, int n) {
        this->data(ndof, this->topology.mirror_halo, j, i, n) = val;
      });
}

void Field::set_bnd(real val) {
  set_top_bnd(val);
  set_bottom_bnd(val);
}

void Field::set_bnd(int ndof, real val) {
  set_top_bnd(ndof, val);
  set_bottom_bnd(ndof, val);
}

void Field::zero() {
  parallel_for(
      "Field zero",
      SimpleBounds<5>(this->total_dofs,
                      this->_nz + 2 * this->topology.mirror_halo,
                      this->topology.n_cells_y + 2 * this->topology.halosize_y,
                      this->topology.n_cells_x + 2 * this->topology.halosize_x,
                      this->topology.nens),
      YAKL_CLASS_LAMBDA(int ndof, int k, int j, int i, int n) {
        this->data(ndof, k, j, i, n) = 0.0;
      });
}

void Field::zero(int ndof) {
  parallel_for(
      "Field zero ndof",
      SimpleBounds<4>(this->_nz + 2 * this->topology.mirror_halo,
                      this->topology.n_cells_y + 2 * this->topology.halosize_y,
                      this->topology.n_cells_x + 2 * this->topology.halosize_x,
                      this->topology.nens),
      YAKL_CLASS_LAMBDA(int k, int j, int i, int n) {
        this->data(ndof, k, j, i, n) = 0.0;
      });
}

// copies data from f into self
void Field::copy(const Field &f, const FIELDOP_EXTENT extent) {

  const auto [is, js, ks] = this->offsets(extent);
  parallel_for(
      "Field copy", this->bounds(extent),
      YAKL_CLASS_LAMBDA(int ndof, int k, int j, int i, int n) {
        this->data(ndof, k + ks, j + js, i + is, n) =
            f.data(ndof, k + ks, j + js, i + is, n);
      });
}

// Computes w (self) = alpha x + y
void Field::waxpy(real alpha, const Field &x, const Field &y,
                  const FIELDOP_EXTENT extent) {

  const auto [is, js, ks] = this->offsets(extent);
  parallel_for(
      "Field waxpy", this->bounds(extent),
      YAKL_CLASS_LAMBDA(int ndof, int k, int j, int i, int n) {
        this->data(ndof, k + ks, j + js, i + is, n) =
            alpha * x.data(ndof, k + ks, j + js, i + is, n) +
            y.data(ndof, k + ks, j + js, i + is, n);
      });
}

// Computes w (self) = alpha x + beta * y
void Field::waxpby(real alpha, real beta, const Field &x, const Field &y,
                   const FIELDOP_EXTENT extent) {

  const auto [is, js, ks] = this->offsets(extent);
  parallel_for(
      "Field waxpby", this->bounds(extent),
      YAKL_CLASS_LAMBDA(int ndof, int k, int j, int i, int n) {
        this->data(ndof, k + ks, j + js, i + is, n) =
            alpha * x.data(ndof, k + ks, j + js, i + is, n) +
            beta * y.data(ndof, k + ks, j + js, i + is, n);
      });
}

// Computes w (self) = alpha x + beta * y + gamma * z
void Field::waxpbypcz(real alpha, real beta, real gamma, const Field &x,
                      const Field &y, const Field &z,
                      const FIELDOP_EXTENT extent) {

  const auto [is, js, ks] = this->offsets(extent);
  parallel_for(
      "Field waxpbypcz", this->bounds(extent),
      YAKL_CLASS_LAMBDA(int ndof, int k, int j, int i, int n) {
        this->data(ndof, k + ks, j + js, i + is, n) =
            alpha * x.data(ndof, k + ks, j + js, i + is, n) +
            beta * y.data(ndof, k + ks, j + js, i + is, n) +
            gamma * z.data(ndof, k + ks, j + js, i + is, n);
      });
}

void Field::exchange() { m_exchange->exchange_data(data); }
