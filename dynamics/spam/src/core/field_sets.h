#pragma once

#include "common.h"
#include "exchange.h"
#include "fields.h"
#include "topology.h"
#include <initializer_list>

template <uint num_fields> class ExchangeSet {

public:
  std::array<Exchange, num_fields> exchanges_arr;
  bool is_initialized;

  ExchangeSet();
  ExchangeSet(const ExchangeSet<num_fields> &exch) = delete;
  ExchangeSet &operator=(const ExchangeSet<num_fields> &exch) = delete;
  void initialize(const std::array<FieldDescription, num_fields> &desc_arr);
  void initialize(const ExchangeSet<num_fields> &es);
  void printinfo();
};

template <uint num_fields> class FieldSet {

public:
  std::array<Field, num_fields> fields_arr;
  std::string baseName;
  bool is_initialized;

  FieldSet();
  // FieldSet( const FieldSet<num_fields> &vs) = delete;
  // FieldSet& operator=( const FieldSet<num_fields> &vs) = delete;
  void printinfo();
  void initialize(const std::string name,
                  const std::array<FieldDescription, num_fields> &desc_arr,
                  ExchangeSet<num_fields> &exchange_set);
  void initialize(const FieldSet<num_fields> &vs, const std::string name);
  void copy(const FieldSet<num_fields> &vs);
  void waxpy(real alpha, const FieldSet<num_fields> &x,
             const FieldSet<num_fields> &y);
  void waxpby(real alpha, real beta, const FieldSet<num_fields> &x,
              const FieldSet<num_fields> &y);
  void waxpbypcz(real alpha, real beta, real gamma,
                 const FieldSet<num_fields> &x, const FieldSet<num_fields> &y,
                 const FieldSet<num_fields> &z);
  void exchange();

  //void exchange(const std::initializer_list<int> &indices);
  template <int num_exchanges>
  void exchange(const int(&indices)[num_exchanges]);
};

template <uint num_fields> ExchangeSet<num_fields>::ExchangeSet() {
  this->is_initialized = false;
}

template <uint num_fields>
void ExchangeSet<num_fields>::initialize(
    const std::array<FieldDescription, num_fields> &desc_arr) {
  for (int i = 0; i < num_fields; i++) {
    this->exchanges_arr[i].initialize(desc_arr[i].topology, desc_arr[i].basedof,
                                      desc_arr[i].extdof, desc_arr[i].ndofs);
  }
  this->is_initialized = true;
}

template <uint num_fields> void ExchangeSet<num_fields>::printinfo() {
  std::cout << "exchange set info\n" << std::flush;
  for (int i = 0; i < num_fields; i++) {
    this->exchanges_arr[i].printinfo();
  }
}

// EVENTUALLY THIS SHOULD BE MORE CLEVER IE PACK ALL THE FIELDS AT ONCE, ETC.
// template <uint num_fields>
// void ExchangeSet<num_fields>::exchange_variable_set(FieldSet<num_fields> &vs)
// {
//  for (int i = 0; i < num_fields; i++) {
//    this->exchanges_arr[i].exchange_field(vs.fields_arr[i].data);
//  }
//}

template <uint num_fields>
void ExchangeSet<num_fields>::initialize(const ExchangeSet<num_fields> &es) {
  for (int i = 0; i < num_fields; i++) {
    this->exchanges_arr[i].initialize(es.exchanges_arr[i]);
  }
  this->is_initialized = true;
}

template <uint num_fields> FieldSet<num_fields>::FieldSet() {
  this->is_initialized = false;
}

template <uint num_fields> void FieldSet<num_fields>::printinfo() {
  std::cout << "field set info " << this->baseName << "\n" << std::flush;
  for (int i = 0; i < num_fields; i++) {
    this->fields_arr[i].printinfo();
  }
}

template <uint num_fields>
void FieldSet<num_fields>::initialize(
    const std::string name,
    const std::array<FieldDescription, num_fields> &desc_arr,
    ExchangeSet<num_fields> &exchange_set) {
  this->baseName = name;
  for (int i = 0; i < num_fields; i++) {
    this->fields_arr[i].initialize(
        desc_arr[i].topology, &exchange_set.exchanges_arr[i], desc_arr[i].name,
        desc_arr[i].basedof, desc_arr[i].extdof, desc_arr[i].ndofs);
  }
  this->is_initialized = true;
}

// creates a new FieldSet that has the same fields as vs
template <uint num_fields>
void FieldSet<num_fields>::initialize(
    const FieldSet<num_fields> &vs,
    const std::string name) //,  EXCHANGE_TYPE exch_type)
{
  this->baseName = name;
  for (int i = 0; i < num_fields; i++) {
    this->fields_arr[i].initialize(vs.fields_arr[i], vs.fields_arr[i].name);
  }
  this->is_initialized = true;
}

// copies data from vs into self
template <uint num_fields>
void FieldSet<num_fields>::copy(const FieldSet<num_fields> &vs) {
  for (int i = 0; i < num_fields; i++) {
    this->fields_arr[i].copy(vs.fields_arr[i]);
  }
}

// Computes w (self) = alpha x + y
template <uint num_fields>
void FieldSet<num_fields>::waxpy(real alpha, const FieldSet<num_fields> &x,
                                 const FieldSet<num_fields> &y) {
  //for (int i = 0; i < num_fields; i++) {
  //  this->fields_arr[i].waxpy(alpha, x.fields_arr[i], y.fields_arr[i]);
  //}

  const auto topology = fields_arr[0].topology;
  const int is = topology.is;
  const int js = topology.js;
  const int ks = topology.ks;

  SArray<real5d, 1, num_fields> this_data;
  SArray<real5d, 1, num_fields> x_data;
  SArray<real5d, 1, num_fields> y_data;
  SArray<int, 1, num_fields> ndofs;

  int max_nz = 0;
  for (int i = 0; i < num_fields; ++i) {
    ndofs(i) = this->fields_arr[i].total_dofs;
    max_nz = std::max(max_nz, fields_arr[i]._nz);
    this_data(i) = this->fields_arr[i].data;
    x_data(i) = x.fields_arr[i].data;
    y_data(i) = y.fields_arr[i].data;
  }

  parallel_for(
      "Field Set waxpy",
      SimpleBounds<4>(max_nz, topology.n_cells_y,
                      topology.n_cells_x, topology.nens),
      YAKL_LAMBDA(int k, int j, int i, int n) {
        for (int f = 0; f < num_fields; ++f) {
          for (int l = 0; l < ndofs(f); ++l) {
            this_data(f)(l, k + ks, j + js, i + is, n) =
                alpha * x_data(f)(l, k + ks, j + js, i + is, n) + y_data(f)(l, k + ks, j + js, i + is, n);
          }
        }
      });
}

// Computes w (self) = alpha x + beta y
template <uint num_fields>
void FieldSet<num_fields>::waxpby(real alpha, real beta,
                                  const FieldSet<num_fields> &x,
                                  const FieldSet<num_fields> &y) {
  //for (int i = 0; i < num_fields; i++) {
  //  this->fields_arr[i].waxpby(alpha, beta, x.fields_arr[i], y.fields_arr[i]);
  //}

  const auto topology = fields_arr[0].topology;
  const int is = topology.is;
  const int js = topology.js;
  const int ks = topology.ks;

  SArray<real5d, 1, num_fields> this_data;
  SArray<real5d, 1, num_fields> x_data;
  SArray<real5d, 1, num_fields> y_data;
  SArray<int, 1, num_fields> ndofs;

  int max_nz = 0;
  for (int i = 0; i < num_fields; ++i) {
    ndofs(i) = this->fields_arr[i].total_dofs;
    max_nz = std::max(max_nz, fields_arr[i]._nz);
    this_data(i) = this->fields_arr[i].data;
    x_data(i) = x.fields_arr[i].data;
    y_data(i) = y.fields_arr[i].data;
  }

  parallel_for(
      "Field Set waxpby",
      SimpleBounds<4>(max_nz, topology.n_cells_y,
                      topology.n_cells_x, topology.nens),
      YAKL_LAMBDA(int k, int j, int i, int n) {
        for (int f = 0; f < num_fields; ++f) {
          for (int l = 0; l < ndofs(f); ++l) {
            this_data(f)(l, k + ks, j + js, i + is, n) =
                alpha * x_data(f)(l, k + ks, j + js, i + is, n) + beta * y_data(f)(l, k + ks, j + js, i + is, n);
          }
        }
      });
}

template <uint num_fields>
void FieldSet<num_fields>::waxpbypcz(real alpha, real beta, real gamma,
                                     const FieldSet<num_fields> &x,
                                     const FieldSet<num_fields> &y,
                                     const FieldSet<num_fields> &z) {
  //for (int i = 0; i < num_fields; i++) {
  //  this->fields_arr[i].waxpbypcz(alpha, beta, gamma, x.fields_arr[i],
  //                                y.fields_arr[i], z.fields_arr[i]);
  //}

  const auto topology = fields_arr[0].topology;
  const int is = topology.is;
  const int js = topology.js;
  const int ks = topology.ks;

  SArray<real5d, 1, num_fields> this_data;
  SArray<real5d, 1, num_fields> x_data;
  SArray<real5d, 1, num_fields> y_data;
  SArray<real5d, 1, num_fields> z_data;
  SArray<int, 1, num_fields> ndofs;

  int max_nz = 0;
  for (int i = 0; i < num_fields; ++i) {
    ndofs(i) = this->fields_arr[i].total_dofs;
    max_nz = std::max(max_nz, fields_arr[i]._nz);
    this_data(i) = this->fields_arr[i].data;
    x_data(i) = x.fields_arr[i].data;
    y_data(i) = y.fields_arr[i].data;
    z_data(i) = z.fields_arr[i].data;
  }

  parallel_for(
      "Field Set waxpbypcz",
      SimpleBounds<4>(max_nz, topology.n_cells_y,
                      topology.n_cells_x, topology.nens),
      YAKL_LAMBDA(int k, int j, int i, int n) {
        for (int f = 0; f < num_fields; ++f) {
          for (int l = 0; l < ndofs(f); ++l) {
            this_data(f)(l, k + ks, j + js, i + is, n) =
                alpha * x_data(f)(l, k + ks, j + js, i + is, n) +
                beta * y_data(f)(l, k + ks, j + js, i + is, n) +
                gamma * z_data(f)(l, k + ks, j + js, i + is, n);
          }
        }
      });
}

template <uint num_fields> void FieldSet<num_fields>::exchange() {
  //for (auto &f : fields_arr) {
  //  f.exchange();
  //}
  
  const auto topology = fields_arr[0].topology;
  const int is = topology.is;
  const int js = topology.js;
  const int ks = topology.ks;
  
  const int halosize_x = topology.halosize_x;
  const int n_cells_x = topology.n_cells_x;
  const int n_cells_y = topology.n_cells_y;

  SArray<real5d, 1, num_fields> this_data;
  SArray<int, 1, num_fields> ndofs;
  SArray<int, 1, num_fields> nzs;
  SArray<int, 1, num_fields> koffset;

  int max_nz = 0;
  for (int i = 0; i < num_fields; ++i) {
    ndofs(i) = this->fields_arr[i].total_dofs;
    max_nz = std::max(max_nz, fields_arr[i]._nz);
    nzs(i) = fields_arr[i]._nz;
    koffset(i) = fields_arr[i].extdof == 1 ? 0 : 1;
    this_data(i) = this->fields_arr[i].data;
  }

  parallel_for(
      "Field Set exchange x",
      SimpleBounds<4>(topology.halosize_x, max_nz,
                      topology.n_cells_y, topology.nens),
      YAKL_LAMBDA(int ii, int k, int j, int n) {
        for (int f = 0; f < num_fields; ++f) {
          for (int l = 0; l < ndofs(f); ++l) {
            this_data(f)(l, k + ks, j + js, ii + is - halosize_x, n) =
                this_data(f)(l, k + ks, j + js, ii + is + n_cells_x - halosize_x, n);
            this_data(f)(l, k + ks, j + js, ii + is + n_cells_x, n) =
                this_data(f)(l, k + ks, j + js, ii + is, n);
          }
        }
   });
  
  //for (auto &f : fields_arr) {
  //  f.exchange_direct();
  //}
  //for (auto &f : fields_arr) {
  //  f.exchange_mirror();
  //}

    parallel_for(
        "Field Set exchange mirror",
        SimpleBounds<4>(topology.mirror_halo,
                        topology.n_cells_x + 2 * halosize_x, topology.n_cells_y,
                        topology.nens),
        YAKL_LAMBDA(int kk, int i, int j, int n) {
          for (int f = 0; f < num_fields; ++f) {
            const int koff = koffset(f);
            const int _nz = nzs(f);
            for (int l = 0; l < ndofs(f); ++l) {
              // Mirror Top
              this_data(f)(l, _nz + ks + kk, j + js, i, n) =
                  this_data(f)(l, _nz + ks - kk - 1 - koff, j + js, i, n);
              // Mirror Bottom
              this_data(f)(l, ks - kk - 1, j + js, i, n) =
                  this_data(f)(l, ks + kk + koff, j + js, i, n);
            }
          }
        });
}

// EVENTUALLY THIS SHOULD BE MORE CLEVER IE PACK ALL THE FIELDS AT ONCE, ETC.
template <uint num_fields>
template <int num_exchanges>
//void FieldSet<num_fields>::exchange(const std::initializer_list<int> &indices) {
void FieldSet<num_fields>::exchange(const int(&indices)[num_exchanges]) {
  //for (int i : indices) {
  //  fields_arr[i].exchange();
  //}
  
  const auto topology = fields_arr[0].topology;
  const int is = topology.is;
  const int js = topology.js;
  const int ks = topology.ks;
  
  const int halosize_x = topology.halosize_x;
  const int n_cells_x = topology.n_cells_x;
  const int n_cells_y = topology.n_cells_y;

  SArray<real5d, 1, num_exchanges> this_data;
  SArray<int, 1, num_exchanges> ndofs;
  SArray<int, 1, num_exchanges> nzs;
  SArray<int, 1, num_exchanges> koffset;

  int max_nz = 0;
  for (int i = 0; i < num_exchanges; ++i) {
    const int j = indices[i];
    ndofs(i) = this->fields_arr[j].total_dofs;
    max_nz = std::max(max_nz, fields_arr[j]._nz);
    nzs(i) = fields_arr[j]._nz;
    koffset(i) = fields_arr[j].extdof == 1 ? 0 : 1;
    this_data(i) = this->fields_arr[j].data;
  }

  parallel_for(
      "Field Set exchange x",
      SimpleBounds<4>(topology.halosize_x, max_nz,
                      topology.n_cells_y, topology.nens),
      YAKL_LAMBDA(int ii, int k, int j, int n) {
        for (int f = 0; f < num_exchanges; ++f) {
          for (int l = 0; l < ndofs(f); ++l) {
            this_data(f)(l, k + ks, j + js, ii + is - halosize_x, n) =
                this_data(f)(l, k + ks, j + js, ii + is + n_cells_x - halosize_x, n);
            this_data(f)(l, k + ks, j + js, ii + is + n_cells_x, n) =
                this_data(f)(l, k + ks, j + js, ii + is, n);
          }
        }
   });

    parallel_for(
        "Field Set exchange mirror",
        SimpleBounds<4>(topology.mirror_halo,
                        topology.n_cells_x + 2 * halosize_x, topology.n_cells_y,
                        topology.nens),
        YAKL_LAMBDA(int kk, int i, int j, int n) {
          for (int f = 0; f < num_exchanges; ++f) {
            const int koff = koffset(f);
            const int _nz = nzs(f);
            for (int l = 0; l < ndofs(f); ++l) {
              // Mirror Top
              this_data(f)(l, _nz + ks + kk, j + js, i, n) =
                  this_data(f)(l, _nz + ks - kk - 1 - koff, j + js, i, n);
              // Mirror Bottom
              this_data(f)(l, ks - kk - 1, j + js, i, n) =
                  this_data(f)(l, ks + kk + koff, j + js, i, n);
            }
          }
        });
}
