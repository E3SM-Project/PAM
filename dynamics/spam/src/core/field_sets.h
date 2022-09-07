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
  void copy(const FieldSet<num_fields> &vs,
            const FIELDOP_EXTENT extent = FIELDOP_EXTENT::WITHOUT_HALOS);
  void waxpy(real alpha, const FieldSet<num_fields> &x,
             const FieldSet<num_fields> &y,
             const FIELDOP_EXTENT extent = FIELDOP_EXTENT::WITHOUT_HALOS);
  void waxpby(real alpha, real beta, const FieldSet<num_fields> &x,
              const FieldSet<num_fields> &y,
              const FIELDOP_EXTENT extent = FIELDOP_EXTENT::WITHOUT_HALOS);
  void waxpbypcz(real alpha, real beta, real gamma,
                 const FieldSet<num_fields> &x, const FieldSet<num_fields> &y,
                 const FieldSet<num_fields> &z,
                 const FIELDOP_EXTENT extent = FIELDOP_EXTENT::WITHOUT_HALOS);
  void exchange();
  void exchange(const std::initializer_list<int> &indices);
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
void FieldSet<num_fields>::copy(const FieldSet<num_fields> &vs,
                                const FIELDOP_EXTENT extent) {
  for (int i = 0; i < num_fields; i++) {
    this->fields_arr[i].copy(vs.fields_arr[i], extent);
  }
}

// Computes w (self) = alpha x + y
template <uint num_fields>
void FieldSet<num_fields>::waxpy(real alpha, const FieldSet<num_fields> &x,
                                 const FieldSet<num_fields> &y,
                                 const FIELDOP_EXTENT extent) {
  for (int i = 0; i < num_fields; i++) {
    this->fields_arr[i].waxpy(alpha, x.fields_arr[i], y.fields_arr[i], extent);
  }
}

// Computes w (self) = alpha x + beta y
template <uint num_fields>
void FieldSet<num_fields>::waxpby(real alpha, real beta,
                                  const FieldSet<num_fields> &x,
                                  const FieldSet<num_fields> &y,
                                  const FIELDOP_EXTENT extent) {
  for (int i = 0; i < num_fields; i++) {
    this->fields_arr[i].waxpby(alpha, beta, x.fields_arr[i], y.fields_arr[i],
                               extent);
  }
}

// Computes w (self) = alpha x + beta * y + gamma * z
template <uint num_fields>
void FieldSet<num_fields>::waxpbypcz(real alpha, real beta, real gamma,
                                     const FieldSet<num_fields> &x,
                                     const FieldSet<num_fields> &y,
                                     const FieldSet<num_fields> &z,
                                     const FIELDOP_EXTENT extent) {
  for (int i = 0; i < num_fields; i++) {
    this->fields_arr[i].waxpbypcz(alpha, beta, gamma, x.fields_arr[i],
                                  y.fields_arr[i], z.fields_arr[i], extent);
  }
}

template <uint num_fields> void FieldSet<num_fields>::exchange() {
  for (auto &f : fields_arr) {
    f.exchange();
  }
}

// EVENTUALLY THIS SHOULD BE MORE CLEVER IE PACK ALL THE FIELDS AT ONCE, ETC.
template <uint num_fields>
void FieldSet<num_fields>::exchange(const std::initializer_list<int> &indices) {
  for (int i : indices) {
    fields_arr[i].exchange();
  }
}
