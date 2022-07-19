#pragma once

#include "common.h"
#include "exchange.h"
#include "fields.h"
#include "topology.h"

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
                  const std::array<FieldDescription, num_fields> &desc_arr);
  void initialize(const FieldSet<num_fields> &vs, const std::string name);
  void copy(const FieldSet<num_fields> &vs);
  void waxpy(real alpha, const FieldSet<num_fields> &x,
             const FieldSet<num_fields> &y);
  void waxpbypcz(real alpha, real beta, real gamma,
                 const FieldSet<num_fields> &x, const FieldSet<num_fields> &y,
                 const FieldSet<num_fields> &z);
};

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
  void exchange_variable_set(FieldSet<num_fields> &vs);
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
template <uint num_fields>
void ExchangeSet<num_fields>::exchange_variable_set(FieldSet<num_fields> &vs) {
  for (int i = 0; i < num_fields; i++) {
    this->exchanges_arr[i].exchange_field(vs.fields_arr[i]);
  }
}

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
    const std::array<FieldDescription, num_fields> &desc_arr) {
  this->baseName = name;
  for (int i = 0; i < num_fields; i++) {
    this->fields_arr[i].initialize(desc_arr[i].topology, desc_arr[i].name,
                                   desc_arr[i].basedof, desc_arr[i].extdof,
                                   desc_arr[i].ndofs);
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
  for (int i = 0; i < num_fields; i++) {
    this->fields_arr[i].waxpy(alpha, x.fields_arr[i], y.fields_arr[i]);
  }
}

// Computes w (self) = alpha x + beta * y + gamma * z
template <uint num_fields>
void FieldSet<num_fields>::waxpbypcz(real alpha, real beta, real gamma,
                                     const FieldSet<num_fields> &x,
                                     const FieldSet<num_fields> &y,
                                     const FieldSet<num_fields> &z) {
  for (int i = 0; i < num_fields; i++) {
    this->fields_arr[i].waxpbypcz(alpha, beta, gamma, x.fields_arr[i],
                                  y.fields_arr[i], z.fields_arr[i]);
  }
}
