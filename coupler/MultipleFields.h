
#pragma once

#include "pam_const.h"

// Aggregate multiple fields into a single field that makes it easier to operate
// on them together inside the same kernel. used mostly for tracers
template <int MAX_FIELDS, class T>
class MultipleFields {
public:
  yakl::SArray<T,1,MAX_FIELDS> fields;
  int num_fields;

  YAKL_INLINE MultipleFields() { num_fields = 0; }

  YAKL_INLINE MultipleFields(MultipleFields const &rhs) {
    this->num_fields = rhs.num_fields;
    for (int i=0; i < num_fields; i++) {
      this->fields(i) = rhs.fields(i);
    }
  }

  YAKL_INLINE MultipleFields & operator=(MultipleFields const &rhs) {
    this->num_fields = rhs.num_fields;
    for (int i=0; i < num_fields; i++) {
      this->fields(i) = rhs.fields(i);
    }
    return *this;
  }

  YAKL_INLINE MultipleFields(MultipleFields &&rhs) {
    this->num_fields = rhs.num_fields;
    for (int i=0; i < num_fields; i++) {
      this->fields(i) = rhs.fields(i);
    }
  }

  YAKL_INLINE MultipleFields& operator=(MultipleFields &&rhs) {
    this->num_fields = rhs.num_fields;
    for (int i=0; i < num_fields; i++) {
      this->fields(i) = rhs.fields(i);
    }
    return *this;
  }

  YAKL_INLINE void add_field( T &field ) {
    this->fields(num_fields) = field;
    num_fields++;
  }

  YAKL_INLINE real &operator() (int tr, int k, int j, int i, int iens) const {
    return this->fields(tr)(k,j,i,iens);
  }
};



template <class T, int N>
using MultiField = MultipleFields< max_fields , Array<T,N,memDevice,styleC> >;



