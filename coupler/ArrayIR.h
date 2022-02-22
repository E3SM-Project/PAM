
#pragma once

#include <iostream>
#include <string>
#include <vector>

// The ArrayIR class holds metadata for a contiguous array in either host or device memory
// For safety, it should be assumed that device memory cannot be accessed on the host and
// vice versa. Any assumptions beyond this are up to the user if it is known a priori that,
// for instance, memory was allocated with managed or unified memory.
// This is primarily to communicate arrays in an agnostic format between YAKL Array objects
// and Kokkos View objects.
// 
// It should be assumed that the data pointed to is contiguous, dimensioned by the given
// dimensions with the left-most dimension varying the slowest and the rightmost varying
// the fastest (i.e., C-style ordering).
// 
// The user is strongly encouraged to use routines that directly translate the ArrayIR
// metadata to and from YAKL Array objects and Kokkos View objects to minimize the
// possibility of wrongly transcribing the memory space or dimensions.
// 
// There is no allocation or deallocation in this class. It merely holds metadata for
// existing allocated Array or View objects (or, really, any array-like object that exists
// in a context of potentially separate memory spaces).
// 
// Also for safety, it is expected that metadata is set only via the constructor and then is
// read-only from there. There is no reason for this data to change after creation.

namespace pam {

  int constexpr IRMemHost   = 1;
  int constexpr IRMemDevice = 2;

  template <class T, int N, int memSpace>
  struct ArrayIR {
  public:
    T *          my_data;
    unsigned int dims[N];
    std::string  label;

    ArrayIR () { this->my_data = nullptr; };
    ~ArrayIR() { this->my_data = nullptr; };

    ArrayIR(T *my_data , std::vector<unsigned int> dims , std::string label = "") {
      // Check and set data pointer
      if (my_data == nullptr) {
        std::cerr << "Error: ArrayIR constructor called with data == nullptr." << std::endl;
        throw "";
      }
      this->my_data = my_data;

      // Check and set dimensions
      if (N != dims.size() ) {
        std::cerr << "Error: ArrayIR constructor dimension sizes do not match the template parameter N." << std::endl;
        throw "";
      }
      for (int i=0; i < dims.size(); i++) {
        if (dims[i] <= 0) {
          std::cerr << "Error: ArrayIR constructor dimension is <= 0." << std::endl;
          throw "";
        }
        this->dims[i] = dims[i];
      }

      // Set label
      this->label = label;
    }


    ArrayIR            (ArrayIR<T,N,memSpace> const &rhs) {
      this->my_data = rhs.my_data;
      for (int i=0; i < N; i++) { this->dims[i] = rhs.dims[i]; }
      this->label = rhs.label;
    }
    ArrayIR & operator=(ArrayIR<T,N,memSpace> const &rhs) {
      if (this == &rhs) return *this;
      this->my_data = rhs.my_data;
      for (int i=0; i < N; i++) { this->dims[i] = rhs.dims[i]; }
      this->label = rhs.label;
      return *this;
    }
    ArrayIR            (ArrayIR<T,N,memSpace> &&rhs) {
      this->my_data = rhs.my_data;
      for (int i=0; i < N; i++) { this->dims[i] = rhs.dims[i]; }
      this->label = rhs.label;
    }
    ArrayIR& operator= (ArrayIR<T,N,memSpace> &&rhs) {
      if (this == &rhs) return *this;
      this->my_data = rhs.my_data;
      for (int i=0; i < N; i++) { this->dims[i] = rhs.dims[i]; }
      this->label = rhs.label;
      return *this;
    }

    T * data    () const { return this->my_data; }
    T * get_data() const { return this->my_data; }

    std::vector<unsigned int> get_dims() const {
      std::vector<unsigned int> ret(N);
      for (int i=0; i < N; i++) { ret[i] = this->dims[i]; }
      return ret;
    }

    std::string get_label() const { return this->label; }

    bool initialized() const { return this->my_data != nullptr; }
  };

}


