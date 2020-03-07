#ifndef _VARIABLESET_H_
#define _VARIABLESET_H_


#include "common.h"
#include "exchange.h"
#include "fields.h"
#include "topology.h"
#include <array>
#include <string>


template<uint ndims, uint num_fields> class VariableSet {

public:

  std::array<Field<ndims>, num_fields> fields_arr;
  std::string baseName;
  bool is_initialized;

  VariableSet();
  VariableSet( const VariableSet<ndims,num_fields> &vs) = delete;
  VariableSet& operator=( const VariableSet<ndims,num_fields> &vs) = delete;
  void printinfo();
  void initialize(const std::string name, std::array<std::string, num_fields> &names_arr, std::array<const Topology<ndims> *, num_fields> &topo_arr, SArray<int, num_fields, 4> &ndofs_arr);
  void initialize(const VariableSet<ndims, num_fields> &vs, const std::string name);
  void copy(const VariableSet<ndims, num_fields> &vs);
  void waxpy(real alpha, const VariableSet<ndims, num_fields> &x, const VariableSet<ndims, num_fields> &y);
};


template<uint ndims, uint num_fields> class ExchangeSet {

public:
  std::array<Exchange<ndims>, num_fields> exchanges_arr;
  bool is_initialized;

  ExchangeSet();
  ExchangeSet( const ExchangeSet<ndims,num_fields> &exch) = delete;
  ExchangeSet& operator=( const ExchangeSet<ndims,num_fields> &exch) = delete;
  void initialize(std::array<const Topology<ndims> *, num_fields> &topo_arr, SArray<int, num_fields, 4> &ndofs_arr);
  void initialize(const ExchangeSet<ndims, num_fields> &es);
  void printinfo();
  void exchange_variable_set(VariableSet<ndims, num_fields> &vs);
};



template<uint ndims, uint num_fields> ExchangeSet<ndims,num_fields>::ExchangeSet()
{
  this->is_initialized = false;
  std::cout << "CREATED EXCHANGE SET\n";
}


template<uint ndims, uint num_fields> void ExchangeSet<ndims,num_fields>::initialize(std::array<const Topology<ndims> *, num_fields> &topo_arr, SArray<int, num_fields, 4> &ndofs_arr)
{
  for (int i=0; i<num_fields; i++)
  {
    this->exchanges_arr[i].initialize(*topo_arr[i], ndofs_arr(i,0), ndofs_arr(i,1), ndofs_arr(i,2), ndofs_arr(i,3));
  }
  this->is_initialized = true;
}

template<uint ndims, uint num_fields> void ExchangeSet<ndims,num_fields>::printinfo()
{
  std::cout << "exchange set info\n" << std::flush;
  for (int i=0; i<num_fields; i++)
  {
    this->exchanges_arr[i].printinfo();
  }
}

//EVENTUALLY THIS SHOULD BE MORE CLEVER IE PACK ALL THE FIELDS AT ONCE, ETC.
template<uint ndims, uint num_fields> void ExchangeSet<ndims,num_fields>::exchange_variable_set(VariableSet<ndims, num_fields> &vs)
{
  for (int i=0; i<num_fields; i++)
  {
    this->exchanges_arr[i].pack(vs.fields_arr[i]);
    this->exchanges_arr[i].exchange();
    this->exchanges_arr[i].unpack(vs.fields_arr[i]);
  }
}

template<uint ndims, uint num_fields> void ExchangeSet<ndims,num_fields>::initialize(const ExchangeSet<ndims, num_fields> &es)
{
  for (int i=0; i<num_fields; i++)
  {
    this->exchanges_arr[i].initialize(es.exchanges_arr[i]);
  }
  this->is_initialized = true;
}




    template<uint ndims, uint num_fields> VariableSet<ndims,num_fields>::VariableSet()
    {
      this->is_initialized = false;
      std::cout << "CREATED VARIABLE SET\n";
    }

  template<uint ndims, uint num_fields> void VariableSet<ndims,num_fields>::printinfo()
  {
  std::cout << "variable set info " << this->baseName << "\n" << std::flush;
  for (int i=0; i<num_fields; i++)
  {
    this->fields_arr[i].printinfo();
  }
  }


  template<uint ndims, uint num_fields> void VariableSet<ndims,num_fields>::initialize(const std::string name, std::array<std::string, num_fields> &names_arr, std::array<const Topology<ndims> *, num_fields> &topo_arr, SArray<int, num_fields, 4> &ndofs_arr)
  {
    this->baseName = name;
    for (int i=0; i<num_fields; i++)
    {
      this->fields_arr[i].initialize(*topo_arr[i], names_arr[i], ndofs_arr(i,0), ndofs_arr(i,1), ndofs_arr(i,2), ndofs_arr(i,3));
    }
    this->is_initialized = true;
  }

  // creates a new VariableSet that has the same fields as vs
  template<uint ndims, uint num_fields> void VariableSet<ndims,num_fields>::initialize(const VariableSet<ndims, num_fields> &vs, const std::string name) //,  EXCHANGE_TYPE exch_type)
  {
    this->baseName = name;
    for (int i=0; i<num_fields; i++)
    {
      this->fields_arr[i].initialize(vs.fields_arr[i], vs.fields_arr[i].name);
    }
    this->is_initialized = true;
  }

  // copies data from vs into self
  template<uint ndims, uint num_fields> void VariableSet<ndims,num_fields>::copy(const VariableSet<ndims, num_fields> &vs)
  {
    for (int i=0; i<num_fields; i++)
    {
      this->fields_arr[i].copy(vs.fields_arr[i]);
    }
  }

  // Computes w (self) = alpha x + y
  template<uint ndims, uint num_fields> void VariableSet<ndims,num_fields>::waxpy(real alpha, const VariableSet<ndims, num_fields> &x, const VariableSet<ndims, num_fields> &y)
  {
    for (int i=0; i<num_fields; i++)
    {
      this->fields_arr[i].waxpy(alpha, x.fields_arr[i], y.fields_arr[i]);
    }
  }




#endif
