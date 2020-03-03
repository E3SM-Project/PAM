#ifndef _VARIABLESET_H_
#define _VARIABLESET_H_


#include "common.h"
#include "exchange.h"
#include "field.h"
#include "topology.h"
#include <array>
#include <string>

int constexpr EXCHANGE_TYPE_NONE = 0;
int constexpr EXCHANGE_TYPE_CLONE = 0;
int constexpr EXCHANGE_TYPE_SHARED = 0;

template<int num_fields> class VariableSet {

public:

  std::array<Field,num_fields> fields_arr;
  std::array<Exchange,num_fields> exchanges_arr;
  std::string baseName;
  bool has_exchange;

  void copy(VariableSet &vs);
  void waxpy(real alpha, VariableSet &x, VariableSet &y);
  void exchange();

  void clone(VariableSet &vs, std::string name, int exchange_type = EXCHANGE_TYPE_NONE);

  void initialize(std::string name, bool exchange);

};

#endif
