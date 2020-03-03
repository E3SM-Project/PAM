#ifndef _VARIABLESET_H_
#define _VARIABLESET_H_


#include "common.h"
#include "exchange.h"
#include "field.h"
#include "topology.h"
#include "STDLIB"

class VariableSet {

public:

  int num_fields;
  std::vector<Field> fields_arr;
  std::vector<Exchange> exchanges_arr;
  std::string baseName;

  void clone(VariableSet &vs, std::string name);
  void copy(VariableSet &vs);
  void waxpy(real alpha, VariableSet &x, VariableSet &y);
  void exchange();

  void create(Topology &topology, bool create_exchange = true);

};
#endif
