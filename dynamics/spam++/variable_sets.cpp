#include "variable_set.h"


  void VariableSet::initialize(std::string name, bool exchange = true){
    baseName = name;
    has_exchange = exchange;
  }

  // creates a new VariableSet that has the same fields as vs

  void VariableSet::clone(VariableSet &vs, std::string name, int exchange_type = EXCHANGE_TYPE_NONE) {
    baseName = name;

    if (exchange_type == EXCHANGE_TYPE_NONE) {
      has_exchange = false;
    }
    else {
      has_exchange = true;
    }

    for (int i=0; i<num_fields; i++)
    {
      fields_arr[i].clone(vs.fields_arr[i], vs.fields_arr[i].name);
      if (exchange_type == EXCHANGE_TYPE_SHARED && vs.has_exchange) { exchanges_arr[i] = vs.exchanges_arr[i];}
      if (exchange_type == EXCHANGE_TYPE_CLONED && vs.has_exchange) { exchanges_arr[i].clone(vs.exchanges_arr[i]); }
    }
  }

  // copies data from vs into self
  void VariableSet::copy(VariableSet &vs) {
    for (int i=0; i<num_fields; i++)
    {
      fields_arr[i].copy(vs.fields_arr[i]);
    }
  }

  // Computes w (self) = alpha x + y
  void VariableSet::waxpy(real alpha, VariableSet &x, VariableSet &y) {
    for (int i=0; i<num_fields; i++)
    {
      fields_arr[i].waxpy(alpha, x.fields_arr[i], y.fields_arr[i]);
    }
  }

// EVENTUALLY THIS SHOULD BE MORE CLEVER IE PACK ALL THE FIELDS AT ONCE, ETC.
  void exchange() {
    for (int i=0; i<num_fields; i++)
    {
      exchanges_arr[i].pack(fields_arr[i]);
      exchanges_arr[i].exchange();
      exchanges_arr[i].unpack(fields_arr[i]);
    }
  }
