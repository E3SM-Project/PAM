#include "variable_set.h"

  // creates a new VariableSet vs that is a deep clone of self, without copying data over
  // SHARES the exchanges with vs
  void VariableSet::clone(VariableSet &vs, std::string name) {
    baseName = name;
    num_fields = vs.num_fields;
    //CREATE FIELDS_ARR AND EXCHANGES_ARR

    for (int i=0; i<num_fields; i++)
    {
      fields_arr[i].clone(vs.fields_arr[i], vs.fields_arr[i].name);
      exchanges_arr[i] = vs.exchanges_arr[i];
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
      exchanges_arr[i].exchange.pack(fields_arr[i]);
      exchanges_arr[i].exchange.exchange();
      exchanges_arr[i].exchange.unpack(fields_arr[i]);
    }
  }
