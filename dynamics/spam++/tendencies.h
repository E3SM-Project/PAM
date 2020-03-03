
#ifndef _TENDENCIES_H_
#define _TENDENCIES_H_

class Tendencies {

  void compute_rhs(VariableSet &const_vars, VariableSet &x, VariableSet &diagnostic_vars, VariableSet &xtend, Topology &topology);
};

#endif
