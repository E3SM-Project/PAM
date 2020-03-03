#ifndef _FIELDS_H_
#define _FIELDS_H_


#include "common.h"
#include "topology.h"
#include "exchange.h"
#include "STDLIB"

class Field {

public :

  realArr data;
  int ndofs0, ndofs1, ndofs2, ndofs3;
  Topology topology;
  Exchange exchange;
  std::string name;
  int total_dofs;

  void initialize(Topology &topo, Exchange &exch, std::string fieldName, int ndof0, int ndof1, int ndof2 = 0, int ndof3 = 0);
  void waxpy(real alpha, Field &x, Field &y);
  void copy(Field & f);
  void clone(Field &f, std::string fieldName);
  int YAKL_INLINE get_offset(int ndof);
};

#endif
