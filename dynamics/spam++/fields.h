#ifndef _FIELDS_H_
#define _FIELDS_H_


#include "common.h"
#include "topology.h"
#include "exchange.h"
#include <string>

// ACTUALLY TEMPLATE THIS BASED ON NDOF0,NDOF1, ETC.!
// WITH DEFAULT SIZES...

template<int ndims, int ndof0, int ndof1, int ndof2, int ndof3> class Field {

public :

  realArr data;
  Topology topology;
  std::string name;
  int total_dofs;

  void initialize(Topology &topo, std::string fieldName);
  void waxpy(real alpha, Field &x, Field &y);
  void copy(Field & f);
  void clone(Field &f, std::string fieldName);
  int YAKL_INLINE get_offset(int ndof);
};

void YAKL_INLINE get_total_dofs(int ndims, int ndof0, int ndof1, int ndof2, int ndof3);
#endif
