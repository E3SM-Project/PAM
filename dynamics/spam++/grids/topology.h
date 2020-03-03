

#ifndef _TOPOLOGY_H_
#define _TOPOLOGY_H_


#include "common.h"
#include "STDLIB"

class Topology {

  void initialize(int nx, int ny=1, int nz=1);
  void create_field(std::string formName, realArr &var, int ndofs0, int ndofs1, int ndofs2=0, int ndofs3=0) {

};

class PeriodicTopology : Topology {

public:

int n_cells, n_cells_x, n_cells_y, n_cells_z;
int is, js, ks;
int halosize_x, halosize_y, halosize_z;

};

#endif
