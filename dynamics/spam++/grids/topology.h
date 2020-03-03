

#ifndef _TOPOLOGY_H_
#define _TOPOLOGY_H_


#include "common.h"
#include "fields.h"
#include <string>

template<int ndims, int halosize_x, int halosize_y, int halosize_z> class Topology {
public:

  void create_arr(std::string formName, realArr &var, int ndofs0, int ndofs1, int ndofs2, int ndofs3) {

};

class PeriodicTopology : Topology {


public:

  void initialize(int nx, int ny, int nz);

int n_cells, n_cells_x, n_cells_y, n_cells_z;
int is, js, ks;

};

#endif
