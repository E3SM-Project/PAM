

#include "topology.h"

// THIS IS BROKEN FOR PARALLEL
  void PeriodicTopology::initialize(int nx, int ny=1, int nz=1) {

    n_cells_x = nx;
    n_cells_y = ny;
    n_cells_z = nz;

    halosize_x = maxhalosize_x;
    halosize_y = maxhalosize_y;
    halosize_z = maxhalosize_z;

    if (ndims == 1) {
      n_cells_y = 1;
      n_cells_z = 1;
      halosize_y = 0;
      halosize_z = 0;
    }
    if (ndims == 2) {
      n_cells_z = 1;
      halosize_z = 0;
    }

    n_cells = n_cells_x * n_cells_y * n_cells_z;

    is = halosize_x;
    js = halosize_y;
    ks = halosize_z;

  }

  void PeriodicTopology::create_field(std::string formName, realArr &var, int ndofs0, int ndofs1, int ndofs2=0, int ndofs3=0) {
    int ndofs_per_cell;
    if (ndims == 1) { ndofs_per_cell = ndofs0 + ndofs1; }
    if (ndims == 2) { ndofs_per_cell = ndofs0 + 2*ndofs1 + ndofs2; }
    if (ndims == 3) { ndofs_per_cell = ndofs0 + 3*ndofs1 + 3*ndofs2 + ndofs3; }
    var = realArr(formName, ndofs_per_cell, n_cells_z+haloxsize_z, n_cells_y+haloxsize_y, n_cells_x+haloxsize_x);
}
