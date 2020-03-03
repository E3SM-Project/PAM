

#include "topology.h"

// THIS IS BROKEN FOR PARALLEL
  void PeriodicTopology::initialize(int nx, int ny, int nz) {

    n_cells_x = nx;
    n_cells_y = ny;
    n_cells_z = nz;

    if (ndims == 1) {
      n_cells_y = 1;
      n_cells_z = 1;
    }
    if (ndims == 2) {
      n_cells_z = 1;
    }

    n_cells = n_cells_x * n_cells_y * n_cells_z;

    is = halosize_x;
    js = halosize_y;
    ks = halosize_z;

  }

  void PeriodicTopology::create_arr(std::string formName, realArr &var, int ndofs0, int ndofs1, int ndofs2, int ndofs3) {
    var = realArr(formName, get_total_dofs(ndims, ndofs0, ndofs1, ndofs2, ndofs3), n_cells_z+haloxsize_z, n_cells_y+haloxsize_y, n_cells_x+haloxsize_x);
}
