

#ifndef _TOPOLOGY_H_
#define _TOPOLOGY_H_


#include "common.h"
#include <string>


template<uint ndims> class Topology {


public:

int n_cells, n_cells_with_halo, n_cells_x, n_cells_y, n_cells_z;
int is, js, ks;
int halosize_x, halosize_y, halosize_z;
bool is_initialized;
Topology<ndims>();
Topology<ndims>( const Topology<ndims> &topo) = delete;
Topology<ndims>& operator=( const Topology<ndims> &topo) = delete;
void initialize(int nx, int ny, int nz, int halox, int haloy, int haloz);
void printinfo();

};



  template<uint ndims> Topology<ndims>::Topology()
  {
    this->is_initialized = false;
    std::cout << "CREATED TOPOLOGY\n";
  }


// THIS IS BROKEN FOR PARALLEL
  template<uint ndims> void Topology<ndims>::initialize(int nx, int ny, int nz, int halox, int haloy, int haloz)
  {
    this->n_cells_x = nx;
    this->n_cells_y = ny;
    this->n_cells_z = nz;
    this->halosize_x = halox;
    this->halosize_y = haloy;
    this->halosize_z = haloz;

    if (ndims == 1) {
      this->n_cells_y = 1;
      this->n_cells_z = 1;
      this->halosize_y = 0;
      this->halosize_z = 0;
    }
    if (ndims == 2) {
      this->n_cells_z = 1;
      this->halosize_z = 0;
    }

    this->n_cells = this->n_cells_x * this->n_cells_y * this->n_cells_z;
    this->n_cells_with_halo = (this->n_cells_x + 2*this->halosize_x) * (this->n_cells_y + 2*this->halosize_y) * (this->n_cells_z + 2*this->halosize_z);

    this->is = this->halosize_x;
    this->js = this->halosize_y;
    this->ks = this->halosize_z;

    this->is_initialized = true;
  }

  template<uint ndims> void Topology<ndims>::printinfo()
  {
   std::cout << "topology info\n" << std::flush;
   std::cout << "nx " << this->n_cells_x << " ny " << this->n_cells_y << " nz " << this->n_cells_z << "\n" << std::flush;
   std::cout << "halosize_x " << this->halosize_x << " halosize_y " << this->halosize_y << " halosize_z " << this->halosize_z << "\n" << std::flush;
   std::cout << "is " << this->is << " js " << this->js << " ks " << this->ks << "\n" << std::flush;
 }


#endif
