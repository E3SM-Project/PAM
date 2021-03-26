

#ifndef _TOPOLOGY_H_
#define _TOPOLOGY_H_


#include "common.h"
#include <string>
#include "parallel.h"


class Topology {


public:

int n_cells, n_cells_with_halo, n_cells_x, n_cells_y, n_cells_z;
int is, js, ks;
int halosize_x, halosize_y, halosize_z;
int i_beg, j_beg, k_beg;
int i_end, j_end, k_end;
SArray<int,2> x_neigh;
SArray<int,2> y_neigh;
SArray<int,2> z_neigh;
int ll_neigh, ur_neigh, ul_neigh, lr_neigh;
int nprocx, nprocy, nprocz;
int nx_glob, ny_glob, nz_glob;

bool is_initialized;
bool primal;

Topology();
Topology( const Topology &topo) = delete;
Topology& operator=( const Topology &topo) = delete;
void initialize(Parallel &par, bool isprimal);
void printinfo();

};



  Topology::Topology()
  {
    this->is_initialized = false;
    std::cout << "CREATED TOPOLOGY\n";
  }


  void Topology::initialize(Parallel &par, bool isprimal)
  {
    this->n_cells_x = par.nx;
    this->n_cells_y = par.ny;
    this->n_cells_z = par.nz;
    this->halosize_x = par.halox;
    this->halosize_y = par.haloy;
    this->halosize_z = par.haloz;
    this->i_beg = par.i_beg;
    this->i_end = par.i_end;
    this->j_beg = par.j_beg;
    this->j_end = par.j_end;
    this->k_beg = par.k_beg;
    this->k_end = par.k_end;
    this->x_neigh(0) = par.x_neigh(0);
    this->x_neigh(1) = par.x_neigh(1);
    this->y_neigh(0) = par.y_neigh(0);
    this->y_neigh(1) = par.y_neigh(1);
    this->z_neigh(0) = par.z_neigh(0);
    this->z_neigh(1) = par.z_neigh(1);
    this->nprocx = par.nprocx;
    this->nprocy = par.nprocy;
    this->nprocz = par.nprocz;
    this->nx_glob = par.nx_glob;
    this->ny_glob = par.ny_glob;
    this->nz_glob = par.nz_glob;
    
    this->primal = isprimal;
    
    if (ndims == 1) {
      this->n_cells_y = 1;
      this->n_cells_z = 1;
      this->halosize_y = 0;
      this->halosize_z = 0;
      this->j_beg = 0;
      this->j_end = 0;
      this->k_beg = 0;
      this->k_end = 0;
      this->y_neigh(0) = -1;
      this->y_neigh(1) = -1;
      this->z_neigh(0) = -1;
      this->z_neigh(1) = -1;
      this->nprocy = 1;
      this->nprocz = 1;
      this->ny_glob = 1;
      this->nz_glob = 1;
    }

    if (ndims == 2) {
      this->n_cells_z = 1;
      this->halosize_z = 0;
      this->k_beg = 0;
      this->k_end = 0;
      this->z_neigh(0) = -1;
      this->z_neigh(1) = -1;
      this->nprocz = 1;
      this->nz_glob = 1;
      this->ll_neigh = par.ll_neigh;
      this->lr_neigh = par.lr_neigh;
      this->ul_neigh = par.ul_neigh;
      this->ur_neigh = par.ur_neigh;
    }

    this->n_cells = this->n_cells_x * this->n_cells_y * this->n_cells_z;
    this->n_cells_with_halo = (this->n_cells_x + 2*this->halosize_x) * (this->n_cells_y + 2*this->halosize_y) * (this->n_cells_z + 2*this->halosize_z);

    this->is = this->halosize_x;
    this->js = this->halosize_y;
    this->ks = this->halosize_z;

    this->is_initialized = true;
  }

  void Topology::printinfo()
  {
  if (this->primal)
   {std::cout << "topology info: type primal\n" << std::flush;
   std::cout << "nx " << this->n_cells_x << " ny " << this->n_cells_y << " nz " << this->n_cells_z << "\n" << std::flush;
   std::cout << "halosize_x " << this->halosize_x << " halosize_y " << this->halosize_y << " halosize_z " << this->halosize_z << "\n" << std::flush;
   std::cout << "is " << this->is << " js " << this->js << " ks " << this->ks << "\n" << std::flush;
 }
  else
  {std::cout << "topology info: type dual\n" << std::flush;}
   std::cout << "nx " << this->n_cells_x << " ny " << this->n_cells_y << " nz " << this->n_cells_z << "\n" << std::flush;
   std::cout << "halosize_x " << this->halosize_x << " halosize_y " << this->halosize_y << " halosize_z " << this->halosize_z << "\n" << std::flush;
   std::cout << "is " << this->is << " js " << this->js << " ks " << this->ks << "\n" << std::flush;
 }


#endif
