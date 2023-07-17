

#pragma once

#include "common.h"
#include "parallel.h"

class Topology {
public:
  // int n_cells_layers, n_cells_layers_with_halo;
  // int n_cells_interfaces, n_cells_interfaces_internal,
  // n_cells_interfaces_with_halo;
  int n_cells_x, n_cells_y;
  int is, js;
  int halosize_x, halosize_y;
  int i_beg, j_beg;
  int i_end, j_end;
  SArray<int, 1, 2> x_neigh;
  SArray<int, 1, 2> y_neigh;
  int ll_neigh, ur_neigh, ul_neigh, lr_neigh;
  int nprocx, nprocy;
  int nx_glob, ny_glob, nl, ni;
  int ks, mirror_halo;
  int nens;

  BND_TYPE xbnd, ybnd;

  bool is_initialized;
  bool primal;

  Topology();
  void initialize(Parallel &par, bool isprimal);
  void printinfo() const;
};

Topology::Topology() { this->is_initialized = false; }

void Topology::initialize(Parallel &par, bool isprimal) {
  this->n_cells_x = par.nx;
  this->n_cells_y = par.ny;
  this->halosize_x = par.halox;
  this->halosize_y = par.haloy;
  this->i_beg = par.i_beg;
  this->i_end = par.i_end;
  this->j_beg = par.j_beg;
  this->j_end = par.j_end;
  this->x_neigh(0) = par.x_neigh(0);
  this->x_neigh(1) = par.x_neigh(1);
  this->y_neigh(0) = par.y_neigh(0);
  this->y_neigh(1) = par.y_neigh(1);
  this->nprocx = par.nprocx;
  this->nprocy = par.nprocy;
  this->nx_glob = par.nx_glob;
  this->ny_glob = par.ny_glob;
  this->nens = par.nens;

  this->xbnd = par.xbnd;
  this->ybnd = par.ybnd;

  this->primal = isprimal;

  this->nl = par.nz;
  this->ni = par.nz + 1;
#ifdef PAMC_EXTRUDED
  if (this->primal) {
    this->nl = par.nz - 1;
    this->ni = par.nz;
  }
#endif

  if (ndims == 1) {
    this->n_cells_y = 1;
    this->halosize_y = 0;
    this->j_beg = 0;
    this->j_end = 0;
    this->y_neigh(0) = -1;
    this->y_neigh(1) = -1;
    this->nprocy = 1;
    this->ny_glob = 1;
    this->ybnd = BND_TYPE::NONE;
  }

  if (ndims == 2) {
    this->ll_neigh = par.ll_neigh;
    this->lr_neigh = par.lr_neigh;
    this->ul_neigh = par.ul_neigh;
    this->ur_neigh = par.ur_neigh;
  }

  //   this->n_cells_layers = this->n_cells_x * this->n_cells_y * this->nl;
  //   this->n_cells_layers_with_halo = (this->n_cells_x + 2*this->halosize_x) *
  //   (this->n_cells_y + 2*this->halosize_y) * (this->nl +
  //   2*this->mirror_halo); this->n_cells_interfaces = this->n_cells_x *
  //   this->n_cells_y * this->ni; this->n_cells_interfaces_internal =
  //   this->n_cells_x * this->n_cells_y * (this->ni-2);
  //   this->n_cells_interfaces_with_halo = (this->n_cells_x +
  //   2*this->halosize_x) * (this->n_cells_y + 2*this->halosize_y) * (this->ni
  //   + 2*this->mirror_halo);

  this->is = this->halosize_x;
  this->js = this->halosize_y;

#ifdef PAMC_EXTRUDED
  this->mirror_halo = mirroringhalo;
  this->ks = this->mirror_halo;
#else
  this->mirror_halo = 0;
  this->ks = 0;
#endif

  this->is_initialized = true;
}

void Topology::printinfo() const {
  if (this->primal) {
    std::cout << "topology info: type primal\n" << std::flush;
  } else {
    std::cout << "topology info: type dual\n" << std::flush;
  }

  std::cout << "nx " << this->n_cells_x << " ny " << this->n_cells_y << " nl "
            << this->nl << " ni " << this->ni << "\n"
            << std::flush;
  std::cout << "halosize_x " << this->halosize_x << " halosize_y "
            << this->halosize_y << " mirroring_halo " << this->mirror_halo
            << "\n"
            << std::flush;
  std::cout << "is " << this->is << " js " << this->js << " ks " << this->ks
            << "\n"
            << std::flush;
}
