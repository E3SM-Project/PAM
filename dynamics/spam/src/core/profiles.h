#pragma once

#include "common.h"
#include "topology.h"

class Profile {

public:
  real3d data;
  Topology topology;
  std::string name;
  int total_dofs; // total number of dofs in a "horiz" slice ie at each layer or
                  // interface
  // basedof = 0,1,2 in 2D or 0,1 in 1D ie vertices,edges,cells or
  // vertices,cells extdof = 0,1 (0 = interfaces, 1=layers)
  int basedof, extdof, ndofs;
  int _nz; //, _nloop, _nloop_halo;

  bool is_initialized;
  Profile();
  // Profile( const Profile &f) = delete;
  // Profile& operator=( const Profile &f) = delete;
  void printinfo();
  void initialize(const Topology &topo, const std::string profName, int bdof,
                  int edof, int nd);
};

Profile::Profile() { this->is_initialized = false; }

void Profile::printinfo() {
  std::cout << "profile info " << this->name << "\n" << std::flush;
  std::cout << this->data << std::flush;
  std::cout << "total_dofs " << this->total_dofs << " basedof " << this->basedof
            << " extdof " << this->extdof << " ndofs " << this->ndofs << "\n"
            << std::flush;
}

void Profile::initialize(const Topology &topo, const std::string profName,
                         int bdof, int edof, int nd) {

  this->topology = topo;
  this->name = profName;
  this->basedof = bdof;
  this->extdof = edof;
  this->ndofs = nd;

  if (this->extdof == 0) {
    this->_nz = this->topology.ni;
  }
  if (this->extdof == 1) {
    this->_nz = this->topology.nl;
  }

  this->total_dofs = this->ndofs;
  if (ndims == 2 && this->basedof == 1) {
    this->total_dofs = 2 * this->ndofs;
  } // 2 edges per cell in 2D

  this->data = real3d(this->name.c_str(), this->total_dofs,
                      this->_nz, // + 2 * this->topology.mirror_halo,
                      this->topology.nens);

  this->is_initialized = true;
}
