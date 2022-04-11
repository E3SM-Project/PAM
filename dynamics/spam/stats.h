#pragma once

#include "common.h"
#include "topology.h"
#include "geometry.h"
#include "params.h"
#include "parallel.h"

class Stat
{
public:
  real3d data;
  std::string name;
  int ndofs, size, nens;

void initialize(std::string statName, int statndof, int statsize, int statnens, int masterproc)
{
  name = statName;
  ndofs = statndof;
  size = statsize;
  nens = statnens;

  if (masterproc)
    {
    data = real3d(name.c_str(), ndofs, size, nens);
  }
}

};


template <uint nprog, uint nconst, uint nstats> class Stats
{
public:
  std::array<Stat,nstats> stats_arr;
  MPI_Request Req [nstats];
  MPI_Status  Status[nstats];
  int ierr;
  int masterproc;
  const Topology *primal_topology;
  const Topology *dual_topology;
  Geometry<1,1,1> *primal_geometry;
  Geometry<1,1,1> *dual_geometry;
  int nens;
  int statsize;

  void initialize(Parameters &params, Parallel &par, const Topology &primal_topo, const Topology &dual_topo, Geometry<1,1,1> &primal_geom, Geometry<1,1,1> &dual_geom)
  {
    this->primal_topology = &primal_topo;
    this->dual_topology = &dual_topo;
    this->primal_geometry = &primal_geom;
    this->dual_geometry = &dual_geom;
    this->nens = params.nens;
    this->masterproc = par.masterproc;
    this->statsize = params.Nsteps/params.Nstat + 1;
  }
  
  virtual void compute( VariableSet<nprog> &progvars,  VariableSet<nconst> &constvars, int tind) {};

};