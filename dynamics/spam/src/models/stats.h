#pragma once

#include "common.h"
#include "geometry.h"
#include "parallel.h"
#include "topology.h"

class Stat {
public:
  realHost3d data;
  std::string name;
  int ndofs, size, nens;

  void initialize(std::string statName, int statndof, int statsize,
                  int statnens, int masterproc) {
    name = statName;
    ndofs = statndof;
    size = statsize;
    nens = statnens;

    if (masterproc) {
      data = realHost3d(name.c_str(), ndofs, size, nens);
    }
  }
};

class Stats {
public:
  std::array<Stat, nstats> stats_arr;
  MPI_Request Req[nstats];
  MPI_Status Status[nstats];
  int ierr;
  int masterproc;
  Geometry<Straight> primal_geometry;
  Geometry<Twisted> dual_geometry;
  Equations *equations;
  int nens;
  int statsize;

  void initialize(ModelParameters &params, Parallel &par,
                  const Geometry<Straight> &primal_geom,
                  const Geometry<Twisted> &dual_geom, Equations &eqs) {
    this->primal_geometry = primal_geom;
    this->dual_geometry = dual_geom;
    this->equations = &eqs;
    this->nens = params.nens;
    this->masterproc = par.masterproc;
    this->statsize = params.statSize + 1;
  }

  virtual void compute(FieldSet<nprognostic> &progvars,
                       FieldSet<nconstant> &constvars, int tind){};
};
