#pragma once

#include "common.h"
#include "geometry.h"
#include "parallel.h"
#include "topology.h"
#include "pam_coupler.h" //Has DataManager and pam_const
using pam::PamCoupler;

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
  int nens;
  int statsize;

  void initialize(PamCoupler &coupler, ModelParameters &params, Parallel &par,
                  const Geometry<Straight> &primal_geom,
                  const Geometry<Twisted> &dual_geom) {

    std::string inFile = coupler.get_option<std::string>("standalone_input_file");
    YAML::Node config = YAML::LoadFile(inFile);        
    real simTime = config["simTime" ].as<real>(0.0_fp);
    real gcm_physics_dt = config["gcm_physics_dt" ].as<real>();
    real dycore_stat_freq = config["dycore_stat_freq" ].as<real>();
    if (simTime <= 0) {simTime = gcm_physics_dt;}
     
    this->primal_geometry = primal_geom;
    this->dual_geometry = dual_geom;
    this->nens = params.nens;
    this->masterproc = par.masterproc;
    this->statsize = simTime / dycore_stat_freq + 1;
  }

  virtual void compute(FieldSet<nprognostic> &progvars,
                       FieldSet<nconstant> &constvars, int tind){};
};
