#pragma once

#include "common.h"
#include "parallel.h"

class Parameters {
public:
  int nx_glob = -1;
  int ny_glob = -1;
  int nz_dual = -1;
  int nens = -1;

  int simSteps = -1;
  int Nsteps = -1;
  int Nout = -1;
  real dtcrm = -1.;
  real dtphys = -1.;
  int crm_per_phys = -1;
  int Nstat = -1;
  std::string outputName;
  std::string tstype;
  real si_tolerance = -1;

  real xlen, ylen;
  real xc, yc;

  int masterproc;
  bool inner_mpi;
};

void readParamsFile(std::string inFile, Parameters &params, Parallel &par,
                    int nz) {
  // Read config file
  YAML::Node config = YAML::LoadFile(inFile);

  params.inner_mpi = config["inner_mpi"].as<bool>(false);
  params.nx_glob = config["crm_nx"].as<int>();
  params.ny_glob = config["crm_ny"].as<int>();
  par.nprocx = config["nprocx"].as<int>();
  par.nprocy = config["nprocy"].as<int>();
  params.nens = config["nens"].as<int>();
  params.dtphys = config["dtphys"].as<real>();
  params.crm_per_phys = config["crm_per_phys"].as<int>(0);
  params.Nout = config["outSteps"].as<int>(0);
  params.Nstat = config["statSteps"].as<int>(0);
  params.simSteps = config["simSteps"].as<int>(0);
  params.tstype = config["tstype"].as<std::string>();
  params.si_tolerance = config["si_tolerance"].as<real>(1e-8);
  params.outputName = config["dycore_out_prefix"].as<std::string>("output");
  params.nz_dual = nz;
  params.Nsteps = params.simSteps * params.crm_per_phys;

  params.dtcrm = params.dtphys / params.crm_per_phys;
}

void read_params_coupler(Parameters &params, Parallel &par,
                         pam::PamCoupler &coupler) {

  params.inner_mpi = false;
  params.nx_glob = coupler.get_nx();
  params.ny_glob = coupler.get_ny();
  par.nprocx = 1;
  par.nprocy = 1;
  params.nens = coupler.get_nens();
  params.dtphys = coupler.get_option<real>("crm_dt");
  params.crm_per_phys = 1;

  params.Nout = 1;
  params.Nstat = 1;
  params.simSteps = 1;
  params.tstype = "ssprk3";
  params.si_tolerance = 1e-8;
  params.outputName = "pamc_output";
  params.nz_dual = coupler.get_nz();
  params.Nsteps = params.simSteps * params.crm_per_phys;

  params.dtcrm = params.dtphys / params.crm_per_phys;
}
