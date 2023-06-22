#pragma once

#include "common.h"
#include "pam_coupler.h"
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
  int si_monitor_convergence;
  int si_verbosity_level;
  int si_max_iters;

  real tanh_upwind_coeff = -1;

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
  params.si_monitor_convergence = config["si_monitor_convergence"].as<int>(2);
  params.si_verbosity_level =
      config["si_verbosity_level"].as<int>(params.si_monitor_convergence);
  params.si_max_iters = config["si_max_iters"].as<int>(
      params.si_monitor_convergence > 1 ? 50 : 5);
  params.tanh_upwind_coeff = config["tanh_upwind_coeff"].as<real>(250);
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
#ifdef _MAN
  params.tstype = "ssprk3";
#else
  params.tstype = "si";
#endif
  params.si_tolerance = 1e-8;
  params.si_monitor_convergence = 0;
  params.si_verbosity_level = 0;
  params.si_max_iters = 3;
  params.tanh_upwind_coeff = 250;
  params.outputName = "pamc_output";
  params.nz_dual = coupler.get_nz();
  params.Nsteps = params.simSteps * params.crm_per_phys;

  params.dtcrm = params.dtphys / params.crm_per_phys;
}

void finalize_parallel(Parameters &params, Parallel &par) {
  int ierr = MPI_Comm_rank(MPI_COMM_WORLD, &par.actualrank);
  par.masterproc = par.actualrank == 0;
  params.masterproc = par.masterproc;
  if (params.inner_mpi) {
    ierr = MPI_Comm_size(MPI_COMM_WORLD, &par.nranks);
    ierr = MPI_Comm_rank(MPI_COMM_WORLD, &par.myrank);
  } else {
    par.nranks = 1;
    par.myrank = 0;
  }

  // if (!(par.nprocx * par.nprocy == par.nranks)) {endrun("Error: nranks !=
  // nprocx * nprocy");}

  // Get my process grid IDs
  par.py = floor(par.myrank / par.nprocx);
  par.px = par.myrank - par.nprocx * par.py;

  // Get my beginning and ending global indices; and domain sizes
  double nper;
  nper = ((double)params.nx_glob) / par.nprocx;
  par.i_beg = (int)round(nper * par.px);
  par.i_end = (int)round(nper * (par.px + 1)) - 1;
  par.nx = par.i_end - par.i_beg + 1;
  par.nx_glob = params.nx_glob;

  if (ndims >= 2) {
    nper = ((double)params.ny_glob) / par.nprocy;
    par.j_beg = (int)round(nper * par.py);
    par.j_end = (int)round(nper * (par.py + 1)) - 1;
    par.ny = par.j_end - par.j_beg + 1;
    par.ny_glob = params.ny_glob;
  }

  par.nz = params.nz_dual;
  par.nens = params.nens;

  // Determine my neighbors
  // x-dir
  par.x_neigh(0) = ij_to_l(wrap(par.px - 1, par.nprocx), par.py, par.nprocx);
  par.x_neigh(1) = ij_to_l(wrap(par.px + 1, par.nprocx), par.py, par.nprocx);

  // y-dir
  par.y_neigh(0) = -1;
  par.y_neigh(1) = -1;
  if (ndims == 2) {
    par.y_neigh(0) = ij_to_l(par.px, wrap(par.py - 1, par.nprocy), par.nprocx);
    par.y_neigh(1) = ij_to_l(par.px, wrap(par.py + 1, par.nprocy), par.nprocx);

    par.ll_neigh = ij_to_l(wrap(par.px - 1, par.nprocx),
                           wrap(par.py - 1, par.nprocy), par.nprocx);
    par.lr_neigh = ij_to_l(wrap(par.px + 1, par.nprocx),
                           wrap(par.py - 1, par.nprocy), par.nprocx);
    par.ur_neigh = ij_to_l(wrap(par.px + 1, par.nprocx),
                           wrap(par.py + 1, par.nprocy), par.nprocx);
    par.ul_neigh = ij_to_l(wrap(par.px - 1, par.nprocx),
                           wrap(par.py + 1, par.nprocy), par.nprocx);
  }

  // set boundaries
  par.xbnd = BND_TYPE::PERIODIC;
  par.ybnd = BND_TYPE::PERIODIC;

  // set halos
  par.halox = maxhalosize;
  par.haloy = maxhalosize;

// Debug output for the parallel decomposition
#ifdef PAM_DEBUG
  ierr = MPI_Barrier(MPI_COMM_WORLD);
  for (int rr = 0; rr < par.nranks; rr++) {
    if (rr == par.myrank) {
      std::cout << "Hello! My rank is: " << par.myrank << "\n";
      std::cout << "My proc grid ID is: " << par.px << " , " << par.py << "\n";
      std::cout << "I have: " << par.nx << " x " << par.ny << " grid points."
                << "\n";
      std::cout << "I start at index: " << par.i_beg << " x " << par.j_beg
                << "\n";
      std::cout << "I end at index: " << par.i_end << " x " << par.j_end
                << "\n";
      std::cout << "My x neighbors are: " << par.x_neigh(0) << " "
                << par.x_neigh(1) << "\n";
      std::cout << "My y neighbors are: " << par.y_neigh(0) << " "
                << par.y_neigh(1) << "\n";
      std::cout << "My corner neighbors are: " << par.ll_neigh << " "
                << par.ul_neigh << " " << par.ur_neigh << " " << par.lr_neigh
                << "\n";
    }
    ierr = MPI_Barrier(MPI_COMM_WORLD);
  }
  ierr = MPI_Barrier(MPI_COMM_WORLD);
#endif
}

void check_and_print_parameters(const Parameters &params, const Parallel &par) {
  // Check time stepping params
  if (not(params.dtphys > 0.0_fp) or not(params.simSteps > 0) or
      not(params.crm_per_phys > 0) or not(params.Nout > 0) or
      not(params.Nstat > 0)) {
    endrun("spam++ must use step-based time control logic ie set simSteps >0, "
           "dtphys>0, crm_per_phys >0, outSteps >0, statSteps >0");
  }

  // Print out the values
  if (par.masterproc) {
    std::cout << "nx:         " << params.nx_glob << "\n";
    std::cout << "ny:         " << params.ny_glob << "\n";
    std::cout << "nl dual:         " << params.nz_dual << "\n";
    std::cout << "ni dual:         " << params.nz_dual + 1 << "\n";
    std::cout << "nens:         " << params.nens << "\n";

    std::cout << "halox:      " << par.halox << "\n";
    std::cout << "haloy:      " << par.haloy << "\n";

    std::cout << "dtcrm:         " << params.dtcrm << "\n";
    std::cout << "dtphys:         " << params.dtphys << "\n";
    std::cout << "Nsteps:     " << params.Nsteps << "\n";
    std::cout << "simSteps:     " << params.simSteps << "\n";
    std::cout << "crm per phys:     " << params.crm_per_phys << "\n";
    std::cout << "Nout:       " << params.Nout << "\n";
    std::cout << "outputName: " << params.outputName << "\n";

    std::cout << "nranks:     " << par.nranks << "\n";
    std::cout << "nprocx:     " << par.nprocx << "\n";
    std::cout << "nprocy:     " << par.nprocy << "\n";

    std::cout << "xlen:       " << params.xlen << "\n";
    std::cout << "ylen:       " << params.ylen << "\n";
    std::cout << "xc:         " << params.xc << "\n";
    std::cout << "yc:         " << params.yc << "\n";
  }
};
