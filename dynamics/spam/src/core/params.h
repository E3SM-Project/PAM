#pragma once

#include "common.h"
#include "pam_coupler.h"
#include "parallel.h"

namespace pamc {

class Parameters {
public:
  int nx_glob = -1;
  int ny_glob = -1;
  int nz_dual = -1;
  int nens = -1;

  real out_freq = -1.;
  real stat_freq = -1.;
  real dtcrm = -1.;
  real sim_time = -1.;
  real dt_crm_phys = -1.;
  int crm_per_phys = -1;
  int statSize = -1;
  std::string outputName;
  std::string tstype;
  real si_tolerance = -1;
  int si_monitor_convergence;
  int si_verbosity_level;
  int si_max_iters;
  int si_nquad;
  bool si_two_point_discrete_gradient;

  real tanh_upwind_coeff = -1;

  real xlen, ylen;
  real xc, yc;

  int masterproc;
  bool inner_mpi;

  bool couple_wind = true;
  // solve a system to exactly invert the velocity averaging done
  // during conversion to coupler state when coupling winds
  bool couple_wind_exact_inverse = false;
  bool clip_negative_densities = true;

  bool clip_vertical_velocities = false;
  bool adjust_crm_per_phys_using_vert_cfl = false;
  real target_cfl = 0.7;
  real max_w = 50.0;
};

void readParamsFile(std::string inFile, Parameters &params, Parallel &par,
                    int nz) {
#ifdef PAM_STANDALONE
  // Read config file
  YAML::Node config = YAML::LoadFile(inFile);

  params.inner_mpi = config["inner_mpi"].as<bool>(false);
  params.nx_glob = config["crm_nx"].as<int>();
  params.ny_glob = config["crm_ny"].as<int>();
  par.nprocx = config["nprocx"].as<int>();
  par.nprocy = config["nprocy"].as<int>();
  params.nens = config["nens"].as<int>();
  params.sim_time = config["sim_time"].as<real>();
  params.dt_crm_phys = config["dt_crm_phys"].as<real>();
  params.crm_per_phys = config["crm_per_phys"].as<int>();
  params.out_freq = config["out_freq"].as<real>(-1.);
  params.stat_freq = config["stat_freq"].as<real>(-1.);
  params.tstype = config["tstype"].as<std::string>();
  params.si_tolerance = config["si_tolerance"].as<real>(1e-8);
  params.si_monitor_convergence = config["si_monitor_convergence"].as<int>(2);
  params.si_verbosity_level =
      config["si_verbosity_level"].as<int>(params.si_monitor_convergence);
  params.si_max_iters = config["si_max_iters"].as<int>(
      params.si_monitor_convergence > 1 ? 50 : 5);
  params.si_nquad = config["si_nquad"].as<int>(4);
  params.si_two_point_discrete_gradient =
      config["si_two_point_discrete_gradient"].as<bool>(false);
  params.tanh_upwind_coeff = config["tanh_upwind_coeff"].as<real>(250);
  params.outputName = config["dycore_out_prefix"].as<std::string>("output");
  params.nz_dual = nz;

  params.couple_wind = config["couple_wind"].as<bool>(true);
  params.couple_wind_exact_inverse =
      config["couple_wind_exact_inverse"].as<bool>(false);

  params.clip_negative_densities =
      config["clip_negative_densities"].as<bool>(true);

  params.clip_vertical_velocities =
      config["clip_vertical_velocities"].as<bool>(false);

  params.adjust_crm_per_phys_using_vert_cfl =
      config["adjust_crm_per_phys_using_vert_cfl"].as<bool>(false);

  params.target_cfl = config["max_w"].as<real>(0.7);
  params.max_w = config["max_w"].as<real>(50.0);

  // ADD A CHECK HERE THAT TOTAL TIME IS EXACTLY DIVISIBLE BY STAT_FREQ
  if (params.stat_freq >= 0.) {
    params.statSize = params.sim_time / params.stat_freq;
  } else {
    params.statSize = 0;
  }

  params.dtcrm = params.dt_crm_phys / params.crm_per_phys;
#endif
}

void read_params_coupler(Parameters &params, Parallel &par,
                         pam::PamCoupler &coupler) {

  params.inner_mpi = false;
  params.nx_glob = coupler.get_nx();
  params.ny_glob = coupler.get_ny();
  par.nprocx = 1;
  par.nprocy = 1;
  params.nens = coupler.get_nens();
  // THESE NAMES IN COUPLER DONT CORRESPOND WITH CONFIG FILE NAMES
  // FIX THIS
  params.dt_crm_phys = coupler.get_option<real>("crm_dt");
  if (coupler.option_exists("crm_dyn_per_phys")) {
    params.crm_per_phys = coupler.get_option<int>("crm_dyn_per_phys");
  } else {
    params.crm_per_phys = 1;
  }

  // Allow coupler options to override defaults
  if (coupler.option_exists("spam_couple_wind_exact_inverse")) {
    params.couple_wind_exact_inverse = coupler.get_option<bool>("spam_couple_wind_exact_inverse");
  }
  if (coupler.option_exists("spam_clip_negative_densities")) {
    params.clip_negative_densities = coupler.get_option<bool>("spam_clip_negative_densities");
  }
  if (coupler.option_exists("spam_clip_vertical_velocities")) {
    params.clip_vertical_velocities = coupler.get_option<bool>("spam_clip_vertical_velocities");
  }
  if (coupler.option_exists("spam_adjust_crm_per_phys_using_vert_cfl")) {
    params.adjust_crm_per_phys_using_vert_cfl = coupler.get_option<bool>("spam_adjust_crm_per_phys_using_vert_cfl");
  }
  if (coupler.option_exists("spam_target_cfl")) {
    params.target_cfl = coupler.get_option<real>("spam_target_cfl");
  }
  if (coupler.option_exists("spam_max_w")) {
    params.max_w = coupler.get_option<real>("spam_max_w");
  }

#ifdef PAMC_MAN
  params.tstype = "ssprk3";
#else
  params.tstype = "si";
#endif
  params.si_tolerance = 1e-8;
  params.si_monitor_convergence = 0;
  params.si_verbosity_level = 0;
  params.si_max_iters = 3;
  params.si_nquad = 2;
  params.si_two_point_discrete_gradient = false;
  params.tanh_upwind_coeff = 250;
  params.outputName = "pamc_output";
  params.nz_dual = coupler.get_nz();
  params.statSize = 0;
  params.dtcrm = params.dt_crm_phys / params.crm_per_phys;
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

void check_and_print_parameters(const Parameters &params, const Parallel &par,
                                bool verbose = false) {

  // Print out the values
  if (par.masterproc and verbose) {
    std::cout << "nx:         " << params.nx_glob << "\n";
    std::cout << "ny:         " << params.ny_glob << "\n";
    std::cout << "nl dual:         " << params.nz_dual << "\n";
    std::cout << "ni dual:         " << params.nz_dual + 1 << "\n";
    std::cout << "nens:         " << params.nens << "\n";

    std::cout << "halox:      " << par.halox << "\n";
    std::cout << "haloy:      " << par.haloy << "\n";

    std::cout << "sim_time:         " << params.sim_time << "\n";
    std::cout << "dtcrm:         " << params.dtcrm << "\n";
    std::cout << "dt_crm_phys:         " << params.dt_crm_phys << "\n";
    std::cout << "statSize:     " << params.statSize << "\n";
    std::cout << "crm per phys:     " << params.crm_per_phys << "\n";
    std::cout << "out_freq:       " << params.out_freq << "\n";
    std::cout << "stat_freq:       " << params.stat_freq << "\n";
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
} // namespace pamc
