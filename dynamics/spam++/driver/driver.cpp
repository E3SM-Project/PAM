

#include "common.h"
#include "RKSimple.h"
#include "variable_sets.h"
//#include "util.h"
#include "fileio.h"
#include "topology.h"
#include "geometry.h"
#include <cmath>
#include <iostream>
#include "mpi.h"
#include "params.h"
#include "model.h"

// *********************** //

int main(int argc, char** argv) {

  yakl::init();

  {
    Stats<number_of_dims, nprognostic, nconstant, nstats> stats;
    VariableSet<number_of_dims, nprognostic> prognostic_vars;
    VariableSet<number_of_dims, nconstant> constant_vars;
    VariableSet<number_of_dims, ndiagnostic> diagnostic_vars;
    VariableSet<number_of_dims, nauxiliary> auxiliary_vars;
    ExchangeSet<number_of_dims, nprognostic> prog_exchange;
    ExchangeSet<number_of_dims, nconstant> const_exchange;
    ExchangeSet<number_of_dims, nauxiliary> aux_exchange;
    FileIO<number_of_dims, nprognostic, nconstant, ndiagnostic, nstats> io;
    Tendencies<number_of_dims, nprognostic, nconstant, nauxiliary> tendencies;
    Diagnostics<number_of_dims, nprognostic, nconstant, ndiagnostic> diagnostics;
    Topology<number_of_dims> topology;
    ModelParameters params;
    Parallel par;

// SETTING THIS STUFF AT COMPILE TIME DOESN'T WORK BECAUSE OF SCOPING. IS THERE A WAY TO PROMOTE VARIABLES OUT OF IF STATEMENTS?

    //if (time_type == TIME_TYPE::KGRK)
    //{
      RKSimpleTimeIntegrator<number_of_dims, nprognostic, nconstant, nauxiliary, n_time_stages> tint;
    //}

    //if (geom_type == GEOM_TYPE::UNIFORM_RECT)
    //{
    UniformRectangularGeometry<number_of_dims,ic_quad_pts,ic_quad_pts,ic_quad_pts> ic_geometry;
    UniformRectangularGeometry<number_of_dims,1,1,1> tendencies_geometry;
    //}

    // Initialize MPI
    int ierr = MPI_Init( &argc , &argv );
    ierr = MPI_Comm_size(MPI_COMM_WORLD,&par.nranks);
    ierr = MPI_Comm_rank(MPI_COMM_WORLD,&par.myrank);

    // Read the parameters
    // Default input file is "input.txt" unless the user passes in another file
    std::cout << "reading parameters\n" << std::flush;
    std::string inFile = "input.txt";
    if (argc > 1) inFile = argv[1];
    readParamsFile<number_of_dims>(inFile, params, par);
    set_model_specific_params<number_of_dims>(inFile, params);
    std::cout << "read parameters\n" << std::flush;

    // Initialize the grid
    std::cout << "start init topo/geom\n" << std::flush;
    topology.initialize(par);
    ic_geometry.initialize(topology, params);
    tendencies_geometry.initialize(topology, params);
    std::cout << "finish init topo/geom\n" << std::flush;

    // Allocate the wind variable and the advected quantities
    std::cout << "start init variable sets and exchange sets\n" << std::flush;
// PROPERLY THESE CAN BE OF SIZE NDIMS+2...
    SArray<int, nprognostic, 4> prog_ndofs_arr;
    SArray<int, nconstant, 4> const_ndofs_arr;
    SArray<int, ndiagnostic, 4> diag_ndofs_arr;
    SArray<int, nauxiliary, 4> aux_ndofs_arr;
    std::array<std::string, nprognostic> prog_names_arr;
    std::array<std::string, nconstant> const_names_arr;
    std::array<std::string, ndiagnostic> diag_names_arr;
    std::array<std::string, nauxiliary> aux_names_arr;
    std::array<const Topology<number_of_dims> *, nprognostic> prog_topo_arr;
    std::array<const Topology<number_of_dims> *, nconstant> const_topo_arr;
    std::array<const Topology<number_of_dims> *, ndiagnostic> diag_topo_arr;
    std::array<const Topology<number_of_dims> *, nauxiliary> aux_topo_arr;
    initialize_variables<number_of_dims, nprognostic, nconstant, nauxiliary, ndiagnostic>(topology,
        prog_ndofs_arr, const_ndofs_arr, aux_ndofs_arr, diag_ndofs_arr,
        prog_names_arr, const_names_arr, aux_names_arr, diag_names_arr,
        prog_topo_arr, const_topo_arr, aux_topo_arr, diag_topo_arr);

    prognostic_vars.initialize("x", prog_names_arr, prog_topo_arr, prog_ndofs_arr);
    constant_vars.initialize("cons", const_names_arr, const_topo_arr, const_ndofs_arr);
    diagnostic_vars.initialize("diag", diag_names_arr, diag_topo_arr, diag_ndofs_arr);
    auxiliary_vars.initialize("aux", aux_names_arr, aux_topo_arr, aux_ndofs_arr);
    prog_exchange.initialize(prog_topo_arr, prog_ndofs_arr);
    const_exchange.initialize(const_topo_arr, const_ndofs_arr);
    aux_exchange.initialize(aux_topo_arr, aux_ndofs_arr);
    std::cout << "finish init variable sets\n" << std::flush;

    // Initialize the statistics
    stats.initialize(params, par, topology, tendencies_geometry);

    // Create the outputter
    std::cout << "start io init\n" << std::flush;
    io.initialize(params.outputName, topology, par, prognostic_vars, constant_vars, diagnostic_vars, stats);
    std::cout << "end io init\n" << std::flush;

    // set the initial conditions
    std::cout << "start set initial conditions\n" << std::flush;
    set_initial_conditions<nprognostic, nconstant, ic_quad_pts, ic_quad_pts, ic_quad_pts>(params, prognostic_vars, constant_vars, ic_geometry);
    // Do a boundary exchange
    prog_exchange.exchange_variable_set(prognostic_vars);
    const_exchange.exchange_variable_set(constant_vars);
    std::cout << "end set initial conditions\n" << std::flush;

    // Initialize the time stepper
    std::cout << "start ts init\n" << std::flush;
    tendencies.initialize(topology, tendencies_geometry, aux_exchange);
    diagnostics.initialize(topology, tendencies_geometry);
    tint.initialize(tendencies, prognostic_vars, constant_vars, auxiliary_vars, prog_exchange);
// SPECIFIC TO RK SCHEMES...
// SHOULD PROBABLY BE PART OF TINT INITIALIZE?
    set_stage_coefficients<n_time_stages>(time_type, tint.stage_coeffs);
    std::cout << "end ts init\n" << std::flush;

    // Output the initial model state
    std::cout << "start initial output\n" << std::flush;
    diagnostics.compute_diag(constant_vars, prognostic_vars, diagnostic_vars);
    io.outputInit(params.etime);
    std::cout << "end initial output\n" << std::flush;

    // Time stepping loop
    std::cout << "start timestepping loop\n" << std::flush;
    stats.compute(prognostic_vars, constant_vars, 0);
    for (uint nstep = 1; nstep<params.Nsteps+1; nstep++) {

      yakl::fence();
      tint.stepForward(params.dt);
      yakl::fence();

      params.etime += params.dt;
      if (nstep%params.Nout == 0) {
        std::cout << "step " << nstep << " time " << params.etime << "\n";
        diagnostics.compute_diag(constant_vars, prognostic_vars, diagnostic_vars);
        io.output(nstep, params.etime);
      }

      if (nstep%params.Nstat == 0)
      {
        stats.compute(prognostic_vars, constant_vars, nstep/params.Nstat);
      }

    }
    io.outputStats(stats);

    std::cout << "end timestepping loop\n" << std::flush;

    std::cout << "start io close\n" << std::flush;
    io.close();
    std::cout << "end io close\n" << std::flush;

   }

  int ierr = MPI_Finalize();
}
