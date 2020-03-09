

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

    VariableSet<number_of_dims, nprognostic> prognostic_vars;
    VariableSet<number_of_dims, nconstant> constant_vars;
    VariableSet<number_of_dims, ndiagnostic> diagnostic_vars;
    ExchangeSet<number_of_dims, nprognostic> prog_exchange;
    ExchangeSet<number_of_dims, nconstant> const_exchange;
    ExchangeSet<number_of_dims, ndiagnostic> diag_exchange;
    FileIO<number_of_dims, nprognostic, nconstant, ndiagnostic> io;
    Tendencies<number_of_dims, nprognostic, nconstant, ndiagnostic> tendencies;
    Topology<number_of_dims> topology;
    Parameters params;

// EVENTUALLY THESE TYPES SHOULD BE SETTABLE AT COMPILE TIME...
    RKSimpleTimeIntegrator<number_of_dims, nprognostic, nconstant, ndiagnostic, n_time_stages> tint;
    UniformRectangularGeometry<number_of_dims,ic_quad_pts,ic_quad_pts,ic_quad_pts> ic_geometry;
    UniformRectangularGeometry<number_of_dims,1,1,1> tendencies_geometry;

    // Initialize MPI
    int ierr = MPI_Init( &argc , &argv );
    ierr = MPI_Comm_size(MPI_COMM_WORLD,&params.nranks);
    ierr = MPI_Comm_rank(MPI_COMM_WORLD,&params.myrank);

    // Determine if I'm the master process
    if (params.myrank == 0) { params.masterproc = 1;}
    else { params.masterproc = 0; }

    // Read the parameters
    // Default input file is "input.txt" unless the user passes in another file
    std::cout << "reading parameters\n" << std::flush;
    std::string inFile = "input.txt";
    if (argc > 1) inFile = argv[1];
    set_model_specific_params(params);
    readParamsFile(inFile, params);
    std::cout << "read parameters\n" << std::flush;

    // Initialize the grid
    std::cout << "start init topo/geom\n" << std::flush;
    topology.initialize(params.nx, params.ny, params.nz, maxhalosize, maxhalosize, maxhalosize);
// THIS IS SPECIFIC TO UNIFORM RECT GEOMETRY...
    ic_geometry.initialize(topology, params.xlen/params.nx, params.ylen/params.ny, params.zlen/params.nz, params.xlen, params.ylen, params.zlen, params.xc, params.yc, params.zc);
    tendencies_geometry.initialize(topology, params.xlen/params.nx, params.ylen/params.ny, params.zlen/params.nz, params.xlen, params.ylen, params.zlen, params.xc, params.yc, params.zc);
    std::cout << "finish init topo/geom\n" << std::flush;

    // Allocate the wind variable and the advected quantities
    std::cout << "start init variable sets and exchange sets\n" << std::flush;
// PROPERLY THESE CAN BE OF SIZE NDIMS+2...
    SArray<int, nprognostic, 4> prog_ndofs_arr;
    SArray<int, nconstant, 4> const_ndofs_arr;
    SArray<int, ndiagnostic, 4> diag_ndofs_arr;
    std::array<std::string, nprognostic> prog_names_arr;
    std::array<std::string, nconstant> const_names_arr;
    std::array<std::string, ndiagnostic> diag_names_arr;
    std::array<const Topology<number_of_dims> *, nprognostic> prog_topo_arr;
    std::array<const Topology<number_of_dims> *, nconstant> const_topo_arr;
    std::array<const Topology<number_of_dims> *, ndiagnostic> diag_topo_arr;
    initialize_variables<number_of_dims, nprognostic, nconstant, ndiagnostic>(topology, prog_ndofs_arr, const_ndofs_arr, diag_ndofs_arr, prog_names_arr, const_names_arr, diag_names_arr, prog_topo_arr, const_topo_arr, diag_topo_arr);

    prognostic_vars.initialize("x", prog_names_arr, prog_topo_arr, prog_ndofs_arr);
    constant_vars.initialize("cons", const_names_arr, const_topo_arr, const_ndofs_arr);
    diagnostic_vars.initialize("diag", diag_names_arr, diag_topo_arr, diag_ndofs_arr);
    prog_exchange.initialize(prog_topo_arr, prog_ndofs_arr);
    const_exchange.initialize(const_topo_arr, const_ndofs_arr);
    diag_exchange.initialize(diag_topo_arr, diag_ndofs_arr);
    std::cout << "finish init variable sets\n" << std::flush;

    // Create the outputter
    std::cout << "start io init\n" << std::flush;
    io.initialize(params.outputName, topology, prognostic_vars, constant_vars);
    std::cout << "end io init\n" << std::flush;

    // set the initial conditions
    std::cout << "start set initial conditions\n" << std::flush;
    set_initial_conditions<nprognostic, nconstant, ndiagnostic, ic_quad_pts, ic_quad_pts, ic_quad_pts>(prognostic_vars, constant_vars, ic_geometry);
    // Do a boundary exchange
    prog_exchange.exchange_variable_set(prognostic_vars);
    const_exchange.exchange_variable_set(constant_vars);
    std::cout << "end set initial conditions\n" << std::flush;

    // Initialize the time stepper
    std::cout << "start ts init\n" << std::flush;
    tendencies.initialize(topology, tendencies_geometry, diag_exchange);
    tint.initialize(tendencies, prognostic_vars, constant_vars, diagnostic_vars, prog_exchange);
    set_stage_coefficients<n_time_stages>(time_type, tint.stage_coeffs);
    std::cout << "end ts init\n" << std::flush;

    // Output the initial model state
    std::cout << "start initial output\n" << std::flush;
    io.outputInit(params.etime);
    std::cout << "end initial output\n" << std::flush;

    // Time stepping loop
    std::cout << "start timestepping loop\n" << std::flush;
    for (uint nstep = 1; nstep<params.Nsteps+1; nstep++) {

      yakl::fence();
      tint.stepForward(params.dt);
      yakl::fence();

      params.etime += params.dt;
      if (nstep%params.Nout == 0) {
// UNSAFE IN PARALLEL
        std::cout << "step " << nstep << " time " << params.etime << "\n";
        io.output(nstep, params.etime);
      }

    }

    std::cout << "end timestepping loop\n" << std::flush;

    std::cout << "start io close\n" << std::flush;
    io.close();
    std::cout << "end io close\n" << std::flush;

   }

  int ierr = MPI_Finalize();
}
