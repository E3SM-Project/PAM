

#include "common.h"
#include "RKSimple.h"
#include "SSPRK.h"
#include "variable_sets.h"
#include "fileio.h"
#include "topology.h"
#include "geometry.h"
#include <cmath>
#include <iostream>
#include "mpi.h"
#include "params.h"
#include "model.h"
#include "initial_conditions.h"

// *********************** //

int main(int argc, char** argv) {

  yakl::init();

  {
    Stats<nprognostic, nconstant, nstats> stats;
    VariableSet<nprognostic> prognostic_vars;
    VariableSet<nconstant> constant_vars;
    VariableSet<ndiagnostic> diagnostic_vars;
    VariableSet<nauxiliary> auxiliary_vars;
    ExchangeSet<nprognostic> prog_exchange;
    ExchangeSet<nconstant> const_exchange;
    ExchangeSet<nauxiliary> aux_exchange;
    FileIO<nprognostic, nconstant, ndiagnostic, nstats> io;
    Tendencies<nprognostic, nconstant, nauxiliary> tendencies;
    Diagnostics<nprognostic, nconstant, ndiagnostic> diagnostics;
    Topology primal_topology;
    Topology dual_topology;
    ModelParameters params;
    Parallel par;
    
#if _TIME_TYPE==0
      RKSimpleTimeIntegrator<nprognostic, nconstant, nauxiliary, 4> tint;
#endif
#if _TIME_TYPE==1
      SSPKKTimeIntegrator<nprognostic, nconstant, nauxiliary, 3> tint;
#endif
#if _TIME_TYPE==2
      SSPKKTimeIntegrator<nprognostic, nconstant, nauxiliary, 2> tint;
#endif


#if _GEOM_TYPE==0   
    UniformRectangularStraightGeometry<ndims,ic_quad_pts,ic_quad_pts,ic_quad_pts> ic_primal_geometry;
    UniformRectangularTwistedGeometry<ndims,ic_quad_pts,ic_quad_pts,ic_quad_pts> ic_dual_geometry;
    UniformRectangularStraightGeometry<ndims,1,1,1> tendencies_primal_geometry;
    UniformRectangularTwistedGeometry<ndims,1,1,1> tendencies_dual_geometry;
#endif

    // Initialize MPI
    int ierr = MPI_Init( &argc , &argv );
    ierr = MPI_Comm_size(MPI_COMM_WORLD,&par.nranks);
    ierr = MPI_Comm_rank(MPI_COMM_WORLD,&par.myrank);

    // Read the parameters
    // Default input file is "input.txt" unless the user passes in another file
    std::cout << "reading parameters\n" << std::flush;
    std::string inFile = "input.txt";
    if (argc > 1) inFile = argv[1];
    readParamsFile(inFile, params, par);
    set_model_specific_params(inFile, params);
    set_ic_specific_params(inFile, params);
    std::cout << "read parameters\n" << std::flush;

    // Initialize the grid
    std::cout << "start init topo/geom\n" << std::flush;
    primal_topology.initialize(par,true);
    ic_primal_geometry.initialize(primal_topology, params);
    tendencies_primal_geometry.initialize(primal_topology, params);
//SHOULD REALLY GENERATE THIS FROM PRIMAL STUFF!
//THIS IS BROKEN FOR EXTRUDED MESHES...
    dual_topology.initialize(par,false);
    ic_dual_geometry.initialize(dual_topology, params);
    tendencies_dual_geometry.initialize(dual_topology, params);
    std::cout << "finish init topo/geom\n" << std::flush;

    // Allocate the variables
    std::cout << "start init variable sets and exchange sets\n" << std::flush;
    //this gives basedof, extdof and ndofs
    SArray<int, nprognostic, 3> prog_dofs_arr;
    SArray<int, nconstant, 3> const_dofs_arr;
    SArray<int, ndiagnostic, 3> diag_dofs_arr;
    SArray<int, nauxiliary, 3> aux_dofs_arr;
    std::array<std::string, nprognostic> prog_names_arr;
    std::array<std::string, nconstant> const_names_arr;
    std::array<std::string, ndiagnostic> diag_names_arr;
    std::array<std::string, nauxiliary> aux_names_arr;
    std::array<const Topology *, nprognostic> prog_topo_arr;
    std::array<const Topology *, nconstant> const_topo_arr;
    std::array<const Topology *, ndiagnostic> diag_topo_arr;
    std::array<const Topology *, nauxiliary> aux_topo_arr;
    initialize_variables<nprognostic, nconstant, nauxiliary, ndiagnostic>(primal_topology, dual_topology,
        prog_dofs_arr, const_dofs_arr, aux_dofs_arr, diag_dofs_arr,
        prog_names_arr, const_names_arr, aux_names_arr, diag_names_arr,
        prog_topo_arr, const_topo_arr, aux_topo_arr, diag_topo_arr);

    prognostic_vars.initialize("x", prog_names_arr, prog_topo_arr, prog_dofs_arr);
    constant_vars.initialize("cons", const_names_arr, const_topo_arr, const_dofs_arr);
    diagnostic_vars.initialize("diag", diag_names_arr, diag_topo_arr, diag_dofs_arr);
    auxiliary_vars.initialize("aux", aux_names_arr, aux_topo_arr, aux_dofs_arr);
    prog_exchange.initialize(prog_topo_arr, prog_dofs_arr);
    const_exchange.initialize(const_topo_arr, const_dofs_arr);
    aux_exchange.initialize(aux_topo_arr, aux_dofs_arr);
    std::cout << "finish init variable sets\n" << std::flush;

    // Initialize the statistics
    stats.initialize(params, par, primal_topology, dual_topology, tendencies_primal_geometry, tendencies_dual_geometry);

    // Create the outputter
    std::cout << "start io init\n" << std::flush;
    io.initialize(params.outputName, primal_topology, dual_topology, par, prognostic_vars, constant_vars, diagnostic_vars, stats);
    std::cout << "end io init\n" << std::flush;

    // set the initial conditions
    std::cout << "start set initial conditions\n" << std::flush;
    set_initial_conditions<nprognostic, nconstant, ic_quad_pts, ic_quad_pts, ic_quad_pts>(params, prognostic_vars, constant_vars, ic_primal_geometry, ic_dual_geometry);
    // Do a boundary exchange
    prog_exchange.exchange_variable_set(prognostic_vars);
    const_exchange.exchange_variable_set(constant_vars);
    std::cout << "end set initial conditions\n" << std::flush;

    // Initialize the time stepper, tendencies, diagnostics
    std::cout << "start ts init\n" << std::flush;
    tendencies.initialize(params, primal_topology, dual_topology, tendencies_primal_geometry, tendencies_dual_geometry, aux_exchange, const_exchange);
    tendencies.compute_constants(constant_vars, prognostic_vars);
    diagnostics.initialize(primal_topology, dual_topology, tendencies_primal_geometry, tendencies_dual_geometry);
    tint.initialize(tendencies, prognostic_vars, constant_vars, auxiliary_vars, prog_exchange);
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
