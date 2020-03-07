

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

// HOW DO WE MAKE THIS UNIVERSAL? VIA PRE-PROCESSOR FLAGS, OF COURSE...
#include "advection.h"

// WHERE SHOULD THESE REALLY LIVE?
// PROBABLY IN A PARAMETERS CLASS, WHICH SPECIFIC TO EACH MODEL...
// parameters
uint nx, ny, nz;
uint Nsteps, Nout;
real dt, etime;
std::string outputName;

// *********** COMPILE TIME CONSTANTS ************** //
// EVENTUALLY THESE SHOULD ALL BE COMPILE TIME/PRE-PROCESSOR FLAGS...
// WITH REASONABLE DEFAULTS!

// Number of Dimensions
uint constexpr number_of_dims = 1;

// Spatial order of accuracy for the model
uint constexpr reconstruction_order_x = 2;
uint constexpr reconstruction_order_y = 2;
uint constexpr reconstruction_order_z = 2;
uint constexpr differential_order_x = 2;
uint constexpr differential_order_y = 2;
uint constexpr differential_order_z = 2;

// Reconstruction type
RECONSTRUCTION_TYPE reconstruction_type_x = RECONSTRUCTION_TYPE::CFV;
RECONSTRUCTION_TYPE reconstruction_type_y = RECONSTRUCTION_TYPE::CFV;
RECONSTRUCTION_TYPE reconstruction_type_z = RECONSTRUCTION_TYPE::CFV;

// SHOULD MAYBE ALSO COLLAPSE IN 1D/2D?
// Halo sizes
uint maxhalosize_x = mymax(reconstruction_order_x,differential_order_x)/2; // IS THIS ALWAYS CORRECT?
uint maxhalosize_y = mymax(reconstruction_order_y,differential_order_y)/2; // IS THIS ALWAYS CORRECT?
uint maxhalosize_z = mymax(reconstruction_order_z,differential_order_z)/2; // IS THIS ALWAYS CORRECT?

// MAKE SURE NQUADPTS COLLAPSES PROPERLY IN 1D/2D

// initial condition quadrature pts
uint constexpr ic_quad_pts_x = 3;
uint constexpr ic_quad_pts_y = 3;
uint constexpr ic_quad_pts_z = 3;

// Time scheme
TIME_TYPE time_type = TIME_TYPE::KGRK;
uint constexpr n_time_stages = 4;

// *********************** //

// FIX THIS TO ACTUALLY READ A FILE, EVENTUALLY
void readParamsFile(std::string inFile) {
  // Topology
  nx = 20;
  ny = 100;
  nz = 10;

  // Time stepping
  dt = 0.5;
  Nsteps = 500;
  Nout = 5;

  outputName = "output.nc";

  // EVENTUALLY THIS SHOULD BE SET BY THE INITIAL CONDITION CHOICE!
  xlen = 50.0;
  ylen = 50.0;
  zlen = 50.0;
  xc = 25.0;
  yc = 25.0;
  zc = 25.0;
  etime = 0.0;
};


class Parallel {

public:

  int nranks;
  int myrank;
  int masterproc;
};

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
    AdvectionTendencies<number_of_dims, nprognostic, nconstant, ndiagnostic> tendencies;
    Topology<number_of_dims> topology;
    Parallel par;

// EVENTUALLY THESE TYPES SHOULD BE SETTABLE AT COMPILE TIME...
    RKSimpleTimeIntegrator<number_of_dims, nprognostic, nconstant, ndiagnostic, n_time_stages> tint;
    UniformRectangularGeometry<number_of_dims,ic_quad_pts_x,ic_quad_pts_y,ic_quad_pts_z> geometry;

    // Initialize MPI
    int ierr = MPI_Init( &argc , &argv );
    ierr = MPI_Comm_size(MPI_COMM_WORLD,&par.nranks);
    ierr = MPI_Comm_rank(MPI_COMM_WORLD,&par.myrank);

    // Determine if I'm the master process
    if (par.myrank == 0) { par.masterproc = 1;}
    else { par.masterproc = 0; }

    // Read the parameters
    // Default input file is "input.txt" unless the user passes in another file
    std::cout << "reading parameters\n" << std::flush;
    std::string inFile = "input.txt";
    if (argc > 1) inFile = argv[1];
    readParamsFile(inFile);
    std::cout << "read parameters\n" << std::flush;

    // Initialize the grid
    std::cout << "start init topo/geom\n" << std::flush;
    topology.initialize(nx, ny, nz, maxhalosize_x, maxhalosize_y, maxhalosize_z);
    geometry.initialize(topology, xlen/nx, ylen/ny, zlen/nz, xlen, ylen, zlen, xc, yc, zc);
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
    io.initialize(outputName, topology, prognostic_vars, constant_vars);
    std::cout << "end io init\n" << std::flush;

    // set the initial conditions
    set_initial_conditions<number_of_dims, nprognostic, nconstant, ndiagnostic, ic_quad_pts_x,ic_quad_pts_y,ic_quad_pts_z>(prognostic_vars, constant_vars, geometry);
    // Do a boundary exchange
    prog_exchange.exchange_variable_set(prognostic_vars);
    const_exchange.exchange_variable_set(constant_vars);

    // Initialize the time stepper
    std::cout << "start ts init\n" << std::flush;
    tendencies.initialize(topology, diag_exchange);
    tint.initialize(tendencies, prognostic_vars, constant_vars, diagnostic_vars, prog_exchange);
    set_stage_coefficients<n_time_stages>(time_type, tint.stage_coeffs);
    std::cout << "end ts init\n" << std::flush;

    // Output the initial model state
    std::cout << "start initial output\n" << std::flush;
    io.outputInit(etime);
    std::cout << "end initial output\n" << std::flush;

    // Time stepping loop
    std::cout << "start timestepping loop\n" << std::flush;
    for (uint nstep = 1; nstep<Nsteps+1; nstep++) {

      yakl::fence();
      tint.stepForward(dt);
      yakl::fence();

      etime += dt;
      if (nstep%Nout == 0) {
// UNSAFE IN PARALLEL
        std::cout << "step " << nstep << " time " << etime << "\n";
        io.output(nstep, etime);
      }

    }

    std::cout << "end timestepping loop\n" << std::flush;

    std::cout << "start io close\n" << std::flush;
    io.close();
    std::cout << "end io close\n" << std::flush;

   }

  int ierr = MPI_Finalize();
}
