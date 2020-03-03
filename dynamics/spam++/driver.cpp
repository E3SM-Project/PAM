

#include "driver.h"


// *********************** //

int main(int argc, char** argv) {

  yakl::init();

  {

    VariableSet<nprognostic> prognostic_vars;
    VariableSet<nconstant> constant_vars;
    VariableSet<ndiagnostic> diagnostic_vars;
    FileIO io;
    Tendencies tendencies;

// EVENTUALLY THESE TYPES SHOULD BE SETTABLE AT COMPILE TIME...
    RKSimpleTimeIntegrator<n_time_stages> tint;
    // MAKE SURE IT COLLAPSES PROPERLY IN 1D/2D
    PeriodicTopology<ndims,maxhalosize_x,maxhalosize_y,maxhalosize_z> topology;
    UniformRectangularGeometry<ndims,ic_quad_pts_x,ic_quad_pts_y,ic_quad_pts_z> geometry;


    // Read the parameters
    // Default input file is "input.txt" unless the user passes in another file
    std::string inFile = "input.txt";
    if (argc > 1) inFile = argv[1];
    readParamsFile(inFile);

    // Initialize the grid
    topology.initialize(nx, ny, nz);
    geometry.initialize(topology, xlen, xc, ylen, yc, zlen, zc);

    // Allocate the wind variable and the advected quantities
    initialize_variables(prognostic_vars, constant_vars, diagnostic_vars);

    // Create the outputter
    io.initialize();

    // set the initial conditions
    set_initial_conditions(prognostic_vars, constant_vars, diagnostic_vars);

    // Initialize the time stepper
    tendencies.initialize();
    tint.initialize(tendencies, prognostic_vars, topology, constant_vars, diagnostic_vars);
//PROPER WAY OF TREATING THIS?
    tinit.stage_coeffs = { 0.25, 1./3., 1./2., 1.};

    // Output the initial model state
    io.outputInit(prognostic_vars, constant_vars);

    // Time stepping loop
    time = 0.0;
    for (int nstep = 0; nstep<Nsteps; nstep++) {

      yakl::fence();
      tint.stepForward(params);
      yakl::fence();

      time += dt;
      if (nstep%Nout == 0) {
// UNSAFE IN PARALLEL
        std::cout << "step" << nstep << "time" << time << "\n";
        io.output(prognostic_vars, nstep, time);
      }

    }
    io.close();
   }

}
