

#include "stdlib.h"
#include <iostream>
#include <string>

#include "common.h"
#include "RKSimple.h"
#include "variable_set.h"
#include "topology.h"
#include "uniform_geometry.h"
#include "tendencies.h"
#include "finitevolume.h"
#include "divergence.h"
#include "driver.h"

// HOW DO WE MAKE THIS UNIVERSAL? VIA PRE-PROCESSOR FLAGS, OF COURSE...
#include "advection.h"

// *********************** //

int main(int argc, char** argv) {

  yakl::init();

  {
    PrognosticVars      prognostic_vars;
    ConstantVars      constant_vars;
    DiagnosticVars      diagnostic_vars;
    FileIO         io;
    ModelTendencies tendencies;

// EVENTUALLY THIS SHOULD BE SETTABLE AT COMPILE TIME...
    RKSimpleTimeIntegrator tint;
    PeriodicTopology    topology;
    UniformRectangularGeometry geometry;

    // Read the parameters
    // Default input file is "input.txt" unless the user passes in another file
    std::string inFile = "input.txt";
    if (argc > 1) inFile = argv[1];
    readParamsFile(inFile);

    // Initialize the grid
    topology.initialize(nx, ny, nz);
    geometry.initialize(topology, xlen, xc, ylen, yc, zlen, zc);

    // Allocate the wind variable and the advected quantities
    prognostic_vars.create(topology, true);
    constant_vars.create(topology, true);
    diagnostic_vars.create(topology, true);

    // Initialize the output
    io.initialize();

    // set the initial conditions
    set_initial_conditions(constant_vars, prognostic_vars);

    // Initialize the time stepper
    tint.initialize(tendencies, prognostic_vars, topology, constant_vars, diagnostic_vars, dt);

    // Output the initial model state
    io.outputInit(constant_vars, prognostic_vars);

    // Time stepping loop
    real time = 0.0;
    for (int nstep = 0; nstep<Nsteps; nstep++) {

      yakl::fence();
      tint.stepForward(params);
      yakl::fence();

      time += dt;
      if (nstep%Nout == 0) {
// UNSAFE IN PARALLEL
        std::cout << "step" << nstep << "time" << time << "\n";
        io.output(prognostic_vars);
      }

    }
   }

}
