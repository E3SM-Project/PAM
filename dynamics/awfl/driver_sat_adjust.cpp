
#include "const.h"
#include "Spatial_euler3d_cons_expl_cart_fv_Agrid.h"
#include "Temporal_ader.h"
#include "Profiles.h"
#include "PhysicsSaturationAdjustment.h"

// Define the Spatial operator based on constants from the Temporal operator
typedef Spatial_euler3d_cons_expl_cart_fv_Agrid<nTimeDerivs,timeAvg,nAder> Spatial;

// Define the Temporal operator based on the Spatial operator
typedef Temporal_ader<Spatial> Temporal;

typedef PhysicsSaturationAdjustment Physics;

int main(int argc, char** argv) {
  yakl::init();
  {

    if (argc <= 1) { endrun("ERROR: Must pass the input YAML filename as a parameter"); }
    std::string inFile(argv[1]);
    YAML::Node config = YAML::LoadFile(inFile);
    if ( !config            ) { endrun("ERROR: Invalid YAML input file"); }
    if ( !config["simTime"] ) { endrun("ERROR: no simTime entry"); }
    if ( !config["outFreq"] ) { endrun("ERROR: no outFreq entry"); }
    real simTime = config["simTime"].as<real>();
    real outFreq = config["outFreq"].as<real>();
    int numOut = 0;

    // Create the model and the physics
    Temporal model;
    Physics physics;

    // Initialize the model and the physics
    model        .init( inFile , physics.numTracers );
    physics      .init( inFile );

    // Define the tracers
    physics.addTracers(model.spaceOp);

    // Initialize the dry state
    Spatial::StateArr state = model.spaceOp.createStateArr();
    model.spaceOp.initState(state, physics);

    // Initialize the tracers
    Spatial::TracerArr tracers = model.spaceOp.createTracerArr();
    physics.initTracers(model.spaceOp,tracers);

    // Adjust the model state to account for moisture
    model.spaceOp.adjustStateForMoisture(state,tracers,physics);

    real etime = 0;

    model.spaceOp.output( state , tracers , physics , etime );
    
    while (etime < simTime) {
      real dt = model.spaceOp.computeTimeStep(0.8, state, tracers, physics);
      if (etime + dt > simTime) { dt = simTime - etime; }
      model.timeStep( state , tracers , physics , dt );
      etime += dt;
      if (etime / outFreq >= numOut+1) {
        std::cout << "Etime , dt: " << etime << " , " << dt << "\n";
        model.spaceOp.output( state , tracers , physics , etime );
        numOut++;
      }
    }

    std::cout << "Elapsed Time: " << etime << "\n";

    model.finalize( state , tracers );

  }
  yakl::finalize();
}



