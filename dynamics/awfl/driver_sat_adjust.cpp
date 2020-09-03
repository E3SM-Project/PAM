
#include "const.h"
#include "Spatial_euler3d_cons_expl_cart_fv_Agrid.h"
#include "Temporal_ader.h"
#include "Profiles.h"
#include "Microphysics_saturation_adjustment.h"
#include "DataManager.h"

// Define the Spatial operator based on constants from the Temporal operator
typedef Spatial_euler3d_cons_expl_cart_fv_Agrid<nTimeDerivs,timeAvg,nAder> Spatial;

// Define the Temporal operator based on the Spatial operator
typedef Temporal_ader<Spatial> Dycore;

typedef Microphysics_saturation_adjustment Microphysics;

int main(int argc, char** argv) {
  yakl::init();
  {

    DataManager dm;

    if (argc <= 1) { endrun("ERROR: Must pass the input YAML filename as a parameter"); }
    std::string inFile(argv[1]);
    YAML::Node config = YAML::LoadFile(inFile);
    if ( !config            ) { endrun("ERROR: Invalid YAML input file"); }
    if ( !config["simTime"] ) { endrun("ERROR: no simTime entry"); }
    if ( !config["outFreq"] ) { endrun("ERROR: no outFreq entry"); }
    real simTime = config["simTime"].as<real>();
    real outFreq = config["outFreq"].as<real>();
    int numOut = 0;

    // Create the dycore and the microphysics
    Dycore       dycore;
    Microphysics micro;

    // Initialize the dycore and the microphysics
    dycore.init( inFile , micro.num_tracers , dm );
    micro .init( inFile , dycore.spaceOp    , dm );

    // Initialize the dry state
    dycore.spaceOp.init_state( dm , micro );

    // Initialize the tracers
    micro.init_tracers( dycore.spaceOp , dm );

    // Adjust the dycore state to account for moisture
    dycore.spaceOp.adjust_state_for_moisture( dm , micro );

    real etime = 0;

    dycore.spaceOp.output( dm , micro , etime );
    
    while (etime < simTime) {
      real dt = dycore.spaceOp.computeTimeStep( 0.8 , dm , micro );
      if (etime + dt > simTime) { dt = simTime - etime; }
      // dycore.timeStep( dm , micro , dt );
      // etime += dt;
      // if (etime / outFreq >= numOut+1) {
        std::cout << "Etime , dt: " << etime << " , " << dt << "\n";
      //   dycore.spaceOp.output( dm , micro , etime );
      //   numOut++;
      // }
    }

    std::cout << "Elapsed Time: " << etime << "\n";

    dycore.finalize( dm );

  }
  yakl::finalize();
}



