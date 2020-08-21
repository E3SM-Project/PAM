
#include "const.h"
#include "Spatial_euler3d_cons_expl_cart_fv_Agrid.h"
#include "Temporal_ssprk3.h"

// Define the Spatial operator based on constants from the Temporal operator
typedef Spatial_euler3d_cons_expl_cart_fv_Agrid<nTimeDerivs,timeAvg,nAder> Spatial;
// Define the Temporal operator based on the Spatial operator
typedef Temporal_ssprk3<Spatial> Temporal;

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

    Temporal model;

    model.spaceOp.addTracer("theta","replica of theta");

    model.init(inFile);

    Spatial::StateArr  state   = model.spaceOp.createStateArr ();
    Spatial::TracerArr tracers = model.spaceOp.createTracerArr();

    model.spaceOp.initState  (state  );
    model.spaceOp.initTracers(tracers);

    real etime = 0;

    model.spaceOp.output( state , tracers , etime );
    
    while (etime < simTime) {
      real dt = model.spaceOp.computeTimeStep(0.8,state);
      if (etime + dt > simTime) { dt = simTime - etime; }
      model.timeStep( state , tracers , dt );
      etime += dt;
      if (etime / outFreq >= numOut+1) {
        std::cout << "Etime , dt: " << etime << " , " << dt << "\n";
        model.spaceOp.output( state , tracers , etime );
        numOut++;
      }
    }

    std::cout << "Elapsed Time: " << etime << "\n";

    model.finalize( state , tracers );

  }
  yakl::finalize();
}



