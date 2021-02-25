
#include "const.h"
#include "Spatial_euler3d_cons_hevi1_cart_fv_Agrid.h"
#include "Temporal_ader.h"
#include "Profiles.h"
#include "Microphysics_kessler.h"
#include "DataManager.h"

// Define the Spatial operator based on constants from the Temporal operatora header file
typedef Spatial_operator<nTimeDerivs,timeAvg,nAder> Spatial;

// Define the Temporal operator based on the Spatial operator
typedef Temporal_operator<Spatial> Dycore;

typedef Microphysics Microphysics;


void allocate_coupler_state( std::string inFile , DataManager &dm);


int main(int argc, char** argv) {
  yakl::init();
  {

    DataManager dm;

    if (argc <= 1) { endrun("ERROR: Must pass the input YAML filename as a parameter"); }
    std::string inFile(argv[1]);
    YAML::Node config = YAML::LoadFile(inFile);
    if ( !config            ) { endrun("ERROR: Invalid YAML input file"); }
    real simTime = config["simTime"].as<real>();
    real outFreq = config["outFreq"].as<real>();
    int numOut = 0;

    allocate_coupler_state( inFile , dm );

    // Create the dycore and the microphysics
    Dycore       dycore;
    Microphysics micro;

    // Initialize the dycore and the microphysics
    dycore.init( inFile , micro.num_tracers , dm );
    micro .init( inFile , dycore            , dm );

    // Initialize the dry state
    dycore.init_state_and_tracers( dm , micro );

    real etime = 0;

    dycore.output( dm , micro , etime );
    
    while (etime < simTime) {
      real dt = dycore.compute_time_step( 0.8 , dm , micro );
      if (etime + dt > simTime) { dt = simTime - etime; }
      dycore.timeStep( dm , micro , dt );
      micro.timeStep( dm , dt );
      etime += dt;
      if (etime / outFreq >= numOut+1) {
        std::cout << "Etime , dt: " << etime << " , " << dt << "\n";
        dycore.output( dm , micro , etime );
        numOut++;
      }
    }

    dycore.timeStep( dm , micro , 1. );

    std::cout << "Elapsed Time: " << etime << "\n";

    dycore.finalize( dm );

  }
  yakl::finalize();
}



void allocate_coupler_state( std::string inFile , DataManager &dm) {
  YAML::Node config = YAML::LoadFile(inFile);
  if ( !config            ) { endrun("ERROR: Invalid YAML input file"); }
  real nx = config["nx"].as<real>();
  real ny = config["ny"].as<real>();
  std::string vcoords_file = config["vcoords"].as<std::string>();
  yakl::SimpleNetCDF nc;
  nc.open(vcoords_file);
  nz = nc.getDimSize("num_interfaces") - 1;
  nc.close();
  real nens = config["nens"].as<real>();
  dm.register_and_allocate<real>( "density_dry"  , "dry density"               , {nz,ny,nx,nens} , {"z","y","x"} );
  dm.register_and_allocate<real>( "uvel"         , "x-direction velocity"      , {nz,ny,nx,nens} , {"z","y","x"} );
  dm.register_and_allocate<real>( "vvel"         , "y-direction velocity"      , {nz,ny,nx,nens} , {"z","y","x"} );
  dm.register_and_allocate<real>( "wvel"         , "z-direction velocity"      , {nz,ny,nx,nens} , {"z","y","x"} );
  dm.register_and_allocate<real>( "theta_dry"    , "dry potential temperature" , {nz,ny,nx,nens} , {"z","y","x"} );
  dm.register_and_allocate<real>( "pressure_dry" , "dry pressure"              , {nz,ny,nx,nens} , {"z","y","x"} );
  yakl::memset( dm.get_collapsed("density_dry" ) , 0._fp );
  yakl::memset( dm.get_collapsed("uvel"        ) , 0._fp );
  yakl::memset( dm.get_collapsed("vvel"        ) , 0._fp );
  yakl::memset( dm.get_collapsed("wvel"        ) , 0._fp );
  yakl::memset( dm.get_collapsed("theta_dry"   ) , 0._fp );
  yakl::memset( dm.get_collapsed("pressure_dry") , 0._fp );
}



