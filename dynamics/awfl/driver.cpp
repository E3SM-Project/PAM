
#include "const.h"
#include "Spatial_cons_expl_fv_Agrid.h"
#include "Temporal_ssprk3.h"
#include "Profiles.h"
#include "Microphysics_kessler.h"
#include "DataManager.h"

// Define the Spatial operator based on constants from the Temporal operatora header file
typedef Spatial_operator<nTimeDerivs,timeAvg,nAder> Spatial;

// Define the Temporal operator based on the Spatial operator
typedef Temporal_operator<Spatial> Dycore;

typedef Microphysics Microphysics;


template <class MICRO>
void allocate_coupler_state( std::string inFile , DataManager &dm , MICRO &micro );


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

    // Create the dycore and the microphysics
    Dycore       dycore;
    Microphysics micro;

    allocate_coupler_state( inFile , dm , micro );

    // Initialize the dycore and the microphysics
    dycore.init( inFile , micro.num_tracers , dm );
    micro .init( inFile , dycore            , dm );

    // Initialize the dry state
    dycore.init_state_and_tracers( dm , micro );

    real etime = 0;

    dycore.output( dm , micro , etime );
    
    while (etime < simTime) {
      real dt = dycore.compute_time_step( 0.4 , dm , micro );
      if (etime + dt > simTime) { dt = simTime - etime; }
      dycore.timeStep( dm , micro , dt );
      dycore.convert_dynamics_to_coupler_state( dm , micro );
      micro.timeStep( dm , dt );
      dycore.convert_coupler_state_to_dynamics( dm , micro );
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



template <class MICRO>
void allocate_coupler_state( std::string inFile , DataManager &dm , MICRO &micro ) {
  YAML::Node config = YAML::LoadFile(inFile);
  if ( !config            ) { endrun("ERROR: Invalid YAML input file"); }
  int nx = config["nx"].as<real>();
  int ny = config["ny"].as<real>();
  std::string vcoords_file = config["vcoords"].as<std::string>();
  yakl::SimpleNetCDF nc;
  nc.open(vcoords_file);
  int nz = nc.getDimSize("num_interfaces") - 1;
  nc.close();
  int nens = config["nens"].as<real>();

  dm.register_and_allocate<real>( "density_dry"  , "dry density"          , {nz,ny,nx,nens} , {"z","y","x","nens"} );
  dm.register_and_allocate<real>( "uvel"         , "x-direction velocity" , {nz,ny,nx,nens} , {"z","y","x","nens"} );
  dm.register_and_allocate<real>( "vvel"         , "y-direction velocity" , {nz,ny,nx,nens} , {"z","y","x","nens"} );
  dm.register_and_allocate<real>( "wvel"         , "z-direction velocity" , {nz,ny,nx,nens} , {"z","y","x","nens"} );
  dm.register_and_allocate<real>( "temp"         , "temperature"          , {nz,ny,nx,nens} , {"z","y","x","nens"} );
  dm.register_and_allocate<real>( "pressure_dry" , "dry pressure"         , {nz,ny,nx,nens} , {"z","y","x","nens"} );

  auto density_dry  = dm.get_collapsed<real>("density_dry" );
  auto uvel         = dm.get_collapsed<real>("uvel"        );
  auto vvel         = dm.get_collapsed<real>("vvel"        );
  auto wvel         = dm.get_collapsed<real>("wvel"        );
  auto temp         = dm.get_collapsed<real>("temp"        );
  auto pressure_dry = dm.get_collapsed<real>("pressure_dry");

  parallel_for( Bounds<1>(nz*ny*nx*nens) , YAKL_LAMBDA (int i) {
    density_dry (i) = 0;
    uvel        (i) = 0;
    vvel        (i) = 0;
    wvel        (i) = 0;
    temp        (i) = 0;
    pressure_dry(i) = 0;
  });

}



