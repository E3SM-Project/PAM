
#include "pam_const.h"
#include "pam_coupler.h"
#include "Dycore.h"
#include "Microphysics.h"
#include "DataManager.h"


int main(int argc, char** argv) {
  yakl::init();
  {
    yakl::timer_start("main");

    DataManager dm;

    if (argc <= 1) { endrun("ERROR: Must pass the input YAML filename as a parameter"); }
    std::string inFile(argv[1]);
    YAML::Node config = YAML::LoadFile(inFile);
    if ( !config            ) { endrun("ERROR: Invalid YAML input file"); }
    real simTime   = config["simTime"].as<real>();
    real outFreq   = config["outFreq"].as<real>();
    int  nx        = config["nx"     ].as<int>();
    int  ny        = config["ny"     ].as<int>();
    int  nens      = config["nens"   ].as<int>();
    real xlen      = config["xlen"   ].as<real>();
    real ylen      = config["ylen"   ].as<real>();
    real dtphys_in = config["dtphys" ].as<real>();

    // Store vertical coordinates
    std::string vcoords_file = config["vcoords"].as<std::string>();
    yakl::SimpleNetCDF nc;
    nc.open(vcoords_file);
    int nz = nc.getDimSize("num_interfaces") - 1;
    real1d zint_in("zint_in",nz+1);
    nc.read(zint_in,"vertical_interfaces");
    nc.close();

    dm.register_and_allocate<real>( "vertical_interface_height" , "vertical_interface_height" , {nz+1,nens} , {"zp1","nens"} );
    auto zint = dm.get<real,2>("vertical_interface_height");
    parallel_for( "driver.cpp 1" , Bounds<2>(nz+1,nens) , YAKL_LAMBDA (int k, int iens) {
      zint(k,iens) = zint_in(k);
    });

    dm.register_and_allocate<real>( "vertical_midpoint_height" , "vertical_midpoint_heignt" , {nz,nens} , {"z","nens"} );
    auto zmid = dm.get<real,2>("vertical_midpoint_height");
    parallel_for("driver.cpp 2"  , Bounds<2>(nz,nens) , YAKL_LAMBDA (int k, int iens) {
      zmid(k,iens) = 0.5_fp*(zint_in(k) + zint_in(k+1));
    });

    int numOut = 0;

    // Create the dycore and the microphysics
    Dycore       dycore;
    Microphysics micro;

    allocate_coupler_state( nz , ny , nx , nens , dm );

    // Initialize the dycore and the microphysics
    dycore.init( inFile , ny , nx , nens , xlen , ylen , micro.get_num_tracers() , dm );
    micro .init( inFile , ny , nx , nens , dycore , dm );

    #ifdef PAM_STANDALONE
      std::cout << "Dycore: " << dycore.dycore_name() << std::endl;
      std::cout << "Micro : " << micro .micro_name() << std::endl;
    #endif

    // Initialize the dry state
    dycore.init_state_and_tracers( dm , micro );

    real etime = 0;

    if (outFreq >= 0) dycore.output( dm , micro , etime );

    real dtphys = dtphys_in;
    while (etime < simTime) {
      if (dtphys_in == 0.) { dtphys = dycore.compute_time_step(dm, micro); }
      if (etime + dtphys > simTime) { dtphys = simTime - etime; }

      yakl::timer_start("micro");
      micro.timeStep( dm , dtphys );
      yakl::timer_stop("micro");

      yakl::timer_start("dycore");
      dycore.timeStep( dm , micro , dtphys );
      yakl::timer_stop("dycore");

      etime += dtphys;
      if (outFreq >= 0. && etime / outFreq >= numOut+1) {
        std::cout << "Etime , dtphys: " << etime << " , " << dtphys << "\n";
        yakl::timer_start("output");
        dycore.output( dm , micro , etime );
        yakl::timer_stop("output");
        numOut++;
      }
    }

    std::cout << "Elapsed Time: " << etime << "\n";

    dycore.finalize( dm );

    yakl::timer_stop("main");
  }
  yakl::finalize();
}
