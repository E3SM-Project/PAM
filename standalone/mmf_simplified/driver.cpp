
#include "pam_const.h"
#include "pam_coupler.h"
#include "Dycore.h"
#include "Microphysics.h"
#include "DataManager.h"
#include "mmf_interface.h"
#include "perturb_temperature.h"


int main(int argc, char** argv) {
  yakl::init();
  {
    yakl::timer_start("main");

    // Creates one PamCoupler object in static scope (no parameters assumes 1 thread)
    mmf_interface::init();

    if (argc <= 1) { endrun("ERROR: Must pass the input YAML filename as a parameter"); }
    std::string inFile(argv[1]);
    YAML::Node config = YAML::LoadFile(inFile);
    if ( !config            ) { endrun("ERROR: Invalid YAML input file"); }
    real simTime        = config["simTime"    ].as<real>();
    int  crm_nx         = config["crm_nx"     ].as<int>();
    int  crm_ny         = config["crm_ny"     ].as<int>();
    int  nens           = config["nens"       ].as<int>();
    real xlen           = config["xlen"       ].as<real>();
    real ylen           = config["ylen"       ].as<real>();
    real dt_gcm         = config["dt_gcm"     ].as<real>();
    real dt_crm_phys    = config["dt_crm_phys"].as<real>();
    real outFreq        = config["outFreq"    ].as<real>();
    std::string coldata = config["column_data"].as<std::string>();

    // Store vertical coordinates
    std::string vcoords_file = config["vcoords"].as<std::string>();
    yakl::SimpleNetCDF nc;
    nc.open(vcoords_file);
    int crm_nz = nc.getDimSize("num_interfaces") - 1;
    real1d zint_in;
    nc.read(zint_in,"vertical_interfaces");
    nc.close();

    // Allocates the coupler state (density_dry, uvel, vvel, wvel, temp, vert grid, hydro background) for thread 0
    mmf_interface::allocate_coupler_state( crm_nz , crm_ny , crm_nx , nens );

    // Create the dycore and the microphysics
    Dycore       dycore;
    Microphysics micro;

    // Set physical constants for coupler at thread 0 using microphysics data
    mmf_interface::set_phys_constants( micro.R_d , micro.R_v , micro.cp_d , micro.cp_v , micro.grav , micro.p0 );

    // Set the vertical grid in the coupler
    mmf_interface::set_grid( xlen , ylen , zint_in );

    int numOut = 0;

    // This is for the dycore to pull out to determine how to do idealized test cases
    mmf_interface::set_option( "standalone_input_file" , inFile );

    ///////////////////////////////////////////////////////////////////////////////
    // This is the end of the mmf_interface code. Using coupler directly from here
    ///////////////////////////////////////////////////////////////////////////////
    auto &coupler = mmf_interface::get_coupler();

    // Initialize the microphysics (registers the microphysics's tracers)
    micro .init( coupler );
    // Before calling dycore.init(coupler), be sure ALL TRACERS HAVE BEEN REGISTERED ALREADY in the coupler
    dycore.init( coupler );

    #ifdef PAM_STANDALONE
      std::cout << "Dycore: " << dycore.dycore_name() << std::endl;
      std::cout << "Micro : " << micro .micro_name() << std::endl;
    #endif

    // Initialize CRM data using supercell background state
    nc.open(coldata,yakl::NETCDF_MODE_READ);
    real1d input_zmid, input_rho_d, input_uvel, input_vvel, input_wvel, input_temp, input_rho_v, input_rho_c;
    // These read calls will allocate the arrays internally
    nc.read(input_zmid  ,"z"           );
    nc.read(input_rho_d ,"density_dry" );
    nc.read(input_uvel  ,"uvel"        );
    nc.read(input_vvel  ,"vvel"        );
    nc.read(input_wvel  ,"wvel"        );
    nc.read(input_temp  ,"temp"        );
    nc.read(input_rho_v ,"water_vapor" );
    nc.read(input_rho_c ,"cloud_liquid");
    nc.close();

    int input_nz = input_zmid.dimension[0];

    auto crm_zmid = coupler.dm.get<real,2>("vertical_midpoint_height");

    real1d gcm_rho_d("gcm_rho_d",crm_nz);
    real1d gcm_uvel ("gcm_uvel ",crm_nz);
    real1d gcm_vvel ("gcm_vvel ",crm_nz);
    real1d gcm_wvel ("gcm_wvel ",crm_nz);
    real1d gcm_temp ("gcm_temp ",crm_nz);
    real1d gcm_rho_v("gcm_rho_v",crm_nz);
    real1d gcm_rho_c("gcm_rho_c",crm_nz);

    // Linear interpolation from col data to PAM model data
    parallel_for( crm_nz , YAKL_LAMBDA (int k) {
      int iens = 0;

      int  ind1, ind2;
      real d1  , d2  ;
      if (crm_zmid(k,iens) < input_zmid(0)) {
        ind1 = 0;
        ind2 = 0;
        d1   = 1;
        d2   = 1;
      } else if (crm_zmid(k,iens) > input_zmid(input_nz-1)) {
        ind1 = input_nz-1;
        ind2 = input_nz-1;
        d1   = 1;
        d2   = 1;
      } else {
        for (int kk=0; kk < input_nz; kk++) {
          if (crm_zmid(k,iens) > input_zmid(kk)) {
            ind1 = kk;
          } else {
            break;
          }
        }
        ind2 = ind1+1;
        d1 = crm_zmid(k,iens) - input_zmid(ind1);
        d2 = input_zmid(ind2) - crm_zmid(k,iens);
      }

      real w1 = 1._fp - d1 / (d1 + d2);
      real w2 = 1._fp - d2 / (d1 + d2);
      
      gcm_rho_d(k) = w1 * input_rho_d(ind1) + w2 * input_rho_d(ind2);
      gcm_uvel (k) = w1 * input_uvel (ind1) + w2 * input_uvel (ind2);
      gcm_vvel (k) = w1 * input_vvel (ind1) + w2 * input_vvel (ind2);
      gcm_wvel (k) = w1 * input_wvel (ind1) + w2 * input_wvel (ind2);
      gcm_temp (k) = w1 * input_temp (ind1) + w2 * input_temp (ind2);
      gcm_rho_v(k) = w1 * input_rho_v(ind1) + w2 * input_rho_v(ind2);
      gcm_rho_c(k) = w1 * input_rho_c(ind1) + w2 * input_rho_c(ind2);
    });

    auto rho_d = coupler.dm.get<real,4>("density_dry" );
    auto uvel  = coupler.dm.get<real,4>("uvel"        );
    auto vvel  = coupler.dm.get<real,4>("vvel"        );
    auto wvel  = coupler.dm.get<real,4>("wvel"        );
    auto temp  = coupler.dm.get<real,4>("temp"        );
    auto rho_v = coupler.dm.get<real,4>("water_vapor" );
    auto rho_c = coupler.dm.get<real,4>("cloud_liquid");

    // Set CRM data to uniform column
    parallel_for( Bounds<4>(crm_nz,crm_ny,crm_nx,nens) , YAKL_LAMBDA (int k, int j, int i, int iens) {
      rho_d(k,j,i,iens) = gcm_rho_d(k);
      uvel (k,j,i,iens) = gcm_uvel (k);
      vvel (k,j,i,iens) = gcm_vvel (k);
      wvel (k,j,i,iens) = gcm_wvel (k);
      temp (k,j,i,iens) = gcm_temp (k);
      rho_v(k,j,i,iens) = gcm_rho_v(k);
      rho_c(k,j,i,iens) = gcm_rho_c(k);
    });

    perturb_temperature( coupler , 0 );

    // Now that we have an initial state, define hydrostasis for each ensemble member
    coupler.update_hydrostasis( coupler.compute_pressure_array() );

    real etime = 0;

    // if (outFreq >= 0) dycore.output( coupler , etime );

    // real dtphys = dtphys_in;
    // while (etime < simTime) {
    //   if (dtphys_in == 0.) { dtphys = dycore.compute_time_step(coupler); }
    //   if (etime + dtphys > simTime) { dtphys = simTime - etime; }

    //   yakl::timer_start("micro");
    //   micro.timeStep( coupler , dtphys );
    //   yakl::timer_stop("micro");

    //   yakl::timer_start("dycore");
    //   dycore.timeStep( coupler , dtphys );
    //   yakl::timer_stop("dycore");

    //   etime += dtphys;
    //   if (outFreq >= 0. && etime / outFreq >= numOut+1) {
    //     std::cout << "Etime , dtphys: " << etime << " , " << dtphys << "\n";
    //     yakl::timer_start("output");
    //     dycore.output( coupler , etime );
    //     yakl::timer_stop("output");
    //     numOut++;
    //   }
    // }

    // std::cout << "Elapsed Time: " << etime << "\n";

    // dycore.finalize( coupler );

    mmf_interface::finalize();

    yakl::timer_stop("main");
  }
  yakl::finalize();
}
