
#include "pam_coupler.h"
#include "Dycore.h"
#include "Microphysics.h"
#include "SGS.h"
#include "mmf_interface.h"
#include "perturb_temperature.h"
#include "gcm_forcing.h"
#include "sponge_layer.h"
#include "saturation_adjustment.h"
#include "output.h"


int main(int argc, char** argv) {
  MPI_Init( &argc , &argv );
  yakl::init();
  {
    using yakl::intrinsics::abs;
    using yakl::intrinsics::maxval;
    yakl::timer_start("main");

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
    real out_freq       = config["out_freq"   ].as<real>();

    int nranks;
    int myrank;
    MPI_Comm_size( MPI_COMM_WORLD , &nranks );
    MPI_Comm_rank( MPI_COMM_WORLD , &myrank );
    bool mainproc = (myrank == 0);

    auto &coupler = mmf_interface::get_coupler();

    // This is for the dycore to pull out to determine how to do idealized test cases
    coupler.set_option<std::string>( "standalone_input_file" , inFile );

    // Store vertical coordinates
    std::string vcoords_file = config["vcoords"].as<std::string>();
    yakl::SimpleNetCDF nc;
    nc.open(vcoords_file);
    int crm_nz = nc.getDimSize("num_interfaces") - 1;
    real1d zint_in;
    nc.read(zint_in,"vertical_interfaces");
    nc.close();

    // Allocates the coupler state (density_dry, uvel, vvel, wvel, temp, vert grid, hydro background) for thread 0
    coupler.allocate_coupler_state( crm_nz , crm_ny , crm_nx , nens );

    // NORMALLY THIS WOULD BE DONE INSIDE THE CRM, BUT WE'RE USING CONSTANTS DEFINED BY THE CRM MICRO SCHEME
    // Create the dycore and the microphysics
    Dycore       dycore;
    Microphysics micro;
    SGS          sgs;

    // Set physical constants for coupler at thread 0 using microphysics data
    coupler.set_phys_constants( micro.R_d , micro.R_v , micro.cp_d , micro.cp_v , micro.grav , micro.p0 );

    // Set the vertical grid in the coupler
    coupler.set_grid( xlen , ylen , zint_in );

    auto &dm = coupler.get_data_manager_readwrite();

    dm.register_and_allocate<real>("gcm_density_dry","GCM column dry density"     ,{crm_nz,nens});
    dm.register_and_allocate<real>("gcm_uvel"       ,"GCM column u-velocity"      ,{crm_nz,nens});
    dm.register_and_allocate<real>("gcm_vvel"       ,"GCM column v-velocity"      ,{crm_nz,nens});
    dm.register_and_allocate<real>("gcm_wvel"       ,"GCM column w-velocity"      ,{crm_nz,nens});
    dm.register_and_allocate<real>("gcm_temp"       ,"GCM column temperature"     ,{crm_nz,nens});
    dm.register_and_allocate<real>("gcm_water_vapor","GCM column water vapor mass",{crm_nz,nens});

    micro .init( coupler );
    sgs   .init( coupler );
    dycore.init( coupler );  // dycore should set idealized conditions here

    #ifdef PAM_STANDALONE
      if (mainproc) {
        std::cout << "Dycore: " << dycore.dycore_name() << std::endl;
        std::cout << "Micro : " << micro .micro_name () << std::endl;
        std::cout << "SGS   : " << sgs   .sgs_name   () << std::endl;
        std::cout << "\n";
        std::cout << "crm_nx:   " << crm_nx << "\n";
        std::cout << "crm_ny:   " << crm_ny << "\n";
        std::cout << "crm_nz:   " << crm_nz << "\n";
        std::cout << "xlen (m): " << xlen << "\n";
        std::cout << "ylen (m): " << ylen << "\n";
        std::cout << "Vertical interface heights: ";
        auto zint_host = zint_in.createHostCopy();
        for (int k=0; k < crm_nz+1; k++) {
          std::cout << zint_host(k) << "  ";
        }
        std::cout << "\n\n";
      }
    #endif

    auto gcm_rho_d = dm.get<real,2>("gcm_density_dry");
    auto gcm_uvel  = dm.get<real,2>("gcm_uvel"       );
    auto gcm_vvel  = dm.get<real,2>("gcm_vvel"       );
    auto gcm_wvel  = dm.get<real,2>("gcm_wvel"       );
    auto gcm_temp  = dm.get<real,2>("gcm_temp"       );
    auto gcm_rho_v = dm.get<real,2>("gcm_water_vapor");

    auto rho_d = dm.get<real const,4>("density_dry" );
    auto uvel  = dm.get<real const,4>("uvel"        );
    auto vvel  = dm.get<real const,4>("vvel"        );
    auto wvel  = dm.get<real const,4>("wvel"        );
    auto temp  = dm.get<real const,4>("temp"        );
    auto rho_v = dm.get<real const,4>("water_vapor" );

    // Compute a column to force the model with by averaging the columns at init
    parallel_for( Bounds<4>(crm_nz,crm_ny,crm_nx,nens) , YAKL_LAMBDA (int k, int j, int i, int iens) {
      gcm_rho_d(k,iens) = 0;
      gcm_uvel (k,iens) = 0;
      gcm_vvel (k,iens) = 0;
      gcm_wvel (k,iens) = 0;
      gcm_temp (k,iens) = 0;
      gcm_rho_v(k,iens) = 0;
    });
    real r_nx_ny = 1._fp / (crm_nx*crm_ny);  // Avoid costly divisions
    parallel_for( Bounds<4>(crm_nz,crm_ny,crm_nx,nens) , YAKL_LAMBDA (int k, int j, int i, int iens) {
      yakl::atomicAdd( gcm_rho_d(k,iens) , rho_d(k,j,i,iens) * r_nx_ny );
      yakl::atomicAdd( gcm_uvel (k,iens) , uvel (k,j,i,iens) * r_nx_ny );
      yakl::atomicAdd( gcm_vvel (k,iens) , vvel (k,j,i,iens) * r_nx_ny );
      yakl::atomicAdd( gcm_wvel (k,iens) , wvel (k,j,i,iens) * r_nx_ny );
      yakl::atomicAdd( gcm_temp (k,iens) , temp (k,j,i,iens) * r_nx_ny );
      yakl::atomicAdd( gcm_rho_v(k,iens) , rho_v(k,j,i,iens) * r_nx_ny );
    });

    modules::perturb_temperature( coupler , 0 );

    // Now that we have an initial state, define hydrostasis for each ensemble member
    coupler.update_hydrostasis( coupler.compute_pressure_array() );

    coupler.add_mmf_function( "apply_gcm_forcing_tendencies" , modules::apply_gcm_forcing_tendencies );
    coupler.add_mmf_function( "dycore" , [&] (pam::PamCoupler &coupler, real dt) { dycore.timeStep(coupler,dt); } );
    coupler.add_mmf_function( "sgs"    , [&] (pam::PamCoupler &coupler, real dt) { sgs   .timeStep(coupler,dt); } );
    coupler.add_mmf_function( "micro"  , [&] (pam::PamCoupler &coupler, real dt) { micro .timeStep(coupler,dt); } );
    coupler.add_mmf_function( "sponge_layer"                 , modules::sponge_layer                 );
    // coupler.add_dycore_function( "saturation_adjustment" , saturation_adjustment );

    real etime_gcm = 0;
    int  num_out = 0;

    // Output the initial state
    if (out_freq >= 0. ) output( coupler , etime_gcm );

    yakl::timer_start("main_loop");
    while (etime_gcm < simTime) {
      if (etime_gcm + dt_gcm > simTime) { dt_gcm = simTime - etime_gcm; }

      modules::compute_gcm_forcing_tendencies( coupler , dt_gcm );

      real etime_crm = 0;
      real simTime_crm = dt_gcm;
      real dt_crm = dt_crm_phys;
      while (etime_crm < simTime_crm) {
        if (dt_crm == 0.) { dt_crm = dycore.compute_time_step(coupler); }
        if (etime_crm + dt_crm > simTime_crm) { dt_crm = simTime_crm - etime_crm; }

        coupler.run_mmf_function( "apply_gcm_forcing_tendencies" , dt_crm );
        coupler.run_mmf_function( "dycore"                       , dt_crm );
        coupler.run_mmf_function( "sponge_layer"                 , dt_crm );
        coupler.run_mmf_function( "sgs"                          , dt_crm );
        coupler.run_mmf_function( "micro"                        , dt_crm );

        etime_crm += dt_crm;
        etime_gcm += dt_crm;
        if (out_freq >= 0. && etime_gcm / out_freq >= num_out+1) {
          yakl::timer_start("output");
          output( coupler , etime_gcm );
          yakl::timer_stop("output");
          real maxw = maxval(abs(dm.get_collapsed<real const>("wvel")));
          if (mainproc) {
            std::cout << "Etime , dtphys, maxw: " << etime_gcm << " , " 
                                                  << dt_crm    << " , "
                                                  << std::setw(10) << maxw << std::endl;
          }
          num_out++;
        }
      }
    }
    yakl::timer_stop("main_loop");

    if (mainproc) {
      std::cout << "Elapsed Time: " << etime_gcm << "\n";
    }

    dycore.finalize( coupler );

    mmf_interface::finalize();

    yakl::timer_stop("main");
  }
  yakl::finalize();
  MPI_Finalize();
}







