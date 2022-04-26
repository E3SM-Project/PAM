
#include "pam_coupler.h"
#include "Dycore.h"
#include "Microphysics.h"
#include "SGS.h"
#include "mmf_interface.h"
#include "perturb_temperature.h"
#include "gcm_forcing.h"
#include "sponge_layer.h"
#include "saturation_adjustment.h"
#include "broadcast_initial_gcm_column.h"
#include "output.h"
#include "supercell_init.h"


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
    auto simTime        = config["simTime"    ].as<real>();
    auto crm_nx         = config["crm_nx"     ].as<int>();
    auto crm_ny         = config["crm_ny"     ].as<int>();
    auto nens           = config["nens"       ].as<int>();
    auto xlen           = config["xlen"       ].as<real>();
    auto ylen           = config["ylen"       ].as<real>();
    auto dt_gcm         = config["dt_gcm"     ].as<real>();
    auto dt_crm_phys    = config["dt_crm_phys"].as<real>();
    auto out_freq       = config["out_freq"   ].as<real>();
    auto out_prefix     = config["out_prefix" ].as<std::string>();

    int nranks;
    int myrank;
    MPI_Comm_size( MPI_COMM_WORLD , &nranks );
    MPI_Comm_rank( MPI_COMM_WORLD , &myrank );
    bool mainproc = (myrank == 0);

    auto &coupler = mmf_interface::get_coupler();

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

    // Compute a supercell initial column
    real1d rho_d_col("rho_d_col",crm_nz);
    real1d uvel_col ("uvel_col" ,crm_nz);
    real1d vvel_col ("vvel_col" ,crm_nz);
    real1d wvel_col ("wvel_col" ,crm_nz);
    real1d temp_col ("temp_col" ,crm_nz);
    real1d rho_v_col("rho_v_col",crm_nz);

    supercell_init( zint_in, rho_d_col, uvel_col, vvel_col, wvel_col, temp_col, rho_v_col, coupler.get_R_d() ,
                    coupler.get_R_v(), coupler.get_grav() );

    // Set the GCM column data for each ensemble to the supercell initial state
    auto &dm = coupler.get_data_manager_readwrite();

    auto gcm_rho_d = dm.get<real,2>("gcm_density_dry");
    auto gcm_uvel  = dm.get<real,2>("gcm_uvel"       );
    auto gcm_vvel  = dm.get<real,2>("gcm_vvel"       );
    auto gcm_wvel  = dm.get<real,2>("gcm_wvel"       );
    auto gcm_temp  = dm.get<real,2>("gcm_temp"       );
    auto gcm_rho_v = dm.get<real,2>("gcm_water_vapor");

    parallel_for( Bounds<2>(crm_nz,nens) , YAKL_LAMBDA (int k, int iens) {
      gcm_rho_d(k,iens) = rho_d_col(k);
      gcm_uvel (k,iens) = uvel_col (k);
      gcm_vvel (k,iens) = vvel_col (k);
      gcm_wvel (k,iens) = wvel_col (k);
      gcm_temp (k,iens) = temp_col (k);
      gcm_rho_v(k,iens) = rho_v_col(k);
    });

    // NORMALLY THIS WOULD BE DONE INSIDE THE CRM, BUT WE'RE USING CONSTANTS DEFINED BY THE CRM MICRO SCHEME
    // Create the dycore and the microphysics
    Dycore       dycore;
    Microphysics micro;
    SGS          sgs;

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

    // Set physical constants for coupler at thread 0 using microphysics data
    coupler.set_phys_constants( micro.R_d , micro.R_v , micro.cp_d , micro.cp_v , micro.grav , micro.p0 );

    // Set the vertical grid in the coupler
    coupler.set_grid( xlen , ylen , zint_in );

    micro .init( coupler );
    sgs   .init( coupler );
    dycore.init( coupler );

    // Initialize the CRM internal state from the initial GCM column and random temperature perturbations
    modules::broadcast_initial_gcm_column( coupler );

    // Now that we have an initial state, define hydrostasis for each ensemble member
    coupler.update_hydrostasis();

    modules::perturb_temperature( coupler , 0 );

    coupler.add_mmf_function( "apply_gcm_forcing_tendencies" , modules::apply_gcm_forcing_tendencies );
    coupler.add_mmf_function( "dycore" , [&] (pam::PamCoupler &coupler, real dt) { dycore.timeStep(coupler,dt); } );
    coupler.add_mmf_function( "sgs"    , [&] (pam::PamCoupler &coupler, real dt) { sgs   .timeStep(coupler,dt); } );
    coupler.add_mmf_function( "micro"  , [&] (pam::PamCoupler &coupler, real dt) { micro .timeStep(coupler,dt); } );
    coupler.add_mmf_function( "sponge_layer"                 , modules::sponge_layer                 );
    // coupler.add_dycore_function( "saturation_adjustment" , saturation_adjustment );

    real etime_gcm = 0;
    int  num_out = 0;

    // Output the initial state
    if (out_freq >= 0. ) output( coupler , out_prefix , etime_gcm );

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
          output( coupler , out_prefix , etime_gcm );
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







