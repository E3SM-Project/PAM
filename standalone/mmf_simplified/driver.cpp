
#include "pam_const.h"
#include "pam_coupler.h"
#include "Dycore.h"
#include "Microphysics.h"
#include "SGS.h"
#include "DataManager.h"
#include "mmf_interface.h"
#include "perturb_temperature.h"
#include "gcm_forcing.h"
#include "gcm_density_forcing.h"
#include "sponge_layer.h"
#include "saturation_adjustment.h"


void output( pam::PamCoupler const &coupler , real etime );


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
    real out_freq       = config["outFreq"    ].as<real>();
    std::string coldata = config["column_data"].as<std::string>();
    bool advect_tke     = config["advect_tke" ].as<bool>();

    int nranks;
    int myrank;
    MPI_Comm_size( MPI_COMM_WORLD , &nranks );
    MPI_Comm_rank( MPI_COMM_WORLD , &myrank );
    bool mainproc = (myrank == 0);

    auto &coupler = mmf_interface::get_coupler();

    coupler.set_option<bool>("advect_tke",advect_tke);

    // How to apply the density forcing: loose (force it like sam does); strict (enforce it strictly every time step)
    if (config["density_forcing"]) {
      coupler.set_option<std::string>("density_forcing",config["density_forcing"].as<std::string>());
    } else {
      coupler.set_option<std::string>("density_forcing","loose");
    }

    // Apply forcing at dynamics time step?  If false, it's applied at the CRM physics time step
    if (config["forcing_at_dycore_time_step"]) {
      coupler.set_option<bool>( "forcing_at_dycore_time_step" , config["forcing_at_dycore_time_step"].as<bool>() );
    } else {
      coupler.set_option<bool>( "forcing_at_dycore_time_step" , false );
    }

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

    coupler.dm.register_and_allocate<real>("gcm_density_dry","GCM column dry density"     ,{crm_nz,nens});
    coupler.dm.register_and_allocate<real>("gcm_uvel"       ,"GCM column u-velocity"      ,{crm_nz,nens});
    coupler.dm.register_and_allocate<real>("gcm_vvel"       ,"GCM column v-velocity"      ,{crm_nz,nens});
    coupler.dm.register_and_allocate<real>("gcm_wvel"       ,"GCM column w-velocity"      ,{crm_nz,nens});
    coupler.dm.register_and_allocate<real>("gcm_temp"       ,"GCM column temperature"     ,{crm_nz,nens});
    coupler.dm.register_and_allocate<real>("gcm_water_vapor","GCM column water vapor mass",{crm_nz,nens});

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

    auto gcm_rho_d = coupler.dm.get<real,2>("gcm_density_dry");
    auto gcm_uvel  = coupler.dm.get<real,2>("gcm_uvel"       );
    auto gcm_vvel  = coupler.dm.get<real,2>("gcm_vvel"       );
    auto gcm_wvel  = coupler.dm.get<real,2>("gcm_wvel"       );
    auto gcm_temp  = coupler.dm.get<real,2>("gcm_temp"       );
    auto gcm_rho_v = coupler.dm.get<real,2>("gcm_water_vapor");

    auto rho_d = coupler.dm.get<real const,4>("density_dry" );
    auto uvel  = coupler.dm.get<real const,4>("uvel"        );
    auto vvel  = coupler.dm.get<real const,4>("vvel"        );
    auto wvel  = coupler.dm.get<real const,4>("wvel"        );
    auto temp  = coupler.dm.get<real const,4>("temp"        );
    auto rho_v = coupler.dm.get<real const,4>("water_vapor" );

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

    coupler.add_mmf_function( "sgs"    , [&] (PamCoupler &coupler, real dt) { sgs   .timeStep(coupler,dt); } );
    coupler.add_mmf_function( "micro"  , [&] (PamCoupler &coupler, real dt) { micro .timeStep(coupler,dt); } );
    coupler.add_mmf_function( "dycore" , [&] (PamCoupler &coupler, real dt) { dycore.timeStep(coupler,dt); } );

    bool forcing_at_dycore_time_step = coupler.get_option<bool>("forcing_at_dycore_time_step");

    if (forcing_at_dycore_time_step) {
      if ( coupler.get_option<std::string>("density_forcing") == "strict" ) {
        coupler.add_dycore_function( "gcm_density_forcing" , modules::gcm_density_forcing );
      }
      coupler.add_dycore_function( "apply_gcm_forcing_tendencies" , modules::apply_gcm_forcing_tendencies );
      coupler.add_dycore_function( "sponge_layer"                 , modules::sponge_layer                 );
    } else {
      if ( coupler.get_option<std::string>("density_forcing") == "strict" ) {
        coupler.add_mmf_function( "gcm_density_forcing" , modules::gcm_density_forcing );
      }
      coupler.add_mmf_function( "apply_gcm_forcing_tendencies" , modules::apply_gcm_forcing_tendencies );
      coupler.add_mmf_function( "sponge_layer"                 , modules::sponge_layer                 );
    }

    // coupler.add_dycore_function( "saturation_adjustment" , saturation_adjustment );

    real etime_gcm = 0;
    int  num_out = 0;

    // Output the initial state
    output( coupler , etime_gcm );

    while (etime_gcm < simTime) {
      if (etime_gcm + dt_gcm > simTime) { dt_gcm = simTime - etime_gcm; }

      modules::compute_gcm_forcing_tendencies( coupler , dt_gcm );

      real etime_crm = 0;
      real simTime_crm = dt_gcm;
      real dt_crm = dt_crm_phys;
      while (etime_crm < simTime_crm) {
        if (dt_crm == 0.) { dt_crm = dycore.compute_time_step(coupler); }
        if (etime_crm + dt_crm > simTime_crm) { dt_crm = simTime_crm - etime_crm; }

        coupler.run_mmf_function( "sgs"    , dt_crm );
        coupler.run_mmf_function( "micro"  , dt_crm );
        coupler.run_mmf_function( "dycore" , dt_crm );
        if (! forcing_at_dycore_time_step) {
          coupler.run_mmf_function( "apply_gcm_forcing_tendencies" , dt_crm );
          coupler.run_mmf_function( "sponge_layer"                 , dt_crm );
        }

        etime_crm += dt_crm;
        etime_gcm += dt_crm;
        if (out_freq >= 0. && etime_gcm / out_freq >= num_out+1) {
          yakl::timer_start("output");
          output( coupler , etime_gcm );
          yakl::timer_stop("output");
          real maxw = maxval(abs(coupler.dm.get_collapsed<real const>("wvel")));
          if (mainproc) {
            std::cout << "Etime , dtphys, maxw: " << etime_gcm << " , " 
                                                  << dt_crm    << " , "
                                                  << std::setw(10) << maxw << "\n";
          }
          num_out++;
        }
      }
    }

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


void output( pam::PamCoupler const &coupler , real etime ) {
  int nranks;
  int myrank;
  MPI_Comm_size( MPI_COMM_WORLD , &nranks );
  MPI_Comm_rank( MPI_COMM_WORLD , &myrank );

  std::string inFile = coupler.get_option<std::string>( "standalone_input_file" );
  YAML::Node config = YAML::LoadFile(inFile);
  auto out_prefix = config["out_prefix"].as<std::string>();

  auto dx = coupler.get_dx();
  auto dy = coupler.get_dy();
  auto nx = coupler.get_nx();
  auto ny = coupler.get_ny();
  auto nz = coupler.get_nz();

  MPI_Barrier(MPI_COMM_WORLD);
  for (int rr=0; rr < nranks; rr++) {
    MPI_Barrier(MPI_COMM_WORLD);
    if (rr == myrank) {



      std::string fname = out_prefix + std::string("_") + std::to_string(myrank) + std::string(".nc");

      yakl::SimpleNetCDF nc;
      int ulIndex = 0; // Unlimited dimension index to place this data at
      // Create or open the file
      if (etime == 0.) {
        nc.create(fname);

        // x-coordinate
        real1d xloc("xloc",nx);
        parallel_for( "Spatial.h output 1" , nx , YAKL_LAMBDA (int i) { xloc(i) = (i+0.5)*dx; });
        nc.write(xloc.createHostCopy(),"x",{"x"});

        // y-coordinate
        real1d yloc("yloc",ny);
        parallel_for( "Spatial.h output 2" , ny , YAKL_LAMBDA (int i) { yloc(i) = (i+0.5)*dy; });
        nc.write(yloc.createHostCopy(),"y",{"y"});

        // z-coordinate
        auto zint = coupler.dm.get<real const,2>("vertical_interface_height");
        real1d zmid("zmid",nz);
        parallel_for( "Spatial.h output 3" , nz , YAKL_LAMBDA (int i) {
          zmid(i) = ( zint(i,0) + zint(i+1,0) ) / 2;
        });
        nc.write(zmid.createHostCopy(),"z",{"z"});

        // Create time variable
        nc.write1(0._fp,"t",0,"t");
      } else {
        nc.open(fname,yakl::NETCDF_MODE_WRITE);
        ulIndex = nc.getDimSize("t");

        // Write the elapsed time
        nc.write1(etime,"t",ulIndex,"t");
      }

      std::vector<std::string> tracer_names = coupler.get_tracer_names();
      int num_tracers = coupler.get_num_tracers();
      // Create MultiField of all state and tracer full variables, since we're doing the same operation on each
      pam::MultiField<real const,4> fields;
      fields.add_field( coupler.dm.get<real const,4>("density_dry") );
      fields.add_field( coupler.dm.get<real const,4>("uvel"       ) );
      fields.add_field( coupler.dm.get<real const,4>("vvel"       ) );
      fields.add_field( coupler.dm.get<real const,4>("wvel"       ) );
      fields.add_field( coupler.dm.get<real const,4>("temp"       ) );
      for (int tr=0; tr < num_tracers; tr++) {
        fields.add_field( coupler.dm.get<real const,4>(tracer_names[tr]) );
      }

      // First, write out standard coupler state
      real3d data("data",nz,ny,nx);
      parallel_for( SimpleBounds<3>(nz,ny,nx) , YAKL_LAMBDA (int k, int j, int i) { data(k,j,i) = fields(0,k,j,i,0); });
      nc.write1(data,"density"    ,{"z","y","x"},ulIndex,"t");
      parallel_for( SimpleBounds<3>(nz,ny,nx) , YAKL_LAMBDA (int k, int j, int i) { data(k,j,i) = fields(1,k,j,i,0); });
      nc.write1(data,"uvel"       ,{"z","y","x"},ulIndex,"t");
      parallel_for( SimpleBounds<3>(nz,ny,nx) , YAKL_LAMBDA (int k, int j, int i) { data(k,j,i) = fields(2,k,j,i,0); });
      nc.write1(data,"vvel"       ,{"z","y","x"},ulIndex,"t");
      parallel_for( SimpleBounds<3>(nz,ny,nx) , YAKL_LAMBDA (int k, int j, int i) { data(k,j,i) = fields(3,k,j,i,0); });
      nc.write1(data,"wvel"       ,{"z","y","x"},ulIndex,"t");
      parallel_for( SimpleBounds<3>(nz,ny,nx) , YAKL_LAMBDA (int k, int j, int i) { data(k,j,i) = fields(4,k,j,i,0); });
      nc.write1(data,"temperature",{"z","y","x"},ulIndex,"t");
      for (int tr=0; tr < num_tracers; tr++) {
        parallel_for( SimpleBounds<3>(nz,ny,nx) , YAKL_LAMBDA (int k, int j, int i) {
          // data(k,j,i) = fields(5+tr,k,j,i,0) / (state(idR,hs+k,hs+j,hs+i,0) + hyDensCells(k,0));
          data(k,j,i) = fields(5+tr,k,j,i,0);
        });
        nc.write1(data,tracer_names[tr],{"z","y","x"},ulIndex,"t");
      }

      // Close the file
      nc.close();



    }
    MPI_Barrier(MPI_COMM_WORLD);
  }
  MPI_Barrier(MPI_COMM_WORLD);
}



