
#pragma once

#include "awfl_const.h"
#include "DataManager.h"

extern "C" void wsm6( double *theta , double *qv , double *qc , double *qr , double *qi , double *qs ,
                      double *qg , double *rho_dry , double *exner , double *pressure , double *dz ,
                      double &dt , double &grav , double &cp_d , double &cp_v , double &R_d , double &R_v ,
                      double &svpt0 , double &ep_1 , double &ep_2 , double &qmin , double &xls , double &xlv ,
                      double &xlf , double &rhoair0 , double &rhowater , double &cliq , double &cice ,
                      double &psat , double *rainnc , double *rainncv , double *snownc , double *snowncv ,
                      double *sr , double *graupelnc , double *graupelncv ,
                      int &ids , int &ide , int &jds , int &jde , int &kds , int &kde ,
                      int &ims , int &ime , int &jms , int &jme , int &kms , int &kme ,
                      int &its , int &ite , int &jts , int &jte , int &kts , int &kte );


extern "C" void wsm6init(double &rhoair0 , double &rhowater , double &rhosnow , double &cliq , double &cp_v,
                         int &hail );


class Microphysics {
public:
  int static constexpr num_tracers = 6;

  struct Constants {
    real R_d               ;
    real cp_d              ;
    real cv_d              ;
    real gamma_d           ;
    real kappa_d           ;
    real R_v               ;
    real cp_v              ;
    real cv_v              ;
    real p0                ;
  };

  int tracer_index_vapor;

  Constants constants;

  real grav ;
  real svpt0;      // constant for saturation vapor pressure calculation (K)
  real ep_1;       // constant for virtual temperature (r_v/r_d - 1) (dimensionless)
  real ep_2;       // constant for specific humidity calculation (dimensionless)
  real qmin;
  real xls ;       // latent heat of sublimation of water at 0^oC (J kg^-1) 
  real xlv ;       // latent heat of vaporization of water at 0^oC (J kg^-1)
  real xlf ;       // latent heat of fusion of water at 0^oC (J kg^-1)
  real rhoair0  ;  // density of dry air at 0^oC and 1000mb pressure (kg m^-3)
  real rhowater ;  // density of liquid water at 0^oC (kg m^-3)
  real rhosnow  ;  // density of snow (kg m^-3)
  real cliq ;      // specific heat of liquid water at 0^oC
  real cice ;      // specific heat of ice at 0^oC
  real psat ;

  // TODO: Change this to type int instead of real
  SArray<int,1,num_tracers> tracer_IDs; // tracer index for microphysics tracers

  int static constexpr ID_V = 0;  // Local index for water vapor
  int static constexpr ID_C = 1;  // Local index for cloud water
  int static constexpr ID_R = 2;  // Local index for rain water
  int static constexpr ID_I = 3;  // Local index for cloud ice
  int static constexpr ID_S = 4;  // Local index for snow
  int static constexpr ID_G = 5;  // Local index for graupel



  Microphysics() {
    constants.R_d         = 287.;
    constants.cp_d        = 7.*constants.R_d/2.;
    constants.cv_d        = constants.cp_d - constants.R_d;
    constants.gamma_d     = constants.cp_d / constants.cv_d;
    constants.kappa_d     = constants.R_d  / constants.cp_d;
    constants.R_v         = 461.6;
    constants.cp_v        = 4.*constants.R_v;
    constants.cv_v        = constants.R_v - constants.cp_v;
    constants.p0          = 1.e5;
    grav = 9.81;
    svpt0=273.15;      
    ep_1=constants.R_v/constants.R_d-1.; 
    ep_2=constants.R_d/constants.R_v;   
    qmin=1.e-15;
    xls = 2.85E6;      
    xlv = 2.5E6;       
    xlf = 3.50E5;      
    rhoair0  = 1.28 ;  
    rhowater = 1000.;  
    rhosnow  = 100.;
    cliq = 4190.;      
    cice = 2106.;      
    psat = 610.78;
  }



  int get_num_tracers() const {
    return num_tracers;
  }



  int get_water_vapor_index() const {
    return tracer_index_vapor;
  }



  template <class DC>
  void init(std::string infile , int ny, int nx, int nens , DC &dycore , DataManager &dm) {
    int nz = dm.get_dimension_size("z");
    // Register tracers in the dycore
    //                                        name           description      positive   adds mass
    tracer_IDs(ID_V) = dycore.add_tracer(dm , "water_vapor" , "Water Vapor" , true     , true);
    tracer_IDs(ID_C) = dycore.add_tracer(dm , "cloud_water" , "Cloud Water" , true     , true);
    tracer_IDs(ID_R) = dycore.add_tracer(dm , "rain_water"  , "Rain Water"  , true     , true);
    tracer_IDs(ID_I) = dycore.add_tracer(dm , "cloud_ice"   , "Cloud Ice"   , true     , true);
    tracer_IDs(ID_S) = dycore.add_tracer(dm , "snow"        , "Snow"        , true     , true);
    tracer_IDs(ID_G) = dycore.add_tracer(dm , "graupel"     , "Graupel"     , true     , true);

    // Register and allocate the tracers in the DataManager
    dm.register_and_allocate<real>( "water_vapor"   , "Water Vapor" , {nz,ny,nx,nens} , {"z","y","x","nens"} );
    dm.register_and_allocate<real>( "cloud_water"   , "Cloud Water" , {nz,ny,nx,nens} , {"z","y","x","nens"} );
    dm.register_and_allocate<real>( "rain_water"    , "Rain Water"  , {nz,ny,nx,nens} , {"z","y","x","nens"} );
    dm.register_and_allocate<real>( "cloud_ice"     , "Cloud Ice"   , {nz,ny,nx,nens} , {"z","y","x","nens"} );
    dm.register_and_allocate<real>( "snow"          , "Snow"        , {nz,ny,nx,nens} , {"z","y","x","nens"} );
    dm.register_and_allocate<real>( "graupel"       , "Graupel"     , {nz,ny,nx,nens} , {"z","y","x","nens"} );

    tracer_index_vapor = tracer_IDs(ID_V);

    // Register and allocation non-tracer quantities used by the microphysics
    dm.register_and_allocate<real>( "rainnc"     , "Grid scale precipitation (mm)"                    , {ny,nx,nens} , {"y","x","nens"} );
    dm.register_and_allocate<real>( "rainncv"    , "One time step grid scale precipitation (mm/step)" , {ny,nx,nens} , {"y","x","nens"} );
    dm.register_and_allocate<real>( "snownc"     , "Grid scale snow and ice (mm)"                     , {ny,nx,nens} , {"y","x","nens"} );
    dm.register_and_allocate<real>( "snowncv"    , "One time step grid scale snow and ice (mm/step)"  , {ny,nx,nens} , {"y","x","nens"} );
    dm.register_and_allocate<real>( "sr"         , "One time step mass ratio of snow to total precip" , {ny,nx,nens} , {"y","x","nens"} );
    dm.register_and_allocate<real>( "graupelnc"  , "Grid scale graupel (mm)"                          , {ny,nx,nens} , {"y","x","nens"} );
    dm.register_and_allocate<real>( "graupelncv" , "One time step grid scale graupel (mm/step)"       , {ny,nx,nens} , {"y","x","nens"} );

    auto qv = dm.get<real,4>("water_vapor");
    auto qc = dm.get<real,4>("cloud_water");
    auto qr = dm.get<real,4>("rain_water");
    auto qi = dm.get<real,4>("cloud_ice");
    auto qs = dm.get<real,4>("snow");
    auto qg = dm.get<real,4>("graupel");
    auto rainnc     = dm.get<real,3>("rainnc"    );
    auto rainncv    = dm.get<real,3>("rainncv"   );
    auto snownc     = dm.get<real,3>("snownc"    );
    auto snowncv    = dm.get<real,3>("snowncv"   );
    auto sr         = dm.get<real,3>("sr"        );
    auto graupelnc  = dm.get<real,3>("graupelnc" );
    auto graupelncv = dm.get<real,3>("graupelncv");
    parallel_for( SimpleBounds<4>(nz,ny,nx,nens) , YAKL_LAMBDA (int k, int j, int i, int iens) {
      qv(k,j,i,iens) = 0;
      qc(k,j,i,iens) = 0;
      qr(k,j,i,iens) = 0;
      qi(k,j,i,iens) = 0;
      qs(k,j,i,iens) = 0;
      qg(k,j,i,iens) = 0;
      if (k == 0) {
        rainnc    (j,i,iens) = 0;
        rainncv   (j,i,iens) = 0;
        snownc    (j,i,iens) = 0;
        snowncv   (j,i,iens) = 0;
        sr        (j,i,iens) = 0;
        graupelnc (j,i,iens) = 0;
        graupelncv(j,i,iens) = 0;        
      }
    });

    int hail = 0;
    wsm6init(rhoair0,rhowater,rhosnow,cliq,constants.cp_v, hail );
  }



  real compute_total_mass( DataManager &dm ) {
    auto rho_v = dm.get<real,4>("water_vapor");
    auto rho_c = dm.get<real,4>("cloud_water");
    auto rho_r = dm.get<real,4>("rain_water");
    auto rho_i = dm.get<real,4>("cloud_ice");
    auto rho_s = dm.get<real,4>("snow");
    auto rho_g = dm.get<real,4>("graupel");
    auto zint  = dm.get<real,2>("vertical_interface_height");
    int nz   = dm.get_dimension_size("z");
    int ny   = dm.get_dimension_size("y");
    int nx   = dm.get_dimension_size("x");
    int nens = dm.get_dimension_size("nens");
    real4d tmp("tmp",nz,ny,nx,nens);
    parallel_for( SimpleBounds<4>(nz,ny,nx,nens) , YAKL_LAMBDA (int k, int j, int i, int iens) {
      real dz = (zint(k+1,iens) - zint(k,iens));
      tmp(k,j,i,iens) = ( rho_v(k,j,i,iens) + rho_c(k,j,i,iens) + rho_r(k,j,i,iens) +
                          rho_i(k,j,i,iens) + rho_s(k,j,i,iens) + rho_g(k,j,i,iens) ) * dz;
    });
    return yakl::intrinsics::sum(tmp);
  }



  void timeStep( DataManager &dm , real dt ) {
    auto rho_v        = dm.get_lev_col<real>("water_vapor");
    auto rho_c        = dm.get_lev_col<real>("cloud_water");
    auto rho_r        = dm.get_lev_col<real>("rain_water");
    auto rho_i        = dm.get_lev_col<real>("cloud_ice");
    auto rho_s        = dm.get_lev_col<real>("snow");
    auto rho_g        = dm.get_lev_col<real>("graupel");

    auto rho_dry      = dm.get_lev_col<real>("density_dry");
    auto temp         = dm.get_lev_col<real>("temp");
    auto pressure_dry = dm.get_lev_col<real>("pressure_dry");

    int nz   = dm.get_dimension_size("z"   );
    int ny   = dm.get_dimension_size("y"   );
    int nx   = dm.get_dimension_size("x"   );
    int nens = dm.get_dimension_size("nens");
    int ncol = ny*nx*nens;

    #ifdef PAM_DEBUG
      validate_array_positive(rho_v);
      validate_array_positive(rho_c);
      validate_array_positive(rho_r);
      validate_array_positive(rho_i);
      validate_array_positive(rho_s);
      validate_array_positive(rho_g);
      real mass_init = compute_total_mass( dm );
    #endif

    // These are inputs to kessler(...)
    real2d qv      ("qv"      ,nz,ncol);
    real2d qc      ("qc"      ,nz,ncol);
    real2d qr      ("qr"      ,nz,ncol);
    real2d qi      ("qi"      ,nz,ncol);
    real2d qs      ("qs"      ,nz,ncol);
    real2d qg      ("qg"      ,nz,ncol);
    real2d pressure("pressure",nz,ncol);
    real2d theta   ("theta"   ,nz,ncol);
    real2d exner   ("exner"   ,nz,ncol);
    auto zint_in = dm.get<real,2>("vertical_interface_height");

    // We have to broadcast the midpoint heights to all columns within a CRM to avoid the microphysics needing
    // to know about the difference between nx,ny and nens
    real2d dz("dz",nz,ny*nx*nens);
    parallel_for( Bounds<4>(nz,ny,nx,nens) , YAKL_LAMBDA (int k, int j, int i, int iens) {
      dz(k,j*nx*nens + i*nens + iens) = zint_in(k+1,iens) - zint_in(k,iens);
    });

    // Force constants into local scope
    real gamma_d = this->constants.gamma_d;
    real R_d     = this->constants.R_d;
    real R_v     = this->constants.R_v;
    real cp_d    = this->constants.cp_d;
    real p0      = this->constants.p0;

    // Save initial state, and compute inputs for kessler(...)
    parallel_for( Bounds<2>(nz,ncol) , YAKL_LAMBDA (int k, int i) {
      qv      (k,i) = rho_v(k,i) / rho_dry(k,i);
      qc      (k,i) = rho_c(k,i) / rho_dry(k,i);
      qr      (k,i) = rho_r(k,i) / rho_dry(k,i);
      qi      (k,i) = rho_i(k,i) / rho_dry(k,i);
      qs      (k,i) = rho_s(k,i) / rho_dry(k,i);
      qg      (k,i) = rho_g(k,i) / rho_dry(k,i);
      pressure(k,i) = pressure_dry(k,i) + R_v * rho_v(k,i) * temp(k,i);
      exner   (k,i) = pow( pressure(k,i) / p0 , R_d / cp_d );
      theta   (k,i) = temp(k,i) / exner(k,i);
    });

    auto rainnc     = dm.get_collapsed<real>("rainnc"    );
    auto rainncv    = dm.get_collapsed<real>("rainncv"   );
    auto snownc     = dm.get_collapsed<real>("snownc"    );
    auto snowncv    = dm.get_collapsed<real>("snowncv"   );
    auto sr         = dm.get_collapsed<real>("sr"        );
    auto graupelnc  = dm.get_collapsed<real>("graupelnc" );
    auto graupelncv = dm.get_collapsed<real>("graupelncv");

    auto theta_host        = theta       .createHostCopy();   
    auto qv_host           = qv          .createHostCopy();          
    auto qc_host           = qc          .createHostCopy();          
    auto qr_host           = qr          .createHostCopy();          
    auto qi_host           = qi          .createHostCopy();          
    auto qs_host           = qs          .createHostCopy();          
    auto qg_host           = qg          .createHostCopy();          
    auto rho_dry_host      = rho_dry     .createHostCopy();     
    auto exner_host        = exner       .createHostCopy();   
    auto pressure_host     = pressure    .createHostCopy();
    auto dz_host           = dz          .createHostCopy();          
    auto rainnc_host       = rainnc      .createHostCopy();      
    auto rainncv_host      = rainncv     .createHostCopy();     
    auto snownc_host       = snownc      .createHostCopy();      
    auto snowncv_host      = snowncv     .createHostCopy();     
    auto sr_host           = sr          .createHostCopy();          
    auto graupelnc_host    = graupelnc   .createHostCopy();   
    auto graupelncv_host   = graupelncv  .createHostCopy();  

    int ids = 1; int ide = ncol; int jds = 1; int jde = 1; int kds = 1; int kde = nz;
    int ims = 1; int ime = ncol; int jms = 1; int jme = 1; int kms = 1; int kme = nz;
    int its = 1; int ite = ncol; int jts = 1; int jte = 1; int kts = 1; int kte = nz;

    wsm6( theta_host.data() , qv_host.data() , qc_host.data() , qr_host.data() , qi_host.data() , qs_host.data() ,
          qg_host.data() , rho_dry_host.data() , exner_host.data() , pressure_host.data() , dz_host.data() ,
          dt , grav , constants.cp_d , constants.cp_v , constants.R_d , constants.R_v , svpt0 , ep_1 , ep_2 , qmin ,
          xls , xlv , xlf , rhoair0 , rhowater , cliq , cice , psat , rainnc_host.data() , rainncv_host.data() ,
          snownc_host.data() , snowncv_host.data() , sr_host.data() , graupelnc_host.data() , graupelncv_host.data() ,
          ids , ide , jds , jde , kds , kde ,
          ims , ime , jms , jme , kms , kme ,
          its , ite , jts , jte , kts , kte );

    theta_host       .deep_copy_to( theta );   
    qv_host          .deep_copy_to( qv );          
    qc_host          .deep_copy_to( qc );          
    qr_host          .deep_copy_to( qr );          
    qi_host          .deep_copy_to( qi );          
    qs_host          .deep_copy_to( qs );          
    qg_host          .deep_copy_to( qg );          
    rho_dry_host     .deep_copy_to( rho_dry );     
    exner_host       .deep_copy_to( exner );   
    pressure_host    .deep_copy_to( pressure );
    dz_host          .deep_copy_to( dz );          
    rainnc_host      .deep_copy_to( rainnc );      
    rainncv_host     .deep_copy_to( rainncv );     
    snownc_host      .deep_copy_to( snownc );      
    snowncv_host     .deep_copy_to( snowncv );     
    sr_host          .deep_copy_to( sr );          
    graupelnc_host   .deep_copy_to( graupelnc );   
    graupelncv_host  .deep_copy_to( graupelncv );  

    parallel_for( Bounds<2>(nz,ncol) , YAKL_LAMBDA (int k, int i) {
      rho_v(k,i) = qv(k,i)*rho_dry(k,i);
      rho_c(k,i) = qc(k,i)*rho_dry(k,i);
      rho_r(k,i) = qr(k,i)*rho_dry(k,i);
      rho_i(k,i) = qi(k,i)*rho_dry(k,i);
      rho_s(k,i) = qs(k,i)*rho_dry(k,i);
      rho_r(k,i) = qg(k,i)*rho_dry(k,i);
      temp (k,i) = theta(k,i) * exner(k,i);
    });

    #ifdef PAM_DEBUG
      validate_array_positive(rho_v);
      validate_array_positive(rho_c);
      validate_array_positive(rho_r);
      validate_array_positive(rho_i);
      validate_array_positive(rho_s);
      validate_array_positive(rho_g);
      real mass_final = compute_total_mass( dm );
      real reldiff = abs(mass_final - mass_init) / ( abs(mass_init) + 1.e-20 );
      if ( reldiff > 1.e-13 ) {
        std::cout << "Microphysics mass change is too large: " << reldiff << std::endl;
        // endrun("ERROR: mass not conserved by kessler microphysics");
      }
    #endif

  }



  void output(DataManager &dm, yakl::SimpleNetCDF &nc, int ulIndex, int iens) const {
    auto rainnc     = dm.get<real,3>("rainnc"    );
    auto rainncv    = dm.get<real,3>("rainncv"   );
    auto snownc     = dm.get<real,3>("snownc"    );
    auto snowncv    = dm.get<real,3>("snowncv"   );
    auto sr         = dm.get<real,3>("sr"        );
    auto graupelnc  = dm.get<real,3>("graupelnc" );
    auto graupelncv = dm.get<real,3>("graupelncv");
    int nx = dm.get_dimension_size("x");
    int ny = dm.get_dimension_size("y");
    real2d data("data",ny,nx);

    parallel_for( Bounds<2>(ny,nx) , YAKL_LAMBDA (int j, int i) { data(j,i) = rainnc(j,i,iens); });
    nc.write1(data.createHostCopy(),"rainnc",{"y","x"},ulIndex,"t");

    parallel_for( Bounds<2>(ny,nx) , YAKL_LAMBDA (int j, int i) { data(j,i) = rainncv(j,i,iens); });
    nc.write1(data.createHostCopy(),"rainncv",{"y","x"},ulIndex,"t");

    parallel_for( Bounds<2>(ny,nx) , YAKL_LAMBDA (int j, int i) { data(j,i) = snownc(j,i,iens); });
    nc.write1(data.createHostCopy(),"snownc",{"y","x"},ulIndex,"t");

    parallel_for( Bounds<2>(ny,nx) , YAKL_LAMBDA (int j, int i) { data(j,i) = snowncv(j,i,iens); });
    nc.write1(data.createHostCopy(),"snowncv",{"y","x"},ulIndex,"t");

    parallel_for( Bounds<2>(ny,nx) , YAKL_LAMBDA (int j, int i) { data(j,i) = sr(j,i,iens); });
    nc.write1(data.createHostCopy(),"sr",{"y","x"},ulIndex,"t");

    parallel_for( Bounds<2>(ny,nx) , YAKL_LAMBDA (int j, int i) { data(j,i) = graupelnc(j,i,iens); });
    nc.write1(data.createHostCopy(),"graupelnc",{"y","x"},ulIndex,"t");

    parallel_for( Bounds<2>(ny,nx) , YAKL_LAMBDA (int j, int i) { data(j,i) = graupelncv(j,i,iens); });
    nc.write1(data.createHostCopy(),"graupelncv",{"y","x"},ulIndex,"t");
  }



  std::string micro_name() const {
    return "wsm6";
  }




};



