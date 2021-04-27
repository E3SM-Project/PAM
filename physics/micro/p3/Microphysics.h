
#pragma once

#include "awfl_const.h"
#include "DataManager.h"


extern "C"
void p3_main_fortran(double *qc , double *nc , double *qr , double *nr , double *th_atm , double *qv ,
                     double &dt , double *qi , double *qm , double *ni , double *bm , double *pres ,
                     double *dz , double *nc_nuceat_tend , double *nccn_prescribed , double *ni_activated ,
                     double *inv_qc_relvar , int &it , double *precip_liq_surf , double *precip_ice_surf ,
                     int &its , int &ite , int &kts , int &kte , double *diag_eff_radius_qc ,
                     double *diag_eff_radius_qi , double *rho_qi , bool &do_predict_nc , 
                     bool &do_prescribed_CCN ,double *dpres , double *exner , double *qv2qi_depos_tend ,
                     double *precip_total_tend , double *nevapr , double *qr_evap_tend ,
                     double *precip_liq_flux , double *precip_ice_flux , double *cld_frac_r ,
                     double *cld_frac_l , double *cld_frac_i , double *p3_tend_out , double *mu_c ,
                     double *lamc , double *liq_ice_exchange , double *vap_liq_exchange , 
                     double *vap_ice_exchange , double *qv_prev , double *t_prev , double *col_location ,
                     double *elapsed_s );



extern "C"
void micro_p3_utils_init_fortran(real &cpair , real &rair , real &rh2o , real &rhoh2o , real &mwh2o ,
                                 real &mwdry , real &gravit , real &latvap , real &latice , real &cpliq ,
                                 real &tmelt , real &pi , int &iulog , bool &masterproc );



extern "C"
void p3_init_fortran(char const *lookup_file_dir , int &dir_len , char const *version_p3 , int &ver_len );



class Microphysics {
public:
  // Doesn't actually have to be static or constexpr. Could be assigned in the constructor
  int static constexpr num_tracers = 9;

  // You should set these int he constructor
  struct Constants {
    real R_d    ;
    real cp_d   ;
    real cv_d   ;
    real gamma_d;
    real kappa_d;
    real R_v    ;
    real cp_v   ;
    real cv_v   ;
    real p0     ;
  };

  real grav;
  real cp_l;

  // This must be set during init() so we can return it in the get_water_vapor_index function
  int tracer_index_vapor;

  bool first_step;

  Constants constants;

  SArray<int,1,num_tracers> tracer_IDs; // tracer index for microphysics tracers

  // Indices for all of your tracer quantities
  int static constexpr ID_C  = 0;  // Local index for Cloud Water Mass  
  int static constexpr ID_NC = 1;  // Local index for Cloud Water Number
  int static constexpr ID_R  = 2;  // Local index for Rain Water Mass   
  int static constexpr ID_NR = 3;  // Local index for Rain Water Number 
  int static constexpr ID_I  = 4;  // Local index for Ice Mass          
  int static constexpr ID_M  = 5;  // Local index for Ice Number        
  int static constexpr ID_NI = 6;  // Local index for Ice-Rime Mass     
  int static constexpr ID_BM = 7;  // Local index for Ice-Rime Volume   
  int static constexpr ID_V  = 8;  // Local index for Water Vapor       



  // TODO: Make sure the constants vibe with P3
  // Set constants and likely num_tracers as well, and anything else you can do immediately
  Microphysics() {
    constants.R_d     = 287.;
    constants.cp_d    = 1003.;
    constants.cv_d    = constants.cp_d - constants.R_d;
    constants.gamma_d = constants.cp_d / constants.cv_d;
    constants.kappa_d = constants.R_d  / constants.cp_d;
    constants.R_v     = 461.;
    constants.cp_v    = 1859;
    constants.cv_v    = constants.R_v - constants.cp_v;
    constants.p0      = 1.e5;
    grav              = 9.81;
    first_step        = true;
    cp_l              = 4218.;
  }



  // This must return the correct # of tracers **BEFORE** init(...) is called
  int get_num_tracers() const {
    return num_tracers;
  }



  // This must return the correct index of water vapor **AFTER** init(...) is called
  int get_water_vapor_index() const {
    return tracer_index_vapor;
  }



  // Can do whatever you want, but mainly for registering tracers and allocating data
  // and storing the water vapor tracer index
  template <class DC>
  void init(std::string infile , int ny, int nx, int nens , DC &dycore , DataManager &dm) {
    int nz = dm.get_dimension_size("z");

    // Register tracers in the dycore
    //                                        name                 description            positive   adds mass
    tracer_IDs(ID_C ) = dycore.add_tracer(dm , "cloud_water"     , "Cloud Water Mass"   , true     , true );
    tracer_IDs(ID_NC) = dycore.add_tracer(dm , "cloud_water_num" , "Cloud Water Number" , true     , false);
    tracer_IDs(ID_R ) = dycore.add_tracer(dm , "rain"            , "Rain Water Mass"    , true     , true );
    tracer_IDs(ID_NR) = dycore.add_tracer(dm , "rain_num"        , "Rain Water Number"  , true     , false);
    tracer_IDs(ID_I ) = dycore.add_tracer(dm , "ice"             , "Ice Mass"           , true     , true );
    tracer_IDs(ID_NI) = dycore.add_tracer(dm , "ice_num"         , "Ice Number"         , true     , false);
    tracer_IDs(ID_M ) = dycore.add_tracer(dm , "ice_rime"        , "Ice-Rime Mass"      , true     , false);
    tracer_IDs(ID_BM) = dycore.add_tracer(dm , "ice_rime_vol"    , "Ice-Rime Volume"    , true     , false);
    tracer_IDs(ID_V ) = dycore.add_tracer(dm , "water_vapor"     , "Water Vapor"        , true     , true );

    // Register and allocate the tracers in the DataManager
    dm.register_and_allocate<real>( "cloud_water"     , "Cloud Water Mass"   , {nz,ny,nx,nens} , {"z","y","x","nens"} );
    dm.register_and_allocate<real>( "cloud_water_num" , "Cloud Water Number" , {nz,ny,nx,nens} , {"z","y","x","nens"} );
    dm.register_and_allocate<real>( "rain"            , "Rain Water Mass"    , {nz,ny,nx,nens} , {"z","y","x","nens"} );
    dm.register_and_allocate<real>( "rain_num"        , "Rain Water Number"  , {nz,ny,nx,nens} , {"z","y","x","nens"} );
    dm.register_and_allocate<real>( "ice"             , "Ice Mass"           , {nz,ny,nx,nens} , {"z","y","x","nens"} );
    dm.register_and_allocate<real>( "ice_num"         , "Ice Number"         , {nz,ny,nx,nens} , {"z","y","x","nens"} );
    dm.register_and_allocate<real>( "ice_rime"        , "Ice-Rime Mass"      , {nz,ny,nx,nens} , {"z","y","x","nens"} );
    dm.register_and_allocate<real>( "ice_rime_vol"    , "Ice-Rime Volume"    , {nz,ny,nx,nens} , {"z","y","x","nens"} );
    dm.register_and_allocate<real>( "water_vapor"     , "Water Vapor"        , {nz,ny,nx,nens} , {"z","y","x","nens"} );

    tracer_index_vapor = tracer_IDs(ID_V);

    // Register and allocation non-tracer quantities used by the microphysics
    int p3_nout = 49;
    dm.register_and_allocate<real>( "precip_liq_surf"    , "precipitation rate, liquid       m s-1"              , {             ny,nx,nens} , {                "y","x","nens"} );
    dm.register_and_allocate<real>( "precip_ice_surf"    , "precipitation rate, solid        m s-1"              , {             ny,nx,nens} , {                "y","x","nens"} );
    dm.register_and_allocate<real>( "diag_eff_radius_qc" , "effective radius, cloud          m"                  , {        nz  ,ny,nx,nens} , {          "z"  ,"y","x","nens"} );
    dm.register_and_allocate<real>( "diag_eff_radius_qi" , "effective radius, ice            m"                  , {        nz  ,ny,nx,nens} , {          "z"  ,"y","x","nens"} );
    dm.register_and_allocate<real>( "bulk_qi"            , "bulk density of ice              kg m-3"             , {        nz  ,ny,nx,nens} , {          "z"  ,"y","x","nens"} );
    dm.register_and_allocate<real>( "mu_c"               , "Size distribution shape parameter for radiation"     , {        nz  ,ny,nx,nens} , {          "z"  ,"y","x","nens"} );
    dm.register_and_allocate<real>( "lamc"               , "Size distribution slope parameter for radiation"     , {        nz  ,ny,nx,nens} , {          "z"  ,"y","x","nens"} );
    dm.register_and_allocate<real>( "qv2qi_depos_tend"   , "qitend due to deposition/sublimation"                , {        nz  ,ny,nx,nens} , {          "z"  ,"y","x","nens"} );
    dm.register_and_allocate<real>( "precip_total_tend"  , "Total precipitation (rain + snow)"                   , {        nz  ,ny,nx,nens} , {          "z"  ,"y","x","nens"} );
    dm.register_and_allocate<real>( "nevapr"             , "evaporation of total precipitation (rain + snow)"    , {        nz  ,ny,nx,nens} , {          "z"  ,"y","x","nens"} );
    dm.register_and_allocate<real>( "qr_evap_tend"       , "evaporation of rain"                                 , {        nz  ,ny,nx,nens} , {          "z"  ,"y","x","nens"} );
    dm.register_and_allocate<real>( "precip_liq_flux"    , "grid-box average rain flux (kg m^-2 s^-1) pverp"     , {        nz+1,ny,nx,nens} , {          "zp1","y","x","nens"} );
    dm.register_and_allocate<real>( "precip_ice_flux"    , "grid-box average ice/snow flux (kg m^-2 s^-1) pverp" , {        nz+1,ny,nx,nens} , {          "zp1","y","x","nens"} );
    dm.register_and_allocate<real>( "liq_ice_exchange"   , "sum of liq-ice phase change tendenices"              , {        nz  ,ny,nx,nens} , {          "z"  ,"y","x","nens"} );
    dm.register_and_allocate<real>( "vap_liq_exchange"   , "sum of vap-liq phase change tendenices"              , {        nz  ,ny,nx,nens} , {          "z"  ,"y","x","nens"} );
    dm.register_and_allocate<real>( "vap_ice_exchange"   , "sum of vap-ice phase change tendenices"              , {        nz  ,ny,nx,nens} , {          "z"  ,"y","x","nens"} );
    dm.register_and_allocate<real>( "p3_tend_out"        , "micro physics tendencies"                            , {p3_nout,nz  ,ny,nx,nens} , {"p3_nout","z"  ,"y","x","nens"} );
    dm.register_and_allocate<real>( "qv_prev"            , "qv from the previous step"                           , {        nz  ,ny,nx,nens} , {          "z"  ,"y","x","nens"} );
    dm.register_and_allocate<real>( "t_prev"             , "Temperature from the previous step"                  , {        nz  ,ny,nx,nens} , {          "z"  ,"y","x","nens"} );

    auto cloud_water         = dm.get<real,4>( "cloud_water"        );
    auto cloud_water_num     = dm.get<real,4>( "cloud_water_num"    );
    auto rain                = dm.get<real,4>( "rain"               );
    auto rain_num            = dm.get<real,4>( "rain_num"           );
    auto ice                 = dm.get<real,4>( "ice"                );
    auto ice_num             = dm.get<real,4>( "ice_num"            );
    auto ice_rime            = dm.get<real,4>( "ice_rime"           );
    auto ice_rime_vol        = dm.get<real,4>( "ice_rime_vol"       );
    auto water_vapor         = dm.get<real,4>( "water_vapor"        );
    auto precip_liq_surf     = dm.get<real,3>( "precip_liq_surf"    );
    auto precip_ice_surf     = dm.get<real,3>( "precip_ice_surf"    );
    auto diag_eff_radius_qc  = dm.get<real,4>( "diag_eff_radius_qc" );
    auto diag_eff_radius_qi  = dm.get<real,4>( "diag_eff_radius_qi" );
    auto bulk_qi             = dm.get<real,4>( "bulk_qi"            );
    auto mu_c                = dm.get<real,4>( "mu_c"               );
    auto lamc                = dm.get<real,4>( "lamc"               );
    auto qv2qi_depos_tend    = dm.get<real,4>( "qv2qi_depos_tend"   );
    auto precip_total_tend   = dm.get<real,4>( "precip_total_tend"  );
    auto nevapr              = dm.get<real,4>( "nevapr"             );
    auto qr_evap_tend        = dm.get<real,4>( "qr_evap_tend"       );
    auto precip_liq_flux     = dm.get<real,4>( "precip_liq_flux"    );
    auto precip_ice_flux     = dm.get<real,4>( "precip_ice_flux"    );
    auto liq_ice_exchange    = dm.get<real,4>( "liq_ice_exchange"   );
    auto vap_liq_exchange    = dm.get<real,4>( "vap_liq_exchange"   );
    auto vap_ice_exchange    = dm.get<real,4>( "vap_ice_exchange"   );
    auto p3_tend_out         = dm.get<real,5>( "p3_tend_out"        );
    auto qv_prev             = dm.get<real,4>( "qv_prev"            );
    auto t_prev              = dm.get<real,4>( "t_prev"             );

    parallel_for( SimpleBounds<4>(nz,ny,nx,nens) , YAKL_LAMBDA (int k, int j, int i, int iens) {
      cloud_water       (k,j,i,iens) = 0;
      cloud_water_num   (k,j,i,iens) = 0;
      rain              (k,j,i,iens) = 0;
      rain_num          (k,j,i,iens) = 0;
      ice               (k,j,i,iens) = 0;
      ice_num           (k,j,i,iens) = 0;
      ice_rime          (k,j,i,iens) = 0;
      ice_rime_vol      (k,j,i,iens) = 0;
      water_vapor       (k,j,i,iens) = 0;
      precip_liq_surf   (  j,i,iens) = 0;
      precip_ice_surf   (  j,i,iens) = 0;
      diag_eff_radius_qc(k,j,i,iens) = 0;
      diag_eff_radius_qi(k,j,i,iens) = 0;
      bulk_qi           (k,j,i,iens) = 0;
      mu_c              (k,j,i,iens) = 0;
      lamc              (k,j,i,iens) = 0;
      qv2qi_depos_tend  (k,j,i,iens) = 0;
      precip_total_tend (k,j,i,iens) = 0;
      nevapr            (k,j,i,iens) = 0;
      qr_evap_tend      (k,j,i,iens) = 0;
      precip_liq_flux   (k,j,i,iens) = 0;
      precip_ice_flux   (k,j,i,iens) = 0;
      liq_ice_exchange  (k,j,i,iens) = 0;
      vap_liq_exchange  (k,j,i,iens) = 0;
      vap_ice_exchange  (k,j,i,iens) = 0;
      qv_prev           (k,j,i,iens) = 0;
      t_prev            (k,j,i,iens) = 0;

      if (k == nz-1) {
        precip_liq_flux(nz,j,i,iens) = 0;
        precip_ice_flux(nz,j,i,iens) = 0;
      }

      for (int l=0; l < p3_nout; l++) {
        p3_tend_out(l,k,j,i,iens) = 0;
      }
    });

    real rhoh2o = 1000.;
    real mwdry  = 28.966;
    real mwh2o  = 0.622 * mwdry;
    real latvap = 2.5E6;
    real latice = 3.50E5;
    real tmelt  = 273.15;
    real pi     = M_PI;
    int  iulog  = 1;
    bool masterproc = true;
    micro_p3_utils_init_fortran( constants.cp_d , constants.R_d , constants.R_v , rhoh2o , mwh2o , mwdry ,
                                 grav , latvap , latice, cp_l , tmelt , pi , iulog , masterproc );

    std::string dir = "../../physics/micro/p3";
    std::string ver = "4";
    int dir_len = dir.length();
    int ver_len = ver.length();
    p3_init_fortran( dir.c_str() , dir_len , ver.c_str() , ver_len );
  }



  void timeStep( DataManager &dm , real dt ) {

    // Get the dimensions sizes
    int nz   = dm.get_dimension_size("z"   );
    int ny   = dm.get_dimension_size("y"   );
    int nx   = dm.get_dimension_size("x"   );
    int nens = dm.get_dimension_size("nens");
    int ncol = ny*nx*nens;

    // Get tracers dimensioned as (nz,ny*nx*nens)
    auto rho_c  = dm.get_lev_col<real>("cloud_water"    );
    auto rho_nc = dm.get_lev_col<real>("cloud_water_num");
    auto rho_r  = dm.get_lev_col<real>("rain"           );
    auto rho_nr = dm.get_lev_col<real>("rain_num"       );
    auto rho_i  = dm.get_lev_col<real>("ice"            );
    auto rho_ni = dm.get_lev_col<real>("ice_num"        );
    auto rho_m  = dm.get_lev_col<real>("ice_rime"       );
    auto rho_bm = dm.get_lev_col<real>("ice_rime_vol"   );
    auto rho_v  = dm.get_lev_col<real>("water_vapor"    );

    // Get coupler state
    auto rho_dry      = dm.get_lev_col<real>("density_dry");
    auto temp         = dm.get_lev_col<real>("temp");
    auto pressure_dry = dm.get_lev_col<real>("pressure_dry");

    // Calculate the grid spacing
    auto zint_in = dm.get<real,2>("vertical_interface_height");
    real2d dz("dz",nz,ny*nx*nens);
    parallel_for( Bounds<4>(nz,ny,nx,nens) , YAKL_LAMBDA (int k, int j, int i, int iens) {
      dz(k,j*nx*nens + i*nens + iens) = zint_in(k+1,iens) - zint_in(k,iens);
    });

    // Get everything from the DataManager that's not a tracer but is persistent across multiple micro calls
    auto precip_liq_surf    = dm.get_collapsed<real>("precip_liq_surf"   );
    auto precip_ice_surf    = dm.get_collapsed<real>("precip_ice_surf"   );
    auto diag_eff_radius_qc = dm.get_lev_col  <real>("diag_eff_radius_qc");
    auto diag_eff_radius_qi = dm.get_lev_col  <real>("diag_eff_radius_qi");
    auto bulk_qi            = dm.get_lev_col  <real>("bulk_qi"           );
    auto mu_c               = dm.get_lev_col  <real>("mu_c"              );
    auto lamc               = dm.get_lev_col  <real>("lamc"              );
    auto qv2qi_depos_tend   = dm.get_lev_col  <real>("qv2qi_depos_tend"  );
    auto precip_total_tend  = dm.get_lev_col  <real>("precip_total_tend" );
    auto nevapr             = dm.get_lev_col  <real>("nevapr"            );
    auto qr_evap_tend       = dm.get_lev_col  <real>("qr_evap_tend"      );
    auto liq_ice_exchange   = dm.get_lev_col  <real>("liq_ice_exchange"  );
    auto vap_liq_exchange   = dm.get_lev_col  <real>("vap_liq_exchange"  );
    auto vap_ice_exchange   = dm.get_lev_col  <real>("vap_ice_exchange"  );
    auto qv_prev            = dm.get_lev_col  <real>("qv_prev"           );
    auto t_prev             = dm.get_lev_col  <real>("t_prev"            );

    auto precip_liq_flux_dm = dm.get<real,4>("precip_liq_flux");
    auto precip_ice_flux_dm = dm.get<real,4>("precip_ice_flux");
    real2d precip_liq_flux("precip_liq_flux",nz+1,ny*nx*nens);
    real2d precip_ice_flux("precip_ice_flux",nz+1,ny*nx*nens);
    parallel_for( Bounds<4>(nz+1,ny,nx,nens) , YAKL_LAMBDA (int k, int j, int i, int iens) {
      precip_liq_flux(k,nx*nens*j + nens*i + iens) = precip_liq_flux_dm(k,j,i,iens);
      precip_ice_flux(k,nx*nens*j + nens*i + iens) = precip_ice_flux_dm(k,j,i,iens);
    });

    auto p3_tend_out_dm = dm.get<real,5>("p3_tend_out");
    real3d p3_tend_out("p3_tend_out",49,nz,ny*nx*nens);
    parallel_for( Bounds<5>(49,nz,ny,nx,nens) , YAKL_LAMBDA (int l, int k, int j, int i, int iens) {
      p3_tend_out(l,k,nx*nens*j + nens*i + iens) = p3_tend_out_dm(l,k,j,i,iens);
    });

    // These are inputs to p3
    real2d qc             ("qc"             ,nz,ncol);
    real2d nc             ("nc"             ,nz,ncol);
    real2d qr             ("qr"             ,nz,ncol);
    real2d nr             ("nr"             ,nz,ncol);
    real2d qi             ("qi"             ,nz,ncol);
    real2d ni             ("ni"             ,nz,ncol);
    real2d qm             ("qm"             ,nz,ncol);
    real2d bm             ("bm"             ,nz,ncol);
    real2d qv             ("qv"             ,nz,ncol);
    real2d pressure       ("pressure"       ,nz,ncol);
    real2d theta          ("theta"          ,nz,ncol);
    real2d exner          ("exner"          ,nz,ncol);
    real2d inv_exner      ("inv_exner"      ,nz,ncol);
    real2d dpres          ("dpres"          ,nz,ncol);
    real2d nc_nuceat_tend ("nc_nuceat_tend" ,nz,ncol);
    real2d nccn_prescribed("nccn_prescribed",nz,ncol);
    real2d ni_activated   ("ni_activated"   ,nz,ncol);
    real2d cld_frac_i     ("cld_frac_i"     ,nz,ncol);
    real2d cld_frac_l     ("cld_frac_l"     ,nz,ncol);
    real2d cld_frac_r     ("cld_frac_r"     ,nz,ncol);
    real2d inv_qc_relvar  ("inv_qc_relvar"  ,nz,ncol);
    real2d col_location   ("col_location"   ,3 ,ncol);

    //////////////////////////////////////////////////////////////////////////////
    // Compute quantities needed for inputs to P3
    //////////////////////////////////////////////////////////////////////////////
    // Force constants into local scope
    real gamma_d = this->constants.gamma_d;
    real R_d     = this->constants.R_d;
    real R_v     = this->constants.R_v;
    real cp_d    = this->constants.cp_d;
    real cp_v    = this->constants.cp_v;
    real cp_l    = this->cp_l;
    real p0      = this->constants.p0;

    YAKL_SCOPE( first_step , this->first_step );

    // Save initial state, and compute inputs for kessler(...)
    parallel_for( Bounds<2>(nz,ncol) , YAKL_LAMBDA (int k, int i) {
      // Compute total density
      real rho = rho_dry(k,i) + rho_c(k,i) + rho_r(k,i) + rho_i(k,i) + rho_v(k,i);

      // P3 doesn't do saturation adjustment, so we need to do that ahead of time
      compute_adjusted_state(rho, rho_dry(k,i) , rho_v(k,i) , rho_c(k,i) , temp(k,i),
                             R_v , cp_d , cp_v , cp_l);

      // Compute quantities for P3
      qc       (k,i) = rho_c (k,i) / rho_dry(k,i);
      nc       (k,i) = rho_nc(k,i) / rho_dry(k,i);
      qr       (k,i) = rho_r (k,i) / rho_dry(k,i);
      nr       (k,i) = rho_nr(k,i) / rho_dry(k,i);
      qi       (k,i) = rho_i (k,i) / rho_dry(k,i);
      ni       (k,i) = rho_ni(k,i) / rho_dry(k,i);
      qm       (k,i) = rho_m (k,i) / rho_dry(k,i);
      bm       (k,i) = rho_bm(k,i) / rho_dry(k,i);
      qv       (k,i) = rho_v (k,i) / rho_dry(k,i);
      pressure (k,i) = pressure_dry(k,i) + R_v * rho_v(k,i) * temp(k,i);
      exner    (k,i) = pow( pressure(k,i) / p0 , R_d / cp_d );
      inv_exner(k,i) = 1. / exner(k,i);
      theta    (k,i) = temp(k,i) / exner(k,i);
      // P3 uses dpres to calculate density via the hydrostatic assumption.
      // So we just reverse this to compute dpres to give true density
      dpres(k,i) = rho * grav * dz(k,i);
      // nc_nuceat_tend, nccn_prescribed, and ni_activated are not used
      nc_nuceat_tend (k,i) = 0;
      nccn_prescribed(k,i) = 0;
      ni_activated   (k,i) = 0;
      // Assume cloud fracton is always 1
      cld_frac_l(k,i) = 1;
      cld_frac_i(k,i) = 1;
      cld_frac_r(k,i) = 1;
      // inv_qc_relvar is always set to one
      inv_qc_relvar(k,i) = 1;
      // col_location is for debugging only, and it will be ignored for now
      if (k < 3) { col_location(k,i) = 1; }

      if (first_step) {
        qv_prev(k,i) = qv  (k,i);
        t_prev (k,i) = temp(k,i);
      }
    });

    first_step = false;

    int it = 1;
    int its = 1;
    int ite = ncol;
    int kts = 1;
    int kte = nz;
    bool do_predict_nc = false;
    bool do_prescribed_CCN = false;
    double elapsed_s;

    auto qc_host                 = qc                .createHostCopy();
    auto nc_host                 = nc                .createHostCopy();
    auto qr_host                 = qr                .createHostCopy();
    auto nr_host                 = nr                .createHostCopy();
    auto theta_host              = theta             .createHostCopy();
    auto qv_host                 = qv                .createHostCopy();
    auto qi_host                 = qi                .createHostCopy();
    auto qm_host                 = qm                .createHostCopy();
    auto ni_host                 = ni                .createHostCopy();
    auto bm_host                 = bm                .createHostCopy();
    auto pressure_host           = pressure          .createHostCopy();
    auto dz_host                 = dz                .createHostCopy();
    auto nc_nuceat_tend_host     = nc_nuceat_tend    .createHostCopy();
    auto nccn_prescribed_host    = nccn_prescribed   .createHostCopy();
    auto ni_activated_host       = ni_activated      .createHostCopy();
    auto inv_qc_relvar_host      = inv_qc_relvar     .createHostCopy();
    auto precip_liq_surf_host    = precip_liq_surf   .createHostCopy();
    auto precip_ice_surf_host    = precip_ice_surf   .createHostCopy();
    auto diag_eff_radius_qc_host = diag_eff_radius_qc.createHostCopy();
    auto diag_eff_radius_qi_host = diag_eff_radius_qi.createHostCopy();
    auto bulk_qi_host            = bulk_qi           .createHostCopy();
    auto dpres_host              = dpres             .createHostCopy();
    auto inv_exner_host          = inv_exner         .createHostCopy();
    auto qv2qi_depos_tend_host   = qv2qi_depos_tend  .createHostCopy();
    auto precip_total_tend_host  = precip_total_tend .createHostCopy();
    auto nevapr_host             = nevapr            .createHostCopy();
    auto qr_evap_tend_host       = qr_evap_tend      .createHostCopy();
    auto precip_liq_flux_host    = precip_liq_flux   .createHostCopy();
    auto precip_ice_flux_host    = precip_ice_flux   .createHostCopy();
    auto cld_frac_r_host         = cld_frac_r        .createHostCopy();
    auto cld_frac_l_host         = cld_frac_l        .createHostCopy();
    auto cld_frac_i_host         = cld_frac_i        .createHostCopy();
    auto p3_tend_out_host        = p3_tend_out       .createHostCopy();
    auto mu_c_host               = mu_c              .createHostCopy();
    auto lamc_host               = lamc              .createHostCopy();
    auto liq_ice_exchange_host   = liq_ice_exchange  .createHostCopy();
    auto vap_liq_exchange_host   = vap_liq_exchange  .createHostCopy();
    auto vap_ice_exchange_host   = vap_ice_exchange  .createHostCopy();
    auto qv_prev_host            = qv_prev           .createHostCopy();
    auto t_prev_host             = t_prev            .createHostCopy();
    auto col_location_host       = col_location      .createHostCopy();

    p3_main_fortran(qc_host.data() , nc_host.data() , qr_host.data() , nr_host.data() , theta_host.data() , qv_host.data() , dt , qi_host.data() ,
                    qm_host.data() , ni_host.data() , bm_host.data() , pressure_host.data() , dz_host.data() , nc_nuceat_tend_host.data() ,
                    nccn_prescribed_host.data() , ni_activated_host.data() , inv_qc_relvar_host.data() , it ,
                    precip_liq_surf_host.data() , precip_ice_surf_host.data() , its , ite , kts , kte ,
                    diag_eff_radius_qc_host.data() , diag_eff_radius_qi_host.data() , bulk_qi_host.data() , do_predict_nc ,
                    do_prescribed_CCN , dpres_host.data() , inv_exner_host.data() , qv2qi_depos_tend_host.data() ,
                    precip_total_tend_host.data() , nevapr_host.data() , qr_evap_tend_host.data() , precip_liq_flux_host.data() ,
                    precip_ice_flux_host.data() , cld_frac_r_host.data() , cld_frac_l_host.data() , cld_frac_i_host.data() , p3_tend_out_host.data() , mu_c_host.data() ,
                    lamc_host.data() , liq_ice_exchange_host.data() , vap_liq_exchange_host.data() , vap_ice_exchange_host.data() , 
                    qv_prev_host.data() , t_prev_host.data() , col_location_host.data() , &elapsed_s );

    qc_host                .deep_copy_to( qc                 );
    nc_host                .deep_copy_to( nc                 );
    qr_host                .deep_copy_to( qr                 );
    nr_host                .deep_copy_to( nr                 );
    theta_host             .deep_copy_to( theta              );
    qv_host                .deep_copy_to( qv                 );
    qi_host                .deep_copy_to( qi                 );
    qm_host                .deep_copy_to( qm                 );
    ni_host                .deep_copy_to( ni                 );
    bm_host                .deep_copy_to( bm                 );
    pressure_host          .deep_copy_to( pressure           );
    dz_host                .deep_copy_to( dz                 );
    nc_nuceat_tend_host    .deep_copy_to( nc_nuceat_tend     );
    nccn_prescribed_host   .deep_copy_to( nccn_prescribed    );
    ni_activated_host      .deep_copy_to( ni_activated       );
    inv_qc_relvar_host     .deep_copy_to( inv_qc_relvar      );
    precip_liq_surf_host   .deep_copy_to( precip_liq_surf    );
    precip_ice_surf_host   .deep_copy_to( precip_ice_surf    );
    diag_eff_radius_qc_host.deep_copy_to( diag_eff_radius_qc );
    diag_eff_radius_qi_host.deep_copy_to( diag_eff_radius_qi );
    bulk_qi_host           .deep_copy_to( bulk_qi            );
    dpres_host             .deep_copy_to( dpres              );
    inv_exner_host         .deep_copy_to( inv_exner          );
    qv2qi_depos_tend_host  .deep_copy_to( qv2qi_depos_tend   );
    precip_total_tend_host .deep_copy_to( precip_total_tend  );
    nevapr_host            .deep_copy_to( nevapr             );
    qr_evap_tend_host      .deep_copy_to( qr_evap_tend       );
    precip_liq_flux_host   .deep_copy_to( precip_liq_flux    );
    precip_ice_flux_host   .deep_copy_to( precip_ice_flux    );
    cld_frac_r_host        .deep_copy_to( cld_frac_r         );
    cld_frac_l_host        .deep_copy_to( cld_frac_l         );
    cld_frac_i_host        .deep_copy_to( cld_frac_i         );
    p3_tend_out_host       .deep_copy_to( p3_tend_out        );
    mu_c_host              .deep_copy_to( mu_c               );
    lamc_host              .deep_copy_to( lamc               );
    liq_ice_exchange_host  .deep_copy_to( liq_ice_exchange   );
    vap_liq_exchange_host  .deep_copy_to( vap_liq_exchange   );
    vap_ice_exchange_host  .deep_copy_to( vap_ice_exchange   );
    qv_prev_host           .deep_copy_to( qv_prev            );
    t_prev_host            .deep_copy_to( t_prev             );
    col_location_host      .deep_copy_to( col_location       );
                    
    ///////////////////////////////////////////////////////////////////////////////
    // Convert P3 outputs into dynamics coupler state and tracer masses
    ///////////////////////////////////////////////////////////////////////////////
    parallel_for( Bounds<2>(nz,ncol) , YAKL_LAMBDA (int k, int i) {
      rho_c  (k,i) = max( qc(k,i)*rho_dry(k,i) , 0._fp );
      rho_nc (k,i) = max( nc(k,i)*rho_dry(k,i) , 0._fp );
      rho_r  (k,i) = max( qr(k,i)*rho_dry(k,i) , 0._fp );
      rho_nr (k,i) = max( nr(k,i)*rho_dry(k,i) , 0._fp );
      rho_i  (k,i) = max( qi(k,i)*rho_dry(k,i) , 0._fp );
      rho_ni (k,i) = max( ni(k,i)*rho_dry(k,i) , 0._fp );
      rho_m  (k,i) = max( qm(k,i)*rho_dry(k,i) , 0._fp );
      rho_bm (k,i) = max( bm(k,i)*rho_dry(k,i) , 0._fp );
      rho_v  (k,i) = max( qv(k,i)*rho_dry(k,i) , 0._fp );
      temp   (k,i) = theta(k,i) * exner(k,i);
      // Save qv and temperature for the next call to p3_main
      qv_prev(k,i) = max( qv(k,i) , 0._fp );
      t_prev (k,i) = temp(k,i);
    });

  }



  // Returns saturation vapor pressure
  YAKL_INLINE real saturation_vapor_pressure(real temp) const {
    real tc = temp - 273.15;
    return 610.94 * exp( 17.625*tc / (243.04+tc) );
  }



  YAKL_INLINE real latent_heat_condensation(real temp) const {
    real tc = temp - 273.15;
    return (2500.8 - 2.36*tc + 0.0016*tc*tc - 0.00006*tc*tc*tc)*1000;
  }



  YAKL_INLINE real cp_moist(real rho_d, real rho_v, real rho_c, real cp_d, real cp_v, real cp_l) const {
    // For the moist specific heat, ignore other species than water vapor and cloud droplets
    real rho = rho_d + rho_v + rho_c;
    return rho_d / rho * cp_d + rho_v / rho * cp_v + rho_c / rho * cp_l;
  }



  // Compute an instantaneous adjustment of sub or super saturation
  YAKL_INLINE void compute_adjusted_state(real rho, real rho_d , real &rho_v , real &rho_c , real &temp,
                                          real R_v , real cp_d , real cp_v , real cp_l) const {
    // Define a tolerance for convergence
    real tol = 1.e-6;

    // Saturation vapor pressure at this temperature
    real svp = saturation_vapor_pressure( temp );

    // Vapor pressure at this temperature
    real pv = rho_v * R_v * temp;

    // If we're super-saturated, we need to condense until saturation is reached
    if        (pv > svp) {
      ////////////////////////////////////////////////////////
      // Bisection method
      ////////////////////////////////////////////////////////
      // Set bounds on how much mass to condense
      real cond1  = 0;     // Minimum amount we can condense out
      real cond2 = rho_v;  // Maximum amount we can condense out

      bool keep_iterating = true;
      while (keep_iterating) {
        real rho_cond = (cond1 + cond2) / 2;                    // How much water vapor to condense for this iteration
        real rv_loc = max( 0._fp , rho_v - rho_cond );          // New vapor density
        real rc_loc = max( 0._fp , rho_c + rho_cond );          // New cloud liquid density
        real Lv = latent_heat_condensation(temp);               // Compute latent heat of condensation
        real cp = cp_moist(rho_d,rv_loc,rc_loc,cp_d,cp_v,cp_l); // New moist specific heat at constant pressure
        real temp_loc = temp + rho_cond*Lv/(rho*cp);            // New temperature after condensation
        real svp_loc = saturation_vapor_pressure(temp_loc);     // New saturation vapor pressure after condensation
        real pv_loc = rv_loc * R_v * temp_loc;                  // New vapor pressure after condensation
        // If we're supersaturated still, we need to condense out more water vapor
        // otherwise, we need to condense out less water vapor
        if (pv_loc > svp_loc) {
          cond1 = rho_cond;
        } else {
          cond2 = rho_cond;
        }
        // If we've converged, then we can stop iterating
        if (abs(cond2-cond1) <= tol) {
          rho_v = rv_loc;
          rho_c = rc_loc;
          temp  = temp_loc;
          keep_iterating = false;
        }
      }

    // If we are unsaturated and have cloud liquid
    } else if (pv < svp && rho_c > 0) {
      // If there's cloud, evaporate enough to achieve saturation
      // or all of it if there isn't enough to reach saturation
      ////////////////////////////////////////////////////////
      // Bisection method
      ////////////////////////////////////////////////////////
      // Set bounds on how much mass to evaporate
      real evap1 = 0;     // minimum amount we can evaporate
      real evap2 = rho_c; // maximum amount we can evaporate

      bool keep_iterating = true;
      while (keep_iterating) {
        real rho_evap = (evap1 + evap2) / 2;                    // How much water vapor to evapense
        real rv_loc = max( 0._fp , rho_v + rho_evap );          // New vapor density
        real rc_loc = max( 0._fp , rho_c - rho_evap );          // New cloud liquid density
        real Lv = latent_heat_condensation(temp);               // Compute latent heat of condensation for water
        real cp = cp_moist(rho_d,rv_loc,rc_loc,cp_d,cp_v,cp_l); // New moist specific heat
        real temp_loc = temp - rho_evap*Lv/(rho*cp);            // New temperature after evaporation
        real svp_loc = saturation_vapor_pressure(temp_loc);     // New saturation vapor pressure after evaporation
        real pv_loc = rv_loc * R_v * temp_loc;                  // New vapor pressure after evaporation
        // If we're unsaturated still, we need to evaporate out more water vapor
        // otherwise, we need to evaporate out less water vapor
        if (pv_loc < svp_loc) {
          evap1 = rho_evap;
        } else {
          evap2 = rho_evap;
        }
        // If we've converged, then we can stop iterating
        if (abs(evap2-evap1) <= tol) {
          rho_v = rv_loc;
          rho_c = rc_loc;
          temp  = temp_loc;
          keep_iterating = false;
        }
      }
    }
  }




  // These are outputs that are not tracer mass. Tracer mass is handled by the dycore instead
  // This assumes the NetCDF handler "nc" is already open and will be closed later
  // This is for a single ensemble index
  void output(DataManager &dm, yakl::SimpleNetCDF &nc, int ulIndex, int iens) const {
  }



  std::string micro_name() const {
    return "p3";
  }



};



