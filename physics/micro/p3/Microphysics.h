
#pragma once

#include "pam_coupler.h"
// #include <stdio.h>

#include "scream_cxx_interface_p3.h"

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
  // You should set these in the constructor
  real R_d    ;
  real cp_d   ;
  real cv_d   ;
  real gamma_d;
  real kappa_d;
  real R_v    ;
  real cp_v   ;
  real cv_v   ;
  real p0     ;

  real grav;
  real cp_l;

  bool sgs_shoc;

  bool first_step;

  real etime;

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



  // Set constants and likely num_tracers as well, and anything else you can do immediately
  Microphysics() {
    R_d        = 287.042;
    cp_d       = 1004.64;
    cv_d       = cp_d - R_d;
    gamma_d    = cp_d / cv_d;
    kappa_d    = R_d  / cp_d;
    R_v        = 461.505;
    cp_v       = 1859;
    cv_v       = R_v - cp_v;
    p0         = 1.e5;
    grav       = 9.80616;
    first_step = true;
    cp_l       = 4188.0;
    sgs_shoc   = false;
  }



  // This must return the correct # of tracers **BEFORE** init(...) is called
  static int constexpr get_num_tracers() {
    return 9;
  }



  // Can do whatever you want, but mainly for registering tracers and allocating data
  void init(pam::PamCoupler &coupler) {
    using yakl::c::parallel_for;
    using yakl::c::SimpleBounds;

    int nx   = coupler.get_nx  ();
    int ny   = coupler.get_ny  ();
    int nz   = coupler.get_nz  ();
    int nens = coupler.get_nens();

    // Register tracers in the coupler
    //                 name                description            positive   adds mass
    coupler.add_tracer("cloud_water"     , "Cloud Water Mass"   , true     , true );
    coupler.add_tracer("cloud_water_num" , "Cloud Water Number" , true     , false);
    coupler.add_tracer("rain"            , "Rain Water Mass"    , true     , true );
    coupler.add_tracer("rain_num"        , "Rain Water Number"  , true     , false);
    coupler.add_tracer("ice"             , "Ice Mass"           , true     , true );
    coupler.add_tracer("ice_num"         , "Ice Number"         , true     , false);
    coupler.add_tracer("ice_rime"        , "Ice-Rime Mass"      , true     , false);
    coupler.add_tracer("ice_rime_vol"    , "Ice-Rime Volume"    , true     , false);
    coupler.add_tracer("water_vapor"     , "Water Vapor"        , true     , true );

    auto &dm = coupler.get_data_manager_device_readwrite();

    dm.register_and_allocate<real>("qv_prev","qv from prev step"         ,{nz,ny,nx,nens},{"z","y","x","nens"});
    dm.register_and_allocate<real>("t_prev" ,"Temperature from prev step",{nz,ny,nx,nens},{"z","y","x","nens"});

    dm.register_and_allocate<real>("precip_liq_surf_out","liq surface precipitation rate",{ny,nx,nens},{"y","x","nens"});
    dm.register_and_allocate<real>("precip_ice_surf_out","ice surface precipitation rate",{ny,nx,nens},{"y","x","nens"});

    dm.register_and_allocate<real>("liq_ice_exchange_out","p3 liq to ice phase change tendency",{nz,ny,nx,nens},{"z","y","x","nens"});
    dm.register_and_allocate<real>("vap_liq_exchange_out","p3 vap to liq phase change tendency",{nz,ny,nx,nens},{"z","y","x","nens"});
    dm.register_and_allocate<real>("vap_ice_exchange_out","p3 vap to ice phase change tendency",{nz,ny,nx,nens},{"z","y","x","nens"});

    auto cloud_water     = dm.get<real,4>( "cloud_water"     );
    auto cloud_water_num = dm.get<real,4>( "cloud_water_num" );
    auto rain            = dm.get<real,4>( "rain"            );
    auto rain_num        = dm.get<real,4>( "rain_num"        );
    auto ice             = dm.get<real,4>( "ice"             );
    auto ice_num         = dm.get<real,4>( "ice_num"         );
    auto ice_rime        = dm.get<real,4>( "ice_rime"        );
    auto ice_rime_vol    = dm.get<real,4>( "ice_rime_vol"    );
    auto water_vapor     = dm.get<real,4>( "water_vapor"     );
    auto qv_prev         = dm.get<real,4>( "qv_prev"         );
    auto t_prev          = dm.get<real,4>( "t_prev"          );

    parallel_for( "micro zero" , SimpleBounds<4>(nz,ny,nx,nens) ,
                  YAKL_LAMBDA (int k, int j, int i, int iens) {
      cloud_water    (k,j,i,iens) = 0;
      cloud_water_num(k,j,i,iens) = 0;
      rain           (k,j,i,iens) = 0;
      rain_num       (k,j,i,iens) = 0;
      ice            (k,j,i,iens) = 0;
      ice_num        (k,j,i,iens) = 0;
      ice_rime       (k,j,i,iens) = 0;
      ice_rime_vol   (k,j,i,iens) = 0;
      water_vapor    (k,j,i,iens) = 0;
      qv_prev        (k,j,i,iens) = 0;
      t_prev         (k,j,i,iens) = 0;
    });

    real rhoh2o = 1000.;
    real mwdry  = 28.966;
    real mwh2o  = 18.016;
    real latvap = 2501000.0;
    real latice = 333700.0;
    real tmelt  = 273.15;
    real pi     = 3.14159265;
    int  iulog  = 1;
    #ifndef P3_CXX
      bool masterproc = true;
      micro_p3_utils_init_fortran( cp_d , R_d , R_v , rhoh2o , mwh2o , mwdry ,
                                   grav , latvap , latice, cp_l , tmelt , pi , iulog , masterproc );
    #endif

    #ifndef P3_CXX
      std::string dir;
      if (coupler.option_exists("p3_lookup_data_path")) {
        dir = coupler.get_option<std::string>("p3_lookup_data_path");
      } else {
        dir = "../../../physics/micro/p3/tables"; // default for PAM standalone
      };
      std::string ver = "4.1.1";
      int dir_len = dir.length();
      int ver_len = ver.length();
      p3_init_fortran( dir.c_str() , dir_len , ver.c_str() , ver_len );
    #endif

    coupler.set_option<std::string>("micro","p3");
    coupler.set_option<real>("latvap",latvap);
    coupler.set_option<real>("latice",latice);
    coupler.set_option<real>("R_d" ,R_d  );
    coupler.set_option<real>("R_v" ,R_v  );
    coupler.set_option<real>("cp_d",cp_d );
    coupler.set_option<real>("cp_v",cp_v );
    coupler.set_option<real>("grav",grav );
    coupler.set_option<real>("p0"  ,p0   );

    etime = 0;
  }



  void timeStep( pam::PamCoupler &coupler ) {
    using yakl::c::parallel_for;
    using yakl::c::SimpleBounds;

    real dt = coupler.get_option<real>("crm_dt");

    if (first_step) {
      if (coupler.get_option<std::string>("sgs") == "shoc") sgs_shoc = true;
    }
    auto &dm = coupler.get_data_manager_device_readwrite();

    // Get the dimensions sizes
    int nz   = coupler.get_nz  ();
    int ny   = coupler.get_ny  ();
    int nx   = coupler.get_nx  ();
    int nens = coupler.get_nens();
    int ncol = ny*nx*nens;

    auto zint_in = dm.get<real,2>("vertical_interface_height");

    real crm_dx = coupler.get_xlen() / nx;
    real crm_dy = ny == 1 ? crm_dx : coupler.get_ylen() / ny;

    #ifdef PAM_DEBUG
      real mass0;
      {
        auto rho_v = dm.get<real,4>( "water_vapor" );
        auto rho_c = dm.get<real,4>( "cloud_water" );
        auto rho_r = dm.get<real,4>( "rain"        );
        auto rho_i = dm.get<real,4>( "ice"         );
        real4d mass4d("mass4d",nz,ny,nx,nens);
        parallel_for( SimpleBounds<4>(nz,ny,nx,nens) , YAKL_LAMBDA (int k, int j, int i, int iens) {
          mass4d(k,j,i,iens) = (rho_v(k,j,i,iens) + rho_c(k,j,i,iens) + rho_r(k,j,i,iens) + rho_i(k,j,i,iens)) *
                               crm_dx * crm_dy * (zint_in(k+1,iens) - zint_in(k,iens));
        });
        mass0 = yakl::intrinsics::sum(mass4d);
      }
    #endif

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
    auto rho_dry = dm.get_lev_col<real>("density_dry");
    auto temp    = dm.get_lev_col<real>("temp"       );

    // Calculate the grid spacing
    real2d dz("dz",nz,ny*nx*nens);
    parallel_for( "micro dz" , SimpleBounds<4>(nz,ny,nx,nens) ,
                  YAKL_LAMBDA (int k, int j, int i, int iens) {
      dz(k,j*nx*nens + i*nens + iens) = zint_in(k+1,iens) - zint_in(k,iens);
    });

    // Get everything from the DataManager that's not a tracer but is persistent across multiple micro calls
    auto qv_prev = dm.get_lev_col<real>("qv_prev");
    auto t_prev  = dm.get_lev_col<real>("t_prev" );

    // Allocates inputs and outputs
    real2d qc                ( "qc"                 ,           nz   , ncol );
    real2d nc                ( "nc"                 ,           nz   , ncol );
    real2d qr                ( "qr"                 ,           nz   , ncol );
    real2d nr                ( "nr"                 ,           nz   , ncol );
    real2d qi                ( "qi"                 ,           nz   , ncol );
    real2d ni                ( "ni"                 ,           nz   , ncol );
    real2d qm                ( "qm"                 ,           nz   , ncol );
    real2d bm                ( "bm"                 ,           nz   , ncol );
    real2d qv                ( "qv"                 ,           nz   , ncol );
    real2d pressure_dry      ( "pressure_dry"       ,           nz   , ncol );
    real2d theta             ( "theta"              ,           nz   , ncol );
    real2d exner             ( "exner"              ,           nz   , ncol );
    real2d inv_exner         ( "inv_exner"          ,           nz   , ncol );
    real2d dpres_dry         ( "dpres_dry"          ,           nz   , ncol );
    real2d nc_nuceat_tend    ( "nc_nuceat_tend"     ,           nz   , ncol );
    real2d nccn_prescribed   ( "nccn_prescribed"    ,           nz   , ncol );
    real2d ni_activated      ( "ni_activated"       ,           nz   , ncol );
    real2d cld_frac_i        ( "cld_frac_i"         ,           nz   , ncol );
    real2d cld_frac_l        ( "cld_frac_l"         ,           nz   , ncol );
    real2d cld_frac_r        ( "cld_frac_r"         ,           nz   , ncol );
    real2d inv_qc_relvar     ( "inv_qc_relvar"      ,           nz   , ncol );
    real2d col_location      ( "col_location"       ,           3    , ncol );
    real1d precip_liq_surf   ( "precip_liq_surf"    ,                  ncol );
    real1d precip_ice_surf   ( "precip_ice_surf"    ,                  ncol );
    real2d diag_eff_radius_qc( "diag_eff_radius_qc" ,           nz   , ncol );
    real2d diag_eff_radius_qi( "diag_eff_radius_qi" ,           nz   , ncol );
    real2d bulk_qi           ( "bulk_qi"            ,           nz   , ncol );
    real2d mu_c              ( "mu_c"               ,           nz   , ncol );
    real2d lamc              ( "lamc"               ,           nz   , ncol );
    real2d qv2qi_depos_tend  ( "qv2qi_depos_tend"   ,           nz   , ncol );
    real2d precip_total_tend ( "precip_total_tend"  ,           nz   , ncol );
    real2d nevapr            ( "nevapr"             ,           nz   , ncol );
    real2d qr_evap_tend      ( "qr_evap_tend"       ,           nz   , ncol );
    real2d precip_liq_flux   ( "precip_liq_flux"    ,           nz+1 , ncol );
    real2d precip_ice_flux   ( "precip_ice_flux"    ,           nz+1 , ncol );
    real2d liq_ice_exchange  ( "liq_ice_exchange"   ,           nz   , ncol );
    real2d vap_liq_exchange  ( "vap_liq_exchange"   ,           nz   , ncol );
    real2d vap_ice_exchange  ( "vap_ice_exchange"   ,           nz   , ncol );
    #ifndef P3_CXX
    int p3_nout = 49;
    real3d p3_tend_out       ( "p3_tend_out"        , p3_nout , nz   , ncol );
    #endif

    //////////////////////////////////////////////////////////////////////////////
    // Compute quantities needed for inputs to P3
    //////////////////////////////////////////////////////////////////////////////
    // Force constants into local scope
    real R_d     = this->R_d;
    real R_v     = this->R_v;
    real cp_d    = this->cp_d;
    real cp_v    = this->cp_v;
    real cp_l    = this->cp_l;
    real p0      = this->p0;

    YAKL_SCOPE( first_step , this->first_step );
    YAKL_SCOPE( sgs_shoc   , this->sgs_shoc   );
    YAKL_SCOPE( grav       , this->grav       );

    // Save initial state, and compute inputs for p3(...)
    parallel_for( "micro adjust preprocess" , SimpleBounds<2>(nz,ncol) , YAKL_LAMBDA (int k, int i) {

      // P3 doesn't do saturation adjustment, so we need to do that ahead of time
      // If we're using SHOC, then it does saturation adjustment, so no need to do it here
      if (! sgs_shoc) {
        // Compute total density
        real rho = rho_dry(k,i) + rho_c(k,i) + rho_r(k,i) + rho_i(k,i) + rho_v(k,i);
        compute_adjusted_state(rho, rho_dry(k,i) , rho_v(k,i) , rho_c(k,i) , temp(k,i),
                               R_v , cp_d , cp_v , cp_l);
      }

      // Compute quantities for P3
      qc          (k,i) = rho_c (k,i) / rho_dry(k,i);
      nc          (k,i) = rho_nc(k,i) / rho_dry(k,i);
      qr          (k,i) = rho_r (k,i) / rho_dry(k,i);
      nr          (k,i) = rho_nr(k,i) / rho_dry(k,i);
      qi          (k,i) = rho_i (k,i) / rho_dry(k,i);
      ni          (k,i) = rho_ni(k,i) / rho_dry(k,i);
      qm          (k,i) = rho_m (k,i) / rho_dry(k,i);
      bm          (k,i) = rho_bm(k,i) / rho_dry(k,i);
      qv          (k,i) = rho_v (k,i) / rho_dry(k,i);
      real pressure     = R_d*rho_dry(k,i)*temp(k,i) + R_v*rho_v(k,i)*temp(k,i);
      pressure_dry(k,i) = R_d*rho_dry(k,i)*temp(k,i);
      exner       (k,i) = pow( pressure / p0 , R_d / cp_d );
      inv_exner   (k,i) = 1. / exner(k,i);
      theta       (k,i) = temp(k,i) / exner(k,i);
      // P3 uses dpres to calculate density via the hydrostatic assumption.
      // So we just reverse this to compute dpres to give true density
      dpres_dry(k,i) = rho_dry(k,i) * grav * dz(k,i);
      // nc_nuceat_tend, nccn_prescribed, and ni_activated are not used
      nccn_prescribed(k,i) = 1e3;
      nc_nuceat_tend (k,i) = 1.0;
      ni_activated   (k,i) = 1.0;
      // col_location is for debugging only, and it will be ignored for now
      if (k < 3) { col_location(k,i) = 1; }

      if (first_step) {
        qv_prev(k,i) = qv  (k,i);
        t_prev (k,i) = temp(k,i);
      }
    });

    if (sgs_shoc) {
      inv_qc_relvar = dm.get_lev_col<real>("inv_qc_relvar");
      auto cld_frac = dm.get_lev_col<real>("cldfrac");
      get_cloud_fraction( cld_frac , qc , qr , qi , cld_frac_i , cld_frac_l , cld_frac_r );
      // parallel_for( SimpleBounds<2>(nz,ncol) , YAKL_LAMBDA (int k, int i) {
      //   cld_frac_l(k,i) = 1;
      //   cld_frac_i(k,i) = 1;
      //   cld_frac_r(k,i) = 1;
      // });
    } else {
      parallel_for( SimpleBounds<2>(nz,ncol) , YAKL_LAMBDA (int k, int i) {
        // Assume cloud fracton is always 1
        cld_frac_l(k,i) = 1;
        cld_frac_i(k,i) = 1;
        cld_frac_r(k,i) = 1;
        // inv_qc_relvar is always set to one
        inv_qc_relvar(k,i) = 1;
      });
    }
    double elapsed_s;
    int it, its, ite, kts, kte;
    bool do_predict_nc = false;
    bool do_prescribed_CCN = false;

    #ifdef P3_CXX

      it  = 0;
      its = 0;
      ite = ncol-1;
      kts = 0;
      kte = nz-1;

      // Create room for transposed variables (only 2-D variables need to be transposed)
      auto transposed_qc                 = qc                .createDeviceCopy().reshape(qc                .extent(1),qc                .extent(0)); // inout
      auto transposed_nc                 = nc                .createDeviceCopy().reshape(nc                .extent(1),nc                .extent(0)); // inout
      auto transposed_qr                 = qr                .createDeviceCopy().reshape(qr                .extent(1),qr                .extent(0)); // inout
      auto transposed_nr                 = nr                .createDeviceCopy().reshape(nr                .extent(1),nr                .extent(0)); // inout
      auto transposed_theta              = theta             .createDeviceCopy().reshape(theta             .extent(1),theta             .extent(0)); // inout
      auto transposed_qv                 = qv                .createDeviceCopy().reshape(qv                .extent(1),qv                .extent(0)); // inout
      auto transposed_qi                 = qi                .createDeviceCopy().reshape(qi                .extent(1),qi                .extent(0)); // inout
      auto transposed_qm                 = qm                .createDeviceCopy().reshape(qm                .extent(1),qm                .extent(0)); // inout
      auto transposed_ni                 = ni                .createDeviceCopy().reshape(ni                .extent(1),ni                .extent(0)); // inout
      auto transposed_bm                 = bm                .createDeviceCopy().reshape(bm                .extent(1),bm                .extent(0)); // inout
      auto transposed_pressure_dry       = pressure_dry      .createDeviceCopy().reshape(pressure_dry      .extent(1),pressure_dry      .extent(0)); // in
      auto transposed_dz                 = dz                .createDeviceCopy().reshape(dz                .extent(1),dz                .extent(0)); // in
      auto transposed_nc_nuceat_tend     = nc_nuceat_tend    .createDeviceCopy().reshape(nc_nuceat_tend    .extent(1),nc_nuceat_tend    .extent(0)); // in
      auto transposed_nccn_prescribed    = nccn_prescribed   .createDeviceCopy().reshape(nccn_prescribed   .extent(1),nccn_prescribed   .extent(0)); // in
      auto transposed_ni_activated       = ni_activated      .createDeviceCopy().reshape(ni_activated      .extent(1),ni_activated      .extent(0)); // in
      auto transposed_inv_qc_relvar      = inv_qc_relvar     .createDeviceCopy().reshape(inv_qc_relvar     .extent(1),inv_qc_relvar     .extent(0)); // in
      auto transposed_dpres_dry          = dpres_dry         .createDeviceCopy().reshape(dpres_dry         .extent(1),dpres_dry         .extent(0)); // in
      auto transposed_inv_exner          = inv_exner         .createDeviceCopy().reshape(inv_exner         .extent(1),inv_exner         .extent(0)); // in
      auto transposed_cld_frac_r         = cld_frac_r        .createDeviceCopy().reshape(cld_frac_r        .extent(1),cld_frac_r        .extent(0)); // in
      auto transposed_cld_frac_l         = cld_frac_l        .createDeviceCopy().reshape(cld_frac_l        .extent(1),cld_frac_l        .extent(0)); // in
      auto transposed_cld_frac_i         = cld_frac_i        .createDeviceCopy().reshape(cld_frac_i        .extent(1),cld_frac_i        .extent(0)); // in
      auto transposed_qv_prev            = qv_prev           .createDeviceCopy().reshape(qv_prev           .extent(1),qv_prev           .extent(0)); // in
      auto transposed_t_prev             = t_prev            .createDeviceCopy().reshape(t_prev            .extent(1),t_prev            .extent(0)); // in
      auto transposed_col_location       = col_location      .createDeviceCopy().reshape(col_location      .extent(1),col_location      .extent(0)); // in
      auto transposed_diag_eff_radius_qc = diag_eff_radius_qc.createDeviceCopy().reshape(diag_eff_radius_qc.extent(1),diag_eff_radius_qc.extent(0)); //   out
      auto transposed_diag_eff_radius_qi = diag_eff_radius_qi.createDeviceCopy().reshape(diag_eff_radius_qi.extent(1),diag_eff_radius_qi.extent(0)); //   out
      auto transposed_bulk_qi            = bulk_qi           .createDeviceCopy().reshape(bulk_qi           .extent(1),bulk_qi           .extent(0)); //   out
      auto transposed_qv2qi_depos_tend   = qv2qi_depos_tend  .createDeviceCopy().reshape(qv2qi_depos_tend  .extent(1),qv2qi_depos_tend  .extent(0)); //   out
      auto transposed_precip_liq_flux    = precip_liq_flux   .createDeviceCopy().reshape(precip_liq_flux   .extent(1),precip_liq_flux   .extent(0)); //   out
      auto transposed_precip_ice_flux    = precip_ice_flux   .createDeviceCopy().reshape(precip_ice_flux   .extent(1),precip_ice_flux   .extent(0)); //   out
      auto transposed_liq_ice_exchange   = liq_ice_exchange  .createDeviceCopy().reshape(liq_ice_exchange  .extent(1),liq_ice_exchange  .extent(0)); //   out
      auto transposed_vap_liq_exchange   = vap_liq_exchange  .createDeviceCopy().reshape(vap_liq_exchange  .extent(1),vap_liq_exchange  .extent(0)); //   out
      auto transposed_vap_ice_exchange   = vap_ice_exchange  .createDeviceCopy().reshape(vap_ice_exchange  .extent(1),vap_ice_exchange  .extent(0)); //   out

      // For in and inout variables, copy transposed data (One kernel for efficiency)
      parallel_for( SimpleBounds<2>(nz,ncol) , YAKL_LAMBDA (int k, int i) {
        int k_p3 = nz-1-k;
        transposed_qc             (i,k_p3) = qc             (k,i); // inout
        transposed_nc             (i,k_p3) = nc             (k,i); // inout
        transposed_qr             (i,k_p3) = qr             (k,i); // inout
        transposed_nr             (i,k_p3) = nr             (k,i); // inout
        transposed_theta          (i,k_p3) = theta          (k,i); // inout
        transposed_qv             (i,k_p3) = qv             (k,i); // inout
        transposed_qi             (i,k_p3) = qi             (k,i); // inout
        transposed_qm             (i,k_p3) = qm             (k,i); // inout
        transposed_ni             (i,k_p3) = ni             (k,i); // inout
        transposed_bm             (i,k_p3) = bm             (k,i); // inout
        transposed_pressure_dry   (i,k_p3) = pressure_dry   (k,i); // in
        transposed_dz             (i,k_p3) = dz             (k,i); // in
        transposed_nc_nuceat_tend (i,k_p3) = nc_nuceat_tend (k,i); // in
        transposed_nccn_prescribed(i,k_p3) = nccn_prescribed(k,i); // in
        transposed_ni_activated   (i,k_p3) = ni_activated   (k,i); // in
        transposed_inv_qc_relvar  (i,k_p3) = inv_qc_relvar  (k,i); // in
        transposed_dpres_dry      (i,k_p3) = dpres_dry      (k,i); // in
        transposed_inv_exner      (i,k_p3) = inv_exner      (k,i); // in
        transposed_cld_frac_r     (i,k_p3) = cld_frac_r     (k,i); // in
        transposed_cld_frac_l     (i,k_p3) = cld_frac_l     (k,i); // in
        transposed_cld_frac_i     (i,k_p3) = cld_frac_i     (k,i); // in
        transposed_qv_prev        (i,k_p3) = qv_prev        (k,i); // in
        transposed_t_prev         (i,k_p3) = t_prev         (k,i); // in
        if (k < 3) transposed_col_location(i,k) = col_location(k,i); // in
      });

      // This has fewer parameters than the Fortran call because they got rid of some in the c++ port of scream
      pam::p3_main_cxx( transposed_qc                .create_ArrayIR() , // inout
                        transposed_nc                .create_ArrayIR() , // inout
                        transposed_qr                .create_ArrayIR() , // inout
                        transposed_nr                .create_ArrayIR() , // inout
                        transposed_theta             .create_ArrayIR() , // inout
                        transposed_qv                .create_ArrayIR() , // inout
                        dt                                             , // inout
                        transposed_qi                .create_ArrayIR() , // inout
                        transposed_qm                .create_ArrayIR() , // inout
                        transposed_ni                .create_ArrayIR() , // inout
                        transposed_bm                .create_ArrayIR() , // inout
                        transposed_pressure_dry      .create_ArrayIR() , // in
                        transposed_dz                .create_ArrayIR() , // in
                        transposed_nc_nuceat_tend    .create_ArrayIR() , // in
                        transposed_nccn_prescribed   .create_ArrayIR() , // in
                        transposed_ni_activated      .create_ArrayIR() , // in
                        transposed_inv_qc_relvar     .create_ArrayIR() , // in
                        it                                             , // in
                        precip_liq_surf              .create_ArrayIR() , //   out
                        precip_ice_surf              .create_ArrayIR() , //   out
                        its                                            , // in
                        ite                                            , // in
                        kts                                            , // in
                        kte                                            , // in
                        transposed_diag_eff_radius_qc.create_ArrayIR() , //   out
                        transposed_diag_eff_radius_qi.create_ArrayIR() , //   out
                        transposed_bulk_qi           .create_ArrayIR() , //   out
                        do_predict_nc                                  , // in
                        do_prescribed_CCN                              , // in
                        transposed_dpres_dry         .create_ArrayIR() , // in
                        transposed_inv_exner         .create_ArrayIR() , // in
                        transposed_qv2qi_depos_tend  .create_ArrayIR() , //   out
                        transposed_precip_liq_flux   .create_ArrayIR() , //   out
                        transposed_precip_ice_flux   .create_ArrayIR() , //   out
                        transposed_cld_frac_r        .create_ArrayIR() , // in
                        transposed_cld_frac_l        .create_ArrayIR() , // in
                        transposed_cld_frac_i        .create_ArrayIR() , // in
                        transposed_liq_ice_exchange  .create_ArrayIR() , //   out
                        transposed_vap_liq_exchange  .create_ArrayIR() , //   out
                        transposed_vap_ice_exchange  .create_ArrayIR() , //   out
                        transposed_qv_prev           .create_ArrayIR() , // in
                        transposed_t_prev            .create_ArrayIR() , // in
                        transposed_col_location      .create_ArrayIR() , // in
                       &elapsed_s                                      );//   out {

      // For inout and out variables, copy transposed data (One kernel for efficiency)
      parallel_for( SimpleBounds<2>(nz+1,ncol) , YAKL_LAMBDA (int k, int i) {
        int k_p3 = (nz+1)-1-k;
        precip_liq_flux   (k,i) = transposed_precip_liq_flux   (i,k_p3); //   out
        precip_ice_flux   (k,i) = transposed_precip_ice_flux   (i,k_p3); //   out
        if (k < nz) {
          k_p3 = nz-1-k;
          qc                (k,i) = transposed_qc                (i,k_p3); // inout
          nc                (k,i) = transposed_nc                (i,k_p3); // inout
          qr                (k,i) = transposed_qr                (i,k_p3); // inout
          nr                (k,i) = transposed_nr                (i,k_p3); // inout
          theta             (k,i) = transposed_theta             (i,k_p3); // inout
          qv                (k,i) = transposed_qv                (i,k_p3); // inout
          qi                (k,i) = transposed_qi                (i,k_p3); // inout
          qm                (k,i) = transposed_qm                (i,k_p3); // inout
          ni                (k,i) = transposed_ni                (i,k_p3); // inout
          bm                (k,i) = transposed_bm                (i,k_p3); // inout
          diag_eff_radius_qc(k,i) = transposed_diag_eff_radius_qc(i,k_p3); //   out
          diag_eff_radius_qi(k,i) = transposed_diag_eff_radius_qi(i,k_p3); //   out
          bulk_qi           (k,i) = transposed_bulk_qi           (i,k_p3); //   out
          qv2qi_depos_tend  (k,i) = transposed_qv2qi_depos_tend  (i,k_p3); //   out
          liq_ice_exchange  (k,i) = transposed_liq_ice_exchange  (i,k_p3); //   out
          vap_liq_exchange  (k,i) = transposed_vap_liq_exchange  (i,k_p3); //   out
          vap_ice_exchange  (k,i) = transposed_vap_ice_exchange  (i,k_p3); //   out
        }
      });

    #else

      it  = 1;
      its = 1;
      ite = ncol;
      kts = 1;
      kte = nz;

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
      auto pressure_dry_host       = pressure_dry      .createHostCopy();
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
      auto dpres_dry_host          = dpres_dry         .createHostCopy();
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

      p3_main_fortran(qc_host.data() , nc_host.data() , qr_host.data() , nr_host.data() , theta_host.data() ,
                      qv_host.data() , dt , qi_host.data() , qm_host.data() , ni_host.data() , bm_host.data() ,
                      pressure_dry_host.data() , dz_host.data() , nc_nuceat_tend_host.data() ,
                      nccn_prescribed_host.data() , ni_activated_host.data() , inv_qc_relvar_host.data() , it ,
                      precip_liq_surf_host.data() , precip_ice_surf_host.data() , its , ite , kts , kte ,
                      diag_eff_radius_qc_host.data() , diag_eff_radius_qi_host.data() , bulk_qi_host.data() ,
                      do_predict_nc , do_prescribed_CCN , dpres_dry_host.data() , inv_exner_host.data() ,
                      qv2qi_depos_tend_host.data() , precip_total_tend_host.data() , nevapr_host.data() ,
                      qr_evap_tend_host.data() , precip_liq_flux_host.data() , precip_ice_flux_host.data() ,
                      cld_frac_r_host.data() , cld_frac_l_host.data() , cld_frac_i_host.data() ,
                      p3_tend_out_host.data() , mu_c_host.data() , lamc_host.data() , liq_ice_exchange_host.data() ,
                      vap_liq_exchange_host.data() , vap_ice_exchange_host.data() , qv_prev_host.data() ,
                      t_prev_host.data() , col_location_host.data() , &elapsed_s );

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
      dpres_dry_host         .deep_copy_to( dpres_dry          );
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
    #endif
    
                    
    ////////////////////////////////////////////////////////////////////////////
    // P3 postprocessing
    ////////////////////////////////////////////////////////////////////////////
    auto liq_ice_exchange_out = dm.get_lev_col<real>("liq_ice_exchange_out");
    auto vap_liq_exchange_out = dm.get_lev_col<real>("vap_liq_exchange_out");
    auto vap_ice_exchange_out = dm.get_lev_col<real>("vap_ice_exchange_out");
    parallel_for( "micro post process" , SimpleBounds<2>(nz,ncol) , YAKL_LAMBDA (int k, int i) {
      // Convert P3 outputs into dynamics coupler state and tracer masses
      rho_c  (k,i) = std::max( qc(k,i)*rho_dry(k,i) , 0._fp );
      rho_nc (k,i) = std::max( nc(k,i)*rho_dry(k,i) , 0._fp );
      rho_r  (k,i) = std::max( qr(k,i)*rho_dry(k,i) , 0._fp );
      rho_nr (k,i) = std::max( nr(k,i)*rho_dry(k,i) , 0._fp );
      rho_i  (k,i) = std::max( qi(k,i)*rho_dry(k,i) , 0._fp );
      rho_ni (k,i) = std::max( ni(k,i)*rho_dry(k,i) , 0._fp );
      rho_m  (k,i) = std::max( qm(k,i)*rho_dry(k,i) , 0._fp );
      rho_bm (k,i) = std::max( bm(k,i)*rho_dry(k,i) , 0._fp );
      rho_v  (k,i) = std::max( qv(k,i)*rho_dry(k,i) , 0._fp );
      // While micro changes total pressure, thus changing exner, the definition
      // of theta depends on the old exner pressure, so we'll use old exner here
      temp   (k,i) = theta(k,i) * exner(k,i);
      // Save qv and temperature for the next call to p3_main
      qv_prev(k,i) = std::max( qv(k,i) , 0._fp );
      t_prev (k,i) = temp(k,i);
      // copy diagnostic quantities to data manager
      liq_ice_exchange_out(k,i) = liq_ice_exchange(k,i);
      vap_liq_exchange_out(k,i) = vap_liq_exchange(k,i);
      vap_ice_exchange_out(k,i) = vap_ice_exchange(k,i);
    });

    // output precipitation rates to be aggregated
    auto precip_liq_surf_out = dm.get<real,3>( "precip_liq_surf_out" );
    auto precip_ice_surf_out = dm.get<real,3>( "precip_ice_surf_out" );
    parallel_for( SimpleBounds<4>(nz,ny,nx,nens) , YAKL_LAMBDA (int k, int j, int i, int iens) {
      int icol = j*nx*nens + i*nens + iens;
      precip_liq_surf_out(j,i,iens) = precip_liq_surf(icol);
      precip_ice_surf_out(j,i,iens) = precip_ice_surf(icol);
    });

    #ifdef PAM_DEBUG
      real mass;
      {
        auto rho_v = dm.get<real,4>( "water_vapor" );
        auto rho_c = dm.get<real,4>( "cloud_water" );
        auto rho_r = dm.get<real,4>( "rain"        );
        auto rho_i = dm.get<real,4>( "ice"         );
        real4d mass4d("mass4d",nz,ny,nx,nens);
        real3d sfc_precip_mass3d("sfc_precip_mass3d",ny,nx,nens);
        parallel_for( SimpleBounds<4>(nz,ny,nx,nens) , YAKL_LAMBDA (int k, int j, int i, int iens) {
          int icol = j*nx*nens + i*nens + iens;
          mass4d(k,j,i,iens) = (rho_v(k,j,i,iens) + rho_c(k,j,i,iens) + rho_r(k,j,i,iens) + rho_i(k,j,i,iens)) *
                               crm_dx * crm_dy * (zint_in(k+1,iens) - zint_in(k,iens));
          sfc_precip_mass3d(j,i,iens) = dt*crm_dx*crm_dy*( precip_liq_surf(icol) + precip_ice_surf(icol) );
        });
        mass = yakl::intrinsics::sum(mass4d) + yakl::intrinsics::sum(sfc_precip_mass3d);
      }
      if ( abs(mass-mass0)/mass0 > 1.e-13 ) std::cout << "WARNING: P3 mass isn't conserved to machine precision\n";
    #endif

    first_step = false;
    etime += dt;
  }



  // Returns saturation vapor pressure
  YAKL_INLINE static real saturation_vapor_pressure(real temp) {
    real tc = temp - 273.15;
    return 610.94 * exp( 17.625*tc / (243.04+tc) );
  }



  YAKL_INLINE static real latent_heat_condensation(real temp) {
    real tc = temp - 273.15;
    return (2500.8 - 2.36*tc + 0.0016*tc*tc - 0.00006*tc*tc*tc)*1000;
  }



  YAKL_INLINE static real cp_moist(real rho_d, real rho_v, real rho_c, real cp_d, real cp_v, real cp_l) {
    // For the moist specific heat, ignore other species than water vapor and cloud droplets
    real rho = rho_d + rho_v + rho_c;
    return rho_d / rho * cp_d  +  rho_v / rho * cp_v  +  rho_c / rho * cp_l;
  }



  // Compute an instantaneous adjustment of sub or super saturation
  YAKL_INLINE static void compute_adjusted_state(real rho, real rho_d , real &rho_v , real &rho_c , real &temp,
                                                 real R_v , real cp_d , real cp_v , real cp_l) {
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
        real rv_loc = std::max( 0._fp , rho_v - rho_cond );          // New vapor density
        real rc_loc = std::max( 0._fp , rho_c + rho_cond );          // New cloud liquid density
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
        real rv_loc = std::max( 0._fp , rho_v + rho_evap );          // New vapor density
        real rc_loc = std::max( 0._fp , rho_c - rho_evap );          // New cloud liquid density
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



  void get_cloud_fraction( realConst2d cld_frac_in , realConst2d qc , realConst2d qr , realConst2d qi ,
                           real2d const &cld_frac_i , real2d const &cld_frac_l , real2d const &cld_frac_r ) {
    using yakl::c::SimpleBounds;

    int nz   = cld_frac_in.dimension[0];
    int ncol = cld_frac_in.dimension[1];

    real constexpr mincld = 0.0001;

    parallel_for( SimpleBounds<2>(nz,ncol) , YAKL_LAMBDA (int k, int i) {
      cld_frac_i(k,i) = std::max(cld_frac_in(k,i), mincld);
      cld_frac_l(k,i) = std::max(cld_frac_in(k,i), mincld);
      cld_frac_r(k,i) = std::max(cld_frac_in(k,i), mincld);
    });

    // precipitation fraction 
    // max overlap is the max cloud fraction in all layers above which are
    // connected to a given layer by a continuous band of precip mass. If
    // there's no precip mass falling into a cell, it's precip frac is equal
    // to the cloud frac, which is probably ~zero.
    // Cycle through the layers from top to bottom and determine if the rain 
    // fraction needs to be updated to match cloud fraction in the layer above.
    parallel_for( ncol , YAKL_LAMBDA (int i) {
      for (int k=nz-2; k >= 0; k--) {
        cld_frac_r(k,i) = std::max( cld_frac_in(k+1,i) , cld_frac_r(k,i) );
      }
    });
  }

  void finalize(pam::PamCoupler &coupler) {
  }

  std::string micro_name() const {
    return "p3";
  }



};



