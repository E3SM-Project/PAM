
#pragma once

//#include "awfl_const.h"
#include "DataManager.h"
#include "pam_coupler.h"
#include "call_p3_from_pam.h"
#include "pam_scream_routines.h"

#define P3_FORTRAN
// #define P3_DEBUG

using pam::PamCoupler;


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


int static constexpr num_tracers_micro = 9;

class Microphysics {
public:
  // Doesn't actually have to be static or constexpr. Could be assigned in the constructor
  int static constexpr num_tracers = 9;

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
  YAKL_INLINE static int get_num_tracers() {
    return num_tracers;
  }



  // Can do whatever you want, but mainly for registering tracers and allocating data
  void init(PamCoupler &coupler) {
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

    coupler.dm.register_and_allocate<real>("qv_prev","qv from prev step"         ,{nz,ny,nx,nens},{"z","y","x","nens"});
    coupler.dm.register_and_allocate<real>("t_prev" ,"Temperature from prev step",{nz,ny,nx,nens},{"z","y","x","nens"});

    auto cloud_water     = coupler.dm.get<real,4>( "cloud_water"     );
    auto cloud_water_num = coupler.dm.get<real,4>( "cloud_water_num" );
    auto rain            = coupler.dm.get<real,4>( "rain"            );
    auto rain_num        = coupler.dm.get<real,4>( "rain_num"        );
    auto ice             = coupler.dm.get<real,4>( "ice"             );
    auto ice_num         = coupler.dm.get<real,4>( "ice_num"         );
    auto ice_rime        = coupler.dm.get<real,4>( "ice_rime"        );
    auto ice_rime_vol    = coupler.dm.get<real,4>( "ice_rime_vol"    );
    auto water_vapor     = coupler.dm.get<real,4>( "water_vapor"     );
    auto qv_prev         = coupler.dm.get<real,4>( "qv_prev"         );
    auto t_prev          = coupler.dm.get<real,4>( "t_prev"          );

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
    bool masterproc = true;
    micro_p3_utils_init_fortran( cp_d , R_d , R_v , rhoh2o , mwh2o , mwdry ,
                                 grav , latvap , latice, cp_l , tmelt , pi , iulog , masterproc );

    std::string dir = "../../../physics/micro/p3";
    std::string ver = "4.1.1";
    int dir_len = dir.length();
    int ver_len = ver.length();
    p3_init_fortran( dir.c_str() , dir_len , ver.c_str() , ver_len );

    coupler.set_option<std::string>("micro","p3");

    etime = 0;
  }



  void timeStep( PamCoupler &coupler , real dt ) {
    if (first_step) {
      if (coupler.get_option<std::string>("sgs") == "shoc") sgs_shoc = true;
    }

    // Get the dimensions sizes
    int nz   = coupler.dm.get_dimension_size("z"   );
    int ny   = coupler.dm.get_dimension_size("y"   );
    int nx   = coupler.dm.get_dimension_size("x"   );
    int nens = coupler.dm.get_dimension_size("nens");
    int ncol = ny*nx*nens;

    auto zint_in = coupler.dm.get<real,2>("vertical_interface_height");

    real crm_dx = coupler.get_xlen() / nx;
    real crm_dy = ny == 1 ? crm_dx : coupler.get_ylen() / ny;

    #ifdef PAM_DEBUG
      real mass0;
      {
        auto rho_v = coupler.dm.get<real,4>( "water_vapor" );
        auto rho_c = coupler.dm.get<real,4>( "cloud_water" );
        auto rho_r = coupler.dm.get<real,4>( "rain"        );
        auto rho_i = coupler.dm.get<real,4>( "ice"         );
        real4d mass4d("mass4d",nz,ny,nx,nens);
        parallel_for( SimpleBounds<4>(nz,ny,nx,nens) , YAKL_LAMBDA (int k, int j, int i, int iens) {
          mass4d(k,j,i,iens) = (rho_v(k,j,i,iens) + rho_c(k,j,i,iens) + rho_r(k,j,i,iens) + rho_i(k,j,i,iens)) *
                               crm_dx * crm_dy * (zint_in(k+1,iens) - zint_in(k,iens));
        });
        mass0 = yakl::intrinsics::sum(mass4d);
      }
    #endif

    // Get tracers dimensioned as (nz,ny*nx*nens)
    auto rho_c  = coupler.dm.get_lev_col<real>("cloud_water"    );
    auto rho_nc = coupler.dm.get_lev_col<real>("cloud_water_num");
    auto rho_r  = coupler.dm.get_lev_col<real>("rain"           );
    auto rho_nr = coupler.dm.get_lev_col<real>("rain_num"       );
    auto rho_i  = coupler.dm.get_lev_col<real>("ice"            );
    auto rho_ni = coupler.dm.get_lev_col<real>("ice_num"        );
    auto rho_m  = coupler.dm.get_lev_col<real>("ice_rime"       );
    auto rho_bm = coupler.dm.get_lev_col<real>("ice_rime_vol"   );
    auto rho_v  = coupler.dm.get_lev_col<real>("water_vapor"    );

    // Get coupler state
    auto rho_dry = coupler.dm.get_lev_col<real>("density_dry");
    auto temp    = coupler.dm.get_lev_col<real>("temp"       );

    // Calculate the grid spacing
    real2d dz("dz",nz,ny*nx*nens);
    parallel_for( "micro dz" , SimpleBounds<4>(nz,ny,nx,nens) ,
                  YAKL_LAMBDA (int k, int j, int i, int iens) {
      dz(k,j*nx*nens + i*nens + iens) = zint_in(k+1,iens) - zint_in(k,iens);
    });

    // Get everything from the DataManager that's not a tracer but is persistent across multiple micro calls
    auto qv_prev = coupler.dm.get_lev_col<real>("qv_prev");
    auto t_prev  = coupler.dm.get_lev_col<real>("t_prev" );

    // Allocates inputs and outputs
    int p3_nout = 49;
    real2d qc                ( "qc"                 ,           nz   , ncol );
    real2d nc                ( "nc"                 ,           nz   , ncol );
    real2d qr                ( "qr"                 ,           nz   , ncol );
    real2d nr                ( "nr"                 ,           nz   , ncol );
    real2d qi                ( "qi"                 ,           nz   , ncol );
    real2d ni                ( "ni"                 ,           nz   , ncol );
    real2d qm                ( "qm"                 ,           nz   , ncol );
    real2d bm                ( "bm"                 ,           nz   , ncol );
    real2d qv                ( "qv"                 ,           nz   , ncol );
    real2d pressure          ( "pressure"           ,           nz   , ncol );
    real2d theta             ( "theta"              ,           nz   , ncol );
    real2d exner             ( "exner"              ,           nz   , ncol );
    real2d inv_exner         ( "inv_exner"          ,           nz   , ncol );
    real2d dpres             ( "dpres"              ,           nz   , ncol );
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
    real3d p3_tend_out       ( "p3_tend_out"        , p3_nout , nz   , ncol );

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
      // Compute total density
      real rho = rho_dry(k,i) + rho_c(k,i) + rho_r(k,i) + rho_i(k,i) + rho_v(k,i);

      // P3 doesn't do saturation adjustment, so we need to do that ahead of time
      // If we're using SHOC, then it does saturation adjustment, so no need to do it here
      if (! sgs_shoc) {
        compute_adjusted_state(rho, rho_dry(k,i) , rho_v(k,i) , rho_c(k,i) , temp(k,i),
                               R_v , cp_d , cp_v , cp_l);
      }

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
      pressure (k,i) = R_d * rho_dry(k,i) * temp(k,i) + R_v * rho_v(k,i) * temp(k,i);
      exner    (k,i) = pow( pressure(k,i) / p0 , R_d / cp_d );
      inv_exner(k,i) = 1. / exner(k,i);
      theta    (k,i) = temp(k,i) / exner(k,i);
      // P3 uses dpres to calculate density via the hydrostatic assumption.
      // So we just reverse this to compute dpres to give true density
      dpres(k,i) = rho_dry(k,i) * grav * dz(k,i);
      // nc_nuceat_tend, nccn_prescribed, and ni_activated are not used
      nc_nuceat_tend (k,i) = 0;
      nccn_prescribed(k,i) = 0;
      ni_activated   (k,i) = 0;
      // col_location is for debugging only, and it will be ignored for now
      if (k < 3) { col_location(k,i) = 1; }

      if (first_step) {
        qv_prev(k,i) = qv  (k,i);
        t_prev (k,i) = temp(k,i);
      }
    });

    if (sgs_shoc) {
      auto ast      = coupler.dm.get_lev_col<real>("cldfrac");
      inv_qc_relvar = coupler.dm.get_lev_col<real>("relvar" );
      get_cloud_fraction( ast , qc , qr , qi , cld_frac_i , cld_frac_l , cld_frac_r );
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
    int its, ite, kts, kte;
    int it = 1;
    bool do_predict_nc = false;
    bool do_prescribed_CCN = false;


    #ifdef P3_DEBUG
    if (etime > 10) {
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

      p3_main_fortran(qc_host.data() , nc_host.data() , qr_host.data() , nr_host.data() , theta_host.data() ,
                      qv_host.data() , dt , qi_host.data() , qm_host.data() , ni_host.data() , bm_host.data() ,
                      pressure_host.data() , dz_host.data() , nc_nuceat_tend_host.data() ,
                      nccn_prescribed_host.data() , ni_activated_host.data() , inv_qc_relvar_host.data() , it ,
                      precip_liq_surf_host.data() , precip_ice_surf_host.data() , its , ite , kts , kte ,
                      diag_eff_radius_qc_host.data() , diag_eff_radius_qi_host.data() , bulk_qi_host.data() ,
                      do_predict_nc , do_prescribed_CCN , dpres_host.data() , inv_exner_host.data() ,
                      qv2qi_depos_tend_host.data() , precip_total_tend_host.data() , nevapr_host.data() ,
                      qr_evap_tend_host.data() , precip_liq_flux_host.data() , precip_ice_flux_host.data() ,
                      cld_frac_r_host.data() , cld_frac_l_host.data() , cld_frac_i_host.data() ,
                      p3_tend_out_host.data() , mu_c_host.data() , lamc_host.data() , liq_ice_exchange_host.data() ,
                      vap_liq_exchange_host.data() , vap_ice_exchange_host.data() , qv_prev_host.data() ,
                      t_prev_host.data() , col_location_host.data() , &elapsed_s );

      its = 0;
      ite = ncol-1;
      kts = 0;
      kte = nz-1;
      pam::call_p3_main_from_pam(dt , it , its , ite , kts , kte , do_predict_nc , do_prescribed_CCN ,
                                 elapsed_s ,
                                 pam::yakl_array_to_arrayIR( qc                 ) ,
                                 pam::yakl_array_to_arrayIR( nc                 ) ,
                                 pam::yakl_array_to_arrayIR( qr                 ) ,
                                 pam::yakl_array_to_arrayIR( nr                 ) ,
                                 pam::yakl_array_to_arrayIR( theta              ) ,
                                 pam::yakl_array_to_arrayIR( qv                 ) ,
                                 pam::yakl_array_to_arrayIR( qi                 ) ,
                                 pam::yakl_array_to_arrayIR( qm                 ) ,
                                 pam::yakl_array_to_arrayIR( ni                 ) ,
                                 pam::yakl_array_to_arrayIR( bm                 ) ,
                                 pam::yakl_array_to_arrayIR( pressure           ) ,
                                 pam::yakl_array_to_arrayIR( dz                 ) ,
                                 pam::yakl_array_to_arrayIR( nc_nuceat_tend     ) ,
                                 pam::yakl_array_to_arrayIR( nccn_prescribed    ) ,
                                 pam::yakl_array_to_arrayIR( ni_activated       ) ,
                                 pam::yakl_array_to_arrayIR( inv_qc_relvar      ) ,
                                 pam::yakl_array_to_arrayIR( precip_liq_surf    ) ,
                                 pam::yakl_array_to_arrayIR( precip_ice_surf    ) ,
                                 pam::yakl_array_to_arrayIR( diag_eff_radius_qc ) ,
                                 pam::yakl_array_to_arrayIR( diag_eff_radius_qi ) ,
                                 pam::yakl_array_to_arrayIR( bulk_qi            ) ,
                                 pam::yakl_array_to_arrayIR( dpres              ) ,
                                 pam::yakl_array_to_arrayIR( inv_exner          ) ,
                                 pam::yakl_array_to_arrayIR( qv2qi_depos_tend   ) ,
                                 pam::yakl_array_to_arrayIR( precip_total_tend  ) ,
                                 pam::yakl_array_to_arrayIR( nevapr             ) ,
                                 pam::yakl_array_to_arrayIR( qr_evap_tend       ) ,
                                 pam::yakl_array_to_arrayIR( precip_liq_flux    ) ,
                                 pam::yakl_array_to_arrayIR( precip_ice_flux    ) ,
                                 pam::yakl_array_to_arrayIR( cld_frac_r         ) ,
                                 pam::yakl_array_to_arrayIR( cld_frac_l         ) ,
                                 pam::yakl_array_to_arrayIR( cld_frac_i         ) ,
                                 pam::yakl_array_to_arrayIR( p3_tend_out        ) ,
                                 pam::yakl_array_to_arrayIR( mu_c               ) ,
                                 pam::yakl_array_to_arrayIR( lamc               ) ,
                                 pam::yakl_array_to_arrayIR( liq_ice_exchange   ) ,
                                 pam::yakl_array_to_arrayIR( vap_liq_exchange   ) ,
                                 pam::yakl_array_to_arrayIR( vap_ice_exchange   ) ,
                                 pam::yakl_array_to_arrayIR( qv_prev            ) ,
                                 pam::yakl_array_to_arrayIR( t_prev             ) ,
                                 pam::yakl_array_to_arrayIR( col_location       ) );

      using yakl::intrinsics::abs;
      using yakl::intrinsics::sum;
      using std::abs;
      #define DEBUG_PRINT(var,var_host) std::cout << std::setw(20) << #var": " << std::setw(20) << sum(abs(var)) << " , " \
                                                                               << std::setw(20) << sum(abs(var_host)) << " , " \
                                                                               << std::setw(20) << abs(sum(abs(var)) - sum(abs(var_host))) / (sum(abs(var_host))+1.e-50) << std::endl
      DEBUG_PRINT(qc                ,qc_host                );
      DEBUG_PRINT(nc                ,nc_host                );
      DEBUG_PRINT(qr                ,qr_host                );
      DEBUG_PRINT(nr                ,nr_host                );
      DEBUG_PRINT(theta             ,theta_host             );
      DEBUG_PRINT(qv                ,qv_host                );
      DEBUG_PRINT(qi                ,qi_host                );
      DEBUG_PRINT(qm                ,qm_host                );
      DEBUG_PRINT(ni                ,ni_host                );
      DEBUG_PRINT(bm                ,bm_host                );
      DEBUG_PRINT(pressure          ,pressure_host          );
      DEBUG_PRINT(dz                ,dz_host                );
      DEBUG_PRINT(nc_nuceat_tend    ,nc_nuceat_tend_host    );
      DEBUG_PRINT(nccn_prescribed   ,nccn_prescribed_host   );
      DEBUG_PRINT(ni_activated      ,ni_activated_host      );
      DEBUG_PRINT(inv_qc_relvar     ,inv_qc_relvar_host     );
      DEBUG_PRINT(precip_liq_surf   ,precip_liq_surf_host   );
      DEBUG_PRINT(precip_ice_surf   ,precip_ice_surf_host   );
      DEBUG_PRINT(diag_eff_radius_qc,diag_eff_radius_qc_host);
      DEBUG_PRINT(diag_eff_radius_qi,diag_eff_radius_qi_host);
      DEBUG_PRINT(bulk_qi           ,bulk_qi_host           );
      DEBUG_PRINT(dpres             ,dpres_host             );
      DEBUG_PRINT(inv_exner         ,inv_exner_host         );
      DEBUG_PRINT(qv2qi_depos_tend  ,qv2qi_depos_tend_host  );
      // DEBUG_PRINT(precip_total_tend ,precip_total_tend_host );
      // DEBUG_PRINT(nevapr            ,nevapr_host            );
      // DEBUG_PRINT(qr_evap_tend      ,qr_evap_tend_host      );
      DEBUG_PRINT(precip_liq_flux   ,precip_liq_flux_host   );
      DEBUG_PRINT(precip_ice_flux   ,precip_ice_flux_host   );
      DEBUG_PRINT(cld_frac_r        ,cld_frac_r_host        );
      DEBUG_PRINT(cld_frac_l        ,cld_frac_l_host        );
      DEBUG_PRINT(cld_frac_i        ,cld_frac_i_host        );
      // DEBUG_PRINT(p3_tend_out       ,p3_tend_out_host       );
      // DEBUG_PRINT(mu_c              ,mu_c_host              );
      // DEBUG_PRINT(lamc              ,lamc_host              );
      DEBUG_PRINT(liq_ice_exchange  ,liq_ice_exchange_host  );
      DEBUG_PRINT(vap_liq_exchange  ,vap_liq_exchange_host  );
      DEBUG_PRINT(vap_ice_exchange  ,vap_ice_exchange_host  );
      DEBUG_PRINT(qv_prev           ,qv_prev_host           );
      DEBUG_PRINT(t_prev            ,t_prev_host            );
      DEBUG_PRINT(col_location      ,col_location_host      );
      abort();
    }
    #endif


    #ifdef P3_FORTRAN

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

      p3_main_fortran(qc_host.data() , nc_host.data() , qr_host.data() , nr_host.data() , theta_host.data() ,
                      qv_host.data() , dt , qi_host.data() , qm_host.data() , ni_host.data() , bm_host.data() ,
                      pressure_host.data() , dz_host.data() , nc_nuceat_tend_host.data() ,
                      nccn_prescribed_host.data() , ni_activated_host.data() , inv_qc_relvar_host.data() , it ,
                      precip_liq_surf_host.data() , precip_ice_surf_host.data() , its , ite , kts , kte ,
                      diag_eff_radius_qc_host.data() , diag_eff_radius_qi_host.data() , bulk_qi_host.data() ,
                      do_predict_nc , do_prescribed_CCN , dpres_host.data() , inv_exner_host.data() ,
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
    
    #else

      its = 0;
      ite = ncol-1;
      kts = 0;
      kte = nz-1;
      pam::call_p3_main_from_pam(dt , it , its , ite , kts , kte , do_predict_nc , do_prescribed_CCN ,
                                 elapsed_s ,
                                 pam::yakl_array_to_arrayIR( qc                 ) ,
                                 pam::yakl_array_to_arrayIR( nc                 ) ,
                                 pam::yakl_array_to_arrayIR( qr                 ) ,
                                 pam::yakl_array_to_arrayIR( nr                 ) ,
                                 pam::yakl_array_to_arrayIR( theta              ) ,
                                 pam::yakl_array_to_arrayIR( qv                 ) ,
                                 pam::yakl_array_to_arrayIR( qi                 ) ,
                                 pam::yakl_array_to_arrayIR( qm                 ) ,
                                 pam::yakl_array_to_arrayIR( ni                 ) ,
                                 pam::yakl_array_to_arrayIR( bm                 ) ,
                                 pam::yakl_array_to_arrayIR( pressure           ) ,
                                 pam::yakl_array_to_arrayIR( dz                 ) ,
                                 pam::yakl_array_to_arrayIR( nc_nuceat_tend     ) ,
                                 pam::yakl_array_to_arrayIR( nccn_prescribed    ) ,
                                 pam::yakl_array_to_arrayIR( ni_activated       ) ,
                                 pam::yakl_array_to_arrayIR( inv_qc_relvar      ) ,
                                 pam::yakl_array_to_arrayIR( precip_liq_surf    ) ,
                                 pam::yakl_array_to_arrayIR( precip_ice_surf    ) ,
                                 pam::yakl_array_to_arrayIR( diag_eff_radius_qc ) ,
                                 pam::yakl_array_to_arrayIR( diag_eff_radius_qi ) ,
                                 pam::yakl_array_to_arrayIR( bulk_qi            ) ,
                                 pam::yakl_array_to_arrayIR( dpres              ) ,
                                 pam::yakl_array_to_arrayIR( inv_exner          ) ,
                                 pam::yakl_array_to_arrayIR( qv2qi_depos_tend   ) ,
                                 pam::yakl_array_to_arrayIR( precip_total_tend  ) ,
                                 pam::yakl_array_to_arrayIR( nevapr             ) ,
                                 pam::yakl_array_to_arrayIR( qr_evap_tend       ) ,
                                 pam::yakl_array_to_arrayIR( precip_liq_flux    ) ,
                                 pam::yakl_array_to_arrayIR( precip_ice_flux    ) ,
                                 pam::yakl_array_to_arrayIR( cld_frac_r         ) ,
                                 pam::yakl_array_to_arrayIR( cld_frac_l         ) ,
                                 pam::yakl_array_to_arrayIR( cld_frac_i         ) ,
                                 pam::yakl_array_to_arrayIR( p3_tend_out        ) ,
                                 pam::yakl_array_to_arrayIR( mu_c               ) ,
                                 pam::yakl_array_to_arrayIR( lamc               ) ,
                                 pam::yakl_array_to_arrayIR( liq_ice_exchange   ) ,
                                 pam::yakl_array_to_arrayIR( vap_liq_exchange   ) ,
                                 pam::yakl_array_to_arrayIR( vap_ice_exchange   ) ,
                                 pam::yakl_array_to_arrayIR( qv_prev            ) ,
                                 pam::yakl_array_to_arrayIR( t_prev             ) ,
                                 pam::yakl_array_to_arrayIR( col_location       ) );

    #endif
                    
    ///////////////////////////////////////////////////////////////////////////////
    // Convert P3 outputs into dynamics coupler state and tracer masses
    ///////////////////////////////////////////////////////////////////////////////
    parallel_for( "micro post process" , SimpleBounds<2>(nz,ncol) , YAKL_LAMBDA (int k, int i) {
      rho_c  (k,i) = yakl::max( qc(k,i)*rho_dry(k,i) , 0._fp );
      rho_nc (k,i) = yakl::max( nc(k,i)*rho_dry(k,i) , 0._fp );
      rho_r  (k,i) = yakl::max( qr(k,i)*rho_dry(k,i) , 0._fp );
      rho_nr (k,i) = yakl::max( nr(k,i)*rho_dry(k,i) , 0._fp );
      rho_i  (k,i) = yakl::max( qi(k,i)*rho_dry(k,i) , 0._fp );
      rho_ni (k,i) = yakl::max( ni(k,i)*rho_dry(k,i) , 0._fp );
      rho_m  (k,i) = yakl::max( qm(k,i)*rho_dry(k,i) , 0._fp );
      rho_bm (k,i) = yakl::max( bm(k,i)*rho_dry(k,i) , 0._fp );
      rho_v  (k,i) = yakl::max( qv(k,i)*rho_dry(k,i) , 0._fp );
      // While micro changes total pressure, thus changing exner, the definition
      // of theta depends on the old exner pressure, so we'll use old exner here
      temp   (k,i) = theta(k,i) * exner(k,i);
      // Save qv and temperature for the next call to p3_main
      qv_prev(k,i) = yakl::max( qv(k,i) , 0._fp );
      t_prev (k,i) = temp(k,i);
    });

    #ifdef PAM_DEBUG
      real mass;
      {
        auto rho_v = coupler.dm.get<real,4>( "water_vapor" );
        auto rho_c = coupler.dm.get<real,4>( "cloud_water" );
        auto rho_r = coupler.dm.get<real,4>( "rain"        );
        auto rho_i = coupler.dm.get<real,4>( "ice"         );
        real4d mass4d("mass4d",nz,ny,nx,nens);
        real3d sfc_precip_mass3d("sfc_precip_mass3d",ny,nx,nens);
        parallel_for( SimpleBounds<4>(nz,ny,nx,nens) , YAKL_LAMBDA (int k, int j, int i, int iens) {
          int icol = j*nx*nens + i*nens + iens;
          mass4d(k,j,i,iens) = (rho_v(k,j,i,iens) + rho_c(k,j,i,iens) + rho_r(k,j,i,iens) + rho_i(k,j,i,iens)) *
                               crm_dx * crm_dy * (zint_in(k+1,iens) - zint_in(k,iens));
          sfc_precip_mass3d(j,i,iens) = dt*crm_dx*crm_dy*( precip_liq_surf(icol)*1000. + precip_ice_surf(icol)*1000. );
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
        real rv_loc = yakl::max( 0._fp , rho_v - rho_cond );          // New vapor density
        real rc_loc = yakl::max( 0._fp , rho_c + rho_cond );          // New cloud liquid density
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
        real rv_loc = yakl::max( 0._fp , rho_v + rho_evap );          // New vapor density
        real rc_loc = yakl::max( 0._fp , rho_c - rho_evap );          // New cloud liquid density
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



  void get_cloud_fraction( realConst2d ast , realConst2d qc , realConst2d qr , realConst2d qi ,
                           real2d const &cld_frac_i , real2d const &cld_frac_l , real2d const &cld_frac_r ) {
    int nz   = ast.dimension[0];
    int ncol = ast.dimension[1];

    real constexpr mincld = 0.0001;
    real constexpr qsmall = 1.e-14;

    parallel_for( Bounds<2>(nz,ncol) , YAKL_LAMBDA (int k, int i) {
      cld_frac_i(k,i) = yakl::max(ast(k,i), mincld);
      cld_frac_l(k,i) = yakl::max(ast(k,i), mincld);
      cld_frac_r(k,i) = yakl::max(ast(k,i), mincld);
    });

    // precipitation fraction 
    // max overlap is the max cloud fraction in all layers above which are
    // connected to this one by a continuous band of precip mass. If
    // there's no precip mass falling into a cell, it's precip frac is equal
    // to the cloud frac, which is probably ~zero.
    // IF rain or ice mix ratios are smaller than threshold,
    // then leave cld_frac_r as cloud fraction at current level
    parallel_for( ncol , YAKL_LAMBDA (int i) {
      for (int k=nz-2; k >= 0; k--) {
        if ( qr(k+1,i) >= qsmall || qi(k+1,i) >= qsmall ) {
          cld_frac_r(k,i) = yakl::max( cld_frac_r(k+1,i) , cld_frac_r(k,i) );
        }
      }
    });
  }



  std::string micro_name() const {
    return "p3";
  }



};



