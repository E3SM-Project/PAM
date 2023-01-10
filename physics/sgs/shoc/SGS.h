
#pragma once

#include "pam_coupler.h"

#define SHOC_CXX

#ifdef SHOC_CXX
 #include "shoc_main_cxx.h"
#endif

extern "C"
void shoc_init_fortran(int &nlev, double &gravit, double &rair, double &rh2o, double &cpair, double &zvir, double &latvap,
                       double &latice, double &karman, double *pref_mid, int &nbot_shoc, int &ntop_shoc);


extern "C"
void shoc_main_fortran(int &shcol, int &nlev, int &nlevi, double &dtime, int &nadv, real *host_dx, double *host_dy,
                       double *thv, double *zt_grid, double *zi_grid, double *pres, double *presi, double *pdel,
                       double *wthl_sfc, double *wqw_sfc, double *uw_sfc, double *vw_sfc, double *wtracer_sfc,
                       int &num_qtracers, double *w_field, double *exner, double *phis, double *host_dse, double *tke,
                       double *thetal, double *qw, double *u_wind, double *v_wind, double *qtracers, double *wthv_sec,
                       double *tkh, double *tk, double *shoc_ql, double *shoc_cldfrac, double *pblh, double *shoc_mix,
                       double *isotropy, double *w_sec, double *thl_sec, double *qw_sec, double *qwthl_sec,
                       double *wthl_sec, double *wqw_sec, double *wtke_sec, double *uw_sec, double *vw_sec, double *w3,
                       double *wqls_sec, double *brunt, double *shoc_ql2 );


class SGS {
public:

  // You should set these in the constructor
  real R_d          ;
  real cp_d         ;
  real cv_d         ;
  real gamma_d      ;
  real kappa_d      ;
  real R_v          ;
  real cp_v         ;
  real cv_v         ;
  real p0           ;
  bool micro_kessler;
  bool micro_p3     ;
  real latvap       ;
  real latice       ;
  real karman       ;

  real grav;
  real cp_l;

  real etime;

  int npbl;

  bool first_step;

  // Indices for all of your tracer quantities
  int static constexpr ID_TKE  = 0;  // Local index for Turbulent Kinetic Energy (m^2/s^2)



  // Set constants and likely num_tracers as well, and anything else you can do immediately
  SGS() {
    R_d           = 287.042;
    cp_d          = 1004.64;
    cv_d          = cp_d - R_d;
    gamma_d       = cp_d / cv_d;
    kappa_d       = R_d  / cp_d;
    R_v           = 461.505;
    cp_v          = 1859;
    cv_v          = R_v - cp_v;
    p0            = 1.e5;
    grav          = 9.80616;
    first_step    = true;
    cp_l          = 4218.;
    micro_kessler = false;
    micro_p3      = false;
    latvap        = 2501000.0;
    latice        = 333700.0;
    karman        = 0.4;
    npbl          = -1;
    etime         = 0;
  }



  // This must return the correct # of tracers **BEFORE** init(...) is called
  static int constexpr get_num_tracers() {
    return 1;
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
    //                 name    description                              positive   adds mass
    coupler.add_tracer("tke" , "Turbulent Kinetic Energy (m^2/s^2)"   , true     , false );

    auto &dm = coupler.get_data_manager_readwrite();

    // Register and allocation non-tracer quantities used by the microphysics
    dm.register_and_allocate<real>( "wthv_sec" , "Buoyancy flux [K m/s]"                , {nz,ny,nx,nens} , {"z","y","x","nens"} );
    dm.register_and_allocate<real>( "tk"       , "Eddy coefficient for momentum [m2/s]" , {nz,ny,nx,nens} , {"z","y","x","nens"} );
    dm.register_and_allocate<real>( "tkh"      , "Eddy coefficent for heat [m2/s]"      , {nz,ny,nx,nens} , {"z","y","x","nens"} );
    dm.register_and_allocate<real>( "cldfrac"  , "Cloud fraction [-]"                   , {nz,ny,nx,nens} , {"z","y","x","nens"} );
    dm.register_and_allocate<real>( "relvar"   , "Relative cloud water variance"        , {nz,ny,nx,nens} , {"z","y","x","nens"} );

    auto tke      = dm.get<real,4>( "tke"      );
    auto wthv_sec = dm.get<real,4>( "wthv_sec" );
    auto tk       = dm.get<real,4>( "tk"       );
    auto tkh      = dm.get<real,4>( "tkh"      );
    auto cldfrac  = dm.get<real,4>( "cldfrac"  );
    auto relvar   = dm.get<real,4>( "relvar"   );

    parallel_for( "sgs zero" , SimpleBounds<4>(nz,ny,nx,nens) , YAKL_LAMBDA (int k, int j, int i, int iens) {
      tke     (k,j,i,iens) = 0;
      wthv_sec(k,j,i,iens) = 0;
      tk      (k,j,i,iens) = 0;
      tkh     (k,j,i,iens) = 0;
      cldfrac (k,j,i,iens) = 0;
      relvar  (k,j,i,iens) = 0;
    });

    coupler.set_option<std::string>("sgs","shoc");
  }



  void timeStep( pam::PamCoupler &coupler , real dt ) {
    using yakl::c::parallel_for;
    using yakl::c::SimpleBounds;

    // Get the dimensions sizes
    int nz   = coupler.get_nz();
    int ny   = coupler.get_ny();
    int nx   = coupler.get_nx();
    int nens = coupler.get_nens();
    int ncol = ny*nx*nens;

    auto &dm = coupler.get_data_manager_readwrite();

    auto pres_mid_tmp = coupler.compute_pressure_array();

    auto pres_mid = pres_mid_tmp.reshape<2>({nz  ,ncol});

    auto zint_pam = dm.get<real,2>("vertical_interface_height");
    auto zmid_pam = dm.get<real,2>("vertical_midpoint_height" );

    real crm_dx = coupler.get_xlen() / nx;
    real crm_dy = ny == 1 ? crm_dx : coupler.get_ylen() / ny;

    // SHOC init requires reference pressure, which we do not have available for the init() call
    if (first_step) {
      // Invert the first column in x, z, and ensemble to use as reference pressure for shoc
      real1d pref_shoc("pref_shoc",nz);
      parallel_for( SimpleBounds<1>(nz) , YAKL_LAMBDA (int k) {
        pref_shoc(k) = pres_mid(nz-1-k,0);
      });
      real zvir = R_v / R_d - 1.;
      int kbot, ktop;

      kbot = nz;
      ktop = 1;
      shoc_init_fortran( nz , grav , R_d , R_v , cp_d , zvir , latvap , latice , karman ,
                         pref_shoc.createHostCopy().data() , kbot , ktop );

      // This check is here instead of init because it's not guaranteed the micro has called init before sgs
      if (! coupler.option_exists("micro")) {
        endrun("ERROR: SHOC requires coupler.set_option<std::string>(\"micro\",...) to be set");
      }
      std::string micro_scheme = coupler.get_option<std::string>("micro");
      if      (micro_scheme == "kessler") { micro_kessler = true; }
      else if (micro_scheme == "p3"     ) { micro_p3      = true; }
      else { endrun("ERROR: SHOC only meant to run with kessler or p3 microphysics"); }
    }

    #ifdef PAM_DEBUG
      real mass0;
      {
        auto rho_d = dm.get<real,4>( "density_dry" );
        auto rho_v = dm.get<real,4>( "water_vapor" );
        real4d rho_c;
        if      (micro_kessler) { rho_c = dm.get<real,4>( "cloud_liquid" ); }
        else if (micro_p3     ) { rho_c = dm.get<real,4>( "cloud_water"  ); }
        real4d mass4d("mass4d",nz,ny,nx,nens);
        parallel_for( SimpleBounds<4>(nz,ny,nx,nens) , YAKL_LAMBDA (int k, int j, int i, int iens) {
          mass4d(k,j,i,iens) = ( rho_d(k,j,i,iens) + rho_v(k,j,i,iens) + rho_c(k,j,i,iens) ) *
                               crm_dx * crm_dy * (zint_pam(k+1,iens) - zint_pam(k,iens));
        });
        mass0 = yakl::intrinsics::sum(mass4d);
      }
    #endif

    // Get saved SHOC-related variables
    auto tke      = dm.get_lev_col<real>( "tke"      ); // PAM Tracer                 ; don't compute
    auto wthv_sec = dm.get_lev_col<real>( "wthv_sec" ); // Reuse from last SHOC output; don't compute
    auto tk       = dm.get_lev_col<real>( "tk"       ); // Reuse from last SHOC output; don't compute
    auto tkh      = dm.get_lev_col<real>( "tkh"      ); // Reuse from last SHOC output; don't compute
    auto cldfrac  = dm.get_lev_col<real>( "cldfrac"  ); // Reuse from last SHOC output; don't compute
    auto relvar   = dm.get_lev_col<real>( "relvar"   ); // Computed on output for P3
    // Get coupler state
    auto rho_d        = dm.get_lev_col<real>( "density_dry"  );
    auto uvel         = dm.get_lev_col<real>( "uvel"         );
    auto vvel         = dm.get_lev_col<real>( "vvel"         );
    auto wvel         = dm.get_lev_col<real>( "wvel"         );
    auto temp         = dm.get_lev_col<real>( "temp"         );
    auto rho_v        = dm.get_lev_col<real>( "water_vapor"  ); // Water vapor mass

    // TODO: What do we do about rho_dry and wvel ???

    // Grab cloud liquid tracer, and determine what other tracers to diffuse in SHOC
    // TODO: If we add new MMF modules that add other tracers that need to be diffused, then this situation
    //       must be handled
    pam::MultiField<real,2> qtracers_pam;  // Extra tracers for SHOC to diffuse
    real2d rho_c;                     // Cloud liquid mass (name differs for different micro schemes)
    if        (micro_kessler) {
      rho_c = dm.get_lev_col<real>( "cloud_liquid" );
      qtracers_pam.add_field( dm.get_lev_col<real>("precip_liquid") );
    } else if (micro_p3     ) {
      rho_c = dm.get_lev_col<real>( "cloud_water"  );
      qtracers_pam.add_field( dm.get_lev_col<real>("cloud_water_num") );
      qtracers_pam.add_field( dm.get_lev_col<real>("rain"           ) );
      qtracers_pam.add_field( dm.get_lev_col<real>("rain_num"       ) );
      qtracers_pam.add_field( dm.get_lev_col<real>("ice"            ) );
      qtracers_pam.add_field( dm.get_lev_col<real>("ice_num"        ) );
      qtracers_pam.add_field( dm.get_lev_col<real>("ice_rime"       ) );
      qtracers_pam.add_field( dm.get_lev_col<real>("ice_rime_vol"   ) );
    }

    int num_qtracers = qtracers_pam.get_num_fields();

    real4d zint_tmp("zint_tmp",nz+1,ny,nx,nens);
    real4d zmid_tmp("zint_tmp",nz  ,ny,nx,nens);

    parallel_for( SimpleBounds<4>(nz+1,ny,nx,nens) , YAKL_LAMBDA (int k, int j, int i, int iens) {
      zint_tmp(k,j,i,iens) = zint_pam(k,iens);
      if (k < nz) zmid_tmp(k,j,i,iens) = zmid_pam(k,iens);
    });

    auto zint = zint_tmp.reshape<2>({nz+1,ncol});
    auto zmid = zmid_tmp.reshape<2>({nz  ,ncol});

    // Create variables for SHOC call (these are all inverted for vertical indices
    real1d shoc_host_dx    ("shoc_host_dx"    ,                  ncol); // grid spacing of host model in x direction [m]
    real1d shoc_host_dy    ("shoc_host_dy"    ,                  ncol); // grid spacing of host model in y direction [m]
    real2d shoc_zt_grid    ("shoc_zt_grid"    ,             nz  ,ncol); // heights, for thermo grid [m]
    real2d shoc_zi_grid    ("shoc_zi_grid"    ,             nz+1,ncol); // heights, for interface grid [m]
    real2d shoc_pres       ("shoc_pres"       ,             nz  ,ncol); // pressure levels on thermo grid [Pa]
    real2d shoc_presi      ("shoc_presi"      ,             nz+1,ncol); // pressure levels on interface grid [Pa]
    real2d shoc_pdel       ("shoc_pdel"       ,             nz  ,ncol); // Differences in pressure levels [Pa]
    real2d shoc_thv        ("shoc_thv"        ,             nz  ,ncol); // virtual potential temperature [K]
    real2d shoc_w_field    ("shoc_w_field"    ,             nz  ,ncol); // large scale vertical velocity [m/s]
    real1d shoc_wthl_sfc   ("shoc_wthl_sfc"   ,                  ncol); // Surface sensible heat flux [K m/s]
    real1d shoc_wqw_sfc    ("shoc_wqw_sfc"    ,                  ncol); // Surface latent heat flux [kg/kg m/s]
    real1d shoc_uw_sfc     ("shoc_uw_sfc"     ,                  ncol); // Surface momentum flux (u-direction) [m2/s2]
    real1d shoc_vw_sfc     ("shoc_vw_sfc"     ,                  ncol); // Surface momentum flux (v-direction) [m2/s2]
    real2d shoc_wtracer_sfc("shoc_wtracer_sfc",num_qtracers,     ncol); // Surface flux for tracers [varies]
    real2d shoc_exner      ("shoc_exner"      ,             nz  ,ncol); // Exner function [-]
    real2d shoc_inv_exner  ("shoc_inv_exner"  ,             nz  ,ncol); // 1/Exner [-]
    real1d shoc_phis       ("shoc_phis"       ,                  ncol); // Host model surface geopotential height
    real2d shoc_host_dse   ("shoc_host_dse"   ,             nz  ,ncol); // dry static energy [J/kg];  dse = Cp*T + g*z + phis
    real2d shoc_tke        ("shoc_tke"        ,             nz  ,ncol); // turbulent kinetic energy [m2/s2]
    real2d shoc_thetal     ("shoc_thetal"     ,             nz  ,ncol); // liquid water potential temperature [K]
    real2d shoc_qw         ("shoc_qw"         ,             nz  ,ncol); // total water mixing ratio [kg/kg]
    real2d shoc_u_wind     ("shoc_u_wind"     ,             nz  ,ncol); // u wind component [m/s]
    real2d shoc_v_wind     ("shoc_v_wind"     ,             nz  ,ncol); // v wind component [m/s]
    real2d shoc_wthv_sec   ("shoc_wthv_sec"   ,             nz  ,ncol); // buoyancy flux [K m/s]
    real3d shoc_qtracers   ("shoc_qtracers"   ,num_qtracers,nz  ,ncol); // tracers [varies]
    real2d shoc_tk         ("shoc_tk"         ,             nz  ,ncol); // eddy coefficient for momentum [m2/s]
    real2d shoc_tkh        ("shoc_tkh"        ,             nz  ,ncol); // eddy coefficent for heat [m2/s]
    real2d shoc_cldfrac    ("shoc_cldfrac"    ,             nz  ,ncol); // Cloud fraction [-]
    real2d shoc_ql         ("shoc_ql"         ,             nz  ,ncol); // cloud liquid mixing ratio [kg/kg]
    real1d shoc_pblh       ("shoc_pblh"       ,                  ncol); // OUT: planetary boundary layer depth [m]
    real2d shoc_ql2        ("shoc_ql2"        ,             nz  ,ncol); // OUT: cloud liquid mixing ratio variance [kg^2/kg^2]
    real2d shoc_mix        ("shoc_mix"        ,             nz  ,ncol); // OUT: Turbulent length scale [m]
    real2d shoc_w_sec      ("shoc_w_sec"      ,             nz  ,ncol); // OUT: vertical velocity variance [m2/s2]
    real2d shoc_thl_sec    ("shoc_thl_sec"    ,             nz+1,ncol); // OUT: temperature variance [K^2]
    real2d shoc_qw_sec     ("shoc_qw_sec"     ,             nz+1,ncol); // OUT: moisture variance [kg2/kg2]
    real2d shoc_qwthl_sec  ("shoc_qwthl_sec"  ,             nz+1,ncol); // OUT: temp moisture covariance [K kg/kg]
    real2d shoc_wthl_sec   ("shoc_wthl_sec"   ,             nz+1,ncol); // OUT: vertical heat flux [K m/s]
    real2d shoc_wqw_sec    ("shoc_wqw_sec"    ,             nz+1,ncol); // OUT: vertical moisture flux [K m/s]
    real2d shoc_wtke_sec   ("shoc_wtke_sec"   ,             nz+1,ncol); // OUT: vertical tke flux [m3/s3]
    real2d shoc_uw_sec     ("shoc_uw_sec"     ,             nz+1,ncol); // OUT: vertical zonal momentum flux [m2/s2]
    real2d shoc_vw_sec     ("shoc_vw_sec"     ,             nz+1,ncol); // OUT: vertical meridional momentum flux [m2/s2]
    real2d shoc_w3         ("shoc_w3"         ,             nz+1,ncol); // OUT: third moment vertical velocity [m3/s3]
    real2d shoc_wqls_sec   ("shoc_wqls_sec"   ,             nz  ,ncol); // OUT: liquid water flux [kg/kg m/s]
    real2d shoc_brunt      ("shoc_brunt"      ,             nz  ,ncol); // OUT: brunt vaisala frequency [s-1]
    real2d shoc_isotropy   ("shoc_isotropy"   ,             nz  ,ncol); // OUT: return to isotropic timescale [s]

    real p0     = this->p0    ;
    real grav   = this->grav  ;
    real R_d    = this->R_d   ;
    real R_v    = this->R_v   ;
    real cp_d   = this->cp_d  ;
    real latvap = this->latvap;

    // Compute inputs for SHOC (reordering the vertical dimension)
    parallel_for( SimpleBounds<2>(nz+1,ncol) , YAKL_LAMBDA (int k, int i) {
      if (k == 0) {
        shoc_host_dx    (i) = crm_dx;
        shoc_host_dy    (i) = crm_dy;
        shoc_wthl_sfc   (i) = 0;
        shoc_wqw_sfc    (i) = 0;
        shoc_uw_sfc     (i) = 0;
        shoc_vw_sfc     (i) = 0;
        shoc_phis       (i) = zint(0,i) * grav;
        for (int tr = 0; tr < num_qtracers; tr++) {
          shoc_wtracer_sfc(tr,i) = 0;
        }
      }
      if (k < nz) {
        int k_shoc = nz-1-k;
        real z       = zmid    (k,i);
        real dz      = zint(k+1,i) - zint(k,i);
        real press   = pres_mid(k,i);
        real t       = temp    (k,i);
        real qv      = rho_v   (k,i) / rho_d(k,i);
        real ql      = rho_c   (k,i) / rho_d(k,i);
        real exner   = pow( press / p0 , R_d / cp_d );
        real theta   = t / exner;
        // https://glossary.ametsoc.org/wiki/Virtual_potential_temperature
        real theta_v = theta * (1 + 0.61_fp * qv);
        // https://glossary.ametsoc.org/wiki/Liquid_water_potential_temperature
        // According to update_host_dse, the simplified version is used here
        real theta_l = theta - (latvap/cp_d) * ql;
        // dry static energy = Cp*T + g*z + phis
        real dse     = cp_d * t + grav * z;
        shoc_ql       (k_shoc,i) = ql;
        shoc_qw       (k_shoc,i) = qv + ql;
        shoc_zt_grid  (k_shoc,i) = z;
        shoc_pres     (k_shoc,i) = press;
        shoc_pdel     (k_shoc,i) = grav * rho_d(k,i) * dz;
        shoc_thv      (k_shoc,i) = theta_v;
        shoc_w_field  (k_shoc,i) = wvel(k,i);
        shoc_exner    (k_shoc,i) = exner;
        shoc_inv_exner(k_shoc,i) = 1._fp / exner;
        shoc_host_dse (k_shoc,i) = dse;
        // TKE is a tracer, so it's stored as mass-weighted. SHOC doesn't want mass-weighted
        shoc_tke      (k_shoc,i) = tke(k,i) / rho_d(k,i);  
        shoc_thetal   (k_shoc,i) = theta_l;
        shoc_u_wind   (k_shoc,i) = uvel(k,i);
        shoc_v_wind   (k_shoc,i) = vvel(k,i);
        shoc_wthv_sec (k_shoc,i) = wthv_sec(k,i);
        for (int tr=0; tr < num_qtracers; tr++) {
          shoc_qtracers(tr,k_shoc,i) = qtracers_pam(tr,k,i) / rho_d(k,i);
        }
        shoc_tk     (k_shoc,i) = tk     (k,i);
        shoc_tkh    (k_shoc,i) = tkh    (k,i);
        shoc_cldfrac(k_shoc,i) = cldfrac(k,i);
      }
      int k_shoc = nz-k;
      shoc_zi_grid(k_shoc,i) = zint    (k,i);
      real pres_int;
      if      (k == 0 ) {
        pres_int = pres_mid(k  ,i) + grav*rho_d(k  ,i)*(zint(k+1,i)-zint(k  ,i))/2;
      } else if (k == nz) {
        pres_int = pres_mid(k-1,i) - grav*rho_d(k-1,i)*(zint(k  ,i)-zint(k-1,i))/2;
      } else {
        pres_int = 0.5_fp * ( pres_mid(k-1,i) - grav*rho_d(k-1,i)*(zint(k  ,i)-zint(k-1,i))/2 +
                              pres_mid(k  ,i) + grav*rho_d(k  ,i)*(zint(k+1,i)-zint(k  ,i))/2 ); }
      shoc_presi  (k_shoc,i) = pres_int;
    });

    int nadv = 1;

#if defined(SHOC_FORTRAN)
    auto shoc_host_dx_host     = shoc_host_dx    .createHostCopy();
    auto shoc_host_dy_host     = shoc_host_dy    .createHostCopy();
    auto shoc_zt_grid_host     = shoc_zt_grid    .createHostCopy();
    auto shoc_zi_grid_host     = shoc_zi_grid    .createHostCopy();
    auto shoc_pres_host        = shoc_pres       .createHostCopy();
    auto shoc_presi_host       = shoc_presi      .createHostCopy();
    auto shoc_pdel_host        = shoc_pdel       .createHostCopy();
    auto shoc_thv_host         = shoc_thv        .createHostCopy();
    auto shoc_w_field_host     = shoc_w_field    .createHostCopy();
    auto shoc_wthl_sfc_host    = shoc_wthl_sfc   .createHostCopy();
    auto shoc_wqw_sfc_host     = shoc_wqw_sfc    .createHostCopy();
    auto shoc_uw_sfc_host      = shoc_uw_sfc     .createHostCopy();
    auto shoc_vw_sfc_host      = shoc_vw_sfc     .createHostCopy();
    auto shoc_wtracer_sfc_host = shoc_wtracer_sfc.createHostCopy();
    auto shoc_exner_host       = shoc_exner      .createHostCopy();
    auto shoc_inv_exner_host   = shoc_inv_exner  .createHostCopy();
    auto shoc_phis_host        = shoc_phis       .createHostCopy();
    auto shoc_host_dse_host    = shoc_host_dse   .createHostCopy();
    auto shoc_tke_host         = shoc_tke        .createHostCopy();
    auto shoc_thetal_host      = shoc_thetal     .createHostCopy();
    auto shoc_qw_host          = shoc_qw         .createHostCopy();
    auto shoc_u_wind_host      = shoc_u_wind     .createHostCopy();
    auto shoc_v_wind_host      = shoc_v_wind     .createHostCopy();
    auto shoc_wthv_sec_host    = shoc_wthv_sec   .createHostCopy();
    auto shoc_qtracers_host    = shoc_qtracers   .createHostCopy();
    auto shoc_tk_host          = shoc_tk         .createHostCopy();
    auto shoc_tkh_host         = shoc_tkh        .createHostCopy();
    auto shoc_cldfrac_host     = shoc_cldfrac    .createHostCopy();
    auto shoc_ql_host          = shoc_ql         .createHostCopy();
    auto shoc_pblh_host        = shoc_pblh       .createHostCopy();
    auto shoc_ql2_host         = shoc_ql2        .createHostCopy();
    auto shoc_mix_host         = shoc_mix        .createHostCopy();
    auto shoc_w_sec_host       = shoc_w_sec      .createHostCopy();
    auto shoc_thl_sec_host     = shoc_thl_sec    .createHostCopy();
    auto shoc_qw_sec_host      = shoc_qw_sec     .createHostCopy();
    auto shoc_qwthl_sec_host   = shoc_qwthl_sec  .createHostCopy();
    auto shoc_wthl_sec_host    = shoc_wthl_sec   .createHostCopy();
    auto shoc_wqw_sec_host     = shoc_wqw_sec    .createHostCopy();
    auto shoc_wtke_sec_host    = shoc_wtke_sec   .createHostCopy();
    auto shoc_uw_sec_host      = shoc_uw_sec     .createHostCopy();
    auto shoc_vw_sec_host      = shoc_vw_sec     .createHostCopy();
    auto shoc_w3_host          = shoc_w3         .createHostCopy();
    auto shoc_wqls_sec_host    = shoc_wqls_sec   .createHostCopy();
    auto shoc_brunt_host       = shoc_brunt      .createHostCopy();
    auto shoc_isotropy_host    = shoc_isotropy   .createHostCopy();

    int nzp1 = nz+1;
    // IMPORTANT: SHOC appears to actually want 1/exner for the exner parameter
    shoc_main_fortran( ncol, nz, nzp1, dt, nadv, 
                       shoc_host_dx_host.data(), shoc_host_dy_host.data(), shoc_thv_host.data(), 
                       shoc_zt_grid_host.data(), shoc_zi_grid_host.data(), shoc_pres_host.data(), shoc_presi_host.data(),
                       shoc_pdel_host.data(),
                       shoc_wthl_sfc_host.data(), shoc_wqw_sfc_host.data(), shoc_uw_sfc_host.data(), shoc_vw_sfc_host.data(), 
                       shoc_wtracer_sfc_host.data(), num_qtracers, shoc_w_field_host.data(), 
                       shoc_inv_exner_host.data(), shoc_phis_host.data(), 
                       shoc_host_dse_host.data(), shoc_tke_host.data(), shoc_thetal_host.data(), shoc_qw_host.data(), 
                       shoc_u_wind_host.data(), shoc_v_wind_host.data(), shoc_qtracers_host.data(),
                       shoc_wthv_sec_host.data(), shoc_tkh_host.data(), shoc_tk_host.data(),
                       shoc_ql_host.data(), shoc_cldfrac_host.data(),
                       shoc_pblh_host.data(),
                       shoc_mix_host.data(), shoc_isotropy_host.data(),
                       shoc_w_sec_host.data(), shoc_thl_sec_host.data(), shoc_qw_sec_host.data(), shoc_qwthl_sec_host.data(),
                       shoc_wthl_sec_host.data(), shoc_wqw_sec_host.data(), shoc_wtke_sec_host.data(),
                       shoc_uw_sec_host.data(), shoc_vw_sec_host.data(), shoc_w3_host.data(),
                       shoc_wqls_sec_host.data(), shoc_brunt_host.data(), shoc_ql2_host.data() );

    shoc_host_dx_host    .deep_copy_to(shoc_host_dx    );
    shoc_host_dy_host    .deep_copy_to(shoc_host_dy    );
    shoc_zt_grid_host    .deep_copy_to(shoc_zt_grid    );
    shoc_zi_grid_host    .deep_copy_to(shoc_zi_grid    );
    shoc_pres_host       .deep_copy_to(shoc_pres       );
    shoc_presi_host      .deep_copy_to(shoc_presi      );
    shoc_pdel_host       .deep_copy_to(shoc_pdel       );
    shoc_thv_host        .deep_copy_to(shoc_thv        );
    shoc_w_field_host    .deep_copy_to(shoc_w_field    );
    shoc_wthl_sfc_host   .deep_copy_to(shoc_wthl_sfc   );
    shoc_wqw_sfc_host    .deep_copy_to(shoc_wqw_sfc    );
    shoc_uw_sfc_host     .deep_copy_to(shoc_uw_sfc     );
    shoc_vw_sfc_host     .deep_copy_to(shoc_vw_sfc     );
    shoc_wtracer_sfc_host.deep_copy_to(shoc_wtracer_sfc);
    shoc_exner_host      .deep_copy_to(shoc_exner      );
    shoc_inv_exner_host  .deep_copy_to(shoc_inv_exner  );
    shoc_phis_host       .deep_copy_to(shoc_phis       );
    shoc_host_dse_host   .deep_copy_to(shoc_host_dse   );
    shoc_tke_host        .deep_copy_to(shoc_tke        );
    shoc_thetal_host     .deep_copy_to(shoc_thetal     );
    shoc_qw_host         .deep_copy_to(shoc_qw         );
    shoc_u_wind_host     .deep_copy_to(shoc_u_wind     );
    shoc_v_wind_host     .deep_copy_to(shoc_v_wind     );
    shoc_wthv_sec_host   .deep_copy_to(shoc_wthv_sec   );
    shoc_qtracers_host   .deep_copy_to(shoc_qtracers   );
    shoc_tk_host         .deep_copy_to(shoc_tk         );
    shoc_tkh_host        .deep_copy_to(shoc_tkh        );
    shoc_cldfrac_host    .deep_copy_to(shoc_cldfrac    );
    shoc_ql_host         .deep_copy_to(shoc_ql         );
    shoc_pblh_host       .deep_copy_to(shoc_pblh       );
    shoc_ql2_host        .deep_copy_to(shoc_ql2        );
    shoc_mix_host        .deep_copy_to(shoc_mix        );
    shoc_w_sec_host      .deep_copy_to(shoc_w_sec      );
    shoc_thl_sec_host    .deep_copy_to(shoc_thl_sec    );
    shoc_qw_sec_host     .deep_copy_to(shoc_qw_sec     );
    shoc_qwthl_sec_host  .deep_copy_to(shoc_qwthl_sec  );
    shoc_wthl_sec_host   .deep_copy_to(shoc_wthl_sec   );
    shoc_wqw_sec_host    .deep_copy_to(shoc_wqw_sec    );
    shoc_wtke_sec_host   .deep_copy_to(shoc_wtke_sec   );
    shoc_uw_sec_host     .deep_copy_to(shoc_uw_sec     );
    shoc_vw_sec_host     .deep_copy_to(shoc_vw_sec     );
    shoc_w3_host         .deep_copy_to(shoc_w3         );
    shoc_wqls_sec_host   .deep_copy_to(shoc_wqls_sec   );
    shoc_brunt_host      .deep_copy_to(shoc_brunt      );
    shoc_isotropy_host   .deep_copy_to(shoc_isotropy   );
#elif defined(SHOC_CXX)
    auto shoc_host_dx_ir = MakeIRType<real, 1>(shoc_host_dx.data(), {ncol}, shoc_host_dx.get_memory_space(), "shoc_host_dx");
    auto shoc_host_dy_ir = MakeIRType<real, 1>(shoc_host_dy.data(), {ncol}, shoc_host_dy.get_memory_space(), "shoc_host_dy");
    auto shoc_zt_grid_ir = MakeIRType<real, 2>(shoc_zt_grid.data(), {nz, ncol}, shoc_zt_grid.get_memory_space(), "shoc_zt_grid");
    auto shoc_zi_grid_ir = MakeIRType<real, 2>(shoc_zi_grid.data(), {nz+1, ncol}, shoc_zi_grid.get_memory_space(), "shoc_zi_grid");
    auto shoc_pres_ir    = MakeIRType<real, 2>(shoc_pres.data(),    {nz, ncol}, shoc_pres.get_memory_space(), "shoc_pres");
    auto shoc_presi_ir   = MakeIRType<real, 2>(shoc_presi.data(),   {nz+1, ncol}, shoc_presi.get_memory_space(), "shoc_presi");
    auto shoc_pdel_ir    = MakeIRType<real, 2>(shoc_pdel.data(),    {nz, ncol}, shoc_pdel.get_memory_space(), "shoc_pdel");
    auto shoc_thv_ir     = MakeIRType<real, 2>(shoc_thv.data(),     {nz, ncol}, shoc_thv.get_memory_space(), "shoc_thv");
    auto shoc_w_field_ir = MakeIRType<real, 2>(shoc_w_field.data(), {nz, ncol}, shoc_w_field.get_memory_space(), "shoc_w_field");
    auto shoc_wthl_sfc_ir= MakeIRType<real, 1>(shoc_wthl_sfc.data(),{ncol}, shoc_wthl_sfc.get_memory_space(), "shoc_wthl_sfc");
    auto shoc_wqw_sfc_ir = MakeIRType<real, 1>(shoc_wqw_sfc.data(), {ncol}, shoc_wqw_sfc.get_memory_space(), "shoc_wqw_sfc");
    auto shoc_uw_sfc_ir  = MakeIRType<real, 1>(shoc_uw_sfc.data(),  {ncol}, shoc_uw_sfc.get_memory_space(), "shoc_uw_sfc");
    auto shoc_vw_sfc_ir  = MakeIRType<real, 1>(shoc_vw_sfc.data(),  {ncol}, shoc_vw_sfc.get_memory_space(), "shoc_vw_sfc");

    auto shoc_wtracer_sfc_ir = MakeIRType<real, 2>(shoc_wtracer_sfc.data(), {num_qtracers, ncol}, shoc_wtracer_sfc.get_memory_space(), "shoc_wtracer_sfc");
    auto shoc_exner_ir       = MakeIRType<real, 2>(shoc_exner.data(), {nz, ncol}, shoc_exner.get_memory_space(), "shoc_exner");
    auto shoc_inv_exner_ir   = MakeIRType<real, 2>(shoc_inv_exner.data(), {nz, ncol}, shoc_inv_exner.get_memory_space(), "shoc_inv_exner");
    auto shoc_phis_ir        = MakeIRType<real, 1>(shoc_phis.data(), {ncol}, shoc_phis.get_memory_space(), "shoc_phis");
    auto shoc_host_dse_ir    = MakeIRType<real, 2>(shoc_host_dse.data(), {nz, ncol}, shoc_host_dse.get_memory_space(), "shoc_host_dse");
    auto shoc_tke_ir         = MakeIRType<real, 2>(shoc_tke.data(), {nz, ncol}, shoc_tke.get_memory_space(), "shoc_tke");
    auto shoc_thetal_ir      = MakeIRType<real, 2>(shoc_thetal.data(), {nz, ncol}, shoc_thetal.get_memory_space(), "shoc_thetal");
    auto shoc_qw_ir          = MakeIRType<real, 2>(shoc_qw.data(), {nz, ncol}, shoc_qw.get_memory_space(), "shoc_qw");
    auto shoc_u_wind_ir      = MakeIRType<real, 2>(shoc_u_wind.data(), {nz, ncol}, shoc_u_wind.get_memory_space(), "shoc_u_wind");
    auto shoc_v_wind_ir      = MakeIRType<real, 2>(shoc_v_wind.data(), {nz, ncol}, shoc_v_wind.get_memory_space(), "shoc_v_wind");
    auto shoc_wthv_sec_ir    = MakeIRType<real, 2>(shoc_wthv_sec.data(), {nz, ncol}, shoc_wthv_sec.get_memory_space(), "shoc_wthv_sec");
    auto shoc_qtracers_ir    = MakeIRType<real, 3>(shoc_qtracers.data(), {num_qtracers, nz, ncol}, shoc_qtracers.get_memory_space(), "shoc_qtracers");
    auto shoc_tk_ir          = MakeIRType<real, 2>(shoc_tk.data(), {nz, ncol}, shoc_tk.get_memory_space(), "shoc_tk");
    auto shoc_tkh_ir         = MakeIRType<real, 2>(shoc_tkh.data(), {nz, ncol}, shoc_tkh.get_memory_space(), "shoc_tkh");
    auto shoc_cldfrac_ir     = MakeIRType<real, 2>(shoc_cldfrac.data(), {nz, ncol}, shoc_cldfrac.get_memory_space(), "shoc_cldfrac");
    auto shoc_ql_ir          = MakeIRType<real, 2>(shoc_ql.data(), {nz, ncol}, shoc_ql.get_memory_space(), "shoc_ql");
    auto shoc_pblh_ir        = MakeIRType<real, 1>(shoc_pblh.data(), {ncol}, shoc_pblh.get_memory_space(), "shoc_pblh");
    auto shoc_ql2_ir         = MakeIRType<real, 2>(shoc_ql2.data(), {nz, ncol}, shoc_ql2.get_memory_space(), "shoc_ql2");
    auto shoc_mix_ir         = MakeIRType<real, 2>(shoc_mix.data(), {nz, ncol}, shoc_mix.get_memory_space(), "shoc_mix");
    auto shoc_w_sec_ir       = MakeIRType<real, 2>(shoc_w_sec.data(), {nz, ncol}, shoc_w_sec.get_memory_space(), "shoc_w_sec");
    auto shoc_thl_sec_ir     = MakeIRType<real, 2>(shoc_thl_sec.data(), {nz+1, ncol}, shoc_thl_sec.get_memory_space(), "shoc_thl_sec");
    auto shoc_qw_sec_ir      = MakeIRType<real, 2>(shoc_qw_sec.data(), {nz+1, ncol}, shoc_qw_sec.get_memory_space(), "shoc_qw_sec");
    auto shoc_qwthl_sec_ir   = MakeIRType<real, 2>(shoc_qwthl_sec.data(), {nz+1, ncol}, shoc_qwthl_sec.get_memory_space(), "shoc_qwthl_sec");
    auto shoc_wthl_sec_ir    = MakeIRType<real, 2>(shoc_wthl_sec.data(), {nz+1, ncol}, shoc_wthl_sec.get_memory_space(), "shoc_wthl_sec");
    auto shoc_wqw_sec_ir     = MakeIRType<real, 2>(shoc_wqw_sec.data(), {nz+1, ncol}, shoc_wqw_sec.get_memory_space(), "shoc_wqw_sec");
    auto shoc_wtke_sec_ir    = MakeIRType<real, 2>(shoc_wtke_sec.data(), {nz+1, ncol}, shoc_wtke_sec.get_memory_space(), "shoc_wtke_sec");
    auto shoc_uw_sec_ir      = MakeIRType<real, 2>(shoc_uw_sec.data(), {nz+1, ncol}, shoc_uw_sec.get_memory_space(), "shoc_uw_sec");
    auto shoc_vw_sec_ir      = MakeIRType<real, 2>(shoc_vw_sec.data(), {nz+1, ncol}, shoc_vw_sec.get_memory_space(), "shoc_vw_sec");
    auto shoc_w3_ir          = MakeIRType<real, 2>(shoc_w3.data(), {nz+1, ncol}, shoc_w3.get_memory_space(), "shoc_w3");
    auto shoc_wqls_sec_ir    = MakeIRType<real, 2>(shoc_wqls_sec.data(), {nz, ncol}, shoc_wqls_sec.get_memory_space(), "shoc_wqls_sec");
    auto shoc_brunt_ir       = MakeIRType<real, 2>(shoc_brunt.data(), {nz, ncol}, shoc_brunt.get_memory_space(),  "shoc_brunt");
    auto shoc_isotropy_ir    = MakeIRType<real, 2>(shoc_isotropy.data(), {nz, ncol}, shoc_isotropy.get_memory_space(), "shoc_isotropy");

    int nzp1 = nz+1;
    // IMPORTANT: SHOC appears to actually want 1/exner for the exner parameter
    shoc_main_cxx( ncol, nz, nzp1, dt, nadv,
                   shoc_host_dx_ir, shoc_host_dy_ir, shoc_thv_ir,
                   shoc_zt_grid_ir, shoc_zi_grid_ir, shoc_pres_ir, shoc_presi_ir,
                   shoc_pdel_ir,
                   shoc_wthl_sfc_ir, shoc_wqw_sfc_ir, shoc_uw_sfc_ir, shoc_vw_sfc_ir,
                   shoc_wtracer_sfc_ir, num_qtracers, shoc_w_field_ir,
                   shoc_inv_exner_ir, shoc_phis_ir,
                   shoc_host_dse_ir, shoc_tke_ir, shoc_thetal_ir, shoc_qw_ir,
                   shoc_u_wind_ir, shoc_v_wind_ir, shoc_qtracers_ir,
                   shoc_wthv_sec_ir, shoc_tkh_ir, shoc_tk_ir,
                   shoc_ql_ir, shoc_cldfrac_ir,
                   shoc_pblh_ir,
                   shoc_mix_ir, shoc_isotropy_ir,
                   shoc_w_sec_ir, shoc_thl_sec_ir, shoc_qw_sec_ir, shoc_qwthl_sec_ir,
                   shoc_wthl_sec_ir, shoc_wqw_sec_ir, shoc_wtke_sec_ir,
                   shoc_uw_sec_ir, shoc_vw_sec_ir, shoc_w3_ir,
                   shoc_wqls_sec_ir, shoc_brunt_ir, shoc_ql2_ir);
#endif

    // Process outputs from SHOC (reordering the vertical dimension)
    parallel_for( SimpleBounds<2>(nz,ncol) , YAKL_LAMBDA (int k, int i) {
      int k_shoc = nz-1-k;
      real qw = shoc_qw(k_shoc,i);
      real ql = shoc_ql(k_shoc,i);
      real qv = qw - ql;
      temp    (k,i) = (shoc_thetal(k_shoc,i) + (latvap/cp_d) * shoc_ql(k_shoc,i)) / shoc_inv_exner(k_shoc,i);
      rho_v   (k,i) = qv * rho_d(k,i);
      rho_c   (k,i) = ql * rho_d(k,i);
      uvel    (k,i) = shoc_u_wind(k_shoc,i);
      vvel    (k,i) = shoc_v_wind(k_shoc,i);
      tke     (k,i) = shoc_tke     (k_shoc,i) * rho_d(k,i);
      wthv_sec(k,i) = shoc_wthv_sec(k_shoc,i);
      tk      (k,i) = shoc_tk      (k_shoc,i);
      tkh     (k,i) = shoc_tkh     (k_shoc,i);
      cldfrac (k,i) = std::min(1._fp , shoc_cldfrac (k_shoc,i) );
      for (int tr=0; tr < num_qtracers; tr++) {
        qtracers_pam(tr,k,i) = shoc_qtracers(tr,k_shoc,i) * rho_d(k,i);
      }
      real rcm  = shoc_ql (k_shoc,i);
      real rcm2 = shoc_ql2(k_shoc,i);
      if ( rcm != 0 && rcm2 != 0 ) {
        relvar(k,i) = std::min( 10._fp , std::max( 0.001_fp , rcm*rcm / rcm2 ) );
      } else {
        relvar(k,i) = 1;
      }
    });

    #ifdef PAM_DEBUG
      real mass;
      {
        auto rho_d = dm.get<real,4>( "density_dry" );
        auto rho_v = dm.get<real,4>( "water_vapor" );
        real4d rho_c;
        if      (micro_kessler) { rho_c = dm.get<real,4>( "cloud_liquid" ); }
        else if (micro_p3     ) { rho_c = dm.get<real,4>( "cloud_water"  ); }
        real4d mass4d("mass4d",nz,ny,nx,nens);
        parallel_for( SimpleBounds<4>(nz,ny,nx,nens) , YAKL_LAMBDA (int k, int j, int i, int iens) {
          mass4d(k,j,i,iens) = ( rho_d(k,j,i,iens) + rho_v(k,j,i,iens) + rho_c(k,j,i,iens) ) *
                               crm_dx * crm_dy * (zint_pam(k+1,iens) - zint_pam(k,iens));
        });
        mass = yakl::intrinsics::sum(mass4d);
      }
      std::cout << "Relative mass Change: " << abs(mass-mass0)/mass0 << " : ";
    #endif

    first_step = false;
    etime += dt;
  }


  void finalize(pam::PamCoupler &coupler) {
  }


  std::string sgs_name() const {
    return "shoc";
  }



};



