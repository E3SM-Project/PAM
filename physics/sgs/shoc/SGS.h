
#pragma once

#include "awfl_const.h"
#include "DataManager.h"
#include "pam_coupler.h"
#include "call_shoc_from_pam.h"
#include "pam_scream_routines.h"

#define SHOC_FORTRAN
// #define SHOC_DEBUG

using pam::PamCoupler;


template <class T, int N, int memSpace1, int memSpace2, int style>
yakl::Array<T,N,memHost,style> operator-(yakl::Array<T,N,memSpace1,style> a, yakl::Array<T,N,memSpace2,style> b) {
  int n = a.totElems();
  yakl::Array<T,N,memHost,style> ret;
  ret = a.createHostCopy();
  for (int i=0; i < n; i++) { ret.myData[i] = a.myData[i] - b.myData[i]; }
  return ret;
}


template <class T, int N, int memSpace1, int style>
yakl::Array<T,N,memHost,style> operator+(yakl::Array<T,N,memSpace1,style> a, T b) {
  int n = a.totElems();
  yakl::Array<T,N,memHost,style> ret;
  ret = a.createHostCopy();
  for (int i=0; i < n; i++) { ret.myData[i] = a.myData[i] + b; }
  return ret;
}


template <class T, int N, int memSpace1, int memSpace2, int style>
yakl::Array<T,N,memHost,style> operator/(yakl::Array<T,N,memSpace1,style> a, yakl::Array<T,N,memSpace2,style> b) {
  int n = a.totElems();
  yakl::Array<T,N,memHost,style> ret;
  ret = a.createHostCopy();
  for (int i=0; i < n; i++) { ret.myData[i] = a.myData[i] / b.myData[i]; }
  return ret;
}


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
  // Doesn't actually have to be static or constexpr. Could be assigned in the constructor
  int static constexpr num_tracers = 1;

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



  // TODO: Make sure the constants vibe with the rest of the model physics
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
    grav          = 9.81;
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
  YAKL_INLINE static int get_num_tracers() {
    return num_tracers;
  }



  // Can do whatever you want, but mainly for registering tracers and allocating data
  void init(PamCoupler &coupler) {
    int nx   = coupler.get_nx  ();
    int ny   = coupler.get_ny  ();
    int nz   = coupler.get_nz  ();
    int nens = coupler.get_nens();

    pam::kokkos_initialize();

    if ( coupler.get_option<bool>("advect_tke") ) {
      // Register tracers in the coupler
      //                 name    description                              positive   adds mass
      coupler.add_tracer("tke" , "Turbulent Kinetic Energy (m^2/s^2)"   , true     , false );
      std::cout << "ADVECTING_TKE\n";
    } else {
      coupler.dm.register_and_allocate<real>( "tke" , "Turbulent Kinetic Energy (m^2/s^2)" , {nz,ny,nx,nens} , {"z","y","x","nens"} );
      std::cout << "*NOT* ADVECTING_TKE\n";
    }

    // Register and allocation non-tracer quantities used by the microphysics
    coupler.dm.register_and_allocate<real>( "wthv_sec" , "Buoyancy flux [K m/s]"                , {nz,ny,nx,nens} , {"z","y","x","nens"} );
    coupler.dm.register_and_allocate<real>( "tk"       , "Eddy coefficient for momentum [m2/s]" , {nz,ny,nx,nens} , {"z","y","x","nens"} );
    coupler.dm.register_and_allocate<real>( "tkh"      , "Eddy coefficent for heat [m2/s]"      , {nz,ny,nx,nens} , {"z","y","x","nens"} );
    coupler.dm.register_and_allocate<real>( "cldfrac"  , "Cloud fraction [-]"                   , {nz,ny,nx,nens} , {"z","y","x","nens"} );
    coupler.dm.register_and_allocate<real>( "relvar"   , "Relative cloud water variance"        , {nz,ny,nx,nens} , {"z","y","x","nens"} );

    auto tke      = coupler.dm.get<real,4>( "tke"      );
    auto wthv_sec = coupler.dm.get<real,4>( "wthv_sec" );
    auto tk       = coupler.dm.get<real,4>( "tk"       );
    auto tkh      = coupler.dm.get<real,4>( "tkh"      );
    auto cldfrac  = coupler.dm.get<real,4>( "cldfrac"  );
    auto relvar   = coupler.dm.get<real,4>( "relvar"   );

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



  void timeStep( PamCoupler &coupler , real dt ) {
    // Get the dimensions sizes
    int nz   = coupler.get_nz();
    int ny   = coupler.get_ny();
    int nx   = coupler.get_nx();
    int nens = coupler.get_nens();
    int ncol = ny*nx*nens;

    auto pres_mid_tmp = coupler.compute_pressure_array();
    auto pres_int_tmp = coupler.interp_pressure_interfaces( pres_mid_tmp );

    auto pres_mid = pres_mid_tmp.reshape<2>({nz  ,ncol});
    auto pres_int = pres_int_tmp.reshape<2>({nz+1,ncol});

    // SHOC init requires reference pressure, which we do not have available for the init() call
    if (first_step) {
      // TODO: See if a more appropriate reference pressure should be used
      // pref is only used to limit PBL height to 400mb, so this should be OK
      // Invert the first column in x, z, and ensemble to use as reference pressure for shoc
      real1d pref_shoc("pref_shoc",nz);
      parallel_for( Bounds<1>(nz) , YAKL_LAMBDA (int k) {
        pref_shoc(k) = pres_mid(nz-1-k,0);
      });
      real zvir = R_v / R_d - 1;
      int kbot, ktop;
      
      // #ifdef SHOC_FORTRAN

        kbot = nz;
        ktop = 1;
        shoc_init_fortran( nz , grav , R_d , R_v , cp_d , zvir , latvap , latice , karman ,
                           pref_shoc.createHostCopy().data() , kbot , ktop );

      // #else

        kbot = nz-1;
        ktop = 0;
        this->npbl = pam::call_shoc_init_from_pam( kbot , ktop , pam::yakl_array_to_arrayIR( pref_shoc ) );

      // #endif

      // This check is here instead of init because it's not guaranteed the micro has called init before sgs
      if (! coupler.option_exists("micro")) {
        endrun("ERROR: SHOC requires coupler.set_option<std::string>(\"micro\",...) to be set");
      }
      std::string micro_scheme = coupler.get_option<std::string>("micro");
      if      (micro_scheme == "kessler") { micro_kessler = true; }
      else if (micro_scheme == "p3"     ) { micro_p3      = true; }
      else { endrun("ERROR: SHOC only meant to run with kessler or p3 microphysics"); }
    }

    // Get saved SHOC-related variables
    auto tke      = coupler.dm.get_lev_col<real>( "tke"      ); // PAM Tracer                 ; don't compute
    auto wthv_sec = coupler.dm.get_lev_col<real>( "wthv_sec" ); // Reuse from last SHOC output; don't compute
    auto tk       = coupler.dm.get_lev_col<real>( "tk"       ); // Reuse from last SHOC output; don't compute
    auto tkh      = coupler.dm.get_lev_col<real>( "tkh"      ); // Reuse from last SHOC output; don't compute
    auto cldfrac  = coupler.dm.get_lev_col<real>( "cldfrac"  ); // Reuse from last SHOC output; don't compute
    auto relvar   = coupler.dm.get_lev_col<real>( "relvar"   ); // Computed on output for P3
    // Get coupler state
    auto rho_d        = coupler.dm.get_lev_col<real>( "density_dry"  );
    auto uvel         = coupler.dm.get_lev_col<real>( "uvel"         );
    auto vvel         = coupler.dm.get_lev_col<real>( "vvel"         );
    auto wvel         = coupler.dm.get_lev_col<real>( "wvel"         );
    auto temp         = coupler.dm.get_lev_col<real>( "temp"         );
    auto rho_v        = coupler.dm.get_lev_col<real>( "water_vapor"  ); // Water vapor mass

    // Grab cloud liquid tracer, and determine what other tracers to diffuse in SHOC
    // TODO: Do we add water vapor and cloud liquid to the diffused tracers, or does SHOC do that already?
    //       I'm pretty sure SHOC does it already for qv and qc. qw=qv+qc is diffused, and other handing is applied to ql
    // TODO: If we add new MMF modules that add other tracers that need to be diffused, then this situation
    //       must be handled
    MultiField<real,2> qtracers_pam;  // Extra tracers for SHOC to diffuse
    real2d rho_c;                     // Cloud liquid mass (name differs for different micro schemes)
    if        (micro_kessler) {
      rho_c = coupler.dm.get_lev_col<real>( "cloud_liquid" );
      qtracers_pam.add_field( coupler.dm.get_lev_col<real>("precip_liquid") );
    } else if (micro_p3     ) {
      rho_c = coupler.dm.get_lev_col<real>( "cloud_water"  );
      qtracers_pam.add_field( coupler.dm.get_lev_col<real>("cloud_water_num") );
      qtracers_pam.add_field( coupler.dm.get_lev_col<real>("rain"           ) );
      qtracers_pam.add_field( coupler.dm.get_lev_col<real>("rain_num"       ) );
      qtracers_pam.add_field( coupler.dm.get_lev_col<real>("ice"            ) );
      qtracers_pam.add_field( coupler.dm.get_lev_col<real>("ice_num"        ) );
      qtracers_pam.add_field( coupler.dm.get_lev_col<real>("ice_rime"       ) );
      qtracers_pam.add_field( coupler.dm.get_lev_col<real>("ice_rime_vol"   ) );
    }

    int num_qtracers = qtracers_pam.get_num_fields();

    auto zint_pam = coupler.dm.get<real,2>("vertical_interface_height");
    auto zmid_pam = coupler.dm.get<real,2>("vertical_midpoint_height" );

    real4d zint_tmp("zint_tmp",nz+1,ny,nx,nens);
    real4d zmid_tmp("zint_tmp",nz  ,ny,nx,nens);

    parallel_for( Bounds<4>(nz+1,ny,nx,nens) , YAKL_LAMBDA (int k, int j, int i, int iens) {
      zint_tmp(k,j,i,iens) = zint_pam(k,iens);
      if (k < nz) zmid_tmp(k,j,i,iens) = zmid_pam(k,iens);
    });

    auto zint = zint_tmp.reshape<2>({nz+1,ncol});
    auto zmid = zmid_tmp.reshape<2>({nz  ,ncol});

    real crm_dx = coupler.get_xlen() / nx;
    real crm_dy;
    if (ny == 1) { crm_dy = crm_dx;                  }
    else         { crm_dy = coupler.get_ylen() / ny; }

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
    real cp_d   = this->cp_d  ;
    real latvap = this->latvap;

    // Compute inputs for SHOC (reordering the vertical dimension)
    parallel_for( Bounds<2>(nz+1,ncol) , YAKL_LAMBDA (int k, int i) {
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
        real z       = zmid    (k_shoc,i);
        real press   = pres_mid(k_shoc,i);
        real t       = temp    (k_shoc,i);
        real qv      = rho_v   (k_shoc,i) / rho_d(k_shoc,i);
        real ql      = rho_c   (k_shoc,i) / rho_d(k_shoc,i);
        real exner   = pow( press / p0 , R_d / cp_d );
        real theta   = t / exner;
        // https://glossary.ametsoc.org/wiki/Virtual_potential_temperature
        // TODO: Figure out which of these SHOC wants
        // real theta_v = theta * (1 + 0.61_fp * qv - ql);
        real theta_v = theta * (1 + 0.61_fp * qv);
        // https://glossary.ametsoc.org/wiki/Liquid_water_potential_temperature
        // According to update_host_dse, the simplified version is used here
        real theta_l = theta - (latvap/cp_d) * ql;
        // dry static energy = Cp*T + g*z + phis
        real dse     = cp_d * t + grav * z;
        shoc_zt_grid  (k,i) = z;
        shoc_pres     (k,i) = press;
        shoc_pdel     (k,i) = pres_int(k_shoc,i) - pres_int(k_shoc+1,i);
        shoc_thv      (k,i) = theta_v;
        shoc_w_field  (k,i) = wvel(k_shoc,i);
        shoc_exner    (k,i) = exner;
        shoc_inv_exner(k,i) = 1._fp / exner;
        shoc_host_dse (k,i) = dse;
        shoc_tke      (k,i) = tke(k_shoc,i);
        shoc_thetal   (k,i) = theta_l;
        shoc_qw       (k,i) = ( rho_v(k_shoc,i) + rho_c(k_shoc,i) ) / rho_d(k_shoc,i);
        shoc_u_wind   (k,i) = uvel(k_shoc,i);
        shoc_v_wind   (k,i) = vvel(k_shoc,i);
        shoc_wthv_sec (k,i) = wthv_sec(k_shoc,i);
        for (int tr=0; tr < num_qtracers; tr++) {
          // TODO: I'm pretty sure from the application of vd_shoc_solve in shoc that this can be total mass
          shoc_qtracers(tr,k,i) = qtracers_pam(tr,k_shoc,i);
        }
        shoc_tk     (k,i) = tk     (k_shoc,i);
        shoc_tkh    (k,i) = tkh    (k_shoc,i);
        shoc_cldfrac(k,i) = cldfrac(k_shoc,i);
        shoc_ql     (k,i) = rho_c  (k_shoc,i) / rho_d(k_shoc,i);
      }
      int k_shoc = nz-k;
      shoc_zi_grid(k,i) = zint    (k_shoc,i);
      shoc_presi  (k,i) = pres_int(k_shoc,i);
    });

    int nadv = 1;



    #ifdef SHOC_DEBUG
    if (etime > 5) {

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
                         shoc_zt_grid_host.data(), shoc_zi_grid_host.data(), shoc_pres_host.data(), shoc_presi_host.data(), shoc_pdel_host.data(),
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

      pam::call_shoc_main_from_pam( ncol, nz, nz+1, dt, nadv, num_qtracers, this->npbl,
                                    pam::yakl_array_to_arrayIR( shoc_host_dx     ),
                                    pam::yakl_array_to_arrayIR( shoc_host_dy     ),
                                    pam::yakl_array_to_arrayIR( shoc_thv         ),
                                    pam::yakl_array_to_arrayIR( shoc_zt_grid     ),
                                    pam::yakl_array_to_arrayIR( shoc_zi_grid     ),
                                    pam::yakl_array_to_arrayIR( shoc_pres        ),
                                    pam::yakl_array_to_arrayIR( shoc_presi       ),
                                    pam::yakl_array_to_arrayIR( shoc_pdel        ),
                                    pam::yakl_array_to_arrayIR( shoc_wthl_sfc    ),
                                    pam::yakl_array_to_arrayIR( shoc_wqw_sfc     ),
                                    pam::yakl_array_to_arrayIR( shoc_uw_sfc      ),
                                    pam::yakl_array_to_arrayIR( shoc_vw_sfc      ),
                                    pam::yakl_array_to_arrayIR( shoc_wtracer_sfc ),
                                    pam::yakl_array_to_arrayIR( shoc_w_field     ),
                                    pam::yakl_array_to_arrayIR( shoc_inv_exner   ),
                                    pam::yakl_array_to_arrayIR( shoc_phis        ),
                                    pam::yakl_array_to_arrayIR( shoc_host_dse    ),
                                    pam::yakl_array_to_arrayIR( shoc_tke         ),
                                    pam::yakl_array_to_arrayIR( shoc_thetal      ),
                                    pam::yakl_array_to_arrayIR( shoc_qw          ),
                                    pam::yakl_array_to_arrayIR( shoc_u_wind      ),
                                    pam::yakl_array_to_arrayIR( shoc_v_wind      ),
                                    pam::yakl_array_to_arrayIR( shoc_qtracers    ),
                                    pam::yakl_array_to_arrayIR( shoc_wthv_sec    ),
                                    pam::yakl_array_to_arrayIR( shoc_tkh         ),
                                    pam::yakl_array_to_arrayIR( shoc_tk          ),
                                    pam::yakl_array_to_arrayIR( shoc_ql          ),
                                    pam::yakl_array_to_arrayIR( shoc_cldfrac     ),
                                    pam::yakl_array_to_arrayIR( shoc_pblh        ),
                                    pam::yakl_array_to_arrayIR( shoc_mix         ),
                                    pam::yakl_array_to_arrayIR( shoc_isotropy    ),
                                    pam::yakl_array_to_arrayIR( shoc_w_sec       ),
                                    pam::yakl_array_to_arrayIR( shoc_thl_sec     ),
                                    pam::yakl_array_to_arrayIR( shoc_qw_sec      ),
                                    pam::yakl_array_to_arrayIR( shoc_qwthl_sec   ),
                                    pam::yakl_array_to_arrayIR( shoc_wthl_sec    ),
                                    pam::yakl_array_to_arrayIR( shoc_wqw_sec     ),
                                    pam::yakl_array_to_arrayIR( shoc_wtke_sec    ),
                                    pam::yakl_array_to_arrayIR( shoc_uw_sec      ),
                                    pam::yakl_array_to_arrayIR( shoc_vw_sec      ),
                                    pam::yakl_array_to_arrayIR( shoc_w3          ),
                                    pam::yakl_array_to_arrayIR( shoc_wqls_sec    ),
                                    pam::yakl_array_to_arrayIR( shoc_brunt       ),
                                    pam::yakl_array_to_arrayIR( shoc_ql2         ) );

      using yakl::intrinsics::abs;
      using yakl::intrinsics::sum;
      using std::abs;



      yakl::SimpleNetCDF nc; 
      nc.create("debug_diff.nc");
      nc.write( (shoc_host_dx.createHostCopy()   - shoc_host_dx_host  ) / (abs(shoc_host_dx_host  ) + 1.e-50) , "shoc_host_dx"     , {       "ncol"} );
      nc.write( (shoc_host_dy.createHostCopy()   - shoc_host_dy_host  ) / (abs(shoc_host_dy_host  ) + 1.e-50) , "shoc_host_dy"     , {       "ncol"} );
      nc.write( (shoc_zt_grid.createHostCopy()   - shoc_zt_grid_host  ) / (abs(shoc_zt_grid_host  ) + 1.e-50) , "shoc_zt_grid"     , {"nz"  ,"ncol"} );
      nc.write( (shoc_zi_grid.createHostCopy()   - shoc_zi_grid_host  ) / (abs(shoc_zi_grid_host  ) + 1.e-50) , "shoc_zi_grid"     , {"nzp1","ncol"} );
      nc.write( (shoc_pres.createHostCopy()      - shoc_pres_host     ) / (abs(shoc_pres_host     ) + 1.e-50) , "shoc_pres"        , {"nz"  ,"ncol"} );
      nc.write( (shoc_presi.createHostCopy()     - shoc_presi_host    ) / (abs(shoc_presi_host    ) + 1.e-50) , "shoc_presi"       , {"nzp1","ncol"} );
      nc.write( (shoc_pdel.createHostCopy()      - shoc_pdel_host     ) / (abs(shoc_pdel_host     ) + 1.e-50) , "shoc_pdel"        , {"nz"  ,"ncol"} );
      nc.write( (shoc_thv.createHostCopy()       - shoc_thv_host      ) / (abs(shoc_thv_host      ) + 1.e-50) , "shoc_thv"         , {"nz"  ,"ncol"} );
      nc.write( (shoc_w_field.createHostCopy()   - shoc_w_field_host  ) / (abs(shoc_w_field_host  ) + 1.e-50) , "shoc_w_field"     , {"nz"  ,"ncol"} );
      nc.write( (shoc_wthl_sfc.createHostCopy()  - shoc_wthl_sfc_host ) / (abs(shoc_wthl_sfc_host ) + 1.e-50) , "shoc_wthl_sfc"    , {       "ncol"} );
      nc.write( (shoc_wqw_sfc.createHostCopy()   - shoc_wqw_sfc_host  ) / (abs(shoc_wqw_sfc_host  ) + 1.e-50) , "shoc_wqw_sfc"     , {       "ncol"} );
      nc.write( (shoc_uw_sfc.createHostCopy()    - shoc_uw_sfc_host   ) / (abs(shoc_uw_sfc_host   ) + 1.e-50) , "shoc_uw_sfc"      , {       "ncol"} );
      nc.write( (shoc_vw_sfc.createHostCopy()    - shoc_vw_sfc_host   ) / (abs(shoc_vw_sfc_host   ) + 1.e-50) , "shoc_vw_sfc"      , {       "ncol"} );
      nc.write( (shoc_exner.createHostCopy()     - shoc_exner_host    ) / (abs(shoc_exner_host    ) + 1.e-50) , "shoc_exner"       , {"nz"  ,"ncol"} );
      nc.write( (shoc_inv_exner.createHostCopy() - shoc_inv_exner_host) / (abs(shoc_inv_exner_host) + 1.e-50) , "shoc_inv_exner"   , {"nz"  ,"ncol"} );
      nc.write( (shoc_phis.createHostCopy()      - shoc_phis_host     ) / (abs(shoc_phis_host     ) + 1.e-50) , "shoc_phis"        , {       "ncol"} );
      nc.write( (shoc_host_dse.createHostCopy()  - shoc_host_dse_host ) / (abs(shoc_host_dse_host ) + 1.e-50) , "shoc_host_dse"    , {"nz"  ,"ncol"} );
      nc.write( (shoc_tke.createHostCopy()       - shoc_tke_host      ) / (abs(shoc_tke_host      ) + 1.e-50) , "shoc_tke"         , {"nz"  ,"ncol"} );
      nc.write( (shoc_thetal.createHostCopy()    - shoc_thetal_host   ) / (abs(shoc_thetal_host   ) + 1.e-50) , "shoc_thetal"      , {"nz"  ,"ncol"} );
      nc.write( (shoc_qw.createHostCopy()        - shoc_qw_host       ) / (abs(shoc_qw_host       ) + 1.e-50) , "shoc_qw"          , {"nz"  ,"ncol"} );
      nc.write( (shoc_u_wind.createHostCopy()    - shoc_u_wind_host   ) / (abs(shoc_u_wind_host   ) + 1.e-50) , "shoc_u_wind"      , {"nz"  ,"ncol"} );
      nc.write( (shoc_v_wind.createHostCopy()    - shoc_v_wind_host   ) / (abs(shoc_v_wind_host   ) + 1.e-50) , "shoc_v_wind"      , {"nz"  ,"ncol"} );
      nc.write( (shoc_wthv_sec.createHostCopy()  - shoc_wthv_sec_host ) / (abs(shoc_wthv_sec_host ) + 1.e-50) , "shoc_wthv_sec"    , {"nz"  ,"ncol"} );
      nc.write( (shoc_tk.createHostCopy()        - shoc_tk_host       ) / (abs(shoc_tk_host       ) + 1.e-50) , "shoc_tk"          , {"nz"  ,"ncol"} );
      nc.write( (shoc_cldfrac.createHostCopy()   - shoc_cldfrac_host  ) / (abs(shoc_cldfrac_host  ) + 1.e-50) , "shoc_cldfrac"     , {"nz"  ,"ncol"} );
      nc.write( (shoc_ql.createHostCopy()        - shoc_ql_host       ) / (abs(shoc_ql_host       ) + 1.e-50) , "shoc_ql"          , {"nz"  ,"ncol"} );
      nc.write( (shoc_pblh.createHostCopy()      - shoc_pblh_host     ) / (abs(shoc_pblh_host     ) + 1.e-50) , "shoc_pblh"        , {       "ncol"} );
      nc.write( (shoc_ql2.createHostCopy()       - shoc_ql2_host      ) / (abs(shoc_ql2_host      ) + 1.e-50) , "shoc_ql2"         , {"nz"  ,"ncol"} );
      nc.write( (shoc_mix.createHostCopy()       - shoc_mix_host      ) / (abs(shoc_mix_host      ) + 1.e-50) , "shoc_mix"         , {"nz"  ,"ncol"} );
      nc.write( (shoc_w_sec.createHostCopy()     - shoc_w_sec_host    ) / (abs(shoc_w_sec_host    ) + 1.e-50) , "shoc_w_sec"       , {"nz"  ,"ncol"} );
      nc.write( (shoc_thl_sec.createHostCopy()   - shoc_thl_sec_host  ) / (abs(shoc_thl_sec_host  ) + 1.e-50) , "shoc_thl_sec"     , {"nzp1","ncol"} );
      nc.write( (shoc_qw_sec.createHostCopy()    - shoc_qw_sec_host   ) / (abs(shoc_qw_sec_host   ) + 1.e-50) , "shoc_qw_sec"      , {"nzp1","ncol"} );
      nc.write( (shoc_qwthl_sec.createHostCopy() - shoc_qwthl_sec_host) / (abs(shoc_qwthl_sec_host) + 1.e-50) , "shoc_qwthl_sec"   , {"nzp1","ncol"} );
      nc.write( (shoc_wthl_sec.createHostCopy()  - shoc_wthl_sec_host ) / (abs(shoc_wthl_sec_host ) + 1.e-50) , "shoc_wthl_sec"    , {"nzp1","ncol"} );
      nc.write( (shoc_wqw_sec.createHostCopy()   - shoc_wqw_sec_host  ) / (abs(shoc_wqw_sec_host  ) + 1.e-50) , "shoc_wqw_sec"     , {"nzp1","ncol"} );
      nc.write( (shoc_wtke_sec.createHostCopy()  - shoc_wtke_sec_host ) / (abs(shoc_wtke_sec_host ) + 1.e-50) , "shoc_wtke_sec"    , {"nzp1","ncol"} );
      nc.write( (shoc_uw_sec.createHostCopy()    - shoc_uw_sec_host   ) / (abs(shoc_uw_sec_host   ) + 1.e-50) , "shoc_uw_sec"      , {"nzp1","ncol"} );
      nc.write( (shoc_vw_sec.createHostCopy()    - shoc_vw_sec_host   ) / (abs(shoc_vw_sec_host   ) + 1.e-50) , "shoc_vw_sec"      , {"nzp1","ncol"} );
      nc.write( (shoc_w3.createHostCopy()        - shoc_w3_host       ) / (abs(shoc_w3_host       ) + 1.e-50) , "shoc_w3"          , {"nzp1","ncol"} );
      nc.write( (shoc_wqls_sec.createHostCopy()  - shoc_wqls_sec_host ) / (abs(shoc_wqls_sec_host ) + 1.e-50) , "shoc_wqls_sec"    , {"nz"  ,"ncol"} );
      nc.write( (shoc_brunt.createHostCopy()     - shoc_brunt_host    ) / (abs(shoc_brunt_host    ) + 1.e-50) , "shoc_brunt"       , {"nz"  ,"ncol"} );
      nc.write( (shoc_isotropy.createHostCopy()  - shoc_isotropy_host ) / (abs(shoc_isotropy_host ) + 1.e-50) , "shoc_isotropy"    , {"nz"  ,"ncol"} );
      nc.close();




      nc.create("debug_ekat.nc");
      nc.write( shoc_host_dx.createHostCopy()   , "shoc_host_dx"     , {       "ncol"} );
      nc.write( shoc_host_dy.createHostCopy()   , "shoc_host_dy"     , {       "ncol"} );
      nc.write( shoc_zt_grid.createHostCopy()   , "shoc_zt_grid"     , {"nz"  ,"ncol"} );
      nc.write( shoc_zi_grid.createHostCopy()   , "shoc_zi_grid"     , {"nzp1","ncol"} );
      nc.write( shoc_pres.createHostCopy()      , "shoc_pres"        , {"nz"  ,"ncol"} );
      nc.write( shoc_presi.createHostCopy()     , "shoc_presi"       , {"nzp1","ncol"} );
      nc.write( shoc_pdel.createHostCopy()      , "shoc_pdel"        , {"nz"  ,"ncol"} );
      nc.write( shoc_thv.createHostCopy()       , "shoc_thv"         , {"nz"  ,"ncol"} );
      nc.write( shoc_w_field.createHostCopy()   , "shoc_w_field"     , {"nz"  ,"ncol"} );
      nc.write( shoc_wthl_sfc.createHostCopy()  , "shoc_wthl_sfc"    , {       "ncol"} );
      nc.write( shoc_wqw_sfc.createHostCopy()   , "shoc_wqw_sfc"     , {       "ncol"} );
      nc.write( shoc_uw_sfc.createHostCopy()    , "shoc_uw_sfc"      , {       "ncol"} );
      nc.write( shoc_vw_sfc.createHostCopy()    , "shoc_vw_sfc"      , {       "ncol"} );
      nc.write( shoc_exner.createHostCopy()     , "shoc_exner"       , {"nz"  ,"ncol"} );
      nc.write( shoc_inv_exner.createHostCopy() , "shoc_inv_exner"   , {"nz"  ,"ncol"} );
      nc.write( shoc_phis.createHostCopy()      , "shoc_phis"        , {       "ncol"} );
      nc.write( shoc_host_dse.createHostCopy()  , "shoc_host_dse"    , {"nz"  ,"ncol"} );
      nc.write( shoc_tke.createHostCopy()       , "shoc_tke"         , {"nz"  ,"ncol"} );
      nc.write( shoc_thetal.createHostCopy()    , "shoc_thetal"      , {"nz"  ,"ncol"} );
      nc.write( shoc_qw.createHostCopy()        , "shoc_qw"          , {"nz"  ,"ncol"} );
      nc.write( shoc_u_wind.createHostCopy()    , "shoc_u_wind"      , {"nz"  ,"ncol"} );
      nc.write( shoc_v_wind.createHostCopy()    , "shoc_v_wind"      , {"nz"  ,"ncol"} );
      nc.write( shoc_wthv_sec.createHostCopy()  , "shoc_wthv_sec"    , {"nz"  ,"ncol"} );
      nc.write( shoc_tk.createHostCopy()        , "shoc_tk"          , {"nz"  ,"ncol"} );
      nc.write( shoc_cldfrac.createHostCopy()   , "shoc_cldfrac"     , {"nz"  ,"ncol"} );
      nc.write( shoc_ql.createHostCopy()        , "shoc_ql"          , {"nz"  ,"ncol"} );
      nc.write( shoc_pblh.createHostCopy()      , "shoc_pblh"        , {       "ncol"} );
      nc.write( shoc_ql2.createHostCopy()       , "shoc_ql2"         , {"nz"  ,"ncol"} );
      nc.write( shoc_mix.createHostCopy()       , "shoc_mix"         , {"nz"  ,"ncol"} );
      nc.write( shoc_w_sec.createHostCopy()     , "shoc_w_sec"       , {"nz"  ,"ncol"} );
      nc.write( shoc_thl_sec.createHostCopy()   , "shoc_thl_sec"     , {"nzp1","ncol"} );
      nc.write( shoc_qw_sec.createHostCopy()    , "shoc_qw_sec"      , {"nzp1","ncol"} );
      nc.write( shoc_qwthl_sec.createHostCopy() , "shoc_qwthl_sec"   , {"nzp1","ncol"} );
      nc.write( shoc_wthl_sec.createHostCopy()  , "shoc_wthl_sec"    , {"nzp1","ncol"} );
      nc.write( shoc_wqw_sec.createHostCopy()   , "shoc_wqw_sec"     , {"nzp1","ncol"} );
      nc.write( shoc_wtke_sec.createHostCopy()  , "shoc_wtke_sec"    , {"nzp1","ncol"} );
      nc.write( shoc_uw_sec.createHostCopy()    , "shoc_uw_sec"      , {"nzp1","ncol"} );
      nc.write( shoc_vw_sec.createHostCopy()    , "shoc_vw_sec"      , {"nzp1","ncol"} );
      nc.write( shoc_w3.createHostCopy()        , "shoc_w3"          , {"nzp1","ncol"} );
      nc.write( shoc_wqls_sec.createHostCopy()  , "shoc_wqls_sec"    , {"nz"  ,"ncol"} );
      nc.write( shoc_brunt.createHostCopy()     , "shoc_brunt"       , {"nz"  ,"ncol"} );
      nc.write( shoc_isotropy.createHostCopy()  , "shoc_isotropy"    , {"nz"  ,"ncol"} );
      nc.close();



      nc.create("debug_fortran.nc");
      nc.write( shoc_host_dx_host  , "shoc_host_dx"     , {       "ncol"} );
      nc.write( shoc_host_dy_host  , "shoc_host_dy"     , {       "ncol"} );
      nc.write( shoc_zt_grid_host  , "shoc_zt_grid"     , {"nz"  ,"ncol"} );
      nc.write( shoc_zi_grid_host  , "shoc_zi_grid"     , {"nzp1","ncol"} );
      nc.write( shoc_pres_host     , "shoc_pres"        , {"nz"  ,"ncol"} );
      nc.write( shoc_presi_host    , "shoc_presi"       , {"nzp1","ncol"} );
      nc.write( shoc_pdel_host     , "shoc_pdel"        , {"nz"  ,"ncol"} );
      nc.write( shoc_thv_host      , "shoc_thv"         , {"nz"  ,"ncol"} );
      nc.write( shoc_w_field_host  , "shoc_w_field"     , {"nz"  ,"ncol"} );
      nc.write( shoc_wthl_sfc_host , "shoc_wthl_sfc"    , {       "ncol"} );
      nc.write( shoc_wqw_sfc_host  , "shoc_wqw_sfc"     , {       "ncol"} );
      nc.write( shoc_uw_sfc_host   , "shoc_uw_sfc"      , {       "ncol"} );
      nc.write( shoc_vw_sfc_host   , "shoc_vw_sfc"      , {       "ncol"} );
      nc.write( shoc_exner_host    , "shoc_exner"       , {"nz"  ,"ncol"} );
      nc.write( shoc_inv_exner_host, "shoc_inv_exner"   , {"nz"  ,"ncol"} );
      nc.write( shoc_phis_host     , "shoc_phis"        , {       "ncol"} );
      nc.write( shoc_host_dse_host , "shoc_host_dse"    , {"nz"  ,"ncol"} );
      nc.write( shoc_tke_host      , "shoc_tke"         , {"nz"  ,"ncol"} );
      nc.write( shoc_thetal_host   , "shoc_thetal"      , {"nz"  ,"ncol"} );
      nc.write( shoc_qw_host       , "shoc_qw"          , {"nz"  ,"ncol"} );
      nc.write( shoc_u_wind_host   , "shoc_u_wind"      , {"nz"  ,"ncol"} );
      nc.write( shoc_v_wind_host   , "shoc_v_wind"      , {"nz"  ,"ncol"} );
      nc.write( shoc_wthv_sec_host , "shoc_wthv_sec"    , {"nz"  ,"ncol"} );
      nc.write( shoc_tk_host       , "shoc_tk"          , {"nz"  ,"ncol"} );
      nc.write( shoc_cldfrac_host  , "shoc_cldfrac"     , {"nz"  ,"ncol"} );
      nc.write( shoc_ql_host       , "shoc_ql"          , {"nz"  ,"ncol"} );
      nc.write( shoc_pblh_host     , "shoc_pblh"        , {       "ncol"} );
      nc.write( shoc_ql2_host      , "shoc_ql2"         , {"nz"  ,"ncol"} );
      nc.write( shoc_mix_host      , "shoc_mix"         , {"nz"  ,"ncol"} );
      nc.write( shoc_w_sec_host    , "shoc_w_sec"       , {"nz"  ,"ncol"} );
      nc.write( shoc_thl_sec_host  , "shoc_thl_sec"     , {"nzp1","ncol"} );
      nc.write( shoc_qw_sec_host   , "shoc_qw_sec"      , {"nzp1","ncol"} );
      nc.write( shoc_qwthl_sec_host, "shoc_qwthl_sec"   , {"nzp1","ncol"} );
      nc.write( shoc_wthl_sec_host , "shoc_wthl_sec"    , {"nzp1","ncol"} );
      nc.write( shoc_wqw_sec_host  , "shoc_wqw_sec"     , {"nzp1","ncol"} );
      nc.write( shoc_wtke_sec_host , "shoc_wtke_sec"    , {"nzp1","ncol"} );
      nc.write( shoc_uw_sec_host   , "shoc_uw_sec"      , {"nzp1","ncol"} );
      nc.write( shoc_vw_sec_host   , "shoc_vw_sec"      , {"nzp1","ncol"} );
      nc.write( shoc_w3_host       , "shoc_w3"          , {"nzp1","ncol"} );
      nc.write( shoc_wqls_sec_host , "shoc_wqls_sec"    , {"nz"  ,"ncol"} );
      nc.write( shoc_brunt_host    , "shoc_brunt"       , {"nz"  ,"ncol"} );
      nc.write( shoc_isotropy_host , "shoc_isotropy"    , {"nz"  ,"ncol"} );
      nc.close();




      real adiff_shoc_host_dx     = 0;
      real adiff_shoc_host_dy     = 0;
      real adiff_shoc_zt_grid     = 0;
      real adiff_shoc_zi_grid     = 0;
      real adiff_shoc_pres        = 0;
      real adiff_shoc_presi       = 0;
      real adiff_shoc_pdel        = 0;
      real adiff_shoc_thv         = 0;
      real adiff_shoc_w_field     = 0;
      real adiff_shoc_wthl_sfc    = 0;
      real adiff_shoc_wqw_sfc     = 0;
      real adiff_shoc_uw_sfc      = 0;
      real adiff_shoc_vw_sfc      = 0;
      real adiff_shoc_wtracer_sfc = 0;
      real adiff_shoc_exner       = 0;
      real adiff_shoc_inv_exner   = 0;
      real adiff_shoc_phis        = 0;
      real adiff_shoc_host_dse    = 0;
      real adiff_shoc_tke         = 0;
      real adiff_shoc_thetal      = 0;
      real adiff_shoc_qw          = 0;
      real adiff_shoc_u_wind      = 0;
      real adiff_shoc_v_wind      = 0;
      real adiff_shoc_wthv_sec    = 0;
      real adiff_shoc_qtracers    = 0;
      real adiff_shoc_tk          = 0;
      real adiff_shoc_tkh         = 0;
      real adiff_shoc_cldfrac     = 0;
      real adiff_shoc_ql          = 0;
      real adiff_shoc_pblh        = 0;
      real adiff_shoc_ql2         = 0;
      real adiff_shoc_mix         = 0;
      real adiff_shoc_w_sec       = 0;
      real adiff_shoc_thl_sec     = 0;
      real adiff_shoc_qw_sec      = 0;
      real adiff_shoc_qwthl_sec   = 0;
      real adiff_shoc_wthl_sec    = 0;
      real adiff_shoc_wqw_sec     = 0;
      real adiff_shoc_wtke_sec    = 0;
      real adiff_shoc_uw_sec      = 0;
      real adiff_shoc_vw_sec      = 0;
      real adiff_shoc_w3          = 0;
      real adiff_shoc_wqls_sec    = 0;
      real adiff_shoc_brunt       = 0;
      real adiff_shoc_isotropy    = 0;
      for (int k=0; k < nz+1; k++) {
        for (int i=0; i < ncol; i++) {
          adiff_shoc_thl_sec   += abs(shoc_thl_sec  (k,i) - shoc_thl_sec_host  (k,i));
          adiff_shoc_qw_sec    += abs(shoc_qw_sec   (k,i) - shoc_qw_sec_host   (k,i));
          adiff_shoc_qwthl_sec += abs(shoc_qwthl_sec(k,i) - shoc_qwthl_sec_host(k,i));
          adiff_shoc_wthl_sec  += abs(shoc_wthl_sec (k,i) - shoc_wthl_sec_host (k,i));
          adiff_shoc_wqw_sec   += abs(shoc_wqw_sec  (k,i) - shoc_wqw_sec_host  (k,i));
          adiff_shoc_wtke_sec  += abs(shoc_wtke_sec (k,i) - shoc_wtke_sec_host (k,i));
          adiff_shoc_uw_sec    += abs(shoc_uw_sec   (k,i) - shoc_uw_sec_host   (k,i));
          adiff_shoc_vw_sec    += abs(shoc_vw_sec   (k,i) - shoc_vw_sec_host   (k,i));
          adiff_shoc_w3        += abs(shoc_w3       (k,i) - shoc_w3_host       (k,i));
          adiff_shoc_zi_grid   += abs(shoc_zi_grid  (k,i) - shoc_zi_grid_host  (k,i));
          adiff_shoc_presi     += abs(shoc_presi    (k,i) - shoc_presi_host    (k,i));
          if (k < nz) {
            adiff_shoc_zt_grid   += abs(shoc_zt_grid  (k,i) - shoc_zt_grid_host  (k,i));
            adiff_shoc_pres      += abs(shoc_pres     (k,i) - shoc_pres_host     (k,i));
            adiff_shoc_pdel      += abs(shoc_pdel     (k,i) - shoc_pdel_host     (k,i));
            adiff_shoc_thv       += abs(shoc_thv      (k,i) - shoc_thv_host      (k,i));
            adiff_shoc_w_field   += abs(shoc_w_field  (k,i) - shoc_w_field_host  (k,i));
            adiff_shoc_exner     += abs(shoc_exner    (k,i) - shoc_exner_host    (k,i));
            adiff_shoc_inv_exner += abs(shoc_inv_exner(k,i) - shoc_inv_exner_host(k,i));
            adiff_shoc_host_dse  += abs(shoc_host_dse (k,i) - shoc_host_dse_host (k,i));
            adiff_shoc_tke       += abs(shoc_tke      (k,i) - shoc_tke_host      (k,i));
            adiff_shoc_thetal    += abs(shoc_thetal   (k,i) - shoc_thetal_host   (k,i));
            adiff_shoc_qw        += abs(shoc_qw       (k,i) - shoc_qw_host       (k,i));
            adiff_shoc_u_wind    += abs(shoc_u_wind   (k,i) - shoc_u_wind_host   (k,i));
            adiff_shoc_v_wind    += abs(shoc_v_wind   (k,i) - shoc_v_wind_host   (k,i));
            adiff_shoc_wthv_sec  += abs(shoc_wthv_sec (k,i) - shoc_wthv_sec_host (k,i));
            adiff_shoc_tk        += abs(shoc_tk       (k,i) - shoc_tk_host       (k,i));
            adiff_shoc_tkh       += abs(shoc_tkh      (k,i) - shoc_tkh_host      (k,i));
            adiff_shoc_cldfrac   += abs(shoc_cldfrac  (k,i) - shoc_cldfrac_host  (k,i));
            adiff_shoc_ql        += abs(shoc_ql       (k,i) - shoc_ql_host       (k,i));
            adiff_shoc_ql2       += abs(shoc_ql2      (k,i) - shoc_ql2_host      (k,i));
            adiff_shoc_mix       += abs(shoc_mix      (k,i) - shoc_mix_host      (k,i));
            adiff_shoc_w_sec     += abs(shoc_w_sec    (k,i) - shoc_w_sec_host    (k,i));
            adiff_shoc_wqls_sec  += abs(shoc_wqls_sec (k,i) - shoc_wqls_sec_host (k,i));
            adiff_shoc_brunt     += abs(shoc_brunt    (k,i) - shoc_brunt_host    (k,i));
            adiff_shoc_isotropy  += abs(shoc_isotropy (k,i) - shoc_isotropy_host (k,i));
            for (int tr=0; tr < num_qtracers; tr++) {
              adiff_shoc_qtracers += abs(shoc_qtracers(tr,k,i) - shoc_qtracers_host(tr,k,i));
            }
          }
          if (k == 0) {
            adiff_shoc_host_dx  += abs(shoc_host_dx (i) - shoc_host_dx_host (i));
            adiff_shoc_host_dy  += abs(shoc_host_dy (i) - shoc_host_dy_host (i));
            adiff_shoc_wthl_sfc += abs(shoc_wthl_sfc(i) - shoc_wthl_sfc_host(i));
            adiff_shoc_wqw_sfc  += abs(shoc_wqw_sfc (i) - shoc_wqw_sfc_host (i));
            adiff_shoc_uw_sfc   += abs(shoc_uw_sfc  (i) - shoc_uw_sfc_host  (i));
            adiff_shoc_vw_sfc   += abs(shoc_vw_sfc  (i) - shoc_vw_sfc_host  (i));
            adiff_shoc_phis     += abs(shoc_phis    (i) - shoc_phis_host    (i));
            adiff_shoc_pblh     += abs(shoc_pblh    (i) - shoc_pblh_host    (i));
            for (int tr=0; tr < num_qtracers; tr++) {
              adiff_shoc_wtracer_sfc += abs(shoc_wtracer_sfc(tr,i) - shoc_wtracer_sfc_host(tr,i));
            }
          }
        }
      }
      std::cout << "shoc_host_dx    : " << std::scientific << adiff_shoc_host_dx     /sum(abs(shoc_host_dx_host     )) << std::endl;
      std::cout << "shoc_host_dy    : " << std::scientific << adiff_shoc_host_dy     /sum(abs(shoc_host_dy_host     )) << std::endl;
      std::cout << "shoc_zt_grid    : " << std::scientific << adiff_shoc_zt_grid     /sum(abs(shoc_zt_grid_host     )) << std::endl;
      std::cout << "shoc_zi_grid    : " << std::scientific << adiff_shoc_zi_grid     /sum(abs(shoc_zi_grid_host     )) << std::endl;
      std::cout << "shoc_pres       : " << std::scientific << adiff_shoc_pres        /sum(abs(shoc_pres_host        )) << std::endl;
      std::cout << "shoc_presi      : " << std::scientific << adiff_shoc_presi       /sum(abs(shoc_presi_host       )) << std::endl;
      std::cout << "shoc_pdel       : " << std::scientific << adiff_shoc_pdel        /sum(abs(shoc_pdel_host        )) << std::endl;
      std::cout << "shoc_thv        : " << std::scientific << adiff_shoc_thv         /sum(abs(shoc_thv_host         )) << std::endl;
      std::cout << "shoc_w_field    : " << std::scientific << adiff_shoc_w_field     /sum(abs(shoc_w_field_host     )) << std::endl;
      std::cout << "shoc_wthl_sfc   : " << std::scientific << adiff_shoc_wthl_sfc    /sum(abs(shoc_wthl_sfc_host    )) << std::endl;
      std::cout << "shoc_wqw_sfc    : " << std::scientific << adiff_shoc_wqw_sfc     /sum(abs(shoc_wqw_sfc_host     )) << std::endl;
      std::cout << "shoc_uw_sfc     : " << std::scientific << adiff_shoc_uw_sfc      /sum(abs(shoc_uw_sfc_host      )) << std::endl;
      std::cout << "shoc_vw_sfc     : " << std::scientific << adiff_shoc_vw_sfc      /sum(abs(shoc_vw_sfc_host      )) << std::endl;
      std::cout << "shoc_wtracer_sfc: " << std::scientific << adiff_shoc_wtracer_sfc /sum(abs(shoc_wtracer_sfc_host )) << std::endl;
      std::cout << "shoc_exner      : " << std::scientific << adiff_shoc_exner       /sum(abs(shoc_exner_host       )) << std::endl;
      std::cout << "shoc_inv_exner  : " << std::scientific << adiff_shoc_inv_exner   /sum(abs(shoc_inv_exner_host   )) << std::endl;
      std::cout << "shoc_phis       : " << std::scientific << adiff_shoc_phis        /sum(abs(shoc_phis_host        )) << std::endl;
      std::cout << "shoc_host_dse   : " << std::scientific << adiff_shoc_host_dse    /sum(abs(shoc_host_dse_host    )) << std::endl;
      std::cout << "shoc_tke        : " << std::scientific << adiff_shoc_tke         /sum(abs(shoc_tke_host         )) << std::endl;
      std::cout << "shoc_thetal     : " << std::scientific << adiff_shoc_thetal      /sum(abs(shoc_thetal_host      )) << std::endl;
      std::cout << "shoc_qw         : " << std::scientific << adiff_shoc_qw          /sum(abs(shoc_qw_host          )) << std::endl;
      std::cout << "shoc_u_wind     : " << std::scientific << adiff_shoc_u_wind      /sum(abs(shoc_u_wind_host      )) << std::endl;
      std::cout << "shoc_v_wind     : " << std::scientific << adiff_shoc_v_wind      /sum(abs(shoc_v_wind_host      )) << std::endl;
      std::cout << "shoc_wthv_sec   : " << std::scientific << adiff_shoc_wthv_sec    /sum(abs(shoc_wthv_sec_host    )) << std::endl;
      std::cout << "shoc_qtracers   : " << std::scientific << adiff_shoc_qtracers    /sum(abs(shoc_qtracers_host    )) << std::endl;
      std::cout << "shoc_tk         : " << std::scientific << adiff_shoc_tk          /sum(abs(shoc_tk_host          )) << std::endl;
      std::cout << "shoc_tkh        : " << std::scientific << adiff_shoc_tkh         /sum(abs(shoc_tkh_host         )) << std::endl;
      std::cout << "shoc_cldfrac    : " << std::scientific << adiff_shoc_cldfrac     /sum(abs(shoc_cldfrac_host     )) << std::endl;
      std::cout << "shoc_ql         : " << std::scientific << adiff_shoc_ql          /sum(abs(shoc_ql_host          )) << std::endl;
      std::cout << "shoc_pblh       : " << std::scientific << adiff_shoc_pblh        /sum(abs(shoc_pblh_host        )) << std::endl;
      std::cout << "shoc_ql2        : " << std::scientific << adiff_shoc_ql2         /sum(abs(shoc_ql2_host         )) << std::endl;
      std::cout << "shoc_mix        : " << std::scientific << adiff_shoc_mix         /sum(abs(shoc_mix_host         )) << std::endl;
      std::cout << "shoc_w_sec      : " << std::scientific << adiff_shoc_w_sec       /sum(abs(shoc_w_sec_host       )) << std::endl;
      std::cout << "shoc_thl_sec    : " << std::scientific << adiff_shoc_thl_sec     /sum(abs(shoc_thl_sec_host     )) << std::endl;
      std::cout << "shoc_qw_sec     : " << std::scientific << adiff_shoc_qw_sec      /sum(abs(shoc_qw_sec_host      )) << std::endl;
      std::cout << "shoc_qwthl_sec  : " << std::scientific << adiff_shoc_qwthl_sec   /sum(abs(shoc_qwthl_sec_host   )) << std::endl;
      std::cout << "shoc_wthl_sec   : " << std::scientific << adiff_shoc_wthl_sec    /sum(abs(shoc_wthl_sec_host    )) << std::endl;
      std::cout << "shoc_wqw_sec    : " << std::scientific << adiff_shoc_wqw_sec     /sum(abs(shoc_wqw_sec_host     )) << std::endl;
      std::cout << "shoc_wtke_sec   : " << std::scientific << adiff_shoc_wtke_sec    /sum(abs(shoc_wtke_sec_host    )) << std::endl;
      std::cout << "shoc_uw_sec     : " << std::scientific << adiff_shoc_uw_sec      /sum(abs(shoc_uw_sec_host      )) << std::endl;
      std::cout << "shoc_vw_sec     : " << std::scientific << adiff_shoc_vw_sec      /sum(abs(shoc_vw_sec_host      )) << std::endl;
      std::cout << "shoc_w3         : " << std::scientific << adiff_shoc_w3          /sum(abs(shoc_w3_host          )) << std::endl;
      std::cout << "shoc_wqls_sec   : " << std::scientific << adiff_shoc_wqls_sec    /sum(abs(shoc_wqls_sec_host    )) << std::endl;
      std::cout << "shoc_brunt      : " << std::scientific << adiff_shoc_brunt       /sum(abs(shoc_brunt_host       )) << std::endl;
      std::cout << "shoc_isotropy   : " << std::scientific << adiff_shoc_isotropy    /sum(abs(shoc_isotropy_host    )) << std::endl;
      abort();
    }
    #endif



    #ifdef SHOC_FORTRAN

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
                         shoc_zt_grid_host.data(), shoc_zi_grid_host.data(), shoc_pres_host.data(), shoc_presi_host.data(), shoc_pdel_host.data(),
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

    #else

      pam::call_shoc_main_from_pam( ncol, nz, nz+1, dt, nadv, num_qtracers, this->npbl,
                                    pam::yakl_array_to_arrayIR( shoc_host_dx     ),
                                    pam::yakl_array_to_arrayIR( shoc_host_dy     ),
                                    pam::yakl_array_to_arrayIR( shoc_thv         ),
                                    pam::yakl_array_to_arrayIR( shoc_zt_grid     ),
                                    pam::yakl_array_to_arrayIR( shoc_zi_grid     ),
                                    pam::yakl_array_to_arrayIR( shoc_pres        ),
                                    pam::yakl_array_to_arrayIR( shoc_presi       ),
                                    pam::yakl_array_to_arrayIR( shoc_pdel        ),
                                    pam::yakl_array_to_arrayIR( shoc_wthl_sfc    ),
                                    pam::yakl_array_to_arrayIR( shoc_wqw_sfc     ),
                                    pam::yakl_array_to_arrayIR( shoc_uw_sfc      ),
                                    pam::yakl_array_to_arrayIR( shoc_vw_sfc      ),
                                    pam::yakl_array_to_arrayIR( shoc_wtracer_sfc ),
                                    pam::yakl_array_to_arrayIR( shoc_w_field     ),
                                    pam::yakl_array_to_arrayIR( shoc_inv_exner   ),
                                    pam::yakl_array_to_arrayIR( shoc_phis        ),
                                    pam::yakl_array_to_arrayIR( shoc_host_dse    ),
                                    pam::yakl_array_to_arrayIR( shoc_tke         ),
                                    pam::yakl_array_to_arrayIR( shoc_thetal      ),
                                    pam::yakl_array_to_arrayIR( shoc_qw          ),
                                    pam::yakl_array_to_arrayIR( shoc_u_wind      ),
                                    pam::yakl_array_to_arrayIR( shoc_v_wind      ),
                                    pam::yakl_array_to_arrayIR( shoc_qtracers    ),
                                    pam::yakl_array_to_arrayIR( shoc_wthv_sec    ),
                                    pam::yakl_array_to_arrayIR( shoc_tkh         ),
                                    pam::yakl_array_to_arrayIR( shoc_tk          ),
                                    pam::yakl_array_to_arrayIR( shoc_ql          ),
                                    pam::yakl_array_to_arrayIR( shoc_cldfrac     ),
                                    pam::yakl_array_to_arrayIR( shoc_pblh        ),
                                    pam::yakl_array_to_arrayIR( shoc_mix         ),
                                    pam::yakl_array_to_arrayIR( shoc_isotropy    ),
                                    pam::yakl_array_to_arrayIR( shoc_w_sec       ),
                                    pam::yakl_array_to_arrayIR( shoc_thl_sec     ),
                                    pam::yakl_array_to_arrayIR( shoc_qw_sec      ),
                                    pam::yakl_array_to_arrayIR( shoc_qwthl_sec   ),
                                    pam::yakl_array_to_arrayIR( shoc_wthl_sec    ),
                                    pam::yakl_array_to_arrayIR( shoc_wqw_sec     ),
                                    pam::yakl_array_to_arrayIR( shoc_wtke_sec    ),
                                    pam::yakl_array_to_arrayIR( shoc_uw_sec      ),
                                    pam::yakl_array_to_arrayIR( shoc_vw_sec      ),
                                    pam::yakl_array_to_arrayIR( shoc_w3          ),
                                    pam::yakl_array_to_arrayIR( shoc_wqls_sec    ),
                                    pam::yakl_array_to_arrayIR( shoc_brunt       ),
                                    pam::yakl_array_to_arrayIR( shoc_ql2         ) );

    #endif

    // Process outputs from SHOC (reordering the vertical dimension)
    parallel_for( Bounds<2>(nz,ncol) , YAKL_LAMBDA (int k, int i) {
      int k_shoc = nz-1-k;
      // TODO: What about rho_dry ??
      uvel(k,i) = shoc_u_wind(k_shoc,i);
      vvel(k,i) = shoc_v_wind(k_shoc,i);
      // TODO: What about wvel ??
      temp(k,i) = ( shoc_host_dse(k_shoc,i) - grav * zmid(k,i) ) / cp_d;
      rho_v(k,i) = shoc_qw(k_shoc,i) * rho_d(k,i) - shoc_ql(k_shoc,i) * rho_d(k,i);
      rho_c(k,i) = shoc_ql(k_shoc,i) * rho_d(k,i);
      for (int tr=0; tr < num_qtracers; tr++) {
        qtracers_pam(tr,k,i) = shoc_qtracers(tr,k_shoc,i);
      }
      tke     (k,i) = shoc_tke     (k_shoc,i);
      wthv_sec(k,i) = shoc_wthv_sec(k_shoc,i);
      tk      (k,i) = shoc_tk      (k_shoc,i);
      tkh     (k,i) = shoc_tkh     (k_shoc,i);
      cldfrac (k,i) = min(1._fp , shoc_cldfrac (k_shoc,i) );
      real rcm  = shoc_ql (k_shoc,i);
      real rcm2 = shoc_ql2(k_shoc,i);
      if ( rcm != 0 && rcm2 != 0 ) {
        relvar(k,i) = min( 10._fp , max( 0.001_fp , rcm*rcm / rcm2 ) );
      } else {
        relvar(k,i) = 1;
      }
    });

    first_step = false;
    etime += dt;
  }


  void finalize(PamCoupler &coupler) {
    pam::kokkos_finalize();
  }


  std::string sgs_name() const {
    return "shoc";
  }



};



