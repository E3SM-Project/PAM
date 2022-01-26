
#pragma once

#include "awfl_const.h"
#include "DataManager.h"
#include "pam_coupler.h"

using pam::PamCoupler;


extern "C"
void shoc_init(int &nlev, double &gravit, double &rair, double &rh2o, double &cpair, double &zvir, double &latvap,
               double &latice, double &karman, double *pref_mid, int &nbot_shoc, int &ntop_shoc);


extern "C"
void shoc_main(int &shcol, int &nlev, int &nlevi, double *dtime, int &nadv, real *host_dx, double *host_dy,
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

  bool first_step;

  // Indices for all of your tracer quantities
  int static constexpr ID_TKE  = 0;  // Local index for Turbulent Kinetic Energy (m^2/s^2)



  // TODO: Make sure the constants vibe with P3
  // Set constants and likely num_tracers as well, and anything else you can do immediately
  SGS() {
    R_d           = 287.;
    cp_d          = 1003.;
    cv_d          = cp_d - R_d;
    gamma_d       = cp_d / cv_d;
    kappa_d       = R_d  / cp_d;
    R_v           = 461.;
    cp_v          = 1859;
    cv_v          = R_v - cp_v;
    p0            = 1.e5;
    grav          = 9.81;
    first_step    = true;
    cp_l          = 4218.;
    micro_kessler = false;
    micro_p3      = false;
    latvap        = 2.5E6 ;
    latice        = 3.50E5;
    karman        = 0.4;
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
    //                 name    description                              positive   adds mass
    coupler.add_tracer("tke" , "Turbulent Kinetic Energy (m^2/s^2)"   , true     , false );

    // Register and allocation non-tracer quantities used by the microphysics
    int p3_nout = 49;
    coupler.dm.register_and_allocate<real>( "wthv_sec" , "Buoyancy flux [K m/s]"                , {nz,ny,nx,nens} , {"z","y","x","nens"} );
    coupler.dm.register_and_allocate<real>( "tk"       , "Eddy coefficient for momentum [m2/s]" , {nz,ny,nx,nens} , {"z","y","x","nens"} );
    coupler.dm.register_and_allocate<real>( "tkh"      , "Eddy coefficent for heat [m2/s]"      , {nz,ny,nx,nens} , {"z","y","x","nens"} );
    coupler.dm.register_and_allocate<real>( "cldfrac"  , "Cloud fraction [-]"                   , {nz,ny,nx,nens} , {"z","y","x","nens"} );

    auto tke      = coupler.dm.get<real,4>( "tke"      );
    auto wthv_sec = coupler.dm.get<real,4>( "wthv_sec" );
    auto tk       = coupler.dm.get<real,4>( "tk"       );
    auto tkh      = coupler.dm.get<real,4>( "tkh"      );
    auto cldfrac  = coupler.dm.get<real,4>( "cldfrac"  );

    parallel_for( "sgs zero" , SimpleBounds<4>(nz,ny,nx,nens) , YAKL_LAMBDA (int k, int j, int i, int iens) {
      tke     (k,j,i,iens) = 0;
      wthv_sec(k,j,i,iens) = 0;
      tk      (k,j,i,iens) = 0;
      tkh     (k,j,i,iens) = 0;
      cldfrac (k,j,i,iens) = 0;
    });
  }



  void timeStep( PamCoupler &coupler , real dt , real etime ) {
    // Get the dimensions sizes
    int nz   = coupler.get_nz();
    int ny   = coupler.get_ny();
    int nx   = coupler.get_nx();
    int nens = coupler.get_nens();
    int ncol = ny*nx*nens;

    auto pres_mid = coupler.compute_pressure_array(); 
    auto pres_int = coupler.interp_pressure_interfaces( pres );

    // SHOC init requires reference pressure, which we do not have available for the init() call
    if (first_step) {
      // TODO: See if a more appropriate reference pressure should be used
      // pref is only used to limit PBL height to 400mb, so this should be OK
      // Invert the first column in x, z, and ensemble to use as reference pressure for shoc
      real1d pref_shoc("pref_shoc",nz);
      parallel_for( Bounds<2>(nz) , YAKL_LAMBDA (int k, int j, int i, int iens) {
        pref_shoc(k,j,i,iens) = pres_mid(nz-1-k,0,0,0);
      });
      shoc_init( nz , grav , R_d , R_v , cp_d , R_v / R_d -1 , latvap , latice , karman ,
                 pref_shoc.createHostCopy().data() , nlev , 1 );

      // This check is here instead of init because it's not guaranteed the micro has called init before sgs
      if (! coupler.option_exists("micro")) {
        endrun("ERROR: SHOC requires coupler.set_option<std::string>(\"micro\",...) to be set");
      }
      std::string micro_scheme = coupler.get_option<std::string>("micro");
      if      (micro_scheme == "kessler") { micro_kessler = true; }
      else if (micro_scheme == "p3"     ) { micro_p3      = true; }
      else { endrun("ERROR: SHOC only meant to run with kessler or p3 microphysics"); }

      first_step = false;
    }

    // Get saved SHOC-related variables
    auto tke      = coupler.dm.get_lev_col<real>( "tke"      ); // PAM Tracer                 ; don't compute
    auto wthv_sec = coupler.dm.get_lev_col<real>( "wthv_sec" ); // Reuse from last SHOC output; don't compute
    auto tk       = coupler.dm.get_lev_col<real>( "tk"       ); // Reuse from last SHOC output; don't compute
    auto tkh      = coupler.dm.get_lev_col<real>( "tkh"      ); // Reuse from last SHOC output; don't compute
    auto cldfrac  = coupler.dm.get_lev_col<real>( "cldfrac"  ); // Reuse from last SHOC output; don't compute
    // Get coupler state
    auto rho_d        = coupler.dm.get_lev_col<real>( "density_dry"  );
    auto uvel         = coupler.dm.get_lev_col<real>( "uvel"         );
    auto vvel         = coupler.dm.get_lev_col<real>( "vvel"         );
    auto wvel         = coupler.dm.get_lev_col<real>( "wvel"         );
    auto temp         = coupler.dm.get_lev_col<real>( "temp"         );
    auto rho_v        = coupler.dm.get_lev_col<real>( "water_vapor"  ); // Water vapor mass

    // Grab cloud liquid tracer, and determine what other tracers to diffuse in SHOC
    // TODO: Maybe tracers don't have to be mixing ratios? Dividing number concentration by rho_dry doesn't seem to make much sense
    // TODO: Do we add water vapor and cloud liquid to the diffused tracers, or does SHOC do that already?
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

    auto zint = coupler.dm.get<real,2>("vertical_interface_height");
    auto zmid = coupler.dm.get<real,2>("vertical_midpoint_height" );

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
        real z       = zmid     (k_shoc,i);
        real press   = press_mid(k_shoc,i);
        real temp    = temp     (k_shoc,i);
        real qv      = rho_v    (k_shoc,i) / rho_d(k_shoc,i);
        real ql      = rho_c    (k_shoc,i) / rho_d(k_shoc,i);
        real exner   = pow( press / p0 , R_d / cp_d );
        real theta   = temp / exner;
        // https://glossary.ametsoc.org/wiki/Virtual_potential_temperature
        real theta_v = theta * (1 + 0.61_fp * qv - ql);
        // https://glossary.ametsoc.org/wiki/Liquid_water_potential_temperature
        real theta_l = theta - (theta/temp) * (latvap/cp_d) * ql;
        // dry static energy = Cp*T + g*z + phis
        real dse     = cp_d * temp + grav * z + grav * zint(0,i);
        shoc_zt_grid (k,i) = z;
        shoc_pres    (k,i) = press;
        shoc_pdel    (k,i) = pres_int(k_shoc+1,i) - pres_int(k_shoc,i);
        shoc_thv     (k,i) = theta_v;
        shoc_w_field (k,i) = wvel(k_shoc,i);
        shoc_exner   (k,i) = exner;
        shoc_host_dse(k,i) = dse;
        shoc_tke     (k,i) = tke(k_shoc,i);
        shoc_thetal  (k,i) = theta_l;
        shoc_qw      (k,i) = ( rho_v(k_shoc,i) + rho_c(k_shoc,i) ) / rho_d(k_shoc,i);
        shoc_u_wind  (k,i) = uvel(k_shoc,i);
        shoc_v_wind  (k,i) = vvel(k_shoc,i);
        shoc_wthv_sec(k,i) = wthv_sec(k_shoc,i);
        for (int tr=0; tr < num_qtracers; tr++) {
          shoc_qtracers(tr,k,i) = qtracers_pam(tr,k_shoc,i) / rho_d(k_shoc,i);
        }
        shoc_tk     (k,i) = tk     (k_shoc,i);
        shoc_tkh    (k,i) = tkh    (k_shoc,i);
        shoc_cldfrac(k,i) = cldfrac(k_shoc,i);
        shoc_ql     (k,i) = rho_c  (k_shoc,i) / rho_d(k_shoc,i);
      }
      shoc_zi_grid(k,i) = zint    (nz-k,i);
      shoc_presi  (k,i) = pres_int(nz-k,i);
    });

    shoc_main_fortran( ncol, nlev, nlev+1, dt, 1, 
                       shoc_host_dx, shoc_host_dy,shoc_thv, 
                       shoc_zt_grid,shoc_zi_grid,shoc_pres,shoc_presi,shoc_pdel,
                       shoc_wthl_sfc, shoc_wqw_sfc, shoc_uw_sfc, shoc_vw_sfc, 
                       shoc_wtracer_sfc,shoc_num_qtracers,shoc_w_field, 
                       shoc_exner,shoc_phis, 
                       shoc_host_dse, shoc_tke, shoc_thetal, shoc_qw, 
                       shoc_u_wind, shoc_v_wind,shoc_qtracers,
                       shoc_wthv_sec,shoc_tkh,shoc_tk,
                       shoc_ql,shoc_cldfrac,
                       shoc_pblh,
                       shoc_mix, shoc_isotropy,
                       shoc_w_sec, shoc_thl_sec, shoc_qw_sec, shoc_qwthl_sec,
                       shoc_wthl_sec, shoc_wqw_sec, shoc_wtke_sec,
                       shoc_uw_sec, shoc_vw_sec, shoc_w3,
                       shoc_wqls_sec, shoc_brunt, shoc_ql2 );

  }



  // Returns saturation vapor pressure
  // TODO: Make this the same as P3's
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



  std::string micro_name() const {
    return "p3";
  }



};



