
#pragma once

#include "pam_const.h"
#include "DataManager.h"
#include "vertical_interp.h"
#include "YAKL_netcdf.h"


namespace pam {

  YAKL_INLINE real hydrostatic_pressure( real2d const &hy_params , real z_in , real zbot , real ztop , int iens             ) {
    real z = ( z_in - zbot ) / (ztop - zbot);
    real a0 = hy_params(0,iens);
    real a1 = hy_params(1,iens);
    real a2 = hy_params(2,iens);
    real a3 = hy_params(3,iens);
    real a4 = hy_params(4,iens);
    real lnp = a0 + ( a1 + ( a2 + ( a3 + a4 * z ) * z ) * z ) * z;
    return exp(lnp);
  }


  YAKL_INLINE real hydrostatic_density( real2d const &hy_params , real z_in , real zbot , real ztop , int iens , real grav ) {
    real z = ( z_in - zbot ) / (ztop - zbot);
    real a0 = hy_params(0,iens);
    real a1 = hy_params(1,iens);
    real a2 = hy_params(2,iens);
    real a3 = hy_params(3,iens);
    real a4 = hy_params(4,iens);
    real p = hydrostatic_pressure( hy_params , z_in , zbot , ztop , iens );
    real mult = a1 + (2*a2 + (3*a3 + 4*a4*z) * z) * z;
    real dpdz = mult*p/(ztop-zbot);
    return -dpdz/grav;
  }


  class PamCoupler {
    public:

    real R_d;
    real R_v;
    real grav;

    DataManager dm;


    PamCoupler() {
      R_d  = 287.;
      R_v  = 461.;
      grav = 9.8;
    }



    PamCoupler(real R_d, real R_v) {
      this->R_d = R_d;
      this->R_v = R_v;
    }



    inline void set_gas_constants(real R_d, real R_v) {
      this->R_d = R_d;
      this->R_v = R_v;
    }


    inline void set_vertical_grid(real2d const &zint_in) {
      int nz   = dm.get_dimension_size("z");
      int nens = dm.get_dimension_size("nens");
      auto zint = dm.get<real,2>("vertical_interface_height");
      auto zmid = dm.get<real,2>("vertical_midpoint_height" );
      parallel_for( SimpleBounds<2>(nz+1,nens) , YAKL_LAMBDA (int k, int iens) {
        zint(k,iens) = zint_in(k,iens);
        if (k < nz) zmid(k,iens) = 0.5_fp * (zint_in(k,iens) + zint_in(k+1,iens));
      });
    }


    inline void set_vertical_grid(real1d const &zint_in) {
      int nz   = dm.get_dimension_size("z");
      int nens = dm.get_dimension_size("nens");
      auto zint = dm.get<real,2>("vertical_interface_height");
      auto zmid = dm.get<real,2>("vertical_midpoint_height" );
      parallel_for( SimpleBounds<2>(nz+1,nens) , YAKL_LAMBDA (int k, int iens) {
        zint(k,iens) = zint_in(k);
        if (k < nz) zmid(k,iens) = 0.5_fp * (zint_in(k) + zint_in(k+1));
      });
    }


    inline void allocate_coupler_state( int nz, int ny, int nx, int nens ) {
      dm.register_and_allocate<real>( "density_dry"               , "dry density"               , {nz,ny,nx,nens} , {"z","y","x","nens"} );
      dm.register_and_allocate<real>( "uvel"                      , "x-direction velocity"      , {nz,ny,nx,nens} , {"z","y","x","nens"} );
      dm.register_and_allocate<real>( "vvel"                      , "y-direction velocity"      , {nz,ny,nx,nens} , {"z","y","x","nens"} );
      dm.register_and_allocate<real>( "wvel"                      , "z-direction velocity"      , {nz,ny,nx,nens} , {"z","y","x","nens"} );
      dm.register_and_allocate<real>( "temp"                      , "temperature"               , {nz,ny,nx,nens} , {"z","y","x","nens"} );
      dm.register_and_allocate<real>( "diagnostic_pressure"       , "pressure"                  , {nz,ny,nx,nens} , {"z","y","x","nens"} );
      dm.register_and_allocate<real>( "vertical_interface_height" , "vertical interface height" , {nz+1    ,nens} , {"zp1"      ,"nens"} );
      dm.register_and_allocate<real>( "vertical_midpoint_height"  , "vertical midpoint height"  , {nz      ,nens} , {"z"        ,"nens"} );
      dm.register_and_allocate<real>( "hydrostasis_parameters"    , "hydrostasis parameters"    , {5       ,nens} , {"five"     ,"nens"} );
      dm.register_and_allocate<real>( "hydrostatic_pressure"      , "hydrostasis pressure"      , {nz      ,nens} , {"z"        ,"nens"} );
      dm.register_and_allocate<real>( "hydrostatic_density"       , "hydrostasis density"       , {nz      ,nens} , {"z"        ,"nens"} );

      auto density_dry  = dm.get_collapsed<real>("density_dry"              );
      auto uvel         = dm.get_collapsed<real>("uvel"                     );
      auto vvel         = dm.get_collapsed<real>("vvel"                     );
      auto wvel         = dm.get_collapsed<real>("wvel"                     );
      auto temp         = dm.get_collapsed<real>("temp"                     );
      auto diag_press   = dm.get_collapsed<real>("diagnostic_pressure"      );
      auto zint         = dm.get_collapsed<real>("vertical_interface_height");
      auto zmid         = dm.get_collapsed<real>("vertical_midpoint_height" );
      auto hy_params    = dm.get_collapsed<real>("hydrostasis_parameters"   );
      auto hy_press     = dm.get_collapsed<real>("hydrostatic_pressure"     );
      auto hy_dens      = dm.get_collapsed<real>("hydrostatic_density"      );

      parallel_for( Bounds<1>(nz*ny*nx*nens) , YAKL_LAMBDA (int i) {
        density_dry (i) = 0;
        uvel        (i) = 0;
        vvel        (i) = 0;
        wvel        (i) = 0;
        temp        (i) = 0;
        diag_press  (i) = 0;
        if (i < (nz+1)*nens) zint(i) = 0;
        if (i < (nz  )*nens) {
          zmid    (i) = 0;
          hy_press(i) = 0;
          hy_dens (i) = 0;
        }
        if (i < 5     *nens) hy_params(i) = 0;
      });
    }


    YAKL_INLINE static real compute_pressure( real rho_d, real rho_v, real T, real R_d, real R_v ) {
      return rho_d*R_d*T + rho_v*R_v*T;
    }


    inline void update_diagnostic_pressure( ) {
      YAKL_SCOPE( R_d , this->R_d );
      YAKL_SCOPE( R_v , this->R_v );

      auto dens_dry = dm.get_lev_col<real>("density_dry");
      auto dens_wv  = dm.get_lev_col<real>("water_vapor");
      auto temp     = dm.get_lev_col<real>("temp");
      auto pressure = dm.get_lev_col<real>("diagnostic_pressure");

      int nz   = dens_dry.dimension[0];
      int ncol = dens_dry.dimension[1];

      parallel_for( SimpleBounds<2>(nz,ncol) , YAKL_LAMBDA (int k, int i) {
        real rho_d = dens_dry(k,i);
        real rho_v = dens_wv (k,i);
        real T     = temp    (k,i);
        pressure(k,i) = compute_pressure( rho_d , rho_v , T , R_d , R_v );
      });
    }


    inline real4d interp_pressure_interfaces( ) {
      auto zint      = dm.get<real,2>("vertical_interface_height");
      auto press     = dm.get<real,4>("diagnostic_pressure");
      auto hy_press  = dm.get<real,2>("hydrostatic_pressure");
      auto hy_params = dm.get<real,2>("hydrostasis_parameters");

      int nz   = dm.get_dimension_size("z");
      int ny   = dm.get_dimension_size("y");
      int nx   = dm.get_dimension_size("x");
      int nens = dm.get_dimension_size("nens");

      real4d press_pert("press_pert",nz,ny,nx,nens);

      // Compute pressure perturbation
      parallel_for( SimpleBounds<4>(nz,ny,nx,nens) , YAKL_LAMBDA (int k, int j, int i, int iens) {
        press_pert(k,j,i,iens) = press(k,j,i,iens) - hy_press(k,iens);
      });

      // Interpolate pressure perturbation from cells to edges
      VerticalInterp<5> vert_interp;
      vert_interp.init(zint);
      auto press_edges = vert_interp.cells_to_edges( press_pert ,
                                                     vert_interp.BC_ZERO_GRADIENT ,
                                                     vert_interp.BC_ZERO_GRADIENT );

      // Add hydrostasis at cell edges to get back full pressure
      parallel_for( SimpleBounds<4>(nz+1,ny,nx,nens) , YAKL_LAMBDA (int k, int j, int i, int iens) {
        press_edges(k,j,i,iens) += hydrostatic_pressure( hy_params , zint(k,iens) , zint(0,iens) , zint(nz,iens) , iens );
      });
      
      return press_edges;
    }


    inline void update_hydrostasis( ) {
      update_diagnostic_pressure();

      auto zint      = dm.get<real,2>("vertical_interface_height");
      auto zmid      = dm.get<real,2>("vertical_midpoint_height" );
      auto pressure  = dm.get<real,4>("diagnostic_pressure"      );
      auto hy_params = dm.get<real,2>("hydrostasis_parameters"   );
      auto hy_press  = dm.get<real,2>("hydrostatic_pressure"     );
      auto hy_dens   = dm.get<real,2>("hydrostatic_density"      );

      int nz   = dm.get_dimension_size("z");
      int ny   = dm.get_dimension_size("y");
      int nx   = dm.get_dimension_size("x");
      int nens = dm.get_dimension_size("nens");

      int k0 = 0;
      int k1 = nz/4;
      int k2 = nz/2;
      int k3 = (3*nz)/4;
      int k4 = nz-1;

      parallel_for( "hydro fitting" , SimpleBounds<1>(nens) , YAKL_LAMBDA (int iens) {
        real zbot = zint(0 ,iens);
        real ztop = zint(nz,iens);

        SArray<double,1,5> z;
        z(0) = ( zmid(k0,iens) - zbot ) / (ztop - zbot);
        z(1) = ( zmid(k1,iens) - zbot ) / (ztop - zbot);
        z(2) = ( zmid(k2,iens) - zbot ) / (ztop - zbot);
        z(3) = ( zmid(k3,iens) - zbot ) / (ztop - zbot);
        z(4) = ( zmid(k4,iens) - zbot ) / (ztop - zbot);

        SArray<double,2,5,5> vand;
        for (int j=0; j < 5; j++) {
          for (int i=0; i < 5; i++) {
            vand(j,i) = pow( z(i) , (double) j );
          }
        }

        auto vand_inv = matinv_ge_cr( vand );

        // Fit to just one column, assuming all columns are fairly similar
        // This will only be used for idealized test cases anyway
        // Another function that passes in an ensemble of pressure profiles will be
        //   used for GCM coupling in an MMF setting
        SArray<double,1,5> logp;
        logp(0) = log(pressure(k0,0,0,iens));
        logp(1) = log(pressure(k1,0,0,iens));
        logp(2) = log(pressure(k2,0,0,iens));
        logp(3) = log(pressure(k3,0,0,iens));
        logp(4) = log(pressure(k4,0,0,iens));

        auto params = matmul_cr( vand_inv , logp );

        for (int i=0; i < 5; i++) { hy_params(i,iens) = params(i); }
      });

      YAKL_SCOPE( grav , this->grav );

      // Compute hydrostatic pressure and density as cell averages
      SArray<real,1,9> gll_pts, gll_wts;
      get_gll_points ( gll_pts );
      get_gll_weights( gll_wts );
      parallel_for( "hydro pressure" , SimpleBounds<2>(nz,nens) , YAKL_LAMBDA (int k, int iens) {
        real p = 0;
        real d = 0;
        for (int kk=0; kk < 9; kk++) {
          real dz = zint(k+1,iens) - zint(k,iens);
          real zloc = zint(k,iens) + 0.5_fp*dz + gll_pts(kk)*dz;
          real wt = gll_wts(kk);
          p += hydrostatic_pressure( hy_params , zloc , zint(0,iens) , zint(nz,iens) , iens        ) * wt;
          d += hydrostatic_density ( hy_params , zloc , zint(0,iens) , zint(nz,iens) , iens , grav ) * wt;
        }
        hy_press(k,iens) = p;
        hy_dens (k,iens) = d;
      });

    }



    template <class FP> YAKL_INLINE void get_gll_points(SArray<FP,1,9> &rslt) {
      rslt(0)=-0.50000000000000000000000000000000000000;
      rslt(1)=-0.44987899770573007865617262220916897903;
      rslt(2)=-0.33859313975536887672294271354567122536;
      rslt(3)=-0.18155873191308907935537603435432960651;
      rslt(4)=0.00000000000000000000000000000000000000;
      rslt(5)=0.18155873191308907935537603435432960651;
      rslt(6)=0.33859313975536887672294271354567122536;
      rslt(7)=0.44987899770573007865617262220916897903;
      rslt(8)=0.50000000000000000000000000000000000000;
    }



    template <class FP> YAKL_INLINE void get_gll_weights(SArray<FP,1,9> &rslt) {
      rslt(0)=0.013888888888888888888888888888888888889;
      rslt(1)=0.082747680780402762523169860014604152919;
      rslt(2)=0.13726935625008086764035280928968636297;
      rslt(3)=0.17321425548652317255756576606985914397;
      rslt(4)=0.18575963718820861678004535147392290249;
      rslt(5)=0.17321425548652317255756576606985914397;
      rslt(6)=0.13726935625008086764035280928968636297;
      rslt(7)=0.082747680780402762523169860014604152919;
      rslt(8)=0.013888888888888888888888888888888888889;
    }



  };

}


