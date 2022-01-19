
#pragma once

#include "pam_const.h"
#include "DataManager.h"
#include "vertical_interp.h"
#include "YAKL_netcdf.h"
#include "Notes.h"


namespace pam {



  YAKL_INLINE real hydrostatic_pressure( realConst2d hy_params , real z_in , real zbot , real ztop ,
                                         int iens             ) {
    real z = ( z_in - zbot ) / (ztop - zbot);
    real a0 = hy_params(0,iens);
    real a1 = hy_params(1,iens);
    real a2 = hy_params(2,iens);
    real a3 = hy_params(3,iens);
    real a4 = hy_params(4,iens);
    real a5 = hy_params(5,iens);
    real a6 = hy_params(6,iens);
    real a7 = hy_params(7,iens);
    real a8 = hy_params(8,iens);
    real a9 = hy_params(9,iens);
    real lnp = a0 + ( a1 + ( a2 + ( a3 + ( a4 + ( a5 + ( a6 + ( a7 + ( a8 + a9*z)*z)*z)*z)*z)*z)*z)*z)*z;
    return exp(lnp);
  }



  YAKL_INLINE real hydrostatic_density( realConst2d hy_params , real z_in , real zbot , real ztop ,
                                        int iens , real grav ) {
    real z = ( z_in - zbot ) / (ztop - zbot);
    real a1 = hy_params(1,iens);
    real a2 = hy_params(2,iens);
    real a3 = hy_params(3,iens);
    real a4 = hy_params(4,iens);
    real a5 = hy_params(5,iens);
    real a6 = hy_params(6,iens);
    real a7 = hy_params(7,iens);
    real a8 = hy_params(8,iens);
    real a9 = hy_params(9,iens);
    real p = hydrostatic_pressure( hy_params , z_in , zbot , ztop , iens );
    real mult = a1 + (2*a2 + (3*a3 + (4*a4 + (5*a5 + (6*a6 + (7*a7 + (8*a8 + 9*a9*z)*z)*z)*z)*z)*z)*z)*z;
    real dpdz = mult*p/(ztop-zbot);
    return -dpdz/grav;
  }



  YAKL_INLINE real hydrostatic_density_deriv( realConst2d hy_params , real z_in , real zbot , real ztop ,
                                              int iens , real grav ) {
    real z = ( z_in - zbot ) / (ztop - zbot);
    real a0 = hy_params(0,iens);
    real a1 = hy_params(1,iens);
    real a2 = hy_params(2,iens);
    real a3 = hy_params(3,iens);
    real a4 = hy_params(4,iens);
    real a5 = hy_params(5,iens);
    real a6 = hy_params(6,iens);
    real a7 = hy_params(7,iens);
    real a8 = hy_params(8,iens);
    real a9 = hy_params(9,iens);
    real z2 = z *z;
    real z3 = z2*z;
    real z4 = z3*z;
    real z5 = z4*z;
    real z6 = z5*z;
    real z7 = z6*z;
    real z8 = z7*z;
    real z9 = z8*z;
    real tmp1 = 9*a9*z8 + 8*a8*z7 + 7*a7*z6 + 6*a6*z5 + 5*a5*z4 + 4*a4*z3 + 3*a3*z2 + 2*a2*z + a1;
    real tmp2 = exp(a9*z9 + a8*z8 + a7*z7 + a6*z6 + a5*z5 + a4*z4 + a3*z3 + a2*z2 + a1*z + a0);
    real tmp3 = 2*(36*a9*z7 + 28*a8*z6 + 21*a7*z5 + 15*a6*z4 + 10*a5*z3 + 6*a4*z2 + 3*a3*z + a2);
    real tmp4 = tmp1*tmp1 * tmp2 + tmp3 * tmp2;
    return -tmp4/grav/(ztop-zbot)/(ztop-zbot);
  }



  YAKL_INLINE real compute_pressure( real rho_d, real rho_v, real T, real R_d, real R_v ) {
    return rho_d*R_d*T + rho_v*R_v*T;
  }



  class PamCoupler {
    public:

    Notes notes;

    real R_d;    // Dry air gas constant
    real R_v;    // Water vapor gas constant
    real cp_d;   // Dry air specific heat at constant pressure
    real cp_v;   // Water vapor specific heat at constant pressure
    real grav;   // Acceleration due to gravity (m s^-2): typically 9.81
    real p0;     // Reference pressure (Pa): typically 10^5
    real xlen;   // Domain length in the x-direction in meters
    real ylen;   // Domain length in the y-direction in meters
    real dt_gcm; // Time step of the GCM for this MMF invocation

    DataManager dm;


    PamCoupler() {
      this->R_d  = 287 ;
      this->R_v  = 461 ;
      this->cp_d = 1004;
      this->cp_v = 1859;
      this->grav = 9.81;
      this->p0   = 1.e5;
      this->xlen = -1;
      this->ylen = -1;
    }


    ~PamCoupler() {
      dm.finalize();
      notes.finalize();
    }


    void set_dt_gcm(real dt_gcm) { this->dt_gcm = dt_gcm; }


    real get_dt_gcm() const { return this->dt_gcm; }


    real get_xlen() const { return this->xlen; }


    real get_ylen() const { return this->ylen; }


    int get_nx() const {
      if (dm.find_dimension("x") == -1) return -1;
      return dm.get_dimension_size("x");
    }


    int get_ny() const {
      if (dm.find_dimension("y") == -1) return -1;
      return dm.get_dimension_size("y");
    }


    int get_nz() const {
      if (dm.find_dimension("z") == -1) return -1;
      return dm.get_dimension_size("z");
    }


    int get_nens() const {
      if (dm.find_dimension("nens") == -1) return -1;
      return dm.get_dimension_size("nens");
    }


    int get_ncrms() const {
      return get_nens();
    }


    void add_note( std::string key , std::string value ) {
      notes.add_note(key,value);
    }


    std::string get_note( std::string key ) const {
      return notes.get_note(key);
    }


    bool note_exists( std::string key ) const {
      return notes.note_exists(key);
    }


    void delete_note( std::string key ) {
      notes.delete_note(key);
    }


    inline void set_phys_constants(real R_d, real R_v, real cp_d, real cp_v, real grav=9.81, real p0=1.e5) {
      this->R_d  = R_d ;
      this->R_v  = R_v ;
      this->cp_d = cp_d;
      this->cp_v = cp_v;
      this->grav = grav;
      this->p0   = p0  ;
    }


    inline void set_grid(real xlen, real ylen, real2d const &zint_in) {
      int nz    = get_nz();
      int nens  = get_nens();
      this->xlen = xlen;
      this->ylen = ylen;
      auto zint = dm.get<real,2>("vertical_interface_height");
      auto zmid = dm.get<real,2>("vertical_midpoint_height" );
      parallel_for( "vert grid 1" , SimpleBounds<2>(nz+1,nens) , YAKL_LAMBDA (int k, int iens) {
        zint(k,iens) = zint_in(k,iens);
        if (k < nz) zmid(k,iens) = 0.5_fp * (zint_in(k,iens) + zint_in(k+1,iens));
      });
    }


    inline void set_grid(real xlen, real ylen, real1d const &zint_in) {
      int nz    = get_nz();
      int nens  = get_nens();
      this->xlen = xlen;
      this->ylen = ylen;
      auto zint = dm.get<real,2>("vertical_interface_height");
      auto zmid = dm.get<real,2>("vertical_midpoint_height" );
      parallel_for( "vert grid 2" , SimpleBounds<2>(nz+1,nens) , YAKL_LAMBDA (int k, int iens) {
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
      dm.register_and_allocate<real>( "vertical_interface_height" , "vertical interface height" , {nz+1    ,nens} , {"zp1"      ,"nens"} );
      dm.register_and_allocate<real>( "vertical_midpoint_height"  , "vertical midpoint height"  , {nz      ,nens} , {"z"        ,"nens"} );
      dm.register_and_allocate<real>( "hydrostasis_parameters"    , "hydrostasis parameters"    , {10      ,nens} , {"ten"      ,"nens"} );
      dm.register_and_allocate<real>( "hydrostatic_pressure"      , "hydrostasis pressure"      , {nz      ,nens} , {"z"        ,"nens"} );
      dm.register_and_allocate<real>( "hydrostatic_density"       , "hydrostasis density"       , {nz      ,nens} , {"z"        ,"nens"} );

      auto density_dry  = dm.get_collapsed<real>("density_dry"              );
      auto uvel         = dm.get_collapsed<real>("uvel"                     );
      auto vvel         = dm.get_collapsed<real>("vvel"                     );
      auto wvel         = dm.get_collapsed<real>("wvel"                     );
      auto temp         = dm.get_collapsed<real>("temp"                     );
      auto zint         = dm.get_collapsed<real>("vertical_interface_height");
      auto zmid         = dm.get_collapsed<real>("vertical_midpoint_height" );
      auto hy_params    = dm.get_collapsed<real>("hydrostasis_parameters"   );
      auto hy_press     = dm.get_collapsed<real>("hydrostatic_pressure"     );
      auto hy_dens      = dm.get_collapsed<real>("hydrostatic_density"      );

      parallel_for( "coupler zero" , SimpleBounds<1>(nz*ny*nx*nens) , YAKL_LAMBDA (int i) {
        density_dry (i) = 0;
        uvel        (i) = 0;
        vvel        (i) = 0;
        wvel        (i) = 0;
        temp        (i) = 0;
        if (i < (nz+1)*nens) zint(i) = 0;
        if (i < (nz  )*nens) {
          zmid    (i) = 0;
          hy_press(i) = 0;
          hy_dens (i) = 0;
        }
        if (i < 5     *nens) hy_params(i) = 0;
      });
    }



    inline void update_hydrostasis( realConst4d pressure ) {
      using yakl::intrinsics::matmul_cr;
      using yakl::intrinsics::matinv_ge;

      auto zint      = dm.get<real,2>("vertical_interface_height");
      auto zmid      = dm.get<real,2>("vertical_midpoint_height" );
      auto hy_params = dm.get<real,2>("hydrostasis_parameters"   );
      auto hy_press  = dm.get<real,2>("hydrostatic_pressure"     );
      auto hy_dens   = dm.get<real,2>("hydrostatic_density"      );

      int nz   = get_nz();
      int ny   = get_ny();
      int nx   = get_nx();
      int nens = get_nens();

      int constexpr npts = 10;
      int constexpr npts_tanh = npts - 2;
      SArray<real,1,npts_tanh> tanh_pts;
      real constexpr mu = 1;  // [0,infty) : larger mu means sampling is more clustered toward the domain edges
      for (int ii=0; ii < npts_tanh; ii++) {
        tanh_pts(ii) = (tanh((2*mu*ii)/(npts_tanh-1)-mu)/tanh(mu)+1)*0.5_fp;
      }

      SArray<int,1,npts> z_indices;
      z_indices(0) = 0;
      z_indices(1) = 1;
      for (int ii=1; ii < npts_tanh; ii++) {
        z_indices(ii+1) = tanh_pts(ii) * (nz-1);
      }
      z_indices(npts-2) = nz-2;
      z_indices(npts-1) = nz-1;

      if (z_indices(2)      <= 1   ) z_indices(2     ) = 2;
      if (z_indices(npts-3) >= nz-2) z_indices(npts-3) = nz-3;

      parallel_for( "hydro fitting" , SimpleBounds<1>(nens) , YAKL_LAMBDA (int iens) {
        real zbot = zint(0 ,iens);
        real ztop = zint(nz,iens);

        SArray<double,1,npts> z;
        for (int ii=0; ii < npts; ii++) {
          z(ii) = ( zmid(z_indices(ii),iens) - zbot ) / (ztop - zbot);;
        }

        SArray<double,2,npts,npts> vand;
        for (int j=0; j < npts; j++) {
          for (int i=0; i < npts; i++) {
            vand(j,i) = pow( z(i) , (double) j );
          }
        }

        // TODO: Implement partial pivoting to improve the condition number here
        auto vand_inv = matinv_ge( vand );

        // Fit to just one column, assuming all columns are fairly similar
        // This will only be used for idealized test cases anyway
        // Another function that passes in an ensemble of pressure profiles will be
        //   used for GCM coupling in an MMF setting
        SArray<double,1,npts> logp;
        for (int ii=0; ii < npts; ii++) {
          logp(ii) = log(pressure(z_indices(ii),0,0,iens));
        }

        auto params = matmul_cr( vand_inv , logp );

        for (int i=0; i < npts; i++) { hy_params(i,iens) = params(i); }
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



    real4d compute_pressure_array() const {
      auto dens_dry = dm.get<real const,4>("density_dry");
      auto dens_wv  = dm.get<real const,4>("water_vapor");
      auto temp     = dm.get<real const,4>("temp");

      int nz   = get_nz();
      int ny   = get_ny();
      int nx   = get_nx();
      int nens = get_nens();

      real4d pressure("pressure",nz,ny,nx,nens);

      YAKL_SCOPE( R_d , this->R_d );
      YAKL_SCOPE( R_v , this->R_v );

      parallel_for( "coupler pressure" , SimpleBounds<4>(nz,ny,nx,nens) , YAKL_LAMBDA (int k, int j, int i, int iens) {
        real rho_d = dens_dry(k,j,i,iens);
        real rho_v = dens_wv (k,j,i,iens);
        real T     = temp    (k,j,i,iens);
        pressure(k,j,i,iens) = compute_pressure( rho_d , rho_v , T , R_d , R_v );
      });

      return pressure;
    }



    real4d interp_pressure_interfaces( realConst4d press ) const {
      auto zint      = dm.get<real const,2>("vertical_interface_height");
      auto hy_press  = dm.get<real const,2>("hydrostatic_pressure");
      auto hy_params = dm.get<real const,2>("hydrostasis_parameters");

      int nz   = get_nz();
      int ny   = get_ny();
      int nx   = get_nx();
      int nens = get_nens();

      real4d press_pert("press_pert",nz,ny,nx,nens);

      // Compute pressure perturbation
      parallel_for( "coup press pert" , SimpleBounds<4>(nz,ny,nx,nens) , YAKL_LAMBDA (int k, int j, int i, int iens) {
        press_pert(k,j,i,iens) = press(k,j,i,iens) - hy_press(k,iens);
      });

      // Interpolate pressure perturbation from cells to edges
      VerticalInterp<pam_ord> vert_interp;
      vert_interp.init(zint);
      auto press_edges = vert_interp.cells_to_edges( press_pert ,
                                                     vert_interp.BC_ZERO_GRADIENT ,
                                                     vert_interp.BC_ZERO_GRADIENT );

      // Add hydrostasis at cell edges to get back full pressure
      parallel_for( "coup press edges" , SimpleBounds<4>(nz+1,ny,nx,nens) , YAKL_LAMBDA (int k, int j, int i, int iens) {
        press_edges(k,j,i,iens) += hydrostatic_pressure( hy_params , zint(k,iens) , zint(0,iens) , zint(nz,iens) , iens );
      });
      
      return press_edges;
    }



    real4d interp_density_interfaces( realConst4d dens ) const {
      auto zint      = dm.get<real const,2>("vertical_interface_height");
      auto hy_dens   = dm.get<real const,2>("hydrostatic_density");
      auto hy_params = dm.get<real const,2>("hydrostasis_parameters");

      int nz   = get_nz();
      int ny   = get_ny();
      int nx   = get_nx();
      int nens = get_nens();

      real4d dens_pert("dens_pert",nz,ny,nx,nens);

      YAKL_SCOPE( grav , this->grav );

      // Compute density perturbation
      parallel_for( "coup dens pert" , SimpleBounds<4>(nz,ny,nx,nens) , YAKL_LAMBDA (int k, int j, int i, int iens) {
        dens_pert(k,j,i,iens) = dens(k,j,i,iens) - hy_dens(k,iens);
      });

      // Interpolate density perturbation from cells to edges
      VerticalInterp<pam_ord> vert_interp;
      vert_interp.init(zint);
      auto dens_edges = vert_interp.cells_to_edges( dens_pert ,
                                                    vert_interp.BC_ZERO_GRADIENT ,
                                                    vert_interp.BC_ZERO_GRADIENT );

      // Add hydrostasis at cell edges to get back full density
      parallel_for( "coup dens edges" , SimpleBounds<4>(nz+1,ny,nx,nens) , YAKL_LAMBDA (int k, int j, int i, int iens) {
        dens_edges(k,j,i,iens) += hydrostatic_density( hy_params , zint(k,iens) , zint(0,iens) , zint(nz,iens) , iens , grav );
      });
      
      return dens_edges;
    }



    template <class FP> YAKL_INLINE static void get_gll_points(SArray<FP,1,9> &rslt) {
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



    template <class FP> YAKL_INLINE static void get_gll_weights(SArray<FP,1,9> &rslt) {
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


