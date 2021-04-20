
#pragma once

#include "awfl_const.h"
#include "DataManager.h"

class Microphysics {
public:
  int static constexpr num_tracers = 3;

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

  // TODO: Change this to type int instead of real
  SArray<real,1,num_tracers> tracer_IDs; // tracer index for microphysics tracers

  int static constexpr ID_V = 0;  // Local index for water vapor
  int static constexpr ID_C = 1;  // Local index for cloud liquid
  int static constexpr ID_R = 2;  // Local index for precipitated liquid (rain)



  Microphysics() {
    constants.R_d         = 287.;
    constants.cp_d        = 1003.;
    constants.cv_d        = constants.cp_d - constants.R_d;
    constants.gamma_d     = constants.cp_d / constants.cv_d;
    constants.kappa_d     = constants.R_d  / constants.cp_d;
    constants.R_v         = 461.;
    constants.cp_v        = 1859;
    constants.cv_v        = constants.R_v - constants.cp_v;
    constants.p0          = 1.e5;
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
    //                                        name              description       positive   adds mass
    tracer_IDs(ID_V) = dycore.add_tracer(dm , "water_vapor"   , "Water Vapor"   , true     , true);
    tracer_IDs(ID_C) = dycore.add_tracer(dm , "cloud_liquid"  , "Cloud liquid"  , true     , true);
    tracer_IDs(ID_R) = dycore.add_tracer(dm , "precip_liquid" , "precip_liquid" , true     , true);

    // Register and allocate the tracers in the DataManager
    dm.register_and_allocate<real>( "water_vapor"   , "Water Vapor"   , {nz,ny,nx,nens} , {"z","y","x","nens"} );
    dm.register_and_allocate<real>( "cloud_liquid"  , "Cloud liquid"  , {nz,ny,nx,nens} , {"z","y","x","nens"} );
    dm.register_and_allocate<real>( "precip_liquid" , "precip_liquid" , {nz,ny,nx,nens} , {"z","y","x","nens"} );

    tracer_index_vapor = tracer_IDs(ID_V);

    // Register and allocation non-tracer quantities used by the microphysics
    dm.register_and_allocate<real>( "precl" , "precipitation rate" , {ny,nx,nens} , {"y","x","nens"} );
  }



  void timeStep( DataManager &dm , real dt ) {
    auto rho_v        = dm.get_lev_col<real>("water_vapor");
    auto rho_c        = dm.get_lev_col<real>("cloud_liquid");
    auto rho_r        = dm.get_lev_col<real>("precip_liquid");
    auto rho_dry      = dm.get_lev_col<real>("density_dry");
    auto temp         = dm.get_lev_col<real>("temp");
    auto pressure_dry = dm.get_lev_col<real>("pressure_dry");

    int nz   = dm.get_dimension_size("z"   );
    int ny   = dm.get_dimension_size("y"   );
    int nx   = dm.get_dimension_size("x"   );
    int nens = dm.get_dimension_size("nens");
    int ncol = ny*nx*nens;

    // These are inputs to kessler(...)
    real2d qv          ("qv"          ,nz,ncol);
    real2d qc          ("qc"          ,nz,ncol);
    real2d qr          ("qr"          ,nz,ncol);
    real2d theta_dry   ("theta_dry"   ,nz,ncol);
    real2d exner_dry   ("exner_dry"   ,nz,ncol);
    real2d zmid_in = dm.get<real,2>("vertical_midpoint_height");

    real2d zmid("zmid",nz,ny*nx*nens);
    parallel_for( Bounds<4>(nz,ny,nx,nens) , YAKL_LAMBDA (int k, int j, int i, int iens) {
      zmid(k,j*nx*nens + i*nens + iens) = zmid_in(k,iens);
    });

    // Force constants into local scope
    real gamma_d = this->constants.gamma_d;
    real R_d     = this->constants.R_d;
    real R_v     = this->constants.R_v;
    real cp_d    = this->constants.cp_d;
    real p0      = this->constants.p0;

    // Save initial state, and compute inputs for kessler(...)
    parallel_for( Bounds<2>(nz,ncol) , YAKL_LAMBDA (int k, int i) {
      if (rho_v(k,i) < 0) rho_v(k,i) = 0;
      if (rho_c(k,i) < 0) rho_c(k,i) = 0;
      if (rho_r(k,i) < 0) rho_r(k,i) = 0;
      qv          (k,i) = rho_v(k,i) / rho_dry(k,i);
      qc          (k,i) = rho_c(k,i) / rho_dry(k,i);
      qr          (k,i) = rho_r(k,i) / rho_dry(k,i);
      exner_dry   (k,i) = pow( pressure_dry(k,i) / p0 , R_d / cp_d );
      theta_dry   (k,i) = temp(k,i) / exner_dry(k,i);
    });

    auto precl = dm.get_collapsed<real>("precl");

    kessler(theta_dry, qv, qc, qr, rho_dry, precl, zmid, exner_dry, dt, R_d, cp_d, p0);

    parallel_for( Bounds<2>(nz,ncol) , YAKL_LAMBDA (int k, int i) {
      rho_v(k,i) = qv(k,i)*rho_dry(k,i);
      rho_c(k,i) = qc(k,i)*rho_dry(k,i);
      rho_r(k,i) = qr(k,i)*rho_dry(k,i);
      temp (k,i) = theta_dry(k,i) * exner_dry(k,i);
    });

  }



  ///////////////////////////////////////////////////////////////////////////////
  //
  //  Version:  2.0
  //
  //  Date:  January 22nd, 2015
  //
  //  Change log:
  //  v2 - Added sub-cycling of rain sedimentation so as not to violate
  //       CFL condition.
  //
  //  The KESSLER subroutine implements the Kessler (1969) microphysics
  //  parameterization as described by Soong and Ogura (1973) and Klemp
  //  and Wilhelmson (1978, KW). KESSLER is called at the end of each
  //  time step and makes the final adjustments to the potential
  //  temperature and moisture variables due to microphysical processes
  //  occurring during that time step. KESSLER is called once for each
  //  vertical column of grid cells. Increments are computed and added
  //  into the respective variables. The Kessler scheme contains three
  //  moisture categories: water vapor, cloud water (liquid water that
  //  moves with the flow), and rain water (liquid water that falls
  //  relative to the surrounding air). There  are no ice categories.
  //  
  //  Variables in the column are ordered from the surface to the top.
  //
  //  Parameters:
  //     theta (inout) - dry potential temperature (K)
  //     qv    (inout) - water vapor mixing ratio (gm/gm) (dry mixing ratio)
  //     qc    (inout) - cloud water mixing ratio (gm/gm) (dry mixing ratio)
  //     qr    (inout) - rain  water mixing ratio (gm/gm) (dry mixing ratio)
  //     rho   (in   ) - dry air density (not mean state as in KW) (kg/m^3)
  //     pk    (in   ) - Exner function  (not mean state as in KW) (p/p0)**(R/cp)
  //     dt    (in   ) - time step (s)
  //     z     (in   ) - heights of thermodynamic levels in the grid column (m)
  //     precl (  out) - Precipitation rate (m_water/s)
  //     Rd    (in   ) - Dry air ideal gas constant
  //     cp    (in   ) - Specific heat of dry air at constant pressure
  //     p0    (in   ) - Reference pressure (Pa)
  //
  // Output variables:
  //     Increments are added into t, qv, qc, qr, and precl which are
  //     returned to the routine from which KESSLER was called. To obtain
  //     the total precip qt, after calling the KESSLER routine, compute:
  //
  //       qt = sum over surface grid cells of (precl * cell area)  (kg)
  //       [here, the conversion to kg uses (10^3 kg/m^3)*(10^-3 m/mm) = 1]
  //
  //
  //  Written in Fortran by: Paul Ullrich
  //                         University of California, Davis
  //                         Email: paullrich@ucdavis.edu
  //
  //  Ported to C++ / YAKL by: Matt Norman
  //                           Oak Ridge National Laboratory
  //                           normanmr@ornl.gov
  //                           https://mrnorman.github.io
  //
  //  Based on a code by Joseph Klemp
  //  (National Center for Atmospheric Research)
  //
  //  Reference:
  //
  //    Klemp, J. B., W. C. Skamarock, W. C., and S.-H. Park, 2015:
  //    Idealized Global Nonhydrostatic Atmospheric Test Cases on a Reduced
  //    Radius Sphere. Journal of Advances in Modeling Earth Systems. 
  //    doi:10.1002/2015MS000435
  //
  ///////////////////////////////////////////////////////////////////////////////

  void kessler(real2d &theta, real2d &qv, real2d &qc, real2d &qr, real2d const &rho,
               real1d &precl, real2d const &z, real2d const &pk, real dt,
               real Rd, real cp, real p0) {
    int nz = theta.dimension[0];
    int ncol = theta.dimension[1];

    // Maximum time step size in accordance with CFL condition
    if (dt <= 0) { endrun("kessler.f90 called with nonpositive dt"); }

    real psl    = p0 / 100;  //  pressure at sea level (mb)
    real rhoqr  = 1000._fp;  //  density of liquid water (kg/m^3)
    real lv     = 2.5e6_fp;  //  latent heat of vaporization (J/kg)

    real2d r    ("r"    ,nz  ,ncol);
    real2d rhalf("rhalf",nz  ,ncol);
    real2d pc   ("pc"   ,nz  ,ncol);
    real2d velqr("velqr",nz  ,ncol);
    real2d dt2d ("dt2d" ,nz-1,ncol);

    parallel_for( Bounds<2>(nz,ncol) , YAKL_LAMBDA (int k, int i) {
      r    (k,i) = 0.001_fp * rho(k,i);
      rhalf(k,i) = sqrt( rho(0,i) / rho(k,i) );
      pc   (k,i) = 3.8_fp / ( pow( pk(k,i) , cp/Rd ) * psl );
      // Liquid water terminal velocity (m/s) following KW eq. 2.15
      velqr(k,i) = 36.34_fp * pow( qr(k,i)*r(k,i) , 0.1364_fp ) * rhalf(k,i);
      // Compute maximum stable time step for each cell
      if (k < nz-1) {
        if (velqr(k,i) > 1.e-10_fp) {
          dt2d(k,i) = 0.8_fp * (z(k+1,i)-z(k,i))/velqr(k,i);
        } else {
          dt2d(k,i) = dt;
        }
      }
      // Initialize precip rate to zero
      if (k == 0) {
        precl(i) = 0;
      }
    });

    // Reduce down the minimum time step among the cells
    real dt_max = yakl::intrinsics::minval(dt2d);

    // Number of subcycles
    int rainsplit = ceil(dt / dt_max);
    real dt0 = dt / static_cast<real>(rainsplit);

    real2d sed("sed",nz,ncol);

    // Subcycle through rain process
    for (int nt=0; nt < rainsplit; nt++) {

      // Sedimentation term using upstream differencing
      parallel_for( Bounds<2>(nz,ncol) , YAKL_LAMBDA (int k, int i) {
        if (k == 0) {
          // Precipitation rate (m/s)
          precl(i) = precl(i) + rho(0,i) * qr(0,i) * velqr(0,i) / rhoqr;
        }
        if (k == nz-1) {
          sed(nz-1,i) = -dt0*qr(nz-1,i)*velqr(nz-1,i)/(0.5_fp * (z(nz-1,i)-z(nz-2,i)));
        } else {
          sed(k,i) = dt0 * ( r(k+1,i)*qr(k+1,i)*velqr(k+1,i) - 
                             r(k  ,i)*qr(k  ,i)*velqr(k  ,i) ) / ( r(k,i)*(z(k+1,i)-z(k,i)) );
        }
      });

      // Adjustment terms
      parallel_for( Bounds<2>(nz,ncol) , YAKL_LAMBDA (int k, int i) {
        // Autoconversion and accretion rates following KW eq. 2.13a,b
        real qrprod = qc(k,i) - ( qc(k,i)-dt0*max( 0.001_fp * (qc(k,i)-0.001_fp) , 0._fp ) ) /
                                ( 1 + dt0 * 2.2_fp * pow( qr(k,i) , 0.875_fp ) );
        qc(k,i) = max( qc(k,i)-qrprod , 0._fp );
        qr(k,i) = max( qr(k,i)+qrprod+sed(k,i) , 0._fp );

        // Saturation vapor mixing ratio (gm/gm) following KW eq. 2.11
        real tmp = pk(k,i)*theta(k,i)-36._fp;
        real qvs = pc(k,i)*exp( 17.27_fp * (pk(k,i)*theta(k,i)-273._fp) / tmp );
        real prod = (qv(k,i)-qvs) / (1._fp + qvs*(4093._fp * lv/cp)/(tmp*tmp));

        // Evaporation rate following KW eq. 2.14a,b
        real tmp1 = dt0*( ( ( 1.6_fp + 124.9_fp * pow( r(k,i)*qr(k,i) , 0.2046_fp ) ) *
                            pow( r(k,i)*qr(k,i) , 0.525_fp ) ) /
                          ( 2550000._fp * pc(k,i) / (3.8_fp * qvs)+540000._fp) ) * 
                        ( max(qvs-qv(k,i),0._fp) / (r(k,i)*qvs) );
        real tmp2 = max( -prod-qc(k,i) , 0._fp );
        real tmp3 = qr(k,i);
        real ern = min( tmp1 , min( tmp2 , tmp3 ) );

        // Saturation adjustment following KW eq. 3.10
        theta(k,i)= theta(k,i) + lv / (cp*pk(k,i)) * 
                                 ( max( prod , -qc(k,i) ) - ern );
        qv(k,i) = max( qv(k,i) - max( prod , -qc(k,i) ) + ern , 0._fp );
        qc(k,i) = qc(k,i) + max( prod , -qc(k,i) );
        qr(k,i) = qr(k,i) - ern;

        // Recalculate liquid water terminal velocity
        velqr(k,i)  = 36.34_fp * pow( qr(k,i)*r(k,i) , 0.1364_fp ) * rhalf(k,i);
        if (k == 0 && nt == rainsplit-1) {
          precl(i) = precl(i) / static_cast<real>(rainsplit);
        }
      });

    }

  }



  void output(DataManager &dm, yakl::SimpleNetCDF &nc, int ulIndex, int iens) const {
    auto precl = dm.get<real,3>("precl");
    int nx = dm.get_dimension_size("x");
    int ny = dm.get_dimension_size("y");
    real2d data("data",ny,nx);
    parallel_for( Bounds<2>(ny,nx) , YAKL_LAMBDA (int j, int i) {
      data(j,i) = precl(j,i,iens);
    });
    nc.write1(data.createHostCopy(),"precl",{"y","x"},ulIndex,"t");
  }



  std::string micro_name() const {
    return "kessler";
  }




};



