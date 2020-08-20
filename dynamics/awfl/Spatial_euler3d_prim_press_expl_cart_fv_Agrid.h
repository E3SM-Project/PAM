
#pragma once

#include "const.h"
#include "phys_params.h"
#include "TransformMatrices.h"
#include "WenoLimiter.h"
#include "Profiles.h"


/*********************************************
 ***** Required inside the Spatial class *****
 *********************************************
class Location;
    - Stores the indices of a single location on the grid

class StateArr; // OR
typedef [class] StateArr;
    - Declare a type for the model state

class TendArr; // OR
typedef [class] TendArr;
    - Declare a type for the model tendencies
    - It must have room for the nTimeDerivs dimension for the time integrator

StateArr createStateArr()
    - Create and return a new StateArr object

TendArr createTendArr(int nTimeDerivs)
    - Create and return a new TendArr object

YAKL_INLINE real &get(StateArr const &state , Location const &loc , int splitIndex)
    - Return the state value at the given location

YAKL_INLINE real &get(TendArr const &tend , Location const &loc , int timeDeriv , int splitIndex)
    - Return the tendency value at the given location

int numSplit()
    - Return the number of split components for this operator
    - The temporal operator will iterate through the splittings

real computeTimeStep(real cfl)
    - Return the time step in seconds from the cfl value

void init(int nTimeDerivs, bool timeAvg, std::string inFile)
    - Initialize internal data structures

void initState( StateArr &state )
    - Initialize the state

void computeTendencies( StateArr const &state , TendArr &tend , real dt , int splitIndex)
    - Compute tendency and time derivatives of the tendency if they are requested

template <class F> void applyTendencies( F const &applySingleTendency , int splitIndex )
    - Loop through the domain, and apply tendencies to the state

const char * getSpatialName()
    - Return the name and info for this spatial operator

void output(StateArr const &state, real etime)
    - Output to file
*/


class Spatial {
public:

  int static constexpr hs = (ord-1)/2;
  int static constexpr numState = 6;

  // Stores a single index location
  struct Location {
    int l;
    int k;
    int j;
    int i;
  };

  typedef real4d StateArr;  // Spatial index

  typedef real5d TendArr;   // (time derivative & spatial index)

  real6d fwaves;            // state edge estimates and flux difference splitting waves
  // Hydrostatically balanced values for density, potential temperature, and pressure
  real1d hyDensCells;
  real1d hyPressureCells;
  real2d hyDensGLL;
  real2d hyPressureGLL;
  real2d hyPressureDerivGLL;
  // Matrices to transform DOFs from one form to another
  SArray<real,2,ord,ngll> coefs_to_gll;
  SArray<real,2,ord,ngll> coefs_to_deriv_gll;
  SArray<real,2,ord,ngll> sten_to_gll;
  SArray<real,2,ord,ngll> sten_to_deriv_gll;
  SArray<real,3,ord,ord,ord> wenoRecon;
  SArray<real,1,hs+2> idl;
  real sigma;
  // For ADER spatial derivative computation
  SArray<real,2,ngll,ngll> derivMatrix;
  // For quadrature
  SArray<real,1,ord> gllWts_ord;
  SArray<real,1,ord> gllPts_ord;
  SArray<real,1,ngll> gllWts_ngll;
  SArray<real,1,ngll> gllPts_ngll;

  int static constexpr idR = 0;  // density perturbation
  int static constexpr idU = 1;  // u
  int static constexpr idV = 2;  // v
  int static constexpr idW = 3;  // w
  int static constexpr idT = 4;  // potential temperature
  int static constexpr idP = 5;  // pressure perturbation

  int static constexpr BC_PERIODIC = 0;
  int static constexpr BC_WALL     = 1;

  int static constexpr DATA_SPEC_THERMAL = 1;

  bool sim2d;

  real dx;
  real dy;
  real dz;

  real dtInit;

  bool dimSwitch;

  // Values read from input file
  int         nx;
  int         ny;
  int         nz;
  real        xlen;
  real        ylen;
  real        zlen;
  int         bc_x;
  int         bc_y;
  int         bc_z;
  bool        weno_scalars;
  bool        weno_winds;
  std::string outFile;
  int         dataSpec;


  static_assert(ord%2 == 1,"ERROR: ord must be an odd integer");


  StateArr createStateArr() {
    return StateArr("stateArr",numState,nz+2*hs,ny+2*hs,nx+2*hs);
  }



  TendArr createTendArr() {
    return TendArr("tendArr",numState,nTimeDerivs,nz,ny,nx);
  }



  YAKL_INLINE real &get(StateArr const &state , Location const &loc , int splitIndex) {
    return state(loc.l , hs+loc.k , hs+loc.j , hs+loc.i);
  }



  YAKL_INLINE real &get(TendArr const &tend , Location const &loc , int timeDeriv , int splitIndex) {
    return tend(loc.l , timeDeriv , loc.k , loc.j , loc.i);
  }



  int numSplit() {
    return 4;
  }



  real computeTimeStep(real cfl, StateArr const &state) {
    if (dtInit <= 0) {
      real maxwave = 0;
      real3d dt3d("dt3d",nz,ny,nx);
      parallel_for( Bounds<3>(nz,ny,nx) , YAKL_LAMBDA (int k, int j, int i) {
        real r = state(idR,hs+k,hs+j,hs+i) + hyDensCells(k);
        real u = state(idU,hs+k,hs+j,hs+i);
        real v = state(idV,hs+k,hs+j,hs+i);
        real w = state(idW,hs+k,hs+j,hs+i);
        real t = state(idT,hs+k,hs+j,hs+i);
        real p = state(idP,hs+k,hs+j,hs+i) + hyPressureCells(k);
        real cs = sqrt(GAMMA*p/r);
        real udt = cfl * dx / max( abs(u-cs) , abs(u+cs) );
        real vdt = cfl * dy / max( abs(v-cs) , abs(v+cs) );
        real wdt = cfl * dz / max( abs(w-cs) , abs(w+cs) );
        dt3d(k,j,i) = min( min(udt,vdt) , wdt );
      });
      yakl::ParallelMin<real,memDevice> pmin(nz*nx*ny);
      dtInit = pmin(dt3d.data());
    }
    return dtInit;
  }



  // Initialize crap needed by recon()
  void init(std::string inFile) {
    dtInit = 0;
    dimSwitch = true;

    // Read the input file
    YAML::Node config = YAML::LoadFile(inFile);
    if ( !config                 ) { endrun("ERROR: Invalid YAML input file"); }
    if ( !config["nx"]           ) { endrun("ERROR: No nx in input file"); }
    if ( !config["ny"]           ) { endrun("ERROR: No ny in input file"); }
    if ( !config["nz"]           ) { endrun("ERROR: No nz in input file"); }
    if ( !config["xlen"]         ) { endrun("ERROR: No xlen in input file"); }
    if ( !config["ylen"]         ) { endrun("ERROR: No ylen in input file"); }
    if ( !config["zlen"]         ) { endrun("ERROR: No zlen in input file"); }
    if ( !config["bc_x"]         ) { endrun("ERROR: No bc_x in input file"); }
    if ( !config["bc_y"]         ) { endrun("ERROR: No bc_y in input file"); }
    if ( !config["bc_z"]         ) { endrun("ERROR: No bc_z in input file"); }
    if ( !config["weno_scalars"] ) { endrun("ERROR: No weno_scalars in input file"); }
    if ( !config["weno_winds"]   ) { endrun("ERROR: No weno_winds in input file"); }
    if ( !config["initData"]     ) { endrun("ERROR: No initData in input file"); }
    if ( !config["outFile"]      ) { endrun("ERROR: No outFile in input file"); }

    nx = config["nx"].as<int>();
    ny = config["ny"].as<int>();
    nz = config["nz"].as<int>();

    sim2d = ny == 1;

    xlen = config["xlen"].as<real>();
    ylen = config["ylen"].as<real>();
    zlen = config["zlen"].as<real>();

    weno_scalars = config["weno_scalars"].as<bool>();
    weno_winds   = config["weno_winds"].as<bool>();

    std::string dataStr = config["initData"].as<std::string>();
    if        (dataStr == "thermal") {
      dataSpec = DATA_SPEC_THERMAL;
    } else {
      endrun("ERROR: Invalid dataSpec");
    }

    std::string bc_x_str = config["bc_x"].as<std::string>();
    if        (bc_x_str == "periodic" ) {
      bc_x = BC_PERIODIC;
    } else if (bc_x_str == "wall"     ) {
      bc_x = BC_WALL;
    } else {
      endrun("Invalid bc_x");
    }

    std::string bc_y_str = config["bc_y"].as<std::string>();
    if        (bc_y_str == "periodic" ) {
      bc_y = BC_PERIODIC;
    } else if (bc_y_str == "wall"     ) {
      bc_y = BC_WALL;
    } else {
      endrun("Invalid bc_y");
    }

    std::string bc_z_str = config["bc_z"].as<std::string>();
    if        (bc_z_str == "periodic" ) {
      bc_z = BC_PERIODIC;
    } else if (bc_z_str == "wall"     ) {
      bc_z = BC_WALL;
    } else {
      endrun("Invalid bc_z");
    }

    outFile = config["outFile"].as<std::string>();

    dx = xlen/nx;
    dy = ylen/ny;
    dz = zlen/nz;

    TransformMatrices::weno_sten_to_coefs(this->wenoRecon);

    // Store to_gll and wenoRecon
    {
      SArray<real,2,ord,ord>  g2c;
      SArray<real,2,ord,ord>  s2c;
      SArray<real,2,ord,ngll> c2g_lower;
      SArray<real,2,ord,ord>  c2g;
      SArray<real,2,ord,ord>  c2d;

      TransformMatrices::gll_to_coefs      (g2c      );
      TransformMatrices::sten_to_coefs     (s2c      );
      TransformMatrices::coefs_to_gll_lower(c2g_lower);
      TransformMatrices::coefs_to_gll      (c2g      );
      TransformMatrices::coefs_to_deriv    (c2d      );

      this->coefs_to_gll       = c2g_lower;
      this->coefs_to_deriv_gll = c2g_lower * c2d;
      this->sten_to_gll        = c2g_lower       * s2c;
      this->sten_to_deriv_gll  = c2g_lower * c2d * s2c;

    }
    // Store ader derivMatrix
    {
      SArray<real,2,ngll,ngll> g2c;
      SArray<real,2,ngll,ngll> c2d;
      SArray<real,2,ngll,ngll> c2g;

      TransformMatrices::gll_to_coefs  (g2c);
      TransformMatrices::coefs_to_deriv(c2d);
      TransformMatrices::coefs_to_gll  (c2g);

      this->derivMatrix = c2g * c2d * g2c;
    }
    TransformMatrices::get_gll_points (this->gllPts_ord);
    TransformMatrices::get_gll_weights(this->gllWts_ord);
    TransformMatrices::get_gll_points (this->gllPts_ngll);
    TransformMatrices::get_gll_weights(this->gllWts_ngll);

    weno::wenoSetIdealSigma(this->idl,this->sigma);

    fwaves = real6d("fwaves",numState,2,nTimeDerivs,nz+1,ny+1,nx+1);
    hyDensCells        = real1d("hyDensCells       ",nz     );
    hyPressureCells    = real1d("hyPressureCells   ",nz     );
    hyDensGLL          = real2d("hyDensGLL         ",nz,ngll);
    hyPressureGLL      = real2d("hyPressureGLL     ",nz,ngll);
    hyPressureDerivGLL = real2d("hyPressureDerivGLL",nz,ngll);
  }



  // Initialize the state
  void initState( StateArr &state ) {
    // Setup hydrostatic background state
    parallel_for( Bounds<1>(nz) , YAKL_LAMBDA (int k) {
      // Compute cell averages
      hyDensCells    (k) = 0;
      hyPressureCells(k) = 0;
      for (int kk=0; kk<ord; kk++) {
        real zloc = (k+0.5_fp)*dz + gllPts_ord(kk)*dz;
        if        (dataSpec == DATA_SPEC_THERMAL) {
          // Compute constant theta hydrostatic background state
          real t  = 300;
          real wt = gllWts_ord(kk);
          hyDensCells    (k) += profiles::initConstTheta_density (t,zloc) * wt;
          hyPressureCells(k) += profiles::initConstTheta_pressure(t,zloc) * wt;
        }
      }
      // Compute ngll GLL points
      for (int kk=0; kk<ngll; kk++) {
        real zloc = (k+0.5_fp)*dz + gllPts_ngll(kk)*dz;
        if        (dataSpec == DATA_SPEC_THERMAL) {
          // Compute constant theta hydrostatic background state
          real t = 300;
          hyDensGLL         (k,kk) = profiles::initConstTheta_density      (t,zloc);
          hyPressureGLL     (k,kk) = profiles::initConstTheta_pressure     (t,zloc);
          hyPressureDerivGLL(k,kk) = profiles::initConstTheta_pressureDeriv(t,zloc);
        }
      }
    });

    // Compute the state
    parallel_for( Bounds<3>(nz,ny,nx) , YAKL_LAMBDA (int k, int j, int i) {
      for (int l=0; l < numState; l++) {
        state(l,hs+k,hs+j,hs+i) = 0;
      }
      for (int kk=0; kk<ord; kk++) {
        for (int jj=0; jj<ord; jj++) {
          for (int ii=0; ii<ord; ii++) {
            real zloc = (k+0.5_fp)*dz + gllPts_ord(kk)*dz;
            real yloc;
            if (sim2d) {
              yloc = ylen/2;
            } else {
              yloc = (j+0.5_fp)*dy + gllPts_ord(jj)*dy;
            }
            real xloc = (i+0.5_fp)*dx + gllPts_ord(ii)*dx;

            real wt = gllWts_ord(kk) * gllWts_ord(jj) * gllWts_ord(ii);

            if        (dataSpec == DATA_SPEC_THERMAL) {
              real constexpr t0 = 300;
              real rh = profiles::initConstTheta_density (t0,zloc);
              real ph = profiles::initConstTheta_pressure(t0,zloc);

              // Compute constant theta hydrostatic background state
              real tp = profiles::ellipsoid_linear(xloc, yloc, zloc, xlen/2, ylen/2, 2000, 2000, 2000, 2000, 2 );
              real t  = t0+tp;
              real pp = C0*pow(rh*t,GAMMA) - ph;
              state(idT,hs+k,hs+j,hs+i) += t  * wt;
              state(idP,hs+k,hs+j,hs+i) += pp * wt;
            }
          }
        }
      }
    });
  }



  // Compute state and tendency time derivatives from the state
  void computeTendencies( StateArr &state , TendArr &tend , real dt , int splitIndex ) {
    if (dimSwitch) {
      if        (splitIndex == 0) {
        computeTendenciesX( state , tend , dt );
      } else if (splitIndex == 1) {
        if (sim2d) {
          memset(tend,0._fp);
        } else {
          computeTendenciesY( state , tend , dt );
        }
      } else if (splitIndex == 2) {
        computeTendenciesZ( state , tend , dt );
      } else if (splitIndex == 3) {
        computeTendenciesS( state , tend , dt );
      }
    } else {
      if        (splitIndex == 0) {
        computeTendenciesS( state , tend , dt );
      } else if (splitIndex == 1) {
        computeTendenciesZ( state , tend , dt );
      } else if (splitIndex == 2) {
        if (sim2d) {
          memset(tend,0._fp);
        } else {
          computeTendenciesY( state , tend , dt );
        }
      } else if (splitIndex == 3) {
        computeTendenciesX( state , tend , dt );
      }
    }
    if (splitIndex == numSplit()) dimSwitch = ! dimSwitch;
  } // computeTendencies



  void computeTendenciesX( StateArr &state , TendArr &tend , real dt ) {
    // Populate the halos
    if        (bc_x == BC_PERIODIC) {
      parallel_for( Bounds<3>(nz,ny,hs) , YAKL_LAMBDA(int k, int j, int ii) {
        for (int l=0; l < numState; l++) {
          state(l,hs+k,hs+j,      ii) = state(l,hs+k,hs+j,nx+ii);
          state(l,hs+k,hs+j,hs+nx+ii) = state(l,hs+k,hs+j,hs+ii);
        }
      });
    } else if (bc_x == BC_WALL) {
      parallel_for( Bounds<3>(nz,ny,hs) , YAKL_LAMBDA(int k, int j, int ii) {
        for (int l=0; l < numState; l++) {
          if (l == idU) {
            state(l,hs+k,hs+j,      ii) = 0;
            state(l,hs+k,hs+j,hs+nx+ii) = 0;
          } else {
            state(l,hs+k,hs+j,      ii) = state(l,hs+k,hs+j,hs     );
            state(l,hs+k,hs+j,hs+nx+ii) = state(l,hs+k,hs+j,hs+nx-1);
          }
        }
      });
    }

    // Loop through all cells, reconstruct in x-direction, compute centered tendencies, store cell-edge state estimates
    parallel_for( Bounds<3>(nz,ny,nx) , YAKL_LAMBDA (int k, int j, int i) {
      ////////////////////////////////////////////////////////////////
      // Reconstruct rho, u, v, w, theta, pressure, dudx, and dpdx
      ////////////////////////////////////////////////////////////////
      SArray<real,2,nAder,ngll> r_DTs;
      SArray<real,2,nAder,ngll> u_DTs;
      SArray<real,2,nAder,ngll> du_DTs;
      SArray<real,2,nAder,ngll> v_DTs;
      SArray<real,2,nAder,ngll> dv_DTs;
      SArray<real,2,nAder,ngll> w_DTs;
      SArray<real,2,nAder,ngll> dw_DTs;
      SArray<real,2,nAder,ngll> t_DTs;
      SArray<real,2,nAder,ngll> dt_DTs;
      SArray<real,2,nAder,ngll> p_DTs;
      SArray<real,2,nAder,ngll> dp_DTs;
      {
        SArray<real,1,ord> stencil;

        // Density
        for (int ii=0; ii < ord; ii++) { stencil(ii) = state(idR,hs+k,hs+j,i+ii); }
        reconstruct_gll_values( stencil , r_DTs , weno_scalars );
        for (int ii=0; ii < ngll; ii++) { r_DTs(0,ii) += hyDensCells(k); } // Add hydrostasis back on

        // u values and derivatives
        for (int ii=0; ii < ord; ii++) { stencil(ii) = state(idU,hs+k,hs+j,i+ii); }
        reconstruct_gll_values_and_derivs( stencil , u_DTs , du_DTs , dx , weno_winds );

        // v
        for (int ii=0; ii < ord; ii++) { stencil(ii) = state(idV,hs+k,hs+j,i+ii); }
        reconstruct_gll_values_and_derivs( stencil , v_DTs , dv_DTs , dx , weno_winds );

        // w
        for (int ii=0; ii < ord; ii++) { stencil(ii) = state(idW,hs+k,hs+j,i+ii); }
        reconstruct_gll_values_and_derivs( stencil , w_DTs , dw_DTs , dx , weno_winds );

        // theta
        for (int ii=0; ii < ord; ii++) { stencil(ii) = state(idT,hs+k,hs+j,i+ii); }
        reconstruct_gll_values_and_derivs( stencil , t_DTs , dt_DTs , dx , weno_winds );

        // pressure values and derivatives
        for (int ii=0; ii < ord; ii++) { stencil(ii) = state(idP,hs+k,hs+j,i+ii); }
        reconstruct_gll_values_and_derivs( stencil , p_DTs , dp_DTs , dx , weno_scalars );
        // Add hydrostasis, and then divide by C0 to get (rho*theta)^gamma
        for (int ii=0; ii < ngll; ii++) { p_DTs(0,ii) += hyPressureCells(k); }
      }

      ///////////////////////////////////////////////////////////////
      // Compute other values needed for centered tendencies and DTs
      ///////////////////////////////////////////////////////////////
      SArray<real,2,nAder,ngll> r_u_DTs;
      SArray<real,2,nAder,ngll> u_u_DTs;
      SArray<real,2,nAder,ngll> rr_dp_DTs;
      SArray<real,2,nAder,ngll> u_dv_DTs;
      SArray<real,2,nAder,ngll> u_dw_DTs;
      SArray<real,2,nAder,ngll> u_dt_DTs;
      SArray<real,2,nAder,ngll> u_dp_DTs;
      SArray<real,2,nAder,ngll> p_du_DTs;
      for (int ii=0; ii < ngll; ii++) {
        real r  = r_DTs (0,ii);
        real u  = u_DTs (0,ii);
        real v  = v_DTs (0,ii);
        real w  = w_DTs (0,ii);
        real t  = t_DTs (0,ii);
        real p  = p_DTs (0,ii);
        real du = du_DTs(0,ii);
        real dv = dv_DTs(0,ii);
        real dw = dw_DTs(0,ii);
        real dt = dt_DTs(0,ii);
        real dp = dp_DTs(0,ii);
        r_u_DTs  (0,ii) = r*u;
        u_u_DTs  (0,ii) = u*u;
        rr_dp_DTs(0,ii) = dp/r;
        u_dv_DTs (0,ii) = u*dv;
        u_dw_DTs (0,ii) = u*dw;
        u_dt_DTs (0,ii) = u*dt;
        u_dp_DTs (0,ii) = u*dp;
        p_du_DTs (0,ii) = p*du;
      }

      //////////////////////////////////////////
      // Compute time derivatives if necessary
      //////////////////////////////////////////
      if (nAder > 1) {
        for (int kt=0; kt < nAder-1; kt++) {
          // Compute state at kt+1
          for (int ii=0; ii<ngll; ii++) {
            // Compute d_dx(r*u) and d_dx(u*u)
            real dru_dx  = 0;
            real duu_dx  = 0;
            for (int s=0; s<ngll; s++) {
              dru_dx += derivMatrix(s,ii) * r_u_DTs(kt,s);
              duu_dx += derivMatrix(s,ii) * u_u_DTs(kt,s);
            }
            dru_dx /= dx;
            duu_dx /= dx;
            // Compute state at kt+1
            r_DTs(kt+1,ii) = -( dru_dx                                  ) / (kt+1);
            u_DTs(kt+1,ii) = -( duu_dx/2 + rr_dp_DTs(kt,ii)             ) / (kt+1);
            v_DTs(kt+1,ii) = -( u_dv_DTs(kt,ii)                         ) / (kt+1);
            w_DTs(kt+1,ii) = -( u_dw_DTs(kt,ii)                         ) / (kt+1);
            t_DTs(kt+1,ii) = -( u_dt_DTs(kt,ii)                         ) / (kt+1);
            p_DTs(kt+1,ii) = -( u_dp_DTs(kt,ii) + GAMMA*p_du_DTs(kt,ii) ) / (kt+1);
          }
          if (bc_x == BC_WALL) {
            if (i == nx-1) u_DTs(kt+1,ngll-1) = 0;
            if (i == 0   ) u_DTs(kt+1,0     ) = 0;
          }
          // Compute du, dv, dw, dt, dp, r*u, u*u, dp/r, u*dv, u*dw, u*dt, u*dp, p*du
          for (int ii=0; ii<ngll; ii++) {
            // Differentiate du, dv, dw, dt, and dp at kt+1
            real du_dx = 0;
            real dv_dx = 0;
            real dw_dx = 0;
            real dt_dx = 0;
            real dp_dx = 0;
            for (int s=0; s<ngll; s++) {
              du_dx += derivMatrix(s,ii) * u_DTs(kt+1,s);
              dv_dx += derivMatrix(s,ii) * v_DTs(kt+1,s);
              dw_dx += derivMatrix(s,ii) * w_DTs(kt+1,s);
              dt_dx += derivMatrix(s,ii) * t_DTs(kt+1,s);
              dp_dx += derivMatrix(s,ii) * p_DTs(kt+1,s);
            }
            du_DTs(kt+1,ii) = du_dx / dx;
            dv_DTs(kt+1,ii) = dv_dx / dx;
            dw_DTs(kt+1,ii) = dw_dx / dx;
            dt_DTs(kt+1,ii) = dt_dx / dx;
            dp_DTs(kt+1,ii) = dp_dx / dx;
            // Compute r*u, u*u, dp/r, u*dv, u*dw, u*dt, u*dp, p*du
            real tot_r_u   = 0;
            real tot_u_u   = 0;
            real tot_rr_dp = 0;
            real tot_u_dv  = 0;
            real tot_u_dw  = 0;
            real tot_u_dt  = 0;
            real tot_u_dp  = 0;
            real tot_p_du  = 0;
            rr_dp_DTs(kt+1,ii) = 0;
            for (int rt=0; rt <= kt+1; rt++) {
              tot_r_u   += r_DTs(rt,ii) * u_DTs (kt+1-rt,ii);
              tot_u_u   += u_DTs(rt,ii) * u_DTs (kt+1-rt,ii);
              tot_rr_dp += dp_DTs(rt,ii) - r_DTs(rt,ii) * rr_dp_DTs(kt+1-rt,ii);
              tot_u_dv  += u_DTs(rt,ii) * dv_DTs(kt+1-rt,ii);
              tot_u_dw  += u_DTs(rt,ii) * dw_DTs(kt+1-rt,ii);
              tot_u_dt  += u_DTs(rt,ii) * dt_DTs(kt+1-rt,ii);
              tot_u_dp  += u_DTs(rt,ii) * dp_DTs(kt+1-rt,ii);
              tot_p_du  += p_DTs(rt,ii) * du_DTs(kt+1-rt,ii);
            }
            r_u_DTs  (kt+1,ii) = tot_r_u;
            u_u_DTs  (kt+1,ii) = tot_u_u;
            rr_dp_DTs(kt+1,ii) = tot_rr_dp / r_DTs(0,ii);
            u_dv_DTs (kt+1,ii) = tot_u_dv;
            u_dw_DTs (kt+1,ii) = tot_u_dw;
            u_dt_DTs (kt+1,ii) = tot_u_dt;
            u_dp_DTs (kt+1,ii) = tot_u_dp;
            p_du_DTs (kt+1,ii) = tot_p_du;
          }
        }
      }

      //////////////////////////////////////////
      // Time average if necessary
      //////////////////////////////////////////
      if (timeAvg) {
        // Compute time averages
        for (int ii=0; ii<ngll; ii++) {
          real dtmult = 1;
          real r_tavg     = 0;
          real u_tavg     = 0;
          real v_tavg     = 0;
          real w_tavg     = 0;
          real t_tavg     = 0;
          real p_tavg     = 0;
          real r_u_tavg   = 0;
          real u_u_tavg   = 0;
          real rr_dp_tavg = 0;
          real u_dv_tavg  = 0;
          real u_dw_tavg  = 0;
          real u_dt_tavg  = 0;
          real u_dp_tavg  = 0;
          real p_du_tavg  = 0;
          for (int kt=0; kt<nAder; kt++) {
            r_tavg     += r_DTs    (kt,ii) * dtmult / (kt+1);
            u_tavg     += u_DTs    (kt,ii) * dtmult / (kt+1);
            v_tavg     += v_DTs    (kt,ii) * dtmult / (kt+1);
            w_tavg     += w_DTs    (kt,ii) * dtmult / (kt+1);
            t_tavg     += t_DTs    (kt,ii) * dtmult / (kt+1);
            p_tavg     += p_DTs    (kt,ii) * dtmult / (kt+1);
            r_u_tavg   += r_u_DTs  (kt,ii) * dtmult / (kt+1);
            u_u_tavg   += u_u_DTs  (kt,ii) * dtmult / (kt+1);
            rr_dp_tavg += rr_dp_DTs(kt,ii) * dtmult / (kt+1);
            u_dv_tavg  += u_dv_DTs (kt,ii) * dtmult / (kt+1);
            u_dw_tavg  += u_dw_DTs (kt,ii) * dtmult / (kt+1);
            u_dt_tavg  += u_dt_DTs (kt,ii) * dtmult / (kt+1);
            u_dp_tavg  += u_dp_DTs (kt,ii) * dtmult / (kt+1);
            p_du_tavg  += p_du_DTs (kt,ii) * dtmult / (kt+1);
            dtmult *= dt;
          }
          r_DTs    (0,ii) = r_tavg    ;
          u_DTs    (0,ii) = u_tavg    ;
          v_DTs    (0,ii) = v_tavg    ;
          w_DTs    (0,ii) = w_tavg    ;
          t_DTs    (0,ii) = t_tavg    ;
          p_DTs    (0,ii) = p_tavg    ;
          r_u_DTs  (0,ii) = r_u_tavg  ;
          u_u_DTs  (0,ii) = u_u_tavg  ;
          rr_dp_DTs(0,ii) = rr_dp_tavg;
          u_dv_DTs (0,ii) = u_dv_tavg ;
          u_dw_DTs (0,ii) = u_dw_tavg ;
          u_dt_DTs (0,ii) = u_dt_tavg ;
          u_dp_DTs (0,ii) = u_dp_tavg ;
          p_du_DTs (0,ii) = p_du_tavg ;
        }
      }

      //////////////////////////////////////////
      // Compute centered tendencies
      //////////////////////////////////////////
      // Compute flux-form updates
      real drdt = - ( r_u_DTs(0,ngll-1) - r_u_DTs(0,0) ) / dx;
      real dudt = - ( u_u_DTs(0,ngll-1) - u_u_DTs(0,0) ) / dx / 2;
      real dvdt = 0;
      real dwdt = 0;
      real dtdt = 0;
      real dpdt = 0;
      // for (int ii=0; ii < ngll; ii++) {
      //   real wt = gllWts_ngll(ii);
      //   dudt += -rr_dp_DTs(0,ii) * wt;
      //   if (! sim2d) {
      //     dvdt += -u_dv_DTs(0,ii) * wt;
      //   }
      //   dwdt += -u_dw_DTs(0,ii) * wt;
      //   dtdt += -u_dt_DTs(0,ii) * wt;
      //   dpdt += ( -u_dp_DTs(0,ii) - GAMMA*p_du_DTs(0,ii) ) * wt;
      // }
      real r = state(idR,hs+k,hs+j,hs+i) + hyDensCells(k);
      real u = state(idU,hs+k,hs+j,hs+i);
      real p = state(idP,hs+k,hs+j,hs+i) + hyPressureCells(k);
      dudt += - (1./r) * ( (p_DTs(0,ngll-1)-hyPressureCells(k)) - (p_DTs(0,0)-hyPressureCells(k)) ) / dx;
      if (! sim2d) {
        dvdt += - u * ( v_DTs(0,ngll-1) - v_DTs(0,0) ) / dx;
      }
      dwdt += - u * ( w_DTs(0,ngll-1) - w_DTs(0,0) ) / dx;
      dtdt += - u * ( t_DTs(0,ngll-1) - t_DTs(0,0) ) / dx;
      dpdt += - u *       ( p_DTs(0,ngll-1) - p_DTs(0,0) ) / dx
              - GAMMA*p * ( u_DTs(0,ngll-1) - u_DTs(0,0) ) / dx;

      tend(idR,0,k,j,i) = drdt;
      tend(idU,0,k,j,i) = dudt;
      tend(idV,0,k,j,i) = dvdt;
      tend(idW,0,k,j,i) = dwdt;
      tend(idT,0,k,j,i) = dtdt;
      tend(idP,0,k,j,i) = dpdt;

      //////////////////////////////////////////
      // Store cell edge estimates of the state
      //////////////////////////////////////////
      fwaves(idR,1,0,k,j,i  ) = r_DTs(0,0     );
      fwaves(idU,1,0,k,j,i  ) = u_DTs(0,0     );
      fwaves(idV,1,0,k,j,i  ) = v_DTs(0,0     );
      fwaves(idW,1,0,k,j,i  ) = w_DTs(0,0     );
      fwaves(idT,1,0,k,j,i  ) = t_DTs(0,0     );
      fwaves(idP,1,0,k,j,i  ) = p_DTs(0,0     );

      fwaves(idR,0,0,k,j,i+1) = r_DTs(0,ngll-1);
      fwaves(idU,0,0,k,j,i+1) = u_DTs(0,ngll-1);
      fwaves(idV,0,0,k,j,i+1) = v_DTs(0,ngll-1);
      fwaves(idW,0,0,k,j,i+1) = w_DTs(0,ngll-1);
      fwaves(idT,0,0,k,j,i+1) = t_DTs(0,ngll-1);
      fwaves(idP,0,0,k,j,i+1) = p_DTs(0,ngll-1);
    });

    ////////////////////////////////////////////////
    // BCs for the state edge estimates
    ////////////////////////////////////////////////
    parallel_for( Bounds<4>(numState,nTimeDerivs,nz,ny) , YAKL_LAMBDA (int l, int kt, int k, int j) {
      if        (bc_x == BC_PERIODIC) {
        fwaves(l,0,kt,k,j,0 ) = fwaves(l,0,kt,k,j,nx);
        fwaves(l,1,kt,k,j,nx) = fwaves(l,1,kt,k,j,0 );
      } else if (bc_x == BC_WALL    ) {
        if (l == idU) {
          fwaves(l,0,kt,k,j,0 ) = 0;
          fwaves(l,1,kt,k,j,0 ) = 0;
          fwaves(l,0,kt,k,j,nx) = 0;
          fwaves(l,1,kt,k,j,nx) = 0;
        } else {
          fwaves(l,0,kt,k,j,0 ) = fwaves(l,1,kt,k,j,0 );
          fwaves(l,1,kt,k,j,nx) = fwaves(l,0,kt,k,j,nx);
        }
      }
    });

    //////////////////////////////////////////////////////////
    // Split the flux differences into waves
    //////////////////////////////////////////////////////////
    parallel_for( Bounds<3>(nz,ny,nx+1) , YAKL_LAMBDA (int k, int j, int i) {
      // Compute interface-averaged values
      real r = 0.5_fp * ( fwaves(idR,0,0,k,j,i) + fwaves(idR,1,0,k,j,i) );
      real u = 0.5_fp * ( fwaves(idU,0,0,k,j,i) + fwaves(idU,1,0,k,j,i) );
      real v = 0.5_fp * ( fwaves(idV,0,0,k,j,i) + fwaves(idV,1,0,k,j,i) );
      real w = 0.5_fp * ( fwaves(idW,0,0,k,j,i) + fwaves(idW,1,0,k,j,i) );
      real t = 0.5_fp * ( fwaves(idT,0,0,k,j,i) + fwaves(idT,1,0,k,j,i) );
      real p = 0.5_fp * ( fwaves(idP,0,0,k,j,i) + fwaves(idP,1,0,k,j,i) );
      real cs2 = GAMMA*p/r;
      real cs  = sqrt(cs2);
      for (int kt=0; kt < nTimeDerivs; kt++) {
        // Get left and right state
        real r_L = fwaves(idR,0,kt,k,j,i);
        real u_L = fwaves(idU,0,kt,k,j,i);
        real v_L = fwaves(idV,0,kt,k,j,i);
        real w_L = fwaves(idW,0,kt,k,j,i);
        real t_L = fwaves(idT,0,kt,k,j,i);
        real p_L = fwaves(idP,0,kt,k,j,i);
        real r_R = fwaves(idR,1,kt,k,j,i);
        real u_R = fwaves(idU,1,kt,k,j,i);
        real v_R = fwaves(idV,1,kt,k,j,i);
        real w_R = fwaves(idW,1,kt,k,j,i);
        real t_R = fwaves(idT,1,kt,k,j,i);
        real p_R = fwaves(idP,1,kt,k,j,i);
        // Compute state jump across the interface
        real dr = r_R - r_L;
        real du = u_R - u_L;
        real dv = v_R - v_L;
        real dw = w_R - w_L;
        real dt = t_R - t_L;
        real dp = p_R - p_L;
        // Compute flux jump across the interface
        real df1 = u*dr + r*du;
        real df2 = u*du + dp/r;
        real df3 = u*dv;
        real df4 = u*dw;
        real df5 = u*dt;
        real df6 = u*dp + GAMMA*p*du;
        // Compute characteristics from the jump
        real w1 = -r*df2/(2*cs) + df6/(2*cs2);
        real w2 =  r*df2/(2*cs) + df6/(2*cs2);
        real w3 = df1 - df6/(cs2);
        real w4 = df3;
        real w5 = df4;
        real w6 = df5;
        // set fwaves to zero for this interface
        for (int l=0; l < numState; l++) {
          fwaves(l,0,kt,k,j,i) = 0;
          fwaves(l,1,kt,k,j,i) = 0;
        }
        // Wave 1 (u-cs)
        fwaves(idR,0,kt,k,j,i) = w1;
        fwaves(idU,0,kt,k,j,i) = -cs*w1/r;
        fwaves(idP,0,kt,k,j,i) = cs2*w1;
        // Wave 2 (u+cs)
        fwaves(idR,1,kt,k,j,i) = w2;
        fwaves(idU,1,kt,k,j,i) = cs*w2/r;
        fwaves(idP,1,kt,k,j,i) = cs2*w2;
        // Waves 3-5 (u)
        if (u < 0) {
          fwaves(idR,0,kt,k,j,i) += w3;
          if (! sim2d) {
            fwaves(idV,0,kt,k,j,i) += w4;
          }
          fwaves(idW,0,kt,k,j,i) += w5;
          fwaves(idT,0,kt,k,j,i) += w6;
        } else {
          fwaves(idR,1,kt,k,j,i) += w3;
          if (! sim2d) {
            fwaves(idV,1,kt,k,j,i) += w4;
          }
          fwaves(idW,1,kt,k,j,i) += w5;
          fwaves(idT,1,kt,k,j,i) += w6;
        }
      }
    });

    //////////////////////////////////////////////////////////
    // Add flux waves to the tendencies
    //////////////////////////////////////////////////////////
    parallel_for( Bounds<4>(numState,nz,ny,nx) , YAKL_LAMBDA(int l, int k, int j, int i) {
      for (int kt=0; kt < nTimeDerivs; kt++) {
        tend(l,kt,k,j,i) += - ( fwaves(l,1,kt,k,j,i) + fwaves(l,0,kt,k,j,i+1) ) / dx;
      }
    });
  }



  void computeTendenciesY( StateArr &state , TendArr &tend , real dt ) {
  }



  void computeTendenciesZ( StateArr &state , TendArr &tend , real dt ) {
    // Populate the halos
    if        (bc_z == BC_PERIODIC) {
      parallel_for( Bounds<3>(ny,nx,hs) , YAKL_LAMBDA(int j, int i, int kk) {
        for (int l=0; l < numState; l++) {
          state(l,      kk,hs+j,hs+i) = state(l,nz+kk,hs+j,hs+i);
          state(l,hs+nz+kk,hs+j,hs+i) = state(l,hs+kk,hs+j,hs+i);
        }
      });
    } else if (bc_z == BC_WALL) {
      parallel_for( Bounds<3>(ny,nx,hs) , YAKL_LAMBDA(int j, int i, int kk) {
        for (int l=0; l < numState; l++) {
          if (l == idW) {
            state(l,      kk,hs+j,hs+i) = 0;
            state(l,hs+nz+kk,hs+j,hs+i) = 0;
          } else {
            state(l,      kk,hs+j,hs+i) = state(l,hs     ,hs+j,hs+i);
            state(l,hs+nz+kk,hs+j,hs+i) = state(l,hs+nz-1,hs+j,hs+i);
          }
        }
      });
    }

    // Loop through all cells, reconstruct in x-direction, compute centered tendencies, store cell-edge state estimates
    parallel_for( Bounds<3>(nz,ny,nx) , YAKL_LAMBDA (int k, int j, int i) {
      ////////////////////////////////////////////////////////////////
      // Reconstruct rho, u, v, w, theta, pressure, dwdz, and dpdz
      ////////////////////////////////////////////////////////////////
      SArray<real,2,nAder,ngll> r_DTs;
      SArray<real,2,nAder,ngll> u_DTs;
      SArray<real,2,nAder,ngll> du_DTs;
      SArray<real,2,nAder,ngll> v_DTs;
      SArray<real,2,nAder,ngll> dv_DTs;
      SArray<real,2,nAder,ngll> w_DTs;
      SArray<real,2,nAder,ngll> dw_DTs;
      SArray<real,2,nAder,ngll> t_DTs;
      SArray<real,2,nAder,ngll> dt_DTs;
      SArray<real,2,nAder,ngll> p_DTs;
      SArray<real,2,nAder,ngll> dp_DTs;
      {
        SArray<real,1,ord> stencil;

        // Density
        for (int kk=0; kk < ord; kk++) { stencil(kk) = state(idR,k+kk,hs+j,hs+i); }
        reconstruct_gll_values( stencil , r_DTs , weno_scalars );
        for (int kk=0; kk < ngll; kk++) { r_DTs(0,kk) += hyDensGLL(k,kk); } // Add hydrostasis back on

        // u
        for (int kk=0; kk < ord; kk++) { stencil(kk) = state(idU,k+kk,hs+j,hs+i); }
        reconstruct_gll_values_and_derivs( stencil , u_DTs , du_DTs , dz , weno_winds );

        // v
        for (int kk=0; kk < ord; kk++) { stencil(kk) = state(idV,k+kk,hs+j,hs+i); }
        reconstruct_gll_values_and_derivs( stencil , v_DTs , dv_DTs , dz , weno_winds );

        // w values and derivatives
        for (int kk=0; kk < ord; kk++) { stencil(kk) = state(idW,k+kk,hs+j,hs+i); }
        reconstruct_gll_values_and_derivs( stencil , w_DTs , dw_DTs , dz , weno_winds );
        if (bc_z == BC_WALL) {
          if (k == nz-1) w_DTs(0,ngll-1) = 0;
          if (k == 0   ) w_DTs(0,0     ) = 0;
        }

        // theta
        for (int kk=0; kk < ord; kk++) { stencil(kk) = state(idT,k+kk,hs+j,hs+i); }
        reconstruct_gll_values_and_derivs( stencil , t_DTs , dt_DTs , dz , weno_winds );

        // pressure values and derivatives
        for (int kk=0; kk < ord; kk++) { stencil(kk) = state(idP,k+kk,hs+j,hs+i); }
        reconstruct_gll_values_and_derivs( stencil , p_DTs , dp_DTs , dz , weno_winds );
        // Add hydrostasis, and then divide by C0 to get (rho*theta)^gamma
        for (int kk=0; kk < ngll; kk++) { p_DTs(0,kk) += hyPressureGLL(k,kk); }
      }

      ///////////////////////////////////////////////////////////////
      // Compute other values needed for centered tendencies and DTs
      ///////////////////////////////////////////////////////////////
      SArray<real,2,nAder,ngll> r_w_DTs;
      SArray<real,2,nAder,ngll> w_du_DTs;
      SArray<real,2,nAder,ngll> w_dv_DTs;
      SArray<real,2,nAder,ngll> w_w_DTs;
      SArray<real,2,nAder,ngll> rr_dp_DTs;
      SArray<real,2,nAder,ngll> w_dt_DTs;
      SArray<real,2,nAder,ngll> w_dp_DTs;
      SArray<real,2,nAder,ngll> p_dw_DTs;
      for (int kk=0; kk < ngll; kk++) {
        real r  = r_DTs (0,kk);
        real u  = u_DTs (0,kk);
        real v  = v_DTs (0,kk);
        real w  = w_DTs (0,kk);
        real t  = t_DTs (0,kk);
        real p  = p_DTs (0,kk);
        real du = du_DTs(0,kk);
        real dv = dv_DTs(0,kk);
        real dw = dw_DTs(0,kk);
        real dt = dt_DTs(0,kk);
        real dp = dp_DTs(0,kk);
        r_w_DTs  (0,kk) = r*w;
        w_du_DTs (0,kk) = w*du;
        w_dv_DTs (0,kk) = w*dv;
        w_w_DTs  (0,kk) = w*w;
        rr_dp_DTs(0,kk) = dp/r;
        w_dt_DTs (0,kk) = w*dt;
        w_dp_DTs (0,kk) = w*(dp+hyPressureDerivGLL(k,kk));
        p_dw_DTs (0,kk) = p*dw;
      }

      //////////////////////////////////////////
      // Compute time derivatives if necessary
      //////////////////////////////////////////
      if (nAder > 1) {
        for (int kt=0; kt < nAder-1; kt++) {
          // Compute state at kt+1
          for (int ii=0; ii<ngll; ii++) {
            // Compute d_dz(r*w) and d_dz(w*w)
            real drw_dz  = 0;
            real dww_dz  = 0;
            for (int s=0; s<ngll; s++) {
              drw_dz += derivMatrix(s,ii) * r_w_DTs(kt,s);
              dww_dz += derivMatrix(s,ii) * w_w_DTs(kt,s);
            }
            drw_dz /= dz;
            dww_dz /= dz;
            // Compute state at kt+1
            r_DTs(kt+1,ii) = -( drw_dz                                  ) / (kt+1);
            u_DTs(kt+1,ii) = -( w_du_DTs(kt,ii)                         ) / (kt+1);
            v_DTs(kt+1,ii) = -( w_dv_DTs(kt,ii)                         ) / (kt+1);
            w_DTs(kt+1,ii) = -( dww_dz/2 + rr_dp_DTs(kt,ii)             ) / (kt+1);
            t_DTs(kt+1,ii) = -( w_dt_DTs(kt,ii)                         ) / (kt+1);
            p_DTs(kt+1,ii) = -( w_dp_DTs(kt,ii) + GAMMA*p_dw_DTs(kt,ii) ) / (kt+1);
          }
          if (bc_z == BC_WALL) {
            if (k == nz-1) w_DTs(kt+1,ngll-1) = 0;
            if (k == 0   ) w_DTs(kt+1,0     ) = 0;
          }
          // Compute du, dv, dw, dt, dp, r*w, w*w, dp/r, w*du, w*dv, w*dt, w*dp, p*dw
          for (int ii=0; ii<ngll; ii++) {
            // Differentiate du, dv, dw, dt, and dp at kt+1
            real du_dz = 0;
            real dv_dz = 0;
            real dw_dz = 0;
            real dt_dz = 0;
            real dp_dz = 0;
            for (int s=0; s<ngll; s++) {
              du_dz += derivMatrix(s,ii) * u_DTs(kt+1,s);
              dv_dz += derivMatrix(s,ii) * v_DTs(kt+1,s);
              dw_dz += derivMatrix(s,ii) * w_DTs(kt+1,s);
              dt_dz += derivMatrix(s,ii) * t_DTs(kt+1,s);
              dp_dz += derivMatrix(s,ii) * p_DTs(kt+1,s);
            }
            du_DTs(kt+1,ii) = du_dz / dz;
            dv_DTs(kt+1,ii) = dv_dz / dz;
            dw_DTs(kt+1,ii) = dw_dz / dz;
            dt_DTs(kt+1,ii) = dt_dz / dz;
            dp_DTs(kt+1,ii) = dp_dz / dz;
            // Compute r*w, w*du, w*dv, w*w, dp/r, w*dt, w*dp, p*dw
            real tot_r_w   = 0;
            real tot_w_du  = 0;
            real tot_w_dv  = 0;
            real tot_w_w   = 0;
            real tot_rr_dp = 0;
            real tot_w_dt  = 0;
            real tot_w_dp  = 0;
            real tot_p_dw  = 0;
            rr_dp_DTs(kt+1,ii) = 0;
            for (int rt=0; rt <= kt+1; rt++) {
              tot_r_w   += r_DTs(rt,ii) * w_DTs (kt+1-rt,ii);
              tot_w_du  += w_DTs(rt,ii) * du_DTs(kt+1-rt,ii);
              tot_w_dv  += w_DTs(rt,ii) * dv_DTs(kt+1-rt,ii);
              tot_w_w   += w_DTs(rt,ii) * w_DTs (kt+1-rt,ii);
              tot_rr_dp += dp_DTs(rt,ii) - r_DTs(rt,ii) * rr_dp_DTs(kt+1-rt,ii);
              tot_w_dt  += w_DTs(rt,ii) * dt_DTs(kt+1-rt,ii);
              tot_w_dp  += w_DTs(rt,ii) * dp_DTs(kt+1-rt,ii);
              tot_p_dw  += p_DTs(rt,ii) * dw_DTs(kt+1-rt,ii);
            }
            r_w_DTs  (kt+1,ii) = tot_r_w ;
            w_du_DTs (kt+1,ii) = tot_w_du;
            w_dv_DTs (kt+1,ii) = tot_w_dv;
            w_w_DTs  (kt+1,ii) = tot_w_w ;
            rr_dp_DTs(kt+1,ii) = tot_rr_dp / r_DTs(0,ii);
            w_dt_DTs (kt+1,ii) = tot_w_dt;
            w_dp_DTs (kt+1,ii) = tot_w_dp;
            p_dw_DTs (kt+1,ii) = tot_p_dw;
          }
        }
      }

      //////////////////////////////////////////
      // Time average if necessary
      //////////////////////////////////////////
      if (timeAvg) {
        // Compute time averages
        for (int ii=0; ii<ngll; ii++) {
          real dtmult = 1;
          real r_tavg     = 0;
          real u_tavg     = 0;
          real v_tavg     = 0;
          real w_tavg     = 0;
          real t_tavg     = 0;
          real p_tavg     = 0;
          real r_w_tavg   = 0;
          real w_du_tavg  = 0;
          real w_dv_tavg  = 0;
          real w_w_tavg   = 0;
          real rr_dp_tavg = 0;
          real w_dt_tavg  = 0;
          real w_dp_tavg  = 0;
          real p_dw_tavg  = 0;
          for (int kt=0; kt<nAder; kt++) {
            r_tavg     += r_DTs    (kt,ii) * dtmult / (kt+1);
            u_tavg     += u_DTs    (kt,ii) * dtmult / (kt+1);
            v_tavg     += v_DTs    (kt,ii) * dtmult / (kt+1);
            w_tavg     += w_DTs    (kt,ii) * dtmult / (kt+1);
            t_tavg     += t_DTs    (kt,ii) * dtmult / (kt+1);
            p_tavg     += p_DTs    (kt,ii) * dtmult / (kt+1);
            r_w_tavg   += r_w_DTs  (kt,ii) * dtmult / (kt+1);
            w_du_tavg  += w_du_DTs (kt,ii) * dtmult / (kt+1);
            w_dv_tavg  += w_dv_DTs (kt,ii) * dtmult / (kt+1);
            w_w_tavg   += w_w_DTs  (kt,ii) * dtmult / (kt+1);
            rr_dp_tavg += rr_dp_DTs(kt,ii) * dtmult / (kt+1);
            w_dt_tavg  += w_dt_DTs (kt,ii) * dtmult / (kt+1);
            w_dp_tavg  += w_dp_DTs (kt,ii) * dtmult / (kt+1);
            p_dw_tavg  += p_dw_DTs (kt,ii) * dtmult / (kt+1);
            dtmult *= dt;
          }
          r_DTs    (0,ii) = r_tavg    ;
          u_DTs    (0,ii) = u_tavg    ;
          v_DTs    (0,ii) = v_tavg    ;
          w_DTs    (0,ii) = w_tavg    ;
          t_DTs    (0,ii) = t_tavg    ;
          p_DTs    (0,ii) = p_tavg    ;
          r_w_DTs  (0,ii) = r_w_tavg  ;
          w_du_DTs (0,ii) = w_du_tavg ;
          w_dv_DTs (0,ii) = w_dv_tavg ;
          w_w_DTs  (0,ii) = w_w_tavg  ;
          rr_dp_DTs(0,ii) = rr_dp_tavg;
          w_dt_DTs (0,ii) = w_dt_tavg ;
          w_dp_DTs (0,ii) = w_dp_tavg ;
          p_dw_DTs (0,ii) = p_dw_tavg ;
        }
      }

      //////////////////////////////////////////
      // Compute centered tendencies
      //////////////////////////////////////////
      // Compute flux-form updates
      real drdt = - ( r_w_DTs(0,ngll-1) - r_w_DTs(0,0) ) / dz;
      real dudt = 0;
      real dvdt = 0;
      real dwdt = - ( w_w_DTs(0,ngll-1) - w_w_DTs(0,0) ) / dz / 2;
      real dtdt = 0;
      real dpdt = 0;
      // for (int kk=0; kk < ngll; kk++) {
      //   real wt = gllWts_ngll(kk);
      //   dudt += -w_du_DTs (0,kk) * wt;
      //   if (! sim2d) {
      //     dvdt += -w_dv_DTs (0,kk) * wt;
      //   }
      //   dwdt += -rr_dp_DTs(0,kk) * wt;
      //   dtdt += -w_dt_DTs(0,kk) * wt;
      //   dpdt += ( -w_dp_DTs(0,kk) - GAMMA*p_dw_DTs(0,kk) ) * wt;
      // }
      real r = state(idR,hs+k,hs+j,hs+i) + hyDensCells(k);
      real w = state(idW,hs+k,hs+j,hs+i);
      real p = state(idP,hs+k,hs+j,hs+i) + hyPressureCells(k);
      dudt += - w * ( u_DTs(0,ngll-1) - u_DTs(0,0) ) / dz;
      if (! sim2d) {
        dvdt += - w * ( v_DTs(0,ngll-1) - v_DTs(0,0) ) / dz;
      }
      dwdt += - (1./r) * ( (p_DTs(0,ngll-1)-hyPressureGLL(k,ngll-1)) - (p_DTs(0,0)-hyPressureGLL(k,0)) ) / dz;
      dtdt += - w * ( t_DTs(0,ngll-1) - t_DTs(0,0) ) / dz;
      dpdt += - w *       ( p_DTs(0,ngll-1) - p_DTs(0,0) ) / dz
              - GAMMA*p * ( w_DTs(0,ngll-1) - w_DTs(0,0) ) / dz;

      tend(idR,0,k,j,i) = drdt;
      tend(idU,0,k,j,i) = dudt;
      tend(idV,0,k,j,i) = dvdt;
      tend(idW,0,k,j,i) = dwdt;
      tend(idT,0,k,j,i) = dtdt;
      tend(idP,0,k,j,i) = dpdt;

      //////////////////////////////////////////
      // Store cell edge estimates of the state
      //////////////////////////////////////////
      fwaves(idR,1,0,k  ,j,i) = r_DTs(0,0     );
      fwaves(idU,1,0,k  ,j,i) = u_DTs(0,0     );
      fwaves(idV,1,0,k  ,j,i) = v_DTs(0,0     );
      fwaves(idW,1,0,k  ,j,i) = w_DTs(0,0     );
      fwaves(idT,1,0,k  ,j,i) = t_DTs(0,0     );
      fwaves(idP,1,0,k  ,j,i) = p_DTs(0,0     );

      fwaves(idR,0,0,k+1,j,i) = r_DTs(0,ngll-1);
      fwaves(idU,0,0,k+1,j,i) = u_DTs(0,ngll-1);
      fwaves(idV,0,0,k+1,j,i) = v_DTs(0,ngll-1);
      fwaves(idW,0,0,k+1,j,i) = w_DTs(0,ngll-1);
      fwaves(idT,0,0,k+1,j,i) = t_DTs(0,ngll-1);
      fwaves(idP,0,0,k+1,j,i) = p_DTs(0,ngll-1);

    });

    ////////////////////////////////////////////////
    // BCs for the state edge estimates
    ////////////////////////////////////////////////
    parallel_for( Bounds<4>(numState,nTimeDerivs,ny,nx) , YAKL_LAMBDA (int l, int kt, int j, int i) {
      if        (bc_z == BC_PERIODIC) {
        fwaves(l,0,kt,0 ,j,i) = fwaves(l,0,kt,nz,j,i);
        fwaves(l,1,kt,nz,j,i) = fwaves(l,1,kt,0 ,j,i);
      } else if (bc_z == BC_WALL    ) {
        if (l == idW) {
          fwaves(l,0,kt,0 ,j,i) = 0;
          fwaves(l,1,kt,0 ,j,i) = 0;
          fwaves(l,0,kt,nz,j,i) = 0;
          fwaves(l,1,kt,nz,j,i) = 0;
        } else {
          fwaves(l,0,kt,0 ,j,i) = fwaves(l,1,kt,0 ,j,i);
          fwaves(l,1,kt,nz,j,i) = fwaves(l,0,kt,nz,j,i);
        }
      }
    });

    //////////////////////////////////////////////////////////
    // Split the flux differences into waves
    //////////////////////////////////////////////////////////
    parallel_for( Bounds<3>(nz+1,ny,nx) , YAKL_LAMBDA (int k, int j, int i) {
      // Compute interface-averaged values
      real r = 0.5_fp * ( fwaves(idR,0,0,k,j,i) + fwaves(idR,1,0,k,j,i) );
      real u = 0.5_fp * ( fwaves(idU,0,0,k,j,i) + fwaves(idU,1,0,k,j,i) );
      real v = 0.5_fp * ( fwaves(idV,0,0,k,j,i) + fwaves(idV,1,0,k,j,i) );
      real w = 0.5_fp * ( fwaves(idW,0,0,k,j,i) + fwaves(idW,1,0,k,j,i) );
      real t = 0.5_fp * ( fwaves(idT,0,0,k,j,i) + fwaves(idT,1,0,k,j,i) );
      real p = 0.5_fp * ( fwaves(idP,0,0,k,j,i) + fwaves(idP,1,0,k,j,i) );
      real cs2 = GAMMA*p/r;
      real cs  = sqrt(cs2);
      for (int kt=0; kt < nTimeDerivs; kt++) {
        // Get left and right state
        real r_L = fwaves(idR,0,kt,k,j,i);
        real u_L = fwaves(idU,0,kt,k,j,i);
        real v_L = fwaves(idV,0,kt,k,j,i);
        real w_L = fwaves(idW,0,kt,k,j,i);
        real t_L = fwaves(idT,0,kt,k,j,i);
        real p_L = fwaves(idP,0,kt,k,j,i);
        real r_R = fwaves(idR,1,kt,k,j,i);
        real u_R = fwaves(idU,1,kt,k,j,i);
        real v_R = fwaves(idV,1,kt,k,j,i);
        real w_R = fwaves(idW,1,kt,k,j,i);
        real t_R = fwaves(idT,1,kt,k,j,i);
        real p_R = fwaves(idP,1,kt,k,j,i);
        // Compute state jump across the interface
        real dr = r_R - r_L;
        real du = u_R - u_L;
        real dv = v_R - v_L;
        real dw = w_R - w_L;
        real dt = t_R - t_L;
        real dp = p_R - p_L;
        // Compute flux jump across the interface
        real df1 = w*dr + r*dw;
        real df2 = w*du;
        real df3 = w*dv;
        real df4 = w*dw + dp/r; // dp and dp' have no difference. Hydrostasis is the same for L and R estimates
        real df5 = w*dt;
        real df6 = w*dp + GAMMA*p*dw;
        // Compute characteristics from the jump
        real w1 = -r*df4/(2*cs) + df6/(2*cs2);
        real w2 =  r*df4/(2*cs) + df6/(2*cs2);
        real w3 = df1 - df6/cs2;
        real w4 = df2;
        real w5 = df3;
        real w6 = df5;
        // set fwaves to zero for this interface
        for (int l=0; l < numState; l++) {
          fwaves(l,0,kt,k,j,i) = 0;
          fwaves(l,1,kt,k,j,i) = 0;
        }
        // Wave 1 (w-cs)
        fwaves(idR,0,kt,k,j,i) = w1;
        fwaves(idW,0,kt,k,j,i) = -cs*w1/r;
        fwaves(idP,0,kt,k,j,i) = cs2*w1;
        // Wave 2 (w+cs)
        fwaves(idR,1,kt,k,j,i) = w2;
        fwaves(idW,1,kt,k,j,i) = cs*w2/r;
        fwaves(idP,1,kt,k,j,i) = cs2*w2;
        // Waves 3-6 (w)
        if (w < 0) {
          fwaves(idR,0,kt,k,j,i) += w3;
          fwaves(idU,0,kt,k,j,i) += w4;
          if (! sim2d) {
            fwaves(idV,0,kt,k,j,i) += w5;
          }
          fwaves(idT,0,kt,k,j,i) += w6;
        } else {
          fwaves(idR,1,kt,k,j,i) += w3;
          fwaves(idU,1,kt,k,j,i) += w4;
          if (! sim2d) {
            fwaves(idV,1,kt,k,j,i) += w5;
          }
          fwaves(idT,1,kt,k,j,i) += w6;
        }
      }
    });

    //////////////////////////////////////////////////////////
    // Add flux waves to the tendencies
    //////////////////////////////////////////////////////////
    parallel_for( Bounds<4>(numState,nz,ny,nx) , YAKL_LAMBDA(int l, int k, int j, int i) {
      for (int kt=0; kt < nTimeDerivs; kt++) {
        tend(l,kt,k,j,i) += - ( fwaves(l,1,kt,k,j,i) + fwaves(l,0,kt,k+1,j,i) ) / dz;
      }
    });
  }



  void computeTendenciesS( StateArr &state , TendArr &tend , real dt ) {
    parallel_for( Bounds<3>(nz,ny,nx) , YAKL_LAMBDA(int k, int j, int i) {
      real rp = state(idR,hs+k,hs+j,hs+i);
      real r  = rp + hyDensCells(k);
      tend(idR,0,k,j,i) = 0;
      tend(idU,0,k,j,i) = 0;
      tend(idV,0,k,j,i) = 0;
      tend(idW,0,k,j,i) = -GRAV * rp / r;
      tend(idT,0,k,j,i) = 0;
      tend(idP,0,k,j,i) = 0;
    });
  }



  template <class F> void applyTendencies( F const &applySingleTendency , int splitIndex ) {
    parallel_for( Bounds<4>(numState,nz,ny,nx) , YAKL_LAMBDA (int l, int k, int j, int i) {
      applySingleTendency({l,k,j,i});
    });
  }



  const char * getSpatialName() { return "1-D Uniform Transport with upwind FV on A-grid"; }



  void output(StateArr const &state, real etime) {
    yakl::SimpleNetCDF nc;
    int ulIndex = 0; // Unlimited dimension index to place this data at
    // Create or open the file
    if (etime == 0.) {
      nc.create(outFile);
      // x-coordinate
      real1d xloc("xloc",nx);
      parallel_for( nx , YAKL_LAMBDA (int i) { xloc(i) = (i+0.5)*dx; });
      nc.write(xloc,"x",{"x"});
      // y-coordinate
      real1d yloc("yloc",ny);
      parallel_for( ny , YAKL_LAMBDA (int i) { yloc(i) = (i+0.5)*dy; });
      nc.write(yloc,"y",{"y"});
      // z-coordinate
      real1d zloc("zloc",nz);
      parallel_for( nz , YAKL_LAMBDA (int i) { zloc(i) = (i+0.5)*dz; });
      nc.write(zloc,"z",{"z"});
      // hydrostatic density, theta, and pressure
      nc.write(hyDensCells    ,"hyDens"    ,{"z"});
      nc.write(hyPressureCells,"hyPressure",{"z"});
      // Create time variable
      nc.write1(0._fp,"t",0,"t");
    } else {
      nc.open(outFile,yakl::NETCDF_MODE_WRITE);
      ulIndex = nc.getDimSize("t");
      // Write the elapsed time
      nc.write1(etime,"t",ulIndex,"t");
    }
    real3d data("data",nz,ny,nx);
    // rho'
    parallel_for( Bounds<3>(nz,ny,nx) , YAKL_LAMBDA (int k, int j, int i) { data(k,j,i) = state(idR,hs+k,hs+j,hs+i); });
    nc.write1(data,"dens_pert",{"z","y","x"},ulIndex,"t");
    // u
    parallel_for( Bounds<3>(nz,ny,nx) , YAKL_LAMBDA (int k, int j, int i) { data(k,j,i) = state(idU,hs+k,hs+j,hs+i); });
    nc.write1(data,"u",{"z","y","x"},ulIndex,"t");
    // v
    parallel_for( Bounds<3>(nz,ny,nx) , YAKL_LAMBDA (int k, int j, int i) { data(k,j,i) = state(idV,hs+k,hs+j,hs+i); });
    nc.write1(data,"v",{"z","y","x"},ulIndex,"t");
    // w
    parallel_for( Bounds<3>(nz,ny,nx) , YAKL_LAMBDA (int k, int j, int i) { data(k,j,i) = state(idW,hs+k,hs+j,hs+i); });
    nc.write1(data,"w",{"z","y","x"},ulIndex,"t");
    // theta'
    parallel_for( Bounds<3>(nz,ny,nx) , YAKL_LAMBDA (int k, int j, int i) { data(k,j,i) = state(idT,hs+k,hs+j,hs+i); });
    nc.write1(data,"pot_temp_pert",{"z","y","x"},ulIndex,"t");
    // pressure'
    parallel_for( Bounds<3>(nz,ny,nx) , YAKL_LAMBDA (int k, int j, int i) { data(k,j,i) = state(idP,hs+k,hs+j,hs+i); });
    nc.write1(data,"pressure_pert",{"z","y","x"},ulIndex,"t");
    // Close the file
    nc.close();
  }



  void finalize(StateArr const &state) {}



  // ord stencil values to ngll GLL values and ngll GLL derivatives; store in DTs
  YAKL_INLINE void reconstruct_gll_values_and_derivs( SArray<real,1,ord> const &stencil , SArray<real,2,nAder,ngll> &DTs ,
                                                      SArray<real,2,nAder,ngll> &deriv_DTs, real dx , bool doweno ) {
    if (doweno) {

      // Reconstruct values
      SArray<real,1,ord> wenoCoefs;
      weno::compute_weno_coefs( wenoRecon , stencil , wenoCoefs , idl , sigma );
      // Transform ord weno coefficients into ngll GLL points
      for (int ii=0; ii<ngll; ii++) {
        real tmp       = 0;
        real deriv_tmp = 0;
        for (int s=0; s < ord; s++) {
          real coef = wenoCoefs(s);
          tmp       += coefs_to_gll      (s,ii) * coef;
          deriv_tmp += coefs_to_deriv_gll(s,ii) * coef;
        }
        DTs      (0,ii) = tmp;
        deriv_DTs(0,ii) = deriv_tmp / dx;
      }

    } else {

      // Transform ord stencil cell averages into ngll GLL points
      for (int ii=0; ii<ngll; ii++) {
        real tmp       = 0;
        real deriv_tmp = 0;
        for (int s=0; s < ord; s++) {
          real sten = stencil(s);
          tmp       += sten_to_gll      (s,ii) * sten;
          deriv_tmp += sten_to_deriv_gll(s,ii) * sten;
        }
        DTs      (0,ii) = tmp;
        deriv_DTs(0,ii) = deriv_tmp / dx;
      }

    } // if doweno
  }



  // ord stencil values to ngll GLL values; store in DTs
  YAKL_INLINE void reconstruct_gll_values( SArray<real,1,ord> const stencil , SArray<real,2,nAder,ngll> &DTs , bool doweno ) {
    if (doweno) {

      // Reconstruct values
      SArray<real,1,ord> wenoCoefs;
      weno::compute_weno_coefs( wenoRecon , stencil , wenoCoefs , idl , sigma );
      // Transform ord weno coefficients into ngll GLL points
      for (int ii=0; ii<ngll; ii++) {
        real tmp = 0;
        for (int s=0; s < ord; s++) {
          tmp += coefs_to_gll(s,ii) * wenoCoefs(s);
        }
        DTs(0,ii) = tmp;
      }

    } else {

      // Transform ord stencil cell averages into ngll GLL points
      for (int ii=0; ii<ngll; ii++) {
        real tmp = 0;
        for (int s=0; s < ord; s++) {
          tmp += sten_to_gll(s,ii) * stencil(s);
        }
        DTs(0,ii) = tmp;
      }

    } // if doweno
  }


};
