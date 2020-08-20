
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
  int static constexpr numState = 5;

  // Stores a single index location
  struct Location {
    int l;
    int k;
    int j;
    int i;
  };

  typedef real4d StateArr;  // Spatial index

  typedef real5d TendArr;   // (time derivative & spatial index)

  real3d pressure;
  real6d fwaves;            // state edge estimates and flux difference splitting waves
  // Hydrostatically balanced values for density, potential temperature, and pressure
  real1d hyDensCells;
  real1d hyPressureCells;
  real2d hyDensGLL;
  real2d hyPressureGLL;
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
  int static constexpr idT = 4;  // potential temperature perturbation

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
    return 3;
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
        real p = C0*pow(r*t,GAMMA);
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

    pressure = real3d("pressure",nz+2*hs,ny+2*hs,nx+2*hs);
    fwaves = real6d("fwaves",numState,2,nTimeDerivs,nz+1,ny+1,nx+1);
    hyDensCells          = real1d("hyDensCells       ",nz     );
    hyPressureCells      = real1d("hyPressureCells   ",nz     );
    hyDensGLL            = real2d("hyDensGLL         ",nz,ngll);
    hyPressureGLL        = real2d("hyPressureGLL     ",nz,ngll);
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
          hyDensGLL    (k,kk) = profiles::initConstTheta_density (t,zloc);
          hyPressureGLL(k,kk) = profiles::initConstTheta_pressure(t,zloc);
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
              // Compute constant theta hydrostatic background state
              real tp = profiles::ellipsoid_linear(xloc, yloc, zloc, xlen/2, ylen/2, 2000, 2000, 2000, 2000, 2 );
              state(idT,hs+k,hs+j,hs+i) += (300+tp) * wt;
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
      }
    } else {
      if        (splitIndex == 0) {
        computeTendenciesZ( state , tend , dt );
      } else if (splitIndex == 1) {
        if (sim2d) {
          memset(tend,0._fp);
        } else {
          computeTendenciesY( state , tend , dt );
        }
      } else if (splitIndex == 2) {
        computeTendenciesX( state , tend , dt );
      }
    }
    if (splitIndex == numSplit()) dimSwitch = ! dimSwitch;
  } // computeTendencies



  void computeTendenciesX( StateArr &state , TendArr &tend , real dt ) {
    // Precompute pressure perturbation
    parallel_for( Bounds<3>(nz,ny,nx) , YAKL_LAMBDA(int k, int j, int i) {
      real r = state(idR,hs+k,hs+j,hs+i) + hyDensCells(k);
      real t = state(idT,hs+k,hs+j,hs+i);
      pressure(hs+k,hs+j,hs+i) = C0*pow(r*t,GAMMA) - hyPressureCells(k);
    });

    // Populate the halos
    if        (bc_x == BC_PERIODIC) {
      parallel_for( Bounds<3>(nz,ny,hs) , YAKL_LAMBDA(int k, int j, int ii) {
        for (int l=0; l < numState; l++) {
          state(l,hs+k,hs+j,      ii) = state(l,hs+k,hs+j,nx+ii);
          state(l,hs+k,hs+j,hs+nx+ii) = state(l,hs+k,hs+j,hs+ii);
        }
        pressure(hs+k,hs+j,      ii) = pressure(hs+k,hs+j,nx+ii);
        pressure(hs+k,hs+j,hs+nx+ii) = pressure(hs+k,hs+j,hs+ii);
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
        pressure(hs+k,hs+j,      ii) = pressure(hs+k,hs+j,hs     );
        pressure(hs+k,hs+j,hs+nx+ii) = pressure(hs+k,hs+j,hs+nx-1);
      });
    }

    // Loop through all cells, reconstruct in x-direction, compute centered tendencies, store cell-edge state estimates
    parallel_for( Bounds<3>(nz,ny,nx) , YAKL_LAMBDA (int k, int j, int i) {
      ////////////////////////////////////////////////////////////////
      // Reconstruct rho, u, v, w, theta, pressure, dudx, and dpdx
      ////////////////////////////////////////////////////////////////
      SArray<real,2,nAder,ngll> r_DTs;
      SArray<real,2,nAder,ngll> u_DTs;
      SArray<real,2,nAder,ngll> v_DTs;
      SArray<real,2,nAder,ngll> dv_DTs;
      SArray<real,2,nAder,ngll> w_DTs;
      SArray<real,2,nAder,ngll> dw_DTs;
      SArray<real,2,nAder,ngll> t_DTs;
      SArray<real,2,nAder,ngll> dt_DTs;
      SArray<real,2,nAder,ngll> rtgamma_DTs;
      SArray<real,2,nAder,ngll> dp_DTs;
      {
        SArray<real,1,ord> stencil;

        // Density
        for (int ii=0; ii < ord; ii++) { stencil(ii) = state(idR,hs+k,hs+j,i+ii); }
        reconstruct_gll_values( stencil , r_DTs , weno_scalars );
        for (int ii=0; ii < ngll; ii++) { r_DTs(0,ii) += hyDensCells(k); } // Add hydrostasis back on

        // u values and derivatives
        for (int ii=0; ii < ord; ii++) { stencil(ii) = state(idU,hs+k,hs+j,i+ii); }
        reconstruct_gll_values( stencil , u_DTs , weno_winds );

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
        for (int ii=0; ii < ord; ii++) { stencil(ii) = pressure(hs+k,hs+j,i+ii); }
        reconstruct_gll_values_and_derivs( stencil , rtgamma_DTs , dp_DTs , dx , weno_scalars );
        // Add hydrostasis, and then divide by C0 to get (rho*theta)^gamma
        for (int ii=0; ii < ngll; ii++) { rtgamma_DTs(0,ii) = ( rtgamma_DTs(0,ii) + hyPressureCells(k) ) / C0; }
      }

      ///////////////////////////////////////////////////////////////
      // Compute other values needed for centered tendencies and DTs
      ///////////////////////////////////////////////////////////////
      SArray<real,2,nAder,ngll> r_u_DTs;
      SArray<real,2,nAder,ngll> u_u_DTs;
      SArray<real,2,nAder,ngll> r_t_DTs;
      SArray<real,2,nAder,ngll> rr_dp_DTs;
      SArray<real,2,nAder,ngll> u_dv_DTs;
      SArray<real,2,nAder,ngll> u_dw_DTs;
      SArray<real,2,nAder,ngll> u_dt_DTs;
      for (int ii=0; ii < ngll; ii++) {
        real r  = r_DTs (0,ii);
        real u  = u_DTs (0,ii);
        real v  = v_DTs (0,ii);
        real w  = w_DTs (0,ii);
        real t  = t_DTs (0,ii);
        real dv = dv_DTs(0,ii);
        real dw = dw_DTs(0,ii);
        real dt = dt_DTs(0,ii);
        real dp = dp_DTs(0,ii);
        r_u_DTs  (0,ii) = r*u;
        u_u_DTs  (0,ii) = u*u;
        r_t_DTs  (0,ii) = r*t;
        rr_dp_DTs(0,ii) = dp/r;
        u_dv_DTs (0,ii) = u*dv;
        u_dw_DTs (0,ii) = u*dw;
        u_dt_DTs (0,ii) = u*dt;
      }

      //////////////////////////////////////////
      // Compute time derivatives if necessary
      //////////////////////////////////////////
      if (nAder > 1) {
      }

      //////////////////////////////////////////
      // Time average if necessary
      //////////////////////////////////////////
      if (timeAvg) {
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
      for (int ii=0; ii < ngll; ii++) {
        real wt = gllWts_ngll(ii);
        dudt += -rr_dp_DTs(0,ii) * wt;
        if (! sim2d) {
          dvdt += -u_dv_DTs (0,ii) * wt;
        }
        dwdt += -u_dw_DTs (0,ii) * wt;
        dtdt += -u_dt_DTs (0,ii) * wt;
      }
      tend(idR,0,k,j,i) = drdt;
      tend(idU,0,k,j,i) = dudt;
      tend(idV,0,k,j,i) = dvdt;
      tend(idW,0,k,j,i) = dwdt;
      tend(idT,0,k,j,i) = dtdt;

      //////////////////////////////////////////
      // Store cell edge estimates of the state
      //////////////////////////////////////////
      fwaves(idR,1,0,k,j,i  ) = r_DTs(0,0     );
      fwaves(idU,1,0,k,j,i  ) = u_DTs(0,0     );
      fwaves(idV,1,0,k,j,i  ) = v_DTs(0,0     );
      fwaves(idW,1,0,k,j,i  ) = w_DTs(0,0     );
      fwaves(idT,1,0,k,j,i  ) = t_DTs(0,0     );

      fwaves(idR,0,0,k,j,i+1) = r_DTs(0,ngll-1);
      fwaves(idU,0,0,k,j,i+1) = u_DTs(0,ngll-1);
      fwaves(idV,0,0,k,j,i+1) = v_DTs(0,ngll-1);
      fwaves(idW,0,0,k,j,i+1) = w_DTs(0,ngll-1);
      fwaves(idT,0,0,k,j,i+1) = t_DTs(0,ngll-1);

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
      real p = C0*pow(r*t,GAMMA);
      real cs2 = GAMMA*p/r;
      real cs  = sqrt(cs2);
      for (int kt=0; kt < nTimeDerivs; kt++) {
        // Get left and right state
        real r_L = fwaves(idR,0,kt,k,j,i);
        real u_L = fwaves(idU,0,kt,k,j,i);
        real v_L = fwaves(idV,0,kt,k,j,i);
        real w_L = fwaves(idW,0,kt,k,j,i);
        real t_L = fwaves(idT,0,kt,k,j,i);
        real p_L = C0*pow(r_L*t_L,GAMMA);
        real r_R = fwaves(idR,1,kt,k,j,i);
        real u_R = fwaves(idU,1,kt,k,j,i);
        real v_R = fwaves(idV,1,kt,k,j,i);
        real w_R = fwaves(idW,1,kt,k,j,i);
        real t_R = fwaves(idT,1,kt,k,j,i);
        real p_R = C0*pow(r_R*t_R,GAMMA);
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
        // Compute characteristics from the jump
        real w1 = 0.5_fp*df1 - r*df2/(2*cs) + r*df5/(2*t);
        real w2 = 0.5_fp*df1 + r*df2/(2*cs) + r*df5/(2*t);
        real w3 = -r*df5/t;
        real w4 = df3;
        real w5 = df4;
        // set fwaves to zero for this interface
        for (int l=0; l < numState; l++) {
          fwaves(l,0,kt,k,j,i) = 0;
          fwaves(l,1,kt,k,j,i) = 0;
        }
        // Wave 1 (u-cs)
        fwaves(idR,0,kt,k,j,i) = w1;
        fwaves(idU,0,kt,k,j,i) = -cs*w1/r;
        // Wave 2 (u+cs)
        fwaves(idR,1,kt,k,j,i) = w2;
        fwaves(idU,1,kt,k,j,i) = cs*w2/r;
        // Waves 3-5 (u)
        if (u < 0) {
          fwaves(idR,0,kt,k,j,i) += w3;
          if (! sim2d) {
            fwaves(idV,0,kt,k,j,i) += w4;
          }
          fwaves(idW,0,kt,k,j,i) += w5;
          fwaves(idT,0,kt,k,j,i) += -t*w3/r;
        } else {
          fwaves(idR,1,kt,k,j,i) += w3;
          if (! sim2d) {
            fwaves(idV,1,kt,k,j,i) += w4;
          }
          fwaves(idW,1,kt,k,j,i) += w5;
          fwaves(idT,1,kt,k,j,i) += -t*w3/r;
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
    // Precompute pressure perturbation
    parallel_for( Bounds<3>(nz,ny,nx) , YAKL_LAMBDA(int k, int j, int i) {
      real r = state(idR,hs+k,hs+j,hs+i) + hyDensCells(k);
      real t = state(idT,hs+k,hs+j,hs+i);
      pressure(hs+k,hs+j,hs+i) = C0*pow(r*t,GAMMA) - hyPressureCells(k);
    });

    // Populate the halos
    if        (bc_y == BC_PERIODIC) {
      parallel_for( Bounds<3>(nz,nx,hs) , YAKL_LAMBDA(int k, int i, int jj) {
        for (int l=0; l < numState; l++) {
          state(l,hs+k,      jj,hs+i) = state(l,hs+k,ny+jj,hs+i);
          state(l,hs+k,hs+ny+jj,hs+i) = state(l,hs+k,hs+jj,hs+i);
        }
        pressure(hs+k,      jj,hs+i) = pressure(hs+k,ny+jj,hs+i);
        pressure(hs+k,hs+ny+jj,hs+i) = pressure(hs+k,hs+jj,hs+i);
      });
    } else if (bc_y == BC_WALL) {
      parallel_for( Bounds<3>(nz,nx,hs) , YAKL_LAMBDA(int k, int i, int jj) {
        for (int l=0; l < numState; l++) {
          if (l == idV) {
            state(l,hs+k,      jj,hs+i) = 0;
            state(l,hs+k,hs+ny+jj,hs+i) = 0;
          } else {
            state(l,hs+k,      jj,hs+i) = state(l,hs+k,hs     ,hs+i);
            state(l,hs+k,hs+ny+jj,hs+i) = state(l,hs+k,hs+ny-1,hs+i);
          }
        }
        pressure(hs+k,      jj,hs+i) = pressure(hs+k,hs     ,hs+i);
        pressure(hs+k,hs+ny+jj,hs+i) = pressure(hs+k,hs+ny-1,hs+i);
      });
    }

    // Loop through all cells, reconstruct in y-direction, compute centered tendencies, store cell-edge state estimates
    parallel_for( Bounds<3>(nz,ny,nx) , YAKL_LAMBDA (int k, int j, int i) {
      ////////////////////////////////////////////////////////////////
      // Reconstruct rho, u, v, w, theta, pressure, dudy, and dpdy
      ////////////////////////////////////////////////////////////////
      SArray<real,2,nAder,ngll> r_DTs;
      SArray<real,2,nAder,ngll> u_DTs;
      SArray<real,2,nAder,ngll> du_DTs;
      SArray<real,2,nAder,ngll> v_DTs;
      SArray<real,2,nAder,ngll> w_DTs;
      SArray<real,2,nAder,ngll> dw_DTs;
      SArray<real,2,nAder,ngll> t_DTs;
      SArray<real,2,nAder,ngll> dt_DTs;
      SArray<real,2,nAder,ngll> rtgamma_DTs;
      SArray<real,2,nAder,ngll> dp_DTs;
      {
        SArray<real,1,ord> stencil;

        // Density
        for (int jj=0; jj < ord; jj++) { stencil(jj) = state(idR,hs+k,j+jj,hs+i); }
        reconstruct_gll_values( stencil , r_DTs , weno_scalars );
        for (int jj=0; jj < ngll; jj++) { r_DTs(0,jj) += hyDensCells(k); } // Add hydrostasis back on

        // u values and derivatives
        for (int jj=0; jj < ord; jj++) { stencil(jj) = state(idU,hs+k,j+jj,hs+i); }
        reconstruct_gll_values_and_derivs( stencil , u_DTs , du_DTs , dy , weno_winds );

        // v
        for (int jj=0; jj < ord; jj++) { stencil(jj) = state(idV,hs+k,j+jj,hs+i); }
        reconstruct_gll_values( stencil , v_DTs , weno_winds );

        // w
        for (int jj=0; jj < ord; jj++) { stencil(jj) = state(idW,hs+k,j+jj,hs+i); }
        reconstruct_gll_values_and_derivs( stencil , w_DTs , dw_DTs , dy , weno_winds );

        // theta
        for (int jj=0; jj < ord; jj++) { stencil(jj) = state(idT,hs+k,j+jj,hs+i); }
        reconstruct_gll_values_and_derivs( stencil , t_DTs , dt_DTs , dy , weno_winds );

        // pressure values and derivatives
        for (int jj=0; jj < ord; jj++) { stencil(jj) = pressure(hs+k,j+jj,hs+i); }
        reconstruct_gll_values_and_derivs( stencil , rtgamma_DTs , dp_DTs , dy , weno_scalars );
        // Add hydrostasis, and then divide by C0 to get (rho*theta)^gamma
        for (int jj=0; jj < ngll; jj++) { rtgamma_DTs(0,jj) = ( rtgamma_DTs(0,jj) + hyPressureCells(k) ) / C0; }
      }

      ///////////////////////////////////////////////////////////////
      // Compute other values needed for centered tendencies and DTs
      ///////////////////////////////////////////////////////////////
      SArray<real,2,nAder,ngll> r_v_DTs;
      SArray<real,2,nAder,ngll> v_du_DTs;
      SArray<real,2,nAder,ngll> v_v_DTs;
      SArray<real,2,nAder,ngll> rr_dp_DTs;
      SArray<real,2,nAder,ngll> v_dw_DTs;
      SArray<real,2,nAder,ngll> v_dt_DTs;
      SArray<real,2,nAder,ngll> r_t_DTs;
      for (int jj=0; jj < ngll; jj++) {
        real r  = r_DTs (0,jj);
        real u  = u_DTs (0,jj);
        real v  = v_DTs (0,jj);
        real w  = w_DTs (0,jj);
        real t  = t_DTs (0,jj);
        real du = du_DTs(0,jj);
        real dw = dw_DTs(0,jj);
        real dt = dt_DTs(0,jj);
        real dp = dp_DTs(0,jj);
        r_v_DTs  (0,jj) = r*v;
        v_du_DTs (0,jj) = v*du;
        v_v_DTs  (0,jj) = v*v;
        rr_dp_DTs(0,jj) = dp/r;
        v_dw_DTs (0,jj) = v*dw;
        v_dt_DTs (0,jj) = v*dt;
        r_t_DTs  (0,jj) = r*t;
      }

      //////////////////////////////////////////
      // Compute time derivatives if necessary
      //////////////////////////////////////////
      if (nAder > 1) {
      }

      //////////////////////////////////////////
      // Time average if necessary
      //////////////////////////////////////////
      if (timeAvg) {
      }

      //////////////////////////////////////////
      // Compute centered tendencies
      //////////////////////////////////////////
      for (int kt=0; kt < nTimeDerivs; kt++) {
        // Compute flux-form updates
        real drdt = - ( r_v_DTs(kt,ngll-1) - r_v_DTs(kt,0) ) / dy;
        real dudt = 0;
        real dvdt = - ( v_v_DTs(kt,ngll-1) - v_v_DTs(kt,0) ) / dy / 2;
        real dwdt = 0;
        real dtdt = 0;
        for (int jj=0; jj < ngll; jj++) {
          real wt = gllWts_ngll(jj);
          dudt += -v_du_DTs (kt,jj) * wt;
          dvdt += -rr_dp_DTs(kt,jj) * wt;
          dwdt += -v_dw_DTs (kt,jj) * wt;
          dtdt += -v_dt_DTs (kt,jj) * wt;
        }
        tend(idR,kt,k,j,i) = drdt;
        tend(idU,kt,k,j,i) = dudt;
        tend(idV,kt,k,j,i) = dvdt;
        tend(idW,kt,k,j,i) = dwdt;
        tend(idT,kt,k,j,i) = dtdt;
      }

      //////////////////////////////////////////
      // Store cell edge estimates of the state
      //////////////////////////////////////////
      for (int kt=0; kt < nTimeDerivs; kt++) {
        fwaves(idR,1,kt,k,j  ,i) = r_DTs(kt,0     );
        fwaves(idU,1,kt,k,j  ,i) = u_DTs(kt,0     );
        fwaves(idV,1,kt,k,j  ,i) = v_DTs(kt,0     );
        fwaves(idW,1,kt,k,j  ,i) = w_DTs(kt,0     );
        fwaves(idT,1,kt,k,j  ,i) = t_DTs(kt,0     );

        fwaves(idR,0,kt,k,j+1,i) = r_DTs(kt,ngll-1);
        fwaves(idU,0,kt,k,j+1,i) = u_DTs(kt,ngll-1);
        fwaves(idV,0,kt,k,j+1,i) = v_DTs(kt,ngll-1);
        fwaves(idW,0,kt,k,j+1,i) = w_DTs(kt,ngll-1);
        fwaves(idT,0,kt,k,j+1,i) = t_DTs(kt,ngll-1);
      }

    });

    ////////////////////////////////////////////////
    // BCs for the state edge estimates
    ////////////////////////////////////////////////
    parallel_for( Bounds<4>(numState,nTimeDerivs,nz,nx) , YAKL_LAMBDA (int l, int kt, int k, int i) {
      if        (bc_y == BC_PERIODIC) {
        fwaves(l,0,kt,k,0 ,i) = fwaves(l,0,kt,k,ny,i);
        fwaves(l,1,kt,k,ny,i) = fwaves(l,1,kt,k,0 ,i);
      } else if (bc_y == BC_WALL    ) {
        if (l == idV) {
          fwaves(l,0,kt,k,0 ,i) = 0;
          fwaves(l,1,kt,k,0 ,i) = 0;
          fwaves(l,0,kt,k,ny,i) = 0;
          fwaves(l,1,kt,k,ny,i) = 0;
        } else {
          fwaves(l,0,kt,k,0 ,i) = fwaves(l,1,kt,k,0 ,i);
          fwaves(l,1,kt,k,ny,i) = fwaves(l,0,kt,k,ny,i);
        }
      }
    });

    //////////////////////////////////////////////////////////
    // Split the flux differences into waves
    //////////////////////////////////////////////////////////
    parallel_for( Bounds<3>(nz,ny+1,nx) , YAKL_LAMBDA (int k, int j, int i) {
      // Compute interface-averaged values
      real r = 0.5_fp * ( fwaves(idR,0,0,k,j,i) + fwaves(idR,1,0,k,j,i) );
      real u = 0.5_fp * ( fwaves(idU,0,0,k,j,i) + fwaves(idU,1,0,k,j,i) );
      real v = 0.5_fp * ( fwaves(idV,0,0,k,j,i) + fwaves(idV,1,0,k,j,i) );
      real w = 0.5_fp * ( fwaves(idW,0,0,k,j,i) + fwaves(idW,1,0,k,j,i) );
      real t = 0.5_fp * ( fwaves(idT,0,0,k,j,i) + fwaves(idT,1,0,k,j,i) );
      real p = C0*pow(r*t,GAMMA);
      real cs2 = GAMMA*p/r;
      real cs  = sqrt(cs2);
      for (int kt=0; kt < nTimeDerivs; kt++) {
        // Get left and right state
        real r_L = fwaves(idR,0,kt,k,j,i);
        real u_L = fwaves(idU,0,kt,k,j,i);
        real v_L = fwaves(idV,0,kt,k,j,i);
        real w_L = fwaves(idW,0,kt,k,j,i);
        real t_L = fwaves(idT,0,kt,k,j,i);
        real p_L = C0*pow(r_L*t_L,GAMMA);
        real r_R = fwaves(idR,1,kt,k,j,i);
        real u_R = fwaves(idU,1,kt,k,j,i);
        real v_R = fwaves(idV,1,kt,k,j,i);
        real w_R = fwaves(idW,1,kt,k,j,i);
        real t_R = fwaves(idT,1,kt,k,j,i);
        real p_R = C0*pow(r_R*t_R,GAMMA);
        // Compute state jump across the interface
        real dr = r_R - r_L;
        real du = u_R - u_L;
        real dv = v_R - v_L;
        real dw = w_R - w_L;
        real dt = t_R - t_L;
        real dp = p_R - p_L;
        // Compute flux jump across the interface
        real df1 = v*dr + r*dv;
        real df2 = v*du;
        real df3 = v*dv + dp/r;
        real df4 = v*dw;
        real df5 = v*dt;
        // Compute characteristics from the jump
        real w1 = 0.5_fp*df1 - r*df3/(2*cs) + r*df5/(2*t);
        real w2 = 0.5_fp*df1 + r*df3/(2*cs) + r*df5/(2*t);
        real w3 = -r*df5/t;
        real w4 = df2;
        real w5 = df4;
        // set fwaves to zero for this interface
        for (int l=0; l < numState; l++) {
          fwaves(l,0,kt,k,j,i) = 0;
          fwaves(l,1,kt,k,j,i) = 0;
        }
        // Wave 1 (v-cs)
        fwaves(idR,0,kt,k,j,i) = w1;
        fwaves(idV,0,kt,k,j,i) = -cs*w1/r;
        // Wave 2 (v+cs)
        fwaves(idR,1,kt,k,j,i) = w2;
        fwaves(idV,1,kt,k,j,i) = cs*w2/r;
        // Waves 3-5 (v)
        if (v < 0) {
          fwaves(idR,0,kt,k,j,i) += w3;
          fwaves(idU,0,kt,k,j,i) += w4;
          fwaves(idW,0,kt,k,j,i) += w5;
          fwaves(idT,0,kt,k,j,i) += -t*w3/r;
        } else {
          fwaves(idR,1,kt,k,j,i) += w3;
          fwaves(idU,1,kt,k,j,i) += w4;
          fwaves(idW,1,kt,k,j,i) += w5;
          fwaves(idT,1,kt,k,j,i) += -t*w3/r;
        }
      }
    });

    //////////////////////////////////////////////////////////
    // Add flux waves to the tendencies
    //////////////////////////////////////////////////////////
    parallel_for( Bounds<4>(numState,nz,ny,nx) , YAKL_LAMBDA(int l, int k, int j, int i) {
      for (int kt=0; kt < nTimeDerivs; kt++) {
        tend(l,kt,k,j,i) += - ( fwaves(l,1,kt,k,j,i) + fwaves(l,0,kt,k,j+1,i) ) / dy;
      }
    });
  }



  void computeTendenciesZ( StateArr &state , TendArr &tend , real dt ) {
    // Precompute pressure perturbation
    parallel_for( Bounds<3>(nz,ny,nx) , YAKL_LAMBDA(int k, int j, int i) {
      real r = state(idR,hs+k,hs+j,hs+i) + hyDensCells(k);
      real t = state(idT,hs+k,hs+j,hs+i);
      pressure(hs+k,hs+j,hs+i) = C0*pow(r*t,GAMMA) - hyPressureCells(k);
    });

    // Populate the halos
    if        (bc_z == BC_PERIODIC) {
      parallel_for( Bounds<3>(ny,nx,hs) , YAKL_LAMBDA(int j, int i, int kk) {
        for (int l=0; l < numState; l++) {
          state(l,      kk,hs+j,hs+i) = state(l,nz+kk,hs+j,hs+i);
          state(l,hs+nz+kk,hs+j,hs+i) = state(l,hs+kk,hs+j,hs+i);
        }
        pressure(      kk,hs+j,hs+i) = pressure(nz+kk,hs+j,hs+i);
        pressure(hs+nz+kk,hs+j,hs+i) = pressure(hs+kk,hs+j,hs+i);
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
        pressure(      kk,hs+j,hs+i) = pressure(hs     ,hs+j,hs+i);
        pressure(hs+nz+kk,hs+j,hs+i) = pressure(hs+nz-1,hs+j,hs+i);
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
      SArray<real,2,nAder,ngll> t_DTs;
      SArray<real,2,nAder,ngll> dt_DTs;
      SArray<real,2,nAder,ngll> rtgamma_DTs;
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
        reconstruct_gll_values( stencil , w_DTs , weno_scalars );

        // theta
        for (int kk=0; kk < ord; kk++) { stencil(kk) = state(idT,k+kk,hs+j,hs+i); }
        reconstruct_gll_values_and_derivs( stencil , t_DTs , dt_DTs , dz , weno_winds );

        // pressure values and derivatives
        for (int kk=0; kk < ord; kk++) { stencil(kk) = pressure(k+kk,hs+j,hs+i); }
        reconstruct_gll_values_and_derivs( stencil , rtgamma_DTs , dp_DTs , dz , weno_scalars );
        // Add hydrostasis, and then divide by C0 to get (rho*theta)^gamma
        for (int kk=0; kk < ngll; kk++) { rtgamma_DTs(0,kk) = ( rtgamma_DTs(0,kk) + hyPressureGLL(k,kk) ) / C0; }
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
      SArray<real,2,nAder,ngll> r_t_DTs;
      for (int kk=0; kk < ngll; kk++) {
        real r  = r_DTs (0,kk);
        real u  = u_DTs (0,kk);
        real v  = v_DTs (0,kk);
        real w  = w_DTs (0,kk);
        real t  = t_DTs (0,kk);
        real du = du_DTs(0,kk);
        real dv = dv_DTs(0,kk);
        real dt = dt_DTs(0,kk);
        real dp = dp_DTs(0,kk);
        r_w_DTs  (0,kk) = r*w;
        w_du_DTs (0,kk) = w*du;
        w_dv_DTs (0,kk) = w*dv;
        w_w_DTs  (0,kk) = w*w;
        rr_dp_DTs(0,kk) = dp/r;
        w_dt_DTs (0,kk) = w*dt;
        r_t_DTs  (0,kk) = r*t;
      }

      //////////////////////////////////////////
      // Compute time derivatives if necessary
      //////////////////////////////////////////
      if (nAder > 1) {
      }

      //////////////////////////////////////////
      // Time average if necessary
      //////////////////////////////////////////
      if (timeAvg) {
      }

      //////////////////////////////////////////
      // Compute centered tendencies
      //////////////////////////////////////////
      for (int kt=0; kt < nTimeDerivs; kt++) {
        // Compute flux-form updates
        real drdt = - ( r_w_DTs(kt,ngll-1) - r_w_DTs(kt,0) ) / dz;
        real dudt = 0;
        real dvdt = 0;
        real dwdt = - ( w_w_DTs(kt,ngll-1) - w_w_DTs(kt,0) ) / dz / 2;
        real dtdt = 0;
        for (int kk=0; kk < ngll; kk++) {
          real wt = gllWts_ngll(kk);
          dudt += -w_du_DTs (kt,kk) * wt;
          if (! sim2d) {
            dvdt += -w_dv_DTs (kt,kk) * wt;
          }
          dwdt += -rr_dp_DTs(kt,kk) * wt;
          dtdt += -w_dt_DTs (kt,kk) * wt;
          // Source term: -g*(rho-rho_h)/rho
          if (kt == 0) {
            dwdt += ( -GRAV * (r_DTs(0,kk) - hyDensGLL(k,kk)) / r_DTs(0,kk) ) * wt;
          } else {
            dwdt += -GRAV * wt; // No time derivatives for hydrostasis
          }
        }
        tend(idR,kt,k,j,i) = drdt;
        tend(idU,kt,k,j,i) = dudt;
        tend(idV,kt,k,j,i) = dvdt;
        tend(idW,kt,k,j,i) = dwdt;
        tend(idT,kt,k,j,i) = dtdt;
      }

      //////////////////////////////////////////
      // Store cell edge estimates of the state
      //////////////////////////////////////////
      for (int kt=0; kt < nTimeDerivs; kt++) {
        fwaves(idR,1,kt,k  ,j,i) = r_DTs(kt,0     );
        fwaves(idU,1,kt,k  ,j,i) = u_DTs(kt,0     );
        fwaves(idV,1,kt,k  ,j,i) = v_DTs(kt,0     );
        fwaves(idW,1,kt,k  ,j,i) = w_DTs(kt,0     );
        fwaves(idT,1,kt,k  ,j,i) = t_DTs(kt,0     );

        fwaves(idR,0,kt,k+1,j,i) = r_DTs(kt,ngll-1);
        fwaves(idU,0,kt,k+1,j,i) = u_DTs(kt,ngll-1);
        fwaves(idV,0,kt,k+1,j,i) = v_DTs(kt,ngll-1);
        fwaves(idW,0,kt,k+1,j,i) = w_DTs(kt,ngll-1);
        fwaves(idT,0,kt,k+1,j,i) = t_DTs(kt,ngll-1);
      }

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
      real p = C0*pow(r*t,GAMMA);
      real cs2 = GAMMA*p/r;
      real cs  = sqrt(cs2);
      for (int kt=0; kt < nTimeDerivs; kt++) {
        // Get left and right state
        real r_L = fwaves(idR,0,kt,k,j,i);
        real u_L = fwaves(idU,0,kt,k,j,i);
        real v_L = fwaves(idV,0,kt,k,j,i);
        real w_L = fwaves(idW,0,kt,k,j,i);
        real t_L = fwaves(idT,0,kt,k,j,i);
        real p_L = C0*pow(r_L*t_L,GAMMA);
        real r_R = fwaves(idR,1,kt,k,j,i);
        real u_R = fwaves(idU,1,kt,k,j,i);
        real v_R = fwaves(idV,1,kt,k,j,i);
        real w_R = fwaves(idW,1,kt,k,j,i);
        real t_R = fwaves(idT,1,kt,k,j,i);
        real p_R = C0*pow(r_R*t_R,GAMMA);
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
        // Compute characteristics from the jump
        real w1 = 0.5_fp*df1 - r*df4/(2*cs) + r*df5/(2*t);
        real w2 = 0.5_fp*df1 + r*df4/(2*cs) + r*df5/(2*t);
        real w3 = -r*df5/t;
        real w4 = df2;
        real w5 = df3;
        // set fwaves to zero for this interface
        for (int l=0; l < numState; l++) {
          fwaves(l,0,kt,k,j,i) = 0;
          fwaves(l,1,kt,k,j,i) = 0;
        }
        // Wave 1 (w-cs)
        fwaves(idR,0,kt,k,j,i) = w1;
        fwaves(idW,0,kt,k,j,i) = -cs*w1/r;
        // Wave 2 (w+cs)
        fwaves(idR,1,kt,k,j,i) = w2;
        fwaves(idW,1,kt,k,j,i) = cs*w2/r;
        // Waves 3-5 (w)
        if (w < 0) {
          fwaves(idR,0,kt,k,j,i) += w3;
          fwaves(idU,0,kt,k,j,i) += w4;
          if (! sim2d) {
            fwaves(idV,0,kt,k,j,i) += w5;
          }
          fwaves(idT,0,kt,k,j,i) += -t*w3/r;
        } else {
          fwaves(idR,1,kt,k,j,i) += w3;
          fwaves(idU,1,kt,k,j,i) += w4;
          if (! sim2d) {
            fwaves(idV,1,kt,k,j,i) += w5;
          }
          fwaves(idT,1,kt,k,j,i) += -t*w3/r;
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
    parallel_for( Bounds<3>(nz,ny,nx) , YAKL_LAMBDA (int k, int j, int i) {
      real r = state(idR,hs+k,hs+j,hs+i) + hyDensCells (k);
      real t = state(idT,hs+k,hs+j,hs+i);
      real p = C0*pow(r*t,GAMMA);
      data(k,j,i) = p - hyPressureCells(k);
    });
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
