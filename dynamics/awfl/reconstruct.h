
#pragma once

#include "awfl_const.h"

namespace awfl {



  // Reconstruction data
  struct Recon {
    bool weno_scalars; // Use WENO limiting for scalars?
    bool weno_winds;   // Use WENO limiting for winds?
    // Transformation matrices for various degrees of freedom
    SArray<real,2,ord,ngll>       coefs_to_gll;
    SArray<real,2,ord,ngll>       sten_to_gll;
    SArray<real,2,ord,ord >       sten_to_coefs;
    SArray<real,3,hs+1,hs+1,hs+1> weno_recon_lower;   // WENO reconstruction matrices
    SArray<real,1,hs+2>           idl;                // Ideal weights for WENO
    real                          sigma;              // WENO sigma parameter (handicap high-order TV estimate)
    SArray<real,2,ngll,ngll>      derivMatrix;        // Transform matrix: ngll GLL pts -> ngll GLL derivs
    // Vertical grid and reconstruction matrix information
    real4d vert_sten_to_gll;
    real4d vert_sten_to_coefs;
    real5d vert_weno_recon_lower;

    void allocate_and_initialize(pam::PamCoupler const &coupler) {
      using yakl::c::parallel_for;
      using yakl::c::SimpleBounds;
      using yakl::intrinsics::matmul_cr;
      auto nens           = coupler.get_nens();
      auto nx             = coupler.get_nx();
      auto ny             = coupler.get_ny();
      auto xlen           = coupler.get_xlen();
      auto ylen           = coupler.get_ylen();
      auto nz             = coupler.get_nz();
      auto dz             = coupler.get_data_manager_readonly().get<real const,2>("vertical_cell_dz");
      auto vert_interface = coupler.get_data_manager_readonly().get<real const,2>("vertical_interface_height");

      if (coupler.option_exists("standalone_input_file")) {
        std::string inFile = coupler.get_option<std::string>( "standalone_input_file" );
        YAML::Node config = YAML::LoadFile(inFile);
        weno_scalars = config["weno_scalars"].as<bool>(true);
        weno_winds   = config["weno_winds"  ].as<bool>(true);
      } else {
        weno_scalars = true;
        weno_winds   = true;
      }

      real2d dz_ghost("dz_ghost",nz+2*hs,nens);
      parallel_for( "Spatial.h init 2" , SimpleBounds<2>(nz+2*hs,nens) , YAKL_LAMBDA (int k, int iens) {
        if      (k >= hs && k < hs+nz) { dz_ghost(k,iens) = dz(k-hs,iens); }
        else if (k < hs              ) { dz_ghost(k,iens) = dz(0   ,iens); }
        else if (k >= hs+nz          ) { dz_ghost(k,iens) = dz(nz-1,iens); }
      });

      real2d vert_interface_ghost("vert_interface_ghost",nz+2*hs+1,nens);
      parallel_for( "Spatial.h init 3" , nens , YAKL_LAMBDA (int iens) {
        vert_interface_ghost(0,iens) = vert_interface(0,iens) - hs*dz(0,iens);
        for (int k=1; k < nz+2*hs+1; k++) {
          vert_interface_ghost(k,iens) = vert_interface_ghost(k-1,iens) + dz_ghost(k-1,iens);
        }
      });

      vert_sten_to_gll      = real4d("vert_sten_to_gll"     ,nz,ord,ngll,nens);
      vert_sten_to_coefs    = real4d("vert_sten_to_coefs"   ,nz,ord,ord ,nens);
      vert_weno_recon_lower = real5d("vert_weno_recon_lower",nz,hs+1,hs+1,hs+1,nens);

      YAKL_SCOPE( vert_sten_to_gll      , this->vert_sten_to_gll      );
      YAKL_SCOPE( vert_sten_to_coefs    , this->vert_sten_to_coefs    );
      YAKL_SCOPE( vert_weno_recon_lower , this->vert_weno_recon_lower );

      auto vint_host      = vert_interface_ghost .createHostCopy();
      auto vert_s2g_host  = vert_sten_to_gll     .createHostCopy();
      auto vert_s2c_host  = vert_sten_to_coefs   .createHostCopy();
      auto vert_weno_host = vert_weno_recon_lower.createHostCopy();

      SArray<real,2,ord,ngll> c2g;
      TransformMatrices::coefs_to_gll_lower(c2g);

      for (int k=0; k < nz; k++) {
        for (int iens = 0; iens < nens; iens++) {
          // Store stencil locations
          SArray<double,1,ord+1> locs;
          for (int kk=0; kk < ord+1; kk++) {
            locs(kk) = vint_host(k+kk,iens);
          }

          // Normalize stencil locations
          double zmid = ( locs(hs+1) + locs(hs) ) / 2;
          double dzmid = locs(hs+1) - locs(hs);
          for (int kk=0; kk < ord+1; kk++) {
            locs(kk) = ( locs(kk) - zmid ) / dzmid;
          }

          // Compute reconstruction matrices
          SArray<double,2,ord,ord> s2c_var_in;
          SArray<double,3,hs+1,hs+1,hs+1> weno_recon_lower_var;
          TransformMatrices_variable::sten_to_coefs_variable<ord>(locs,s2c_var_in);
          TransformMatrices_variable::weno_lower_sten_to_coefs<ord>(locs,weno_recon_lower_var);
          SArray<real,2,ord,ord> s2c_var;
          for (int jj=0; jj < ord; jj++) {
            for (int ii=0; ii < ord; ii++) {
              s2c_var(jj,ii) = s2c_var_in(jj,ii);
            }
          }
          auto s2g_var = matmul_cr( c2g , s2c_var );

          // Store reconstruction matrices
          for (int jj=0; jj < ord; jj++) {
            for (int ii=0; ii < ord; ii++) {
              vert_s2c_host(k,jj,ii,iens) = s2c_var(jj,ii);
            }
          }

          for (int jj=0; jj < ord; jj++) {
            for (int ii=0; ii < ngll; ii++) {
              vert_s2g_host(k,jj,ii,iens) = s2g_var(jj,ii);
            }
          }

          for (int kk=0; kk < hs+1; kk++) {
            for (int jj=0; jj < hs+1; jj++) {
              for (int ii=0; ii < hs+1; ii++) {
                vert_weno_host(k,kk,jj,ii,iens) = weno_recon_lower_var(kk,jj,ii);
              }
            }
          }

        }
      }

      vert_s2g_host .deep_copy_to(vert_sten_to_gll     );
      vert_s2c_host .deep_copy_to(vert_sten_to_coefs   );
      vert_weno_host.deep_copy_to(vert_weno_recon_lower);

      // Compute the grid spacing in each dimension
      auto dx = xlen/nx;
      auto dy = ylen/ny;

      // Store the WENO reconstruction matrices
      TransformMatrices::weno_lower_sten_to_coefs(weno_recon_lower);

      // Block exists to avoid name mangling stufff
      {
        SArray<real,2,ord,ord>  s2c;        // Converts ord stencil cell averages to ord coefficients
        SArray<real,2,ord,ngll> c2g_lower;  // Converts ord coefficients to ngll GLL points
        SArray<real,2,ord,ord>  c2d;        // Converts ord coefficients to order differentiated coefficients

        TransformMatrices::sten_to_coefs     (s2c      );
        TransformMatrices::coefs_to_gll_lower(c2g_lower);

        coefs_to_gll       = c2g_lower;
        sten_to_coefs      = s2c;
        sten_to_gll        = matmul_cr( c2g_lower , s2c );
      }
      // Store ader derivMatrix
      {
        SArray<real,2,ngll,ngll> g2c;  // Converts ngll GLL points to ngll coefficients
        SArray<real,2,ngll,ngll> c2d;  // Converts ngll coefficients to ngll differentiated coefficients
        SArray<real,2,ngll,ngll> c2g;  // Converts ngll coefficients to ngll GLL points

        TransformMatrices::gll_to_coefs  (g2c);
        TransformMatrices::coefs_to_deriv(c2d);
        TransformMatrices::coefs_to_gll  (c2g);

        derivMatrix = matmul_cr( c2g , matmul_cr( c2d , g2c ) );
      }

      // Store WENO ideal weights and sigma value
      weno::wenoSetIdealSigma<ord>(idl,sigma);
    }
  };



  YAKL_INLINE void reconstruct_gll_values( SArray<real,1,ord> const &stencil                     ,
                                           SArray<real,1,ngll> &gll                              ,
                                           SArray<real,2,ord,ngll> const &coefs_to_gll           ,
                                           SArray<real,2,ord,ngll> const &sten_to_gll            ,
                                           SArray<real,2,ord,ord>  const &sten_to_coefs          ,
                                           SArray<real,3,hs+1,hs+1,hs+1> const &weno_recon_lower ,
                                           SArray<real,1,hs+2> const &idl                        ,
                                           real sigma, bool doweno ) {
    if (doweno) {

      // Reconstruct values
      SArray<real,1,ord> wenoCoefs;
      weno::compute_weno_coefs<ord>( weno_recon_lower , sten_to_coefs , stencil , wenoCoefs , idl , sigma );
      // Transform ord weno coefficients into ngll GLL points
      for (int ii=0; ii<ngll; ii++) {
        real tmp = 0;
        for (int s=0; s < ord; s++) {
          tmp += coefs_to_gll(s,ii) * wenoCoefs(s);
        }
        gll(ii) = tmp;
      }

    } else {

      // Transform ord stencil cell averages into ngll GLL points
      for (int ii=0; ii<ngll; ii++) {
        real tmp = 0;
        for (int s=0; s < ord; s++) {
          tmp += sten_to_gll(s,ii) * stencil(s);
        }
        gll(ii) = tmp;
      }

    } // if doweno
  }



}


