#pragma once

#include "common.h"
#include "TransformMatrices.h"
#include "TransformMatrices_variable.h"



  template <int ord>
  YAKL_INLINE void map_weights2( SArray<real,1,(ord-1)/2+2> const &idl , SArray<real,1,(ord-1)/2+2> &wts ) {
    int constexpr hs = (ord-1)/2;
    // Map the weights for quicker convergence. WARNING: Ideal weights must be (0,1) before mapping
    for (int i=0; i<hs+2; i++) {
      wts(i) = wts(i) * ( idl(i) + idl(i)*idl(i) - 3._fp*idl(i)*wts(i) + wts(i)*wts(i) ) /
                        ( idl(i)*idl(i) + wts(i) * ( 1._fp - 2._fp * idl(i) ) );
    }
  }


  template <int ord>
  YAKL_INLINE void convexify2( SArray<real,1,(ord-1)/2+2> &wts ) {
    int constexpr hs = (ord-1)/2;
    real sum = 0._fp;
    real const eps = 1.0e-20;
    for (int i=0; i<hs+2; i++) { sum += wts(i); }
    for (int i=0; i<hs+2; i++) { wts(i) /= (sum + eps); }
  }

  template <int ord>
  YAKL_INLINE void compute_weno_coefs2( SArray<real,3,(ord-1)/2+1,(ord-1)/2+1,(ord-1)/2+1> const &recon_lo ,
                                       SArray<real,2,ord,ord> const & recon_hi ,
                                       SArray<real,1,ord> const &u ,
                                       SArray<real,1,ord> &aw ,
                                       SArray<real,1,(ord-1)/2+2> const &idl ,
                                       real const sigma ) {
    int constexpr hs = (ord-1)/2;
    SArray<real,2,hs+1,hs+1> a_lo;
    SArray<real,1,ord> a_hi;
    real const eps = 1.0e-20;

    // Compute three quadratic polynomials (left, center, and right) and the high-order polynomial
    for(int i=0; i<hs+1; i++) {
      for (int ii=0; ii<hs+1; ii++) {
        real tmp = 0;
        for (int s=0; s<hs+1; s++) {
          tmp += recon_lo(i,s,ii) * u(i+s);
        }
        a_lo(i,ii) = tmp;
      }
    }
    for (int ii=0; ii<ord; ii++) {
      real tmp = 0;
      for (int s=0; s<ord; s++) {
        tmp += recon_hi(s,ii) * u(s);
      }
      a_hi(ii) = tmp;
    }

    // Compute "bridge" polynomial
    for (int i=0; i<hs+1; i++) {
      for (int ii=0; ii<hs+1; ii++) {
        a_hi(ii) -= idl(i)*a_lo(i,ii);
      }
    }
    for (int ii=0; ii<ord; ii++) {
      a_hi(ii) /= idl(hs+1);
    }

    SArray<real,1,hs+1> lotmp;
    SArray<real,1,hs+2> tv;

    // Compute total variation of all candidate polynomials
    for (int i=0; i<hs+1; i++) {
      for (int ii=0; ii<hs+1; ii++) {
        lotmp(ii) = a_lo(i,ii);
      }
      tv(i) = TransformMatrices::coefs_to_tv(lotmp);
    }
    tv(hs+1) = TransformMatrices::coefs_to_tv(a_hi);

    real lo_avg;

    // Reduce the bridge polynomial TV to something closer to the other TV values
    lo_avg = 0._fp;
    for (int i=0; i<hs+1; i++) {
      lo_avg += tv(i);
    }
    lo_avg /= hs+1;
    tv(hs+1) = lo_avg + ( tv(hs+1) - lo_avg ) * sigma;

    SArray<real,1,hs+2> wts;

    // WENO weights are proportional to the inverse of TV**2 and then re-confexified
    for (int i=0; i<hs+2; i++) {
      wts(i) = idl(i) / ( tv(i)*tv(i) + eps );
    }
    convexify2<ord>(wts);

    // Map WENO weights for sharper fronts and less sensitivity to "eps"
    map_weights2<ord>(idl,wts);
    convexify2<ord>(wts);

    // WENO polynomial is the weighted sum of candidate polynomials using WENO weights instead of ideal weights
    for (int i=0; i < ord; i++) {
      aw(i) = wts(hs+1) * a_hi(i);
    }
    for (int i=0; i<hs+1; i++) {
      for (int ii=0; ii<hs+1; ii++) {
        aw(ii) += wts(i) * a_lo(i,ii);
      }
    }
  }



 template <uint ndofs, uint ord, uint hs = (ord - 1) / 2>
 void YAKL_INLINE weno_func_vert(SArray<real, 2, ndofs, 2> &edgerecon,
                            SArray<real, 2, ndofs, ord> const &dens,
                            //ADD MORE STUFF HERE!
                            SArray<real,2,ord,2> const &coefs_to_gll           ,
                            SArray<real,2,ord,2> const &sten_to_gll            ,
                            SArray<real,2,ord,ord>  const &sten_to_coefs          ,
                            SArray<real,3,hs+1,hs+1,hs+1> const &weno_recon_lower ,
                            SArray<real, 1, hs + 2> const &wenoIdl,
                            real wenoSigma) {

   SArray<real, 1, ord> stencil;
   SArray<real, 1, 2> gllPts;
   SArray<real,1,ord> wenoCoefs;

   for (int l = 0; l < ndofs; l++) {

       for (int ii = 0; ii < ord; ii++) {
         stencil(ii) = dens(l, ii);
       }

       // Reconstruct values
       compute_weno_coefs2<ord>( weno_recon_lower , sten_to_coefs , stencil , wenoCoefs , wenoIdl , wenoSigma );
       // Transform ord weno coefficients into ngll GLL points
       for (int ii=0; ii<2; ii++) {
         real tmp = 0;
         for (int s=0; s < ord; s++) {
           tmp += coefs_to_gll(s,ii) * wenoCoefs(s);
         }
         gllPts(ii) = tmp;
       }

       edgerecon(l, 0) = gllPts(0);
       edgerecon(l, 1) = gllPts(1);
     }
   }




   template <class T, uint ord, uint hs=(ord-1)/2>
   void create_variable_WENO(
     real4d coefs_to_gll_arr,
     real4d sten_to_gll_arr,
     real4d sten_to_coefs_arr,
     real5d weno_recon_lower_arr,
   real &wenoSigma,
   SArray<real, 1, hs+2> &wenoIdl,
   const Geometry<T> &geom
   )
   {
     using yakl::intrinsics::matmul_cr;

     const auto &topo = geom.topology;

     wenoSetIdealSigma<ord>(wenoIdl,wenoSigma);

     parallel_for("Compute Vertically Variable WENO Func Arrays", SimpleBounds<2>(topo.nl,topo.nens),
         YAKL_LAMBDA(int k, int n) {

           SArray<real,2,ord,2> coefs_to_gll;
           SArray<real,2,ord,2> sten_to_gll;
           SArray<real,2,ord,ord> sten_to_coefs;
           SArray<real,3,hs+1,hs+1,hs+1> weno_recon_lower;

           TransformMatrices::coefs_to_gll_lower(coefs_to_gll);

         // Store stencil locations
         SArray<double,1,ord+1> locs;
         for (int kk=0; kk < ord+1; kk++) {
           locs(kk) = geom.zint(k+kk,n);
         }

        //std::cout << "locs at k = " << k << "\n";
        //for (int kk=0; kk < ord+1; kk++) {std::cout << locs(kk) << "\n";}

           // Normalize stencil locations
           real zmid = ( locs(hs+1) + locs(hs) ) / 2;
           real dzmid = locs(hs+1) - locs(hs);
           for (int kk=0; kk < ord+1; kk++) {
             locs(kk) = ( locs(kk) - zmid ) / dzmid;
           }

           // Compute reconstruction matrices
           SArray<real,2,ord,ord> s2c_var_in;
           SArray<real,3,hs+1,hs+1,hs+1> weno_recon_lower_var;
           TransformMatrices_variable::sten_to_coefs_variable<ord>(locs,s2c_var_in);
           TransformMatrices_variable::weno_lower_sten_to_coefs<ord>(locs,weno_recon_lower);
           SArray<real,2,ord,ord> s2c_var;
           for (int jj=0; jj < ord; jj++) {
             for (int ii=0; ii < ord; ii++) {
               sten_to_coefs(jj,ii) = s2c_var_in(jj,ii);
             }
           }
           sten_to_gll = matmul_cr( coefs_to_gll , sten_to_coefs );


     for (int h=0;h<ord;h++){
       for (int g=0;g<2;g++){
         coefs_to_gll_arr(k,h,g,n) = coefs_to_gll(h,g);
         sten_to_gll_arr(k,h,g,n) = sten_to_gll(h,g);
       }}

     for (int h1=0;h1<ord;h1++){
       for (int h2=0;h2<ord;h2++){
         sten_to_coefs_arr(k,h1,h2,n) = sten_to_coefs(h1,h2);
       }}


     for (int h1=0;h1<hs+1;h1++){
       for (int h2=0;h2<hs+1;h2++){
         for (int h3=0;h3<hs+1;h3++){
           weno_recon_lower_arr(k,h1,h2,h3,n) = weno_recon_lower(h1,h2,h3);
       }}}

       });

   }
