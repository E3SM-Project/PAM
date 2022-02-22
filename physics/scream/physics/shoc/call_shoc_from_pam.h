
#pragma once

#include "ArrayIR.h"
#include "shoc_functions.hpp"

namespace pam {

  template <class T>
  inline Kokkos::View<T*,Kokkos::LayoutRight,ekat::HostDevice>
  arrayIR_to_kokkos_view( ArrayIR<T,1,IRMemHost> const &arrayIR ) {
    if (! arrayIR.initialized()) {
      yakl::yakl_throw("Error: converting non-initialized ArrayIR object into Array is not allowed");
    }
    return Kokkos::View<T*,Kokkos::LayoutRight,ekat::HostDevice>( arrayIR.get_data()  ,
                                                                  arrayIR.get_dims()[0] );
  }

  template <class T>
  inline Kokkos::View<T*,Kokkos::LayoutRight,ekat::DefaultDevice>
  arrayIR_to_kokkos_view( ArrayIR<T,1,IRMemDevice> const &arrayIR ) {
    if (! arrayIR.initialized()) {
      yakl::yakl_throw("Error: converting non-initialized ArrayIR object into Array is not allowed");
    }
    return Kokkos::View<T*,Kokkos::LayoutRight,ekat::DefaultDevice>( arrayIR.get_data()  ,
                                                                     arrayIR.get_dims()[0] );
  }

  template <class T>
  inline Kokkos::View<T**,Kokkos::LayoutRight,ekat::HostDevice>
  arrayIR_to_kokkos_view( ArrayIR<T,2,IRMemHost> const &arrayIR ) {
    if (! arrayIR.initialized()) {
      yakl::yakl_throw("Error: converting non-initialized ArrayIR object into Array is not allowed");
    }
    return Kokkos::View<T**,Kokkos::LayoutRight,ekat::HostDevice>( arrayIR.get_data()    ,
                                                                   arrayIR.get_dims()[0] ,
                                                                   arrayIR.get_dims()[1] );
  }

  template <class T>
  inline Kokkos::View<T**,Kokkos::LayoutRight,ekat::DefaultDevice>
  arrayIR_to_kokkos_view( ArrayIR<T,2,IRMemDevice> const &arrayIR ) {
    if (! arrayIR.initialized()) {
      yakl::yakl_throw("Error: converting non-initialized ArrayIR object into Array is not allowed");
    }
    return Kokkos::View<T**,Kokkos::LayoutRight,ekat::DefaultDevice>( arrayIR.get_data()    ,
                                                                      arrayIR.get_dims()[0] ,
                                                                      arrayIR.get_dims()[1] );
  }

  template <class T>
  inline Kokkos::View<T***,Kokkos::LayoutRight,ekat::HostDevice>
  arrayIR_to_kokkos_view( ArrayIR<T,3,IRMemHost> const &arrayIR ) {
    if (! arrayIR.initialized()) {
      yakl::yakl_throw("Error: converting non-initialized ArrayIR object into Array is not allowed");
    }
    return Kokkos::View<T***,Kokkos::LayoutRight,ekat::HostDevice>( arrayIR.get_data()    ,
                                                                    arrayIR.get_dims()[0] ,
                                                                    arrayIR.get_dims()[1] ,
                                                                    arrayIR.get_dims()[2]  );
  }

  template <class T>
  inline Kokkos::View<T***,Kokkos::LayoutRight,ekat::DefaultDevice>
  arrayIR_to_kokkos_view( ArrayIR<T,3,IRMemDevice> const &arrayIR ) {
    if (! arrayIR.initialized()) {
      yakl::yakl_throw("Error: converting non-initialized ArrayIR object into Array is not allowed");
    }
    return Kokkos::View<T***,Kokkos::LayoutRight,ekat::DefaultDevice>( arrayIR.get_data()    ,
                                                                       arrayIR.get_dims()[0] ,
                                                                       arrayIR.get_dims()[1] ,
                                                                       arrayIR.get_dims()[2]  );
  }


  // TODO: Change this to only accept IRMemDevice
  template <class real, int memSpace>
  void call_shoc_from_pam( int ncol, int nlev, int nlevi, real dt, int nadv, int num_qtracers,
                           ArrayIR<real,1,memSpace> shoc_host_dx    ,
                           ArrayIR<real,1,memSpace> shoc_host_dy    ,
                           ArrayIR<real,2,memSpace> shoc_thv        ,
                           ArrayIR<real,2,memSpace> shoc_zt_grid    ,
                           ArrayIR<real,2,memSpace> shoc_zi_grid    ,
                           ArrayIR<real,2,memSpace> shoc_pres       ,
                           ArrayIR<real,2,memSpace> shoc_presi      ,
                           ArrayIR<real,2,memSpace> shoc_pdel       ,
                           ArrayIR<real,1,memSpace> shoc_wthl_sfc   ,
                           ArrayIR<real,1,memSpace> shoc_wqw_sfc    ,
                           ArrayIR<real,1,memSpace> shoc_uw_sfc     ,
                           ArrayIR<real,1,memSpace> shoc_vw_sfc     ,
                           ArrayIR<real,2,memSpace> shoc_wtracer_sfc,
                           ArrayIR<real,2,memSpace> shoc_w_field    ,
                           ArrayIR<real,2,memSpace> shoc_inv_exner  ,
                           ArrayIR<real,1,memSpace> shoc_phis       ,
                           ArrayIR<real,2,memSpace> shoc_host_dse   ,
                           ArrayIR<real,2,memSpace> shoc_tke        ,
                           ArrayIR<real,2,memSpace> shoc_thetal     ,
                           ArrayIR<real,2,memSpace> shoc_qw         ,
                           ArrayIR<real,2,memSpace> shoc_u_wind     ,
                           ArrayIR<real,2,memSpace> shoc_v_wind     ,
                           ArrayIR<real,3,memSpace> shoc_qtracers   ,
                           ArrayIR<real,2,memSpace> shoc_wthv_sec   ,
                           ArrayIR<real,2,memSpace> shoc_tkh        ,
                           ArrayIR<real,2,memSpace> shoc_tk         ,
                           ArrayIR<real,2,memSpace> shoc_ql         ,
                           ArrayIR<real,2,memSpace> shoc_cldfrac    ,
                           ArrayIR<real,1,memSpace> shoc_pblh       ,
                           ArrayIR<real,2,memSpace> shoc_mix        ,
                           ArrayIR<real,2,memSpace> shoc_isotropy   ,
                           ArrayIR<real,2,memSpace> shoc_w_sec      ,
                           ArrayIR<real,2,memSpace> shoc_thl_sec    ,
                           ArrayIR<real,2,memSpace> shoc_qw_sec     ,
                           ArrayIR<real,2,memSpace> shoc_qwthl_sec  ,
                           ArrayIR<real,2,memSpace> shoc_wthl_sec   ,
                           ArrayIR<real,2,memSpace> shoc_wqw_sec    ,
                           ArrayIR<real,2,memSpace> shoc_wtke_sec   ,
                           ArrayIR<real,2,memSpace> shoc_uw_sec     ,
                           ArrayIR<real,2,memSpace> shoc_vw_sec     ,
                           ArrayIR<real,2,memSpace> shoc_w3         ,
                           ArrayIR<real,2,memSpace> shoc_wqls_sec   ,
                           ArrayIR<real,2,memSpace> shoc_brunt      ,
                           ArrayIR<real,2,memSpace> shoc_ql2        ) {
    // Create some convenient types using ekat and shoc Functions
    typedef real                                          Scalar           ;
    typedef ekat::DefaultDevice                           Device           ;
    typedef typename scream::shoc::Functions<real,Device> Functions        ;
    typedef typename Functions::Spack                     Spack            ;
    typedef typename Functions::SHOCInput                 SHOCInput        ;
    typedef typename Functions::SHOCInputOutput           SHOCInputOutput  ;
    typedef typename Functions::SHOCOutput                SHOCOutput       ;
    typedef typename Functions::SHOCHistoryOutput         SHOCHistoryOutput;

    static_assert( SCREAM_SMALL_PACK_SIZE == 1 ,
                   "ERROR: PAM's shoc integration isn't setup to deal with SCREAM_SMALL_PACK_SIZE > 1" );

    // Transform the PAM ArrayIR metadata and data pointers into Kokkos View objects of the appropriate dimensions
    auto host_dx     = arrayIR_to_kokkos_view( shoc_host_dx     );
    auto host_dy     = arrayIR_to_kokkos_view( shoc_host_dy     );
    auto thv         = arrayIR_to_kokkos_view( shoc_thv         );
    auto zt_grid     = arrayIR_to_kokkos_view( shoc_zt_grid     );
    auto zi_grid     = arrayIR_to_kokkos_view( shoc_zi_grid     );
    auto pres        = arrayIR_to_kokkos_view( shoc_pres        );
    auto presi       = arrayIR_to_kokkos_view( shoc_presi       );
    auto pdel        = arrayIR_to_kokkos_view( shoc_pdel        );
    auto wthl_sfc    = arrayIR_to_kokkos_view( shoc_wthl_sfc    );
    auto wqw_sfc     = arrayIR_to_kokkos_view( shoc_wqw_sfc     );
    auto uw_sfc      = arrayIR_to_kokkos_view( shoc_uw_sfc      );
    auto vw_sfc      = arrayIR_to_kokkos_view( shoc_vw_sfc      );
    auto wtracer_sfc = arrayIR_to_kokkos_view( shoc_wtracer_sfc );
    auto w_field     = arrayIR_to_kokkos_view( shoc_w_field     );
    auto inv_exner   = arrayIR_to_kokkos_view( shoc_inv_exner   );
    auto phis        = arrayIR_to_kokkos_view( shoc_phis        );
    auto host_dse    = arrayIR_to_kokkos_view( shoc_host_dse    );
    auto tke         = arrayIR_to_kokkos_view( shoc_tke         );
    auto thetal      = arrayIR_to_kokkos_view( shoc_thetal      );
    auto qw          = arrayIR_to_kokkos_view( shoc_qw          );
    auto u_wind      = arrayIR_to_kokkos_view( shoc_u_wind      );
    auto v_wind      = arrayIR_to_kokkos_view( shoc_v_wind      );
    auto qtracers    = arrayIR_to_kokkos_view( shoc_qtracers    );
    auto wthv_sec    = arrayIR_to_kokkos_view( shoc_wthv_sec    );
    auto tkh         = arrayIR_to_kokkos_view( shoc_tkh         );
    auto tk          = arrayIR_to_kokkos_view( shoc_tk          );
    auto ql          = arrayIR_to_kokkos_view( shoc_ql          );
    auto cldfrac     = arrayIR_to_kokkos_view( shoc_cldfrac     );
    auto pblh        = arrayIR_to_kokkos_view( shoc_pblh        );
    auto mix         = arrayIR_to_kokkos_view( shoc_mix         );
    auto isotropy    = arrayIR_to_kokkos_view( shoc_isotropy    );
    auto w_sec       = arrayIR_to_kokkos_view( shoc_w_sec       );
    auto thl_sec     = arrayIR_to_kokkos_view( shoc_thl_sec     );
    auto qw_sec      = arrayIR_to_kokkos_view( shoc_qw_sec      );
    auto qwthl_sec   = arrayIR_to_kokkos_view( shoc_qwthl_sec   );
    auto wthl_sec    = arrayIR_to_kokkos_view( shoc_wthl_sec    );
    auto wqw_sec     = arrayIR_to_kokkos_view( shoc_wqw_sec     );
    auto wtke_sec    = arrayIR_to_kokkos_view( shoc_wtke_sec    );
    auto uw_sec      = arrayIR_to_kokkos_view( shoc_uw_sec      );
    auto vw_sec      = arrayIR_to_kokkos_view( shoc_vw_sec      );
    auto w3          = arrayIR_to_kokkos_view( shoc_w3          );
    auto wqls_sec    = arrayIR_to_kokkos_view( shoc_wqls_sec    );
    auto brunt       = arrayIR_to_kokkos_view( shoc_brunt       );
    auto ql2         = arrayIR_to_kokkos_view( shoc_ql2         );

    // There appears to be no documentation in SHOC's C++ code about what the dimension ordering should be
    // for variables within the SHOC* structs. I found the following in scream's shoc_functions_f90.cpp file
    // which is the best I can gather at this point. It's pretty cumbersome having to sort through and match
    // things up.
    // std::vector<int> dim1_2d_sizes = {shcol, shcol, shcol, shcol, shcol,
    //                                   shcol, shcol, shcol, shcol, shcol,
    //                                   shcol, shcol, shcol, shcol, shcol,
    //                                   shcol, shcol, shcol, shcol,
    //                                   shcol, shcol, shcol, shcol, shcol,
    //                                   shcol, shcol, shcol, shcol, shcol,
    //                                   shcol, shcol, shcol, shcol, shcol};
    // std::vector<int> dim2_2d_sizes = {nlev,  nlevi, nlev,         nlevi, nlev,
    //                                   nlev,  nlev,  num_qtracers, nlev,  nlev,
    //                                   nlev,  nlev,  nlev,         nlev,  nlev,
    //                                   nlev,  nlev,  nlev,  nlev,
    //                                   nlev,  nlev,  nlev,         nlevi, nlevi,
    //                                   nlevi, nlevi, nlevi,        nlevi, nlevi,
    //                                   nlevi, nlevi, nlev,         nlev,  nlev};
    // std::vector<const Real*> ptr_array_1d = {host_dx, host_dy, wthl_sfc, wqw_sfc,
    //                                          uw_sfc,  vw_sfc,  phis};
    // std::vector<const Real*> ptr_array_2d = {zt_grid,   zi_grid,  pres,        presi,        pdel,
    //                                          thv,       w_field,  wtracer_sfc, inv_exner,        host_dse,
    //                                          tke,       thetal,   qw,          u_wind,       v_wind,
    //                                          wthv_sec,  tk,       shoc_cldfrac, shoc_ql,
    //                                          shoc_ql2,  shoc_mix, w_sec,       thl_sec,      qw_sec,
    //                                          qwthl_sec, wthl_sec, wqw_sec,     wtke_sec,     uw_sec,
    //                                          vw_sec,    w3,       wqls_sec,    brunt,        isotropy};
    // view_3d horiz_wind_d("horiz_wind",shcol,2,nlev_packs);
    // view_3d qtracers_cxx_d("qtracers",shcol,num_qtracers,nlev_packs);

    // Create the Views that will be used in SHOC main
    typename ekat::KokkosTypes<Device>::view_1d<Scalar> SHOCInput_dx         ( "SHOCInput_dx         " , ncol                );
    typename ekat::KokkosTypes<Device>::view_1d<Scalar> SHOCInput_dy         ( "SHOCInput_dy         " , ncol                );
    typename ekat::KokkosTypes<Device>::view_2d<Spack > SHOCInput_zt_grid    ( "SHOCInput_zt_grid    " , ncol , nlev         );
    typename ekat::KokkosTypes<Device>::view_2d<Spack > SHOCInput_zi_grid    ( "SHOCInput_zi_grid    " , ncol , nlevi        );
    typename ekat::KokkosTypes<Device>::view_2d<Spack > SHOCInput_pres       ( "SHOCInput_pres       " , ncol , nlev         );
    typename ekat::KokkosTypes<Device>::view_2d<Spack > SHOCInput_presi      ( "SHOCInput_presi      " , ncol , nlevi        );
    typename ekat::KokkosTypes<Device>::view_2d<Spack > SHOCInput_pdel       ( "SHOCInput_pdel       " , ncol , nlev         );
    typename ekat::KokkosTypes<Device>::view_2d<Spack > SHOCInput_thv        ( "SHOCInput_thv        " , ncol , nlev         );
    typename ekat::KokkosTypes<Device>::view_2d<Spack > SHOCInput_w_field    ( "SHOCInput_w_field    " , ncol , nlev         );
    typename ekat::KokkosTypes<Device>::view_1d<Scalar> SHOCInput_wthl_sfc   ( "SHOCInput_wthl_sfc   " , ncol                );
    typename ekat::KokkosTypes<Device>::view_1d<Scalar> SHOCInput_wqw_sfc    ( "SHOCInput_wqw_sfc    " , ncol                );
    typename ekat::KokkosTypes<Device>::view_1d<Scalar> SHOCInput_uw_sfc     ( "SHOCInput_uw_sfc     " , ncol                );
    typename ekat::KokkosTypes<Device>::view_1d<Scalar> SHOCInput_vw_sfc     ( "SHOCInput_vw_sfc     " , ncol                );
    typename ekat::KokkosTypes<Device>::view_2d<Spack > SHOCInput_wtracer_sfc( "SHOCInput_wtracer_sfc" , ncol , num_qtracers );
    typename ekat::KokkosTypes<Device>::view_2d<Spack > SHOCInput_inv_exner  ( "SHOCInput_inv_exner  " , ncol , nlev         );
    typename ekat::KokkosTypes<Device>::view_1d<Scalar> SHOCInput_phis       ( "SHOCInput_phis       " , ncol                );

    typename ekat::KokkosTypes<Device>::view_2d<Spack>  SHOCInputOutput_host_dse    ("SHOCInputOutput_host_dse    " , ncol                , nlev );
    typename ekat::KokkosTypes<Device>::view_2d<Spack>  SHOCInputOutput_tke         ("SHOCInputOutput_tke         " , ncol                , nlev );
    typename ekat::KokkosTypes<Device>::view_2d<Spack>  SHOCInputOutput_thetal      ("SHOCInputOutput_thetal      " , ncol                , nlev );
    typename ekat::KokkosTypes<Device>::view_2d<Spack>  SHOCInputOutput_qw          ("SHOCInputOutput_qw          " , ncol                , nlev );
    typename ekat::KokkosTypes<Device>::view_3d<Spack>  SHOCInputOutput_horiz_wind  ("SHOCInputOutput_horiz_wind  " , ncol ,            2 , nlev );
    typename ekat::KokkosTypes<Device>::view_2d<Spack>  SHOCInputOutput_wthv_sec    ("SHOCInputOutput_wthv_sec    " , ncol                , nlev );
    typename ekat::KokkosTypes<Device>::view_3d<Spack>  SHOCInputOutput_qtracers    ("SHOCInputOutput_qtracers    " , ncol , num_qtracers , nlev );
    typename ekat::KokkosTypes<Device>::view_2d<Spack>  SHOCInputOutput_tk          ("SHOCInputOutput_tk          " , ncol                , nlev );
    typename ekat::KokkosTypes<Device>::view_2d<Spack>  SHOCInputOutput_shoc_cldfrac("SHOCInputOutput_shoc_cldfrac" , ncol                , nlev );
    typename ekat::KokkosTypes<Device>::view_2d<Spack>  SHOCInputOutput_shoc_ql     ("SHOCInputOutput_shoc_ql     " , ncol                , nlev );

    typename ekat::KokkosTypes<Device>::view_1d<Scalar> SHOCOutput_pblh    ("SHOCOutput_pblh    " , ncol        );
    typename ekat::KokkosTypes<Device>::view_2d<Spack>  SHOCOutput_shoc_ql2("SHOCOutput_shoc_ql2" , ncol , nlev );

    typename ekat::KokkosTypes<Device>::view_2d<Spack>  SHOCHistoryOutput_shoc_mix ("SHOCHistoryOutput_shoc_mix " , ncol , nlev  );
    typename ekat::KokkosTypes<Device>::view_2d<Spack>  SHOCHistoryOutput_w_sec    ("SHOCHistoryOutput_w_sec    " , ncol , nlev  );
    typename ekat::KokkosTypes<Device>::view_2d<Spack>  SHOCHistoryOutput_thl_sec  ("SHOCHistoryOutput_thl_sec  " , ncol , nlevi );
    typename ekat::KokkosTypes<Device>::view_2d<Spack>  SHOCHistoryOutput_qw_sec   ("SHOCHistoryOutput_qw_sec   " , ncol , nlevi );
    typename ekat::KokkosTypes<Device>::view_2d<Spack>  SHOCHistoryOutput_qwthl_sec("SHOCHistoryOutput_qwthl_sec" , ncol , nlevi );
    typename ekat::KokkosTypes<Device>::view_2d<Spack>  SHOCHistoryOutput_wthl_sec ("SHOCHistoryOutput_wthl_sec " , ncol , nlevi );
    typename ekat::KokkosTypes<Device>::view_2d<Spack>  SHOCHistoryOutput_wqw_sec  ("SHOCHistoryOutput_wqw_sec  " , ncol , nlevi );
    typename ekat::KokkosTypes<Device>::view_2d<Spack>  SHOCHistoryOutput_wtke_sec ("SHOCHistoryOutput_wtke_sec " , ncol , nlevi );
    typename ekat::KokkosTypes<Device>::view_2d<Spack>  SHOCHistoryOutput_uw_sec   ("SHOCHistoryOutput_uw_sec   " , ncol , nlevi );
    typename ekat::KokkosTypes<Device>::view_2d<Spack>  SHOCHistoryOutput_vw_sec   ("SHOCHistoryOutput_vw_sec   " , ncol , nlevi );
    typename ekat::KokkosTypes<Device>::view_2d<Spack>  SHOCHistoryOutput_w3       ("SHOCHistoryOutput_w3       " , ncol , nlevi );
    typename ekat::KokkosTypes<Device>::view_2d<Spack>  SHOCHistoryOutput_wqls_sec ("SHOCHistoryOutput_wqls_sec " , ncol , nlev  );
    typename ekat::KokkosTypes<Device>::view_2d<Spack>  SHOCHistoryOutput_brunt    ("SHOCHistoryOutput_brunt    " , ncol , nlev  );
    typename ekat::KokkosTypes<Device>::view_2d<Spack>  SHOCHistoryOutput_isotropy ("SHOCHistoryOutput_isotropy " , ncol , nlev  );

    // Initialize the inputs
    Kokkos::parallel_for( Kokkos::MDRangePolicy<Kokkos::Rank<2>>({0,0},{ncol,nlevi}) , KOKKOS_LAMBDA (int i, int k) {
      SHOCInput_zi_grid(i,k) = zi_grid(k,i);
      SHOCInput_presi  (i,k) = presi  (k,i);
      if (k < nlev) {
        SHOCInput_zt_grid           (i  ,k) = zt_grid     (k,i);
        SHOCInput_pres              (i  ,k) = pres        (k,i);
        SHOCInput_pdel              (i  ,k) = pdel        (k,i);
        SHOCInput_thv               (i  ,k) = thv         (k,i);
        SHOCInput_w_field           (i  ,k) = w_field     (k,i);
        SHOCInput_inv_exner         (i  ,k) = inv_exner   (k,i);
        SHOCInputOutput_host_dse    (i  ,k) = host_dse    (k,i);
        SHOCInputOutput_tke         (i  ,k) = tke         (k,i);
        SHOCInputOutput_thetal      (i  ,k) = thetal      (k,i);
        SHOCInputOutput_qw          (i  ,k) = qw          (k,i);
        SHOCInputOutput_horiz_wind  (i,0,k) = u_wind      (k,i);
        SHOCInputOutput_horiz_wind  (i,1,k) = v_wind      (k,i);
        SHOCInputOutput_wthv_sec    (i  ,k) = wthv_sec    (k,i);
        SHOCInputOutput_tk          (i  ,k) = tk          (k,i);
        SHOCInputOutput_shoc_cldfrac(i  ,k) = cldfrac     (k,i);
        SHOCInputOutput_shoc_ql     (i  ,k) = ql          (k,i);
        for (int tr = 0; tr < num_qtracers; tr++) {
          SHOCInputOutput_qtracers(i,tr,k) = qtracers(tr,k,i);
        }
      }
      if (k == 0) {
        SHOCInput_dx      (i) = host_dx (i);
        SHOCInput_dy      (i) = host_dy (i);
        SHOCInput_wthl_sfc(i) = wthl_sfc(i);
        SHOCInput_wqw_sfc (i) = wqw_sfc (i);
        SHOCInput_uw_sfc  (i) = uw_sfc  (i);
        SHOCInput_vw_sfc  (i) = vw_sfc  (i);
        SHOCInput_phis    (i) = phis    (i);
        for (int tr = 0; tr < num_qtracers; tr++) {
          SHOCInput_wtracer_sfc(i,tr) = wtracer_sfc(tr,i);
        }
      }
    });

    // Create the structs used by SHOC main
    SHOCInput         shoc_in   ;
    SHOCInputOutput   shoc_inout;
    SHOCOutput        shoc_out  ;
    SHOCHistoryOutput shoc_hist ;

    // Assign the Views to the structs used by SHOC main
    shoc_in.dx          = SHOCInput_dx         ;
    shoc_in.dy          = SHOCInput_dy         ;
    shoc_in.zt_grid     = SHOCInput_zt_grid    ;
    shoc_in.zi_grid     = SHOCInput_zi_grid    ;
    shoc_in.pres        = SHOCInput_pres       ;
    shoc_in.presi       = SHOCInput_presi      ;
    shoc_in.pdel        = SHOCInput_pdel       ;
    shoc_in.thv         = SHOCInput_thv        ;
    shoc_in.w_field     = SHOCInput_w_field    ;
    shoc_in.wthl_sfc    = SHOCInput_wthl_sfc   ;
    shoc_in.wqw_sfc     = SHOCInput_wqw_sfc    ;
    shoc_in.uw_sfc      = SHOCInput_uw_sfc     ;
    shoc_in.vw_sfc      = SHOCInput_vw_sfc     ;
    shoc_in.wtracer_sfc = SHOCInput_wtracer_sfc;
    shoc_in.inv_exner   = SHOCInput_inv_exner  ;
    shoc_in.phis        = SHOCInput_phis       ;

    shoc_inout.host_dse     = SHOCInputOutput_host_dse    ;
    shoc_inout.tke          = SHOCInputOutput_tke         ;
    shoc_inout.thetal       = SHOCInputOutput_thetal      ;
    shoc_inout.qw           = SHOCInputOutput_qw          ;
    shoc_inout.horiz_wind   = SHOCInputOutput_horiz_wind  ;
    shoc_inout.wthv_sec     = SHOCInputOutput_wthv_sec    ;
    shoc_inout.qtracers     = SHOCInputOutput_qtracers    ;
    shoc_inout.tk           = SHOCInputOutput_tk          ;
    shoc_inout.shoc_cldfrac = SHOCInputOutput_shoc_cldfrac;
    shoc_inout.shoc_ql      = SHOCInputOutput_shoc_ql     ;

    shoc_out.pblh     = SHOCOutput_pblh    ;
    shoc_out.shoc_ql2 = SHOCOutput_shoc_ql2;

    shoc_hist.shoc_mix  = SHOCHistoryOutput_shoc_mix ;
    shoc_hist.w_sec     = SHOCHistoryOutput_w_sec    ;
    shoc_hist.thl_sec   = SHOCHistoryOutput_thl_sec  ;
    shoc_hist.qw_sec    = SHOCHistoryOutput_qw_sec   ;
    shoc_hist.qwthl_sec = SHOCHistoryOutput_qwthl_sec;
    shoc_hist.wthl_sec  = SHOCHistoryOutput_wthl_sec ;
    shoc_hist.wqw_sec   = SHOCHistoryOutput_wqw_sec  ;
    shoc_hist.wtke_sec  = SHOCHistoryOutput_wtke_sec ;
    shoc_hist.uw_sec    = SHOCHistoryOutput_uw_sec   ;
    shoc_hist.vw_sec    = SHOCHistoryOutput_vw_sec   ;
    shoc_hist.w3        = SHOCHistoryOutput_w3       ;
    shoc_hist.wqls_sec  = SHOCHistoryOutput_wqls_sec ;
    shoc_hist.brunt     = SHOCHistoryOutput_brunt    ;
    shoc_hist.isotropy  = SHOCHistoryOutput_isotropy ;

    // Kokkos::parallel_for( MDRangePolicy<> )

  }

}


