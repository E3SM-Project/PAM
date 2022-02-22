
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


  template <class real, int memSpace>
  void call_shoc_from_pam( int ncol, int nz, int nzp1, real dt, int nadv, int num_qtracers,
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
    typedef real                                          Scalar           ;
    typedef ekat::DefaultDevice                           Device           ;
    typedef typename scream::shoc::Functions<real,Device> Functions        ;
    typedef typename Functions::Spack                     Spack            ;
    typedef typename Functions::SHOCInput                 SHOCInput        ;
    typedef typename Functions::SHOCInputOutput           SHOCInputOutput  ;
    typedef typename Functions::SHOCOutput                SHOCOutput       ;
    typedef typename Functions::SHOCHistoryOutput         SHOCHistoryOutput;

    int constexpr pack_size = SCREAM_SMALL_PACK_SIZE;

    SHOCInput         shoc_in   ;
    SHOCInputOutput   shoc_inout;
    SHOCOutput        shoc_out  ;
    SHOCHistoryOutput shoc_hist ;

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

    int npacks_z = (nz-1) / pack_size + 1;  // ceiling operation for the number of packs in the z-direction

    /////////////////////////////////////////////////////////////////////
    // INPUTS
    /////////////////////////////////////////////////////////////////////
    typename ekat::KokkosTypes<Device>::view_1d<Scalar> SHOCInput_dx         ( "dx"                         , ncol            );
    typename ekat::KokkosTypes<Device>::view_1d<Scalar> SHOCInput_dy         ( "dy"                         , ncol            );
    typename ekat::KokkosTypes<Device>::view_2d<Spack > SHOCInput_zt_grid    ( "zt_grid"                    , ncol , npacks_z );
    typename ekat::KokkosTypes<Device>::view_2d<Spack > SHOCInput_zi_grid    ( "zi_grid"                    , ncol , npacks_z );
    typename ekat::KokkosTypes<Device>::view_2d<Spack > SHOCInput_pres       ( "pres"                       , ncol , npacks_z );
    typename ekat::KokkosTypes<Device>::view_2d<Spack > SHOCInput_presi      ( "presi"                      , ncol , npacks_z );
    typename ekat::KokkosTypes<Device>::view_2d<Spack > SHOCInput_pdel       ( "pdel"                       , ncol , npacks_z );
    typename ekat::KokkosTypes<Device>::view_2d<Spack > SHOCInput_thv        ( "thv"                        , ncol , npacks_z );
    typename ekat::KokkosTypes<Device>::view_2d<Spack > SHOCInput_w_field    ( "w_field"                    , ncol , npacks_z );
    typename ekat::KokkosTypes<Device>::view_1d<Scalar> SHOCInput_wthl_sfc   ( "wthl_sfc"                   , ncol            );
    typename ekat::KokkosTypes<Device>::view_1d<Scalar> SHOCInput_wqw_sfc    ( "wqw_sfc"                    , ncol            );
    typename ekat::KokkosTypes<Device>::view_1d<Scalar> SHOCInput_uw_sfc     ( "uw_sfc"                     , ncol            );
    typename ekat::KokkosTypes<Device>::view_1d<Scalar> SHOCInput_vw_sfc     ( "vw_sfc"                     , ncol            );
    typename ekat::KokkosTypes<Device>::view_2d<Spack > SHOCInput_wtracer_sfc( "wtracer_sfc" , num_qtracers , ncol            );  // TODO: Why is this a Spack?
    typename ekat::KokkosTypes<Device>::view_2d<Spack > SHOCInput_inv_exner  ( "inv_exner"                  , ncol , npacks_z );
    typename ekat::KokkosTypes<Device>::view_1d<Scalar> SHOCInput_phis       ( "phis"                       , ncol            );






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

    // Kokkos::parallel_for( MDRangePolicy<> )

  }

}


