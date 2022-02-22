
#pragma once

#include "ArrayIR.h"
#include "shoc_functions.hpp"

namespace pam {

  typedef Kokkos::HostSpace DeviceSpace;

  template <class T>
  inline Kokkos::View<T*,Kokkos::LayoutRight,Kokkos::HostSpace>
  arrayIR_to_kokkos_view( ArrayIR<T,1,IRMemHost> const &arrayIR ) {
    if (! arrayIR.initialized()) {
      yakl::yakl_throw("Error: converting non-initialized ArrayIR object into Array is not allowed");
    }
    return Kokkos::View<T*,Kokkos::LayoutRight,Kokkos::HostSpace>( arrayIR.get_label() ,
                                                                   arrayIR.get_data()  ,
                                                                   arrayIR.get_dims()[0] );
  }

  template <class T>
  inline Kokkos::View<T*,Kokkos::LayoutRight,DeviceSpace>
  arrayIR_to_kokkos_view( ArrayIR<T,1,IRMemDevice> const &arrayIR ) {
    if (! arrayIR.initialized()) {
      yakl::yakl_throw("Error: converting non-initialized ArrayIR object into Array is not allowed");
    }
    return Kokkos::View<T*,Kokkos::LayoutRight,DeviceSpace>( arrayIR.get_label() ,
                                                             arrayIR.get_data()  ,
                                                             arrayIR.get_dims()[0] );
  }

  template <class T>
  inline Kokkos::View<T**,Kokkos::LayoutRight,Kokkos::HostSpace>
  arrayIR_to_kokkos_view( ArrayIR<T,2,IRMemHost> const &arrayIR ) {
    if (! arrayIR.initialized()) {
      yakl::yakl_throw("Error: converting non-initialized ArrayIR object into Array is not allowed");
    }
    return Kokkos::View<T**,Kokkos::LayoutRight,Kokkos::HostSpace>( arrayIR.get_label()   ,
                                                                    arrayIR.get_data()    ,
                                                                    arrayIR.get_dims()[0] ,
                                                                    arrayIR.get_dims()[1] );
  }

  template <class T>
  inline Kokkos::View<T**,Kokkos::LayoutRight,DeviceSpace>
  arrayIR_to_kokkos_view( ArrayIR<T,2,IRMemDevice> const &arrayIR ) {
    if (! arrayIR.initialized()) {
      yakl::yakl_throw("Error: converting non-initialized ArrayIR object into Array is not allowed");
    }
    return Kokkos::View<T**,Kokkos::LayoutRight,DeviceSpace>( arrayIR.get_label()   ,
                                                              arrayIR.get_data()    ,
                                                              arrayIR.get_dims()[0] ,
                                                              arrayIR.get_dims()[1] );
  }

  template <class T>
  inline Kokkos::View<T***,Kokkos::LayoutRight,Kokkos::HostSpace>
  arrayIR_to_kokkos_view( ArrayIR<T,3,IRMemHost> const &arrayIR ) {
    if (! arrayIR.initialized()) {
      yakl::yakl_throw("Error: converting non-initialized ArrayIR object into Array is not allowed");
    }
    return Kokkos::View<T***,Kokkos::LayoutRight,Kokkos::HostSpace>( arrayIR.get_label()   ,
                                                                     arrayIR.get_data()    ,
                                                                     arrayIR.get_dims()[0] ,
                                                                     arrayIR.get_dims()[1] ,
                                                                     arrayIR.get_dims()[2]  );
  }

  template <class T>
  inline Kokkos::View<T***,Kokkos::LayoutRight,DeviceSpace>
  arrayIR_to_kokkos_view( ArrayIR<T,3,IRMemDevice> const &arrayIR ) {
    if (! arrayIR.initialized()) {
      yakl::yakl_throw("Error: converting non-initialized ArrayIR object into Array is not allowed");
    }
    return Kokkos::View<T***,Kokkos::LayoutRight,DeviceSpace>( arrayIR.get_label()   ,
                                                               arrayIR.get_data()    ,
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
    int constexpr pack_size = SCREAM_SMALL_PACK_SIZE;
    typedef real                                          Scalar           ;
    typedef Kokkos::HostSpace                             Device           ;
    typedef typename scream::shoc::Functions<real,Device> Functions        ;
    typedef typename Functions::Spack                     Spack            ;
    typedef typename Functions::SHOCInput                 SHOCInput        ;
    typedef typename Functions::SHOCInputOutput           SHOCInputOutput  ;
    typedef typename Functions::SHOCOutput                SHOCOutput       ;
    typedef typename Functions::SHOCHistoryOutput         SHOCHistoryOutput;

    SHOCInput         shoc_in   ;
    SHOCInputOutput   shoc_inout;
    SHOCOutput        shoc_out  ;
    SHOCHistoryOutput shoc_hist ;

    int npacks_z = (nz-1) / pack_size + 1;  // ceiling operation for the number of packs in the z-direction

    /////////////////////////////////////////////////////////////////////
    // INPUTS
    /////////////////////////////////////////////////////////////////////
    // Functions::view_1d<const Scalar> dx         ( "dx"                         , ncol            );
    // Functions::view_1d<const Scalar> dy         ( "dy"                         , ncol            );
    // Functions::view_2d<const Spack>  zt_grid    ( "zt_grid"                    , ncol , npacks_z );
    // Functions::view_2d<const Spack>  zi_grid    ( "zi_grid"                    , ncol , npacks_z );
    // Functions::view_2d<const Spack>  pres       ( "pres"                       , ncol , npacks_z );
    // Functions::view_2d<const Spack>  presi      ( "presi"                      , ncol , npacks_z );
    // Functions::view_2d<const Spack>  pdel       ( "pdel"                       , ncol , npacks_z );
    // Functions::view_2d<const Spack>  thv        ( "thv"                        , ncol , npacks_z );
    // Functions::view_2d<const Spack>  w_field    ( "w_field"                    , ncol , npacks_z );
    // Functions::view_1d<const Scalar> wthl_sfc   ( "wthl_sfc"                   , ncol            );
    // Functions::view_1d<const Scalar> wqw_sfc    ( "wqw_sfc"                    , ncol            );
    // Functions::view_1d<const Scalar> uw_sfc     ( "uw_sfc"                     , ncol            );
    // Functions::view_1d<const Scalar> vw_sfc     ( "vw_sfc"                     , ncol            );
    // Functions::view_2d<const Spack>  wtracer_sfc( "wtracer_sfc" , num_qtracers , ncol            );  // TODO: Why is this a Spack?
    // Functions::view_2d<const Spack>  inv_exner  ( "inv_exner"                  , ncol , npacks_z );
    // Functions::view_1d<const Scalar> phis       ( "phis"                       , ncol            );


    // Kokkos::parallel_for( MDRangePolicy<> )

  }

}


