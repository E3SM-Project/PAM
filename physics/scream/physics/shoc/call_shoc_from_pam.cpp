
#include "call_shoc_from_pam.h"
#include "shoc_functions.hpp"
#include "ArrayIR_kokkos.h"

namespace pam {

  // Computes and returns nbpl
  int call_shoc_init_from_pam (int nbot_shoc, int ntop_shoc, ArrayIR<double,1,IRMemDevice> input_pref_mid) {
    typedef ekat::DefaultDevice                             Device;
    typedef typename scream::shoc::Functions<double,Device> SHOC  ;
    typedef typename SHOC::Spack                            Spack ;

    int nlev = input_pref_mid.get_dims()[0];

    auto pref_mid = arrayIR_to_kokkos_view( input_pref_mid );

    typename ekat::KokkosTypes<Device>::view_1d<Spack> pref_spack("pref_spack",nlev);

    Kokkos::parallel_for( Kokkos::RangePolicy<>(0,nlev) , KOKKOS_LAMBDA (int k) {
      pref_spack(k) = pref_mid(k);
    });

    return SHOC::shoc_init( nbot_shoc , ntop_shoc , pref_spack );
  }


  void call_shoc_main_from_pam( int ncol, int nlev, int nlevi, double dtime, int nadv, int num_qtracers, int npbl,
                                ArrayIR<double,1,IRMemDevice> input_host_dx    ,
                                ArrayIR<double,1,IRMemDevice> input_host_dy    ,
                                ArrayIR<double,2,IRMemDevice> input_thv        ,
                                ArrayIR<double,2,IRMemDevice> input_zt_grid    ,
                                ArrayIR<double,2,IRMemDevice> input_zi_grid    ,
                                ArrayIR<double,2,IRMemDevice> input_pres       ,
                                ArrayIR<double,2,IRMemDevice> input_presi      ,
                                ArrayIR<double,2,IRMemDevice> input_pdel       ,
                                ArrayIR<double,1,IRMemDevice> input_wthl_sfc   ,
                                ArrayIR<double,1,IRMemDevice> input_wqw_sfc    ,
                                ArrayIR<double,1,IRMemDevice> input_uw_sfc     ,
                                ArrayIR<double,1,IRMemDevice> input_vw_sfc     ,
                                ArrayIR<double,2,IRMemDevice> input_wtracer_sfc,
                                ArrayIR<double,2,IRMemDevice> input_w_field    ,
                                ArrayIR<double,2,IRMemDevice> input_inv_exner  ,
                                ArrayIR<double,1,IRMemDevice> input_phis       ,
                                ArrayIR<double,2,IRMemDevice> input_host_dse   ,
                                ArrayIR<double,2,IRMemDevice> input_tke        ,
                                ArrayIR<double,2,IRMemDevice> input_thetal     ,
                                ArrayIR<double,2,IRMemDevice> input_qw         ,
                                ArrayIR<double,2,IRMemDevice> input_u_wind     ,
                                ArrayIR<double,2,IRMemDevice> input_v_wind     ,
                                ArrayIR<double,3,IRMemDevice> input_qtracers   ,
                                ArrayIR<double,2,IRMemDevice> input_wthv_sec   ,
                                ArrayIR<double,2,IRMemDevice> input_tkh        ,
                                ArrayIR<double,2,IRMemDevice> input_tk         ,
                                ArrayIR<double,2,IRMemDevice> input_ql         ,
                                ArrayIR<double,2,IRMemDevice> input_cldfrac    ,
                                ArrayIR<double,1,IRMemDevice> input_pblh       ,
                                ArrayIR<double,2,IRMemDevice> input_mix        ,
                                ArrayIR<double,2,IRMemDevice> input_isotropy   ,
                                ArrayIR<double,2,IRMemDevice> input_w_sec      ,
                                ArrayIR<double,2,IRMemDevice> input_thl_sec    ,
                                ArrayIR<double,2,IRMemDevice> input_qw_sec     ,
                                ArrayIR<double,2,IRMemDevice> input_qwthl_sec  ,
                                ArrayIR<double,2,IRMemDevice> input_wthl_sec   ,
                                ArrayIR<double,2,IRMemDevice> input_wqw_sec    ,
                                ArrayIR<double,2,IRMemDevice> input_wtke_sec   ,
                                ArrayIR<double,2,IRMemDevice> input_uw_sec     ,
                                ArrayIR<double,2,IRMemDevice> input_vw_sec     ,
                                ArrayIR<double,2,IRMemDevice> input_w3         ,
                                ArrayIR<double,2,IRMemDevice> input_wqls_sec   ,
                                ArrayIR<double,2,IRMemDevice> input_brunt      ,
                                ArrayIR<double,2,IRMemDevice> input_ql2        ) {
    // Create some convenient types using ekat and shoc Functions
    typedef double                                          Scalar           ;
    typedef ekat::DefaultDevice                             Device           ;
    typedef typename scream::shoc::Functions<double,Device> SHOC             ;
    typedef typename SHOC::Spack                            Spack            ;
    typedef typename SHOC::SHOCInput                        SHOCInput        ;
    typedef typename SHOC::SHOCInputOutput                  SHOCInputOutput  ;
    typedef typename SHOC::SHOCOutput                       SHOCOutput       ;
    typedef typename SHOC::SHOCHistoryOutput                SHOCHistoryOutput;
    typedef typename SHOC::WorkspaceMgr                     WorkspaceMgr     ;

    static_assert( SCREAM_SMALL_PACK_SIZE == 1 ,
                   "ERROR: PAM's shoc integration isn't setup to deal with SCREAM_SMALL_PACK_SIZE > 1" );

    // Transform the PAM ArrayIR metadata and data pointers into Kokkos View objects of the appropriate dimensions
    auto host_dx     = arrayIR_to_kokkos_view( input_host_dx     );
    auto host_dy     = arrayIR_to_kokkos_view( input_host_dy     );
    auto thv         = arrayIR_to_kokkos_view( input_thv         );
    auto zt_grid     = arrayIR_to_kokkos_view( input_zt_grid     );
    auto zi_grid     = arrayIR_to_kokkos_view( input_zi_grid     );
    auto pres        = arrayIR_to_kokkos_view( input_pres        );
    auto presi       = arrayIR_to_kokkos_view( input_presi       );
    auto pdel        = arrayIR_to_kokkos_view( input_pdel        );
    auto wthl_sfc    = arrayIR_to_kokkos_view( input_wthl_sfc    );
    auto wqw_sfc     = arrayIR_to_kokkos_view( input_wqw_sfc     );
    auto uw_sfc      = arrayIR_to_kokkos_view( input_uw_sfc      );
    auto vw_sfc      = arrayIR_to_kokkos_view( input_vw_sfc      );
    auto wtracer_sfc = arrayIR_to_kokkos_view( input_wtracer_sfc );
    auto w_field     = arrayIR_to_kokkos_view( input_w_field     );
    auto inv_exner   = arrayIR_to_kokkos_view( input_inv_exner   );
    auto phis        = arrayIR_to_kokkos_view( input_phis        );
    auto host_dse    = arrayIR_to_kokkos_view( input_host_dse    );
    auto tke         = arrayIR_to_kokkos_view( input_tke         );
    auto thetal      = arrayIR_to_kokkos_view( input_thetal      );
    auto qw          = arrayIR_to_kokkos_view( input_qw          );
    auto u_wind      = arrayIR_to_kokkos_view( input_u_wind      );
    auto v_wind      = arrayIR_to_kokkos_view( input_v_wind      );
    auto qtracers    = arrayIR_to_kokkos_view( input_qtracers    );
    auto wthv_sec    = arrayIR_to_kokkos_view( input_wthv_sec    );
    auto tkh         = arrayIR_to_kokkos_view( input_tkh         );
    auto tk          = arrayIR_to_kokkos_view( input_tk          );
    auto ql          = arrayIR_to_kokkos_view( input_ql          );
    auto cldfrac     = arrayIR_to_kokkos_view( input_cldfrac     );
    auto pblh        = arrayIR_to_kokkos_view( input_pblh        );
    auto mix         = arrayIR_to_kokkos_view( input_mix         );
    auto isotropy    = arrayIR_to_kokkos_view( input_isotropy    );
    auto w_sec       = arrayIR_to_kokkos_view( input_w_sec       );
    auto thl_sec     = arrayIR_to_kokkos_view( input_thl_sec     );
    auto qw_sec      = arrayIR_to_kokkos_view( input_qw_sec      );
    auto qwthl_sec   = arrayIR_to_kokkos_view( input_qwthl_sec   );
    auto wthl_sec    = arrayIR_to_kokkos_view( input_wthl_sec    );
    auto wqw_sec     = arrayIR_to_kokkos_view( input_wqw_sec     );
    auto wtke_sec    = arrayIR_to_kokkos_view( input_wtke_sec    );
    auto uw_sec      = arrayIR_to_kokkos_view( input_uw_sec      );
    auto vw_sec      = arrayIR_to_kokkos_view( input_vw_sec      );
    auto w3          = arrayIR_to_kokkos_view( input_w3          );
    auto wqls_sec    = arrayIR_to_kokkos_view( input_wqls_sec    );
    auto brunt       = arrayIR_to_kokkos_view( input_brunt       );
    auto ql2         = arrayIR_to_kokkos_view( input_ql2         );

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
        SHOCInput_zt_grid           (i  ,k) = zt_grid  (k,i);
        SHOCInput_pres              (i  ,k) = pres     (k,i);
        SHOCInput_pdel              (i  ,k) = pdel     (k,i);
        SHOCInput_thv               (i  ,k) = thv      (k,i);
        SHOCInput_w_field           (i  ,k) = w_field  (k,i);
        SHOCInput_inv_exner         (i  ,k) = inv_exner(k,i);
        SHOCInputOutput_host_dse    (i  ,k) = host_dse (k,i);
        SHOCInputOutput_tke         (i  ,k) = tke      (k,i);
        SHOCInputOutput_thetal      (i  ,k) = thetal   (k,i);
        SHOCInputOutput_qw          (i  ,k) = qw       (k,i);
        SHOCInputOutput_horiz_wind  (i,0,k) = u_wind   (k,i);
        SHOCInputOutput_horiz_wind  (i,1,k) = v_wind   (k,i);
        SHOCInputOutput_wthv_sec    (i  ,k) = wthv_sec (k,i);
        SHOCInputOutput_tk          (i  ,k) = tk       (k,i);
        SHOCInputOutput_shoc_cldfrac(i  ,k) = cldfrac  (k,i);
        SHOCInputOutput_shoc_ql     (i  ,k) = ql       (k,i);
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

    Kokkos::parallel_for( Kokkos::MDRangePolicy<Kokkos::Rank<2>>({0,0},{ncol,nlevi}) , KOKKOS_LAMBDA (int i, int k) {
      SHOCHistoryOutput_thl_sec  (i,k) = 0;
      SHOCHistoryOutput_qw_sec   (i,k) = 0;
      SHOCHistoryOutput_qwthl_sec(i,k) = 0;
      SHOCHistoryOutput_wthl_sec (i,k) = 0;
      SHOCHistoryOutput_wqw_sec  (i,k) = 0;
      SHOCHistoryOutput_wtke_sec (i,k) = 0;
      SHOCHistoryOutput_uw_sec   (i,k) = 0;
      SHOCHistoryOutput_vw_sec   (i,k) = 0;
      SHOCHistoryOutput_w3       (i,k) = 0;
      if (k < nlev) {
        SHOCOutput_shoc_ql2         (i  ,k) = 0;
        SHOCHistoryOutput_shoc_mix  (i  ,k) = 0;
        SHOCHistoryOutput_w_sec     (i  ,k) = 0;
        SHOCHistoryOutput_wqls_sec  (i  ,k) = 0;
        SHOCHistoryOutput_brunt     (i  ,k) = 0;
        SHOCHistoryOutput_isotropy  (i  ,k) = 0;
        for (int tr = 0; tr < num_qtracers; tr++) {
          SHOCInputOutput_qtracers(i,tr,k) = 0;
        }
      }
      if (k == 0) {
        SHOCOutput_pblh(i) = 0;
      }
    });

    Kokkos::fence();

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

    // Call shoc_main
    auto const nlev_packs   = ekat::npack<Spack>(nlev);
    auto const policy       = ekat::ExeSpaceUtils<ekat::KokkosTypes<Device>::ExeSpace>::get_default_team_policy(ncol,nlev_packs);
    auto const nlevi_packs  = ekat::npack<Spack>(nlevi);
    int  const n_wind_slots = ekat::npack<Spack>(2)*Spack::n;
    int  const n_trac_slots = ekat::npack<Spack>(num_qtracers+3)*Spack::n;
    WorkspaceMgr workspace_mgr(nlevi_packs, 13+(n_wind_slots+n_trac_slots), policy);

    // static Int shoc_main(
    //   const Int&               shcol,                // Number of SHOC columns in the array
    //   const Int&               nlev,                 // Number of levels
    //   const Int&               nlevi,                // Number of levels on interface grid
    //   const Int&               npbl,                 // Maximum number of levels in pbl from surface
    //   const Int&               nadv,                 // Number of times to loop SHOC
    //   const Int&               num_q_tracers,        // Number of tracers
    //   const Scalar&            dtime,                // SHOC timestep [s]
    //   const WorkspaceMgr&      workspace_mgr,        // WorkspaceManager for local variables
    //   const SHOCInput&         shoc_input,           // Input
    //   const SHOCInputOutput&   shoc_input_output,    // Input/Output
    //   const SHOCOutput&        shoc_output,          // Output
    //   const SHOCHistoryOutput& shoc_history_output); // Output (diagnostic)
    const auto elapsed_microsec = SHOC::shoc_main(ncol, nlev, nlevi, npbl, nadv, num_qtracers, dtime,
                                                  workspace_mgr, shoc_in, shoc_inout, shoc_out, shoc_hist);

    Kokkos::fence();

    // Copy the outputs
    Kokkos::parallel_for( Kokkos::MDRangePolicy<Kokkos::Rank<2>>({0,0},{ncol,nlevi}) , KOKKOS_LAMBDA (int i, int k) {
      thl_sec  (k,i) = SHOCHistoryOutput_thl_sec  (i,k)[0];
      qw_sec   (k,i) = SHOCHistoryOutput_qw_sec   (i,k)[0];
      qwthl_sec(k,i) = SHOCHistoryOutput_qwthl_sec(i,k)[0];
      wthl_sec (k,i) = SHOCHistoryOutput_wthl_sec (i,k)[0];
      wqw_sec  (k,i) = SHOCHistoryOutput_wqw_sec  (i,k)[0];
      wtke_sec (k,i) = SHOCHistoryOutput_wtke_sec (i,k)[0];
      uw_sec   (k,i) = SHOCHistoryOutput_uw_sec   (i,k)[0];
      vw_sec   (k,i) = SHOCHistoryOutput_vw_sec   (i,k)[0];
      w3       (k,i) = SHOCHistoryOutput_w3       (i,k)[0];
      if (k < nlev) {
        host_dse    (k,i) = SHOCInputOutput_host_dse    (i  ,k)[0];
        tke         (k,i) = SHOCInputOutput_tke         (i  ,k)[0];
        thetal      (k,i) = SHOCInputOutput_thetal      (i  ,k)[0];
        qw          (k,i) = SHOCInputOutput_qw          (i  ,k)[0];
        u_wind      (k,i) = SHOCInputOutput_horiz_wind  (i,0,k)[0];
        v_wind      (k,i) = SHOCInputOutput_horiz_wind  (i,1,k)[0];
        wthv_sec    (k,i) = SHOCInputOutput_wthv_sec    (i  ,k)[0];
        tk          (k,i) = SHOCInputOutput_tk          (i  ,k)[0];
        cldfrac     (k,i) = SHOCInputOutput_shoc_cldfrac(i  ,k)[0];
        ql          (k,i) = SHOCInputOutput_shoc_ql     (i  ,k)[0];
        ql2         (k,i) = SHOCOutput_shoc_ql2         (i  ,k)[0];
        mix         (k,i) = SHOCHistoryOutput_shoc_mix  (i  ,k)[0];
        w_sec       (k,i) = SHOCHistoryOutput_w_sec     (i  ,k)[0];
        wqls_sec    (k,i) = SHOCHistoryOutput_wqls_sec  (i  ,k)[0];
        brunt       (k,i) = SHOCHistoryOutput_brunt     (i  ,k)[0];
        isotropy    (k,i) = SHOCHistoryOutput_isotropy  (i  ,k)[0];
        for (int tr = 0; tr < num_qtracers; tr++) {
          qtracers(tr,k,i) = SHOCInputOutput_qtracers(i,tr,k)[0];
        }
      }
      if (k == 0) {
        pblh(i) = SHOCOutput_pblh(i);
      }
    });

    Kokkos::fence();

  }

}


