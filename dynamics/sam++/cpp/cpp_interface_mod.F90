
module cpp_interface_mod
  use params, only: crm_rknd, crm_iknd, crm_lknd
  use iso_c_binding
  implicit none

  interface


    subroutine abcoefs() bind(C,name="abcoefs")
    end subroutine abcoefs


    subroutine adams() bind(C,name="adams")
    end subroutine adams

    
    subroutine wrap_arrays( u, v, w, t, p, tabs, qv, qcl, qpl,                                 &
                            qci, qpi, tke2, tk2, dudt, dvdt, dwdt, misc, fluxbu,               &
                            fluxbv, fluxbt, fluxbq, fluxtu, fluxtv, fluxtt, fluxtq,            &
                            fzero, precsfc, precssfc, t0, q0, qv0, tabs0, tv0,                 &
                            u0, v0, tg0, qg0, ug0, vg0, p0, tke0, t01, q01,                    &
                            qp0, qn0, prespot, rho, rhow, bet, gamaz, wsub,                    &
                            qtend, ttend, utend, vtend, sstxy, fcory, fcorzy, latitude,        &
                            longitude, prec_xy, pw_xy, cw_xy, iw_xy, cld_xy, u200_xy,          &
                            usfc_xy, v200_xy, vsfc_xy, w500_xy, w_max, u_max, twsb,            &
                            precflux, uwle, uwsb, vwle, vwsb, tkelediss, tdiff, tlat,          &
                            tlatqi, qifall, qpfall, total_water_evap, total_water_prec, CF3D,  &
                            u850_xy, v850_xy, psfc_xy, swvp_xy, cloudtopheight, echotopheight, &
                            cloudtoptemp, fcorz, fcor, longitude0, latitude0, z0, uhl,         &
                            vhl, taux0, tauy0, z, pres, zi, presi, adz, adzw,                  &
                            dt3, dz ) bind(C,name="wrap_arrays")
      use params, only: crm_rknd, crm_iknd, crm_lknd
      implicit none
      real(crm_rknd), dimension(*) :: u, v, w, t, p, tabs, qv, qcl, qpl,                                 &
                                      qci, qpi, tke2, tk2, dudt, dvdt, dwdt, misc, fluxbu,               &
                                      fluxbv, fluxbt, fluxbq, fluxtu, fluxtv, fluxtt, fluxtq,            &
                                      fzero, precsfc, precssfc, t0, q0, qv0, tabs0, tv0,                 &
                                      u0, v0, tg0, qg0, ug0, vg0, p0, tke0, t01, q01,                    &
                                      qp0, qn0, prespot, rho, rhow, bet, gamaz, wsub,                    &
                                      qtend, ttend, utend, vtend, sstxy, fcory, fcorzy, latitude,        &
                                      longitude, prec_xy, pw_xy, cw_xy, iw_xy, cld_xy, u200_xy,          &
                                      usfc_xy, v200_xy, vsfc_xy, w500_xy, w_max, u_max, twsb,            &
                                      precflux, uwle, uwsb, vwle, vwsb, tkelediss, tdiff, tlat,          &
                                      tlatqi, qifall, qpfall, total_water_evap, total_water_prec, CF3D,  &
                                      u850_xy, v850_xy, psfc_xy, swvp_xy, cloudtopheight, echotopheight, &
                                      cloudtoptemp, fcorz, fcor, longitude0, latitude0, z0, uhl,         &
                                      vhl, taux0, tauy0, z, pres, zi, presi, adz, adzw,                  &
                                      dt3, dz
    end subroutine wrap_arrays


  end interface

end module cpp_interface_mod
