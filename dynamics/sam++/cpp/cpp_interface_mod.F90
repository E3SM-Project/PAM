
module cpp_interface_mod
  use iso_c_binding
  implicit none

  interface

    ! extern "C" void abcoefs(int na, int nb, int nc, int nstep, real const &dt3, real &at, real &bt, real &ct);
    subroutine abcoefs(na, nb, nc, nstep, dt3, at, bt, ct) bind(C,name="abcoefs")
      use iso_c_binding
      use params, only: crm_rknd, crm_iknd, crm_lknd
      implicit none
      integer(crm_iknd), value, intent(in   ) :: na, nb, nc, nstep
      real   (crm_rknd)       , intent(in   ) :: dt3(3)
      real   (crm_rknd)       , intent(  out) :: at, bt, ct
    end subroutine abcoefs

  end interface

end module cpp_interface_mod
