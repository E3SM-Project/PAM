
module cpp_interface_mod
  use iso_c_binding
  implicit none

  interface

    ! extern "C" void abcoefs(real *dt3);
    subroutine abcoefs(dt3) bind(C,name="abcoefs")
      use iso_c_binding
      use params, only: crm_rknd, crm_iknd, crm_lknd
      implicit none
      real(crm_rknd), intent(in   ) :: dt3(3)
    end subroutine abcoefs

  end interface

end module cpp_interface_mod
