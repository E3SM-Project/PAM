
module supercell_init_mod
  use iso_c_binding
  implicit none
  
  interface
    subroutine supercell_init(vert_interface, rho_d_col, uvel_col, vvel_col, wvel_col, temp_col, rho_v_col, &
                              Rd, Rv, grav, nz) bind(C,name="supercell_init_fortran")
      use iso_c_binding
      real(c_double) :: vert_interface(*)
      real(c_double) :: rho_d_col     (*)
      real(c_double) :: uvel_col      (*)
      real(c_double) :: vvel_col      (*)
      real(c_double) :: wvel_col      (*)
      real(c_double) :: temp_col      (*)
      real(c_double) :: rho_v_col     (*)
      real(c_double) :: Rd, Rv, grav
      integer(c_int) :: nz
    end subroutine
  end interface

end module

