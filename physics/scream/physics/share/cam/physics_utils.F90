! physics_utils: common elements for SCREAM physics models.
module physics_utils
  use iso_c_binding, only: c_double, c_float, c_bool
  implicit none
  private
  save
  integer,parameter,public :: rtype8 = c_double ! 8 byte real, compatible with c type double
  integer,parameter,public :: btype  = c_bool ! boolean type, compatible with c
  integer,parameter,public :: itype = selected_int_kind (13) ! 8 byte integer
  integer,parameter,public :: rtype = c_double ! 8 byte real, compatible with c type double
end module physics_utils

