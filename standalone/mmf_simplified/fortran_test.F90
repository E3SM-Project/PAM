
program fortran_test
  use iso_c_binding
  use pam_fortran_interface
  use gator_mod, only: gator_init, gator_finalize
  implicit none
  character(len=maxlen) :: str
  logical(c_bool), pointer, contiguous :: blah(:)
  real(c_double), pointer, contiguous :: blah2(:,:)
  real(8) :: myarray(50,72,40)
  real(8) :: num
  logical :: exists
  call gator_init()
  call pam_register_dimension('myvar',100)
  call pam_set_option('chicken','livers')
  call pam_get_option('chicken',str)
  write(*,*) trim(str)
  call pam_set_option('number',0.1D0)
  call pam_get_option('number',num)
  write(*,*) num
  call pam_option_exists('chicken',exists)
  write(*,*) exists
  call pam_remove_option('chicken')
  call pam_option_exists('chicken',exists)
  write(*,*) exists
  call pam_create_array_logical("blah",[10],"blah array")
  call pam_get_array("blah",blah)
  call pam_create_array_double("blah2",[10,12])
  call pam_get_array("blah2",blah2)
  call pam_mirror_array_readonly("snuffulufagus",myarray)
  call pam_mirror_array_readwrite("snuffulufagus2",myarray)
  call pam_make_readonly("blah2")
  blah2 = 2.
  write(*,*) blah2
  call pam_finalize()
  call gator_finalize()
end program

