
program fortran_test
  use iso_c_binding
  use mmf_fortran_interface
  use gator_mod, only: gator_init, gator_finalize
  implicit none
  character(len=maxlen) :: str
  logical(c_bool), pointer, contiguous :: blah(:)
  real(c_double), pointer, contiguous :: blah2(:,:)
  real(8) :: num
  logical :: exists
  call gator_init()
  call mmf_register_dimension('myvar',100)
  call mmf_set_option('chicken','livers')
  call mmf_get_option('chicken',str)
  write(*,*) trim(str)
  call mmf_set_option('number',0.1D0)
  call mmf_get_option('number',num)
  write(*,*) num
  call mmf_option_exists('chicken',exists)
  write(*,*) exists
  call mmf_remove_option('chicken')
  call mmf_option_exists('chicken',exists)
  write(*,*) exists
  call mmf_create_array_logical("blah","blah array",[10])
  call mmf_get_array("blah",blah)
  call mmf_create_array_double("blah2","blah2 array",[10,12])
  call mmf_get_array("blah2",blah2)
  blah2 = 2.
  write(*,*) blah2
  call mmf_finalize()
  call gator_finalize()
end program

