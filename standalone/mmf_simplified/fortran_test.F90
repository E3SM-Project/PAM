
program fortran_test
  use mmf_fortran_interface
  implicit none
  character(len=maxlen) :: str
  real(8) :: num
  call mmf_register_dimension('myvar',100)
  call mmf_set_option('chicken','livers')
  call mmf_get_option('chicken',str)
  write(*,*) trim(str)
  call mmf_set_option('number',0.1D0)
  call mmf_get_option('number',num)
  write(*,*) num
end program

