
module mmf_fortran_interface
  use iso_c_binding
  implicit none

  integer, parameter :: maxlen = 256

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! Interfaces for Fortran-facing routines
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  interface mmf_set_option
    module procedure mmf_set_option_logical
    module procedure mmf_set_option_integer
    module procedure mmf_set_option_string
    module procedure mmf_set_option_float
    module procedure mmf_set_option_double
  end interface

  interface mmf_get_option
    module procedure mmf_get_option_logical
    module procedure mmf_get_option_integer
    module procedure mmf_get_option_string
    module procedure mmf_get_option_float
    module procedure mmf_get_option_double
  end interface



  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! Interfaces for extern "C" routines
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  interface
    subroutine mmf_finalize() bind(C,name='mmf_interface_finalize')
    end subroutine
  end interface


  interface
    subroutine mmf_set_option_logical_c(key,val) bind(C,name='mmf_interface_set_option_bool')
      use iso_c_binding
      implicit none
      character(kind=c_char) :: key(*)
      logical(c_bool), value :: val
    end subroutine
    subroutine mmf_set_option_integer_c(key,val) bind(C,name='mmf_interface_set_option_int')
      use iso_c_binding
      implicit none
      character(kind=c_char) :: key(*)
      integer(c_int), value  :: val
    end subroutine
    subroutine mmf_set_option_string_c(key,val) bind(C,name='mmf_interface_set_option_string')
      use iso_c_binding
      implicit none
      character(kind=c_char) :: key(*)
      character(kind=c_char) :: val(*)
    end subroutine
    subroutine mmf_set_option_float_c(key,val) bind(C,name='mmf_interface_set_option_float')
      use iso_c_binding
      implicit none
      character(kind=c_char) :: key(*)
      real(c_float), value   :: val
    end subroutine
    subroutine mmf_set_option_double_c(key,val) bind(C,name='mmf_interface_set_option_double')
      use iso_c_binding
      implicit none
      character(kind=c_char) :: key(*)
      real(c_double), value  :: val
    end subroutine
  end interface


  interface
    subroutine mmf_get_option_logical_c(key,val) bind(C,name='mmf_interface_get_option_bool')
      use iso_c_binding
      implicit none
      character(kind=c_char) :: key(*)
      logical(c_bool)        :: val
    end subroutine
    subroutine mmf_get_option_integer_c(key,val) bind(C,name='mmf_interface_get_option_int')
      use iso_c_binding
      implicit none
      character(kind=c_char) :: key(*)
      integer(c_int)         :: val
    end subroutine
    subroutine mmf_get_option_stringlen_c(key,len) bind(C,name='mmf_interface_get_option_stringlen')
      use iso_c_binding
      implicit none
      character(kind=c_char) :: key(*)
      integer(c_int)         :: len
    end subroutine
    subroutine mmf_get_option_string_c(key,val) bind(C,name='mmf_interface_get_option_string')
      use iso_c_binding
      implicit none
      character(kind=c_char) :: key(*)
      character(kind=c_char) :: val(*)
    end subroutine
    subroutine mmf_get_option_float_c(key,val) bind(C,name='mmf_interface_get_option_float')
      use iso_c_binding
      implicit none
      character(kind=c_char) :: key(*)
      real(c_float)          :: val
    end subroutine
    subroutine mmf_get_option_double_c(key,val) bind(C,name='mmf_interface_get_option_double')
      use iso_c_binding
      implicit none
      character(kind=c_char) :: key(*)
      real(c_double)         :: val
    end subroutine
  end interface


  interface
    subroutine mmf_option_exists_c(key,exists) bind(C,name='mmf_interface_option_exists')
      use iso_c_binding
      implicit none
      character(kind=c_char) :: key(*)
      logical(c_bool)        :: exists
    end subroutine
  end interface


  interface
    subroutine mmf_remove_option_c(key) bind(C,name='mmf_interface_remove_option')
      use iso_c_binding
      implicit none
      character(kind=c_char) :: key(*)
    end subroutine
  end interface


  interface
    subroutine mmf_create_array_logical_c(key,desc,dims,ndims) bind(C,name='mmf_interface_create_array_bool')
      use iso_c_binding
      implicit none
      character(kind=c_char) :: key(*), desc(*)
      integer(c_int)         :: dims(*)
      integer(c_int), value  :: ndims
    end subroutine
  end interface
  interface
    subroutine mmf_create_array_int_c(key,desc,dims,ndims) bind(C,name='mmf_interface_create_array_int')
      use iso_c_binding
      implicit none
      character(kind=c_char) :: key(*), desc(*)
      integer(c_int)         :: dims(*)
      integer(c_int), value  :: ndims
    end subroutine
  end interface
  interface
    subroutine mmf_create_array_float_c(key,desc,dims,ndims) bind(C,name='mmf_interface_create_array_float')
      use iso_c_binding
      implicit none
      character(kind=c_char) :: key(*), desc(*)
      integer(c_int)         :: dims(*)
      integer(c_int), value  :: ndims
    end subroutine
  end interface
  interface
    subroutine mmf_create_array_double_c(key,desc,dims,ndims) bind(C,name='mmf_interface_create_array_double')
      use iso_c_binding
      implicit none
      character(kind=c_char) :: key(*), desc(*)
      integer(c_int)         :: dims(*)
      integer(c_int), value  :: ndims
    end subroutine
  end interface


  interface
    subroutine mmf_get_array_bool_c(label,ptr,dims,ndims)
      use iso_c_binding
      implicit none
      character(kind=c_char) :: label(*)
      logical(c_bool)        :: ptr(*)
      integer(c_int)         :: dims(*)
      integer(c_int), value  :: ndims
    end subroutine
    subroutine mmf_get_array_int_c(label,ptr,dims,ndims)
      use iso_c_binding
      implicit none
      character(kind=c_char) :: label(*)
      integer(c_int)         :: ptr(*)
      integer(c_int)         :: dims(*)
      integer(c_int), value  :: ndims
    end subroutine
    subroutine mmf_get_array_float_c(label,ptr,dims,ndims)
      use iso_c_binding
      implicit none
      character(kind=c_char) :: label(*)
      real(c_float)          :: ptr(*)
      integer(c_int)         :: dims(*)
      integer(c_int), value  :: ndims
    end subroutine
    subroutine mmf_get_array_double_c(label,ptr,dims,ndims)
      use iso_c_binding
      implicit none
      character(kind=c_char) :: label(*)
      real(c_double)         :: ptr(*)
      integer(c_int)         :: dims(*)
      integer(c_int), value  :: ndims
    end subroutine
  end interface


  interface
    subroutine mmf_destroy_array_c(key) bind(C,name='mmf_interface_destroy_array')
      use iso_c_binding
      implicit none
      character(kind=c_char) :: key(*)
    end subroutine
  end interface


  interface
    subroutine mmf_array_exists_c(key,exists) bind(C,name='mmf_interface_array_exists')
      use iso_c_binding
      implicit none
      character(kind=c_char) :: key(*)
      logical(c_bool)        :: exists
    end subroutine
  end interface


  interface
    subroutine mmf_register_dimension_c(key,len) bind(C,name='mmf_interface_register_dimension')
      use iso_c_binding
      implicit none
      character(kind=c_char) :: key(*)
      integer(c_int), value  :: len
    end subroutine
  end interface

contains

  function string_f2c(fort,len)
    implicit none
    character(len=*), intent(in) :: fort
    integer         , intent(in) :: len
    character(kind=c_char)       :: string_f2c(len+1)
    integer :: i
    do i = 1 , len
      string_f2c(i) = fort(i:i)
    enddo
    string_f2c(len+1) = char(0)
  end function
  function string_c2f(c,len)
    implicit none
    character(c_char), intent(in) :: c(*)
    integer(c_int)   , intent(in) :: len
    character(len=len)            :: string_c2f
    integer :: i
    do i = 1 , len
      string_c2f(i:i) = c(i)
    enddo
  end function


  subroutine mmf_register_dimension(key,len)
    implicit none
    character(len=*), intent(in) :: key
    integer         , intent(in) :: len
    call mmf_register_dimension_c( string_f2c(key,len_trim(key)) , int(len,c_int) )
  end subroutine


  subroutine mmf_set_option_logical(key,val)
    implicit none
    character(len=*), intent(in) :: key
    logical         , intent(in) :: val
    call mmf_set_option_logical_c( string_f2c(key,len_trim(key)) , logical(val,c_bool) )
  end subroutine
  subroutine mmf_set_option_integer(key,val)
    implicit none
    character(len=*), intent(in) :: key
    integer         , intent(in) :: val
    call mmf_set_option_integer_c( string_f2c(key,len_trim(key)) , int(val,c_int) )
  end subroutine
  subroutine mmf_set_option_string(key,val)
    implicit none
    character(len=*), intent(in) :: key
    character(len=*), intent(in) :: val
    call mmf_set_option_string_c( string_f2c(key,len_trim(key)) , string_f2c(val,len_trim(val)) )
  end subroutine
  subroutine mmf_set_option_float(key,val)
    implicit none
    character(len=*), intent(in) :: key
    real(4)         , intent(in) :: val
    call mmf_set_option_float_c( string_f2c(key,len_trim(key)) , real(val,c_float) )
  end subroutine
  subroutine mmf_set_option_double(key,val)
    implicit none
    character(len=*), intent(in) :: key
    real(8)         , intent(in) :: val
    call mmf_set_option_double_c( string_f2c(key,len_trim(key)) , real(val,c_double) )
  end subroutine


  subroutine mmf_get_option_logical(key,val)
    implicit none
    character(len=*), intent(in   ) :: key
    logical         , intent(  out) :: val
    logical(c_bool) :: val_c
    call mmf_get_option_logical_c( string_f2c(key,len_trim(key)) , val_c )
    val = val_c
  end subroutine
  subroutine mmf_get_option_integer(key,val)
    implicit none
    character(len=*), intent(in   ) :: key
    integer         , intent(  out) :: val
    integer(c_int) :: val_c
    call mmf_get_option_integer_c( string_f2c(key,len_trim(key)) , val_c )
    val = val_c
  end subroutine
  subroutine mmf_get_option_string(key,val)
    implicit none
    character(len=*)     , intent(in   ) :: key
    character(len=maxlen), intent(  out) :: val
    character(kind=c_char), allocatable :: val_c(:)
    integer(c_int) :: len
    call mmf_get_option_stringlen_c( string_f2c(key,len_trim(key)) , len )
    allocate(val_c(len))
    call mmf_get_option_string_c( string_f2c(key,len_trim(key)) , val_c )
    val = string_c2f( val_c , len )
    deallocate(val_c)
  end subroutine
  subroutine mmf_get_option_float(key,val)
    implicit none
    character(len=*), intent(in   ) :: key
    real(4)         , intent(  out) :: val
    real(c_float) :: val_c
    call mmf_get_option_float_c( string_f2c(key,len_trim(key)) , val_c )
    val = val_c
  end subroutine
  subroutine mmf_get_option_double(key,val)
    implicit none
    character(len=*), intent(in   ) :: key
    real(8)         , intent(  out) :: val
    real(c_double) :: val_c
    call mmf_get_option_double_c( string_f2c(key,len_trim(key)) , val_c )
    val = val_c
  end subroutine


endmodule mmf_fortran_interface

