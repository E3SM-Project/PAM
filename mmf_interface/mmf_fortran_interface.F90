
module mmf_fortran_interface
  implicit none

  interface
    subroutine mmf_finalize() bind(C,name='mmf_interface_finalize')
    end subroutine
  end interface


  interface mmf_set_option
    subroutine mmf_set_option_logical(key,val) bind(C,name='mmf_interface_set_option_bool')
      use iso_c_binding
      implicit none
      character(kind=c_char) :: key(*)
      logical(c_bool), value :: val
    end subroutine
    subroutine mmf_set_option_integer(key,val) bind(C,name='mmf_interface_set_option_int')
      use iso_c_binding
      implicit none
      character(kind=c_char) :: key(*)
      integer(c_int), value  :: val
    end subroutine
    subroutine mmf_set_option_string(key,val) bind(C,name='mmf_interface_set_option_string')
      use iso_c_binding
      implicit none
      character(kind=c_char) :: key(*)
      character(kind=c_char) :: val(*)
    end subroutine
    subroutine mmf_set_option_float(key,val) bind(C,name='mmf_interface_set_option_float')
      use iso_c_binding
      implicit none
      character(kind=c_char) :: key(*)
      real(c_float), value   :: val
    end subroutine
    subroutine mmf_set_option_double(key,val) bind(C,name='mmf_interface_set_option_double')
      use iso_c_binding
      implicit none
      character(kind=c_char) :: key(*)
      real(c_double), value  :: val
    end subroutine
  end interface


  interface mmf_get_option
    subroutine mmf_get_option_logical(key,val) bind(C,name='mmf_interface_get_option_bool')
      use iso_c_binding
      implicit none
      character(kind=c_char) :: key(*)
      logical(c_bool)        :: val
    end subroutine
    subroutine mmf_get_option_integer(key,val) bind(C,name='mmf_interface_get_option_int')
      use iso_c_binding
      implicit none
      character(kind=c_char) :: key(*)
      integer(c_int)         :: val
    end subroutine
    subroutine mmf_get_option_string(key,val) bind(C,name='mmf_interface_get_option_string')
      use iso_c_binding
      implicit none
      character(kind=c_char) :: key(*)
      character(kind=c_char) :: val(*)
    end subroutine
    subroutine mmf_get_option_float(key,val) bind(C,name='mmf_interface_get_option_float')
      use iso_c_binding
      implicit none
      character(kind=c_char) :: key(*)
      real(c_float)          :: val
    end subroutine
    subroutine mmf_get_option_double(key,val) bind(C,name='mmf_interface_get_option_double')
      use iso_c_binding
      implicit none
      character(kind=c_char) :: key(*)
      real(c_double)         :: val
    end subroutine
  end interface


  interface
    subroutine mmf_option_exists(key,exists) bind(C,name='mmf_interface_option_exists')
      use iso_c_binding
      implicit none
      character(kind=c_char) :: key(*)
      logical(c_bool)        :: exists
    end subroutine
  end interface


  interface
    subroutine mmf_remove_option(key) bind(C,name='mmf_interface_remove_option')
      use iso_c_binding
      implicit none
      character(kind=c_char) :: key(*)
    end subroutine
  end interface


  interface
    subroutine mmf_create_array_logical(key,desc,dims,ndims) bind(C,name='mmf_interface_create_array_bool')
      use iso_c_binding
      implicit none
      character(kind=c_char) :: key(*), desc(*)
      integer(c_int)         :: dims(*)
      integer(c_int), value  :: ndims
    end subroutine
  end interface
  interface
    subroutine mmf_create_array_logical(key,desc,dims,ndims) bind(C,name='mmf_interface_create_array_int')
      use iso_c_binding
      implicit none
      character(kind=c_char) :: key(*), desc(*)
      integer(c_int)         :: dims(*)
      integer(c_int), value  :: ndims
    end subroutine
  end interface
  interface
    subroutine mmf_create_array_logical(key,desc,dims,ndims) bind(C,name='mmf_interface_create_array_float')
      use iso_c_binding
      implicit none
      character(kind=c_char) :: key(*), desc(*)
      integer(c_int)         :: dims(*)
      integer(c_int), value  :: ndims
    end subroutine
  end interface
  interface
    subroutine mmf_create_array_logical(key,desc,dims,ndims) bind(C,name='mmf_interface_create_array_double')
      use iso_c_binding
      implicit none
      character(kind=c_char) :: key(*), desc(*)
      integer(c_int)         :: dims(*)
      integer(c_int), value  :: ndims
    end subroutine
  end interface


  interface
    subroutine mmf_destroy_array(key) bind(C,name='mmf_interface_destroy_array')
      use iso_c_binding
      implicit none
      character(kind=c_char) :: key(*)
    end subroutine
  end interface


  interface
    subroutine mmf_array_exists(key,exists) bind(C,name='mmf_interface_array_exists')
      use iso_c_binding
      implicit none
      character(kind=c_char) :: key(*)
      logical(c_bool)        :: exists
    end subroutine
  end interface


  interface
    subroutine mmf_register_dimension(key,len) bind(C,name='mmf_interface_register_dimension')
      use iso_c_binding
      implicit none
      character(kind=c_char) :: key(*)
      integer(c_int), value  :: len
    end subroutine
  end interface

endmodule mmf_fortran_interface

