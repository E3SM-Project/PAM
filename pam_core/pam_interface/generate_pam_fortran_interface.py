
import os, stat

if (os.path.isfile("pam_fortran_interface.F90")) :
  os.chmod("pam_fortran_interface.F90", stat.S_IRUSR | stat.S_IWUSR)

f = open("pam_fortran_interface.F90", "w")

option_types = ["logical","integer","string","float","double"]
array_types = ["logical","integer","float","double"]
array_iso_c_map = {"logical":"logical(c_bool)" , "integer":"integer(c_int)" , "float":"real(c_float)" , "double":"real(c_double)"}

f.write('''

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!  WARNING: THIS FILE IS GENERATED. DO NOT EDIT DIRECTLY
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module pam_fortran_interface
  use iso_c_binding
  implicit none 
  integer, parameter :: maxlen = 256

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! Interfaces for Fortran-facing routines
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  interface pam_finalize
    module procedure pam_finalize
  end interface

  interface pam_register_dimension
    module procedure pam_register_dimension
  end interface

  interface pam_option_exists
    module procedure pam_option_exists
  end interface

  interface pam_remove_option
    module procedure pam_remove_option
  end interface

  interface pam_destroy_array
    module procedure pam_destroy_array
  end interface

  interface pam_array_exists
    module procedure pam_array_exists
  end interface

  interface pam_set_option\n''')

for tp in option_types :
  f.write(f"    module procedure pam_set_option_{tp}\n")

f.write('''  end interface

  interface pam_get_option\n''')

for tp in option_types :
  f.write(f"    module procedure pam_get_option_{tp}\n")
f.write('''  end interface

  interface pam_get_array\n''')

for tp in array_types :
  for d in range(1,8) :
    f.write(f"    module procedure pam_get_array_{tp}_{d}d\n")

f.write('''  end interface

  interface pam_mirror_array_readonly\n''')

for tp in array_types :
  for d in range(1,8) :
    f.write(f"    module procedure pam_mirror_array_readonly_{tp}_{d}d\n")

f.write('''  end interface

  interface pam_mirror_array_readwrite\n''')

for tp in array_types :
  for d in range(1,8) :
    f.write(f"    module procedure pam_mirror_array_readwrite_{tp}_{d}d\n")

f.write('''  end interface\n\n''')

f.write('''
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! Interfaces for extern "C" routines
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  interface
    subroutine pam_finalize_c() &
               bind(C,name='pam_interface_finalize')
    end subroutine
  end interface

  interface
    subroutine pam_set_option_logical_c(key,val) &
               bind(C,name='pam_interface_set_option_bool')
      use iso_c_binding
      implicit none
      character(kind=c_char) :: key(*)
      logical(c_bool), value :: val
    end subroutine
    subroutine pam_set_option_integer_c(key,val) &
               bind(C,name='pam_interface_set_option_int')
      use iso_c_binding
      implicit none
      character(kind=c_char) :: key(*)
      integer(c_int), value  :: val
    end subroutine
    subroutine pam_set_option_string_c(key,val) &
               bind(C,name='pam_interface_set_option_string')
      use iso_c_binding
      implicit none
      character(kind=c_char) :: key(*)
      character(kind=c_char) :: val(*)
    end subroutine
    subroutine pam_set_option_float_c(key,val) &
               bind(C,name='pam_interface_set_option_float')
      use iso_c_binding
      implicit none
      character(kind=c_char) :: key(*)
      real(c_float), value   :: val
    end subroutine
    subroutine pam_set_option_double_c(key,val) &
               bind(C,name='pam_interface_set_option_double')
      use iso_c_binding
      implicit none
      character(kind=c_char) :: key(*)
      real(c_double), value  :: val
    end subroutine
  end interface

  interface
    subroutine pam_get_option_logical_c(key,val) &
               bind(C,name='pam_interface_get_option_bool')
      use iso_c_binding
      implicit none
      character(kind=c_char) :: key(*)
      logical(c_bool)        :: val
    end subroutine
    subroutine pam_get_option_integer_c(key,val) &
               bind(C,name='pam_interface_get_option_int')
      use iso_c_binding
      implicit none
      character(kind=c_char) :: key(*)
      integer(c_int)         :: val
    end subroutine
    subroutine pam_get_option_stringlen_c(key,len) &
               bind(C,name='pam_interface_get_option_stringlen')
      use iso_c_binding
      implicit none
      character(kind=c_char) :: key(*)
      integer(c_int)         :: len
    end subroutine
    subroutine pam_get_option_string_c(key,val) &
               bind(C,name='pam_interface_get_option_string')
      use iso_c_binding
      implicit none
      character(kind=c_char) :: key(*)
      character(kind=c_char) :: val(*)
    end subroutine
    subroutine pam_get_option_float_c(key,val) &
               bind(C,name='pam_interface_get_option_float')
      use iso_c_binding
      implicit none
      character(kind=c_char) :: key(*)
      real(c_float)          :: val
    end subroutine
    subroutine pam_get_option_double_c(key,val) &
               bind(C,name='pam_interface_get_option_double')
      use iso_c_binding
      implicit none
      character(kind=c_char) :: key(*)
      real(c_double)         :: val
    end subroutine
  end interface

  interface
    subroutine pam_option_exists_c(key,exists) &
               bind(C,name='pam_interface_option_exists')
      use iso_c_binding
      implicit none
      character(kind=c_char) :: key(*)
      logical(c_bool)        :: exists
    end subroutine
  end interface

  interface
    subroutine pam_remove_option_c(key) &
               bind(C,name='pam_interface_remove_option')
      use iso_c_binding
      implicit none
      character(kind=c_char) :: key(*)
    end subroutine
    subroutine pam_make_readonly_c(key) &
               bind(C,name='pam_interface_make_readonly')
      use iso_c_binding
      implicit none
      character(kind=c_char) :: key(*)
    end subroutine
  end interface

  interface
    subroutine pam_create_array_logical_c(key,desc,dims,ndims) &
               bind(C,name='pam_interface_create_array_bool')
      use iso_c_binding
      implicit none
      character(kind=c_char) :: key(*), desc(*)
      integer(c_int)         :: dims(*)
      integer(c_int), value  :: ndims
    end subroutine
  end interface
  interface
    subroutine pam_create_array_integer_c(key,desc,dims,ndims) &
               bind(C,name='pam_interface_create_array_int')
      use iso_c_binding
      implicit none
      character(kind=c_char) :: key(*), desc(*)
      integer(c_int)         :: dims(*)
      integer(c_int), value  :: ndims
    end subroutine
  end interface
  interface
    subroutine pam_create_array_float_c(key,desc,dims,ndims) &
               bind(C,name='pam_interface_create_array_float')
      use iso_c_binding
      implicit none
      character(kind=c_char) :: key(*), desc(*)
      integer(c_int)         :: dims(*)
      integer(c_int), value  :: ndims
    end subroutine
  end interface
  interface
    subroutine pam_create_array_double_c(key,desc,dims,ndims) &
               bind(C,name='pam_interface_create_array_double')
      use iso_c_binding
      implicit none
      character(kind=c_char) :: key(*), desc(*)
      integer(c_int)         :: dims(*)
      integer(c_int), value  :: ndims
    end subroutine
  end interface

  interface
    subroutine pam_mirror_array_readonly_logical_c(key,desc,dims,ndims,arr) &
               bind(C,name='pam_interface_mirror_array_readonly_bool')
      use iso_c_binding
      implicit none
      character(kind=c_char) :: key(*), desc(*)
      integer(c_int)         :: dims(*)
      integer(c_int), value  :: ndims
      logical(c_bool)        :: arr(*)
    end subroutine
  end interface
  interface
    subroutine pam_mirror_array_readonly_integer_c(key,desc,dims,ndims,arr) &
               bind(C,name='pam_interface_mirror_array_readonly_int')
      use iso_c_binding
      implicit none
      character(kind=c_char) :: key(*), desc(*)
      integer(c_int)         :: dims(*)
      integer(c_int), value  :: ndims
      integer(c_int)         :: arr(*)
    end subroutine
  end interface
  interface
    subroutine pam_mirror_array_readonly_float_c(key,desc,dims,ndims,arr) &
               bind(C,name='pam_interface_mirror_array_readonly_float')
      use iso_c_binding
      implicit none
      character(kind=c_char) :: key(*), desc(*)
      integer(c_int)         :: dims(*)
      integer(c_int), value  :: ndims
      real(c_float)          :: arr(*)
    end subroutine
  end interface
  interface
    subroutine pam_mirror_array_readonly_double_c(key,desc,dims,ndims,arr) &
               bind(C,name='pam_interface_mirror_array_readonly_double')
      use iso_c_binding
      implicit none
      character(kind=c_char) :: key(*), desc(*)
      integer(c_int)         :: dims(*)
      integer(c_int), value  :: ndims
      real(c_double)         :: arr(*)
    end subroutine
  end interface

  interface
    subroutine pam_mirror_array_readwrite_logical_c(key,desc,dims,ndims,arr) &
               bind(C,name='pam_interface_mirror_array_readwrite_bool')
      use iso_c_binding
      implicit none
      character(kind=c_char) :: key(*), desc(*)
      integer(c_int)         :: dims(*)
      integer(c_int), value  :: ndims
      logical(c_bool)        :: arr(*)
    end subroutine
  end interface
  interface
    subroutine pam_mirror_array_readwrite_integer_c(key,desc,dims,ndims,arr) &
               bind(C,name='pam_interface_mirror_array_readwrite_int')
      use iso_c_binding
      implicit none
      character(kind=c_char) :: key(*), desc(*)
      integer(c_int)         :: dims(*)
      integer(c_int), value  :: ndims
      integer(c_int)         :: arr(*)
    end subroutine
  end interface
  interface
    subroutine pam_mirror_array_readwrite_float_c(key,desc,dims,ndims,arr) &
               bind(C,name='pam_interface_mirror_array_readwrite_float')
      use iso_c_binding
      implicit none
      character(kind=c_char) :: key(*), desc(*)
      integer(c_int)         :: dims(*)
      integer(c_int), value  :: ndims
      real(c_float)          :: arr(*)
    end subroutine
  end interface
  interface
    subroutine pam_mirror_array_readwrite_double_c(key,desc,dims,ndims,arr) &
               bind(C,name='pam_interface_mirror_array_readwrite_double')
      use iso_c_binding
      implicit none
      character(kind=c_char) :: key(*), desc(*)
      integer(c_int)         :: dims(*)
      integer(c_int), value  :: ndims
      real(c_double)         :: arr(*)
    end subroutine
  end interface

  interface
    subroutine pam_get_array_logical_c(label,ptr,dims,ndims) &
               bind(C,name='pam_interface_get_array_bool')
      use iso_c_binding
      implicit none
      character(kind=c_char) :: label(*)
      type(c_ptr)            :: ptr
      integer(c_int)         :: dims(*)
      integer(c_int), value  :: ndims
    end subroutine
    subroutine pam_get_array_integer_c(label,ptr,dims,ndims) &
               bind(C,name='pam_interface_get_array_int')
      use iso_c_binding
      implicit none
      character(kind=c_char) :: label(*)
      type(c_ptr)            :: ptr
      integer(c_int)         :: dims(*)
      integer(c_int), value  :: ndims
    end subroutine
    subroutine pam_get_array_float_c(label,ptr,dims,ndims) &
               bind(C,name='pam_interface_get_array_float')
      use iso_c_binding
      implicit none
      character(kind=c_char) :: label(*)
      type(c_ptr)            :: ptr
      integer(c_int)         :: dims(*)
      integer(c_int), value  :: ndims
    end subroutine
    subroutine pam_get_array_double_c(label,ptr,dims,ndims) &
               bind(C,name='pam_interface_get_array_double')
      use iso_c_binding
      implicit none
      character(kind=c_char) :: label(*)
      type(c_ptr)            :: ptr
      integer(c_int)         :: dims(*)
      integer(c_int), value  :: ndims
    end subroutine
  end interface

  interface
    subroutine pam_destroy_array_c(key) &
               bind(C,name='pam_interface_destroy_array')
      use iso_c_binding
      implicit none
      character(kind=c_char) :: key(*)
    end subroutine
  end interface

  interface
    subroutine pam_array_exists_c(key,exists) &
               bind(C,name='pam_interface_array_exists')
      use iso_c_binding
      implicit none
      character(kind=c_char) :: key(*)
      logical(c_bool)        :: exists
    end subroutine
  end interface

  interface
    subroutine pam_register_dimension_c(key,len) &
               bind(C,name='pam_interface_register_dimension')
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


  subroutine pam_finalize()
    call pam_finalize_c()
  end subroutine


  subroutine pam_register_dimension(key,len)
    implicit none
    character(len=*), intent(in) :: key
    integer         , intent(in) :: len
    call pam_register_dimension_c( string_f2c(key,len_trim(key)) , int(len,c_int) )
  end subroutine


  subroutine pam_option_exists(key,exists)
    implicit none
    character(len=*), intent(in   ) :: key
    logical         , intent(  out) :: exists
    logical(c_bool) :: exists_c
    call pam_option_exists_c( string_f2c(key,len_trim(key)) , exists_c )
    exists = exists_c
  end subroutine


  subroutine pam_remove_option(key)
    implicit none
    character(len=*), intent(in   ) :: key
    call pam_remove_option_c( string_f2c(key,len_trim(key)) )
  end subroutine


  subroutine pam_make_readonly(key)
    implicit none
    character(len=*), intent(in   ) :: key
    call pam_make_readonly_c( string_f2c(key,len_trim(key)) )
  end subroutine


  subroutine pam_destroy_array(key)
    implicit none
    character(len=*), intent(in) :: key
    call pam_destroy_array_c( string_f2c(key,len_trim(key)) )
  end subroutine


  subroutine pam_array_exists(key,exists)
    implicit none
    character(len=*), intent(in   ) :: key
    logical         , intent(  out) :: exists
    logical(c_bool) :: exists_c
    call pam_array_exists_c( string_f2c(key,len_trim(key)) , exists_c )
    exists = exists_c
  end subroutine


  subroutine pam_set_option_logical(key,val)
    implicit none
    character(len=*), intent(in) :: key
    logical         , intent(in) :: val
    call pam_set_option_logical_c( string_f2c(key,len_trim(key)) , logical(val,c_bool) )
  end subroutine
  subroutine pam_set_option_integer(key,val)
    implicit none
    character(len=*), intent(in) :: key
    integer         , intent(in) :: val
    call pam_set_option_integer_c( string_f2c(key,len_trim(key)) , int(val,c_int) )
  end subroutine
  subroutine pam_set_option_string(key,val)
    implicit none
    character(len=*), intent(in) :: key
    character(len=*), intent(in) :: val
    call pam_set_option_string_c( string_f2c(key,len_trim(key)) , string_f2c(val,len_trim(val)) )
  end subroutine
  subroutine pam_set_option_float(key,val)
    implicit none
    character(len=*), intent(in) :: key
    real(4)         , intent(in) :: val
    call pam_set_option_float_c( string_f2c(key,len_trim(key)) , real(val,c_float) )
  end subroutine
  subroutine pam_set_option_double(key,val)
    implicit none
    character(len=*), intent(in) :: key
    real(8)         , intent(in) :: val
    call pam_set_option_double_c( string_f2c(key,len_trim(key)) , real(val,c_double) )
  end subroutine


  subroutine pam_get_option_logical(key,val)
    implicit none
    character(len=*), intent(in   ) :: key
    logical         , intent(  out) :: val
    logical(c_bool) :: val_c
    call pam_get_option_logical_c( string_f2c(key,len_trim(key)) , val_c )
    val = val_c
  end subroutine
  subroutine pam_get_option_integer(key,val)
    implicit none
    character(len=*), intent(in   ) :: key
    integer         , intent(  out) :: val
    integer(c_int) :: val_c
    call pam_get_option_integer_c( string_f2c(key,len_trim(key)) , val_c )
    val = val_c
  end subroutine
  subroutine pam_get_option_string(key,val)
    implicit none
    character(len=*)     , intent(in   ) :: key
    character(len=maxlen), intent(  out) :: val
    character(kind=c_char), allocatable :: val_c(:)
    integer(c_int) :: len
    call pam_get_option_stringlen_c( string_f2c(key,len_trim(key)) , len )
    allocate(val_c(len))
    call pam_get_option_string_c( string_f2c(key,len_trim(key)) , val_c )
    val = string_c2f( val_c , len )
    deallocate(val_c)
  end subroutine
  subroutine pam_get_option_float(key,val)
    implicit none
    character(len=*), intent(in   ) :: key
    real(4)         , intent(  out) :: val
    real(c_float) :: val_c
    call pam_get_option_float_c( string_f2c(key,len_trim(key)) , val_c )
    val = val_c
  end subroutine
  subroutine pam_get_option_double(key,val)
    implicit none
    character(len=*), intent(in   ) :: key
    real(8)         , intent(  out) :: val
    real(c_double) :: val_c
    call pam_get_option_double_c( string_f2c(key,len_trim(key)) , val_c )
    val = val_c
  end subroutine\n\n''')

for tp in array_types :
  f.write(f'''  subroutine pam_create_array_{tp}(key,dims,desc_in)
    implicit none
    character(len=*), intent(in) :: key
    integer         , intent(in) :: dims(:)
    character(len=*), intent(in), optional :: desc_in
    integer(c_int) :: dims_c(size(dims))
    character(len=maxlen) :: desc = ""
    if (present(desc_in)) desc = trim(desc_in)
    dims_c = dims(ubound(dims,1):lbound(dims,1):-1)
    call pam_create_array_{tp}_c( string_f2c(key,len_trim(key)) , &
                                  string_f2c(desc,len_trim(desc)) , &
                                  dims_c , size(dims_c) )
  end subroutine\n''')

f.write("\n\n")

for tp in array_types :
  iso_c  = array_iso_c_map[tp]
  for d in range(1,8) :
    f.write(f'''  subroutine pam_get_array_{tp}_{d}d(label,ptr)
    implicit none
    character(len=*)                    , intent(in   ) :: label
    {iso_c}, pointer, contiguous                :: ptr(:{",:"*(d-1)})
    integer(c_int) :: dims({d})
    type(c_ptr)    :: ptr_c
    call pam_get_array_{tp}_c( string_f2c(label,len_trim(label)) , ptr_c , dims , {d} )
    dims = dims(ubound(dims,1):lbound(dims,1):-1)
    call c_f_pointer( ptr_c , ptr , [dims(1)''')
    for i in range(d-1) :
      f.write(f",dims({i+2})")
    f.write('''] )
  end subroutine\n''')

f.write("\n\n")

for tp in array_types :
  iso_c  = array_iso_c_map[tp]
  for d in range(1,8) :
    f.write(f'''  subroutine pam_mirror_array_readonly_{tp}_{d}d(label,arr,desc_in)
    implicit none
    character(len=*), intent(in) :: label
    {iso_c}, intent(in) :: arr(:{",:"*(d-1)})
    character(len=*), intent(in), optional :: desc_in
    integer(c_int) :: dims_c({d})
    character(len=maxlen) :: desc = ""
    if (present(desc_in)) desc = trim(desc_in)
    dims_c = shape(arr)
    dims_c = dims_c(ubound(dims_c,1):lbound(dims_c,1):-1)
    call pam_mirror_array_readonly_{tp}_c( string_f2c(label,len_trim(label)) , &
                                           string_f2c(desc,len_trim(desc)) , dims_c , {d} , arr )
  end subroutine\n''')

f.write("\n\n")

for tp in array_types :
  iso_c  = array_iso_c_map[tp]
  for d in range(1,8) :
    f.write(f'''  subroutine pam_mirror_array_readwrite_{tp}_{d}d(label,arr,desc_in)
    implicit none
    character(len=*), intent(in) :: label
    {iso_c}, intent(in) :: arr(:{",:"*(d-1)})
    character(len=*), intent(in), optional :: desc_in
    integer(c_int) :: dims_c({d})
    character(len=maxlen) :: desc = ""
    if (present(desc_in)) desc = trim(desc_in)
    dims_c = shape(arr)
    dims_c = dims_c(ubound(dims_c,1):lbound(dims_c,1):-1)
    call pam_mirror_array_readwrite_{tp}_c( string_f2c(label,len_trim(label)) , &
                                            string_f2c(desc,len_trim(desc)) , dims_c , {d} , arr )
  end subroutine\n''')

f.write("\n\n")

f.write('''\nend module pam_fortran_interface\n\n''')

f.close()
os.chmod("pam_fortran_interface.F90", stat.S_IRUSR | stat.S_IRGRP | stat.S_IROTH)


