!> This module provides utilities for random number XXX
module utils_m
  use kinds_m
  use iso_fortran_env, only : error_unit
  use iso_c_binding
  implicit none

  private
  public :: is_one, statar
  public :: span, catch_error, assert, init_random_seed
  public :: random_n_number, random_n
  public :: strupper, strlower
  public :: fstr_cstr, cstr_fstr
  public :: file_exists, file_open, freeunit

  real(dp), parameter :: sqrt12= sqrt(12._dp)

  interface span
     module procedure span_d, span_r
  end interface span

  interface statar
     module procedure statar_s, statar_d
  end interface statar

  interface random_n_number
     module procedure random_n_number_1d,random_n_number_2d,random_n_number_3d
     module procedure random_n_number_1s,random_n_number_2s,random_n_number_3s
     module procedure random_n_number_1c,random_n_number_2c,random_n_number_3c
  end interface random_n_number

  interface file_exists
     module procedure file_unit_exists
     module procedure file_name_exists
  end interface file_exists

  interface file_open
     module procedure file_open_by_unit
     module procedure file_open_by_name
  end interface file_open
contains
  !> Aborts the program with nonzero exit code
  !> The statement "stop msg" will return 0 exit code when compiled using
  !> gfortran. catch_error() uses the statement "stop 1" which returns an exit code
  !> 1 and a print statement to print the message.
  !> call catch_error("Invalid argument")
  subroutine catch_error (msg, istop)
    character(len=*) :: msg ! Message to print on stdout
    integer, intent(in), optional :: istop
    write(error_unit,'(a)') 'Error: '//trim(msg)
    if (present(istop)) then
       if (abs(istop)>0) &
            call exit(1)
    endif
  end subroutine catch_error

  !> If condition == .false., it aborts the program.
  !> call assert(a == 5)
  subroutine assert (condition, msg, istop)
    logical, intent(in) :: condition
    character(len=*), optional :: msg ! Message to print on stdout
    integer, intent(in), optional :: istop
    if (.not.condition) then
       if(present(msg)) then
          call catch_error("Assert failed: "//trim(msg),istop)
       else
          call catch_error("Assert failed.",istop)
       endif
    end if
  end subroutine assert

  pure function is_one (opt) result (lopt)
!!$   return .TRUE. if integer opt present and == 1, else .FALSE.
    integer, optional,intent(in) :: opt
    logical :: lopt
    lopt= .false.
    if (present(opt)) then
       if (opt==1) then
          lopt= .true.
       end if
    endif
    return
  end function is_one

  function span_r (r1, r0, n) result (s)
    real(sp), dimension(n) :: s
    real(sp), intent(in) :: r1
    real(sp), intent(in) :: r0
    integer :: n

    real(sp) :: dr
    integer :: i
    dr= (r0-r1)/max(1,n-1)
    s= r1+[(i-1,i=1,n)]*dr
  end function span_r

  function span_d (r1, r0, n) result (s)
    real(dp), dimension(n) :: s
    real(dp), intent(in) :: r1
    real(dp), intent(in) :: r0
    integer :: n

    real(dp) :: dr
    integer :: i
    dr= (r0-r1)/max(1,n-1)
    s= r1+[(i-1,i=1,n)]*dr
  end function span_d

  subroutine statar_s (a,verbose,n,mn,mx,avg,rms)
    real(sp), intent(in)            :: a(:)
    integer,  optional, intent(in)  :: verbose
    integer,  intent(out)           :: n
    real(sp), optional, intent(out) :: mn,mx,avg,rms
    n= size(a)
    mn= minval(a)
    mx= maxval(a)
    avg= sum(a)/n
    rms= sum((a-avg)**2)/n
    if (is_one(verbose)) then
       print *,'N#  ',n
       print *,'min ',mn
       print *,'Max ',mx
       print *,'avg ',avg
       print *,'rms ',rms
    endif
  end subroutine statar_s
  subroutine statar_d (a,verbose,n,mn,mx,avg,rms)
    real(dp), intent(in)            :: a(:)
    integer,  optional, intent(in)  :: verbose
    integer,  intent(out)           :: n
    real(dp), optional, intent(out) :: mn,mx,avg,rms
    n= size(a)
    mn= minval(a)
    mx= maxval(a)
    avg= sum(a)/n
    rms= sum((a-avg)**2)/n
    if (is_one(verbose)) then
       print *,'N#  ',n
       print *,'min ',mn
       print *,'Max ',mx
       print *,'avg ',avg
       print *,'rms ',rms
    endif
  end subroutine statar_d

  subroutine init_random_seed()
    integer :: i, n, clock
    integer, allocatable :: seed(:)
    call random_seed(size=n)
    allocate(seed(n))
    call system_clock(count=clock)
    seed = clock + 37 * [ (i - 1, i = 1, n) ]
    call random_seed(put=seed)
  end subroutine init_random_seed

  ! generate a random normal deviate using the polar method.
  ! reference: marsaglia,g. & bray,t.a. 'a convenient method for generating
  ! normal variables', siam rev., vol.6, 260-264, 1964.
  ! this subroutine was taken from: http://jblevins.org/mirror/amiller/rnorm.f90
  function random_n () result (r)
    real(dp) :: r
    real(dp) :: rr(12)
    call random_number(rr)
    r= sum(rr)-6
  end function random_n
  subroutine random_n_number_1d (x)
    real(dp), intent(out), dimension(:) :: x
    real(dp) :: xx(12,size(x))
    call random_number(xx)
    x= sum(xx,1)-6
  end subroutine random_n_number_1d
  subroutine random_n_number_2d (x)
    real(dp), intent(out), dimension(:,:) :: x
    real(dp) :: xx(12,size(x,1),size(x,2))
    call random_number(xx)
    x= sum(xx,1)-6
  end subroutine random_n_number_2d
  subroutine random_n_number_3d (x)
    real(dp), intent(out), dimension(:,:,:) :: x
    real(dp) :: xx(12,size(x,1),size(x,2),size(x,3))
    call random_number(xx)
    x= sum(xx,1)-6
  end subroutine random_n_number_3d
  subroutine random_n_number_1s (x)
    real(sp), intent(out), dimension(:) :: x
    real(sp) :: xx(12,size(x))
    call random_number(xx)
    x= sum(xx,1)-6
  end subroutine random_n_number_1s
  subroutine random_n_number_2s (x)
    real(sp), intent(out), dimension(:,:) :: x
    real(sp) :: xx(12,size(x,1),size(x,2))
    call random_number(xx)
    x= sum(xx,1)-6
  end subroutine random_n_number_2s
  subroutine random_n_number_3s (x)
    real(sp), intent(out), dimension(:,:,:) :: x
    real(sp) :: xx(12,size(x,1),size(x,2),size(x,3))
    call random_number(xx)
    x= sum(xx,1)-6
  end subroutine random_n_number_3s
  subroutine random_n_number_1c (x)
    complex(vp), intent(out), dimension(:) :: x
    real(vp) :: xx(12,2,size(x))
    call random_number(xx)
    x= cmplx(sum(xx(:,1,:),1)-6,sum(xx(:,2,:),1)-6)
  end subroutine random_n_number_1c
  subroutine random_n_number_2c (x)
    complex(vp), intent(out), dimension(:,:) :: x
    real(vp) :: xx(12,2,size(x,1),size(x,2))
    call random_number(xx)
    x= cmplx(sum(xx(:,1,:,:),1)-6,sum(xx(:,2,:,:),1)-6)
  end subroutine random_n_number_2c
  subroutine random_n_number_3c (x)
    complex(vp), intent(out), dimension(:,:,:) :: x
    real(vp) :: xx(12,2,size(x,1),size(x,2),size(x,3))
    call random_number(xx)
    x= cmplx(sum(xx(:,1,:,:,:),1)-6,sum(xx(:,2,:,:,:),1)-6)
  end subroutine random_n_number_3c

  function strlower (string) result (tolower_result)
    character (len=*), intent(in) :: string
    character (len=len(string)) :: tolower_result

    integer         :: i,ii

    do i= 1,len(string)
       ii= iachar(string(i:i))
       select case (ii)
       case (65:90)
          tolower_result(i:i)= achar(ii+32)
       case default
          tolower_result(i:i)= string(i:i)
       end select
    enddo
  end function strlower

  function strupper (string) result (toupper_result)
    character (len=*), intent(in) :: string
    character (len=len(string)) :: toupper_result

    integer         :: i,ii

    do i= 1,len(string)
       ii= iachar(string(i:i))
       select case (ii)
       case (97:122)
          toupper_result(i:i)= achar(ii-32)
       case default
          toupper_result(i:i)= string(i:i)
       end select
    enddo
  end function strupper

  !  type(c_ptr), target, intent(in) :: cstring
  !  character(c_char), pointer :: c_str(:)
  !  character(len=:), allocatable :: f_str
  !  call c_f_pointer(cstring, c_str, [slen])
  !  call fstr_cstr(f_str,c_str)
  subroutine fstr_cstr (f_str,c_str)
    use, intrinsic :: iso_c_binding, only: c_char, c_null_char
    character(kind=c_char, len=1), intent(in) :: c_str(*)
    character(len=:), allocatable, intent(out) :: f_str

    integer :: i, nc
    i= 1
    do while (c_str(i)/=c_null_char)
       i= i+1
    enddo
    nc= i-1
    if (.not. allocated(f_str)) then
       allocate(character(len=nc)::f_str)
    else
       if (len(f_str)<nc) &
            print*,'string too short'
    endif
    f_str = transfer(c_str(1:nc),f_str)
  end subroutine fstr_cstr

  ! did not make c_str allocatable as it must be fixed length in C API
  subroutine cstr_fstr (c_str, f_str)
    use, intrinsic :: iso_c_binding, only: c_char, c_null_char
    character(len=*), intent(in) :: f_str
    character(kind=c_char, len=1), intent(out) :: c_str(*)

    integer :: n,i
    n= len(f_str)
    do i= 1,n
       c_str(i)= f_str(i:i)
    enddo
    c_str(i)= c_null_char
  end subroutine cstr_fstr

  function freeunit () result (f)
    integer :: f
    integer :: umx= 200
    logical :: openstatus

    f= 10
    inquire(f,opened=openstatus)
    do while (openstatus .and. f<umx)
       f= f+1
       inquire(unit=f,opened=openstatus)
    end do
    if (f==umx) &
         f= -1
  end function freeunit

  function file_unit_exists (fh) result (ex)
    integer, intent(in) :: fh
    logical :: ex

    inquire(unit=fh,exist=ex)
  end function file_unit_exists

  function file_name_exists (filename) result (ex)
    character(len=*), intent(in) :: filename
    logical :: ex
    integer :: trimlen

    trimlen= len_trim(filename)
    inquire(file=filename(1:trimlen),exist= ex)
  end function file_name_exists

  function file_open_by_unit (fh) result (is_open)
    integer, intent(in) :: fh
    logical :: is_open

    inquire(unit=fh,opened=is_open)
  end function file_open_by_unit

  function file_open_by_name (filename) result (is_open)
    character(len=*), intent(in) :: filename
    logical :: is_open
    integer :: trimlen

    trimlen= len_trim(filename)
    inquire(file=filename(1:trimlen),opened=is_open)
  end function file_open_by_name
end module utils_m
