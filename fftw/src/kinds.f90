!> This module provides basic kinds of floats, complex, integers, chars
module kinds_m

  use iso_fortran_env
  use iso_c_binding, only: c_float, c_double, c_long_double
  use iso_c_binding, only: c_short, c_int, c_long, c_long_long, c_signed_char, c_char

  !> basic kinds of floats, complex, integers, chars
  !!
  !! 'dp' double precision type.
  !! 'sp' single precision type.
  !! 'gp' non-signal/geometry real data precision
  !! 'vp' precision for real signal data

  implicit none
  private
  public dp,sp,qp,gp,vp
  public i1,i2,i4,i8,b1

  integer, parameter :: qp= c_long_double       ! quad   precision
  integer, parameter :: dp= c_double            ! double precision
  integer, parameter :: sp= c_float             ! single precision

  !$if($SINGLE_PRECISION)
  !integer, parameter :: gp = sp                 ! non-signal/geometry precision
  !$else
  integer, parameter :: gp = dp                 ! default double
  !$endif
  integer, parameter :: vp= sp                  ! signal precision real[complex](kind=sp)

  integer, parameter :: i1= c_signed_char
  integer, parameter :: i2= c_short
  integer, parameter :: i4= c_int
  !$if($BIT32)
  !integer, parameter :: i8= c_long_long         ! 32-bit - short:2, int:4, long:4, long_long:8
  !$else
  integer, parameter :: i8= c_long
  !$endif

  integer, parameter :: b1= c_char
end module kinds_m
