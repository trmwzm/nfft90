module fftw_m
  use, intrinsic :: iso_c_binding
  use kinds_m, only : vp
  use utils_m, only : file_exists, file_open

  implicit none
  private
  include 'fftw3.f03'

  public :: alloc_fftw, dealloc_fftw, fftgood, get_wisdom_file, import_wisdom_file
  public :: freqndx, ndxfreq, export_wisdom_string

  character(len=128) wisdomfile
  character(len=128) wisdompath
  character(len=:), allocatable :: wisdomstring

  type, public :: plan_pair
     type(C_PTR)::  anl, syn  ! or forward/backward, or +1/-1 from sign of exp
  end type plan_pair

  integer, save :: plan_mode = FFTW_MEASURE

  type, public :: fftw
     type(c_ptr), dimension(2) :: p = [C_NULL_PTR,C_NULL_PTR]  ! plan for/backward, if any dirs is +1/-1
     integer :: rank = 0
     integer, dimension(:), allocatable :: dims
     integer, dimension(:), allocatable :: dirs ! output exponential factor signs, or 0 if no FT
     logical :: inplace = .true.
  contains
     procedure :: plan
     procedure :: exec
     procedure :: dealloc
  end type fftw

  interface alloc_fftw
     module procedure alloc_fftw_1d_c, alloc_fftw_2d_c, alloc_fftw_3d_c
     module procedure alloc_fftw_1d_r, alloc_fftw_2d_r, alloc_fftw_3d_r
  end interface alloc_fftw

  interface dealloc_fftw
     module procedure dealloc_fftw_c1, dealloc_fftw_r1,dealloc_fftw_c2
     module procedure dealloc_fftw_r2, dealloc_fftw_c3, dealloc_fftw_r3
  end interface dealloc_fftw

  interface
     integer(C_INT) function strlen (s) bind (C, name='strlen')
       import
       type(C_PTR), value :: s
     end function strlen
     subroutine free_pointer (p) bind (C, name='free')
       import
       type(C_PTR), value :: p
     end subroutine free_pointer
     function fopen (path,mode) bind (c,name='fopen')
       import
       character(c_char) path(*),mode(*)
       type(c_ptr) :: fopen
     end function fopen
     integer(c_int) function fClose (file) bind (c, NAME='fclose')
       import
       type(c_ptr), value :: file
     end function fClose
  end interface
contains
  subroutine alloc_fftw_1d_c (x, n1)
    complex(vp), intent(out), pointer :: x(:)
    integer, intent(in) :: n1
    type(c_ptr) :: p
    p= fftwf_alloc_complex(int(n1, c_size_t))
    call c_f_pointer(p, x, [n1])
  end subroutine alloc_fftw_1d_c
  subroutine alloc_fftw_2d_c (x, n1, n2)
    complex(vp), intent(out), pointer :: x(:, :)
    integer, intent(in) :: n1, n2
    type(c_ptr) :: p
    p= fftwf_alloc_complex(int(n1*n2, c_size_t))
    call c_f_pointer(p, x, [n1, n2])
  end subroutine alloc_fftw_2d_c
  subroutine alloc_fftw_3d_c (x, n1, n2, n3)
    complex(vp), intent(out), pointer :: x(:, :, :)
    integer, intent(in) :: n1, n2, n3
    type(c_ptr) :: p
    p= fftwf_alloc_complex(int(n1*n2*n3, c_size_t))
    call c_f_pointer(p, x, [n1, n2, n3])
  end subroutine alloc_fftw_3d_c

  subroutine dealloc_fftw_c1 (x)
    complex(vp), intent(in), target :: x(:)
    call fftwf_free(c_loc(x))
  end subroutine dealloc_fftw_c1
  subroutine dealloc_fftw_c2 (x)
    complex(vp), intent(in), target :: x(:,:)
    call fftwf_free(c_loc(x))
  end subroutine dealloc_fftw_c2
  subroutine dealloc_fftw_c3 (x)
    complex(vp), intent(in), target :: x(:,:,:)
    call fftwf_free(c_loc(x))
  end subroutine dealloc_fftw_c3

  subroutine alloc_fftw_1d_r (x, n1)
    real(vp), intent(out), pointer :: x(:)
    integer, intent(in) :: n1
    type(c_ptr) :: p
    p= fftwf_alloc_real(int(n1, c_size_t))
    call c_f_pointer(p, x, [n1])
  end subroutine alloc_fftw_1d_r
  subroutine alloc_fftw_2d_r (x, n1, n2)
    real(vp), intent(out), pointer :: x(:, :)
    integer, intent(in) :: n1, n2
    type(c_ptr) :: p
    p= fftwf_alloc_real(int(n1*n2, c_size_t))
    call c_f_pointer(p, x, [n1, n2])
  end subroutine alloc_fftw_2d_r
  subroutine alloc_fftw_3d_r (x, n1, n2, n3)
    real(vp), intent(out), pointer :: x(:, :, :)
    integer, intent(in) :: n1, n2, n3
    type(c_ptr) :: p
    p= fftwf_alloc_real(int(n1*n2*n3, c_size_t))
    call c_f_pointer(p, x, [n1, n2, n3])
  end subroutine alloc_fftw_3d_r

  subroutine dealloc_fftw_r1 (x)
    real(vp), intent(in), target :: x(:)
    call fftwf_free(c_loc(x))
  end subroutine dealloc_fftw_r1
  subroutine dealloc_fftw_r2 (x)
    real(vp), intent(in), target :: x(:,:)
    call fftwf_free(c_loc(x))
  end subroutine dealloc_fftw_r2
  subroutine dealloc_fftw_r3 (x)
    real(vp), intent(in), target :: x(:,:,:)
    call fftwf_free(c_loc(x))
  end subroutine dealloc_fftw_r3

  subroutine plan (this, ashp, a, b, dirs)
    class(fftw), intent(inout) :: this
    integer, dimension(:), intent(in) :: ashp
    complex(vp), dimension(*), intent(out) :: a
    complex(vp), dimension(*), intent(out), optional :: b
    integer, dimension(size(ashp)), intent(in), optional :: dirs

    integer, dimension(:), allocatable :: stds
    integer, dimension(:), allocatable :: dimlp,stdlp
    integer, dimension(:), allocatable :: dimft,stdft
    logical, dimension(:), allocatable :: mskdir
    integer :: idir,ift,i,nft,nlp
    integer, dimension(size(ashp)) :: dirs0

   if (.not.present(dirs)) then
       dirs0= 1
    else
       dirs0= dirs
    endif

    this%dims= ashp
    this%rank= size(this%dims)
    this%dirs= dirs0

    if (this%rank/=size(this%dirs)) &
       stop "incorrect dirs length"

    allocate(stds(this%rank))
    stds= 1
    do i= 2,this%rank
       stds(i)= stds(i-1)*this%dims(i-1)
    end do

    do idir= 1,-1,-2
       ift= min(2,2-idir)
       mskdir= this%dirs==idir ! where FT
       nft= count(mskdir)
       this%p(ift)= c_null_ptr
       if (nft>0) then
          nlp= this%rank-nft
          dimft= pack(this%dims,mskdir)
          stdft= pack(stds,mskdir)
          dimlp= pack(this%dims,.not.mskdir)
          stdlp= pack(stds,.not.mskdir)
          if (size(dimlp)==0) then
             dimlp= [0]
             stdlp= [0]
          endif
          if (present(b)) then
             this%inplace= .false.
             call plan_wrkr(this,a,nft,dimft,stdft,nlp,dimlp,stdlp,idir,b)
          else
             this%inplace= .true.
             call plan_wrkr(this,a,nft,dimft,stdft,nlp,dimlp,stdlp,idir)
          endif
       endif
       if (allocated(dimft)) &
            deallocate(dimft,stdft,dimlp,stdlp)
    enddo
    ! call update wisdom file
    deallocate(stds,mskdir)
  end subroutine plan
  subroutine plan_wrkr (this, a, nft, dimft, stdft, nlp, dimlp, stdlp, idir, b)
    class(fftw), intent(inout) :: this
    complex(vp), dimension(*), intent(out) :: a
    integer, intent(in) :: nft, nlp, idir
    integer, dimension(:), allocatable, intent(in) :: dimft,stdft
    integer, dimension(:), allocatable, intent(in) :: dimlp,stdlp
    complex(vp), dimension(*), intent(out), optional :: b

    integer :: ift, isgn, i
    type(fftwf_iodim) :: idft(nft), idlp(nlp)

    ift= min(2,2-idir)
    isgn= -idir

    do i= 1,nft
       idft(i)%n= dimft(i)
       idft(i)%is= stdft(i)
       idft(i)%os= stdft(i)
    enddo
    do i= 1,nlp
       idlp(i)%n= dimlp(i)
       idlp(i)%is= stdlp(i)
       idlp(i)%os= stdlp(i)
    enddo

    if (.not.present(b)) then
       this%p(ift)= fftwf_plan_guru_dft(nft,idft,nlp,idlp,a,a,isgn,plan_mode)
    else
       this%p(ift)= fftwf_plan_guru_dft(nft,idft,nlp,idlp,a,b,isgn,plan_mode)
    endif
  end subroutine plan_wrkr

  ! exec: need pointer && allocatable versions
  subroutine exec(this, a, b)
    class(fftw), intent(in) :: this
    complex(vp), dimension(*), intent(inout) :: a
    complex(vp), dimension(*), intent(out), optional :: b

    integer :: i

    if (present(b)) then
       if (this%inplace) &
            stop "planned in place"
       do i= 1,2
          if (c_associated(this%p(i))) &
               call fftwf_execute_dft(this%p(i),a,b)
        enddo
    else
       if (.not.this%inplace) &
            stop "planned out of place"
       do i= 1,2
          if (c_associated(this%p(i))) &
               call fftwf_execute_dft(this%p(i),a,a)
       enddo
    end if
  end subroutine exec

  !clean up
  subroutine dealloc(this)
    class(fftw), intent(inout) :: this
    integer :: i
    do i= 1,2
       if (c_associated(this%p(i))) &
            call fftwf_destroy_plan(this%p(i))
    enddo
    if (allocated(this%dims)) then
       deallocate(this%dims)
       deallocate(this%dirs)
    endif
    this%p(1)= c_null_ptr
    this%p(2)= c_null_ptr
  end subroutine dealloc

  subroutine import_wisdom_string ()
    integer :: i
     ! Read FFTW Plans
    i= fftwf_import_wisdom_from_string(trim(wisdomstring)//C_NULL_CHAR )
  end subroutine import_wisdom_string

  subroutine export_wisdom_string
    character(C_CHAR), dimension(:), pointer :: s
    integer(C_SIZE_T) :: n, slen
    type(C_PTR) :: p
    p= fftwf_export_wisdom_to_string()
    if (.not.c_associated(p)) &
         stop 'error exporting wisdom'
    slen= strlen(p)
    call c_f_pointer(p,s,[slen+1])
    if (allocated(wisdomstring)) &
         deallocate(wisdomstring)
    allocate(character(slen) :: wisdomstring)
    do n= 1,slen
      wisdomstring(n:n)= s(n)
    end do
    print*,wisdomstring
    call free_pointer(p)
  end subroutine export_wisdom_string

  pure function f_c_string (f_string) result (c_string)
    character(len=*), intent(in) :: f_string
    character(len=1,kind=c_char) :: c_string(len_trim(f_string)+1)
    integer                      :: n, i

    n= len_trim(f_string)
    do i= 1,n
       c_string(i)= f_string(i:i)
    end do
    c_string(n+1)= c_null_char
  end function f_c_string

  subroutine import_wisdom_file ()
    integer(c_int) :: is
    type(c_ptr) :: stream
    stream= fopen(f_c_string(wisdomfile),f_c_string('r'))
    if (.not.c_associated(stream)) &
         stop "error opening wisdom file"
    is= fftwf_import_wisdom_from_file(stream)
    if (is<0) &
         stop "error reading wisdom file"
    is= fclose(stream)
  end subroutine import_wisdom_file

  subroutine get_wisdom_file ()
    wisdomfile= ''
    call getenv("FFTW_WISDOM",wisdompath)
    if (len_trim(wisdompath)>0) then
       wisdomfile= trim(wisdompath) // '/wisdom_f'
       if (file_exists(wisdomfile)) then
          if (.not.file_open(wisdomfile)) then
             return
          endif
       endif
       wisdomfile= ''
    endif
  end subroutine get_wisdom_file

  integer function fftgood(in, factors)
    integer, intent(in) :: in
    integer, intent(in), optional :: factors(:)

    integer :: nfactors
    integer :: nondiv
    integer, allocatable :: factors0(:), exponents(:)
    if (.not.present(factors)) then
       factors0= [2,3,5]
    else
       factors0= factors
    endif
    nfactors= size(factors0,dim=1)
    allocate(exponents(1:nfactors))
    fftgood= in
    do
      call get_exponents(fftgood,nfactors,factors0,exponents,nondiv)
      if(nondiv==1) &
           exit
      fftgood= fftgood+1
    end do
    deallocate(exponents,factors0)
  end function fftgood

  subroutine get_exponents(num, nfactors, factors, exponents, nondiv)
    integer, intent(in)  :: num
    integer, intent(in)  :: nfactors
    integer, intent(in)  :: factors(nfactors)
    integer, intent(out) :: exponents(nfactors)
    integer, intent(out) :: nondiv

    integer :: ifactor
    nondiv= num
    do ifactor= 1,nfactors
      exponents(ifactor)= 0
      do
        if(mod(nondiv,factors(ifactor))/=0) &
             exit
        nondiv= nondiv/factors(ifactor)
        exponents(ifactor)= exponents(ifactor)+1
      end do
    end do
  end subroutine get_exponents

  elemental function freqndx(ii, nn)
    integer, intent(in) :: ii,nn
    integer :: freqndx

    ! index to frequency
    if(ii<=nn/2 + 1) then
       freqndx= ii-1
    else
       freqndx= ii-nn-1
    end if
  end function freqndx

  elemental function ndxfreq(ii, nn)
    integer, intent(in) :: ii,nn
    integer :: ndxfreq

    ! frequency to index
    if(ii>=0) then
       ndxfreq= ii+1
    else
       ndxfreq= ii+nn+1
    end if
  end function ndxfreq
end module
