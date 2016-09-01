program test_fftw2
  use kinds_m, only: vp, gp
  use fftw_m, only: fftw, alloc_fftw, dealloc_fftw, get_wisdom_file,import_wisdom_file
  use fftw_m, only: export_wisdom_string
  use utils_m, only : statar
  implicit none

  type(fftw) :: ft,fti,ft2
  complex(vp), pointer :: z(:),w(:)
  complex(vp), pointer :: z2(:,:),z22(:,:)
  complex(vp), allocatable :: za(:),wa(:)
  complex(vp), allocatable :: z2a(:,:),z22a(:,:)
  real(vp), allocatable :: xr(:), xi(:)
  real(vp), allocatable :: xr2(:,:), xi2(:,:)
  integer :: n= 256
  integer :: nn
  real(vp) :: mn,mx,avg,rms
  real(vp) :: t1, t2, t3

  call get_wisdom_file()
  call import_wisdom_file()
  !call export_wisdom_string()

  call alloc_fftw(z,n)
  allocate(xr(n),xi(n))

  call random_number(xr)
  call random_number(xi)
  xr= xr-0.5
  xi= xi-0.5

  call ft%plan(shape(z),z,dirs=[1])
  call fti%plan(shape(z),z,dirs=[-1])

  z= cmplx(xr,xi)
  call ft%exec(z)
  call fti%exec(z)
  xr= abs(z/n-cmplx(xr,xi))
  print *, '0: fftw 1-D forw/backward'
  call statar(reshape(xr,[n]),1,nn,mn,mx,avg,rms)
  print *,''

  call ft%dealloc()
  call fti%dealloc()    ! ------- done ----------

  call alloc_fftw(z2,n,n)
  allocate(xr2(n,n),xi2(n,n))

  call random_number(xr2)
  call random_number(xi2)
  xr2= xr2-0.5
  xi2= xi2-0.5

  call ft%plan(shape(z2),z2,dirs=[1,1])
  call fti%plan(shape(z2),z2,dirs=[-1,-1])

  z2= cmplx(xr2,xi2)
  call ft%exec(z2)
  call fti%exec(z2)
  xr2= abs(z2/(n*n)-cmplx(xr2,xi2))
  print *, '1: fftw 2-D forw/backward'
  call statar(reshape(xr2,[n*n]),1,nn,mn,mx,avg,rms)
  print *,''

  call ft%dealloc()
  call fti%dealloc()    ! ------- done ----------

  call random_number(xr2)
  call random_number(xi2)
  xr2= xr2-0.5
  xi2= xi2-0.5

  call ft%plan(shape(z2),z2,dirs=[1,0])
  call fti%plan(shape(z2),z2,dirs=[-1,0])

  z2= cmplx(xr2,xi2)
  call ft%exec(z2)
  call fti%exec(z2)

  xr2= abs(z2/n-cmplx(xr2,xi2))
  print *, '2: fftw 2-D 1st dim forw/backward'
  call statar(reshape(xr2,[size(xr2)]),1,nn,mn,mx,avg,rms)
  print *,''

  call ft%dealloc()
  call fti%dealloc()    ! ------- done ----------

  call random_number(xr2)
  call random_number(xi2)
  xr2= xr2-0.5
  xi2= xi2-0.5

  call ft%plan(shape(z2),z2,dirs=[0,1])
  call fti%plan(shape(z2),z2,dirs=[0,-1])

  z2= cmplx(xr2,xi2)
  call ft%exec(z2)
  call fti%exec(z2)
  xr2= abs(z2/n-cmplx(xr2,xi2))
  print *, '3: fftw 2-D 2nd dim forw/backward'
  call statar(reshape(xr2,[size(xr2)]),1,nn,mn,mx,avg,rms)
  print *,''

  call ft%dealloc()
  call fti%dealloc()    ! ------- done ----------

  call random_number(xr2)
  call random_number(xi2)
  xr2= xr2-0.5
  xi2= xi2-0.5

  call ft%plan(shape(z2),z2,dirs=[1,1])
  call ft2%plan(shape(z2),z2,dirs=[-1,0])
  call fti%plan(shape(z2),z2,dirs=[0,-1])

  z2= cmplx(xr2,xi2)
  call ft%exec(z2)
  call ft2%exec(z2)
  call fti%exec(z2)
  xr2= abs(z2/(n*n)-cmplx(xr2,xi2))
  print *, '4: fftw 2-D forw/ 1st then 2nd backward'
  call statar(reshape(xr2,[size(xr2)]),1,nn,mn,mx,avg,rms)
  print *,''

  call ft%dealloc()
  call ft2%dealloc()
  call fti%dealloc()    ! ------- done ----------

  call alloc_fftw(z22,n,n)
  z22= 0

  call random_number(xr2)
  call random_number(xi2)
  xr2= xr2-0.5
  xi2= xi2-0.5

  call ft%plan(shape(z2),z2,z22,dirs=[1,1])
  call fti%plan(shape(z2),z2,z22,dirs=[-1,-1])

  z2= cmplx(xr2,xi2)
  call ft%exec(z2,z22)
  call fti%exec(z22,z2)
  xr2= abs(z2/(n*n)-cmplx(xr2,xi2))
  print *, '5: fftw 2-D forw/backward out-of-place'
  call statar(reshape(xr2,[size(xr2)]),1,nn,mn,mx,avg,rms)
  print *,''

  call ft%dealloc()
  call fti%dealloc()    ! ------- done ----------

  allocate(za(n))

  call random_number(xr)
  call random_number(xi)
  xr= xr-0.5
  xi= xi-0.5

  call ft%plan(shape(za),za,dirs=[1])
  call fti%plan(shape(za),za,dirs=[-1])

  za= cmplx(xr,xi)
  call ft%exec(za)
  call fti%exec(za)
  xr= abs(za/n-cmplx(xr,xi))
  print *, '6: fftw 1-D forw/backward, -- allocated'
  call statar(reshape(xr,[n]),1,nn,mn,mx,avg,rms)
  print *,''

  deallocate(za,xr,xi)
  call ft%dealloc()
  call fti%dealloc()    ! ------- done ----------

  allocate(z2a(n,n))

  call random_number(xr2)
  call random_number(xi2)
  xr2= xr2-0.5
  xi2= xi2-0.5

  call ft%plan(shape(z2a),z2a,dirs=[1,1])
  call fti%plan(shape(z2a),z2a,dirs=[-1,-1])

  z2a= cmplx(xr2,xi2)
  call ft%exec(z2a)
  call fti%exec(z2a)
  xr2= abs(z2a/(n*n)-cmplx(xr2,xi2))
  print *, '7: fftw 2-D forw/backward, -- allocated'
  call statar(reshape(xr2,[n*n]),1,nn,mn,mx,avg,rms)
  print *,''

  call ft%dealloc()
  call fti%dealloc()    ! ------- done ----------

  call random_number(xr2)
  call random_number(xi2)
  xr2= xr2-0.5
  xi2= xi2-0.5

  call ft%plan(shape(z2a),z2a,dirs=[1,0])
  call fti%plan(shape(z2a),z2a,dirs=[-1,0])

  z2a= cmplx(xr2,xi2)
  call ft%exec(z2a)
  call fti%exec(z2a)

  xr2= abs(z2a/n-cmplx(xr2,xi2))
  print *, '8: fftw 2-D 1st dim forw/backward, -- allocated'
  call statar(reshape(xr2,[size(xr2)]),1,nn,mn,mx,avg,rms)
  print *,''

  call ft%dealloc()
  call fti%dealloc()    ! ------- done ----------

  call random_number(xr2)
  call random_number(xi2)
  xr2= xr2-0.5
  xi2= xi2-0.5

  call ft%plan(shape(z2a),z2a,dirs=[0,1])
  call fti%plan(shape(z2a),z2a,dirs=[0,-1])

  z2a= cmplx(xr2,xi2)
  call ft%exec(z2a)
  call fti%exec(z2a)
  xr2= abs(z2a/n-cmplx(xr2,xi2))
  print *, '9: fftw 2-D 2nd dim forw/backward, -- allocated'
  call statar(reshape(xr2,[size(xr2)]),1,nn,mn,mx,avg,rms)
  print *,''

  call ft%dealloc()
  call fti%dealloc()    ! ------- done ----------

  call random_number(xr2)
  call random_number(xi2)
  xr2= xr2-0.5
  xi2= xi2-0.5

  call ft%plan(shape(z2a),z2a,dirs=[1,1])
  call ft2%plan(shape(z2a),z2a,dirs=[-1,0])
  call fti%plan(shape(z2a),z2a,dirs=[0,-1])

  z2a= cmplx(xr2,xi2)
  call ft%exec(z2a)
  call ft2%exec(z2a)
  call fti%exec(z2a)
  xr2= abs(z2a/(n*n)-cmplx(xr2,xi2))
  print *, '10: fftw 2-D forw/ 1st then 2nd backward, -- allocated'
  call statar(reshape(xr2,[size(xr2)]),1,nn,mn,mx,avg,rms)
  print *,''

  call ft%dealloc()
  call ft2%dealloc()
  call fti%dealloc()    ! ------- done ----------

  allocate(z22a(n,n))
  z22a= 0

  call random_number(xr2)
  call random_number(xi2)
  xr2= xr2-0.5
  xi2= xi2-0.5

  call ft%plan(shape(z2a),z2a,z22a,dirs=[1,1])
  call fti%plan(shape(z2a),z2a,z22a,dirs=[-1,-1])

  z2a= cmplx(xr2,xi2)

  call ft%exec(z2a,z22a)
  call fti%exec(z22a,z2a)

  xr2= abs(z2a/(n*n)-cmplx(xr2,xi2))
  print *, '11: fftw 2-D forw/backward out-of-place, -- allocated'
  call statar(reshape(xr2,[size(xr2)]),1,nn,mn,mx,avg,rms)
  print *,''

  call ft%dealloc()
  call fti%dealloc()    ! ------- done ----------

  call cpu_time(t1)
  call ft%plan(shape(z2),z2,dirs=[1,1])
  call fti%plan(shape(z2),z2,dirs=[-1,-1])
  call cpu_time(t2)

  z2= cmplx(xr2,xi2)
  call ft%exec(z2)
  call fti%exec(z2)
  call cpu_time(t3)

  xr2= abs(z2/(n*n)-cmplx(xr2,xi2))
  print *,''
  print *, '12: fftw 2-D forw/backward'
  print *,'    alloc_fftw:'
  call statar(reshape(xr2,[n*n]),1,nn,mn,mx,avg,rms)
  print *, "Total time:", (t3-t1)*1000, "ms"
  print *, "plan:      ", (t2-t1)*1000, "ms"
  print *, "exec:      ", (t3-t2)*1000, "ms"

  call ft%dealloc()
  call fti%dealloc()    ! ------- done ----------

  call cpu_time(t1)
  call ft%plan(shape(z2a),z2a,dirs=[1,1])
  call fti%plan(shape(z2a),z2a,dirs=[-1,-1])
  call cpu_time(t2)

  z2a= cmplx(xr2,xi2)
  call ft%exec(z2a)
  call fti%exec(z2a)
  call cpu_time(t3)

  xr2= abs(z2a/(n*n)-cmplx(xr2,xi2))
  print *,''
  print *, '13: fftw 2-D forw/backward'
  print *,'    allocate:'
  call statar(reshape(xr2,[n*n]),1,nn,mn,mx,avg,rms)
  print *, "Total time:", (t3-t1)*1000, "ms"
  print *, "plan:      ", (t2-t1)*1000, "ms"
  print *, "exec:      ", (t3-t2)*1000, "ms"

  call ft%dealloc()
  call fti%dealloc()    ! ------- done ----------
  call dealloc_fftw(z)

  call alloc_fftw(z,n)
  call alloc_fftw(w,n)
  allocate(za(n),wa(n))
  allocate(xr(n),xi(n))

  call random_number(xr)
  call random_number(xi)
  xr= xr-0.5
  xi= xi-0.5

  ! in-place
  call ft%plan(shape(z),z,dirs=[1])
  call fti%plan(shape(z),z,dirs=[-1])

  z= cmplx(xr,xi)
  call ft%exec(z)
  call fti%exec(z)
  xr= abs(z/n-cmplx(xr,xi))
  print *, '14: fftw 1-D forw/backward, alloc_fftw, inplace'
  call statar(reshape(xr,[n]),1,nn,mn,mx,avg,rms)
  print *,''

  call ft%dealloc()
  call fti%dealloc()    ! ------- done ----------

  call ft%plan(shape(za),za,dirs=[1])
  call fti%plan(shape(za),za,dirs=[-1])
  print*,'plan done'

  za= cmplx(xr,xi)
  call ft%exec(za)
  call fti%exec(za)
  xr= abs(za/n-cmplx(xr,xi))
  print *, '15: fftw 1-D forw/backward, allocatable, inplace'
  call statar(reshape(xr,[n]),1,nn,mn,mx,avg,rms)
  print *,''

  call ft%dealloc()
  call fti%dealloc()    ! ------- done ----------

  ! out-place
  call ft%plan(shape(z),z,w,dirs=[1])
  call fti%plan(shape(z),w,z,dirs=[-1])

  z= cmplx(xr,xi)
  call ft%exec(z,w)
  call fti%exec(w,z)
  xr= abs(z/n-cmplx(xr,xi))
  print *, '16: fftw 1-D forw/backward, alloc_fftw, outplace'
  call statar(reshape(xr,[n]),1,nn,mn,mx,avg,rms)
  print *,''

  call ft%dealloc()
  call fti%dealloc()    ! ------- done ----------

  call ft%plan(shape(za),za,wa,dirs=[1])
  call fti%plan(shape(za),wa,za,dirs=[-1])

  za= cmplx(xr,xi)
  call ft%exec(za,wa)
  call fti%exec(wa,za)
  xr= abs(za/n-cmplx(xr,xi))
  print *, '17: fftw 1-D forw/backward, allocatable, outplace'
  call statar(reshape(xr,[n]),1,nn,mn,mx,avg,rms)
  print *,''

  call ft%dealloc()
  call fti%dealloc()    ! ------- done ----------

  ! clean up
  call dealloc_fftw(w)
  call dealloc_fftw(z)
  call dealloc_fftw(z2)
  call dealloc_fftw(z22)

  deallocate(xr,xi,za,wa)
  deallocate(z22a,z2a,xr2,xi2)
contains
  
end program test_fftw2
