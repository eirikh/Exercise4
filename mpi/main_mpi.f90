program main_mpi

implicit none
include 'mpif.h'

integer(kind=8) :: i,j,n
integer(kind=8), parameter :: sumsize = 20
integer(kind=8), parameter :: ndim = 2**sumsize 
real(kind=8), parameter :: sinf = (4*datan(1.0d0))**2/6.0D0
real(kind=8):: xerror(sumsize)
real(kind=8):: vector(ndim)
real(kind=8):: xsum(sumsize)
real(kind=4):: wstart, wstop
integer :: ierr

data xsum/ sumsize * 0.0D0/

call second(wstart)

call mpi_init(ierr)

do i = 1,ndim
   vector(i) = 1.0D0/(i*i)
enddo

xsum(1) = vector(1)

do i = 1,sumsize
   n=2**i
   do j=n/2+1,n
      xsum(i) = xsum(i) + vector(j)
   enddo
enddo

do i = 2,sumsize
   xsum(i) = xsum(i) + xsum(i-1)
   xerror(i) = sinf - xsum(i)
enddo

call second(wstop)
write(6,*) 'time',wstop-wstart

!Output
do i = 3,sumsize
   write(*,*) xsum(i), xerror(i), (xsum(i) + xerror(i))
enddo

call mpi_finalize(ierr)

end program