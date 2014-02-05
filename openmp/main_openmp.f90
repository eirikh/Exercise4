program main_openmp

integer(kind=8) :: i,j,n
integer(kind=8), parameter :: sumsize = 20
integer(kind=8), parameter :: ndim = 2**sumsize 
real(kind=8), parameter :: sinf = (4*datan(1.0d0))**2/6.0D0
real(kind=8):: xerror(sumsize)
real(kind=8):: vector(ndim)
real(kind=8):: xsum(sumsize)
real(kind=4)::wstart,wstop
real(kind=8)::omp_get_wtime

data xsum/ sumsize * 0.0D0/

call second(wstart)

!$omp parallel do schedule(dynamic)
do i = 1,ndim
   vector(i) = 1.0D0/(i*i)
enddo
!$omp end parallel do

xsum(1) = vector(1)

!$omp parallel do schedule(dynamic) private(j,n)
do i = 1,sumsize
   n=2**i
   do j=n/2+1,n
      xsum(i) = xsum(i) + vector(j)
   enddo
enddo
!$omp end parallel do

do i = 2,sumsize
   xsum(i) = xsum(i) + xsum(i-1)
   xerror(i) = sinf - xsum(i)
enddo

call second(wstop)
write(6,*) 'time',wstop-wstart

!Output
do i = 3,sumsize
   write(6,*) xsum(i), xerror(i), (xsum(i) + xerror(i))
enddo

end program
