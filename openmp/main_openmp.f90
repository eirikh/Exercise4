program main_openmp

integer(kind=8) :: i,j,n
integer(kind=8), parameter :: sumsize = 20
integer(kind=8), parameter :: ndim = 2**sumsize 
real(kind=8), parameter :: sinf = (4*datan(1.0d0))**2/6.0D0
real(kind=8):: xerror(sumsize)
real(kind=8):: vector(ndim)
real(kind=8):: xsum(sumsize)
real(kind=8)::wstart
real(kind=8)::omp_get_wtime
integer :: omp_get_num_threads
integer :: omp_get_thread_num

data xsum/ sumsize * 0.0D0/

wstart = omp_get_wtime()

!$omp parallel do schedule(dynamic)
do i = 1,ndim
   vector(i) = 1.0D0/(i*i)
enddo
!$omp end parallel do

xsum(1) = vector(1)

!$omp parallel
!$omp barrier
!$omp master
write(6,*) 'threads',omp_get_num_threads()
!$omp end master
!$omp barrier
!$omp do schedule(dynamic) private(j,n)
do i = 1,sumsize
   n=2**i
   do j=n/2+1,n
      xsum(i) = xsum(i) + vector(j)
   enddo
enddo
!$omp end do
!$omp end parallel

do i = 2,sumsize
   xsum(i) = xsum(i) + xsum(i-1)
   xerror(i) = sinf - xsum(i)
enddo

write(6,*) 'time',omp_get_wtime()-wstart

!Output
do i = 3,sumsize
   write(6,*) xsum(i), xerror(i), (xsum(i) + xerror(i))
enddo

end program
