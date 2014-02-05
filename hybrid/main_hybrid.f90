program main_mpi

implicit none
include 'mpif.h'

integer(kind=8) :: i,j,n,chunk,localstop
integer(kind=8), parameter :: sumsize = 14
integer(kind=8), parameter :: ndim = 2**sumsize 
real(kind=8), parameter :: sinf = (4*datan(1.0d0))**2/6.0D0
real(kind=8):: xerror(sumsize)
real(kind=8):: vector(ndim)
real(kind=8):: xsum(sumsize),xtmp = 0
real(kind=4):: wstart, wstop
integer :: omp_get_num_threads, omp_get_thread_num
integer :: ierr,mpisize,mpirank,mpistatus(mpi_status_size), log2

call second(wstart)

call mpi_init(ierr)

call mpi_comm_size(mpi_comm_world, mpisize, ierr)
call mpi_comm_rank(mpi_comm_world, mpirank, ierr)

data xsum/ sumsize * 0.0D0/
chunk=ndim/mpisize

if (mpirank .eq. 0) then
!$omp parallel do schedule(static)
   do i = 1,ndim
      vector(i) = 1.0D0/(i*i)
   enddo
!$omp end parallel do

   xsum(1) = vector(1)

   do i=1,mpisize-1
      call mpi_send(vector(i*chunk+1),chunk,mpi_real8,i,1,mpi_comm_world,ierr)
   enddo
else
   call mpi_recv(vector(1),chunk,mpi_real8,0,1,mpi_comm_world,mpistatus,ierr)
endif


if (mpirank .eq. 0) then
!$omp parallel do schedule(static) private(j,n)
   do i = 1,sumsize-log2(mpisize)
      n=2**i
      do j=n/2+1,n
         xsum(i) = xsum(i) + vector(j)
      enddo
   enddo
!$omp end parallel do
else
   i=sumsize - log2(mpisize) +log2(mpirank) + 1
!$omp parallel do schedule(static) reduction(+:xtmp)
   do j=1,chunk
      xtmp = xtmp + vector(j)
   enddo
!$omp end parallel do
   xsum(i) = xtmp
endif

if (mpirank .eq. 0) then
   call mpi_reduce(mpi_in_place,xsum(1),sumsize,mpi_real8,mpi_sum,0,mpi_comm_world,ierr)
else
   call mpi_reduce(xsum(1),xsum(1),sumsize,mpi_real8,mpi_sum,0,mpi_comm_world,ierr)
endif

if (mpirank .eq. 0) then

   do i = 2,sumsize
      xsum(i) = xsum(i) + xsum(i-1)
      xerror(i) = sinf - xsum(i)
   enddo

   !Output
   call second(wstop)
   write(6,*) 'time',wstop-wstart
   if (mpirank .eq. 0) then
      do i = 3,sumsize
         write(*,*) xsum(i), xerror(i), (xsum(i) + xerror(i))
      enddo
   endif

endif

call mpi_finalize(ierr)

end program

function log2(n) result(i)

integer, intent(in) :: n
integer             :: i,rest

rest = n

i = -1

do while(rest .gt. 0)
   i = i+1
   rest = rest/2
enddo

end function log2
