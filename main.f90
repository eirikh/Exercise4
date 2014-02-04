program main

implicit none

integer :: i,j,n
integer, parameter, kind(=8):: ndim = 16384
real, kind(=8):: vector(ndim)
real, kind(=8):: xsum(12)

!write(*,*) "Enter size of vector v"
!read(*,*) ndim

do i = 1,ndim
   vector(i) = 1.0D0/(i*i)
enddo

do i = 3,14
   n=2**i
   do j=1,n
      xsum(i) = xsum(i) + vector(j)
   enddo
enddo


end program
