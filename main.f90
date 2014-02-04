program main

integer, kind(=8):: ndim
real, kind(=8):: vector(:)
real, kind(=8):: xsum(:)

write(*,*) "Enter size of vector v"
read(*,*) ndim

do i = 1,ndim
   vector(i) = 1.0D0/(i*i)
   xsum(i) = xsum(i) + vector(i)



end program
