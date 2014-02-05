program hybrid_pi
implicit none
include 'mpif.h'
   !! Program: Calculation of PI.
   integer    :: ierr,myrank,ranks,master_rank
   integer(8) :: i,istart,iend 
   real(8)    :: x, mypi,vpi,gsum,lsum,step,tol
   integer(8) :: nsteps
   real(8)    :: starttime
   
   master_rank=0
   gsum = 0.0
   lsum = 0.0
   vpi = 3.141592654
   tol = 0.0000001
   nsteps = 100000000 !!Try different nsteps

   call MPI_INIT(ierr)
   call MPI_COMM_SIZE(MPI_COMM_WORLD,ranks,ierr)
   call MPI_COMM_RANK(MPI_COMM_WORLD,myrank,ierr)

   starttime = MPI_Wtime()

   step = 1.0 / nsteps

   istart = myrank * nsteps / ranks
   iend   = (myrank+1) * nsteps / ranks

   !omp parallel do private(i,x) reduction(+:lsum)
   do i=istart, iend -1
      x = (i + 0.5) * step;
      lsum = lsum + 4.0/(1.0 + x * x);
   end do
   !omp end parallel do

   call MPI_REDUCE(lsum, gsum, 1 , MPI_REAL8, MPI_SUM, master_rank, MPI_COMM_WORLD, ierr)
   
   mypi = step * gsum;

   if (myrank == master_rank) then
      if (abs(mypi-vpi)>tol) then
         write (6,*) 'Error in pi: ',mypi
      else
         write(6,*) 'PI ', mypi
      end if
      write (6,*) 'Runtime : ',MPI_Wtime()-starttime
   end if


   call MPI_FINALIZE(ierr);

end program hybrid_pi

