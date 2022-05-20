PROGRAM Test_MPI_BC
include 'mpif.h'

integer ProcessID, NumProcess, ierror
real :: r2d(3,2),r1d(6)

call MPI_INIT(ierror)
call MPI_COMM_SIZE(MPI_COMM_WORLD, NumProcess, ierror)
call MPI_COMM_RANK(MPI_COMM_WORLD, ProcessID, ierror)

r = 0.0

if( ProcessID==0 )  then
   open (1, file = 'data.dat', status = 'old')

   do i = 1,3  
      read(1,*) r2d(i,1), r2d(i,2)
   end do 
endif

write(*,*)'ProcessID:',ProcessID,'Before value',r2d

r1d = RESHAPE( r2d , (/6/) )

call MPI_BCAST(r1d,6,MPI_REAL,0,MPI_COMM_WORLD,ierror)

r2d = RESHAPE( r1d , (/3,2/) )

write(*,*)'ProcessID:',ProcessID,'After value',r2d

call MPI_FINALIZE(ierror)
END PROGRAM
