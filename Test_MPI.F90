PROGRAM hello_world_mpi
include 'mpif.h'

integer ProcessID, NumProcess, ierror

call MPI_INIT(ierror)
call MPI_COMM_SIZE(MPI_COMM_WORLD, NumProcess, ierror)
call MPI_COMM_RANK(MPI_COMM_WORLD, ProcessID, ierror)

write(*,*)'ProcessID = ',ProcessID, 'NumProcess =',NumProcess

call MPI_FINALIZE(ierror)
END PROGRAM
