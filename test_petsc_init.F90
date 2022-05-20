program test

!#include <petsc/finclude/petsc.h>
!  use petsc

#include <petsc/finclude/petscksp.h>
  use petscksp
  PetscErrorCode :: ierr
  PetscMPIInt :: myrank
  Character(len=100) :: String

  call PetscInitialize(PETSC_NULL_CHARACTER,ierr)
  call MPI_Comm_rank(PETSC_COMM_WORLD,myrank,ierr);  CHKERRQ(ierr)
  call PetscPrintf(PETSC_COMM_WORLD,"Hello World \n",ierr) ;  CHKERRQ(ierr)

  write(String,*)'MyRank: ',myrank,'\n'
  call PetscPrintf(PETSC_COMM_WORLD,String,ierr);  CHKERRQ(ierr)

  call PetscSynchronizedPrintf(PETSC_COMM_WORLD,String,ierr)
  call PetscSynchronizedFlush(PETSC_COMM_WORLD,PETSC_STDOUT,ierr);  CHKERRQ(ierr)

  call PetscFinalize(ierr)

end program test

