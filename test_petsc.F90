program test

!#include <petsc/finclude/petsc.h>
!  use petsc

#include <petsc/finclude/petscksp.h>
  use petscksp
  
  implicit none

  PetscErrorCode :: ierr
  PetscMPIInt :: myrank
  PetscInt :: N,ilow,ihigh,i
  Character(len=100) :: String
  !Integer :: i

  Vec v
  Mat K

  N = 100

  call PetscInitialize(PETSC_NULL_CHARACTER,ierr)
  call MPI_Comm_rank(PETSC_COMM_WORLD,myrank,ierr);  CHKERRQ(ierr)

  write(String,*)'MyRank: ',myrank,'\n'
  call PetscPrintf(PETSC_COMM_SELF,String,ierr);  CHKERRQ(ierr)

  call VecCreate(PETSC_COMM_SELF,v,ierr);  CHKERRQ(ierr)
  call VecSetSizes(v,PETSC_DECIDE,N,ierr);  CHKERRQ(ierr)
  call VecSetType(v,"mpi",ierr);  CHKERRQ(ierr)
  call VecSetFromOptions(v,ierr);  CHKERRQ(ierr)
 
  call VecGetOwnershipRange(v,ilow,ihigh,ierr);  CHKERRQ(ierr)

  do i=ilow,ihigh-1
     write(*,*)i,'ilow',ilow,'ihigh',ihigh
     call VecSetValue(v,i,i,INSERT_VALUES,ierr);  CHKERRQ(ierr)
  end do

  call VecAssemblyBegin(v,ierr);  CHKERRQ(ierr)
  call VecAssemblyEnd(v,ierr);    CHKERRQ(ierr)

  call PetscFinalize(ierr)

end program test

