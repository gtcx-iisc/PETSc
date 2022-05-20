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
  PetscInt :: row(1),col(1),value(1)
  PetscInt, allocatable, dimension(:) :: dnnz,onnz

  Vec v
  Mat K

  N = 100

  allocate(dnnz(0:N-1))
  allocate(onnz(0:N-1))

  do i=0,(N-1)
     dnnz(i) = 1
     onnz(i) = 1
  end do

  call PetscInitialize(PETSC_NULL_CHARACTER,ierr)
  call MPI_Comm_rank(PETSC_COMM_WORLD,myrank,ierr);  CHKERRQ(ierr)

  write(String,*)'MyRank: ',myrank,'\n'
  call PetscPrintf(PETSC_COMM_SELF,String,ierr);  CHKERRQ(ierr)


  !! VECTOR COMPONENTS

  call VecCreate(PETSC_COMM_WORLD,v,ierr);  CHKERRQ(ierr)
  call VecSetSizes(v,PETSC_DECIDE,N,ierr);  CHKERRQ(ierr)
  call VecSetType(v,"mpi",ierr);  CHKERRQ(ierr)
  call VecSetFromOptions(v,ierr);  CHKERRQ(ierr)
 
  call VecGetOwnershipRange(v,ilow,ihigh,ierr);  CHKERRQ(ierr)

  do i=ilow,(ihigh-1)
     write(*,*)i,'ilow',ilow,'ihigh',ihigh
     call VecSetValue(v,i,i,INSERT_VALUES,ierr);  CHKERRQ(ierr)
  end do

  call VecAssemblyBegin(v,ierr);  CHKERRQ(ierr)
  call VecAssemblyEnd(v,ierr);    CHKERRQ(ierr)
  write(*,*)'OK1'

  !! MATRIX COMPONENT

  ! call MatCreate(PETSC_COMM_WWORLD,K,ierr);  CHKERRQ(ierr)
  ! call MatCreateAIJ(PETSC_COMM_WORLD,PETSC_DECIDE,PETSC_DECIDE,N,N, &
  !     1,PETSC_NULL_INTEGER,1,PETSC_NULL_INTEGER, &
  !     K,ierr);  CHKERRQ(ierr)
 
  call MatCreateAIJ(PETSC_COMM_WORLD,N,N,PETSC_DETERMINE,PETSC_DETERMINE, &
      PETSC_DECIDE,dnnz,PETSC_DECIDE,onnz,K,ierr);  CHKERRQ(ierr)

  write(*,*)'OK2'
  ! call MatSetSizes(K,PETSC_DECIDE,PETSC_DECIDE,N,N,ierr);  CHKERRQ(ierr)
  write(*,*)'OK3'
  !call MatSetType(K,"mpi",ierr);   CHKERRQ(ierr)
  call MatSetFromOptions(K,ierr);   CHKERRQ(ierr)
  write(*,*)'OK4'

  call MatGetOwnershipRange(K,ilow,ihigh,ierr);   CHKERRQ(ierr)
  write(*,*)'OK5'


  ! CREATE A DIAGONAL MATRIX
  do i=ilow,(ihigh-1)
     row(1) = i
     col(1) = i
     value(1) = 10*i
     call MatSetValue(K,1,row,1,col,value,INSERT_VALUES,ierr);   CHKERRQ(ierr)
  end do

  write(*,*)'OK6'
  call MatAssemblyBegin(K,MAT_FINAL_ASSEMBLY,ierr);  CHKERRQ(ierr);
  call MatAssemblyEnd(K,MAT_FINAL_ASSEMBLY,ierr);  CHKERRQ(ierr);
  write(*,*)'OK7'

  call MatDestroy(K,ierr);  CHKERRQ(ierr)
  call PetscFinalize(ierr);  CHKERRQ(ierr)

end program test

