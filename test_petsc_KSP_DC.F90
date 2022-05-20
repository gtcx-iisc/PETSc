program main
#include <petsc/finclude/petscvec.h>
#include <petsc/finclude/petscmat.h>
#include <petsc/finclude/petscksp.h>

use petscvec
use petscmat
use petscksp

implicit none

  Vec             x,b
  Mat             A
  KSP             ksp
  PetscInt,parameter ::  m=8
  PetscScalar,parameter ::  two =2.0, one = 1.0
  PetscInt,pointer,dimension(:) ::  dnnz,onnz
  PetscInt    ::  i,rstart,rend,M1,N1
  PetscErrorCode  ierr

  double precision tol

  call PetscInitialize(PETSC_NULL_CHARACTER,ierr)
  
  if (ierr /= 0) then
   print*,'PetscInitialize failed'
   stop
  endif
      
  
  allocate(dnnz(0:m-1))
  allocate(onnz(0:m-1))
  
  do i=0,m-1 
   dnnz(i) = 2
   onnz(i) = 3
  end do
  
  call VecCreate(PETSC_COMM_WORLD,b,ierr);  CHKERRQ(ierr)
  call VecSetSizes(b,PETSC_DECIDE,m,ierr);  CHKERRQ(ierr)
  call VecSetType(b,"mpi",ierr);  CHKERRQ(ierr)
  call VecSetFromOptions(b,ierr);  CHKERRQ(ierr)

  call VecGetOwnershipRange(b,rstart,rend,ierr);  CHKERRQ(ierr)
  call VecDuplicate(b,x,ierr);  CHKERRQ(ierr)

  do i=rstart,rend-1
     call VecSetValue(b,i,one*i*0.5,INSERT_VALUES,ierr);  CHKERRQ(ierr)
  end do

  call VecAssemblyBegin(b,ierr);  CHKERRQ(ierr)
  call VecAssemblyEnd(b,ierr);    CHKERRQ(ierr)

  ! call VecDuplicate(x,b,ierr);  CHKERRQ(ierr)
  call VecView(b,PETSC_VIEWER_STDOUT_WORLD,ierr);  CHKERRQ(ierr)


  call MatCreateAIJ(PETSC_COMM_WORLD,PETSC_DETERMINE,PETSC_DETERMINE,m,m,&
          PETSC_DECIDE,dnnz,PETSC_DECIDE,onnz,A,ierr);  CHKERRQ(ierr)
  call MatSetFromOptions(A,ierr);  CHKERRA(ierr)
  call MatSetUp(A,ierr);  CHKERRA(ierr)
  deallocate(dnnz)
  deallocate(onnz)

  !call MatSetOption(A,MAT_NEW_NONZERO_ALLOCATION_ERR,PETSC_TRUE,&
  !        ierr);  CHKERRQ(ierr)

  call MatSetOption(A,MAT_NEW_NONZERO_ALLOCATION_ERR,PETSC_FALSE,&
           ierr);  CHKERRQ(ierr)
  call MatGetOwnershipRange(A,rstart,rend,ierr);  CHKERRQ(ierr)
  call MatGetSize(A,M1,N1,ierr);  CHKERRQ(ierr)

  do i=rstart,rend-1
   call MatSetValue(A,i,i,one,INSERT_VALUES,ierr);  CHKERRQ(ierr)

   !if (i==0) call MatSetValue(A,i,m-1,two,INSERT_VALUES,ierr);  CHKERRQ(ierr)

   !if (i==m-1) call MatSetValue(A,i,0,two,INSERT_VALUES,ierr);  CHKERRQ(ierr)

   if (i/=m-1) call MatSetValue(A,i,i+1,two,INSERT_VALUES,ierr);  CHKERRQ(ierr)

   if (i/=0)   call MatSetValue(A,i,i-1,two,INSERT_VALUES,ierr);  CHKERRQ(ierr)
   
  end do

  call MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY,ierr);  CHKERRQ(ierr)
  call MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY,ierr);  CHKERRQ(ierr)
  call MatView(A,PETSC_VIEWER_STDOUT_WORLD,ierr);  CHKERRQ(ierr)
  
  call KSPCreate(PETSC_COMM_WORLD,ksp,ierr);  CHKERRQ(ierr)
  call KSPSetType(ksp,KSPCG,ierr);  CHKERRQ(ierr)
  call KSPSetOperators(ksp,A,A,ierr);  CHKERRQ(ierr)
  
  tol = 0.00001
  call KSPSetTolerances(ksp,tol,PETSC_DEFAULT_REAL,&
       PETSC_DEFAULT_REAL,PETSC_DEFAULT_INTEGER,ierr);  CHKERRQ(ierr)
  call KSPSetFromOptions(ksp,ierr)
  call KSPView(ksp,PETSC_VIEWER_STDOUT_WORLD,ierr);  CHKERRQ(ierr)
 
  call KSPSetUp(ksp,ierr);  CHKERRQ(ierr)
  call KSPSolve(ksp,b,x,ierr);  CHKERRQ(ierr)

  call VecView(x,PETSC_VIEWER_STDOUT_WORLD,ierr);  CHKERRQ(ierr)

  call MatDestroy(A,ierr);  CHKERRQ(ierr)
  call PETSCFinalize(ierr);  CHKERRQ(ierr)


end program

!/*TEST
!
!   test:
!      suffix: 1
!      output_file: output/ex4_1.out
!
!   test:
!      suffix: 2
!      nsize: 2
!      output_file: output/ex4_2.out
!
!TEST*/
