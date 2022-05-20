program main
#include <petsc/finclude/petscvec.h>
#include <petsc/finclude/petscmat.h>

use petscvec
use petscmat

implicit none

  Mat             A
  PetscInt,parameter ::  m=8
  PetscScalar,parameter ::  two =2.0, one = 1.0
  PetscInt,pointer,dimension(:) ::  dnnz,onnz
  PetscInt    ::  i,rstart,rend,M1,N1
  PetscErrorCode  ierr

  call PetscInitialize(PETSC_NULL_CHARACTER,ierr)
  
  if (ierr /= 0) then
   print*,'PetscInitialize failed'
   stop
  endif
      
  
  allocate(dnnz(0:m-1))
  allocate(onnz(0:m-1))
  
  do i=0,m-1 
   dnnz(i) = 1
   onnz(i) = 3
  end do
  
  call MatCreateAIJ(PETSC_COMM_WORLD,PETSC_DETERMINE,PETSC_DETERMINE,m,m,&
          PETSC_DECIDE,dnnz,PETSC_DECIDE,onnz,A,ierr);  CHKERRQ(ierr)
  call MatSetFromOptions(A,ierr);  CHKERRA(ierr)
  call MatSetUp(A,ierr);  CHKERRA(ierr)
  deallocate(dnnz)
  deallocate(onnz)

  call MatSetOption(A,MAT_NEW_NONZERO_ALLOCATION_ERR,PETSC_TRUE,&
          ierr);  CHKERRQ(ierr)
  call MatGetOwnershipRange(A,rstart,rend,ierr);  CHKERRQ(ierr)
  call MatGetSize(A,M1,N1,ierr);  CHKERRQ(ierr)

  do i=rstart,rend-1
   call MatSetValue(A,i,i,two,INSERT_VALUES,ierr);  CHKERRQ(ierr)

   if (rend<N1) then
      call MatSetValue(A,i,m-1,one,INSERT_VALUES,ierr);  CHKERRQ(ierr)
      call MatSetValue(A,i,m-2,one,INSERT_VALUES,ierr);  CHKERRQ(ierr)
   end if 

  end do
  call MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY,ierr);  CHKERRQ(ierr)
  call MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY,ierr);  CHKERRQ(ierr)
  call MatView(A,PETSC_VIEWER_STDOUT_WORLD,ierr);  CHKERRQ(ierr)
  call MatDestroy(A,ierr);  CHKERRQ(ierr)
  call PetscFinalize(ierr);  CHKERRQ(ierr)
  
end program

