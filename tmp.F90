program main
#include <petsc/finclude/petscvec.h>
#include <petsc/finclude/petscmat.h>
#include <petsc/finclude/petscksp.h>

use petscvec
use petscmat
use petscksp

implicit none

  Vec             v,b,x_recv
  Mat             A
  KSP             ksp
  PetscInt    ::  m
  PetscScalar ::  ss
  PetscInt,pointer,dimension(:) ::  dnnz,onnz,idx
  PetscScalar,pointer,dimension(:,:) :: blockVal
  PetscInt    ::  i,j,k,rstart,rend,M1,N1
  PetscErrorCode  ierr
  VecScatter     gather
  MPI_Comm       comm

  Integer     :: MaxVertexValency=8
  Integer     :: fxy=10, ftri=11, fni=12, fbv=13
  Integer     :: err,nxy,ntri,nni,nbv,myrank,ii,jj,kk,Flag

  double precision, allocatable, dimension(:) :: ni,phi,x,y,d
  double precision, allocatable, dimension(:,:) :: xy
  Integer, allocatable, dimension(:) :: bv
  Integer, allocatable, dimension(:,:) :: tri

  double precision :: tol, Area

  open(UNIT=fxy ,FILE='xy.out' ,ACTION='READ')
  open(UNIT=ftri,FILE='tri.out',ACTION='READ')
  open(UNIT=fni ,FILE='ni.out' ,ACTION='READ')
  open(UNIT=fbv ,FILE='bv.out' ,ACTION='READ')

  read(fxy,*)nxy
  read(ftri,*)ntri
  read(fni,*)nni
  read(fbv,*)nbv

  allocate(xy(nxy,2))
  allocate(tri(ntri,3))
  allocate(ni(nni))
  allocate(bv(nbv))
  allocate(x(nxy))
  allocate(y(nxy))
  allocate(d(nxy))

  do ii=1,nxy
    read(fxy,*)xy(ii,1),xy(ii,2)
  end do

  x = xy(:,1)
  y = xy(:,2)

  do ii=1,ntri
    read(ftri,*)tri(ii,1),tri(ii,2),tri(ii,3)
  end do

  do ii=1,nni
    read(fni,*)ni(ii)
  end do

  do ii=1,nbv
    read(fbv,*)bv(ii)
  end do

  m = nxy
  comm = MPI_COMM_WORLD

  call PetscInitialize(PETSC_NULL_CHARACTER,ierr)
  
  if (ierr /= 0) then
   print*,'PetscInitialize failed'
   stop
  endif
 
  call MPI_Comm_rank(PETSC_COMM_WORLD,myrank,ierr);  CHKERRQ(ierr)

  if (myrank==0) then
    !write(*,*)xy(:,2)
    !write(*,*)tri(:,1)
    !write(*,*)ni
    !write(*,*)bv
  end if

  allocate(dnnz(0:m-1))
  allocate(onnz(0:m-1))
  allocate(idx(1))
  allocate(blockval(1,1))
  
  do i=0,m-1 
   dnnz(i) = 2
   onnz(i) = MaxVertexValency
  end do
  
  call VecCreateMPI(comm,PETSC_DECIDE,m,b,ierr);  CHKERRQ(ierr)
  call VecSetFromOptions(b,ierr);  CHKERRQ(ierr)

  call VecGetOwnershipRange(b,rstart,rend,ierr);  CHKERRQ(ierr)
  call VecDuplicate(b,v,ierr);  CHKERRQ(ierr)

  call MatCreateAIJ(comm,PETSC_DETERMINE,PETSC_DETERMINE,m,m,&
          PETSC_DECIDE,dnnz,PETSC_DECIDE,onnz,A,ierr);  CHKERRQ(ierr)
  call MatSetFromOptions(A,ierr);  CHKERRA(ierr)
  call MatSetUp(A,ierr);  CHKERRA(ierr)

  deallocate(dnnz)
  deallocate(onnz)

  call MatSetOption(A,MAT_NEW_NONZERO_ALLOCATION_ERR,PETSC_FALSE,&
           ierr);  CHKERRQ(ierr)
  call MatGetOwnershipRange(A,rstart,rend,ierr);  CHKERRQ(ierr)
  call MatGetSize(A,M1,N1,ierr);  CHKERRQ(ierr)

  d = 0.0

  do jj=1,ntri
    do ii=0,2

      i = tri(jj,ii+1)
      j = tri(jj,mod(ii+1,3)+1)
      k = tri(jj,mod(ii+2,3)+1)

      Area = 0.5* ( (x(i)-x(j))*(y(i)-y(k)) - ((x(i)-x(k))*(y(i)-y(j))) )
      d(i) = d(i) + Area/12 * (2*ni(i)+ni(j)+ni(k))

      if( ((i-1)>=rstart) .AND. ((i-1)<rend) )then
        Flag = 1

        do kk=1,nbv
          if (i==bv(kk)) then
            d(i) = 0.0

            !ss = 1.0
            !call MatAssemblyBegin(A,MAT_FLUSH_ASSEMBLY,ierr);  CHKERRQ(ierr)
            !call MatAssemblyEnd(A,MAT_FLUSH_ASSEMBLY,ierr);  CHKERRQ(ierr)
            !call MatSetValue(A,i-1,i-1,ss,INSERT_VALUES,ierr);  CHKERRQ(ierr)
            !call MatAssemblyBegin(A,MAT_FLUSH_ASSEMBLY,ierr);  CHKERRQ(ierr)
            !call MatAssemblyEnd(A,MAT_FLUSH_ASSEMBLY,ierr);  CHKERRQ(ierr)

            Flag = 0
            exit
          end if
        end do

        if (Flag==1) then

          ss = -1/(4*Area) * ( (x(j)-x(k))**2 + (y(j)-y(k))**2 )
          call MatSetValue(A,i-1,i-1,ss,ADD_VALUES,ierr);  CHKERRQ(ierr)

          ss = -1/(4*Area) * ( (x(i)-x(k))*(x(k)-x(j)) + &
                               (y(i)-y(k))*(y(k)-y(j)) )
          call MatSetValue(A,i-1,j-1,ss,ADD_VALUES,ierr);  CHKERRQ(ierr)

          ss = -1/(4*Area) * ( (x(i)-x(j))*(x(j)-x(k)) + &
                               (y(i)-y(j))*(y(j)-y(k)) )
          call MatSetValue(A,i-1,k-1,ss,ADD_VALUES,ierr);  CHKERRQ(ierr)

        end if
      end if

    end do
  end do

  call MatAssemblyBegin(A,MAT_FLUSH_ASSEMBLY,ierr);  CHKERRQ(ierr)
  call MatAssemblyEnd(A,MAT_FLUSH_ASSEMBLY,ierr);  CHKERRQ(ierr)

  do ii=1,nbv
    i = bv(ii)

    if( ((i-1)>=rstart) .AND. ((i-1)<rend) )then
      ss = 1.0
      call MatSetValue(A,i-1,i-1,ss,INSERT_VALUES,ierr);  CHKERRQ(ierr)
    end if
  end do

  call VecGetOwnershipRange(b,rstart,rend,ierr);  CHKERRQ(ierr)

  do ii=1,nbv
    d(bv(ii)) = 0.0
  end do

  do i=1,nxy
    if( ((i-1)>=rstart) .AND. ((i-1)<rend) )then
      ss = d(i)
      call VecSetValue(b,i-1,ss,INSERT_VALUES,ierr);  CHKERRQ(ierr)
    end if
  end do

  call VecAssemblyBegin(b,ierr);  CHKERRQ(ierr)
  call VecAssemblyEnd(b,ierr);    CHKERRQ(ierr)
  !call VecView(b,PETSC_VIEWER_STDOUT_WORLD,ierr);  CHKERRQ(ierr)

  call MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY,ierr);  CHKERRQ(ierr)
  call MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY,ierr);  CHKERRQ(ierr)
  !call MatView(A,PETSC_VIEWER_STDOUT_WORLD,ierr);  CHKERRQ(ierr)
 
  call PetscViewerPushFormat(PETSC_VIEWER_STDOUT_WORLD,&
          PETSC_VIEWER_ASCII_MATLAB,ierr)
  call MatView(A,PETSC_VIEWER_STDOUT_WORLD,ierr)
  call PetscViewerPopFormat(PETSC_VIEWER_STDOUT_WORLD,ierr)

  call KSPCreate(comm,ksp,ierr);  CHKERRQ(ierr)
  call KSPSetOperators(ksp,A,A,ierr);  CHKERRQ(ierr)
  
  tol = 0.0000001
  call KSPSetTolerances(ksp,tol,PETSC_DEFAULT_REAL,&
       PETSC_DEFAULT_REAL,PETSC_DEFAULT_INTEGER,ierr);  CHKERRQ(ierr)
  !call KSPView(ksp,PETSC_VIEWER_STDOUT_WORLD,ierr);  CHKERRQ(ierr)
 
  call KSPSolve(ksp,b,v,ierr);  CHKERRQ(ierr)
  call VecView(v,PETSC_VIEWER_STDOUT_WORLD,ierr);  CHKERRQ(ierr)

  call VecDestroy(v,ierr);  CHKERRQ(ierr)
  call VecDestroy(b,ierr);  CHKERRQ(ierr)
  call MatDestroy(A,ierr);  CHKERRQ(ierr)
  call KSPDestroy(ksp,ierr);  CHKERRQ(ierr)
  call PETSCFinalize(ierr);  CHKERRQ(ierr)

end program

