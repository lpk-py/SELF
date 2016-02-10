MODULE CGSM_2D_FEM_PRECONDITIONER
!  CGSM_2D_FEM_PRECONDITIONER.f90
!
!  This module defines a data-type for a preconditioner to the matrix system
!  arrived at in CGSM_2D_CLASS.f90. The preconditioner can be used to accelerate
!  the convergence of the conjugate gradient method used to solve elliptic or parabolic
!  equations implicitly.
!
!  The preconditioner matrix in this module is derived from the Finite Element approximation
!  using the global quadrature nodes as the grid. To build the preconditioner, only the geometry is needed.
! 
!  
!
!
! Module History 
! 
! o  May 2015
!    
!
! 
! List of Module Procedures
!
!
!
!
!  
! =======================================================================================
USE COMMONDATA 

USE LEGENDRE
USE LAGRANGE_1D_CLASS
USE LAGRANGE_2D_CLASS

USE NODAL_STORAGE_2D_CLASS

USE MAPPEDGEOM_2D_CLASS




IMPLICIT NONE

   TYPE FEM1D_PRECONDITIONER
      integer                 :: nDeg, nElems ! the polynomial degree and the number of elements 
      real(prec), allocatable :: dx(:,:)  ! dimensioned (1:nDeg,1:nElems)
      real(prec), allocatable :: coeffs(:,:) ! dimensioned (-1:1,1:nTotNodes) where nTotNodes = (nDeg+1)*(nElems) - 1

      CONTAINS

      PROCEDURE :: BUILD => BUILD_FEM1D_PRECONDITIONER
      PROCEDURE :: TRASH => TRASH_FEM1D_PRECONDITIONER

      PROCEDURE :: SOLVE => SOLVE_FEM1D_PRECONDITIONER

   END TYPE FEM1D_PRECONDITIONER



CONTAINS

 SUBROUTINE BUILD_FEM1D_PRECONDITIONER( myFEM, theMesh, polyDeg )
 ! S/R BUILD_FEM1D_PRECONDITIONER
 !
 ! 
 ! ====================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS( FEM1D_PRECONDITIONER ), intent(out) :: myFEM
   TYPE( MESH_1D ), intent(in)                :: theMesh
   integer, intent(in)                        :: polyDeg
   ! Local
   integer    :: iNode, iEl, nEl, nTotNodes, iGlob
   real(prec) :: dxL, dxR

      nEl = theMesh % nElems

      ALLOCATE( myFEM % dx(1:polyDeg,1:nEl) )

      nTotNodes = (polyDeg)*(nEl)

      ALLOCATE( myFEM % coeffs(-1:1, 0:nTotNodes) )


      myFEM % coeffs = ZERO

      myFEM % nElems = nEl ! Assign the total number of elements
      myFEM % nDeg = polyDeg ! Assign the polynomial degree


      do iEl = 1, theMesh % nElems ! Loop over all of the elements

         do iNode = 1, polyDeg  ! Loop over the polynomial degree


            myFEM % dx(iNode, iEl) = theMesh % elements(iEl) % geometry % x(iNode) - &
                                     theMesh % elements(iEl) % geometry % x(iNode-1)

            
          
         enddo ! iNode, cycle over the polynomial degree

      enddo ! iEl, cycle over the elements

      ! Here is where we would store the FEM-coefficients, if we did.

      do iEl = 1, theMesh % nElems ! Loop over all of the elements

         do iNode = 1, polyDeg  ! Loop over the polynomial degree


            iGlob = iNode + (polyDeg)*(iEl-1) ! Calculate the global node ID

            dxL = myFEM % dx(iNode,iEl) ! grid spacing to the left

            if( iNode == polyDeg ) then ! we need to use the grid spacing in the next element

               if( iEl /= theMesh % nElems )then ! we have made certain not to be be at the very last point for this calculation
               
                  dxR = myFEM % dx(1,iEl+1)

               endif

            else

               dxR = myFEM % dx(iNode+1, iEl)

            endif


            if( iGlob /= nTotNodes )then

               myFEM % coeffs(-1,iGlob) = -ONE/dxL ! Left coefficient
               
               myFEM % coeffs(0,iGlob) = (ONE/dxL + ONE/dxR) ! central coefficient
   
               myFEM % coeffs(1,iGlob) = -ONE/dxR ! right coefficient 

            endif

         enddo ! iNode, cycle over the polynomial degree

      enddo ! iEl, cycle over the elements

 END SUBROUTINE BUILD_FEM1D_PRECONDITIONER
!
!
!
 SUBROUTINE TRASH_FEM1D_PRECONDITIONER( myFEM )
 ! S/R TRASH_FEM1D_PRECONDITIONER
 !
 ! 
 ! ====================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS( FEM1D_PRECONDITIONER ), intent(inout) :: myFEM


      DEALLOCATE( myFEM % dx )


 END SUBROUTINE TRASH_FEM1D_PRECONDITIONER
!
!
!
 SUBROUTINE SOLVE_FEM1D_PRECONDITIONER( myFEM, rhs, x, sorFac )
 ! S/R SOLVE_FEM1D_PRECONDITIONER
 ! 
 !  The preconditioner matrix is assumed to be the tridiagonal matrix
 !  formed by approximating the stiffness matrix with a linear finite 
 !  element discretization. A single SSOR sweep is done to invert the 
 !  system
 !
 ! ================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS( FEM1D_PRECONDITIONER ), intent(in) :: myFEM
   real(prec), intent(in)                     :: rhs(0:myFEM % nDeg, &
                                                     1:myFEM % nElems )
   real(prec), intent(out)                    :: x(0:myFEM % nDeg, &
                                                   1:myFEM % nElems )
   real(prec), intent(in)                     :: sorFac 
   ! LOCAL 
   integer    :: iEl, iNode, k, thatNode, thatEl
   integer    :: nEl, nDeg, iGlob, iter
   real(prec) :: s, res( 0:myFEM%nDeg*myFEM%nElems ) 

      nEl = myFEM % nElems
      nDeg = myFEM % nDeg   
      res = ZERO
      ! Initialize the solution guess to zero 
      x = ZERO

  do iTer = 1, 100
      ! Do the forward sweep through SSOR
      do iEl = 1, nEl
         do iNode = 1,nDeg

            iGlob = iNode + (nDeg)*(iEl-1) ! Calculate the global node ID
            s = ZERO
            do k = -1,1 ! Loop over the FEM stencil

               ! Get the local node ID and the element ID for the node in the stencil
               thatNode = iNode + k

               if( thatNode > nDeg )then
                  thatNode = 1   
                  thatEl = iEl + 1
               else
                  thatEl = iEl 
               endif
               !print*, myFEM % coeffs(k, iGlob)
               s = s + myFEM % coeffs(k, iGlob)*x(thatNode, thatEl)

            enddo ! k, loop over the FEM stencil

            if( iNode == nDeg .AND. iEl == nEl )then ! the last node

               x(iNode, iEl) = ZERO

            else
      
               x(iNode, iEl) = x(iNode, iEl) + sorFac*( rhs(iNode, iEl) - s )/(myFEM % coeffs(0,iGlob))


            endif

         enddo ! Loop over the local nodes

      enddo ! Loop over the elements

      ! Do the backward sweep through SSOR
      do iEl = nEl, 1, -1
         do iNode = nDeg, 1, -1

            iGlob = iNode + (nDeg)*(iEl-1) ! Calculate the global node ID
            s = ZERO
            do k = -1,1 ! Loop over the FEM stencil

               ! Get the local node ID and the element ID for the node in the stencil
               thatNode = iNode + k

               if( thatNode > nDeg )then
                  thatNode = 1   
                  thatEl = iEl + 1
               else
                  thatEl = iEl 
               endif

               s = s + myFEM % coeffs(k, iGlob)*x(thatNode, thatEl)

            enddo ! k, loop over the FEM stencil

            if( iNode == nDeg .AND. iEl == nEl )then ! the last node

               x(iNode, iEl) = ZERO

            else
      
               x(iNode, iEl) = x(iNode, iEl) + sorFac*( rhs(iNode, iEl) - s )/(myFEM % coeffs(0,iGlob))
               res(iGlob) = abs(rhs(iNode,iEl) - s)
            endif

         enddo ! Loop over the local nodes

      enddo ! Loop over the elements
    

enddo
      ! The solution is left masked  
  


 END SUBROUTINE SOLVE_FEM1D_PRECONDITIONER
!
!
!
END MODULE CGSM_2D_FEM_PRECONDITIONER
