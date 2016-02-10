MODULE DIAGONAL_2D_PRECONDITIONER
!  DIAGONAL_2D_PRECONDITIONER.f90
!
!  This module defines a data-type for a preconditioner to the matrix system
!  arrived at in CGSM_2D_CLASS.f90. The preconditioner can be used to accelerate
!  the convergence of the conjugate gradient method used to solve elliptic or parabolic
!  equations implicitly.
!
!  The preconditioner matrix in this module is derived from diagonal contribution to the laplacian
!  operator for the nodal galerkin spectral method in mapped geometry. It is intended to be used
!  as a preconditioner for elliptic or parabolic problems.
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




IMPLICIT NONE

   TYPE DIAG_PRECONDITIONER
      integer                 :: nDeg ! the polynomial degree and the number of elements 
      real(prec), allocatable :: coeffs(:,:) ! (0:nDeg,0:nDeg)

      CONTAINS

      PROCEDURE :: BUILD => BUILD_DIAG_PRECONDITIONER
      PROCEDURE :: TRASH => TRASH_DIAG_PRECONDITIONER

      PROCEDURE :: SOLVE => SOLVE_DIAG_PRECONDITIONER

   END TYPE DIAG_PRECONDITIONER



CONTAINS

 SUBROUTINE BUILD_DIAG_PRECONDITIONER( myDiag, polyDeg )
 ! S/R BUILD_DIAG_PRECONDITIONER
 !
 ! Here, a blank preconditioner is built.
 ! The module CGSM2D_CLASS.f90 will use calls to the MATRIX_ACTION
 ! routine to build the diagonal coefficients.
 ! 
 ! ====================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS( DIAG_PRECONDITIONER ), intent(out) :: myDiag
   integer, intent(in)                       :: polyDeg
   ! Local
   integer    :: iNode, iEl, nEl, nTotNodes, iGlob
   real(prec) :: dxL, dxR



      ALLOCATE( myDiag % coeffs(0:polyDeg,0:polyDeg) )


      myDiag % coeffs = ONE ! Default to one ( the identity matrix, like using no pre-conditioner at all )

      myDiag % nDeg = polyDeg ! Assign the polynomial degree


      
 END SUBROUTINE BUILD_DIAG_PRECONDITIONER
!
!
!
 SUBROUTINE TRASH_DIAG_PRECONDITIONER( myDiag )
 ! S/R TRASH_DIAG_PRECONDITIONER
 !
 ! 
 ! ====================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS( DIAG_PRECONDITIONER ), intent(inout) :: myDiag


      DEALLOCATE( myDiag % coeffs )


 END SUBROUTINE TRASH_DIAG_PRECONDITIONER
!
!
!
 SUBROUTINE SOLVE_DIAG_PRECONDITIONER( myDiag, rhs, x)
 ! S/R SOLVE_DIAG_PRECONDITIONER
 ! 
 !  The preconditioner is diagonal, so we solve by dividing each entry of rhs by
 !  of the coefficients.
 !
 ! ================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS( DIAG_PRECONDITIONER ), intent(in) :: myDiag
   real(prec), intent(in)                   :: rhs(0:myDiag % nDeg, 0:myDiag % nDeg )
   real(prec), intent(out)                  :: x(0:myDiag % nDeg, 0:myDiag % nDeg)
   ! LOCAL 
   integer    :: iS, iP, nDeg

      nDeg = myDiag % nDeg   
      
      do iP = 0, nDeg
         do iS = 0, nDeg
 
            x(iS,iP) = rhs(iS,iP)/myDiag % coeffs(iS,iP)
 
         enddo
      enddo
  


 END SUBROUTINE SOLVE_DIAG_PRECONDITIONER
!
!
!
END MODULE DIAGONAL_2D_PRECONDITIONER
