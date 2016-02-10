MODULE DIAGONAL_2DSEM_PRECONDITIONER
!  DIAGONAL_2DSEM_PRECONDITIONER.f90
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
! The model settings
USE MODEL_PRECISION
USE CONSTANTS_DICTIONARY
USE MODEL_FLAGS




IMPLICIT NONE

   TYPE DIAG_SEM_PRECONDITIONER
      integer                 :: nDeg, nElems ! the polynomial degree and the number of elements 
      real(prec), allocatable :: coeffs(:,:,:) ! (0:nDeg,0:nDeg)

      CONTAINS

      PROCEDURE :: BUILD => BUILD_DIAG_SEM_PRECONDITIONER
      PROCEDURE :: TRASH => TRASH_DIAG_SEM_PRECONDITIONER

      PROCEDURE :: SOLVE => SOLVE_DIAG_SEM_PRECONDITIONER

   END TYPE DIAG_SEM_PRECONDITIONER



CONTAINS

 SUBROUTINE BUILD_DIAG_SEM_PRECONDITIONER( myDiag, polyDeg, nEl )
 ! S/R BUILD_DIAG_SEM_PRECONDITIONER
 !
 ! Here, a blank preconditioner is built.
 ! The module CGSM2D_CLASS.f90 will use calls to the MATRIX_ACTION
 ! routine to build the diagonal coefficients.
 ! 
 ! ====================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS( DIAG_SEM_PRECONDITIONER ), intent(out) :: myDiag
   integer, intent(in)                           :: polyDeg, nEl
   ! Local
 



      ALLOCATE( myDiag % coeffs(0:polyDeg,0:polyDeg,1:nEl) )


      myDiag % coeffs = ONE ! Default to one ( the identity matrix, like using no pre-conditioner at all )

      myDiag % nDeg = polyDeg ! Assign the polynomial degree

      myDiag % nElems = nEl
      
 END SUBROUTINE BUILD_DIAG_SEM_PRECONDITIONER
!
!
!
 SUBROUTINE TRASH_DIAG_SEM_PRECONDITIONER( myDiag )
 ! S/R TRASH_DIAG_SEM_PRECONDITIONER
 !
 ! 
 ! ====================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS( DIAG_SEM_PRECONDITIONER ), intent(inout) :: myDiag


      DEALLOCATE( myDiag % coeffs )


 END SUBROUTINE TRASH_DIAG_SEM_PRECONDITIONER
!
!
!
 SUBROUTINE SOLVE_DIAG_SEM_PRECONDITIONER( myDiag, rhs, x)
 ! S/R SOLVE_DIAG_SEM_PRECONDITIONER
 ! 
 !  The preconditioner is diagonal, so we solve by dividing each entry of rhs by
 !  of the coefficients.
 !
 ! ================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS( DIAG_SEM_PRECONDITIONER ), intent(in) :: myDiag
   real(prec), intent(in)                       :: rhs(0:myDiag % nDeg, 0:myDiag % nDeg, 1:myDiag % nElems )
   real(prec), intent(out)                      :: x(0:myDiag % nDeg, 0:myDiag % nDeg, 1:myDiag % nElems)
   ! LOCAL 
   integer    :: iS, iP, iEl, nDeg

      nDeg = myDiag % nDeg   
      
      do iEl = 1, myDiag % nElems
         do iP = 0, nDeg
            do iS = 0, nDeg
  
               x(iS,iP,iEl) = rhs(iS,iP,iEl)/myDiag % coeffs(iS,iP,iEl)
 
            enddo
         enddo
      enddo
  


 END SUBROUTINE SOLVE_DIAG_SEM_PRECONDITIONER
!
!
!
END MODULE DIAGONAL_2DSEM_PRECONDITIONER
