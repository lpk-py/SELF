! IterativeSolvers_Class.f90 ( new with v2.1 - 7 Feb. 2016)
! 
! ====================================== LICENSE ================================================= !
!
!  This program is free software; you can redistribute it and/or modify
!  it under the terms of the GNU General Public License as published by
!  the Free Software Foundation; either version 2 of the License, or
!  (at your option) any later version.
!  
!  This program is distributed in the hope that it will be useful,
!  but WITHOUT ANY WARRANTY; without even the implied warranty of
!  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!  GNU General Public License for more details.
!  
!  You should have received a copy of the GNU General Public License
!  along with this program; if not, write to the Free Software
!  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
!  MA 02110-1301, USA.
! 
! ==================================== Module History ============================================ ! 
! 
! o  (ver 2.1) February 2016
!
! ========================================= Logs ================================================= !
! 2016-02-07  Joe  <joe@clay>
!
!   The idea behind this module is to provide a set of iterative solvers that can be used to solve
!   Ax = b, without having to store the matrix "A", and to provide routines that can be reused
!   with the end-user determining the memory layout and data-structure organization that is appropriate
!   for their problem. The solution data, and the right-hand-side are stored in an arbitrary data 
!   structure (of the programmers choosing). The data-structure
!   must be a class that has routines called 
!      (1) MatrixAction -- compute Ax
!      (2) Residual -- computes r = b-Ax
!      (3) PreconditionSolve -- Optionally, a preconditioner can be provided to solve Hz = r, where
!                               H is the preconditioning matrix
!
! To accomplish this, an abstract data-type is provided here with the first two of the above listed
! procedures as "Deferred" procedures --> When the end-user wishes to use this module, these deferred
! procedures must be defined in the extended data type.
!
! 
! Further, the class must have the 
! 
! //////////////////////////////////////////////////////////////////////////////////////////////// !

MODULE IterativeSolvers_Class

USE ModelPrecision
USE ConstantsDictionary

IMPLICIT NONE

  TYPE, ABSTRACT :: IterativeSolver
     INTEGER                 :: maxIters  ! Maximum number of iterates
     INTEGER                 :: nDOF      ! The number of Degrees of Freedom
     REAL(prec)              :: tolerance ! tolerance for convergence
     REAL(prec), ALLOCATABLE :: resi(:)   ! An array for storing the residual at each iterate
     REAL(prec), ALLOCATABLE :: sol(:)
     
     CONTAINS
     
     PROCEDURE :: Initialize => Initialize_IterativeSolver
     PROCEDURE :: Finallize => Finallize_IterativeSolver
     
     PROCEDURE (ActionOfMatrix), DEFERRED :: MatrixAction
     PROCEDURE (BminusAx), DEFERRED       :: Residual
     PROCEDURE (CopyFromType), DEFERRED   :: CopySolution
     
  END TYPE IterativeSolver

  ABSTRACT INTERFACE
     FUNCTION ActionOfMatrix( this )
        USE ModelPrecision
        IMPORT IterativeSolver
        CLASS(IterativeSolver) :: this
        REAL(prec)             :: ActionOfMatrix(1:this % nDOF)
     END FUNCTION ActionOfMatrix
  END INTERFACE
  
  ABSTRACT INTERFACE
     FUNCTION BminusAx( this )
        USE ModelPrecision
        IMPORT IterativeSolver
        CLASS(IterativeSolver) :: this
        REAL(prec)             :: BminusAx(1:this % nDOF)
     END FUNCTION BminusAx
  END INTERFACE

  ABSTRACT INTERFACE
     SUBROUTINE CopyFromType( this, sol )
        USE ModelPrecision
        IMPORT IterativeSolver
        CLASS(IterativeSolver) :: this
        REAL(prec)             :: sol(1:this % nDOF)
     END SUBROUTINE CopyFromType
  END INTERFACE

  TYPE, ABSTRACT, EXTENDS(IterativeSolver) :: ConjugateGradient
     CONTAINS
     PROCEDURE :: Solve => Solve_ConjugateGradient
  END TYPE ConjugateGradient

 CONTAINS
!
!
! 
 SUBROUTINE Initialize_IterativeSolver( this, maxIters, nDOF, tolerance )
 ! S/R Initialize
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS( IterativeSolver ), INTENT(out) :: this
   INTEGER, INTENT(in)                   :: maxIters, nDOF
   REAL(prec), INTENT(in)                :: tolerance
   
      this % maxIters  = maxIters
      this % nDOF      = nDOF
      this % tolerance = tolerance
      
      ALLOCATE( this % sol(1:nDOF), this % resi(0:maxIters) )
      this % sol  = ZERO
      this % resi = ZERO
      
 END SUBROUTINE Initialize_IterativeSolver
!
!
!
 SUBROUTINE Finallize_IterativeSolver( this )
 ! S/R Finallize
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS( IterativeSolver ), INTENT(inout) :: this
      
    DEALLOCATE( this % sol, this % resi )
      
 END SUBROUTINE Finallize_IterativeSolver
!
!
!
!
 SUBROUTINE Solve_ConjugateGradient( this, ioerr )
 !  S/R Solve
 !
 !  This subroutine solves the system Ax = b using the un-preconditioned conjugate gradient method.
 !  The matrix action and residual routines are supplied by a non-abstracted type-extension of
 !  ConjugateGradient. These routines should return an array indexed from 1 to nDOF. Thus,
 !  in addition to a MatrixAction and Residual, the user should map their data-structure to a 1-D 
 !  array.  
 !
 !  On output ioerr is set to an error checking flag. 
 !  If ioerr ==  0, the method converged within the maximum number of iterations.
 !     ioerr == -1, the method did not converge within the maximum number of iterations.
 !     ioerr == -2, something that is not caught by the current construct happened.
 !  
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS( ConjugateGradient ), INTENT(inout) :: this
   INTEGER, INTENT(out)                      :: ioerr
   ! LOCAL
   INTEGER    :: iter
   INTEGER    :: nIt
   REAL(prec) :: TOL
   REAL(prec) :: r(1:this % nDOF )
   REAL(prec) :: v(1:this % nDOF )
   REAL(prec) :: z(1:this % nDOF )
   REAL(prec) :: rNorm
   REAL(prec) :: a, b, num, den, r0
   
      ioerr = -2
      nIt = this % maxIters
      TOL = this % tolerance

      r = this % Residual( )
      r0 = DOT_PRODUCT( r, r ) 
      
      this % resi = ZERO
      this % resi(0) = r0

      ! Apply the preconditioner
      r = z
      num = DOT_PRODUCT( r, z ) ! numerator   r (DOT) (H^(-1)r )
      v = z 

      DO iter = 1,nIt ! Loop over the PCG iterates
 
         ! Compute Ad matrix-vector product
         z = this % MatrixAction( )
         
         ! Compute the search-direction magnitude
         den = DOT_PRODUCT( v, z ) ! denominator

         a = num/den

         ! Update the solution guess
         this % sol = this % sol + a*v

         ! Update the residual
         r = r - a*z
         
         rNorm = DOT_PRODUCT( r, r ) 

         this % resi(iter) = sqrt(rNorm)
         
         IF( sqrt(rNorm)/r0 < TOL ) then
            EXIT
            ioerr = 0
         ENDIF

         ! Apply the preconditioner to the residual
         r = z

         den = num ! r(DOT)[ H^(-1)r ] = r(DOT)d

         ! Calculate the change in the search direction
         num = DOT_PRODUCT( r, z ) 
         
         ! Update the search direction
         b = num/den

         v = z + b*v
         
      ENDDO ! iter, loop over the CG iterates 


      IF( SQRT(rNorm) > TOL )THEN
         PRINT*, 'MODULE IterativeSolvers : ConjugateGradient failed to converge '
         PRINT*, 'Last L-2 residual : ', sqrt(rNorm)
         ioerr=-1
      ENDIF   

 END SUBROUTINE Solve_ConjugateGradient
  
END MODULE IterativeSolvers_Class
