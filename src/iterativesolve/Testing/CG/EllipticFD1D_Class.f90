! EllipticFD1D_Class.f90 ( new with v2.1 - 11 Feb. 2016)
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
! 2016-02-11  Joe  <joe@clay>
!
!   This module provides a simple example that shows how to extend the iterative solvers for a 
!   particular problem. Here, the problem is to solve u_xx = s(x) subject to homogeneous dirichlet
!   boundary conditions. The 1-D elliptic pde is discretized on a uniform grid with centered finite
!   differences. This data-structure provides routines for the matrix-action and residual 
!   calculations that are necessary for using the conjugate gradient solver.
! 
! //////////////////////////////////////////////////////////////////////////////////////////////// !

MODULE EllipticFD1D_Class

USE ModelPrecision
USE ConstantsDictionary
!
USE IterativeSolvers_Class

   ! The type extension defined here (EllipticFD1D) inherits the attributes and routines from the 
   ! IterativeSolvers and ConjugateGradient classes. This includes the attributes
   !
   !  INTEGER                 :: maxIters  ! Maximum number of iterates
   !  INTEGER                 :: nDOF      ! The number of Degrees of Freedom
   !  REAL(prec)              :: tolerance ! tolerance for convergence
   !  REAL(prec), ALLOCATABLE :: resi(:)   ! An array for storing the residual at each iterate
   !  REAL(prec), ALLOCATABLE :: sol(:)
   !
   !  manual constructor and destructor 
   !
   !  PROCEDURE :: Initialize
   !  PROCEDURE :: Finallize
   !
   !  and the solver
   !
   !  PROCEDURE :: Solve
   !
   !  Any procedure can be defined in this module to override the inherited routines. For this 
   !  1-D elliptic problem, this is not necessary. In this module, however, we must define the 
   !  MatrixAction, Residual, and CopySolution routines.
   !
   !  The latter is rather trivial for this problem. The 1-D array storage is already appropriate 
   !  for this problem.
   !
   TYPE, EXTENDS( ConjugateGradient ) :: EllipticFD1D 
      INTEGER                 :: nX
      REAL(prec)              :: dx, diffCoeff
      REAL(prec), ALLOCATABLE :: x(:), localSol(:), flux(:,:)
      
      PROCEDURE :: Initialize => Initialize_EllipticFD1D
      
      
   END TYPE EllipticFD1D
!
!
!==================================================================================================!
!------------------------------- Manual Constructors/Destructors ----------------------------------!
!==================================================================================================!
!
!
 SUBROUTINE Initialize_EllipticFD1D( this, maxIters, tolerance, xL, xR, dx, dCo )
 ! S/R Initialize
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS( EllipticFD1D ), INTENT(out) :: this
   INTEGER, INTENT(in)                :: maxIters
   REAL(prec), INTENT(in)             :: tolerance, xL, xR, dx, dCo
   ! LOCAL
   INTEGER :: nDOF, nX, i
   
      this % maxIters  = maxIters
      this % nX        = (xR-xL-dx)/dx
      this % nDOF      = this % nX + 1
      this % tolerance = tolerance
      
      this % dx = dx
      this % diffCoeff = dCo
      
      nDOF = this % nDOF
      nX   = this % nX
      
      ALLOCATE( this % sol(1:nDOF), this % resi(0:maxIters), this % x(0:nX), this % localSol(0:nX) )
      ALLOCATE( this % flux(0:nX,1:2) )
      this % sol      = ZERO
      this % resi     = ZERO
      this % localSol = ZERO
      this % flux     = ZERO
      
      DO i = 0, nX
         this % x(i) = xL + dx*(REAL(i,prec)) + HALF*dX
      ENDDO
      
 END SUBROUTINE Initialize_EllipticFD1D
!
!
!
 SUBROUTINE Finallize_EllipticFD1D( this )
 ! S/R Finallize
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS( EllipticFD1D ), INTENT(inout) :: this
      
    DEALLOCATE( this % sol, this % resi, this % x, this % localSol, this % flux )
      
 END SUBROUTINE Finallize_EllipticFD1D
!
!
!==================================================================================================!
!----------------------------------------- Type Specific ------------------------------------------!
!==================================================================================================!
!
!
 SUBROUTINE Flux_EllipticFD1D( this )
 !
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS( EllipticFD1D ), INTENT(inout) :: this
   ! LOCAL
   REAL(prec) :: sR, sL, dx, kappa
   INTEGER    :: i, nX
   
    dx = this % dx
    nX = this % nX
    kappa = this % diffCoeff
    
    ! The flux on the left side of the first element is zero (homogeneous neumann)
    this % flux(0,1)  = ZERO
    ! Similarly, we apply zero flux of the right-most boundary of the right-most element.
    this % flux(nX,2) = ZERO
    
    DO i = 0, nX-1
       ! Calculate the flux on the right side of this element
       sL = this % localSol(i)
       sR = this % localSol(i+1)
       this % flux(i,2) = (sR-sL)/dx
       
       ! Global conservation requires that neighboring elements have fluxes that are equal in 
       ! magnitude, but opposite in sign. Because of this, we can set the left hand flux of the next
       ! element equal to the negative of the right hand flux of the current element
       this % flux(i+1,1) = -this % flux(i,2)
    ENDDO
       
       
  
 END SUBROUTINE Flux_EllipticFD1D
!
!
!
 FUNCTION MatrixAction_EllipticFD1D( this )
 !
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   CLASS( EllipticFD1D ) :: this
   REAL(prec)            :: MatrixAction_EllipticFD1D(1:this % nDOF)
   ! LOCAL
   INTEGER    :: i, nX
   REAL(prec) :: dx
   
     dx = this % dx
     nX = this % nX
     
     DO i = 0, nX
     
        MatrixAction_EllipticFD1D(i+1) = this % flux() 
 END FUNCTION MatrixAction_EllipticFD1D
END MODULE EllipticFD1D_Class
