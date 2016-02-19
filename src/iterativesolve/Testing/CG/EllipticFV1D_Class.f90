! EllipticFV1D_Class.f90 ( new with v2.1 - 11 Feb. 2016)
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

MODULE EllipticFV1D_Class

USE ModelPrecision
USE ConstantsDictionary
!
USE IterativeSolvers_Class
USE IterativeSolversParams_Class

IMPLICIT NONE

   ! The type extension defined here (EllipticFV1D) inherits the attributes and routines from the 
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

   
   TYPE, EXTENDS( ConjugateGradient ) :: EllipticFV1D 
      INTEGER                 :: nX
      REAL(prec)              :: dx, diffCoeff
      REAL(prec), ALLOCATABLE :: x(:), localSol(:), flux(:,:)
      REAL(prec), ALLOCATABLE :: source(:)
      
      CONTAINS

      PROCEDURE :: Initialize => Initialize_EllipticFV1D
      PROCEDURE :: Finallize  => Finallize_EllipticFV1D 
   
      PROCEDURE :: CopyTo   => CopyTo_EllipticFV1D
      PROCEDURE :: CopyFrom => CopyFrom_EllipticFV1D
      
      PROCEDURE :: CalculateFlux => CalculateFlux_EllipticFV1D
      PROCEDURE :: MatrixAction  => MatrixAction_EllipticFV1D
      PROCEDURE :: Residual      => Residual_EllipticFV1D

      PROCEDURE :: WriteTecplot => WriteTecplot_EllipticFV1D

   END TYPE EllipticFV1D

 CONTAINS
!
!
!==================================================================================================!
!------------------------------- Manual Constructors/Destructors ----------------------------------!
!==================================================================================================!
!
!
 SUBROUTINE Initialize_EllipticFV1D( this, params )
 ! S/R Initialize
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS( EllipticFV1D ), INTENT(out)         :: this
   TYPE( IterativeSolversParams ), INTENT(in) :: params
   ! LOCAL
   INTEGER :: nDOF, nX, i
   
      this % maxIters  = params % MaximumIterates
      this % dx        = params % dX
      this % nX        = INT( (params % xR - params % xL)/this % dx )
      this % nDOF      = this % nX
      this % tolerance = params % tolerance
      this % diffCoeff = params % kappa
      
      nDOF = this % nDOF
      nX   = this % nX
      
      ALLOCATE( this % sol(1:nDOF), this % resi(0:this % maxIters), this % x(1:nX), this % localSol(0:nX+1) )
      ALLOCATE( this % flux(1:nX+1,1:2), this % source(1:nX) )
      this % sol      = ZERO
      this % resi     = ZERO
      this % localSol = ZERO
      this % flux     = ZERO
      this % source   = ZERO
      
      DO i = 1, nX
         this % x(i) = params % xL + this % dx*(REAL(i-1,prec)) + HALF*this % dX
      ENDDO
      
 END SUBROUTINE Initialize_EllipticFV1D
!
!
!
 SUBROUTINE Finallize_EllipticFV1D( this )
 ! S/R Finallize
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS( EllipticFV1D ), INTENT(inout) :: this
      
    DEALLOCATE( this % sol, this % resi, this % x, this % localSol, this % flux )
      
 END SUBROUTINE Finallize_EllipticFV1D
!
!
!==================================================================================================!
!----------------------------------------- Type Specific ------------------------------------------!
!==================================================================================================!
!
!
 SUBROUTINE CopyFrom_EllipticFV1D( this, sol )
 ! S/R CopyFrom
 !
 !  Copies FROM sol to the Elliptic_FV1D class
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS( EllipticFV1D ), INTENT(inout) :: this
   REAL(prec), INTENT(in)               :: sol(1:this % nDOF)
   
    this % localSol(1:this % nDOF) = sol(1:this % nDOF)
  
 END SUBROUTINE CopyFrom_EllipticFV1D 
!
!
!
 SUBROUTINE CopyTo_EllipticFV1D( this, sol )
 ! S/R CopyTo
 !
 !   Copies TO sol from the Elliptic_FV1D class
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS( EllipticFV1D ), INTENT(in) :: this
   REAL(prec), INTENT(out)           :: sol(1:this % nDOF)
   
    sol(1:this % nDOF) = this % localSol(1:this % nDOF)
  
 END SUBROUTINE CopyTo_EllipticFV1D 
!
!
!
 SUBROUTINE CalculateFlux_EllipticFV1D( this )
 !
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS( EllipticFV1D ), INTENT(inout) :: this
   ! LOCAL
   REAL(prec) :: sR, sL, dx, kappa
   INTEGER    :: i, nX
   
    dx = this % dx
    nX = this % nX
    kappa = this % diffCoeff
    
    ! A homogeneous Dirichlet boundary condition is invoked on the left-most boundary. To enforce
    ! this condition discretely in the finite-volume framework, the ghost-value is set to the negative
    ! of the nearest internal value.
    sL = -this % localSol(1)
    sR = this % localSol(1)
    ! From this prescription of the ghost value, the flux is computed and effectively enforces the
    ! Dirichlet condition at x=xL
    this % flux(1,1)  = -kappa*(sR-sL)/dx

    ! To apply the zero flux condition at x=xR, the last ghost value is set equal to the last 
    ! internal solution value.
    this % localSol(nX+1) = this % localSol(nX)
    
    DO i = 1, nX
       ! Calculate the flux on the right side of this element
       sL = this % localSol(i)
       sR = this % localSol(i+1)
       this % flux(i,2) = kappa*(sR-sL)/dx
       
       ! Global conservation requires that neighboring elements have fluxes that are equal in 
       ! magnitude, but opposite in sign. Because of this, we can set the left hand flux of the next
       ! element equal to the negative of the right hand flux of the current element
       this % flux(i+1,1) = -this % flux(i,2)
    ENDDO
       
       
  
 END SUBROUTINE CalculateFlux_EllipticFV1D
!
!
!
 FUNCTION MatrixAction_EllipticFV1D( this, s )
 !
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   CLASS( EllipticFV1D ) :: this
   REAL(prec)            :: s(1:this % nDOF)
   REAL(prec)            :: MatrixAction_EllipticFV1D(1:this % nDOF)
   ! LOCAL
   INTEGER    :: i, nX
   REAL(prec) :: dx
   
     dx = this % dx
     nX = this % nX
     
     ! In this function we want to return the matrix action applied to the vector "s".
     ! The Flux routine uses the "localSol" storage to calculate the flux. 
     ! Calculating the divergence of the flux gives the matrix-action. Thus, here the input
     ! variable "s" is copied to the localSol storage.
     CALL this % CopyFrom( s )
     ! Then the flux is calculated, enforcing homegeneous Dirichlet on the left boundary, x=xL,  
     ! and homogeneous Neumann on the right boundary, x=xR. See Flux_EllipticFV1D for details on
     ! the boundary condition implementation and flux calculation.
     CALL this % CalculateFlux( )     

     ! The flux divergence is calculated and assigned as the 
     DO i = 1, nX
        MatrixAction_EllipticFV1D(i) = (this % flux(i,2) + this % flux(i,1) )
     ENDDO

 END FUNCTION MatrixAction_EllipticFV1D
!
!
!
 FUNCTION Residual_EllipticFV1D( this, s )
 !
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   CLASS( EllipticFV1D ) :: this
   REAL(prec)            :: s(1:this % nDOF)
   REAL(prec)            :: Residual_EllipticFV1D(1:this % nDOF)
   ! LOCAL
   INTEGER    :: i, nX
   REAL(prec) :: dx
   REAL(prec) :: As(1:this % nDOF)
   
     dx = this % dx
     nX = this % nX
     
     As = this % MatrixAction( s )     

     ! The flux divergence is calculated and assigned as the 
     DO i = 1, nX
        Residual_EllipticFV1D(i) = this % source(i)*dx - As(i)
     ENDDO

 END FUNCTION Residual_EllipticFV1D
!
!
!==================================================================================================!
!----------------------------------------- Type Specific ------------------------------------------!
!==================================================================================================!
!
!
 SUBROUTINE WriteTecplot_EllipticFV1D( this )
 ! S/R WriteTecplot
 ! 
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS( EllipticFV1D ), INTENT(in) :: this
   ! LOCAL
   INTEGER :: i, N, fUnit

     N = this % nX

     OPEN( UNIT = NewUnit(fUnit), &
           FILE = 'FV1DSolution.curve', &
           FORM = 'FORMATTED' )

     WRITE( fUnit, * )'#ellipticSolution'

     DO i = 1, N

        WRITE(fUnit, '(F17.8,1x,F17.8)') this % x(i), this % localSol(i)     

     ENDDO

     CLOSE(fUnit)

     N = this % maxIters
  
     OPEN( UNIT = NewUnit(fUnit), &
           FILE = 'FV1D-CGResidual.curve', &
           FORM = 'FORMATTED' )

     WRITE( fUnit, * )'#CG-Residual'

     DO i = 0, N

        WRITE(fUnit, '(I6,1x,E17.8)') i, this % resi(i)     

     ENDDO

     WRITE( fUnit, * )'#CG-Residual-normalized'

     DO i = 0, N

        WRITE(fUnit, '(I6,1x,E17.8)') i, this % resi(i)/this % resi(0)     

     ENDDO

     CLOSE(fUnit)

 END SUBROUTINE WriteTecplot_EllipticFV1D

END MODULE EllipticFV1D_Class
