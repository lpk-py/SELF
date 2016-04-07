! CGsemElliptic2D_Driver.f90 ( new with v2.1 - 29 March 2016)
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
! o  (ver 2.1) March 2016
!
! ========================================= Logs ================================================= !
! yyyy-mm-dd  Joe  <joe@clay>
!    This module defines the data structure used for implementing the Continuous Galerkin
!    Spectral Element Method. 
!
!    The module  is set up so that we solve d/dx( F ) = Q
!
!    For the heat equation : F = k du/dx, Q = du/dt
!
!    With this "generic flux" (F), and source (Q), it should be a trivial task to solve a different
!    system with CGSEM through the modification of F and Q. Here, the module is set up to solve
!    the heat-equation using the trapezoidal rule for time integration.
!
!   
! 
! //////////////////////////////////////////////////////////////////////////////////////////////// !

PROGRAM CGsemElliptic2D_Driver

USE ModelPrecision

USE QuadMeshClass

USE CGsemParams_Class
USE CGsemElliptic_2D_Class


IMPLICIT NONE

   TYPE( CGsemElliptic_2D ) :: ellipticSolver
   TYPE( CGsemParams )      :: params
   INTEGER :: ioerr

   CALL params % Build( )

   CALL ellipticSolver % Build( params )
   CALL FillInSource( ellipticSolver )
   CALL ellipticSolver % Solve( DirichletBC, ioerr )

   CALL ellipticSolver % ConstructMatrix( )
   CALL ellipticSolver % preconditioner % ConstructPCMatrix( )

   CALL ellipticSolver % WriteTecplot( )
   CALL ellipticSolver % WriteResidual( )

   CALL ellipticSolver % Trash( )

CONTAINS

 SUBROUTINE FillInSource( myCGSEM )
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS( CGsemElliptic_2D ), INTENT(inout) :: myCGSEM
   ! LOCAL
   INTEGER    :: i, j, k, N, nEl
   REAL(prec) :: x, y

      
      N   = myCGSEM % nS
      nEl = myCGSEM % mesh % nElems
 
      DO k = 1, nEl
         DO j = 0, N
            DO i = 0, N

               CALL myCGSEM % mesh % GetPositionAtNode( k, x, y, i, j )
               myCGSEM % source(i,j,k) = exp( -( (x-HALF)**2 + (y-HALF)**2 ) ) !x*y

            ENDDO
         ENDDO
      ENDDO

 END SUBROUTINE FillInSource
!
!
! 
 FUNCTION DirichletBC( x, y )
 !
 ! =============================================================================================== !
   IMPLICIT NONE  
   REAL(prec) :: x, y
   REAL(prec) :: DirichletBC 

      DirichletBC = ZERO

 END FUNCTION DirichletBC

END PROGRAM CGsemElliptic2D_Driver
