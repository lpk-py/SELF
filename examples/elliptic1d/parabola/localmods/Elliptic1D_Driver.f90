! Elliptic1D_Driver.f90 ( new with v2.1 - 29 March 2016)
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

PROGRAM Elliptic1D_Driver

USE ModelPrecision

USE Elliptic1D_Class


IMPLICIT NONE

   TYPE( Elliptic1D ) :: ellipticSolver
   INTEGER :: ioerr


   CALL ellipticSolver % Build( )

   CALL ellipticSolver % Solve( ioerr )

   CALL ellipticSolver % WriteTecplot( )
   CALL ellipticSolver % WriteResidual( )

   CALL ellipticSolver % Trash( )

END PROGRAM Elliptic1D_Driver
