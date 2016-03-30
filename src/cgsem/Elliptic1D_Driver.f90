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
!     This is a test program for the "Elliptic1D" Class. The driver manages the construction, 
!     solving, file I/O, and destruction. The system solved here is
!
!        P_xx = 1               (1a)
!        P(0) = P(1) = 0        (1b)
!
!     The right-hand-side of (1a) can be changed in the subroutine "FillSource_Elliptic1D" in the
!     module "Elliptic1D_Class.f90". The values used for the Dirichlet boundary condition can be
!     set in the subroutine "DirichletConditions_Elliptic1D". Additionally, Neumman boundary 
!     conditions can be set by changing value of the attribute "dirichletMask" and by setting
!     the "boundaryFlux" attribute to the desired value.
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
