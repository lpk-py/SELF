!  Elliptic1D_Driver.f90
!  
!  Copyright 2015 Joseph Schoonover <schoonover.numerics@gmail.com>
!  
!  Elliptic1D_Driver.f90 is part of the Spectral Element Libraries in Fortran (SELF).
!
!  Permission is hereby granted, free of charge, to any person obtaining a copy of this software
!  and associated documentation files (the "Software"), to deal in the Software without restriction,
!  including without limitation the rights to use, copy, modify, merge, publish, distribute, 
!  sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is 
!  furnished to do so, subject to the following conditions:
!
!  The above copyright notice and this permission notice shall be included in all copies or 
!  substantial portions of the Software.
!
!  THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING
!  BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
!  NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM,
!  DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, 
!  OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
!  
! //////////////////////////////////////////////////////////////////////////////////////////////// !
!


PROGRAM Elliptic1D_Driver
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
