! SolutionStorage_Class.f90 ( new with v2.1 - 7 Feb. 2016)
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
!   This module provides an abstract type for handling the storage of solution vectors that can be
!   used in iterative solvers. Residuals, search directions, and other intermediate states in an 
!   iterative solver should be stored in a manner consistent with the problem's natural storage 
!   system.
!   must be a class that has routines called 
!      (1) MatrixAction -- compute Ax
!      (2) Residual -- computes r = b-Ax
!      (3) VectorProduct -- computes the equivalent of a dot-product for two "solution-vectors" (x)
!      (4) UpdateSearchDirection -- takes in an integer (to reference which search direction)
!      (5) SolutionUpdate -- Given a new search direction in the solution, the solution should be 
!                            able to be updated with this information. Because of this, the data-
!                            structure must allow for storage of the search direction(s) (same type as
!                            the solution)
!      (6) PreconditionSolve -- Optionally, a preconditioner can be provided to solve Hz = r, where
!                               H is the preconditioning matrix
!
! To accomplish this, an abstract data-type is provided here with the first four of the above listed
! procedures as "Deferred" procedures --> When the end-user wishes to use this module, these deferred
! procedures must be defined in the extended data type.
!
! The extended data-structure must make use of a "SolutionStorage" type
! 
! Further, the class must have the 
! 
! //////////////////////////////////////////////////////////////////////////////////////////////// !

MODULE IterativeSolvers_Class

USE ModelPrecision

USE CGSEM2D_Class


    TYPE, ABSTRACT SolutionStorage
     
     PROCEDURE (SolDotSol) :: VectorProduct
  END TYPE SolutionStorage


  ABSTRACT INTERFACE
     FUNCTION SolDotSol( this, that )
        IMPORT SolutionStorage
        CLASS(SolutionStorage) :: this, that
        REAL(prec)             :: SolDotSol
     END FUNCTION StateDotState
  END INTERFACE




  
  
END MODULE SolutionStorage_Class
