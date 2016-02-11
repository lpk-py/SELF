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


   TYPE, EXTENDS( ConjugateGradient ) :: EllipticFD1D 

   END TYPE EllipticFD1D

END MODULE EllipticFD1D_Class
