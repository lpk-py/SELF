! Driver_EllipticFV1D.f90 ( new with v2.1 - 11 Feb. 2016)
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
! 2016-02-17  Joe  <joe@clay>
!
!   This program tests the EllipticFV1D_Class module. The module is set-up to solve a problem of the
!   form 
!              u_xx = s
!             u(xL) = 0
!           u_x(xR) = 0
! 
! //////////////////////////////////////////////////////////////////////////////////////////////// !

PROGRAM Driver_EllipticFV1D

USE ModelPrecision
USE ConstantsDictionary
USE Timing 
!
USE IterativeSolvers_Class
USE EllipticFV1D_Class
USE IterativeSolversParams_Class

 IMPLICIT NONE
 

 TYPE( IterativeSolversParams ) :: params
 TYPE( EllipticFV1D )           :: solver
 TYPE( MultiTimers )            :: clocks
 
 INTEGER    :: N, i, ioerr
 REAL(prec) :: x 
 
   CALL params % Build( )

   CALL solver % Initialize( params )

   CALL clocks % Build( )
   CALL clocks % AddTimer( 'Conjugate Gradient Method', 1 )

   ! The RHS of the elliptic equation (the source term) is set up here in addition to the initial 
   ! condition.
   N = solver % nX
   
   DO i = 1, N
      x = solver % x(i)
      solver % source(i) = TheSource( x )
   ENDDO

   ! To get a good measure of the timing, we repeat the solve a specified number of times (params % nTimingCycles)
   DO i = 1, params % nTimingCycles

      solver % localSol = ZERO ! Reset the initial condition to repeat the experiment


      CALL clocks % StartThisTimer( 1 )
      CALL solver % Solve( ioerr )
      CALL clocks % StopThisTimer( 1 )

   ENDDO

   CALL clocks % Write_MultiTimers( )

   ! Should add an "exception handler" - from ioerr, trigger a set of events to occur

   CALL solver % WriteTecplot( )

   CALL solver % Finallize( )

   CALL clocks % Trash( )

CONTAINS

 FUNCTION TheSource( x ) RESULT( s )
   REAL(prec) :: x, s

      s = ONE

 END FUNCTION TheSource

END PROGRAM Driver_EllipticFV1D
