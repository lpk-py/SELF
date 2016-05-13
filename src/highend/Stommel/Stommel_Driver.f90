! Stommel_Driver.f90 ( new with v2.1 - 29 March 2016)
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

!
!   
! 
! //////////////////////////////////////////////////////////////////////////////////////////////// !

PROGRAM Stommel_Driver

USE ModelPrecision

USE QuadMeshClass

USE StommelParams_Class
USE Stommel_Class


IMPLICIT NONE

   TYPE( Stommel )       :: ellipticSolver
   TYPE( StommelParams ) :: params
   INTEGER :: ioerr

   CALL params % Build( )

   CALL ellipticSolver % Build( params )

   CALL SetCoriolisParameter( ellipticSolver )
   CALL SetBathymetry( ellipticSolver )

   ! Now that the bathymetry and the coriolis parameter are set, the potential vorticity gradient
   ! can be calculated and stored in the "Stommel" class. Recall that the potential vorticity in 
   ! the Stommel model is f/h, where "f" is the coriolis parameter and "h" is the fluid depth.
   CALL CalculatePVGradient( ellipticSolver )

   ! Now that the bathymetry and the coriolis parameter are set, the flux coefficent for the
   ! drag operator can be reset. * Recall that the flux coefficient is cDrag/h^2 where "cDrag"
   ! is the linear drag coefficient and "h" is the fluid depth.
   CALL ellipticSolver % ResetFluxCoefficient( )

   CALL SetEkmanPumping( ellipticSolver, params )

   CALL ellipticSolver % Solve( DirichletBC, ioerr )


   CALL ellipticSolver % WriteTecplot( )
   CALL ellipticSolver % WriteResidual( )

   CALL ellipticSolver % Trash( )

CONTAINS

 SUBROUTINE SetCoriolisParameter( myCGSEM, params )
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   TYPE( Stommel ), INTENT(inout)    :: myCGSEM
   TYPE( StommelParams ), INTENT(in) :: params
   ! LOCAL
   INTEGER    :: i, j, k, N, nEl
   REAL(prec) :: x, y, f0, beta

      
      N   = myCGSEM % nS
      nEl = myCGSEM % mesh % nElems
      f0  = params % f0
      beta = params % betaY

      DO k = 1, nEl
         DO j = 0, N
            DO i = 0, N

               CALL myCGSEM % mesh % GetPositionAtNode( k, x, y, i, j ) ! x and y are the physical locations in the mesh
               myCGSEM % fCori(i,j,k) = f0 + beta*y

            ENDDO
         ENDDO
      ENDDO

 END SUBROUTINE SetCoriolisParameter
!
!
!
  SUBROUTINE SetBathymetry( myCGSEM, params )
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   TYPE( Stommel ), INTENT(inout)    :: myCGSEM
   TYPE( StommelParams ), INTENT(in) :: params
   ! LOCAL
   INTEGER    :: i, j, k, N, nEl
   REAL(prec) :: x, y

      
      N   = myCGSEM % nS
      nEl = myCGSEM % mesh % nElems


      DO k = 1, nEl
         DO j = 0, N
            DO i = 0, N

               CALL myCGSEM % mesh % GetPositionAtNode( k, x, y, i, j ) ! x and y are the physical locations in the mesh
               myCGSEM % h(i,j,k) = 1000.0_prec

            ENDDO
         ENDDO
      ENDDO

 END SUBROUTINE SetBathymetry
!
!
!
 SUBROUTINE SetEkmanPumping( myCGSEM, params )
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   TYPE( Stommel ), INTENT(inout)    :: myCGSEM
   TYPE( StommelParams ), INTENT(in) :: params
   ! LOCAL
   INTEGER    :: i, j, k, N, nEl
   REAL(prec) :: x, y, f, Ly

      
      N   = myCGSEM % nS
      nEl = myCGSEM % mesh % nElems
      Ly  = params % Ly

      DO k = 1, nEl
         DO j = 0, N
            DO i = 0, N

               CALL myCGSEM % mesh % GetPositionAtNode( k, x, y, i, j ) ! x and y are the physical locations in the mesh
               f = myCGSEM % fCori(i,j,k)

               myCGSEM % source(i,j,k) = -f*0.001_prec*( sin(pi*y/Ly)  )

            ENDDO
         ENDDO
      ENDDO

 END SUBROUTINE SetEkmanPumping
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
!
!
!
  SUBROUTINE CalculatePVGradient( myCGSEM, params )
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   TYPE( Stommel ), INTENT(inout)    :: myCGSEM
   TYPE( StommelParams ), INTENT(in) :: params
   ! LOCAL
   INTEGER    :: i, j, k, N, nEl
   REAL(prec) :: q(0:myCGSEM % nS, 0:myCGSEM % nP)
   REAL(prec) :: qx(0:myCGSEM % nS, 0:myCGSEM % nP)
   REAL(prec) :: qy(0:myCGSEM % nS, 0:myCGSEM % nP)

      N   = myCGSEM % nS
      nEl = myCGSEM % mesh % nElems

      DO k = 1, nEl

         q = myCGSEM % f(:,:,k)/myCGSEM % h(:,:,k)

         CALL myCGSEM % CalculateGradient( k, q, qx, qy)

         myCGSEM % qX(:,:,k) = qx
         myCGSEM % qY(:,:,k) = qy

      ENDDO


 END SUBROUTINE CalculatePVGradient
!
!
!

END PROGRAM Stommel_Driver
