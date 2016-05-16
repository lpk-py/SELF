! Stommel_Driver.f90
! 
! Copyright 2015 Joseph Schoonover <schoonover.numerics@gmail.com>
! 
! Stommel_Driver.f90 is part of the Spectral Element Libraries in Fortran (SELF).
! 
! Permission is hereby granted, free of charge, to any person obtaining a copy of this software 
! and associated documentation files (the "Software"), to deal in the Software without restriction, 
! including without limitation the rights to use, copy, modify, merge, publish, distribute, 
! sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is 
! furnished to do so, subject to the following conditions: 
! 
! The above copyright notice and this permission notice shall be included in all copies or  
! substantial portions of the Software. 
! 
! THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING 
! BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND 
! NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, 
! DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, 
! OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE. 
!
! //////////////////////////////////////////////////////////////////////////////////////////////// !
  
PROGRAM Stommel_Driver
! ========================================= Logs ================================================= !
!2016-05-13  Joseph Schoonover  schoonover.numerics@gmail.com 
!
! //////////////////////////////////////////////////////////////////////////////////////////////// ! 

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

   CALL SetCoriolisParameter( ellipticSolver, params )
   CALL SetBathymetry( ellipticSolver, params )

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


   CALL ellipticSolver % mesh % WriteTecplot( )
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
               myCGSEM % h(i,j,k) = ONE + 0.1_prec*exp( -( (x-5.0_prec)**2 + (y-5.0_prec)**2 ) )

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
   REAL(prec) :: x, y, f, Ly, A

      
      N   = myCGSEM % nS
      nEl = myCGSEM % mesh % nElems
      Ly  = params % yScale
      A   = ONE
      DO k = 1, nEl
         DO j = 0, N
            DO i = 0, N

               CALL myCGSEM % mesh % GetPositionAtNode( k, x, y, i, j ) ! x and y are the physical locations in the mesh

               myCGSEM % source(i,j,k) = -A*( sin(pi*y/Ly)  )

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
  SUBROUTINE CalculatePVGradient( myCGSEM )
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   TYPE( Stommel ), INTENT(inout)    :: myCGSEM
   ! LOCAL
   INTEGER    :: k, nEl

      nEl = myCGSEM % mesh % nElems

      DO k = 1, nEl

         myCGSEM % Q(:,:,k) = myCGSEM % fCori(:,:,k)/myCGSEM % h(:,:,k)
 

      ENDDO


 END SUBROUTINE CalculatePVGradient
!
!
!

END PROGRAM Stommel_Driver
