! LBSW_InitialConditions.f90
! 
! Copyright 2015 Joseph Schoonover <schoonover.numerics@gmail.com>
! 
! LBSW_InitialConditions.f90 is part of the Spectral Element Libraries in Fortran (SELF).
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
 

PROGRAM LBSW_InitialConditions
! ========================================= Logs ================================================= !
!2016-05-19  Joseph Schoonover  schoonover.numerics@gmail.com 
!
! //////////////////////////////////////////////////////////////////////////////////////////////// ! 

USE ModelPrecision

USE QuadMeshClass

USE LBSW_Class


IMPLICIT NONE

   TYPE( LBSW )         :: shelfwaves
   INTEGER :: ioerr

   CALL params % Build( )

   CALL shelfwaves % Build( params )

   CALL SetBackgroundVelocity( shelfwaves )
   CALL SetBathymetry( shelfwaves )
   CALL shelfwaves % vortInverter % ResetFluxCoefficient( )

   CALL SetInitialVorticity( shelfwaves )
   CALL shelfwaves % vortInverter % Solve( ioerr )
   CALL shelwaves % CalculateVelocity( )

   ! Now that the bathymetry and the background meridional velocity are set, the "PVFactor" can
   ! be set : PVFactor = (background vorticity)*(h_x/h)
   CALL CalculatePVFactor( shelfwaves )

   ! Additionally, the "forcing" can be set where forcing = v*h'_y/h, where h' is the 
   ! shelf steepening component of the bathymetry
   CALL SetForcing( shelfwaves )


   CALL shelfwaves % WriteTecplot( 'LBSW.init' )
   CALL shelfwaves % vortInverter % WriteTecplot( 'Stream.init' )
   CALL shelfwaves % WriteResidual( )

   CALL shelfwaves % Trash( )

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

END PROGRAM LBSW_InitialConditions
