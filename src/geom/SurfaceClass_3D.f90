! SurfaceClass_3D.f90
! 
! Copyright 2015 Joseph Schoonover <schoonover.numerics@gmail.com>
! 
! SurfaceClass_3D.f90 is part of the Spectral Element Libraries in Fortran (SELF).
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
 
 

MODULE SurfaceClass_3D
! ========================================= Logs ================================================= !
!2016-05-11  Joseph Schoonover  schoonover.numerics@gmail.com 
!
! //////////////////////////////////////////////////////////////////////////////////////////////// ! 
! src/common/
USE ModelPrecision
USE ConstantsDictionary
! src/interp/
USE Lagrange_1D_Class
USE Lagrange_2D_Class


   TYPE Surface_3D 
      INTEGER, PRIVATE                  :: nS, nP
      TYPE(Lagrange_2D), PRIVATE        :: interp
      REAL(prec),  ALLOCATABLE, PRIVATE :: x(:,:), y(:,:), z(:,:)
      REAL(prec),  ALLOCATABLE, PRIVATE :: dxds(:,:), dyds(:,:), dzds(:,:)
      REAL(prec),  ALLOCATABLE, PRIVATE :: dxdp(:,:), dydp(:,:), dzdp(:,:)
      
      CONTAINS

      PROCEDURE :: Build => Build_Surface_3D
      PROCEDURE :: Trash => Trash_Surface_3D

      PROCEDURE :: SetNumberOfNodes => SetNumberOfNodes_Surface_3D
      PROCEDURE :: GetNumberOfNodes => GetNumberOfNodes_Surface_3D
      PROCEDURE :: SetInterpNodes => SetInterpNodes_Surface_3D
      PROCEDURE :: GetInterpNodes => GetInterpNodes_Surface_3D
      PROCEDURE :: SetNodes => SetNodes_Surface_3D
      PROCEDURE :: GetNodes => GetNodes_Surface_3D
      PROCEDURE :: GetPositionAtNode => GetPositionAtNode_Surface_3D
      PROCEDURE :: SetSlope => SetSlope_Surface_3D
      PROCEDURE :: GetSlope => GetSlope_Surface_3D

      PROCEDURE :: CalculateMetrics => CalculateMetrics_Surface_3D
      PROCEDURE :: Evaluate => Evaluate_Surface_3D
      PROCEDURE :: EvaluateMetrics => EvaluateMetrics_Surface_3D
      
      
   END TYPE Surface_3D


 INTEGER, PARAMETER, PRIVATE    :: defaultN = 1 ! default number of nodes
 REAL(prec), PARAMETER, PRIVATE :: defaultXY = ZERO
 
 
 CONTAINS
!
!
!==================================================================================================!
!------------------------------- Manual Constructors/Destructors ----------------------------------!
!==================================================================================================!
!
!
 SUBROUTINE Build_Surface_3D( mySurface, x, y, z, s, p )
 ! S/R Build
 !
 ! =============================================================================================== !
 ! DECLARATIONS
  IMPLICIT NONE
  CLASS( Surface_3D ), INTENT(out) :: mySurface
  REAL(prec), INTENT(in), OPTIONAL :: x(0:,0:), y(0:,0:), z(0:,0:), s(0:), p(0:)
  ! LOCAL
  INTEGER :: nS, nP
  
      IF( PRESENT(x) .AND. PRESENT(y) .AND. PRESENT(z) .AND. PRESENT(s) .AND. PRESENT(p) )THEN
         nS = UBOUND(s, DIM=1 )
         nP = UBOUND(p, DIM=1 )
      ELSE
         nS = defaultN
         nP = nS
      ENDIF
      
      CALL mySurface % SetNumberOfNodes( nS, nP )
      
      IF( PRESENT(s) .AND. PRESENT(p) )THEN
         CALL mySurface % interp % Build( nS, nP, s, p )
      ELSE
         CALL mySurface % interp % Build( nS, nP )
      ENDIF

      ALLOCATE( mySurface % x(0:nS,0:nP), mySurface % y(0:nS,0:nP), mySurface % z(0:nS,0:nP) )
      ALLOCATE( mySurface % dxds(0:nS,0:nP), mySurface % dyds(0:nS,0:nP), mySurface % dzds(0:nS,0:nP) )
      ALLOCATE( mySurface % dxdp(0:nS,0:nP), mySurface % dydp(0:nS,0:nP), mySurface % dzdp(0:nS,0:nP) )

      IF( PRESENT(x) .AND. PRESENT(y) .AND. PRESENT(z) )THEN
         CALL mySurface % SetNodes( x, y, z )
      ELSE
         CALL mySurface % SetNodes( )
      ENDIF
      
      CALL mySurface % CalculateMetrics( )
 
 END SUBROUTINE Build_Surface_3D
!
!
!
 SUBROUTINE Trash_Surface_3D( mySurface )
 ! S/R Trash
 !
 ! =============================================================================================== !
 ! DECLARATIONS
  IMPLICIT NONE
  CLASS( Surface_3D ), INTENT(inout)     :: mySurface
  
      CALL mySurface % interp % Trash( )

      DEALLOCATE( mySurface % x, mySurface % y, mySurface % z )

      DEALLOCATE( mySurface % dxds, mySurface % dyds, mySurface % dzds )

 
 END SUBROUTINE Trash_Surface_3D
!
!
!==================================================================================================!
!--------------------------------------- Accessors ------------------------------------------------!
!==================================================================================================!
!
!
 SUBROUTINE SetNumberOfNodes_Surface_3D( mySurface, nS, nP )
 ! S/R SetNumberOfNodes
 !
 ! =============================================================================================== !
 ! DECLARATIONS
  IMPLICIT NONE
  CLASS( Surface_3D ), INTENT(inout) :: mySurface
  INTEGER, INTENT(in)                :: nS, nP
  
      mySurface % nS = nS
      mySurface % nP = nP

 END SUBROUTINE SetNumberOfNodes_Surface_3D
!
!
!
 SUBROUTINE GetNumberOfNodes_Surface_3D( mySurface, nS, nP )
 ! S/R GetNumberOfNodes
 !
 ! =============================================================================================== !
 ! DECLARATIONS
  IMPLICIT NONE
  CLASS( Surface_3D ), INTENT(in) :: mySurface
  INTEGER, INTENT(out)            :: nS, nP
  
      nS = mySurface % nS
      nP = mySurface % nP

 END SUBROUTINE GetNumberOfNodes_Surface_3D
!
!
!
 SUBROUTINE SetInterpNodes_Surface_3D( mySurface, s, p )
 ! S/R SetInterpNodes
 !
 ! =============================================================================================== !
 ! DECLARATIONS
  IMPLICIT NONE
  CLASS( Surface_3D ), INTENT(inout) :: mySurface
  REAL(prec), INTENT(in)             :: s(0:mySurface % nS)
  REAL(prec), INTENT(in)             :: p(0:mySurface % nP)
  
      CALL mySurface % interp % SetNodes( s, p )

 END SUBROUTINE SetInterpNodes_Surface_3D
!
!
!
 SUBROUTINE GetInterpNodes_Surface_3D( mySurface, s, p )
 ! S/R GetInterpNodes
 !
 ! =============================================================================================== !
 ! DECLARATIONS
  IMPLICIT NONE
  CLASS( Surface_3D ), INTENT(in) :: mySurface
  REAL(prec), INTENT(out)         :: s(0:mySurface % nS)
  REAL(prec), INTENT(out)         :: p(0:mySurface % nP)
  
      CALL mySurface % interp % GetNodes( s, p )

 END SUBROUTINE GetInterpNodes_Surface_3D
!
!
!
 SUBROUTINE SetNodes_Surface_3D( mySurface, x, y, z )
 ! S/R SetNodes
 !
 ! =============================================================================================== !
 ! DECLARATIONS
  IMPLICIT NONE
  CLASS( Surface_3D ), INTENT(inout) :: mySurface
  REAL(prec), INTENT(in), OPTIONAL   :: x(0:mySurface % nS, 0:mySurface % nP)
  REAL(prec), INTENT(in), OPTIONAL   :: y(0:mySurface % nS, 0:mySurface % nP)
  REAL(prec), INTENT(in), OPTIONAL   :: z(0:mySurface % nS, 0:mySurface % nP)
  
      IF( PRESENT(x) .AND. PRESENT(y) .AND. PRESENT(z) )THEN
         mySurface % x = x
         mySurface % y = y
         mySurface % z = z
      ELSE
         mySurface % x = defaultXY
         mySurface % y = defaultXY
         mySurface % z = defaultXY
      ENDIF

 END SUBROUTINE SetNodes_Surface_3D
!
!
!
 SUBROUTINE GetNodes_Surface_3D( mySurface, x, y, z )
 ! S/R GetNodes
 !
 ! =============================================================================================== !
 ! DECLARATIONS
  IMPLICIT NONE
  CLASS( Surface_3D ), INTENT(in) :: mySurface
  REAL(prec), INTENT(out)         :: x(0:mySurface % nS, 0:mySurface % nP)
  REAL(prec), INTENT(out)         :: y(0:mySurface % nS, 0:mySurface % nP)
  REAL(prec), INTENT(out)         :: z(0:mySurface % nS, 0:mySurface % nP)
  
      x = mySurface % x
      y = mySurface % y 
      z = mySurface % z

 END SUBROUTINE GetNodes_Surface_3D
!
!
!
 SUBROUTINE GetPositionAtNode_Surface_3D( mySurface, x, y, z, i, j )
 ! S/R GetPositionAtNode
 !
 ! =============================================================================================== !
 ! DECLARATIONS
  IMPLICIT NONE
  CLASS( Surface_3D ), INTENT(in) :: mySurface
  REAL(prec), INTENT(out)         :: x, y ,z
  INTEGER, INTENT(in)             :: i, j
  
      x = mySurface % x(i,j)
      y = mySurface % y(i,j)
      z = mySurface % z(i,j)

 END SUBROUTINE GetPositionAtNode_Surface_3D
!
!
!
 SUBROUTINE SetSlope_Surface_3D( mySurface, dxds, dyds, dzds )
 ! S/R SetSlope
 !
 ! =============================================================================================== !
 ! DECLARATIONS
  IMPLICIT NONE
  CLASS( Surface_3D ), INTENT(inout) :: mySurface
  REAL(prec), INTENT(in)             :: dxds(0:mySurface % nS, 0:mySurface % nP)
  REAL(prec), INTENT(in)             :: dyds(0:mySurface % nS, 0:mySurface % nP)
  REAL(prec), INTENT(in)             :: dzds(0:mySurface % nS, 0:mySurface % nP)
  
      mySurface % dxds = dxds
      mySurface % dyds = dyds
      mySurface % dzds = dzds

 END SUBROUTINE SetSlope_Surface_3D
!
!
!
 SUBROUTINE GetSlope_Surface_3D( mySurface, dxds, dyds, dzds )
 ! S/R GetSlope
 !
 ! =============================================================================================== !
 ! DECLARATIONS
  IMPLICIT NONE
  CLASS( Surface_3D ), INTENT(in) :: mySurface
  REAL(prec), INTENT(out)         :: dxds(0:mySurface % nS, 0:mySurface % nP)
  REAL(prec), INTENT(out)         :: dyds(0:mySurface % nS, 0:mySurface % nP)
  REAL(prec), INTENT(out)         :: dzds(0:mySurface % nS, 0:mySurface % nP)
  
      dxds = mySurface % dxds
      dyds = mySurface % dyds 
      dzds = mySurface % dzds

 END SUBROUTINE GetSlope_Surface_3D
!
!
!==================================================================================================!
!--------------------------------- Type Specific Routines -----------------------------------------!
!==================================================================================================!
!
!
 SUBROUTINE CalculateMetrics_Surface_3D( mySurface )
 ! S/R CalculateMetrics
 !
 ! =============================================================================================== !
 ! DECLARATIONS
  IMPLICIT NONE
  CLASS( Surface_3D ), INTENT(inout) :: mySurface
  ! LOCAL
  REAL(prec), ALLOCATABLE :: dMatS(:,:), dMatP(:,:)
  INTEGER :: iS, iP, nS, nP
  
     CALL mySurface % GetNumberOfNodes( nS, nP )
  
     CALL mySurface % interp % CalculateDerivativeMatrix( dMatS, dMatP )
     
     ! Calculate the s-derivatives
     DO iP = 0, nP           
        mySurface % dxds(0:nS,iP) = MATMUL( dMatS, mySurface % x(0:nS,iP) )
        mySurface % dyds(0:nS,iP) = MATMUL( dMatS, mySurface % y(0:nS,iP) )
        mySurface % dzds(0:nS,iP) = MATMUL( dMatS, mySurface % z(0:nS,iP) )
     ENDDO
     
     ! Calculate the p-derivatives
     DO iS = 0, nS           
        mySurface % dxdp(iS,0:nP) = MATMUL( dMatP, mySurface % x(iS,0:nP) )
        mySurface % dydp(iS,0:nP) = MATMUL( dMatP, mySurface % y(iS,0:nP) )
        mySurface % dzdp(iS,0:nP) = MATMUL( dMatP, mySurface % z(iS,0:nP) )
     ENDDO

     DEALLOCATE( dMatS, dMatP )
     
 END SUBROUTINE CalculateMetrics_Surface_3D
!
!
!
 SUBROUTINE Evaluate_Surface_3D( mySurface, s, p, x, y, z )
 ! S/R Evaluate
 !
 ! =============================================================================================== !
 ! DECLARATIONS
  IMPLICIT NONE
  CLASS( Surface_3D ), INTENT(in) :: mySurface
  REAL(prec), INTENT(in)          :: s, p 
  REAL(prec), INTENT(out)         :: x, y, z
  
      x = mySurface % interp % EvaluateInterpolant( s, p, mySurface % x )
      y = mySurface % interp % EvaluateInterpolant( s, p, mySurface % y )
      z = mySurface % interp % EvaluateInterpolant( s, p, mySurface % z )
      
 END SUBROUTINE Evaluate_Surface_3D
!
!
!
 SUBROUTINE EvaluateMetrics_Surface_3D( mySurface, s, p, dxds, dyds, dzds, dxdp, dydp, dzdp )
 ! S/R EvaluateMetrics
 !
 ! =============================================================================================== !
 ! DECLARATIONS
  IMPLICIT NONE
  CLASS( Surface_3D ), INTENT(in) :: mySurface
  REAL(prec), INTENT(in)          :: s, p 
  REAL(prec), INTENT(out)         :: dxds, dyds, dzds, dxdp, dydp, dzdp
  
      dxds = mySurface % interp % EvaluateInterpolant( s, p, mySurface % dxds )
      dyds = mySurface % interp % EvaluateInterpolant( s, p, mySurface % dyds )
      dzds = mySurface % interp % EvaluateInterpolant( s, p, mySurface % dzds )
      
      dxdp = mySurface % interp % EvaluateInterpolant( s, p, mySurface % dxdp )
      dydp = mySurface % interp % EvaluateInterpolant( s, p, mySurface % dydp )
      dzdp = mySurface % interp % EvaluateInterpolant( s, p, mySurface % dzdp )

 END SUBROUTINE EvaluateMetrics_Surface_3D
!
!
!
END MODULE SurfaceClass_3D
