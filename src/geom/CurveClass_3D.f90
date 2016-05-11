! CurveClass_3D.f90
! 
! Copyright 2015 Joseph Schoonover <schoonover.numerics@gmail.com>
! 
! CurveClass_3D.f90 is part of the Spectral Element Libraries in Fortran (SELF).
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
 
 
MODULE CurveClass_3D
! ========================================= Logs ================================================= !
!2016-05-11  Joseph Schoonover  schoonover.numerics@gmail.com 
!
! //////////////////////////////////////////////////////////////////////////////////////////////// ! 

! src/common/
USE ModelPrecision
USE ConstantsDictionary
! src/interp/
USE Lagrange_1D_Class


   TYPE Curve_3D 
      INTEGER, PRIVATE                  :: nNodes
      TYPE(Lagrange_1D), PRIVATE        :: interp
      REAL(prec),  ALLOCATABLE, PRIVATE :: x(:), y(:), z(:)
      REAL(prec),  ALLOCATABLE, PRIVATE :: dxds(:), dyds(:), dzds(:)

      CONTAINS

      PROCEDURE :: Build => Build_Curve_3D
      PROCEDURE :: Trash => Trash_Curve_3D

      PROCEDURE :: SetNumberOfNodes => SetNumberOfNodes_Curve_3D
      PROCEDURE :: GetNumberOfNodes => GetNumberOfNodes_Curve_3D
      PROCEDURE :: SetInterpNodes => SetInterpNodes_Curve_3D
      PROCEDURE :: GetInterpNodes => GetInterpNodes_Curve_3D
      PROCEDURE :: SetNodes => SetNodes_Curve_3D
      PROCEDURE :: GetNodes => GetNodes_Curve_3D
      PROCEDURE :: SetSlope => SetSlope_Curve_3D
      PROCEDURE :: GetSlope => GetSlope_Curve_3D

      PROCEDURE :: CalculateSlope => CalculateSlope_Curve_3D
      PROCEDURE :: Evaluate => Evaluate_Curve_3D
      PROCEDURE :: EvaluateSlope => EvaluateSlope_Curve_3D
      
      
   END TYPE Curve_3D


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
 SUBROUTINE Build_Curve_3D( myCurve_3D, x, y, z, nodes )
 ! S/R Build
 !
 ! =============================================================================================== !
 ! DECLARATIONS
  IMPLICIT NONE
  CLASS( Curve_3D ), INTENT(out)      :: myCurve_3D
  REAL(prec), INTENT(in), OPTIONAL :: x(0:), y(0:), z(0:), nodes(0:)
  ! LOCAL
  INTEGER :: N
  
      IF( PRESENT(x) .AND. PRESENT(y) .AND. PRESENT(z) .AND. PRESENT(nodes) )THEN
         N = UBOUND(x, DIM=1 )
      ELSE
         N = defaultN
      ENDIF
      
      CALL myCurve_3D % SetNumberOfNodes( N )
      
      IF( PRESENT(nodes) )THEN
         CALL myCurve_3D % interp % Build( N, nodes )
      ELSE
         CALL myCurve_3D % interp % Build( N )
      ENDIF

      ALLOCATE( myCurve_3D % x(0:N), myCurve_3D % y(0:N), myCurve_3D % z(0:N) )

      ALLOCATE( myCurve_3D % dxds(0:N), myCurve_3D % dyds(0:N), myCurve_3D % dzds(0:N) )

      IF( PRESENT(x) .AND. PRESENT(y) )THEN
         CALL myCurve_3D % SetNodes( x, y, z )
      ELSE
         CALL myCurve_3D % SetNodes( )
      ENDIF
      
      CALL myCurve_3D % CalculateSlope( )
 
 END SUBROUTINE Build_Curve_3D
!
!
!
 SUBROUTINE Trash_Curve_3D( myCurve_3D )
 ! S/R Trash
 !
 ! =============================================================================================== !
 ! DECLARATIONS
  IMPLICIT NONE
  CLASS( Curve_3D ), INTENT(inout)     :: myCurve_3D
  
      CALL myCurve_3D % interp % Trash( )

      DEALLOCATE( myCurve_3D % x, myCurve_3D % y, myCurve_3D % z )

      DEALLOCATE( myCurve_3D % dxds, myCurve_3D % dyds, myCurve_3D % dzds )

 
 END SUBROUTINE Trash_Curve_3D
!
!
!==================================================================================================!
!--------------------------------------- Accessors ------------------------------------------------!
!==================================================================================================!
!
!
 SUBROUTINE SetNumberOfNodes_Curve_3D( myCurve_3D, N )
 ! S/R SetNumberOfNodes
 !
 ! =============================================================================================== !
 ! DECLARATIONS
  IMPLICIT NONE
  CLASS( Curve_3D ), INTENT(inout) :: myCurve_3D
  INTEGER, INTENT(in)           :: N
  
      myCurve_3D % nNodes = N

 END SUBROUTINE SetNumberOfNodes_Curve_3D
!
!
!
 SUBROUTINE GetNumberOfNodes_Curve_3D( myCurve_3D, N )
 ! S/R GetNumberOfNodes
 !
 ! =============================================================================================== !
 ! DECLARATIONS
  IMPLICIT NONE
  CLASS( Curve_3D ), INTENT(in) :: myCurve_3D
  INTEGER, INTENT(out)       :: N
  
      N = myCurve_3D % nNodes

 END SUBROUTINE GetNumberOfNodes_Curve_3D
!
!
!
 SUBROUTINE SetInterpNodes_Curve_3D( myCurve_3D, nodes )
 ! S/R SetInterpNodes
 !
 ! =============================================================================================== !
 ! DECLARATIONS
  IMPLICIT NONE
  CLASS( Curve_3D ), INTENT(inout) :: myCurve_3D
  REAL(prec), INTENT(in)        :: nodes(0:myCurve_3D % nNodes)
  
      CALL myCurve_3D % interp % SetNodes( nodes )

 END SUBROUTINE SetInterpNodes_Curve_3D
!
!
!
 SUBROUTINE GetInterpNodes_Curve_3D( myCurve_3D, nodes )
 ! S/R GetInterpNodes
 !
 ! =============================================================================================== !
 ! DECLARATIONS
  IMPLICIT NONE
  CLASS( Curve_3D ), INTENT(in) :: myCurve_3D
  REAL(prec), INTENT(out)    :: nodes(0:myCurve_3D % nNodes)
  
      CALL myCurve_3D % interp % GetNodes( nodes )

 END SUBROUTINE GetInterpNodes_Curve_3D
!
!
!
 SUBROUTINE SetNodes_Curve_3D( myCurve_3D, x, y, z )
 ! S/R SetNodes
 !
 ! =============================================================================================== !
 ! DECLARATIONS
  IMPLICIT NONE
  CLASS( Curve_3D ), INTENT(inout)    :: myCurve_3D
  REAL(prec), INTENT(in), OPTIONAL :: x(0:myCurve_3D % nNodes)
  REAL(prec), INTENT(in), OPTIONAL :: y(0:myCurve_3D % nNodes)
  REAL(prec), INTENT(in), OPTIONAL :: z(0:myCurve_3D % nNodes)
  
      IF( PRESENT(x) .AND. PRESENT(y) .AND. PRESENT(z) )THEN
         myCurve_3D % x = x
         myCurve_3D % y = y
         myCurve_3D % z = z
      ELSE
         myCurve_3D % x = defaultXY
         myCurve_3D % y = defaultXY
         myCurve_3D % z = defaultXY
      ENDIF

 END SUBROUTINE SetNodes_Curve_3D
!
!
!
 SUBROUTINE GetNodes_Curve_3D( myCurve_3D, x, y, z )
 ! S/R GetNodes
 !
 ! =============================================================================================== !
 ! DECLARATIONS
  IMPLICIT NONE
  CLASS( Curve_3D ), INTENT(in) :: myCurve_3D
  REAL(prec), INTENT(out)    :: x(0:myCurve_3D % nNodes)
  REAL(prec), INTENT(out)    :: y(0:myCurve_3D % nNodes)
  REAL(prec), INTENT(out)    :: z(0:myCurve_3D % nNodes)
  
      x = myCurve_3D % x
      y = myCurve_3D % y 
      z = myCurve_3D % z

 END SUBROUTINE GetNodes_Curve_3D
!
!
!
 SUBROUTINE SetSlope_Curve_3D( myCurve_3D, dxds, dyds, dzds )
 ! S/R SetSlope
 !
 ! =============================================================================================== !
 ! DECLARATIONS
  IMPLICIT NONE
  CLASS( Curve_3D ), INTENT(inout) :: myCurve_3D
  REAL(prec), INTENT(in)        :: dxds(0:myCurve_3D % nNodes)
  REAL(prec), INTENT(in)        :: dyds(0:myCurve_3D % nNodes)
  REAL(prec), INTENT(in)        :: dzds(0:myCurve_3D % nNodes)
  
      myCurve_3D % dxds = dxds
      myCurve_3D % dyds = dyds
      myCurve_3D % dzds = dzds

 END SUBROUTINE SetSlope_Curve_3D
!
!
!
 SUBROUTINE GetSlope_Curve_3D( myCurve_3D, dxds, dyds, dzds )
 ! S/R GetSlope
 !
 ! =============================================================================================== !
 ! DECLARATIONS
  IMPLICIT NONE
  CLASS( Curve_3D ), INTENT(in) :: myCurve_3D
  REAL(prec), INTENT(out)    :: dxds(0:myCurve_3D % nNodes)
  REAL(prec), INTENT(out)    :: dyds(0:myCurve_3D % nNodes)
  REAL(prec), INTENT(out)    :: dzds(0:myCurve_3D % nNodes)
  
      dxds = myCurve_3D % dxds
      dyds = myCurve_3D % dyds 
      dzds = myCurve_3D % dzds

 END SUBROUTINE GetSlope_Curve_3D
!
!
!==================================================================================================!
!--------------------------------- Type Specific Routines -----------------------------------------!
!==================================================================================================!
!
!
 SUBROUTINE CalculateSlope_Curve_3D( myCurve_3D )
 ! S/R CalculateSlope
 !
 ! =============================================================================================== !
 ! DECLARATIONS
  IMPLICIT NONE
  CLASS( Curve_3D ), INTENT(inout) :: myCurve_3D
  ! LOCAL
  REAL(prec), ALLOCATABLE :: dMat(:,:)
  
     CALL myCurve_3D % interp % CalculateDerivativeMatrix( dMat )
     
     CALL myCurve_3D % SetSlope( MATMUL(dMat, myCurve_3D % x), &
                              MATMUL(dMat, myCurve_3D % y), &
                              MATMUL(dMat, myCurve_3D % z) )      

     DEALLOCATE( dMat )
     
 END SUBROUTINE CalculateSlope_Curve_3D
!
!
!
 SUBROUTINE Evaluate_Curve_3D( myCurve_3D, s, x, y, z )
 ! S/R Evaluate
 !
 ! =============================================================================================== !
 ! DECLARATIONS
  IMPLICIT NONE
  CLASS( Curve_3D ), INTENT(in) :: myCurve_3D
  REAL(prec), INTENT(in)     :: s 
  REAL(prec), INTENT(out)    :: x, y, z
  
      x = myCurve_3D % interp % EvaluateInterpolant( s, myCurve_3D % x )
      y = myCurve_3D % interp % EvaluateInterpolant( s, myCurve_3D % y )
      z = myCurve_3D % interp % EvaluateInterpolant( s, myCurve_3D % z )
      
 END SUBROUTINE Evaluate_Curve_3D
!
!
!
 SUBROUTINE EvaluateSlope_Curve_3D( myCurve_3D, s, dxds, dyds, dzds )
 ! S/R Evaluate
 !
 ! =============================================================================================== !
 ! DECLARATIONS
  IMPLICIT NONE
  CLASS( Curve_3D ), INTENT(in) :: myCurve_3D
  REAL(prec), INTENT(in)     :: s 
  REAL(prec), INTENT(out)    :: dxds, dyds, dzds
  
      dxds = myCurve_3D % interp % EvaluateInterpolant( s, myCurve_3D % dxds )
      dyds = myCurve_3D % interp % EvaluateInterpolant( s, myCurve_3D % dyds )
      dzds = myCurve_3D % interp % EvaluateInterpolant( s, myCurve_3D % dzds )

 END SUBROUTINE EvaluateSlope_Curve_3D
!
!
!
END MODULE CurveClass_3D
