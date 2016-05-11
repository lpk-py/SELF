! CurveClass_2D.f90
! 
! Copyright 2015 Joseph Schoonover <schoonover.numerics@gmail.com>
! 
! CurveClass_2D.f90 is part of the Spectral Element Libraries in Fortran (SELF).
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
 
 
MODULE CurveClass_2D
! ========================================= Logs ================================================= !
!2016-05-11  Joseph Schoonover  schoonover.numerics@gmail.com 
!
! //////////////////////////////////////////////////////////////////////////////////////////////// ! 
! src/common/
USE ModelPrecision
USE ConstantsDictionary
! src/interp/
USE Lagrange_1D_Class


   TYPE Curve_2D 
      INTEGER, PRIVATE                  :: nNodes
      TYPE(Lagrange_1D), PRIVATE        :: interp
      REAL(prec),  ALLOCATABLE, PRIVATE :: x(:), y(:)
      REAL(prec),  ALLOCATABLE, PRIVATE :: dxds(:), dyds(:)

      CONTAINS

      PROCEDURE :: Build => Build_Curve_2D
      PROCEDURE :: Trash => Trash_Curve_2D

      PROCEDURE :: SetNumberOfNodes => SetNumberOfNodes_Curve_2D
      PROCEDURE :: GetNumberOfNodes => GetNumberOfNodes_Curve_2D
      PROCEDURE :: SetInterpNodes => SetInterpNodes_Curve_2D
      PROCEDURE :: GetInterpNodes => GetInterpNodes_Curve_2D
      PROCEDURE :: SetNodes => SetNodes_Curve_2D
      PROCEDURE :: GetNodes => GetNodes_Curve_2D
      PROCEDURE :: SetSlope => SetSlope_Curve_2D
      PROCEDURE :: GetSlope => GetSlope_Curve_2D

      PROCEDURE :: CalculateSlope => CalculateSlope_Curve_2D
      PROCEDURE :: Evaluate => Evaluate_Curve_2D
      PROCEDURE :: EvaluateSlope => EvaluateSlope_Curve_2D
      
      
   END TYPE Curve_2D


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
 SUBROUTINE Build_Curve_2D( myCurve, x, y, nodes )
 ! S/R Build
 !
 ! =============================================================================================== !
 ! DECLARATIONS
  IMPLICIT NONE
  CLASS( Curve_2D ), INTENT(out)   :: myCurve
  REAL(prec), INTENT(in), OPTIONAL :: x(0:), y(0:), nodes(0:)
  ! LOCAL
  INTEGER :: N
  
      IF( PRESENT(x) .AND. PRESENT(y) .AND. PRESENT(nodes) )THEN
         N = UBOUND(x, DIM=1 )
      ELSE
         N = defaultN
      ENDIF
      
      CALL myCurve % SetNumberOfNodes( N )
      
      IF( PRESENT(nodes) )THEN
         CALL myCurve % interp % Build( N, nodes )
      ELSE
         CALL myCurve % interp % Build( N )
      ENDIF

      ALLOCATE( myCurve % x(0:N), mycurve % y(0:N) )

      ALLOCATE( myCurve % dxds(0:N), mycurve % dyds(0:N) )

      IF( PRESENT(x) .AND. PRESENT(y) )THEN
         CALL myCurve % SetNodes( x, y )
      ELSE
         CALL myCurve % SetNodes( )
      ENDIF
      
      CALL myCurve % CalculateSlope( )
 
 END SUBROUTINE Build_Curve_2D
!
!
!
 SUBROUTINE Trash_Curve_2D( myCurve )
 ! S/R Trash
 !
 ! =============================================================================================== !
 ! DECLARATIONS
  IMPLICIT NONE
  CLASS( Curve_2D ), INTENT(inout) :: myCurve
  
      CALL myCurve % interp % Trash( )

      DEALLOCATE( myCurve % x, mycurve % y )

      DEALLOCATE( myCurve % dxds, mycurve % dyds )

 
 END SUBROUTINE Trash_Curve_2D
!
!
!==================================================================================================!
!--------------------------------------- Accessors ------------------------------------------------!
!==================================================================================================!
!
!
 SUBROUTINE SetNumberOfNodes_Curve_2D( myCurve, N )
 ! S/R SetNumberOfNodes
 !
 ! =============================================================================================== !
 ! DECLARATIONS
  IMPLICIT NONE
  CLASS( Curve_2D ), INTENT(inout) :: myCurve
  INTEGER, INTENT(in)              :: N
  
      myCurve % nNodes = N

 END SUBROUTINE SetNumberOfNodes_Curve_2D
!
!
!
 SUBROUTINE GetNumberOfNodes_Curve_2D( myCurve, N )
 ! S/R GetNumberOfNodes
 !
 ! =============================================================================================== !
 ! DECLARATIONS
  IMPLICIT NONE
  CLASS( Curve_2D ), INTENT(in) :: myCurve
  INTEGER, INTENT(out)          :: N
  
      N = myCurve % nNodes

 END SUBROUTINE GetNumberOfNodes_Curve_2D
!
!
!
 SUBROUTINE SetInterpNodes_Curve_2D( myCurve, nodes )
 ! S/R SetInterpNodes
 !
 ! =============================================================================================== !
 ! DECLARATIONS
  IMPLICIT NONE
  CLASS( Curve_2D ), INTENT(inout) :: myCurve
  REAL(prec), INTENT(in)           :: nodes(0:myCurve % nNodes)
  
      CALL myCurve % interp % SetNodes( nodes )

 END SUBROUTINE SetInterpNodes_Curve_2D
!
!
!
 SUBROUTINE GetInterpNodes_Curve_2D( myCurve, nodes )
 ! S/R GetInterpNodes
 !
 ! =============================================================================================== !
 ! DECLARATIONS
  IMPLICIT NONE
  CLASS( Curve_2D ), INTENT(in) :: myCurve
  REAL(prec), INTENT(out)       :: nodes(0:myCurve % nNodes)
  
      CALL myCurve % interp % GetNodes( nodes )

 END SUBROUTINE GetInterpNodes_Curve_2D
!
!
!
 SUBROUTINE SetNodes_Curve_2D( myCurve, x, y )
 ! S/R SetNodes
 !
 ! =============================================================================================== !
 ! DECLARATIONS
  IMPLICIT NONE
  CLASS( Curve_2D ), INTENT(inout) :: myCurve
  REAL(prec), INTENT(in), OPTIONAL :: x(0:myCurve % nNodes), y(0:myCurve % nNodes)
  
      IF( PRESENT(x) .AND. PRESENT(y) )THEN
         myCurve % x = x
         myCurve % y = y
      ELSE
         myCurve % x = defaultXY
         myCurve % y = defaultXY
      ENDIF

 END SUBROUTINE SetNodes_Curve_2D
!
!
!
 SUBROUTINE GetNodes_Curve_2D( myCurve, x, y )
 ! S/R GetNodes
 !
 ! =============================================================================================== !
 ! DECLARATIONS
  IMPLICIT NONE
  CLASS( Curve_2D ), INTENT(in) :: myCurve
  REAL(prec), INTENT(out)       :: x(0:myCurve % nNodes), y(0:myCurve % nNodes)
  
      x = myCurve % x
      y = myCurve % y 

 END SUBROUTINE GetNodes_Curve_2D
!
!
!
 SUBROUTINE SetSlope_Curve_2D( myCurve, dxds, dyds )
 ! S/R SetSlope
 !
 ! =============================================================================================== !
 ! DECLARATIONS
  IMPLICIT NONE
  CLASS( Curve_2D ), INTENT(inout) :: myCurve
  REAL(prec), INTENT(in)           :: dxds(0:myCurve % nNodes), dyds(0:myCurve % nNodes)
  
      myCurve % dxds = dxds
      myCurve % dyds = dyds

 END SUBROUTINE SetSlope_Curve_2D
!
!
!
 SUBROUTINE GetSlope_Curve_2D( myCurve, dxds, dyds )
 ! S/R GetSlope
 !
 ! =============================================================================================== !
 ! DECLARATIONS
  IMPLICIT NONE
  CLASS( Curve_2D ), INTENT(in) :: myCurve
  REAL(prec), INTENT(out)       :: dxds(0:myCurve % nNodes), dyds(0:myCurve % nNodes)
  
      dxds = myCurve % dxds
      dyds = myCurve % dyds 

 END SUBROUTINE GetSlope_Curve_2D
!
!
!==================================================================================================!
!--------------------------------- Type Specific Routines -----------------------------------------!
!==================================================================================================!
!
!
 SUBROUTINE CalculateSlope_Curve_2D( myCurve )
 ! S/R CalculateSlope
 !
 ! =============================================================================================== !
 ! DECLARATIONS
  IMPLICIT NONE
  CLASS( Curve_2D ), INTENT(inout) :: myCurve
  ! LOCAL
  REAL(prec), ALLOCATABLE :: dMat(:,:)
  
     CALL myCurve % interp % CalculateDerivativeMatrix( dMat )
     
     CALL myCurve % SetSlope( MATMUL(dMat, myCurve % x), MATMUL(dMat, myCurve % y) )      

     DEALLOCATE( dMat )
     
 END SUBROUTINE CalculateSlope_Curve_2D
!
!
!
 SUBROUTINE Evaluate_Curve_2D( myCurve, s, x, y )
 ! S/R Evaluate
 !
 ! =============================================================================================== !
 ! DECLARATIONS
  IMPLICIT NONE
  CLASS( Curve_2D ), INTENT(in) :: myCurve
  REAL(prec), INTENT(in)        :: s 
  REAL(prec), INTENT(out)       :: x, y
  
      x = myCurve % interp % EvaluateInterpolant( s, myCurve % x )
      y = myCurve % interp % EvaluateInterpolant( s, myCurve % y )

 END SUBROUTINE Evaluate_Curve_2D
!
!
!
 SUBROUTINE EvaluateSlope_Curve_2D( myCurve, s, dxds, dyds )
 ! S/R Evaluate
 !
 ! =============================================================================================== !
 ! DECLARATIONS
  IMPLICIT NONE
  CLASS( Curve_2D ), INTENT(in) :: myCurve
  REAL(prec), INTENT(in)        :: s 
  REAL(prec), INTENT(out)       :: dxds, dyds
  
      dxds = myCurve % interp % EvaluateDerivative( s, myCurve % x )
      dyds = myCurve % interp % EvaluateDerivative( s, myCurve % y )

 END SUBROUTINE EvaluateSlope_Curve_2D
!
!
!
END MODULE CurveClass_2D
