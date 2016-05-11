! CurveClass_1D.f90
! 
! Copyright 2015 Joseph Schoonover <schoonover.numerics@gmail.com>
! 
! CurveClass_1D.f90 is part of the Spectral Element Libraries in Fortran (SELF).
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
 
 
MODULE CurveClass_1D
! ========================================= Logs ================================================= !
!2016-05-11  Joseph Schoonover  schoonover.numerics@gmail.com 
!
! //////////////////////////////////////////////////////////////////////////////////////////////// ! 
! src/common/
USE ModelPrecision
USE ConstantsDictionary
! src/interp/
USE Lagrange_1D_Class


   TYPE Curve_1D 
      INTEGER, PRIVATE                  :: nNodes
      TYPE(Lagrange_1D), PRIVATE        :: interp
      REAL(prec),  ALLOCATABLE, PRIVATE :: x(:)
      REAL(prec),  ALLOCATABLE, PRIVATE :: dxds(:)

      CONTAINS

      PROCEDURE :: Build => Build_Curve_1D
      PROCEDURE :: Trash => Trash_Curve_1D

      PROCEDURE :: SetNumberOfNodes => SetNumberOfNodes_Curve_1D
      PROCEDURE :: GetNumberOfNodes => GetNumberOfNodes_Curve_1D
      PROCEDURE :: SetInterpNodes => SetInterpNodes_Curve_1D
      PROCEDURE :: GetInterpNodes => GetInterpNodes_Curve_1D
      PROCEDURE :: SetNodes => SetNodes_Curve_1D
      PROCEDURE :: GetNodes => GetNodes_Curve_1D
      PROCEDURE :: SetSlope => SetSlope_Curve_1D
      PROCEDURE :: GetSlope => GetSlope_Curve_1D

      PROCEDURE :: CalculateSlope => CalculateSlope_Curve_1D
      PROCEDURE :: Evaluate => Evaluate_Curve_1D
      PROCEDURE :: EvaluateSlope => EvaluateSlope_Curve_1D
      
      
   END TYPE Curve_1D


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
 SUBROUTINE Build_Curve_1D( myCurve_1D, x, nodes )
 ! S/R Build
 !
 ! =============================================================================================== !
 ! DECLARATIONS
  IMPLICIT NONE
  CLASS( Curve_1D ), INTENT(out)   :: myCurve_1D
  REAL(prec), INTENT(in), OPTIONAL :: x(0:), nodes(0:)
  ! LOCAL
  INTEGER :: N
  
      IF( PRESENT(x) .AND. PRESENT(nodes) )THEN
         N = UBOUND(x, DIM=1 )
      ELSE
         N = defaultN
      ENDIF
      
      CALL myCurve_1D % SetNumberOfNodes( N )
      
      IF( PRESENT(nodes) )THEN
         CALL myCurve_1D % interp % Build( N, nodes )
      ELSE
         CALL myCurve_1D % interp % Build( N )
      ENDIF

      ALLOCATE( myCurve_1D % x(0:N) )

      ALLOCATE( myCurve_1D % dxds(0:N) )

      IF( PRESENT(x) )THEN
         CALL myCurve_1D % SetNodes( x )
      ELSE
         CALL myCurve_1D % SetNodes( )
      ENDIF
      
      CALL myCurve_1D % CalculateSlope( )
 
 END SUBROUTINE Build_Curve_1D
!
!
!
 SUBROUTINE Trash_Curve_1D( myCurve_1D )
 ! S/R Trash
 !
 ! =============================================================================================== !
 ! DECLARATIONS
  IMPLICIT NONE
  CLASS( Curve_1D ), INTENT(inout) :: myCurve_1D
  
      CALL myCurve_1D % interp % Trash( )

      DEALLOCATE( myCurve_1D % x )

      DEALLOCATE( myCurve_1D % dxds )

 
 END SUBROUTINE Trash_Curve_1D
!
!
!==================================================================================================!
!--------------------------------------- Accessors ------------------------------------------------!
!==================================================================================================!
!
!
 SUBROUTINE SetNumberOfNodes_Curve_1D( myCurve_1D, N )
 ! S/R SetNumberOfNodes
 !
 ! =============================================================================================== !
 ! DECLARATIONS
  IMPLICIT NONE
  CLASS( Curve_1D ), INTENT(inout) :: myCurve_1D
  INTEGER, INTENT(in)              :: N
  
      myCurve_1D % nNodes = N

 END SUBROUTINE SetNumberOfNodes_Curve_1D
!
!
!
 SUBROUTINE GetNumberOfNodes_Curve_1D( myCurve_1D, N )
 ! S/R GetNumberOfNodes
 !
 ! =============================================================================================== !
 ! DECLARATIONS
  IMPLICIT NONE
  CLASS( Curve_1D ), INTENT(in) :: myCurve_1D
  INTEGER, INTENT(out)          :: N
  
      N = myCurve_1D % nNodes

 END SUBROUTINE GetNumberOfNodes_Curve_1D
!
!
!
 SUBROUTINE SetInterpNodes_Curve_1D( myCurve_1D, nodes )
 ! S/R SetInterpNodes
 !
 ! =============================================================================================== !
 ! DECLARATIONS
  IMPLICIT NONE
  CLASS( Curve_1D ), INTENT(inout) :: myCurve_1D
  REAL(prec), INTENT(in)           :: nodes(0:myCurve_1D % nNodes)
  
      CALL myCurve_1D % interp % SetNodes( nodes )

 END SUBROUTINE SetInterpNodes_Curve_1D
!
!
!
 SUBROUTINE GetInterpNodes_Curve_1D( myCurve_1D, nodes )
 ! S/R GetInterpNodes
 !
 ! =============================================================================================== !
 ! DECLARATIONS
  IMPLICIT NONE
  CLASS( Curve_1D ), INTENT(in) :: myCurve_1D
  REAL(prec), INTENT(out)       :: nodes(0:myCurve_1D % nNodes)
  
      CALL myCurve_1D % interp % GetNodes( nodes )

 END SUBROUTINE GetInterpNodes_Curve_1D
!
!
!
 SUBROUTINE SetNodes_Curve_1D( myCurve_1D, x )
 ! S/R SetNodes
 !
 ! =============================================================================================== !
 ! DECLARATIONS
  IMPLICIT NONE
  CLASS( Curve_1D ), INTENT(inout)    :: myCurve_1D
  REAL(prec), INTENT(in), OPTIONAL    :: x(0:myCurve_1D % nNodes)
  
      IF( PRESENT(x) )THEN
         myCurve_1D % x = x
      ELSE
         myCurve_1D % x = defaultXY
      ENDIF

 END SUBROUTINE SetNodes_Curve_1D
!
!
!
 SUBROUTINE GetNodes_Curve_1D( myCurve_1D, x )
 ! S/R GetNodes
 !
 ! =============================================================================================== !
 ! DECLARATIONS
  IMPLICIT NONE
  CLASS( Curve_1D ), INTENT(in) :: myCurve_1D
  REAL(prec), INTENT(out)       :: x(0:myCurve_1D % nNodes)
  
      x = myCurve_1D % x

 END SUBROUTINE GetNodes_Curve_1D
!
!
!
 SUBROUTINE SetSlope_Curve_1D( myCurve_1D, dxds )
 ! S/R SetSlope
 !
 ! =============================================================================================== !
 ! DECLARATIONS
  IMPLICIT NONE
  CLASS( Curve_1D ), INTENT(inout) :: myCurve_1D
  REAL(prec), INTENT(in)           :: dxds(0:myCurve_1D % nNodes)
  
      myCurve_1D % dxds = dxds

 END SUBROUTINE SetSlope_Curve_1D
!
!
!
 SUBROUTINE GetSlope_Curve_1D( myCurve_1D, dxds )
 ! S/R GetSlope
 !
 ! =============================================================================================== !
 ! DECLARATIONS
  IMPLICIT NONE
  CLASS( Curve_1D ), INTENT(in) :: myCurve_1D
  REAL(prec), INTENT(out)       :: dxds(0:myCurve_1D % nNodes)
  
      dxds = myCurve_1D % dxds

 END SUBROUTINE GetSlope_Curve_1D
!
!
!==================================================================================================!
!--------------------------------- Type Specific Routines -----------------------------------------!
!==================================================================================================!
!
!
 SUBROUTINE CalculateSlope_Curve_1D( myCurve_1D )
 ! S/R CalculateSlope
 !
 ! =============================================================================================== !
 ! DECLARATIONS
  IMPLICIT NONE
  CLASS( Curve_1D ), INTENT(inout) :: myCurve_1D
  ! LOCAL
  REAL(prec), ALLOCATABLE :: dMat(:,:)
  
     CALL myCurve_1D % interp % CalculateDerivativeMatrix( dMat )
     
     CALL myCurve_1D % SetSlope( MATMUL(dMat, myCurve_1D % x) )      

     DEALLOCATE( dMat )
     
 END SUBROUTINE CalculateSlope_Curve_1D
!
!
!
 SUBROUTINE Evaluate_Curve_1D( myCurve_1D, s, x )
 ! S/R Evaluate
 !
 ! =============================================================================================== !
 ! DECLARATIONS
  IMPLICIT NONE
  CLASS( Curve_1D ), INTENT(in) :: myCurve_1D
  REAL(prec), INTENT(in)        :: s 
  REAL(prec), INTENT(out)       :: x
  
      x = myCurve_1D % interp % EvaluateInterpolant( s, myCurve_1D % x )

 END SUBROUTINE Evaluate_Curve_1D
!
!
!
 SUBROUTINE EvaluateSlope_Curve_1D( myCurve_1D, s, dxds )
 ! S/R Evaluate
 !
 ! =============================================================================================== !
 ! DECLARATIONS
  IMPLICIT NONE
  CLASS( Curve_1D ), INTENT(in) :: myCurve_1D
  REAL(prec), INTENT(in)        :: s 
  REAL(prec), INTENT(out)       :: dxds
  
      dxds = myCurve_1D % interp % EvaluateInterpolant( s, myCurve_1D % dxds )

 END SUBROUTINE EvaluateSlope_Curve_1D
!
!
!
END MODULE CurveClass_1D
