MODULE GeometrySupportRoutines
! GeometrySupportRoutines.f90 (v2.1 - 12 Dec. 2015) *New with source code version 2.1
! 
! schoonover.numerics@gmail.com
! 
! o (ver 2.1) Dec 2015
!
!  This module contains routines that are useful for computing mappings
!
! 
!  
! ================================================================================================ !
! src/common/
USE ModelPrecision
USE ConstantsDictionary
! src/interp/
USE Lagrange_1D_Class
! src/geom/
USE CurveClass_2D
USE CurveClass_3D



 INTERFACE TransfiniteInterpolation 
    MODULE PROCEDURE :: TransfiniteInterpolation_2D
 END INTERFACE TransfiniteInterpolation
 
 INTERFACE TransfiniteInterpolationMetrics 
    MODULE PROCEDURE :: TransfiniteInterpolationMetrics_2D
 END INTERFACE TransfiniteInterpolationMetrics
 
 INTERFACE Unidirectional
    MODULE PROCEDURE :: Unidirectional_2D
 END INTERFACE Unidirectional 
 
 INTERFACE UnidirectionalDerivative 
    MODULE PROCEDURE :: UnidirectionalDerivative_2D
 END INTERFACE UnidirectionalDerivative 

CONTAINS

 FUNCTION TransfiniteInterpolation_2D( curves, a, b ) RESULT( P )
 ! TransfiniteInterpolation_2D
 !  Takes in the four curves (south, east, north, west) and evaluates the 
 !  bidirectional mapping at xi^1 = a, xi^2 = b. The south and north curves
 !  are assumed to be functions of xi^1, and are located at xi^2 = -1,1 respectively.
 !
 !   An area in 2-D is the assumed geometrical object generated from this evaluation.
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   TYPE( Curve_2D ) :: curves(1:4)
   REAL(prec)       :: a, b
   REAL(prec)       :: P(1:2)
   ! LOCAL
   REAL(prec)  :: P1(1:2), P2(1:2)
   REAL(prec)  :: Psouth(1:2), Pnorth(1:2)
   REAL(prec)  :: leftCurve(1:2), rightCurve(1:2)  
   REAL(prec)  :: c1(1:2), c2(1:2), c3(1:2), c4(1:2)
 

     ! Obtain the corner node locations (Tensor product portion of boolean summation)
     CALL curves(1) % Evaluate( -ONE, c1(1), c1(2) ) ! southwest
     CALL curves(1) % Evaluate( ONE, c2(1), c2(2) ) ! southeast
     CALL curves(3) % Evaluate( -ONE, c3(1), c3(2) ) ! northwest
     CALL curves(3) % Evaluate( ONE, c4(1), c4(2) ) ! northeast
 
   ! Do the unidirectional interpolation between the east and west curves at xi^2 = b
     CALL curves(4) % Evaluate( b, leftCurve(1), leftCurve(2) ) ! west curve
     CALL curves(2) % Evaluate( b, rightCurve(1), rightCurve(2) ) ! east curve

     P1 = Unidirectional_2D( leftCurve, rightCurve, a )

     CALL curves(1) % Evaluate( a, leftCurve(1), leftCurve(2) ) ! south curve
     CALL curves(3) % Evaluate( a, rightCurve(1), rightCurve(2) ) ! north curve

     P2 = Unidirectional( leftCurve, rightCurve, b )

     ! Use the recursive formula
     Psouth = Unidirectional( c1, c2, a )
     Pnorth = Unidirectional( c3, c4, a )
   
     P = P1 + P2 - Unidirectional_2D( Psouth, Pnorth, b ) 
     

 END FUNCTION TransfiniteInterpolation_2D
!
!
!
 FUNCTION TransfiniteInterpolationMetrics_2D( curves, a, b ) RESULT( metricTensor )
 ! FUNCTION  TransfiniteInterpolationMetrics_2D
 !  
 ! ================================================================================
 ! DECLARATIONS
   IMPLICIT NONE
   TYPE( Curve_2D ) :: curves(1:4)
   REAL(prec)       :: a, b
   REAL(prec)       :: metricTensor(1:2,1:2)
   ! LOCAL
   REAL(prec)  :: P1(1:2), P2(1:2), Psouth(1:2), Pnorth(1:2), Pwest(1:2), Peast(1:2)
   REAL(prec)  :: westCurve(1:2), eastCurve(1:2)  
   REAL(prec)  :: southCurve(1:2), northCurve(1:2)
   REAL(prec)  ::  dWest(1:2), dEast(1:2), dSouth(1:2), dNorth(1:2)
   REAL(prec)  :: c1(1:2), c2(1:2), c3(1:2), c4(1:2)
 
 
     CALL curves(4) % Evaluate( b, westCurve(1), westCurve(2) ) ! west curve
     CALL curves(2) % Evaluate( b, eastCurve(1), eastCurve(2) ) ! east curve
     CALL curves(1) % Evaluate( a, southCurve(1), southCurve(2) ) ! south curve
     CALL curves(3) % Evaluate( a, northCurve(1), northCurve(2) ) ! north curve

     ! Obtain the curve derivatives
     CALL curves(4) % EvaluateSlope( b, dWest(1), dWest(2) ) ! west curve
     CALL curves(2) % EvaluateSlope( b, dEast(1), dEast(2) ) ! east curve
     CALL curves(1) % EvaluateSlope( a, dSouth(1), dSouth(2) ) ! south curve
     CALL curves(3) % EvaluateSlope( a, dNorth(1), dNorth(2) ) ! north curve

     ! Obtain the corner node locations (Tensor product portion of boolean summation)
     CALL curves(1) % Evaluate( -ONE, c1(1), c1(2) ) ! southwest
     CALL curves(1) % Evaluate( ONE, c2(1), c2(2) ) ! southeast
     CALL curves(3) % Evaluate( -ONE, c3(1), c3(2) ) ! northwest
     CALL curves(3) % Evaluate( ONE, c4(1), c4(2) ) ! northeast
 
   ! Do the derivative wrt the first computational space variable
     P1 = UnidirectionalDerivative( westCurve, eastCurve )
     Psouth = UnidirectionalDerivative( c1, c2 )
     Pnorth = UnidirectionalDerivative( c3, c4 )

     metricTensor(1:2,1) = P1 + Unidirectional( dSouth - Psouth, dNorth - Pnorth, b ) ! dxds, dyds

  ! Do the derivative wrt the second computational space variable 
     P2 = UnidirectionalDerivative( southCurve, northCurve )
     Pwest = UnidirectionalDerivative( c1, c3 )
     Peast = UnidirectionalDerivative( c2, c4 )

     metricTensor(1:2,2) = P2 + Unidirectional( dWest - Pwest, dEast - Peast, a ) ! dxdp, dydp


 END FUNCTION TransfiniteInterpolationMetrics_2D
!
!
!
 FUNCTION Unidirectional_2D( valLeft, valRight, a ) RESULT( P )
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   REAL(prec) :: valLeft(1:2), valRight(1:2)
   REAL(prec) :: a
   REAL(prec) :: P(1:2)

       P = HALF*( (ONE - a)*valLeft + (ONE + a)*valRight )
    
 END FUNCTION Unidirectional_2D
!
!
!
 FUNCTION Unidirectional_3D( valLeft, valRight, a ) RESULT( P )
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   REAL(prec) :: valLeft(1:3), valRight(1:3)
   REAL(prec) :: a
   REAL(prec) :: P(1:3)

       P = HALF*( (ONE - a)*valLeft + (ONE + a)*valRight )
    
 END FUNCTION Unidirectional_3D
!
!
!
 FUNCTION LinearBlend( a ) RESULT( weights )
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   REAL(prec) :: a
   REAL(prec) :: weights(1:2)

       weights(1) = HALF*(ONE - a)
       weights(2) = HALF*(ONE + a)
    
 END FUNCTION LinearBlend
!
!
!
 FUNCTION UnidirectionalDerivative_2D( valLeft, valRight ) RESULT( P )
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   REAL(prec) :: valLeft(1:2), valRight(1:2)
   REAL(prec) :: P(1:2)

       P = HALF*( valRight - valLeft )
    
 END FUNCTION UnidirectionalDerivative_2D
!
!
!
 FUNCTION UnidirectionalDerivative_3D( valLeft, valRight ) RESULT( P )
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   REAL(prec) :: valLeft(1:3), valRight(1:3)
   REAL(prec) :: P(1:3)

       P = HALF*( valRight - valLeft )
    
 END FUNCTION UnidirectionalDerivative_3D
!
!
!
END Module GeometrySupportRoutines
