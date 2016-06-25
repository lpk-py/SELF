! Lagrange_1D_Class.f90
! 
! Copyright 2015 Joseph Schoonover <schoonover.numerics@gmail.com>
! 
! Lagrange_1D_Class.f90 is part of the Spectral Element Libraries in Fortran (SELF).
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
 
 
MODULE Lagrange_1D_Class
! ========================================= Logs ================================================= !
!2016-05-11  Joseph Schoonover  schoonover.numerics@gmail.com 
!
! //////////////////////////////////////////////////////////////////////////////////////////////// ! 

USE ModelPrecision
USE ConstantsDictionary
USE CommonRoutines

IMPLICIT NONE



   TYPE, PUBLIC :: Lagrange_1D
      INTEGER                 :: nS      ! number of nodal points where we have obersvations
      INTEGER                 :: nSo     ! number of nodal points where we want observations
      REAL(prec), ALLOCATABLE :: s(:)    ! locations where we have obervations
      REAL(prec), ALLOCATABLE :: bWs(:)  ! barycentric weights
      REAL(prec), ALLOCATABLE :: so(:)   ! Locations where we want observations
      REAL(prec), ALLOCATABLE :: Ts(:,:) ! Interpolation matrix to get us from what we have to what
                                         ! we want

      CONTAINS
      
      !-------------!
      ! Constructors/Destructors
      PROCEDURE :: Build => BuildLagrange_1D
      PROCEDURE :: Trash => TrashLagrange_1D
  
      ! Accessors
      PROCEDURE :: Sets => SetsLagrange_1D
      PROCEDURE :: Gets => GetsLagrange_1D
      PROCEDURE :: SetNode => SetNodeLagrange_1D
      PROCEDURE :: GetNode => GetNodeLagrange_1D
      PROCEDURE :: SetbWs => SetbWsLagrange_1D
      PROCEDURE :: GetbWs => GetbWsLagrange_1D
      PROCEDURE :: SetWeight => SetWeightLagrange_1D
      PROCEDURE :: GetWeight => GetWeightLagrange_1D
      PROCEDURE :: SetNumberOfs => SetNumberOfsLagrange_1D
      PROCEDURE :: GetNumberOfs => GetNumberOfsLagrange_1D


      ! Data-structure operations
      PROCEDURE :: CalculatebWs => CalculateBarycentricbWs_1D
      PROCEDURE :: EvaluateInterpolant => EvaluateLagrange_1D
      PROCEDURE :: EvaluateLagrangePolynomial => EvaluateLagrangePolynomial_1D
      PROCEDURE :: EvaluateDerivative => EvaluateDerivative_1D
      PROCEDURE :: CalculateDerivativeMatrix => CalculateDerivativeMatrix_1D
      PROCEDURE :: CalculateInterpolationMatrix => CalculateInterpolationMatrix_1D
      PROCEDURE :: CoarseToFine => CoarseToFine_1D
      
      ! File I/O Routines
      PROCEDURE :: WriteCurve => WriteCurve_Lagrange_1D

    END TYPE Lagrange_1D


 INTEGER, PARAMETER, PUBLIC       :: nDefaults = 1          ! Default number of s
 REAL(PREC),PARAMETER, PUBLIC     :: LagrangeNodeDefault = ZERO ! The default value to set the s
 CHARACTER(17), PARAMETER, PUBLIC :: LagrangeCurveFMT = '(E17.7,2x,E17.7)'
 
 CONTAINS
!
!
!==================================================================================================!
!------------------------------- Manual Constructors/Destructors ----------------------------------!
!==================================================================================================!
!
!
 SUBROUTINE BuildLagrange_1D( myPoly, ns, xIn )
 ! S/R BuildLagrange_1D
 !
 !   This subroutine is the basic constructor. If the number of s are input, THEN the desired 
 !   amount of memory is allocated. If the actual interpolation s are input, THEN the s
 !   are stored and the barycentric bWs are calculated. If the number of s are not given,
 !   THEN the number of s are determined from the size of xIn. In the event that xIn is also 
 !   not given, the default number of s is used (nDefaults) and the s and barycentric
 !   bWs are set to zero.
 !
 !
 ! =============================================================================================== !
 ! DECLARATIONS 
   IMPLICIT NONE
   CLASS(Lagrange_1D), INTENT(inout) :: myPoly
   INTEGER, INTENT(in)               :: ns
   REAL(prec), INTENT(in),optional   :: xIn(:)
   

      ! Set the number of s
      CALL myPoly % SetNumberOfs( ns )
      
      ! Allocate storage
      ALLOCATE( myPoly % s(0:ns), myPoly % bWs(0:ns) )
      
      IF( PRESENT(xIn) )THEN
         CALL myPoly % Sets( xIn )
         CALL myPoly % CalculatebWs( )
      ELSE
         CALL myPoly % Sets( )
         CALL myPoly % SetbWs( )
      ENDIF
 
 END SUBROUTINE BuildLagrange_1D
!
!
!
SUBROUTINE TrashLagrange_1D(myPoly)
 ! S/R TrashLagrange_1D
 !
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS(Lagrange_1D), INTENT(inout) :: myPoly


      myPoly % ns = 0

      DEALLOCATE( myPoly % s, myPoly % bWs )

     

 END SUBROUTINE TrashLagrange_1D
!
!
!==================================================================================================!
!--------------------------------------- Accessors ------------------------------------------------!
!==================================================================================================!
!
!
 SUBROUTINE SetsLagrange_1D( myPoly, xIn )
 ! S/R SetsLagrange_1D
 !  
 !    Sets the attribute "s" to xIn. If xIn is not present, THEN the s are set to the 
 !    default value (LagrangeNodeDefault)
 !
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS(Lagrange_1D), INTENT(inout) :: myPoly
   REAL(prec), INTENT(in), OPTIONAL  :: xIn(0:myPoly % ns)

       IF( PRESENT(xIn) )THEN
          myPoly % s = xIn
       ELSE
          myPoly % s = LagrangeNodeDefault
       ENDIF

 END SUBROUTINE SetsLagrange_1D
!
!
!
 SUBROUTINE GetsLagrange_1D( myPoly, xOut )
 ! S/R GetsLagrange_1D
 !  
 !    Gets the attribute "s" and sets xOut.
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS(Lagrange_1D), INTENT(in) :: myPoly
   REAL(prec), INTENT(out)        :: xOut(0:myPoly % ns)


       xOut = myPoly % s 

 END SUBROUTINE GetsLagrange_1D
!
!
!
 SUBROUTINE SetNodeLagrange_1D( myPoly, i, xIn )
 ! S/R SetNodeLagrange_1D
 !  
 !    Sets the i-th node to xIn
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS(Lagrange_1D), INTENT(inout) :: myPoly
   INTEGER, INTENT(in)               :: i
   REAL(prec), INTENT(in)            :: xIn


      myPoly % s(i) = xIn 

 END SUBROUTINE SetNodeLagrange_1D
!
!
!
 SUBROUTINE GetNodeLagrange_1D( myPoly, i, xOut )
 ! S/R GetNodeLagrange_1D
 !  
 !   Gets ths i-th interpolation node and sets xOut to this value
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS(Lagrange_1D), INTENT(in) :: myPoly
   INTEGER, INTENT(in)            :: i
   REAL(prec), INTENT(out)        :: xOut


      xOut = myPoly % s(i) 

 END SUBROUTINE GetNodeLagrange_1D
!
!
!
 SUBROUTINE SetbWsLagrange_1D( myPoly, wIn )
 ! S/R SetbWsLagrange_1D
 !  
 !    Sets the attribute "bWs" to wIn. If wIn is not present, THEN the barycentric bWs are
 !    set to the default value (LagrangeNodeDefault)
 !
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS(Lagrange_1D), INTENT(inout) :: myPoly
   REAL(prec), INTENT(in), OPTIONAL  :: wIn(0:myPoly % ns)

       IF( PRESENT(wIn) )THEN
          myPoly % bWs = wIn
       ELSE
          myPoly % bWs = LagrangeNodeDefault
       ENDIF

 END SUBROUTINE SetbWsLagrange_1D
!
!
!
 SUBROUTINE GetbWsLagrange_1D( myPoly, wOut )
 ! S/R GetbWsLagrange_1D
 !  
 !    Gets the attribute "bWs" and sets wOut.
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS(Lagrange_1D), INTENT(in) :: myPoly
   REAL(prec), INTENT(out)        :: wOut(0:myPoly % ns)


       wOut = myPoly % bWs 

 END SUBROUTINE GetbWsLagrange_1D
!
!
!
 SUBROUTINE SetWeightLagrange_1D( myPoly, i, wIn )
 ! S/R SetWeightLagrange_1D
 !  
 !    Sets the i-th barycentric weight to wIn
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS(Lagrange_1D), INTENT(inout) :: myPoly
   INTEGER, INTENT(in)               :: i
   REAL(prec), INTENT(in)            :: wIn


       myPoly % bWs(i) = wIn 

 END SUBROUTINE SetWeightLagrange_1D
!
!
!
 SUBROUTINE GetWeightLagrange_1D( myPoly, i, wOut )
 ! S/R GetWeightLagrange_1D
 !  
 !   Gets the i-th barycentric weight and sets wOut to this value
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS(Lagrange_1D), INTENT(in) :: myPoly
   INTEGER, INTENT(in)            :: i
   REAL(prec), INTENT(out)        :: wOut


       wOut = myPoly % bWs(i) 

 END SUBROUTINE GetWeightLagrange_1D
!
!
!
 SUBROUTINE SetNumberOfsLagrange_1D( myPoly, n )
 ! S/R SetNumberOfsLagrange_1D(
 !  
 !    Sets the number of s (ns) to n.
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS(Lagrange_1D), INTENT(inout) :: myPoly
   INTEGER, INTENT(in)             :: n
    
       myPoly % ns = n 

 END SUBROUTINE SetNumberOfsLagrange_1D
!
!
!
 SUBROUTINE GetNumberOfsLagrange_1D( myPoly, n )
 ! S/R GetNumberOfsLagrange_1D
 !  
 !    Sets n to the number of interpolation s (ns)
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS(Lagrange_1D), INTENT(in) :: myPoly
   INTEGER, INTENT(out)           :: n


       n = myPoly % ns 

 END SUBROUTINE GetNumberOfsLagrange_1D
!
!
!==================================================================================================!
!--------------------------------- Type Specific Routines -----------------------------------------!
!==================================================================================================!
!
!
 SUBROUTINE CalculateBarycentricbWs_1D( myPoly )
 ! S/R CalculateBarycentricbWs_1D
 !  
 !    Calculates the barycentric bWs from the interpolation s x(0:nP) and stores them in 
 !    the "bWs" attribute.
 !    
 !
! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS(Lagrange_1D), INTENT(inout) :: myPoly
   
      myPoly % 
 END SUBROUTINE CalculateBarycentricbWs_1D
!
!
!
 FUNCTION EvaluateLagrange_1D( myPoly, xEval, f ) RESULT( inFatX )  
 ! S/R EvaluateLagrange_1D 
 !  
 !    Evaluates the Lagrange interpolating polynomial described by the Lagrange_1D data-structure
 !    "myPoly" and has the nodal values contained in "f" at the location "xEval"
 !
 ! =============================================================================================== !
   IMPLICIT NONE
   CLASS(Lagrange_1D), INTENT(in) :: myPoly
   REAL(prec), INTENT(in)         :: xEval
   REAL(prec), INTENT(in)         :: f(0:myPoly%ns)
   REAL(prec)                     :: inFatX 
   ! LOCAL
   REAL(prec) :: num, den, t, thisX, thisW
   INTEGER    :: jP, nP 
 
     CALL myPoly % GetNumberOfs( nP )
     num = ZERO
     den = ZERO

     DO jP = 0, nP ! loop over the interpolation s

        CALL myPoly % GetNode( jP, thisX )

        IF( AlmostEqual(xEval, thisX) ) THEN !xEval is an interpolation node
          
          inFatX = f(jP)
          RETURN
  
        ELSE ! THEN xEval is not an interpolation node

           CALL myPoly % GetWeight( jP, thisW )
           
           t = thisW/(xEval - thisX)

           num = num + t*f(jP)
           den = den + t

        ENDIF  ! IF the evaluation point is an interpolation node
        

     ENDDO ! jP, loop over the interpolation s


     inFatX = num/den

     RETURN

 END FUNCTION EvaluateLagrange_1D
!
!
!
 FUNCTION EvaluateLagrangePolynomial_1D( myPoly, xEval ) RESULT( lAtX )  
 ! FUNCTION EvaluateLagrangePolynomial_1D ( LAGrange_POLYnomial)
 !  
 !    Evaluates the lagrange interpolating polynomial at the point x.
 !
 ! 
 ! =============================================================================================== !
  IMPLICIT NONE
  CLASS(Lagrange_1D), INTENT(in) :: myPoly
  REAL(prec), INTENT(in)         :: xEval
  REAL(prec)                     :: lAtX(0:myPoly%ns)
  ! LOCAL
  REAL(prec) :: temp1, temp2, thisX, thisW
  INTEGER    :: jP, nP
  LOGICAL    :: xMatchesNode

     CALL myPoly % GetNumberOfs( nP )
     
     xMatchesNode = .FALSE.

     DO jP = 0, nP ! loop over the interpolation s

        CALL myPoly % GetNode( jP, thisX )
        
        lAtX(jP) = ZERO

        IF( AlmostEqual(xEval, thisX) ) THEN !xEval is an interpolation node

          lAtX(jP) = ONE
          
          xMatchesNode = .TRUE.
  
        ENDIF 

     ENDDO ! jP, loop over the interpolation s

     
     IF( xMatchesNode )THEN ! we're done, the evaluation point is one of the s

        RETURN

     ENDIF
     
     ! Otherwise, we calculate the values for the lagrange polynomials

     temp1 = ZERO
     
     DO jP = 0, nP ! loop over the interpolation s

        CALL myPoly % GetWeight( jP, thisW )
        CALL myPoly % GetNode( jP, thisX )

        temp2 = thisW/(xEval - thisX)

        lAtX(jP) = temp2

        temp1 = temp1 + temp2

     ENDDO ! jP, loop over the interpolation s
     
     DO jP = 0, nP ! loop over the interpolation s (again)

        lAtX(jP) = lAtX(jP)/temp1

     ENDDO ! jP, loop over the interpolation s (again)
     

 END FUNCTION EvaluateLagrangePolynomial_1D
!
!
!
 FUNCTION EvaluateDerivative_1D(myPoly, xEval, f) RESULT(dInFdx)  
 ! FUNCTION EvaluateDerivative_1D 
 !  
 !    Evaluates the derivative of the Lagrange interpolating polynomial described by the Lagrange_1D
 !    data-structure "myPoly" and has the nodal values contained in "f" at the location "xEval" 
 ! 
 ! =============================================================================================== !
  IMPLICIT NONE
  CLASS(Lagrange_1D), INTENT(in) :: myPoly
  REAL(prec), INTENT(in)         :: xEval
  REAL(prec), INTENT(in)         :: f(0:myPoly%ns)
  REAL(prec)                     :: dInFdx
  ! LOCAL
  REAL(prec) :: num, den, t, p, thisX, thisW
  INTEGER    :: jP, jS, nP 
  LOGICAL    :: atNode 

     CALL myPoly % GetNumberOfs( nP )

     num = ZERO
     atNode = .FALSE.

     DO jP = 0, nP ! loop over the interpolation s

        CALL myPoly % GetNode( jP, thisX )

        IF( AlmostEqual(xEval, thisX) ) THEN !xEval is an interpolation node
     
           atNode = .TRUE.

           p = f(jP)
           
           CALL myPoly % GetWeight( jP, thisW )
           
           den = -thisW
   
           jS = jP

        ENDIF  ! IF the evaluation point is an interpolation node
        
     ENDDO ! jP, loop over the interpolation s



     IF( atNode ) THEN ! the evaluation location is an interpolating node


        DO jP = 0, nP ! loop over the interpolation s

           IF( .NOT.(jP == jS) ) THEN ! x(jP) is not the interpolating node

             CALL myPoly % GetNode( jP, thisX )
             CALL myPoly % GetWeight( jP, thisW )
             
             num = num + thisW*(p - f(jP))/(xEval - thisX)

           ENDIF

        ENDDO ! jP, loop over the interpolation s


     ELSE ! the evaluation location is not an interpolating node


        den = ZERO

        p = EvaluateLagrange_1D(myPoly, xEval, f) 
        
        DO jP = 0, nP ! loop over the interpolation s

           CALL myPoly % GetNode( jP, thisX )
           CALL myPoly % GetWeight( jP, thisW )

           t = thisW/(xEval - thisX)

           num = num + t*(p - f(jP))/(xEval - thisX)

           den = den + t

        ENDDO ! jP, loop over the interpolation s


     ENDIF ! conditional, IF we're on an interpolating node

     
     dInFdx = num/den

     RETURN

 END FUNCTION EvaluateDerivative_1D
!
!
!
SUBROUTINE CalculateDerivativeMatrix_1D( myPoly, dMat, otherInterpolant )  
 ! S/R CalculateDerivativeMatrix 
 !  Description :
 !  Calculates the derivative matrix "Dmat(:,:)"for calculating the derivative
 !  at the interpolation s "myPoly%s(:)". 
 !  * Note : the diagonal entries are computed using the "Negative Sum Trick" to
 !    reduce roundoff errors. To further reduce the errors, one should compute all
 !    of the off-diagonal components first, sort them from smallest to largest, and
 !    THEN compute the diagonal components. For now, sorting is not done.
 ! 
 !  If another interpolant is provided, then the derivative matrix is created at the other set of 
 !  s. To allow for this option, the derivative matrix is ALLOCATABLE. It is important
 !  that the end user take care to release memory taken up by the derivative matrix.
 ! 
 ! =============================================================================================== !
  IMPLICIT NONE
  CLASS(Lagrange_1D), INTENT(in)          :: myPoly
  REAL(prec), ALLOCATABLE, INTENT(out)    :: dMat(:,:)
  TYPE(Lagrange_1D), INTENT(in), OPTIONAL :: otherInterpolant
  ! LOCAL
  REAL(prec) :: thisX
  REAL(prec) :: x(0:myPoly % ns), w(0:myPoly % ns)
  INTEGER    :: kP, jP, nP, mP
  REAL(prec), ALLOCATABLE :: temp(:) !offDiags(0:myPoly%ns,0:myPoly%ns), temp(0:myPoly%ns)

     CALL myPoly % GetNumberOfs( nP )
     
     IF( PRESENT(otherInterpolant) )THEN ! We will generate a derivative matrix at another set of s
     
        CALL otherInterpolant % GetNumberOfs( mP )
     
        ALLOCATE( dMat(0:mP,0:nP),  temp(0:nP) )
     
        DO kP = 0, mP ! loop over the interpolation s

           DO jP = 0, nP! loop over interpolating polynomial
           
              temp(0:nP) = ZERO
              ! "Cherry Pick" the Lagrange interpolating polynomial
              temp(jP) = ONE
           
              CALL otherInterpolant % GetNode( kP, thisX )
              dMat(kP,jP) = myPoly % EvaluateDerivative( thisX, temp ) 

           ENDDO ! loop over the interpolation s

        ENDDO ! kP, loop over the interpolation s
        
        DEALLOCATE(temp)
     
     ELSE ! We will generate the derivative matrix at the "native" s
     
        ALLOCATE( dMat(0:nP,0:nP) )
        
        CALL myPoly % Gets( x )
        CALL myPoly % GetbWs( w )

        DO kP = 0, nP ! loop over the interpolation s

           dMat(kP,kP) = ZERO
        
           DO jP = 0, nP ! loop over the interpolation s (again)
           
              IF( .NOT. (jP == kP) )THEN

                 dMat(kP,jP) = w(jP)/( w(kP)*( x(kP) - x(jP) ) )

                 dMat(kP,kP) = dMat(kP,kP) - dMat(kP,jP)
      
               ENDIF
        
           ENDDO ! loop over the interpolation s

        ENDDO ! kP, loop over the interpolation s

      ENDIF

 END SUBROUTINE CalculateDerivativeMatrix_1D
!
!
!
 SUBROUTINE CalculateInterpolationMatrix_1D( oldPoly, nNew, zNew, T )  
 ! S/R CalculateInterpolationMatrix 
 !  Description :
 !  Calculates the interpolation matrix between two sets of s for the
 !  interpolant with s "xOld(0:nP)" and barycentric bbWs "bWs(0:nP)"
 !
 ! 
 ! 
 ! =============================================================================================== !
  IMPLICIT NONE
  CLASS(Lagrange_1D), INTENT(in) :: oldPoly
  INTEGER, INTENT(in)            :: nNew
  REAL(prec), INTENT(in)         :: zNew(0:nNew)
  REAL(prec), INTENT(out)        :: T(0:nNew,0:oldPoly % ns)
  ! LOCAL
  REAL(prec) :: temp1, temp2
  REAL(prec) :: x(0:oldPoly % ns), w(0:oldPoly % ns)
  INTEGER :: kP, jP,  oldnP
  LOGICAL :: rowHasMatch 

     CALL oldPoly % GetNumberOfs( oldnP )

     CALL oldPoly % Gets( x )
     CALL oldPoly % GetbWs( w )

     DO kP = 0, nNew ! loop over the new interpolation s

        rowHasMatch = .FALSE.
       
        DO jP = 0, oldnP ! loop over the old interpolation s

           T(kP,jP) = ZERO
           
           IF( AlmostEqual( zNew(kP), x(jP) ) )THEN
      
              rowHasMatch = .TRUE.

              T(kP,jP) = ONE

           ENDIF ! IF the new and old s match

        ENDDO ! jP, loop over the old interpolation s


        IF( .NOT.(rowHasMatch) )THEN ! the interpolation s are not the same


           temp1 = ZERO

           DO jP = 0, oldnP ! loop over the old interpolation s         
              
              temp2 = w(jP)/( zNew(kP) - x(jP) )

              T(kP,jP) = temp2

              temp1 = temp1 + temp2

           ENDDO ! jP, loop over the old interpolation s

           DO jP = 0, oldnP ! loop over the old interpolation s (again)

              T(kP,jP) = T(kP,jP)/temp1

           ENDDO ! jP, loop over the old inteprolation s (again)


        ENDIF ! IF the interpolation s are not the same


     ENDDO ! kP, loop over the new interpolation s

 END SUBROUTINE CalculateInterpolationMatrix_1D
!
!
!
 SUBROUTINE CoarseToFine_1D( oldPoly, fOld, T, nNew, fNew )  
 ! S/R CoarseToFine_1D 
 !  Description :
 !   Uses the matrix "T(0:mP,0:nP)" and the function function values "fOld(0:nP)",
 !   to generate the new function values at a new specified set of mesh points.
 !   The new set of mesh points were the ones sent into S/R CalculateInterpolationMatrix to 
 !   generate the interpolation matrix. This routine is essentially a matrix
 !   multiplication routine.
 !
 !
 ! 
 ! =============================================================================================== !
  IMPLICIT NONE
  CLASS(Lagrange_1D), INTENT(in) :: oldPoly
  INTEGER, INTENT(in)            :: nNew
  REAL(prec), INTENT(in)         :: fOld(0:oldPoly%ns)
  REAL(prec), INTENT(in)         :: T(0:nNew,0:oldPoly%ns)
  REAL(prec), INTENT(out)        :: fNew(0:nNew)

     fNew = MATMUL( T, fOld )

 END SUBROUTINE CoarseToFine_1D
!
!
!==================================================================================================!
!--------------------------------- File I/O Routines ----------------------------------------------!
!==================================================================================================!
!
!
 SUBROUTINE WriteCurve_Lagrange_1D(myPoly, fAts, filename)
 ! S/R WriteCurve_Lagrange_1D
 !
 !
 ! 
 ! =============================================================================================== !
  IMPLICIT NONE
  CLASS(Lagrange_1D), INTENT(in) :: myPoly
  REAL(prec), INTENT(in)         :: fAts(0:myPoly%ns)
  CHARACTER(*), INTENT(in)       :: filename
  ! LOCAL
  INTEGER :: jS, fUnit, nP
  
     CALL myPoly % GetNumberOfs( nP )

     OPEN(UNIT=NewUnit(fUnit),&
          FILE=filename//'.curve',form='FORMATTED')
 
     WRITE(fUnit,*)'#Lagrange1dCurve'
 
     DO jS = 0, nP
        
        WRITE(fUnit,LagrangeCurveFMT) myPoly % s(jS), fAts(jS)

     ENDDO

     close(UNIT=fUnit)

     RETURN
     
 END SUBROUTINE WriteCurve_Lagrange_1D
!
!
!
END MODULE Lagrange_1D_Class
