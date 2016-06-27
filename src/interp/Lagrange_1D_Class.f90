! Lagrange_1D_Class.f90
! 
! Copyright 2015 Joseph Schoonover <schoonover.numerics@gmail.com>, The Florida State University
!
! Lagrange_1D_Class.f90 is part of the Spectral Element Libraries in Fortran (SELF).
! 
! Licensed under the Apache License, Version 2.0 (the "License"); 
! You may obtain a copy of the License at 
!
! http://www.apache.org/licenses/LICENSE-2.0 
!
! Unless required by applicable law or agreed to in writing, software 
! distributed under the License is distributed on an "AS IS" BASIS, 
! WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. 
! See the License for the specific language governing permissions and  
! limitations under the License.
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
      PROCEDURE :: Build => Build_Lagrange_1D
      PROCEDURE :: Trash => Trash_Lagrange_1D
  
      ! Accessors
      PROCEDURE :: SetNodes => SetNodes_Lagrange_1D
      PROCEDURE :: GetNodes => GetNodes_Lagrange_1D
      PROCEDURE :: SetNode => SetNode_Lagrange_1D
      PROCEDURE :: GetNode => GetNode_Lagrange_1D
      PROCEDURE :: SetWeights => SetWeights_Lagrange_1D
      PROCEDURE :: GetWeights => GetWeights_Lagrange_1D
      PROCEDURE :: SetNumberOfNodes => SetNumberOfNodes_Lagrange_1D
      PROCEDURE :: GetNumberOfNodes => GetNumberOfNodes_Lagrange_1D


      ! Data-structure operations
      PROCEDURE :: CalculateBarycentricWeights => CalculateBarycentricWeights_Lagrange_1D
      PROCEDURE :: EvaluateInterpolant => Evaluate_Lagrange_1D
      PROCEDURE :: EvaluateLagrangePolynomial => EvaluateLagrangePolynomial_Lagrange_1D
      PROCEDURE :: EvaluateDerivative => EvaluateDerivative_Lagrange_1D
      PROCEDURE :: CalculateDerivativeMatrix => CalculateDerivativeMatrix_Lagrange_1D
      PROCEDURE :: CalculateInterpolationMatrix => CalculateInterpolationMatrix_Lagrange_1D
      PROCEDURE :: CoarseToFine => CoarseToFine_Lagrange_1D
      
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
 SUBROUTINE Build_Lagrange_1D( myPoly, nS, nSo, s, so )
 ! S/R Build
 !
 !   A manual constructor for the Lagrange_1D class.
 ! 
 !   Usage :
 !      CALL myPoly % Build( nS, nSo, s, so )
 !
 !      Allocates memory and fills in data for the attributes of the Lagrange_1D class. 
 !
 !   Input : 
 !      myPoly       The Lagrange_1D data structure
 !      nS           The number of observations where we have data
 !      nSo          The number of observations that we want
 !      s(0:nS)      The node locations where we have data
 !      so(0:nSo)    The node locations where we want data
 !
 !   Output :
 !      myPoly       The Lagrange_1D data structure with its attributes filled in
 !
 ! =============================================================================================== !
 ! DECLARATIONS 
   IMPLICIT NONE
   CLASS(Lagrange_1D), INTENT(inout) :: myPoly
   INTEGER, INTENT(in)               :: nS, nSo
   REAL(prec), INTENT(in)            :: s(0:nS), so(0:nSo)
   
      ! Set the number of observations (those we have and those we want)
      myPoly % nS  = nS
      myPoly % nSo = nSo
      
      ! Allocate storage
      ALLOCATE( myPoly % s(0:nS), myPoly % bWs(0:nS) )
      ALLOCATE( myPoly % so(0:nSo), myPoly % Ts(0:nSo,0:nS) )
      
      ! Fill in the nodal locations
      myPoly % s  =  s
      myPoly % so = so

      ! and calculate the barycentric weights for quick interpolation.
      CALL myPoly % CalculateBarycentricWeights( )

      ! Using the two nodal locations, we can construct the interpolation matrix. The interpolation
      ! matrix enables quick interpolation.
      CALL myPoly % CalculateInterpolationMatrix( )
 
 END SUBROUTINE Build_Lagrange_1D
!
!
!
SUBROUTINE Trash_Lagrange_1D(myPoly)
 ! S/R Trash
 !
 !   A manual destructor for the Lagrange_1D class.
 ! 
 !   Usage :
 !      CALL myPoly % Trash( )
 !
 !      Deallocates memory for the Lagrange_1D class. 
 !
 !   Input : 
 !      myPoly       The Lagrange_1D data structure
 !
 !   Output :
 !      myPoly       An empty Lagrange_1D data structure.
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS(Lagrange_1D), INTENT(inout) :: myPoly

      DEALLOCATE( myPoly % s, myPoly % bWs )
      DEALLOCATE( myPoly % so, myPoly % Ts )

 END SUBROUTINE Trash_Lagrange_1D
!
!
!==================================================================================================!
!--------------------------------------- Accessors ------------------------------------------------!
!==================================================================================================!
!
!
 SUBROUTINE SetNodes_Lagrange_1D( myPoly, sInput )
 ! S/R SetNodes
 !  
 !   Uses "sInput" to assign the attribute "s" of the Lagrange_1D data structure. This routine is
 !   meant to be a convenience, so users do not have to make reference directly to the "s"
 !   attribute.
 !
 !   Usage :
 !      CALL myPoly % SetNodes( sInput )
 !
 !   Input :
 !      sInput       A REAL array on nodal locations where we have observations
 !
 !   Output :
 !      myPoly       The Lagrange_1D structure with the "s" attribute updated
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS(Lagrange_1D), INTENT(inout) :: myPoly
   REAL(prec), INTENT(in)            :: sInput(0:myPoly % ns)

      myPoly % s = sInput

 END SUBROUTINE SetNodes_Lagrange_1D
!
!
!
 SUBROUTINE GetNodes_Lagrange_1D( myPoly, sOutput )
 ! S/R GetNodes_Lagrange_1D
 !  
 !   Uses the Lagrange_1D data structure to report the locations where we have data ("s") in the
 !   REAL array "sOutput". This routine is meant to be a convenience, so users do not have to make
 !   reference directly to the "s" attribute.
 !
 !   Usage :
 !      CALL myPoly % GetNodes( sOutput )
 !
 !   Input :
 !      myPoly       The Lagrange_1D structure (with "s" assigned)
 !
 !   Output :
 !      sOutput      A REAL array on nodal locations where we have observations
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS(Lagrange_1D), INTENT(in) :: myPoly
   REAL(prec), INTENT(out)        :: sOutput(0:myPoly % ns)

      sOutput = myPoly % s 

 END SUBROUTINE GetNodes_Lagrange_1D
!
!
!
SUBROUTINE SetAlternateNodes_Lagrange_1D( myPoly, sInput )
 ! S/R SetAlternateNodes
 !  
 !   Uses "sInput" to assign the attribute "so" of the Lagrange_1D data structure. Recall that "so"
 !   refers to the locations where we want to have new observations, ie, it is the set of locations
 !   that we will interpolate to.
 !
 !   This routine is meant to be a convenience, so users do not have to make reference directly to 
 !   the "so" attribute.
 !
 !   Usage :
 !      CALL myPoly % SetNodes( sInput )
 !
 !   Input :
 !      sInput       A REAL array on nodal locations where we have observations
 !
 !   Output :
 !      myPoly       The Lagrange_1D structure with the "so" attribute updated
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS(Lagrange_1D), INTENT(inout) :: myPoly
   REAL(prec), INTENT(in)            :: sInput(0:myPoly % ns)

      myPoly % so = sInput

 END SUBROUTINE SetAlternateNodes_Lagrange_1D
!
!
!
 SUBROUTINE GetAlternateNodes_Lagrange_1D( myPoly, sOutput )
 ! S/R GetAlternateNodes
 !  
 !   Uses the Lagrange_1D data structure to report the locations where we want data ("so") in the
 !   REAL array "sOutput". Recall that "so" refers to the locations where we want to have new 
 !   observations, ie, it is the set of locations that we will interpolate to.
 !
 !   This routine is meant to be a convenience, so users do not have to make reference directly to 
 !   the "so" attribute.
 !
 !   Usage :
 !      CALL myPoly % GetNodes( sOutput )
 !
 !   Input :
 !      myPoly       The Lagrange_1D structure (with "so" assigned)
 !
 !   Output :
 !      sOutput      A REAL array on nodal locations where we have observations
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS(Lagrange_1D), INTENT(in) :: myPoly
   REAL(prec), INTENT(out)        :: sOutput(0:myPoly % ns)

      sOutput = myPoly % so 

 END SUBROUTINE GetAlternateNodes_Lagrange_1D
!
!
!
 SUBROUTINE SetWeights_Lagrange_1D( myPoly, wInput )
 ! S/R SetWeights
 !  
 !   Uses "wInput" to assign the attribute "bWs" of the Lagrange_1D data structure. Recall that 
 !   "bWs" refers to barycentric interpolation weights.
 !
 !   This routine is meant to be a convenience, so users do not have to make reference directly to 
 !   the "bWs" attribute.
 !
 !   Usage :
 !      CALL myPoly % SetWeights( wInput )
 !
 !   Input :
 !      wIn          A REAL array on nodal locations where we have observations
 !
 !   Output :
 !      myPoly       The Lagrange_1D structure with the "bWs" attribute updated
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS(Lagrange_1D), INTENT(inout) :: myPoly
   REAL(prec), INTENT(in)            :: wInput(0:myPoly % ns)

      myPoly % bWs = wInput

 END SUBROUTINE SetWeights_Lagrange_1D
!
!
!
 SUBROUTINE GetWeights_Lagrange_1D( myPoly, wOutput )
 ! S/R GetWeights
 !  
 !   Uses the Lagrange_1D data structure to report the barycentric interpolation weights ("bWs") in 
 !   the REAL array "wOutput". Recall that "bWs" refers to barycentric interpolation weights.
 !
 !   This routine is meant to be a convenience, so users do not have to make reference directly to 
 !   the "bWs" attribute.
 !
 !   Usage :
 !      CALL myPoly % GetWeights( wOutput )
 !
 !   Input :
 !      myPoly       The Lagrange_1D structure, with the "bWs" attribute assigned
 !
 !   Output :
 !      wOuput       A REAL array containing the barycentric weights
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS(Lagrange_1D), INTENT(in) :: myPoly
   REAL(prec), INTENT(out)        :: wOutput(0:myPoly % ns)

       wOutput = myPoly % bWs 

 END SUBROUTINE GetWeights_Lagrange_1D
!
!
!
 SUBROUTINE SetNumberOfNodes_Lagrange_1D( myPoly, N )
 ! S/R SetNumberOfNodes
 !  
 !    Sets the "nS" attribute of the Lagrange_1D data structure to N. The "nS" attribute refers
 !    to the number of nodes where we have observations.
 !  
 !
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS(Lagrange_1D), INTENT(inout) :: myPoly
   INTEGER, INTENT(in)             :: n
    
       myPoly % ns = n 

 END SUBROUTINE SetNumberOfNodesLagrange_1D
!
!
!
 SUBROUTINE GetNumberOfNodesLagrange_1D( myPoly, n )
 ! S/R GetNumberOfNodesLagrange_1D
 !  
 !    SetNodes n to the number of interpolation s (ns)
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS(Lagrange_1D), INTENT(in) :: myPoly
   INTEGER, INTENT(out)           :: n


       n = myPoly % ns 

 END SUBROUTINE GetNumberOfNodesLagrange_1D
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
 
     CALL myPoly % GetNumberOfNodes( nP )
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

     CALL myPoly % GetNumberOfNodes( nP )
     
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

     CALL myPoly % GetNumberOfNodes( nP )

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

     CALL myPoly % GetNumberOfNodes( nP )
     
     IF( PRESENT(otherInterpolant) )THEN ! We will generate a derivative matrix at another set of s
     
        CALL otherInterpolant % GetNumberOfNodes( mP )
     
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
        
        CALL myPoly % GetNodes( x )
        CALL myPoly % GetWeights( w )

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
 SUBROUTINE CalculateInterpolationMatrix_1D( oldPoly )  
 ! S/R CalculateInterpolationMatrix 
 !  Description :
 !  Calculates the interpolation matrix between two SetNodes of s for the
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

     CALL oldPoly % GetNumberOfNodes( oldnP )

     CALL oldPoly % GetNodes( x )
     CALL oldPoly % GetWeights( w )

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
  
     CALL myPoly % GetNumberOfNodes( nP )

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
