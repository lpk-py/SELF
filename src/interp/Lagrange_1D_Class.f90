MODULE Lagrange_1D_Class
! Lagrange_1D_Class.f90
! 
! Module History 
! 
! o (v 1.0 - 21 January 2014 )
!
! o (v 2.1 - 17 November 2015)
!
! =======================================================================================
USE ModelPrecision
USE ConstantsDictionary
USE CommonRoutines

IMPLICIT NONE



   TYPE, PUBLIC :: Lagrange_1D
      INTEGER, PRIVATE                 :: nNodes     ! number of nodes
      REAL(prec), PRIVATE, ALLOCATABLE :: nodes(:)   ! Interpolation nodes
      REAL(prec), PRIVATE, ALLOCATABLE :: weights(:) ! barycentric weights

      CONTAINS
      
      !-------------!
      ! Constructors/Destructors
      PROCEDURE :: Build => BuildLagrange_1D
      PROCEDURE :: Trash => TrashLagrange_1D
  
  
      ! Accessors
      PROCEDURE :: SetNodes => SetNodesLagrange_1D
      PROCEDURE :: GetNodes => GetNodesLagrange_1D
      PROCEDURE :: SetNode => SetNodeLagrange_1D
      PROCEDURE :: GetNode => GetNodeLagrange_1D
      PROCEDURE :: SetWeights => SetWeightsLagrange_1D
      PROCEDURE :: GetWeights => GetWeightsLagrange_1D
      PROCEDURE :: SetWeight => SetWeightLagrange_1D
      PROCEDURE :: GetWeight => GetWeightLagrange_1D
      PROCEDURE :: SetNumberOfNodes => SetNumberOfNodesLagrange_1D
      PROCEDURE :: GetNumberOfNodes => GetNumberOfNodesLagrange_1D


      ! Data-structure operations
      PROCEDURE :: CalculateWeights => CalculateBarycentricWeights_1D
      PROCEDURE :: EvaluateInterpolant => EvaluateLagrange_1D
      PROCEDURE :: EvaluateLagrangePolynomial => EvaluateLagrangePolynomial_1D
      PROCEDURE :: EvaluateDerivative => EvaluateDerivative_1D
      PROCEDURE :: CalculateDerivativeMatrix => CalculateDerivativeMatrix_1D
      PROCEDURE :: CalculateInterpolationMatrix => CalculateInterpolationMatrix_1D
      PROCEDURE :: CoarseToFine => CoarseToFine_1D
      
      ! File I/O Routines
      PROCEDURE :: WriteCurve => WriteCurve_Lagrange_1D

    END TYPE Lagrange_1D


 INTEGER, PARAMETER, PUBLIC       :: nDefaultNodes = 1          ! Default number of nodes
 REAL(PREC),PARAMETER, PUBLIC     :: LagrangeNodeDefault = ZERO ! The default value to set the nodes
 CHARACTER(17), PARAMETER, PUBLIC :: LagrangeCurveFMT = '(E17.7,2x,E17.7)'
 
 CONTAINS
!
!
!==================================================================================================!
!------------------------------- Manual Constructors/Destructors ----------------------------------!
!==================================================================================================!
!
!
 SUBROUTINE BuildLagrange_1D( myPoly, nNodes, xIn )
 ! S/R BuildLagrange_1D
 !
 !   This subroutine is the basic constructor. If the number of nodes are input, THEN the desired 
 !   amount of memory is allocated. If the actual interpolation nodes are input, THEN the nodes
 !   are stored and the barycentric weights are calculated. If the number of nodes are not given,
 !   THEN the number of nodes are determined from the size of xIn. In the event that xIn is also 
 !   not given, the default number of nodes is used (nDefaultNodes) and the nodes and barycentric
 !   weights are set to zero.
 !
 !
 ! =============================================================================================== !
 ! DECLARATIONS 
   IMPLICIT NONE
   CLASS(Lagrange_1D), INTENT(inout) :: myPoly
   INTEGER, INTENT(in)               :: nNodes
   REAL(prec), INTENT(in),optional   :: xIn(:)
   

      ! Set the number of nodes
      CALL myPoly % SetNumberOfNodes( nNodes )
      
      ! Allocate storage
      ALLOCATE( myPoly % nodes(0:nNodes), myPoly % weights(0:nNodes) )
      
      IF( PRESENT(xIn) )THEN
         CALL myPoly % SetNodes( xIn )
         CALL myPoly % CalculateWeights( )
      ELSE
         CALL myPoly % SetNodes( )
         CALL myPoly % SetWeights( )
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


      myPoly % nNodes = 0

      DEALLOCATE( myPoly % nodes, myPoly % weights )

     

 END SUBROUTINE TrashLagrange_1D
!
!
!==================================================================================================!
!--------------------------------------- Accessors ------------------------------------------------!
!==================================================================================================!
!
!
 SUBROUTINE SetNodesLagrange_1D( myPoly, xIn )
 ! S/R SetNodesLagrange_1D
 !  
 !    Sets the attribute "nodes" to xIn. If xIn is not present, THEN the nodes are set to the 
 !    default value (LagrangeNodeDefault)
 !
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS(Lagrange_1D), INTENT(inout) :: myPoly
   REAL(prec), INTENT(in), OPTIONAL  :: xIn(0:myPoly % nNodes)

       IF( PRESENT(xIn) )THEN
          myPoly % nodes = xIn
       ELSE
          myPoly % nodes = LagrangeNodeDefault
       ENDIF

 END SUBROUTINE SetNodesLagrange_1D
!
!
!
 SUBROUTINE GetNodesLagrange_1D( myPoly, xOut )
 ! S/R GetNodesLagrange_1D
 !  
 !    Gets the attribute "nodes" and sets xOut.
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS(Lagrange_1D), INTENT(in) :: myPoly
   REAL(prec), INTENT(out)        :: xOut(0:myPoly % nNodes)


       xOut = myPoly % nodes 

 END SUBROUTINE GetNodesLagrange_1D
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


      myPoly % nodes(i) = xIn 

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


      xOut = myPoly % nodes(i) 

 END SUBROUTINE GetNodeLagrange_1D
!
!
!
 SUBROUTINE SetWeightsLagrange_1D( myPoly, wIn )
 ! S/R SetWeightsLagrange_1D
 !  
 !    Sets the attribute "weights" to wIn. If wIn is not present, THEN the barycentric weights are
 !    set to the default value (LagrangeNodeDefault)
 !
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS(Lagrange_1D), INTENT(inout) :: myPoly
   REAL(prec), INTENT(in), OPTIONAL  :: wIn(0:myPoly % nNodes)

       IF( PRESENT(wIn) )THEN
          myPoly % weights = wIn
       ELSE
          myPoly % weights = LagrangeNodeDefault
       ENDIF

 END SUBROUTINE SetWeightsLagrange_1D
!
!
!
 SUBROUTINE GetWeightsLagrange_1D( myPoly, wOut )
 ! S/R GetWeightsLagrange_1D
 !  
 !    Gets the attribute "weights" and sets wOut.
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS(Lagrange_1D), INTENT(in) :: myPoly
   REAL(prec), INTENT(out)        :: wOut(0:myPoly % nNodes)


       wOut = myPoly % weights 

 END SUBROUTINE GetWeightsLagrange_1D
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


       myPoly % weights(i) = wIn 

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


       wOut = myPoly % weights(i) 

 END SUBROUTINE GetWeightLagrange_1D
!
!
!
 SUBROUTINE SetNumberOfNodesLagrange_1D( myPoly, n )
 ! S/R SetNumberOfNodesLagrange_1D(
 !  
 !    Sets the number of nodes (nNodes) to n.
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS(Lagrange_1D), INTENT(inout) :: myPoly
   INTEGER, INTENT(in)             :: n
    
       myPoly % nNodes = n 

 END SUBROUTINE SetNumberOfNodesLagrange_1D
!
!
!
 SUBROUTINE GetNumberOfNodesLagrange_1D( myPoly, n )
 ! S/R GetNumberOfNodesLagrange_1D
 !  
 !    Sets n to the number of interpolation nodes (nNodes)
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS(Lagrange_1D), INTENT(in) :: myPoly
   INTEGER, INTENT(out)           :: n


       n = myPoly % nNodes 

 END SUBROUTINE GetNumberOfNodesLagrange_1D
!
!
!==================================================================================================!
!--------------------------------- Type Specific Routines -----------------------------------------!
!==================================================================================================!
!
!
 SUBROUTINE CalculateBarycentricWeights_1D( myPoly )
 ! S/R CalculateBarycentricWeights_1D
 !  
 !    Calculates the barycentric weights from the interpolation nodes x(0:nP) and stores them in 
 !    the "weights" attribute.
 !    
 !
! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS(Lagrange_1D), INTENT(inout) :: myPoly
   ! LOCAL
   INTEGER :: jP, kP, nP
   
      CALL myPoly % GetNumberOfNodes( nP )

      myPoly % weights(0:nP) = ONE ! initialize the weights to 1

      ! Computes the product w_k = w_k*(x_k - x_j), k /= j
      DO jP = 1,nP ! loop over the interpolation nodes

         DO kP = 0, jP-1 ! loop to perform multiplication for weights

            myPoly % weights(kP) = myPoly % weights(kP)*&
                                  ( myPoly % nodes(kP) - myPoly % nodes(jP) )

            myPoly % weights(jP) = myPoly % weights(jP)*&
                                 ( myPoly % nodes(jP) - myPoly % nodes(kP))

         ENDDO ! kP, mulitplication loop

      ENDDO ! jP, loop over the interpolation nodes


     DO jP = 0, nP
 
        myPoly % weights(jP) = ONE/myPoly % weights(jP)
 
     ENDDO ! jP
    
  
     RETURN

 END SUBROUTINE CalculateBarycentricWeights_1D
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
   REAL(prec), INTENT(in)         :: f(0:myPoly%nNodes)
   REAL(prec)                     :: inFatX 
   ! LOCAL
   REAL(prec) :: num, den, t, thisX, thisW
   INTEGER    :: jP, nP 
 
     CALL myPoly % GetNumberOfNodes( nP )
     num = ZERO
     den = ZERO

     DO jP = 0, nP ! loop over the interpolation nodes

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
        

     ENDDO ! jP, loop over the interpolation nodes


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
  REAL(prec)                     :: lAtX(0:myPoly%nNodes)
  ! LOCAL
  REAL(prec) :: temp1, temp2, thisX, thisW
  INTEGER    :: jP, nP
  LOGICAL    :: xMatchesNode

     CALL myPoly % GetNumberOfNodes( nP )
     
     xMatchesNode = .FALSE.

     DO jP = 0, nP ! loop over the interpolation nodes

        CALL myPoly % GetNode( jP, thisX )
        
        lAtX(jP) = ZERO

        IF( AlmostEqual(xEval, thisX) ) THEN !xEval is an interpolation node

          lAtX(jP) = ONE
          
          xMatchesNode = .TRUE.
  
        ENDIF 

     ENDDO ! jP, loop over the interpolation nodes

     
     IF( xMatchesNode )THEN ! we're done, the evaluation point is one of the nodes

        RETURN

     ENDIF
     
     ! Otherwise, we calculate the values for the lagrange polynomials

     temp1 = ZERO
     
     DO jP = 0, nP ! loop over the interpolation nodes

        CALL myPoly % GetWeight( jP, thisW )
        CALL myPoly % GetNode( jP, thisX )

        temp2 = thisW/(xEval - thisX)

        lAtX(jP) = temp2

        temp1 = temp1 + temp2

     ENDDO ! jP, loop over the interpolation nodes
     
     DO jP = 0, nP ! loop over the interpolation nodes (again)

        lAtX(jP) = lAtX(jP)/temp1

     ENDDO ! jP, loop over the interpolation nodes (again)
     

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
  REAL(prec), INTENT(in)         :: f(0:myPoly%nNodes)
  REAL(prec)                     :: dInFdx
  ! LOCAL
  REAL(prec) :: num, den, t, p, thisX, thisW
  INTEGER    :: jP, jS, nP 
  LOGICAL    :: atNode 

     CALL myPoly % GetNumberOfNodes( nP )

     num = ZERO
     atNode = .FALSE.

     DO jP = 0, nP ! loop over the interpolation nodes

        CALL myPoly % GetNode( jP, thisX )

        IF( AlmostEqual(xEval, thisX) ) THEN !xEval is an interpolation node
     
           atNode = .TRUE.

           p = f(jP)
           
           CALL myPoly % GetWeight( jP, thisW )
           
           den = -thisW
   
           jS = jP

        ENDIF  ! IF the evaluation point is an interpolation node
        
     ENDDO ! jP, loop over the interpolation nodes



     IF( atNode ) THEN ! the evaluation location is an interpolating node


        DO jP = 0, nP ! loop over the interpolation nodes

           IF( .NOT.(jP == jS) ) THEN ! x(jP) is not the interpolating node

             CALL myPoly % GetNode( jP, thisX )
             CALL myPoly % GetWeight( jP, thisW )
             
             num = num + thisW*(p - f(jP))/(xEval - thisX)

           ENDIF

        ENDDO ! jP, loop over the interpolation nodes


     ELSE ! the evaluation location is not an interpolating node


        den = ZERO

        p = EvaluateLagrange_1D(myPoly, xEval, f) 
        
        DO jP = 0, nP ! loop over the interpolation nodes

           CALL myPoly % GetNode( jP, thisX )
           CALL myPoly % GetWeight( jP, thisW )

           t = thisW/(xEval - thisX)

           num = num + t*(p - f(jP))/(xEval - thisX)

           den = den + t

        ENDDO ! jP, loop over the interpolation nodes


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
 !  at the interpolation nodes "myPoly%nodes(:)". 
 !  * Note : the diagonal entries are computed using the "Negative Sum Trick" to
 !    reduce roundoff errors. To further reduce the errors, one should compute all
 !    of the off-diagonal components first, sort them from smallest to largest, and
 !    THEN compute the diagonal components. For now, sorting is not done.
 ! 
 !  If another interpolant is provided, then the derivative matrix is created at the other set of 
 !  nodes. To allow for this option, the derivative matrix is ALLOCATABLE. It is important
 !  that the end user take care to release memory taken up by the derivative matrix.
 ! 
 ! =============================================================================================== !
  IMPLICIT NONE
  CLASS(Lagrange_1D), INTENT(in)          :: myPoly
  REAL(prec), ALLOCATABLE, INTENT(out)    :: dMat(:,:)
  TYPE(Lagrange_1D), INTENT(in), OPTIONAL :: otherInterpolant
  ! LOCAL
  REAL(prec) :: temp1, temp2, thisW, thatW, thisX, thatX
  REAL(prec) :: x(0:myPoly % nNodes), w(0:myPoly % nNodes)
  INTEGER    :: kP, jP, nP, mP
  REAL(prec), ALLOCATABLE :: temp(:) !offDiags(0:myPoly%nNodes,0:myPoly%nNodes), temp(0:myPoly%nNodes)

     CALL myPoly % GetNumberOfNodes( nP )
     
     IF( PRESENT(otherInterpolant) )THEN ! We will generate a derivative matrix at another set of nodes
     
        CALL otherInterpolant % GetNumberOfNodes( mP )
     
        ALLOCATE( dMat(0:mP,0:nP),  temp(0:nP) )
     
        DO kP = 0, mP ! loop over the interpolation nodes

           DO jP = 0, nP! loop over interpolating polynomial
           
              temp(0:nP) = ZERO
              ! "Cherry Pick" the Lagrange interpolating polynomial
              temp(jP) = ONE
           
              CALL otherInterpolant % GetNode( kP, thisX )
              dMat(kP,jP) = myPoly % EvaluateDerivative( thisX, temp ) 

           ENDDO ! loop over the interpolation nodes

        ENDDO ! kP, loop over the interpolation nodes
        
        DEALLOCATE(temp)
     
     ELSE ! We will generate the derivative matrix at the "native" nodes
     
        ALLOCATE( dMat(0:nP,0:nP) )
        
        CALL myPoly % GetNodes( x )
        CALL myPoly % GetWeights( w )

        DO kP = 0, nP ! loop over the interpolation nodes

           dMat(kP,kP) = ZERO
        
           DO jP = 0, nP ! loop over the interpolation nodes (again)
           
              IF( .NOT. (jP == kP) )THEN

                 dMat(kP,jP) = w(jP)/( w(kP)*( x(kP) - x(jP) ) )

                 dMat(kP,kP) = dMat(kP,kP) - dMat(kP,jP)
      
               ENDIF
        
           ENDDO ! loop over the interpolation nodes

        ENDDO ! kP, loop over the interpolation nodes

      ENDIF

 END SUBROUTINE CalculateDerivativeMatrix_1D
!
!
!
 SUBROUTINE CalculateInterpolationMatrix_1D( oldPoly, nNew, zNew, T )  
 ! S/R CalculateInterpolationMatrix 
 !  Description :
 !  Calculates the interpolation matrix between two sets of nodes for the
 !  interpolant with nodes "xOld(0:nP)" and barycentric bweights "weights(0:nP)"
 !
 ! 
 ! 
 ! =============================================================================================== !
  IMPLICIT NONE
  CLASS(Lagrange_1D), INTENT(in) :: oldPoly
  INTEGER, INTENT(in)            :: nNew
  REAL(prec), INTENT(in)         :: zNew(0:nNew)
  REAL(prec), INTENT(out)        :: T(0:nNew,0:oldPoly % nNodes)
  ! LOCAL
  REAL(prec) :: temp1, temp2, newX, oldX, newW, oldW
  REAL(prec) :: x(0:oldPoly % nNodes), w(0:oldPoly % nNodes)
  INTEGER :: kP, jP,  oldnP
  LOGICAL :: rowHasMatch 

     CALL oldPoly % GetNumberOfNodes( oldnP )

     CALL oldPoly % GetNodes( x )
     CALL oldPoly % GetWeights( w )

     DO kP = 0, nNew ! loop over the new interpolation nodes

        rowHasMatch = .FALSE.
       
        DO jP = 0, oldnP ! loop over the old interpolation nodes

           T(kP,jP) = ZERO
           
           IF( AlmostEqual( zNew(kP), x(jP) ) )THEN
      
              rowHasMatch = .TRUE.

              T(kP,jP) = ONE

           ENDIF ! IF the new and old nodes match

        ENDDO ! jP, loop over the old interpolation nodes


        IF( .NOT.(rowHasMatch) )THEN ! the interpolation nodes are not the same


           temp1 = ZERO

           DO jP = 0, oldnP ! loop over the old interpolation nodes         
              
              temp2 = w(jP)/( zNew(kP) - x(jP) )

              T(kP,jP) = temp2

              temp1 = temp1 + temp2

           ENDDO ! jP, loop over the old interpolation nodes

           DO jP = 0, oldnP ! loop over the old interpolation nodes (again)

              T(kP,jP) = T(kP,jP)/temp1

           ENDDO ! jP, loop over the old inteprolation nodes (again)


        ENDIF ! IF the interpolation nodes are not the same


     ENDDO ! kP, loop over the new interpolation nodes

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
  REAL(prec), INTENT(in)         :: fOld(0:oldPoly%nNodes)
  REAL(prec), INTENT(in)         :: T(0:nNew,0:oldPoly%nNodes)
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
 SUBROUTINE WriteCurve_Lagrange_1D(myPoly, fAtNodes, filename)
 ! S/R WriteCurve_Lagrange_1D
 !
 !
 ! 
 ! =============================================================================================== !
  IMPLICIT NONE
  CLASS(Lagrange_1D), INTENT(in) :: myPoly
  REAL(prec), INTENT(in)         :: fAtNodes(0:myPoly%nNodes)
  CHARACTER(*), INTENT(in)       :: filename
  ! LOCAL
  INTEGER :: jS, fUnit, nP
  
     CALL myPoly % GetNumberOfNodes( nP )

     OPEN(UNIT=NewUnit(fUnit),&
          FILE=filename//'.curve',form='FORMATTED')
 
     WRITE(fUnit,*)'#Lagrange1dCurve'
 
     DO jS = 0, nP
        
        WRITE(fUnit,LagrangeCurveFMT) myPoly % nodes(jS), fAtNodes(jS)

     ENDDO

     close(UNIT=fUnit)

     RETURN
     
 END SUBROUTINE WriteCurve_Lagrange_1D
!
!
!
END MODULE Lagrange_1D_Class
