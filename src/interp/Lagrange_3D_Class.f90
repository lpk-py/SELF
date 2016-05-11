! Lagrange_3D_Class.f90
! 
! Copyright 2015 Joseph Schoonover <schoonover.numerics@gmail.com>
! 
! Lagrange_3D_Class.f90 is part of the Spectral Element Libraries in Fortran (SELF).
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
 
 
MODULE Lagrange_3D_Class
! ========================================= Logs ================================================= !
!2016-05-11  Joseph Schoonover  schoonover.numerics@gmail.com 
!
! //////////////////////////////////////////////////////////////////////////////////////////////// ! 

USE ModelPrecision
USE ConstantsDictionary
USE CommonRoutines
USE Lagrange_1D_Class

IMPLICIT NONE




    ! Two-dimensional polynomial is the tensor
    ! product of two one-dimensional polynomials
    TYPE, PUBLIC :: Lagrange_3D
      INTEGER            :: nS, nP, nQ
      TYPE(Lagrange_1D)  :: sInterp
      TYPE(Lagrange_1D)  :: pInterp
      TYPE(Lagrange_1D)  :: qInterp

            CONTAINS
      !-------------!
      PROCEDURE :: Build => BuildLagrange_3D
      PROCEDURE :: Trash => TrashLagrange_3D
  
      PROCEDURE :: SetNodes => SetNodesLagrange_3D
      PROCEDURE :: GetNodes => GetNodesLagrange_3D
      PROCEDURE :: SetNode => SetNodeLagrange_3D
      PROCEDURE :: GetNode => GetNodeLagrange_3D
      PROCEDURE :: SetWeights => SetWeightsLagrange_3D
      PROCEDURE :: GetWeights => GetWeightsLagrange_3D
      PROCEDURE :: SetWeight => SetWeightLagrange_3D
      PROCEDURE :: GetWeight => GetWeightLagrange_3D
      PROCEDURE :: SetNumberOfNodes => SetNumberOfNodesLagrange_3D
      PROCEDURE :: GetNumberOfNodes => GetNumberOfNodesLagrange_3D

      PROCEDURE :: CalculateWeights => CalculateBarycentricWeights_3D
      PROCEDURE :: EvaluateInterpolant => EvaluateLagrange_3D
      PROCEDURE :: EvaluateDerivative => EvaluateDerivative_3D
      PROCEDURE :: CalculateDerivativeMatrix => CalculateDerivativeMatrix_3D
      PROCEDURE :: CalculateInterpolationMatrix => CalculateInterpolationMatrix_3D
      PROCEDURE :: CoarseToFine => CoarseToFine_3D
      PROCEDURE :: WriteTecplot => WriteTecplot_Lagrange_3D
      
    END TYPE Lagrange_3D

 
 CONTAINS
!
!
!==================================================================================================!
!------------------------------- Manual Constructors/Destructors ----------------------------------!
!==================================================================================================!
!
!
 SUBROUTINE BuildLagrange_3D( myPoly, nS, nP, nQ, xIn, yIn, zIn )
 ! S/R BuildLagrange_3D
 !
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS(Lagrange_3D), INTENT(inout) :: myPoly
   INTEGER, INTENT(in)               :: nS, nP, nQ
   REAL(prec), INTENT(in), OPTIONAL  :: xIn(0:nS), yIn(0:nP), zIn(0:nQ)

      CALL myPoly % SetNumberOfNodes( nS, nP, nQ )
      
      IF( PRESENT(xIn) )THEN
         CALL myPoly % sInterp % Build( nS, xIn )
      ELSE
         CALL myPoly % sInterp % Build( nS )
      ENDIF
      
      IF( PRESENT(yIn) )THEN
         CALL myPoly % pInterp % Build( nP, yIn )
      ELSE
         CALL myPoly % pInterp % Build( nP )
      ENDIF
      
      IF( PRESENT(zIn) )THEN
         CALL myPoly % qInterp % Build( nQ, zIn )
      ELSE
         CALL myPoly % qInterp % Build( nQ )
      ENDIF
         


 END SUBROUTINE BuildLagrange_3D
!
!
!
 SUBROUTINE TrashLagrange_3D( myPoly )
 ! S/R TrashLagrange_3D
 ! 
 !    DEALLOCATES memory occupied by myPoly
 !
 ! 
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   CLASS(Lagrange_3D), INTENT(inout) :: myPoly


      CALL myPoly % sInterp % Trash( ) 
      CALL myPoly % pInterp % Trash( )
      CALL myPoly % qInterp % Trash( )

 END SUBROUTINE TrashLagrange_3D
!
!
!==================================================================================================!
!--------------------------------------- Accessors ------------------------------------------------!
!==================================================================================================!
!
!
SUBROUTINE SetNodesLagrange_3D( myPoly, xIn, yIn, zIn )
 ! S/R SetNodesLagrange_3D
 !  
 !
 !
 ! =============================================================================================== !
 ! DECLARATIONS
    IMPLICIT NONE
    CLASS(Lagrange_3D), INTENT(inout) :: myPoly
    REAL(prec), INTENT(in), OPTIONAL  :: xIn(0:myPoly % nS), yIn(0:myPoly % nP), zIn(0:myPoly % nQ)

       IF( PRESENT(xIn) )THEN
          CALL myPoly % sInterp % SetNodes(xIn)
       ELSE
          CALL myPoly % sInterp % SetNodes( )
       ENDIF
       
       IF( PRESENT(yIn) )THEN
          CALL myPoly % pInterp % SetNodes(yIn)
       ELSE
          CALL myPoly % pInterp % SetNodes( )
       ENDIF 
       
       IF( PRESENT(zIn) )THEN
          CALL myPoly % qInterp % SetNodes(zIn)
       ELSE
          CALL myPoly % qInterp % SetNodes( )
       ENDIF 

 END SUBROUTINE SetNodesLagrange_3D
!
!
!
 SUBROUTINE GetNodesLagrange_3D( myPoly, xOut, yOut, zOut )
 ! S/R GetNodesLagrange_3D
 !  
 !
 ! =============================================================================================== !
 ! DECLARATIONS
    IMPLICIT NONE
    CLASS(Lagrange_3D), INTENT(in) :: myPoly
    REAL(prec), INTENT(out)        :: xOut(0:myPoly % nS), yOut(0:myPoly % nP), zOut(0:myPoly % nQ)


       CALL myPoly % sInterp % GetNodes(xOut)
       CALL myPoly % pInterp % GetNodes(yOut)
       CALL myPoly % qInterp % GetNodes(zOut)    

 END SUBROUTINE GetNodesLagrange_3D
!
!
!
 SUBROUTINE SetNodeLagrange_3D( myPoly, iX, iY, iZ, xIn, yIn, zIn )
 ! S/R SetNodeLagrange_3D
 !  
 !
 ! =============================================================================================== !
 ! DECLARATIONS
    IMPLICIT NONE
    CLASS(Lagrange_3D), INTENT(inout) :: myPoly
    INTEGER, INTENT(in)               :: iX, iY, iZ
    REAL(prec), INTENT(out)           :: xIn, yIn, zIn
    

       CALL myPoly % sInterp % SetNode(iX, xIn)
       CALL myPoly % pInterp % SetNode(iY, yIn)
       CALL myPoly % qInterp % SetNode(iZ, zIn)    

 END SUBROUTINE SetNodeLagrange_3D
!
!
!
 SUBROUTINE GetNodeLagrange_3D( myPoly, iX, iY, iZ, xOut, yOut, zOut )
 ! S/R GetNodeLagrange_3D
 !  
 !
 ! =============================================================================================== !
 ! DECLARATIONS
    IMPLICIT NONE
    CLASS(Lagrange_3D), INTENT(in) :: myPoly
    INTEGER, INTENT(in)            :: iX, iY, iZ
    REAL(prec), INTENT(out)        :: xOut, yOut, zOut


       CALL myPoly % sInterp % GetNode(iX, xOut)
       CALL myPoly % pInterp % GetNode(iY, yOut)
       CALL myPoly % qInterp % GetNode(iZ, zOut)    

 END SUBROUTINE GetNodeLagrange_3D
!
!
!
 SUBROUTINE SetWeightsLagrange_3D( myPoly, wsIn, wpIn, wqIn )
 ! S/R SetWeightsLagrange_3D
 !  
 !
 !
 ! =============================================================================================== !
 ! DECLARATIONS
    IMPLICIT NONE
    CLASS(Lagrange_3D), INTENT(inout) :: myPoly
    REAL(prec), INTENT(in), OPTIONAL  :: wsIn(0:myPoly % nS), wpIn(0:myPoly % nP), wqIn(0:myPoly % nQ)

       IF( PRESENT(wsIn) )THEN
          CALL myPoly % sInterp % SetWeights(wsIn)
       ELSE
          CALL myPoly % sInterp % SetWeights( )
       ENDIF
       
       IF( PRESENT(wpIn) )THEN
          CALL myPoly % pInterp % SetWeights(wpIn)
       ELSE
          CALL myPoly % pInterp % SetWeights( )
       ENDIF 
       
       IF( PRESENT(wqIn) )THEN
          CALL myPoly % qInterp % SetWeights(wqIn)
       ELSE
          CALL myPoly % qInterp % SetWeights( )
       ENDIF 

 END SUBROUTINE SetWeightsLagrange_3D
!
!
!
 SUBROUTINE GetWeightsLagrange_3D( myPoly, wsOut, wpOut, wqOut )
 ! S/R GetWeightsLagrange_3D
 !  
 !
 ! =============================================================================================== !
 ! DECLARATIONS
    IMPLICIT NONE
    CLASS(Lagrange_3D), INTENT(in) :: myPoly
    REAL(prec), INTENT(out)        :: wsOut(0:myPoly % nS), wpOut(0:myPoly % nP), wqOut(0:myPoly % nQ)


       CALL myPoly % sInterp % GetWeights(wsOut)
       CALL myPoly % pInterp % GetWeights(wpOut) 
       CALL myPoly % qInterp % GetWeights(wqOut)    

 END SUBROUTINE GetWeightsLagrange_3D
!
!
!
 SUBROUTINE SetWeightLagrange_3D( myPoly, iX, iY,  iZ, wsIn, wpIn, wqIn )
 ! S/R SetWeightLagrange_3D
 !  
 !
 ! =============================================================================================== !
 ! DECLARATIONS
    IMPLICIT NONE
    CLASS(Lagrange_3D), INTENT(inout) :: myPoly
    INTEGER, INTENT(in)               :: iX, iY, iZ
    REAL(prec), INTENT(out)           :: wsIn, wpIn, wqIn


       CALL myPoly % sInterp % SetWeight(iX, wsIn)
       CALL myPoly % pInterp % SetWeight(iY, wpIn)  
       CALL myPoly % qInterp % SetWeight(iZ, wqIn)    

 END SUBROUTINE SetWeightLagrange_3D
!
!
!
 SUBROUTINE GetWeightLagrange_3D( myPoly, iX, iY, iZ, wsOut, wpOut, wqOut )
 ! S/R GetWeightLagrange_3D
 !  
 !
 ! =============================================================================================== !
 ! DECLARATIONS
    IMPLICIT NONE
    CLASS(Lagrange_3D), INTENT(in) :: myPoly
    INTEGER, INTENT(in)            :: iX, iY, iZ
    REAL(prec), INTENT(out)        :: wsOut, wpOut, wqOut


       CALL myPoly % sInterp % GetWeight(iX, wsOut)
       CALL myPoly % pInterp % GetWeight(iY, wpOut)  
       CALL myPoly % qInterp % GetWeight(iZ, wqOut)    

 END SUBROUTINE GetWeightLagrange_3D
!
!
!
 SUBROUTINE SetNumberOfNodesLagrange_3D( myPoly, nS, nP, nQ )
 ! S/R SetNumberOfNodesLagrange_3D
 !  
 !
 ! =============================================================================================== !
 ! DECLARATIONS
    IMPLICIT NONE
    CLASS(Lagrange_3D), INTENT(inout) :: myPoly
    INTEGER, INTENT(in)               :: nS, nP, nQ

       myPoly % nS = nS
       myPoly % nP = nP
       myPoly % nQ = nQ

 END SUBROUTINE SetNumberOfNodesLagrange_3D
!
!
!
 SUBROUTINE GetNumberOfNodesLagrange_3D( myPoly, nS, nP, nQ )
 ! S/R GetNumberOfNodesLagrange_3D
 !  
 !
 ! =============================================================================================== !
 ! DECLARATIONS
    IMPLICIT NONE
    CLASS(Lagrange_3D), INTENT(in) :: myPoly
    INTEGER, INTENT(out)           :: nS, nP, nQ

       nS = myPoly % nS
       nP = myPoly % nP
       nQ = myPoly % nQ

 END SUBROUTINE GetNumberOfNodesLagrange_3D
!
!
!==================================================================================================!
!--------------------------------- Type Specific Routines -----------------------------------------!
!==================================================================================================!
!
!
 SUBROUTINE CalculateBarycentricWeights_3D( myPoly )
 ! S/R CalculateBarycentricWeights_3D
 !  
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS(Lagrange_3D), INTENT(inout) :: myPoly
   
      CALL myPoly % sInterp % CalculateWeights()
      CALL myPoly % pInterp % CalculateWeights()
      CALL myPoly % qInterp % CalculateWeights()

 END SUBROUTINE CalculateBarycentricWeights_3D
!
!
!
 FUNCTION EvaluateLagrange_3D( myPoly, xEval, yEval, zEval, f ) RESULT( inFatXY )  
 ! FUNCTION EvaluateLagrange_3D
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS(Lagrange_3D), INTENT(in) :: myPoly
   REAL(prec), INTENT(in)         :: xEval, yEval, zEval
   REAL(prec), INTENT(in)         :: f(0:myPoly % nS, 0:myPoly % nP, 0:myPoly % nQ)
   REAL(prec)                     :: inFatXY 
 ! LOCAL
   REAL(prec)  :: lS(0:myPoly % nS)
   REAL(prec)  :: lP(0:myPoly % nP)
   REAL(prec)  :: lQ(0:myPoly % nQ)
   INTEGER     :: jS, jP, jQ, nS, nP, nQ
 
     CALL myPoly % GetNumberOfNodes( nS, nP, nQ )
     
     lS = myPoly % sInterp % EvaluateLagrangePolynomial( xEval ) 
     lP = myPoly % pInterp % EvaluateLagrangePolynomial( yEval ) 
     lQ = myPoly % qInterp % EvaluateLagrangePolynomial( zEval ) 

     inFatXY = ZERO

     DO jS = 0,myPoly % nS ! Loop over the s-nodes
        DO jP = 0,myPoly % nP ! Loop over the p-nodes 
           DO jQ = 0,myPoly % nQ ! Loop over the q-nodes
  
              inFatXY = inFatXY + f(jS,jP,jQ)*lS(jS)*lP(jP)*lQ(jQ)

           ENDDO ! jQ, loop over the q-nodes
        ENDDO ! jP, loop over the p-nodes
     ENDDO ! jS, loop over the s-nodes

  

 END FUNCTION EvaluateLagrange_3D
!
!
!
 FUNCTION EvaluateDerivative_3D(myPoly, xEval, yEval, zEval, f) RESULT(gradF)  
 ! FUNCTION EvaluateDerivative_3D 
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS(Lagrange_3D), INTENT(in) :: myPoly
   REAL(prec), INTENT(in)         :: xEval, yEval, zEval
   REAL(prec), INTENT(in)         :: f(0:myPoly % nS, 0:myPoly % nP, 0:myPoly % nQ)
   REAL(prec)                     :: gradF(1:3)
 ! LOCAL
   REAL(prec)  :: lS(0:myPoly % nS)
   REAL(prec)  :: lP(0:myPoly % nP)
   REAL(prec)  :: lQ(0:myPoly % nQ)
   REAL(prec)  :: thisf(0:myPoly % nS), thatf(0:myPoly % nP), theotherf(0:myPoly % nQ)
   INTEGER     :: jS, jP, jQ
   INTEGER     :: nS, nP, nQ
 
     CALL myPoly % GetNumberOfNodes( nS, nP, nQ )
     
     lS = myPoly % sInterp % EvaluateLagrangePolynomial( xEval ) 
     lP = myPoly % pInterp % EvaluateLagrangePolynomial( yEval ) 
     lQ = myPoly % qInterp % EvaluateLagrangePolynomial( zEval ) 


     gradF(1:3) = ZERO

     DO jP = 0, nP ! Loop over p-points

        DO jQ = 0, nQ ! Loop over q-points

           ! Evaluate s-derivative
           thisf = f(0:nS,jP,jQ)
           gradF(1) = gradF(1) + myPoly % sInterp % EvaluateDerivative(xEval, thisf)*lP(jP)*lQ(jQ)

        ENDDO

     ENDDO

   
     DO jS = 0, nS ! Loop over s-points

        DO jQ = 0, nQ ! Loop over q-points

           ! Evaluate p-derivative
           thatf = f(jS,0:nP,jQ)
           gradF(2) = gradF(2) + myPoly % pInterp % EvaluateDerivative(yEval, thatf)*lS(jS)*lQ(jQ)

        ENDDO

     ENDDO


     DO jS = 0, nS ! Loop over s-points

        DO jP = 0, nP ! Loop over p-points

           ! Evaluate q-derivative
           theotherf = f(jS,jP,0:nQ)
           gradF(3) = gradF(3) + myPoly % qInterp % EvaluateDerivative(zEval, theotherf)*lS(jS)*lP(jP)

        ENDDO

     ENDDO


 END FUNCTION EvaluateDerivative_3D
!
!
!
 SUBROUTINE CalculateDerivativeMatrix_3D( myPoly, dMatX, dMatY, dMatZ, otherInterpolant )  
 ! S/R CalculateDerivativeMatrix (POLYnomial Derivative MATRIX)
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS(Lagrange_3D), INTENT(in)          :: myPoly
   REAL(prec), INTENT(out), ALLOCATABLE    :: dMatX(:,:)
   REAL(prec), INTENT(out), ALLOCATABLE    :: dMatY(:,:)
   REAL(prec), INTENT(out), ALLOCATABLE    :: dMatZ(:,:)
   TYPE(Lagrange_3D), INTENT(in), OPTIONAL :: otherInterpolant

     IF( PRESENT(otherInterpolant) )THEN
        CALL myPoly % sInterp % CalculateDerivativeMatrix( dMatX, otherInterpolant % sInterp )
        CALL myPoly % pInterp % CalculateDerivativeMatrix( dMatY, otherInterpolant % pInterp )
        CALL myPoly % qInterp % CalculateDerivativeMatrix( dMatZ, otherInterpolant % qInterp )
     ELSE
        CALL myPoly % sInterp % CalculateDerivativeMatrix( dMatX )
        CALL myPoly % pInterp % CalculateDerivativeMatrix( dMatY )
        CALL myPoly % qInterp % CalculateDerivativeMatrix( dMatZ )
     ENDIF

 END SUBROUTINE CalculateDerivativeMatrix_3D
!
!
!
 SUBROUTINE CalculateInterpolationMatrix_3D( oldPoly, nSnew, nPnew, nQnew, sNew, pNew, qNew, Ts, Tp, Tq )  
 ! S/R CalculateInterpolationMatrix
 !  
 !
 ! =============================================================================================== !
  IMPLICIT NONE
  CLASS(Lagrange_3D), INTENT(in) :: oldPoly
  INTEGER, INTENT(in)            :: nSnew, nPnew, nQnew
  REAL(prec), INTENT(in)         :: sNew(0:nSnew)
  REAL(prec), INTENT(in)         :: pNew(0:nPnew)
  REAL(prec), INTENT(in)         :: qNew(0:nQnew)
  REAL(prec), INTENT(out)        :: Ts(0:,0:)
  REAL(prec), INTENT(out)        :: Tp(0:,0:)
  REAL(prec), INTENT(out)        :: Tq(0:,0:)

     CALL oldPoly % sInterp % CalculateInterpolationMatrix( nSnew, sNew, Ts )  
     CALL oldPoly % pInterp % CalculateInterpolationMatrix( nPnew, pNew, Tp )  
     CALL oldPoly % qInterp % CalculateInterpolationMatrix( nQnew, qNew, Tq )  

 END SUBROUTINE CalculateInterpolationMatrix_3D
!
!
!
 SUBROUTINE CoarseToFine_3D( oldPoly, fOld, Ts, Tp, Tq, nSnew, nPnew, nQnew, fNew ) 
 ! S/R CoarseToFine_3D 
 !  
 !
 ! =============================================================================================== !
  IMPLICIT NONE
  CLASS(Lagrange_3D), INTENT(in) :: oldPoly
  INTEGER, INTENT(in)            :: nSnew, nPnew, nQnew
  REAL(prec), INTENT(in)         :: Ts(0:nSnew,0:oldPoly % nS)
  REAL(prec), INTENT(in)         :: Tp(0:nPnew,0:oldPoly % nP)
  REAL(prec), INTENT(in)         :: Tq(0:nQnew,0:oldPoly % nQ)
  REAL(prec), INTENT(in)         :: fOld(0:oldPoly % nS, 0:oldPoly % nP, 0:oldPoly % nQ)
  REAL(prec), INTENT(out)        :: fNew(0:nSnew, 0:nPnew, 0:nQnew)

  ! LOCAL
  REAL(prec) :: fInt1(0:nSnew, 0:oldPoly % nP, 0:oldPoly % nQ) 
  REAL(prec) :: fInt2(0:nSnew, 0:nPnew, 0:oldPoly % nQ)
  INTEGER :: jS, jP, jQ
  INTEGER :: nSold, nPold, nQold
  

    CALL oldPoly % GetNumberOfNodes( nSold, nPold, nQold )
   

    DO jP  = 0, nPold ! Loop over the old p-points and

       DO jQ = 0, nQold ! Loop over the old q-points
       
       ! interpolate onto the new s-points
          CALL oldPoly % sInterp % CoarseToFine(fOld(0:nSold,jP, jQ), Ts, nSnew, fInt1(0:nSnew,jP, jQ) )    

       ENDDO ! jQ, loop over the old  q-points

    ENDDO ! jP, loop over the old p-points


    DO jS  = 0, nSnew ! Loop over the new s-points and

       DO jQ = 0, nQold ! Loop over the old q-points

       ! interpolate onto the new p-points
       CALL oldPoly % pInterp % CoarseToFine(fInt1(jS,0:nPold, jQ), Tp, nPnew, fInt2(jS,0:nPnew, jQ) )  
      

       ENDDO ! jQ, loop over the old s-points

    ENDDO ! jS, loop over the old s-points

    DO jS  = 0, nSnew ! Loop over the new s-points and

       DO jP = 0, nPnew ! Loop over the new p-points

       ! interpolate onto the new q-points
         CALL oldPoly % qInterp % CoarseToFine(fInt2(jS,jP, 0:nQold), Tq, nQnew, fNew(jS,jP,0:nQnew) )

       ENDDO ! jP, loop over the old p-points

    ENDDO ! jS, loop over the old s-points

  END SUBROUTINE CoarseToFine_3D
!
!
!==================================================================================================!
!--------------------------------- File I/O Routines ----------------------------------------------!
!==================================================================================================!
!
!
 SUBROUTINE WriteTecplot_Lagrange_3D( myPoly, func, funcname, filename )
! S/R WriteTecplot_Lagrange_3D
 !  Description :
 !  
 !
 ! =============================================================================================== !
  IMPLICIT NONE
  CLASS(Lagrange_3D), INTENT(in) :: myPoly
  REAL(prec), INTENT(in)         :: func(0:myPoly % nS, 0:myPoly % nP, 0:myPoly % nQ)
  CHARACTER(*), INTENT(in)       :: funcname, filename  
  ! Local
  INTEGER    :: jS, jP, jQ, nS, nP, nQ, fUnit
  REAL(prec) :: s, p, q
  
    CALL myPoly % GetNumberOfNodes( nS, nP, nQ )
    
    OPEN( UNIT=NewUnit(fUnit), &
          FILE= trim(filename)//'.tec',&
          FORM='formatted' )

    WRITE(fUnit,*) 'VARIABLES = "X", "Y", "Z", "'//trim(funcname)//'"'
    WRITE(fUnit,*) 'ZONE I=',nS+1,', J=', nP+1,', K=', nQ+1,',DATAPACKING=POINT'

    DO jQ = 0, nQ
       DO jP = 0, nP
          DO jS = 0, nS
       
             CALL myPoly % GetNode( jS, jP, jQ, s, p, q )

             WRITE (fUnit,*) s, p, q, func(jS,jP,jQ)

          ENDDO
       ENDDO
    ENDDO

    CLOSE(UNIT=fUnit)

 END SUBROUTINE WriteTecplot_Lagrange_3D
!
!
!
END MODULE Lagrange_3D_Class
