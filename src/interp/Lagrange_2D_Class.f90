! Lagrange_2D_Class.f90
! 
! Copyright 2015 Joseph Schoonover <schoonover.numerics@gmail.com>
! 
! Lagrange_2D_Class.f90 is part of the Spectral Element Libraries in Fortran (SELF).
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
 
 
MODULE Lagrange_2D_Class
! ========================================= Logs ================================================= !
!2016-05-11  Joseph Schoonover  schoonover.numerics@gmail.com 
!
! //////////////////////////////////////////////////////////////////////////////////////////////// ! 
USE ModelPrecision
USE ConstantsDictionary
USE CommonRoutines

IMPLICIT NONE




    ! Two-dimensional polynomial is the tensor
    ! product of two one-dimensional polynomials
    TYPE, PUBLIC ::  Lagrange_2D
      INTEGER                 :: nS, nP           ! Number of nodal points where we have observations
      INTEGER                 :: nSo, nPo         ! Number of nodal points where we want observations
      REAL(prec), ALLOCATABLE :: s(:), p(:)       ! Locations where we have observations
      REAL(prec), ALLOCATABLE :: bWs(:), bWp(:)   ! Barycentric weights in s and p
      REAL(prec), ALLOCATABLE :: so(:), po(:)     ! Locations where we want observations
      REAL(prec), ALLOCATABLE :: Ts(:,:), Tp(:,:) ! Interpolation matrices to get us from what we
                                                  ! have to what we want 

      CONTAINS
      !-------------!
      PROCEDURE :: Build => BuildLagrange_2D
      PROCEDURE :: Trash => TrashLagrange_2D
  
      PROCEDURE :: SetNodes => SetNodesLagrange_2D
      PROCEDURE :: GetNodes => GetNodesLagrange_2D
      PROCEDURE :: SetNode => SetNodeLagrange_2D
      PROCEDURE :: GetNode => GetNodeLagrange_2D
      PROCEDURE :: SetWeights => SetWeightsLagrange_2D
      PROCEDURE :: GetWeights => GetWeightsLagrange_2D
      PROCEDURE :: SetWeight => SetWeightLagrange_2D
      PROCEDURE :: GetWeight => GetWeightLagrange_2D
      PROCEDURE :: SetNumberOfNodes => SetNumberOfNodesLagrange_2D
      PROCEDURE :: GetNumberOfNodes => GetNumberOfNodesLagrange_2D
      
      PROCEDURE :: CalculateWeights => CalculateBarycentricWeights_2D
      PROCEDURE :: EvaluateInterpolant => EvaluateLagrange_2D
      PROCEDURE :: EvaluateDerivative => EvaluateDerivative_2D
      PROCEDURE :: CalculateDerivativeMatrix => CalculateDerivativeMatrix_2D
      PROCEDURE :: CalculateInterpolationMatrix => CalculateInterpolationMatrix_2D
      PROCEDURE :: CoarseToFine => CoarseToFine_2D
      PROCEDURE :: WriteTecplot => WriteTecplot_Lagrange_2D
      
    END TYPE Lagrange_2D




 CONTAINS
!
!
!==================================================================================================!
!------------------------------- Manual Constructors/Destructors ----------------------------------!
!==================================================================================================!
!
!
 SUBROUTINE BuildLagrange_2D( myPoly, nS, nP, nSo, nPo, xIn, yIn )
 ! S/R BuildLagrange_2D
 !
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS(Lagrange_2D), INTENT(inout) :: myPoly
   INTEGER, INTENT(in)               :: nS, nP, nSo, nPo
   REAL(prec), INTENT(in), OPTIONAL  :: xIn(0:nS), yIn(0:nP)

      myPoly % nS  = nS
      myPoly % nP  = nP
      myPoly % nSo = nSo
      myPoly % nPo = nPo
      
      ALLOCATE( myPoly % s(0:nS), myPoly % p(0:nS) )
      myPoly % s = ZERO
      myPoly % p = ZERO

      ALLOCATE( myPoly % bWs(0:nS), myPoly % bWp(0:nP) )

      ALLOCATE( myPoly % so(0:nS), myPoly % po(0:nS) )
      myPoly % so = ZERO
      myPoly % po = ZERO
      
      ! I know that, in this allocate statement, Ts is indexed over the new s-points then the old
      ! s-points and Tp is indexed with the old points first instead. This is because the routine
      ! to map the solution from the points we have to the points we want requires that the second
      ! interpolation matrix is transposed.
      ALLOCATE( myPoly % Ts(0:nSo,0:nS), myPoly % Tp(0:nP,0:nPo) )  
      myPoly % Ts = ZERO
      myPoly % Tp = ZERO

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
         


 END SUBROUTINE BuildLagrange_2D
!
!
!
SUBROUTINE TrashLagrange_2D(myPoly)
 ! S/R TrashLagrange_2D
 ! 
 !    DEALLOCATES memory occupied by myPoly
 !
 ! 
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS(Lagrange_2D), INTENT(inout) :: myPoly


      CALL myPoly % sInterp % Trash( )
      CALL myPoly % pInterp % Trash( )

 END SUBROUTINE TrashLagrange_2D
!
!
!==================================================================================================!
!--------------------------------------- Accessors ------------------------------------------------!
!==================================================================================================!
!
!
 SUBROUTINE SetNodesLagrange_2D( myPoly, xIn, yIn )
 ! S/R SetNodesLagrange_2D
 !  
 !
 !
 ! =============================================================================================== !
 ! DECLARATIONS
    IMPLICIT NONE
    CLASS(Lagrange_2D), INTENT(inout) :: myPoly
    REAL(prec), INTENT(in), OPTIONAL  :: xIn(0:myPoly % nS), yIn(0:myPoly % nP)

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

 END SUBROUTINE SetNodesLagrange_2D
!
!
!
 SUBROUTINE GetNodesLagrange_2D( myPoly, xOut, yOut )
 ! S/R GetNodesLagrange_2D
 !  
 !
 ! =============================================================================================== !
 ! DECLARATIONS
    IMPLICIT NONE
    CLASS(Lagrange_2D), INTENT(in) :: myPoly
    REAL(prec), INTENT(out)        :: xOut(0:myPoly % nS), yOut(0:myPoly % nP)


       CALL myPoly % sInterp % GetNodes(xOut)
       CALL myPoly % pInterp % GetNodes(yOut)    

 END SUBROUTINE GetNodesLagrange_2D
!
!
!
 SUBROUTINE SetNodeLagrange_2D( myPoly, iX, iY, xIn, yIn)
 ! S/R SetNodeLagrange_2D
 !  
 !
 ! =============================================================================================== !
 ! DECLARATIONS
    IMPLICIT NONE
    CLASS(Lagrange_2D), INTENT(inout) :: myPoly
    INTEGER, INTENT(in)               :: iX, iY
    REAL(prec), INTENT(out)           :: xIn, yIn


       CALL myPoly % sInterp % SetNode(iX, xIn)
       CALL myPoly % pInterp % SetNode(iY, yIn)    

 END SUBROUTINE SetNodeLagrange_2D
!
!
!
 SUBROUTINE GetNodeLagrange_2D( myPoly, iX, iY, xOut, yOut)
 ! S/R GetNodeLagrange_2D
 !  
 !
 ! =============================================================================================== !
 ! DECLARATIONS
    IMPLICIT NONE
    CLASS(Lagrange_2D), INTENT(in) :: myPoly
    INTEGER, INTENT(in)            :: iX, iY
    REAL(prec), INTENT(out)        :: xOut, yOut


       CALL myPoly % sInterp % GetNode(iX, xOut)
       CALL myPoly % pInterp % GetNode(iY, yOut)    

 END SUBROUTINE GetNodeLagrange_2D
!
!
!
 SUBROUTINE SetWeightsLagrange_2D( myPoly, wsIn, wpIn )
 ! S/R SetWeightsLagrange_2D
 !  
 !
 !
 ! =============================================================================================== !
 ! DECLARATIONS
    IMPLICIT NONE
    CLASS(Lagrange_2D), INTENT(inout) :: myPoly
    REAL(prec), INTENT(in), OPTIONAL  :: wsIn(0:myPoly % nS), wpIn(0:myPoly % nP)

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

 END SUBROUTINE SetWeightsLagrange_2D
!
!
!
 SUBROUTINE GetWeightsLagrange_2D( myPoly, wsOut, wpOut )
 ! S/R GetWeightsLagrange_2D
 !  
 !
 ! =============================================================================================== !
 ! DECLARATIONS
    IMPLICIT NONE
    CLASS(Lagrange_2D), INTENT(in) :: myPoly
    REAL(prec), INTENT(out)        :: wsOut(0:myPoly % nS), wpOut(0:myPoly % nP)


       CALL myPoly % sInterp % GetWeights(wsOut)
       CALL myPoly % pInterp % GetWeights(wpOut)    

 END SUBROUTINE GetWeightsLagrange_2D
!
!
!
 SUBROUTINE SetWeightLagrange_2D( myPoly, iX, iY, wsIn, wpIn )
 ! S/R SetWeightLagrange_2D
 !  
 !
 ! =============================================================================================== !
 ! DECLARATIONS
    IMPLICIT NONE
    CLASS(Lagrange_2D), INTENT(inout) :: myPoly
    INTEGER, INTENT(in)               :: iX, iY
    REAL(prec), INTENT(out)           :: wsIn, wpIn


       CALL myPoly % sInterp % SetWeight(iX, wsIn)
       CALL myPoly % pInterp % SetWeight(iY, wpIn)    

 END SUBROUTINE SetWeightLagrange_2D
!
!
!
 SUBROUTINE GetWeightLagrange_2D( myPoly, iX, iY, wsOut, wpOut )
 ! S/R GetWeightLagrange_2D
 !  
 !
 ! =============================================================================================== !
 ! DECLARATIONS
    IMPLICIT NONE
    CLASS(Lagrange_2D), INTENT(in) :: myPoly
    INTEGER, INTENT(in)            :: iX, iY
    REAL(prec), INTENT(out)        :: wsOut, wpOut


       CALL myPoly % sInterp % GetWeight(iX, wsOut)
       CALL myPoly % pInterp % GetWeight(iY, wpOut)    

 END SUBROUTINE GetWeightLagrange_2D
!
!
!
 SUBROUTINE SetNumberOfNodesLagrange_2D( myPoly, nS, nP )
 ! S/R SetNumberOfNodesLagrange_2D
 !  
 !
 ! =============================================================================================== !
 ! DECLARATIONS
    IMPLICIT NONE
    CLASS(Lagrange_2D), INTENT(inout) :: myPoly
    INTEGER, INTENT(in)               :: nS, nP

       myPoly % nS = nS
       myPoly % nP = nP

 END SUBROUTINE SetNumberOfNodesLagrange_2D
!
!
!
 SUBROUTINE GetNumberOfNodesLagrange_2D( myPoly, nS, nP )
 ! S/R GetNumberOfNodesLagrange_2D
 !  
 !
 ! =============================================================================================== !
 ! DECLARATIONS
    IMPLICIT NONE
    CLASS(Lagrange_2D), INTENT(in) :: myPoly
    INTEGER, INTENT(out)           :: nS, nP

       nS = myPoly % nS
       nP = myPoly % nP

 END SUBROUTINE GetNumberOfNodesLagrange_2D
!
!
!==================================================================================================!
!--------------------------------- Type Specific Routines -----------------------------------------!
!==================================================================================================!
!
!
 SUBROUTINE CalculateBarycentricWeights_2D(myPoly)
 ! S/R CalculateBarycentricWeights_2D
 !  
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS(Lagrange_2D), INTENT(inout) :: myPoly

   
      CALL myPoly % sInterp % CalculateWeights( )
      CALL myPoly % pInterp % CalculateWeights( )

 END SUBROUTINE CalculateBarycentricWeights_2D
!
!
!
 FUNCTION EvaluateLagrange_2D(myPoly, xEval, yEval, f) RESULT (inFatXY)  
 ! FUNCTION EvaluateLagrange_2D
 !  
 !
 ! =============================================================================================== !
  IMPLICIT NONE
  CLASS(Lagrange_2D), INTENT(in) :: myPoly
  REAL(prec), INTENT(in)         :: xEval, yEval
  REAL(prec), INTENT(in)         :: f(0:myPoly % nS, 0:myPoly % nP)
  REAL(prec)                     :: inFatXY 
  ! LOCAL
  REAL(prec)  :: lS(0:myPoly % nS)
  REAL(prec)  :: lP(0:myPoly % nP)
  INTEGER     :: jS, jP, nS, nP
 
     CALL myPoly % GetNumberOfNodes( nS, nP )
     
     lS = myPoly % sInterp % EvaluateLagrangePolynomial( xEval ) 
     lP = myPoly % pInterp % EvaluateLagrangePolynomial( yEval )

     inFatXY = ZERO

     DO jS = 0, nS ! Loop over the s-nodes
        DO jP = 0, nP ! Loop over the p-nodes 
  
           inFatXY = inFatXY + f(jS,jP)*lS(jS)*lP(jP)

        ENDDO ! jS, loop over the p-nodes
     ENDDO ! jP, loop over the s-nodes

  

 END FUNCTION EvaluateLagrange_2D
!
!
!
 FUNCTION EvaluateDerivative_2D( myPoly, xEval, yEval, f ) RESULT(gradF)  
 ! FUNCTION EvaluateDerivative_2D 
 !  
 !
 ! =============================================================================================== !
 ! DECLARATIONS
  IMPLICIT NONE
  CLASS(Lagrange_2D), INTENT(in)  :: myPoly
  REAL(prec), INTENT(in)          :: xEval, yEval
  REAL(prec), INTENT(in)          :: f(0:myPoly % nS, 0:myPoly % nP)
  REAL(prec)                      :: gradF(1:2)
  ! LOCAL
  REAL(prec)                      :: dfds, dfdp
  REAL(prec)  :: lS(0:myPoly % nS), f1(0:myPoly % nS)
  REAL(prec)  :: lP(0:myPoly % nP), f2(0:myPoly % nP)
  INTEGER     :: jS, jP, nS, nP
 
     CALL myPoly % GetNumberOfNodes( nS, nP )
     
     lS = myPoly % sInterp % EvaluateLagrangePolynomial(xEval) 
     lP = myPoly % pInterp % EvaluateLagrangePolynomial(yEval)

     dfds = ZERO
     dfdp = ZERO

     DO jP = 0, nP ! Loop over p-points

        f1 = f(0:nS,jP)
        ! Evaluate s-derivative
        dfds = dfds + ( myPoly % sInterp % EvaluateDerivative(xEval, f1) )*lP(jP)

     ENDDO

     DO jS = 0, nS ! Loop over s-points

        f2 = f(jS,0:nP)
        ! Evaluate p-derivative
        dfdp = dfdp + ( myPoly % pInterp % EvaluateDerivative(yEval, f2) )*lS(jS)

     ENDDO
     
     gradf(1) = dfds
     gradf(2) = dfdp

 END FUNCTION EvaluateDerivative_2D
!
!
!
 SUBROUTINE CalculateDerivativeMatrix_2D( myPoly, dMatX, dMatY, otherInterpolant )  
 ! S/R CalculateDerivativeMatrix 
 !  
 !
 ! =============================================================================================== !
  IMPLICIT NONE
  CLASS(Lagrange_2D), INTENT(in)          :: myPoly
  REAL(prec), INTENT(out), ALLOCATABLE    :: dMatX(:,:)
  REAL(prec), INTENT(out), ALLOCATABLE    :: dMatY(:,:)
  TYPE(Lagrange_2D), INTENT(in), OPTIONAL :: otherInterpolant
  
     IF( PRESENT(otherInterpolant) )THEN
        CALL myPoly % sInterp % CalculateDerivativeMatrix( dMatX, otherInterpolant % sInterp )
        CALL myPoly % pInterp % CalculateDerivativeMatrix( dMatY, otherInterpolant % pInterp )
     ELSE
        CALL myPoly % sInterp % CalculateDerivativeMatrix( dMatX )
        CALL myPoly % pInterp % CalculateDerivativeMatrix( dMatY )
     ENDIF

 END SUBROUTINE CalculateDerivativeMatrix_2D
!
!
!
 SUBROUTINE CalculateInterpolationMatrix_2D( oldPoly, nSnew, nPnew, sNew, pNew, Ts, Tp)  
 ! S/R CalculateInterpolationMatrix
 !  
 !
 ! =============================================================================================== !
  IMPLICIT NONE
  CLASS(Lagrange_2D), INTENT(in) :: oldPoly
  INTEGER, INTENT(in)            :: nSnew, nPnew
  REAL(prec), INTENT(in)         :: sNew(0:nSnew)
  REAL(prec), INTENT(in)         :: pNew(0:nPnew)
  REAL(prec), INTENT(out)        :: Ts(0:nSnew,0:oldPoly % nS)
  REAL(prec), INTENT(out)        :: Tp(0:,0:oldPoly % nP)


     CALL oldPoly % sInterp % CalculateInterpolationMatrix( nSnew, sNew, Ts )  
     CALL oldPoly % pInterp % CalculateInterpolationMatrix( nPnew, pNew, Tp )  

 END SUBROUTINE CalculateInterpolationMatrix_2D
!
!
!
 SUBROUTINE CoarseToFine_2D( oldPoly, fOld, Ts, Tp, nSnew, nPnew, fNew) 
 ! S/R CoarseToFine_2D 
 !  
 !
 ! =============================================================================================== !
  IMPLICIT NONE
  CLASS(Lagrange_2D), INTENT(in)  :: oldPoly
  INTEGER, INTENT(in)             :: nSnew, nPnew
  REAL(prec), INTENT(in)          :: Ts(0:nSnew,0:oldPoly % nS)
  REAL(prec), INTENT(in)          :: Tp(0:nPnew,0:oldPoly % nP)
  REAL(prec), INTENT(in)          :: fOld(0:oldPoly % nS, 0:oldPoly % nP)
  REAL(prec), INTENT(out)         :: fNew(0:nSnew, 0:nPnew)

  ! LOCAL
  REAL(prec) :: fInt(0:nSnew,0:oldPoly % nP) 
  INTEGER :: j, nSold, nPold

      CALL oldPoly % GetNumberOfNodes( nSold, nPold )
       
      ! Interpolation in 2-D can be broken into two Matrix-Matrix multiplies,
      ! To interpolate from the points (s_i, p_j) to (s_a, p_b), we use the interpolation formula,
      !
      !                    F_{a,b} = \sum_{j} [ \sum_{i} [F_{i,j} l_i(s_a)] l_j(s_b)]            (1)
      !
      ! l_i(s_a) is the a-th row of the i-th column of the "Interpolation Matrix" (Ts). In 2-D,
      ! we can treat the 2-D array of nodal values as a 2-D matrix. The inner sum over i in the
      ! Eq. (1) is matrix-matrix multiply. Let M denote the number of points in the new
      ! grid, and N the number of points in the old grid. (Ts)_{a,i} = l_i(s_a) is an MxN matrix,
      ! and F_{i,j} is NxN. The matrix-matrix multiply
      !
      !                           F*_{a,j} = Ts*F                                                (2)
      !
      ! Yields an MxN set of nodal values that are interpolated onto the new nodes in the 
      ! s-direction, but not in the p-direction. The outer sum in Eq. (1) then can be written
      !
      !        F_{a,b} = \sum_{j}[ F*_{a,j}l_j(s_b) ] = \sum_{j}[F*_{a,j}Tp_{b,j}]               (3)
      !  
      ! Which amounts to a matrix-matrix multiply with the transpose of the interpolation matrix.
      ! This is why, when the interpolation matrices are made, the transpose of the p-Interpolation 
      ! is stored. This illustrates why the algorithm here looks rather simple, with two 
      ! matrix-matrix multiplies.
  
      fInt = MATMUL( Ts, fOld )
      fNew = MATMUL( Tp, fIntT )


  END SUBROUTINE CoarseToFine_2D
!
!
!==================================================================================================!
!--------------------------------- File I/O Routines ----------------------------------------------!
!==================================================================================================!
!
!
 SUBROUTINE WriteTecplot_Lagrange_2D( myPoly, func, funcname, filename )
! S/R WriteTecplot_Lagrange_2D
 !  Description :
 !  
 !
 ! =============================================================================================== !
  IMPLICIT NONE
  CLASS(Lagrange_2D), INTENT(in) :: myPoly
  REAL(prec), INTENT(in)         :: func(0:myPoly % nS, 0:myPoly % nP)
  CHARACTER(*), INTENT(in)       :: funcname, filename  
  ! Local
  INTEGER    :: jS, jP, nS, nP, fUnit
  REAL(prec) :: s, p
  
    CALL myPoly % GetNumberOfNodes( nS, nP )
    
    OPEN( UNIT=NewUnit(fUnit), &
          FILE= trim(filename)//'.tec',&
          FORM='formatted' )

    WRITE(fUnit,*) 'VARIABLES = "X", "Y", "'//trim(funcname)//'"'
    WRITE(fUnit,*) 'ZONE I=',nS+1,', J=', nP+1,',DATAPACKING=POINT'

    DO jP = 0, nP
       DO jS = 0, nS
       
          CALL myPoly % GetNode( jS, jP, s, p )

          WRITE (fUnit,*) s, p, func(jS,jP)

       ENDDO
    ENDDO

    CLOSE(UNIT=fUnit)

 END SUBROUTINE WriteTecplot_Lagrange_2D
!
!
!
END MODULE Lagrange_2D_Class
