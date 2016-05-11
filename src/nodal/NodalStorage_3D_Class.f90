! NodalStorage_3D_Class.f90
! 
! Copyright 2015 Joseph Schoonover <schoonover.numerics@gmail.com>
! 
! NodalStorage_3D_Class.f90 is part of the Spectral Element Libraries in Fortran (SELF).
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
 

MODULE NodalStorage_3D_Class
! ========================================= Logs ================================================= !
!2016-05-11  Joseph Schoonover  schoonover.numerics@gmail.com 
!
! //////////////////////////////////////////////////////////////////////////////////////////////// ! 
USE ConstantsDictionary
USE ModelPrecision
USE ModelFlags
USE Legendre
USE Lagrange_1D_Class
USE Lagrange_3D_Class

IMPLICIT NONE


    TYPE NodalStorage_3D
      INTEGER                          :: nS, nP, nQ
      TYPE(Lagrange_3D)                :: interp       ! Interpolant
      REAL(prec), ALLOCATABLE, PRIVATE :: qWeightS(:)  ! dimensioned (0:nS) : Quadrature weights in the "S" direction
      REAL(prec), ALLOCATABLE, PRIVATE :: qWeightP(:)  ! dimensioned (0:nP) : Quadrature weights in the "P" direction
      REAL(prec), ALLOCATABLE, PRIVATE :: qWeightQ(:)  ! dimensioned (0:nP) : Quadrature weights in the "P" direction
      REAL(prec), ALLOCATABLE, PRIVATE :: dMatS(:,:)   ! dimensioned (0:nS, 0:nS) : Contains the s-derivative matrix
      REAL(prec), ALLOCATABLE, PRIVATE :: dMatP(:,:)   ! dimensioned (0:nP, 0:nP) : Contains the p-derivative matrix
      REAL(prec), ALLOCATABLE, PRIVATE :: dMatQ(:,:)   ! dimensioned (0:nP, 0:nP) : Contains the p-derivative matrix
      REAL(prec), ALLOCATABLE, PRIVATE :: lagWest(:)   ! s-lagrange polynomials at x=-1, dimensioned (0:nS)
      REAL(prec), ALLOCATABLE, PRIVATE :: lagEast(:)   ! s-lagrange polynomials at x=1, dimensioned (0:nS)
      REAL(prec), ALLOCATABLE, PRIVATE :: lagSouth(:)  ! p-lagrange polynomials at x=-1, dimensioned (0:nP)
      REAL(prec), ALLOCATABLE, PRIVATE :: lagNorth(:)  ! p-lagrange polynomials at x=1, dimensioned (0:nP)
      REAL(prec), ALLOCATABLE, PRIVATE :: lagBottom(:) ! q-lagrange polynomials at x=-1, dimensioned (0:nP)
      REAL(prec), ALLOCATABLE, PRIVATE :: lagTop(:)    ! q-lagrange polynomials at x=1, dimensioned (0:nP)

      CONTAINS

      ! Manual Constructors/Destructors
      PROCEDURE :: Build => BuildNodalStorage_3D
      PROCEDURE :: Trash => TrashNodalStorage_3D

      ! Accessors
      PROCEDURE :: SetNumberOfNodes => SetNumberOfNodes_NodalStorage_3D
      PROCEDURE :: GetNumberOfNodes => GetNumberOfNodes_NodalStorage_3D
      PROCEDURE :: SetNodes => SetNodes_NodalStorage_3D
      PROCEDURE :: GetNodes => GetNodes_NodalStorage_3D
      PROCEDURE :: SetQuadratureWeights => SetQuadratureWeights_NodalStorage_3D
      PROCEDURE :: GetQuadratureWeights => GetQuadratureWeights_NodalStorage_3D
      PROCEDURE :: GetDerivativeMatrix => GetDerivativeMatrix_NodalStorage_3D
      PROCEDURE :: SetEasternInterpolants => SetEasternInterpolants_NodalStorage_3D
      PROCEDURE :: GetEasternInterpolants => GetEasternInterpolants_NodalStorage_3D
      PROCEDURE :: SetWesternInterpolants => SetWesternInterpolants_NodalStorage_3D
      PROCEDURE :: GetWesternInterpolants => GetWesternInterpolants_NodalStorage_3D
      PROCEDURE :: SetSouthernInterpolants => SetSouthernInterpolants_NodalStorage_3D
      PROCEDURE :: GetSouthernInterpolants => GetSouthernInterpolants_NodalStorage_3D
      PROCEDURE :: SetNorthernInterpolants => SetNorthernInterpolants_NodalStorage_3D
      PROCEDURE :: GetNorthernInterpolants => GetNorthernInterpolants_NodalStorage_3D
      PROCEDURE :: SetBottomInterpolants => SetBottomInterpolants_NodalStorage_3D
      PROCEDURE :: GetBottomInterpolants => GetBottomInterpolants_NodalStorage_3D
      PROCEDURE :: SetTopInterpolants => SetTopInterpolants_NodalStorage_3D
      PROCEDURE :: GetTopInterpolants => GetTopInterpolants_NodalStorage_3D

      ! Type-Specific
      PROCEDURE :: GenerateMappingMatrix => GenerateMappingMatrix_NodalStorage_3D
      PROCEDURE :: ApplyDerivativeMatrix => ApplyDerivativeMatrix_NodalStorage_3D
      PROCEDURE :: CalculateAtBoundaries => CalculateAtBoundaries_NodalStorage_3D
      !PROCEDURE :: Integrate => Integrate_NodalStorage_3D
      
      
    END TYPE NodalStorage_3D


! Default parameters
 INTEGER, PRIVATE :: defaultN          = 1
 INTEGER, PRIVATE :: defaultQuadrature = GAUSS_LOBATTO
 INTEGER, PRIVATE :: defaultForm       = CG


 CONTAINS
!
!
!==================================================================================================!
!------------------------------- Manual Constructors/Destructors ----------------------------------!
!==================================================================================================!
!
!
 SUBROUTINE BuildNodalStorage_3D( myNodal, N, M, P, quadrature, approxForm  )
 ! S/R BuildNodalStorage_3D
 !  
 !   This constructor allocates space for the NodalStorage_3D class and fills in the quadrature and
 !   derivative matrix. The type of quadrature can be set to GAUSS or GAUSS_LOBATTO; this determines
 !   the nodes that are generated and stored. The approximation form can be CG for Continuous 
 !   Galerkin or DG  for Discontinuous Galerkin. The number of nodes (minus 1) is N. These parameters
 !   are all optional. Therefore, to call this function with the optional parameters, use the 
 !   syntax :
 !
 !   CALL myNodal % Build( N          = numberOfPoints, &
 !                         M          = numberOfPoints, &
 !                         P          = numberOfPoints, &
 !                         quadrature = yourQuadrature, &
 !                         approxForm = yourForm          )
 !
 !   If none of the optional arguments are given, the Build defaults to a linear interpolant through
 !   s=-1 and s=1 and the Continuous Galerkin approximation is used giving the standard derivative
 !   matrix.
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS(NodalStorage_3D), INTENT(out) :: myNodal
   INTEGER, INTENT(in), optional       :: N, M, P, quadrature, approxForm
   !LOCAL
   INTEGER                 :: iP, jP, thisQuadrature, thisN, thisM, thisP, thisForm
   REAL(prec), ALLOCATABLE :: tempS(:), tempWs(:), tempDmatS(:,:)
   REAL(prec), ALLOCATABLE :: tempP(:), tempWp(:), tempDmatP(:,:)
   REAL(prec), ALLOCATABLE :: tempQ(:), tempWq(:), tempDmatQ(:,:)

   ! Check and setup the optional arguments
      IF( PRESENT(N) )THEN
         thisN = N
      ELSE
         thisN = defaultN
      ENDIF
      
      IF( PRESENT(M) )THEN
         thisM = M
      ELSE
         thisM = thisN
      ENDIF
      
      IF( PRESENT(P) )THEN
         thisP = P
      ELSE
         thisP = thisN
      ENDIF
      
      IF( PRESENT(quadrature) )THEN
         thisQuadrature = quadrature
      ELSE
         thisQuadrature = defaultQuadrature
      ENDIF
      
      IF( PRESENT(approxForm) )THEN
         thisForm = approxForm
      ELSE
         thisForm = defaultForm
      ENDIF
      
      ! Allocate space
      CALL myNodal % SetNumberOfNodes( thisN )
      ALLOCATE( tempS(0:thisN), tempWs(0:thisN), tempDmatS(0:thisN,0:thisN) )
      ALLOCATE( tempP(0:thisN), tempWp(0:thisM), tempDmatP(0:thisM,0:thisM) )
      ALLOCATE( tempQ(0:thisP), tempWq(0:thisP), tempDmatQ(0:thisP,0:thisP) )

      ALLOCATE( myNodal % lagWest(0:thisN), myNodal % lagEast(0:thisN) )
      ALLOCATE( myNodal % lagSouth(0:thisM), myNodal % lagNorth(0:thisM) )
      ALLOCATE( myNodal % lagBottom(0:thisP), myNodal % lagTop(0:thisP) )

      ALLOCATE( myNodal % dMatS(0:thisN,0:thisN) )
      ALLOCATE( myNodal % dMatP(0:thisM,0:thisM) )
      ALLOCATE( myNodal % dMatQ(0:thisP,0:thisP) )

      ALLOCATE( myNodal % qWeightS(0:thisN) )
      ALLOCATE( myNodal % qWeightP(0:thisM) )
      ALLOCATE( myNodal % qWeightQ(0:thisP) )

      ! Generate the quadrature
      CALL GenerateLegendreQuadrature( thisN, tempS, tempWs, thisQuadrature )
      CALL GenerateLegendreQuadrature( thisM, tempP, tempWp, thisQuadrature )
      CALL GenerateLegendreQuadrature( thisP, tempQ, tempWq, thisQuadrature )
      
      ! Store the quadrature weights
      CALL myNodal % SetQuadratureWeights( tempWs, tempWp, tempWq )

      ! Build and store the interpolant
      CALL myNodal % interp % Build( thisN, thisM, thisP, tempS, tempP, tempQ )
   
      ! Calculate and store the interpolants evaluated at the endpoints
      myNodal % lagWest = myNodal % interp % sInterp % EvaluateLagrangePolynomial( -ONE )

      myNodal % lagEast = myNodal % interp % sInterp % EvaluateLagrangePolynomial( ONE )
      
      myNodal % lagSouth = myNodal % interp % pInterp % EvaluateLagrangePolynomial( -ONE )

      myNodal % lagNorth = myNodal % interp % pInterp % EvaluateLagrangePolynomial( ONE )
      
      myNodal % lagBottom = myNodal % interp % qInterp % EvaluateLagrangePolynomial( -ONE )

      myNodal % lagTop = myNodal % interp % qInterp % EvaluateLagrangePolynomial( ONE )

      ! Build the derivative matrices
      CALL myNodal % interp % CalculateDerivativeMatrix( tempDmatS, tempDmatP, tempDmatQ )

      IF( thisForm == CG )then ! Continuous Galerkin, store the derivative matrices as is

         myNodal % dMatS = tempDmatS
         myNodal % dMatP = tempDmatP
         myNodal % dMatQ = tempDmatQ

      ELSEIF( thisForm == DG )then
      
         ! For Discontinuous Galerkin, the matrix is transposed and multiplied by a ratio of quadrature
         ! weights.
        DO jP = 0, thisN ! loop over the matrix rows
            DO iP = 0, thisN ! loop over the matrix columns

               myNodal % dMatS(iP,jP) = -tempDmatS(jP,iP)*&
                                         myNodal % qWeightS(jP)/&
                                         myNodal % qWeightS(iP)

            ENDDO
         ENDDO
         
         DO jP = 0, thisM ! loop over the matrix rows
            DO iP = 0, thisM ! loop over the matrix columns

               myNodal % dMatP(iP,jP) = -tempDmatP(jP,iP)*&
                                         myNodal % qWeightP(jP)/&
                                         myNodal % qWeightP(iP)

            ENDDO
         ENDDO
         
         DO jP = 0, thisP ! loop over the matrix rows
            DO iP = 0, thisP ! loop over the matrix columns

               myNodal % dMatQ(iP,jP) = -tempDmatQ(jP,iP)*&
                                         myNodal % qWeightQ(jP)/&
                                         myNodal % qWeightQ(iP)

            ENDDO
         ENDDO
        
      ELSE

         PRINT*,'Module NodalStorage_3D_Class.f90 : S/R BuildNodalStorage_3D : Invalid SEM form. Stopping'
         STOP

      ENDIF

      DEALLOCATE( tempS, tempWs, tempDmatS )
      DEALLOCATE( tempP, tempWp, tempDmatP )
      DEALLOCATE( tempQ, tempWq, tempDmatQ )

 END SUBROUTINE BuildNodalStorage_3D
!
!
!
 SUBROUTINE TrashNodalStorage_3D( myNodal)
 ! S/R TrashNodalStorage_3D
 ! 
 !  
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS(NodalStorage_3D), INTENT(inout) :: myNodal

   CALL myNodal % interp % TRASH( )

   DEALLOCATE( myNodal % lagWest, myNodal % lagEast )
   DEALLOCATE( myNodal % lagSouth, myNodal % lagNorth )
   DEALLOCATE( myNodal % lagBottom, myNodal % lagTop )
 
   DEALLOCATE( myNodal % qWeightS, myNodal % qWeightP, myNodal % qWeightQ )

   DEALLOCATE( myNodal % dMatS, myNodal % dMatP, myNodal % dMatQ )


 END SUBROUTINE TrashNodalStorage_3D
!
!
!==================================================================================================!
!--------------------------------------- Accessors ------------------------------------------------!
!==================================================================================================!
!
!
 SUBROUTINE SetNumberOfNodes_NodalStorage_3D( myNodal, nS, nP, nQ  )
 ! S/R SetNumberOfNodes
 !
 !    
 !
 ! =============================================================================================== !
 ! DECLARATIONS
    IMPLICIT NONE
    CLASS(NodalStorage_3D), INTENT(inout) :: myNodal
    INTEGER, INTENT(in)                   :: nS
    INTEGER, INTENT(in), OPTIONAL         :: nP, nQ

       myNodal % nS = nS
       IF( PRESENT(nP) )THEN
          myNodal % nP = nP
       ELSE
          myNodal % nP = nS
       ENDIF

       IF( PRESENT(nQ) )THEN
          myNodal % nQ = nQ
       ELSE
          myNodal % nQ = nS
       ENDIF
       
 END SUBROUTINE SetNumberOfNodes_NodalStorage_3D
!
!
!
 SUBROUTINE GetNumberOfNodes_NodalStorage_3D( myNodal, nS, nP, nQ  )
 ! S/R GetNumberOfNodes
 !
 !
 ! =============================================================================================== !
 ! DECLARATIONS
    IMPLICIT NONE
    CLASS(NodalStorage_3D), INTENT(in) :: myNodal
    INTEGER, INTENT(out)               :: nS, nP, nQ

      nS = myNodal % nS
      nP = myNodal % nP
      nQ = myNodal % nQ

 END SUBROUTINE GetNumberOfNodes_NodalStorage_3D
!
!
!
 SUBROUTINE SetNodes_NodalStorage_3D( myNodal, sNodes, pNodes, qNodes )
 ! S/R SetQuadratureWeights_NodalStorage_3D
 !
 !    
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS(NodalStorage_3D), INTENT(inout) :: myNodal
   REAL(prec), INTENT(in)                :: sNodes(0:myNodal % nS)
   REAL(prec), INTENT(in)                :: pNodes(0:myNodal % nP)
   REAL(prec), INTENT(in)                :: qNodes(0:myNodal % nQ)
  
       
      CALL myNodal % interp % SetNodes( sNodes, pNodes, qNodes )
       

 END SUBROUTINE SetNodes_NodalStorage_3D
!
!
!
 SUBROUTINE GetNodes_NodalStorage_3D( myNodal, sNodes, pNodes, qNodes  )
 ! S/R GetNodes
 !
 !    
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS(NodalStorage_3D), INTENT(in)  :: myNodal
   REAL(prec), INTENT(out)             :: sNodes(0:myNodal % nS)
   REAL(prec), INTENT(out)             :: pNodes(0:myNodal % nP)
   REAL(prec), INTENT(out)             :: qNodes(0:myNodal % nQ)
  
       CALL myNodal % interp % GetNodes( sNodes, pNodes, qNodes )

 END SUBROUTINE GetNodes_NodalStorage_3D
!
!
!
 SUBROUTINE SetQuadratureWeights_NodalStorage_3D( myNodal, quadWeightS, quadWeightP, quadWeightQ  )
 ! S/R SetNodes
 !
 !    
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS(NodalStorage_3D), INTENT(inout) :: myNodal
   REAL(prec), INTENT(in)                :: quadWeightS(0:myNodal % nS)
   REAL(prec), INTENT(in)                :: quadWeightP(0:myNodal % nP)
   REAL(prec), INTENT(in)                :: quadWeightQ(0:myNodal % nQ)
  
       myNodal % qWeightS = quadWeightS
       myNodal % qWeightP = quadWeightP
       myNodal % qWeightQ = quadWeightQ
       
 END SUBROUTINE SetQuadratureWeights_NodalStorage_3D
!
!
!
 SUBROUTINE GetQuadratureWeights_NodalStorage_3D( myNodal, quadWeightS, quadWeightP, quadWeightQ )
 ! S/R GetQuadratureWeights_NodalStorage_3D
 !
 !    
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS(NodalStorage_3D), INTENT(in)  :: myNodal
   REAL(prec), INTENT(out)             :: quadWeightS(0:myNodal % nS)
   REAL(prec), INTENT(out)             :: quadWeightP(0:myNodal % nP)
   REAL(prec), INTENT(out)             :: quadWeightQ(0:myNodal % nQ)
  
       quadWeightS = myNodal % qWeightS
       quadWeightP = myNodal % qWeightP
       quadWeightQ = myNodal % qWeightQ

 END SUBROUTINE GetQuadratureWeights_NodalStorage_3D
!
!
!
 SUBROUTINE SetDerivativeMatrix_NodalStorage_3D( myNodal, derMatS, derMatP, derMatQ )
 ! S/R SetDerivativeMatrix_NodalStorage_3D
 !
 !    
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS(NodalStorage_3D), INTENT(inout) :: myNodal
   REAL(prec), INTENT(in)                :: derMatS(0:myNodal % nS,0:myNodal % nS)
   REAL(prec), INTENT(in)                :: derMatP(0:myNodal % nP,0:myNodal % nP)
   REAL(prec), INTENT(in)                :: derMatQ(0:myNodal % nQ,0:myNodal % nQ)
    
      myNodal % dMatS = derMatS 
      myNodal % dMatP = derMatP
      myNodal % dMatQ = derMatQ

 END SUBROUTINE SetDerivativeMatrix_NodalStorage_3D
!
!
!
 SUBROUTINE GetDerivativeMatrix_NodalStorage_3D( myNodal, derMatS, derMatP, derMatQ )
 ! S/R GetDerivativeMatrix_NodalStorage_3D
 !
 !    
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS(NodalStorage_3D), INTENT(in) :: myNodal
   REAL(prec), INTENT(out)            :: derMatS(0:myNodal % nS,0:myNodal % nS)
   REAL(prec), INTENT(out)            :: derMatP(0:myNodal % nP,0:myNodal % nP)
   REAL(prec), INTENT(out)            :: derMatQ(0:myNodal % nQ,0:myNodal % nQ)
    
      derMatS = myNodal % dMatS
      derMatP = myNodal % dMatP
      derMatQ = myNodal % dMatQ

 END SUBROUTINE GetDerivativeMatrix_NodalStorage_3D
!
!
!
 SUBROUTINE SetEasternInterpolants_NodalStorage_3D( myNodal, lEast )
 ! S/R SetEasternInterpolants
 !
 !    
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS(NodalStorage_3D), INTENT(inout) :: myNodal
   REAL(prec), INTENT(in)                :: lEast(0:myNodal % nS)
    
      myNodal % lagEast = lEast

 END SUBROUTINE SetEasternInterpolants_NodalStorage_3D
!
!
!
 SUBROUTINE GetEasternInterpolants_NodalStorage_3D( myNodal, lEast )
 ! S/R GetEasternInterpolants
 !
 !    
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS(NodalStorage_3D), INTENT(in) :: myNodal
   REAL(prec), INTENT(out)            :: lEast(0:myNodal % nS)
    
      lEast = myNodal % lagEast

 END SUBROUTINE GetEasternInterpolants_NodalStorage_3D
!
!
!
 SUBROUTINE SetWesternInterpolants_NodalStorage_3D( myNodal, lWest )
 ! S/R SetWesternInterpolants
 !
 !    
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS(NodalStorage_3D), INTENT(inout) :: myNodal
   REAL(prec), INTENT(in)                :: lWest(0:myNodal % nS)
    
      myNodal % lagWest = lWest

 END SUBROUTINE SetWesternInterpolants_NodalStorage_3D
!
!
!
 SUBROUTINE GetWesternInterpolants_NodalStorage_3D( myNodal, lWest )
 ! S/R GetWesternInterpolants
 !
 !    
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS(NodalStorage_3D), INTENT(in) :: myNodal
   REAL(prec), INTENT(out)            :: lWest(0:myNodal % nS)
    
      lWest = myNodal % lagWest

 END SUBROUTINE GetWesternInterpolants_NodalStorage_3D
!
!
!
 SUBROUTINE SetNorthernInterpolants_NodalStorage_3D( myNodal, lNorth )
 ! S/R SetNorthernInterpolants
 !
 !    
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS(NodalStorage_3D), INTENT(inout) :: myNodal
   REAL(prec), INTENT(in)                :: lNorth(0:myNodal % nP)
    
      myNodal % lagNorth = lNorth

 END SUBROUTINE SetNorthernInterpolants_NodalStorage_3D
!
!
!
 SUBROUTINE GetNorthernInterpolants_NodalStorage_3D( myNodal, lNorth )
 ! S/R GetNorthernInterpolants
 !
 !    
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS(NodalStorage_3D), INTENT(in) :: myNodal
   REAL(prec), INTENT(out)            :: lNorth(0:myNodal % nP)
    
      lNorth = myNodal % lagNorth

 END SUBROUTINE GetNorthernInterpolants_NodalStorage_3D
!
!
!
 SUBROUTINE SetSouthernInterpolants_NodalStorage_3D( myNodal, lSouth )
 ! S/R SetSouthernInterpolants
 !
 !    
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS(NodalStorage_3D), INTENT(inout) :: myNodal
   REAL(prec), INTENT(in)                :: lSouth(0:myNodal % nP)
    
      myNodal % lagSouth = lSouth

 END SUBROUTINE SetSouthernInterpolants_NodalStorage_3D
!
!
!
 SUBROUTINE GetSouthernInterpolants_NodalStorage_3D( myNodal, lSouth )
 ! S/R GetSouthernInterpolants
 !
 !    
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS(NodalStorage_3D), INTENT(in) :: myNodal
   REAL(prec), INTENT(out)            :: lSouth(0:myNodal % nP)
    
      lSouth = myNodal % lagSouth

 END SUBROUTINE GetSouthernInterpolants_NodalStorage_3D
!
!
!
 SUBROUTINE SetBottomInterpolants_NodalStorage_3D( myNodal, lBottom )
 ! S/R SetBottomInterpolants
 !
 !    
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS(NodalStorage_3D), INTENT(inout) :: myNodal
   REAL(prec), INTENT(in)                :: lBottom(0:myNodal % nQ)
    
      myNodal % lagBottom = lBottom

 END SUBROUTINE SetBottomInterpolants_NodalStorage_3D
!
!
!
 SUBROUTINE GetBottomInterpolants_NodalStorage_3D( myNodal, lBottom )
 ! S/R GetBottomInterpolants
 !
 !    
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS(NodalStorage_3D), INTENT(in) :: myNodal
   REAL(prec), INTENT(out)            :: lBottom(0:myNodal % nQ)
    
      lBottom = myNodal % lagBottom

 END SUBROUTINE GetBottomInterpolants_NodalStorage_3D
!
!
!
 SUBROUTINE SetTopInterpolants_NodalStorage_3D( myNodal, lTop )
 ! S/R SetTopInterpolants
 !
 !    
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS(NodalStorage_3D), INTENT(inout) :: myNodal
   REAL(prec), INTENT(in)                :: lTop(0:myNodal % nP)
    
      myNodal % lagTop = lTop

 END SUBROUTINE SetTopInterpolants_NodalStorage_3D
!
!
!
 SUBROUTINE GetTopInterpolants_NodalStorage_3D( myNodal, lTop )
 ! S/R GetTopInterpolants
 !
 !    
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS(NodalStorage_3D), INTENT(in) :: myNodal
   REAL(prec), INTENT(out)            :: lTop(0:myNodal % nP)
    
      lTop = myNodal % lagTop

 END SUBROUTINE GetTopInterpolants_NodalStorage_3D
!
!
!==================================================================================================!
!--------------------------- Type Specific Routines -----------------------------------------------!
!==================================================================================================!
!
!
 SUBROUTINE GenerateMappingMatrix_NodalStorage_3D( myNodal, nNew, xNew, plotMatS, plotMatP, plotMatQ )
 ! S/R GenerateMappingMatrix_NodalStorage_3D
 !
 !    This subroutine generates the interpolation matrix to interpolate from the quadrature nodes
 !    to a grid with M+1 points. If the grid points are not given, a grid with uniform spacing is
 !    assumed
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS(NodalStorage_3D), INTENT(in) :: myNodal
   INTEGER, INTENT(in)                :: nNew
   REAL(prec), INTENT(in), OPTIONAL   :: xNew(0:nNew)
   REAL(prec), INTENT(out)            :: plotMatS(0:nNew,0:myNodal % nS)
   REAL(prec), INTENT(out)            :: plotMatP(0:nNew,0:myNodal % nP)
   REAL(prec), INTENT(out)            :: plotMatQ(0:nNew,0:myNodal % nQ)
   ! Local
   INTEGER           :: N, M, P
   REAL(prec)        :: x(0:nNew)
   
   
      CALL myNodal % GetNumberOfNodes( N, M, P )
      
      IF( PRESENT(xNew) )THEN
         x = xNew
      ELSE
         x = UniformPoints( -ONE, ONE, nNew )
      ENDIF
      
      CALL myNodal % interp % CalculateInterpolationMatrix( nNew, nNew, nNew, &
                                                            x, x, x, &
                                                            plotMatS, plotMatP, plotMatQ )
      
   
 END SUBROUTINE GenerateMappingMatrix_NodalStorage_3D
!
!
!
! SUBROUTINE Integrate_NodalStorage_3D( myNodal, f, intF )
 ! S/R Integrate_NodalStorage_3D
 !
 !    This subroutine integrates a function by multiplying its nodal values (f) by the quadrature 
 !    weights. The result is a real value stored in intF
 !
 ! =============================================================================================== !
 ! DECLARATIONS
!   IMPLICIT NONE
!   CLASS(NodalStorage_3D), INTENT(in) :: myNodal
!   REAL(prec), INTENT(in)             :: f(0:myNodal % nS)
!   REAL(prec), INTENT(out)            :: intF
   ! Local
!   INTEGER :: N
   
!      CALL myNodal % GetNumberOfNodes( N )
      
!      intF = DOT_PRODUCT( myNodal % qWeights(0:N), f(0:N) )
   
! END SUBROUTINE Integrate_NodalStorage_3D
!
!
!
 SUBROUTINE ApplyDerivativeMatrix_NodalStorage_3D( myNodal, f, dFds, dFdp, dfdQ )
 ! S/R ApplyDerivativeMatrix_NodalStorage_3D
 !
 !    This subroutine calculates the matrix-vector product between the derivative matrix and the 
 !    array of real values specified in "f". The result is an array of real values stored in "dF"
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS(NodalStorage_3D), INTENT(in) :: myNodal
   REAL(prec), INTENT(in)             :: f(0:myNodal % nS, 0:myNodal % nP, 0:myNodal % nQ)
   REAL(prec), INTENT(out)            :: dFds(0:myNodal % nS, 0:myNodal % nP, 0:myNodal % nQ)
   REAL(prec), INTENT(out)            :: dFdp(0:myNodal % nS, 0:myNodal % nP, 0:myNodal % nQ)
   REAL(prec), INTENT(out)            :: dFdq(0:myNodal % nS, 0:myNodal % nP, 0:myNodal % nQ)
   ! Local
   INTEGER :: N, M, P, i, j
   
      CALL myNodal % GetNumberOfNodes( N, M, P )
      
      DO i = 0, M
         DO j = 0, P
            dFds(0:N,i,j) = MATMUL( myNodal % dMatS(0:N,0:N), f(0:N,i,j) )
         ENDDO
      ENDDO
   
      DO i = 0, N
         DO j = 0, P
            dFdp(i,0:M,j) = MATMUL( myNodal % dMatP(0:M,0:M), f(i,0:M,j) )
         ENDDO
      ENDDO
      
      DO i = 0, N
         DO j = 0, M
            dFdq(i,j,0:P) = MATMUL( myNodal % dMatQ(0:P,0:P), f(i,j,0:P) )
         ENDDO
      ENDDO     
      
 END SUBROUTINE ApplyDerivativeMatrix_NodalStorage_3D
!
!
!
 SUBROUTINE CalculateAtBoundaries_NodalStorage_3D( myNodal, f, fWest, fEast, fSouth, fNorth, fBottom, fTop )
 ! S/R CalculateAtBoundaries_NodalStorage_3D
 !
 !    This subroutine estimates a function at the s=-1 and s=1 whose nodal values are given in the 
 !    real array "f". fWest is the function value estimated at s=-1, and fEast is the function value
 !    estimated at s=1.
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS(NodalStorage_3D), INTENT(in) :: myNodal
   REAL(prec), INTENT(in)             :: f(0:myNodal % nS, 0:myNodal % nP, 0:myNodal % nQ)
   REAL(prec), INTENT(out)            :: fWest(0:myNodal % nP,0:myNodal % nQ)
   REAL(prec), INTENT(out)            :: fEast(0:myNodal % nP,0:myNodal % nQ)
   REAL(prec), INTENT(out)            :: fSouth(0:myNodal % nS,0:myNodal % nQ)
   REAL(prec), INTENT(out)            :: fNorth(0:myNodal % nS,0:myNodal % nQ)
   REAL(prec), INTENT(out)            :: fBottom(0:myNodal % nS,0:myNodal % nP)
   REAL(prec), INTENT(out)            :: fTop(0:myNodal % nS,0:myNodal % nP)
   ! Local
   INTEGER :: N, M, P, i, j
   
      CALL myNodal % GetNumberOfNodes( N, M, P )
      
      DO i = 0,M
         DO j = 0,P
            fWest(i,j) = DOT_PRODUCT( myNodal % lagWest(0:N), f(0:N,i,j) )
            fEast(i,j) = DOT_PRODUCT( myNodal % lagEast(0:N), f(0:N,i,j) )
         ENDDO
      ENDDO
      
      DO i = 0,N
         DO j = 0,P
            fSouth(i,j) = DOT_PRODUCT( myNodal % lagSouth(0:M), f(i,0:M,j) )
            fNorth(i,j) = DOT_PRODUCT( myNodal % lagNorth(0:M), f(i,0:M,j) )
         ENDDO
      ENDDO
      
      DO i = 0,M
         DO j = 0,N
            fBottom(i,j) = DOT_PRODUCT( myNodal % lagBottom(0:P), f(i,j,0:P) )
            fTop(i,j) = DOT_PRODUCT( myNodal % lagTop(0:P), f(i,j,0:P) )
         ENDDO
      ENDDO
      
      
 END SUBROUTINE CalculateAtBoundaries_NodalStorage_3D
!
!
!
END MODULE NodalStorage_3D_Class
