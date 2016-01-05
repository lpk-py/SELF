MODULE NodalStorage_2D_Class
! NodalStorage_2D_Class.f90
! 
! Module History 
! 
! o (v 1.0 - 9 February 2014 )
!
! o (v 2.1 - 29 November 2015)
!
! =======================================================================================
USE ConstantsDictionary
USE ModelPrecision
USE ModelFlags
USE Legendre
USE Lagrange_1D_Class
USE Lagrange_2D_Class

IMPLICIT NONE


    TYPE NodalStorage_2D
      INTEGER                          :: nS, nP
      TYPE(Lagrange_2D)                :: interp       ! Interpolant
      REAL(prec), ALLOCATABLE, PRIVATE :: qWeightS(:)  ! dimensioned (0:nS) : Quadrature weights in the "S" direction
      REAL(prec), ALLOCATABLE, PRIVATE :: qWeightP(:)  ! dimensioned (0:nP) : Quadrature weights in the "P" direction
      REAL(prec), ALLOCATABLE, PRIVATE :: dMatS(:,:)   ! dimensioned (0:nS, 0:nS) : Contains the s-derivative matrix
      REAL(prec), ALLOCATABLE, PRIVATE :: dMatP(:,:)   ! dimensioned (0:nP, 0:nP) : Contains the p-derivative matrix
      REAL(prec), ALLOCATABLE, PRIVATE :: lagWest(:)   ! s-lagrange polynomials at x=-1, dimensioned (0:nS)
      REAL(prec), ALLOCATABLE, PRIVATE :: lagEast(:)   ! s-lagrange polynomials at x=1, dimensioned (0:nS)
      REAL(prec), ALLOCATABLE, PRIVATE :: lagSouth(:)  ! p-lagrange polynomials at x=-1, dimensioned (0:nP)
      REAL(prec), ALLOCATABLE, PRIVATE :: lagNorth(:)  ! p-lagrange polynomials at x=1, dimensioned (0:nP)


      CONTAINS

      ! Manual Constructors/Destructors
      PROCEDURE :: Build => BuildNodalStorage_2D
      PROCEDURE :: Trash => TrashNodalStorage_2D

      ! Accessors
      PROCEDURE :: SetNumberOfNodes => SetNumberOfNodes_NodalStorage_2D
      PROCEDURE :: GetNumberOfNodes => GetNumberOfNodes_NodalStorage_2D
      PROCEDURE :: SetNodes => SetNodes_NodalStorage_2D
      PROCEDURE :: GetNodes => GetNodes_NodalStorage_2D
      PROCEDURE :: SetQuadratureWeights => SetQuadratureWeights_NodalStorage_2D
      PROCEDURE :: GetQuadratureWeights => GetQuadratureWeights_NodalStorage_2D
      PROCEDURE :: GetDerivativeMatrix => GetDerivativeMatrix_NodalStorage_2D
      PROCEDURE :: SetEasternInterpolants => SetEasternInterpolants_NodalStorage_2D
      PROCEDURE :: GetEasternInterpolants => GetEasternInterpolants_NodalStorage_2D
      PROCEDURE :: SetWesternInterpolants => SetWesternInterpolants_NodalStorage_2D
      PROCEDURE :: GetWesternInterpolants => GetWesternInterpolants_NodalStorage_2D
      PROCEDURE :: SetSouthernInterpolants => SetSouthernInterpolants_NodalStorage_2D
      PROCEDURE :: GetSouthernInterpolants => GetSouthernInterpolants_NodalStorage_2D
      PROCEDURE :: SetNorthernInterpolants => SetNorthernInterpolants_NodalStorage_2D
      PROCEDURE :: GetNorthernInterpolants => GetNorthernInterpolants_NodalStorage_2D

      ! Type-Specific
      PROCEDURE :: GenerateMappingMatrix => GenerateMappingMatrix_NodalStorage_2D
      PROCEDURE :: ApplyDerivativeMatrix => ApplyDerivativeMatrix_NodalStorage_2D
      PROCEDURE :: CalculateAtBoundaries => CalculateAtBoundaries_NodalStorage_2D
      !PROCEDURE :: Integrate => Integrate_NodalStorage_2D
      
      
    END TYPE NodalStorage_2D


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
 SUBROUTINE BuildNodalStorage_2D( myNodal, N, M, quadrature, approxForm  )
 ! S/R BuildNodalStorage_2D
 !  
 !   This constructor allocates space for the NodalStorage_2D class and fills in the quadrature and
 !   derivative matrix. The type of quadrature can be set to GAUSS or GAUSS_LOBATTO; this determines
 !   the nodes that are generated and stored. The approximation form can be CG for Continuous 
 !   Galerkin or DG  for Discontinuous Galerkin. The number of nodes (minus 1) is N. These parameters
 !   are all optional. Therefore, to call this function with the optional parameters, use the 
 !   syntax :
 !
 !   CALL myNodal % Build( N          = numberOfPoints, &
 !                         M          = numberOfPoints, &
 !                         quadrature = yourQuadrature, &
 !                         approxForm = yourForm          )
 !
 !   If none of the optional arguments are given, the Build defaults to a linear interpolant through
 !   s=-1 and s=1 and the Continuous Galerkin approximation is used giving the standard derivative
 !   matrix.
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS(NodalStorage_2D), INTENT(out) :: myNodal
   INTEGER, INTENT(in), optional       :: N, M, quadrature, approxForm
   !LOCAL
   INTEGER                 :: iP, jP, thisQuadrature, thisN, thisM, thisForm
   REAL(prec), ALLOCATABLE :: tempS(:), tempWs(:), tempDmatS(:,:)
   REAL(prec), ALLOCATABLE :: tempP(:), tempWp(:), tempDmatP(:,:)

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

      ALLOCATE( myNodal % lagWest(0:thisN), myNodal % lagEast(0:thisN) )
      ALLOCATE( myNodal % lagSouth(0:thisM), myNodal % lagNorth(0:thisM) )

      ALLOCATE( myNodal % dMatS(0:thisN,0:thisN) )
      ALLOCATE( myNodal % dMatP(0:thisM,0:thisM) )

      ALLOCATE( myNodal % qWeightS(0:thisN) )
      ALLOCATE( myNodal % qWeightP(0:thisM) )

      ! Generate the quadrature
      CALL GenerateLegendreQuadrature( thisN, tempS, tempWs, thisQuadrature )
      CALL GenerateLegendreQuadrature( thisM, tempP, tempWp, thisQuadrature )
      
      ! Store the quadrature weights
      CALL myNodal % SetQuadratureWeights( tempWs, tempWp )

      ! Build and store the interpolant
      CALL myNodal % interp % Build( thisN, thisM, tempS, tempP )
   
      ! Calculate and store the interpolants evaluated at the endpoints
      myNodal % lagWest = myNodal % interp % sInterp % EvaluateLagrangePolynomial( -ONE )

      myNodal % lagEast = myNodal % interp % sInterp % EvaluateLagrangePolynomial( ONE )
      
      myNodal % lagSouth = myNodal % interp % pInterp % EvaluateLagrangePolynomial( -ONE )

      myNodal % lagNorth = myNodal % interp % pInterp % EvaluateLagrangePolynomial( ONE )

      ! Build the derivative matrices
      CALL myNodal % interp % CalculateDerivativeMatrix( tempDmatS, tempDmatP )

      IF( thisForm == CG )then ! Continuous Galerkin, store the derivative matrices as is

         myNodal % dMatS = tempDmatS
         myNodal % dMatP = tempDmatP

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
        
      ELSE

         PRINT*,'Module NodalStorage_2D_Class.f90 : S/R BuildNodalStorage_2D : Invalid SEM form. Stopping'
         STOP

      ENDIF

      DEALLOCATE( tempS, tempWs, tempDmatS )
      DEALLOCATE( tempP, tempWp, tempDmatP )

 END SUBROUTINE BuildNodalStorage_2D
!
!
!
 SUBROUTINE TrashNodalStorage_2D( myNodal)
 ! S/R TrashNodalStorage_2D
 ! 
 !  
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS(NodalStorage_2D), INTENT(inout) :: myNodal

   CALL myNodal % interp % TRASH( )

   DEALLOCATE( myNodal % lagWest, myNodal % lagEast )
   DEALLOCATE( myNodal % lagSouth, myNodal % lagNorth )
 
   DEALLOCATE( myNodal % qWeightS, myNodal % qWeightP )

   DEALLOCATE( myNodal % dMatS, myNodal % dMatP )


 END SUBROUTINE TrashNodalStorage_2D
!
!
!==================================================================================================!
!--------------------------------------- Accessors ------------------------------------------------!
!==================================================================================================!
!
!
 SUBROUTINE SetNumberOfNodes_NodalStorage_2D( myNodal, nS, nP  )
 ! S/R SetNumberOfNodes
 !
 !    
 !
 ! =============================================================================================== !
 ! DECLARATIONS
    IMPLICIT NONE
    CLASS(NodalStorage_2D), INTENT(inout) :: myNodal
    INTEGER, INTENT(in)                   :: nS
    INTEGER, INTENT(in), OPTIONAL         :: nP

       myNodal % nS = nS
       IF( PRESENT(nP) )THEN
          myNodal % nP = nP
       ELSE
          myNodal % nP = nS
       ENDIF

 END SUBROUTINE SetNumberOfNodes_NodalStorage_2D
!
!
!
 SUBROUTINE GetNumberOfNodes_NodalStorage_2D( myNodal, nS, nP  )
 ! S/R GetNumberOfNodes
 !
 !
 ! =============================================================================================== !
 ! DECLARATIONS
    IMPLICIT NONE
    CLASS(NodalStorage_2D), INTENT(in) :: myNodal
    INTEGER, INTENT(out)               :: nS, nP

      nS = myNodal % nS
      nP = myNodal % nP

 END SUBROUTINE GetNumberOfNodes_NodalStorage_2D
!
!
!
 SUBROUTINE SetNodes_NodalStorage_2D( myNodal, sNodes, pNodes )
 ! S/R SetQuadratureWeights_NodalStorage_2D
 !
 !    
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS(NodalStorage_2D), INTENT(inout) :: myNodal
   REAL(prec), INTENT(in)                :: sNodes(0:myNodal % nS)
   REAL(prec), INTENT(in), OPTIONAL      :: pNodes(0:myNodal % nP)
  
       IF( PRESENT(pNodes) )THEN
          CALL myNodal % interp % SetNodes( sNodes, pNodes )
       ELSE
          CALL myNodal % interp % SetNodes( sNodes, sNodes )
       ENDIF

 END SUBROUTINE SetNodes_NodalStorage_2D
!
!
!
 SUBROUTINE GetNodes_NodalStorage_2D( myNodal, sNodes, pNodes  )
 ! S/R GetNodes
 !
 !    
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS(NodalStorage_2D), INTENT(in)  :: myNodal
   REAL(prec), INTENT(out)             :: sNodes(0:myNodal % nS)
   REAL(prec), INTENT(out)             :: pNodes(0:myNodal % nP)
  
       CALL myNodal % interp % GetNodes( sNodes, pNodes )

 END SUBROUTINE GetNodes_NodalStorage_2D
!
!
!
 SUBROUTINE SetQuadratureWeights_NodalStorage_2D( myNodal, quadWeightS, quadWeightP  )
 ! S/R SetNodes
 !
 !    
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS(NodalStorage_2D), INTENT(inout) :: myNodal
   REAL(prec), INTENT(in)                :: quadWeightS(0:myNodal % nS)
   REAL(prec), INTENT(in), OPTIONAL      :: quadWeightP(0:myNodal % nP)
  
       myNodal % qWeightS = quadWeightS
       IF( PRESENT(quadWeightP) )THEN
          myNodal % qWeightP = quadWeightP
       ELSE
          myNodal % qWeightP = quadWeightS
       ENDIF

 END SUBROUTINE SetQuadratureWeights_NodalStorage_2D
!
!
!
 SUBROUTINE GetQuadratureWeights_NodalStorage_2D( myNodal, quadWeightS, quadWeightP  )
 ! S/R GetQuadratureWeights_NodalStorage_2D
 !
 !    
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS(NodalStorage_2D), INTENT(in)  :: myNodal
   REAL(prec), INTENT(out)             :: quadWeightS(0:myNodal % nS)
   REAL(prec), INTENT(out)             :: quadWeightP(0:myNodal % nP)
  
       quadWeightS = myNodal % qWeightS
       quadWeightP = myNodal % qWeightP

 END SUBROUTINE GetQuadratureWeights_NodalStorage_2D
!
!
!
 SUBROUTINE SetDerivativeMatrix_NodalStorage_2D( myNodal, derMatS, derMatP )
 ! S/R SetDerivativeMatrix_NodalStorage_2D
 !
 !    
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS(NodalStorage_2D), INTENT(inout) :: myNodal
   REAL(prec), INTENT(in)                :: derMatS(0:myNodal % nS,0:myNodal % nS)
   REAL(prec), INTENT(in)                :: derMatP(0:myNodal % nP,0:myNodal % nP)
    
      myNodal % dMatS = derMatS 
      myNodal % dMatP = derMatP

 END SUBROUTINE SetDerivativeMatrix_NodalStorage_2D
!
!
!
 SUBROUTINE GetDerivativeMatrix_NodalStorage_2D( myNodal, derMatS, derMatP )
 ! S/R GetDerivativeMatrix_NodalStorage_2D
 !
 !    
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS(NodalStorage_2D), INTENT(in) :: myNodal
   REAL(prec), INTENT(out)            :: derMatS(0:myNodal % nS,0:myNodal % nS)
   REAL(prec), INTENT(out)            :: derMatP(0:myNodal % nP,0:myNodal % nP)
    
      derMatS = myNodal % dMatS
      derMatP = myNodal % dMatP

 END SUBROUTINE GetDerivativeMatrix_NodalStorage_2D
!
!
!
 SUBROUTINE SetEasternInterpolants_NodalStorage_2D( myNodal, lEast )
 ! S/R SetEasternInterpolants
 !
 !    
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS(NodalStorage_2D), INTENT(inout) :: myNodal
   REAL(prec), INTENT(in)                :: lEast(0:myNodal % nS)
    
      myNodal % lagEast = lEast

 END SUBROUTINE SetEasternInterpolants_NodalStorage_2D
!
!
!
 SUBROUTINE GetEasternInterpolants_NodalStorage_2D( myNodal, lEast )
 ! S/R GetEasternInterpolants
 !
 !    
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS(NodalStorage_2D), INTENT(in) :: myNodal
   REAL(prec), INTENT(out)            :: lEast(0:myNodal % nS)
    
      lEast = myNodal % lagEast

 END SUBROUTINE GetEasternInterpolants_NodalStorage_2D
!
!
!
 SUBROUTINE SetWesternInterpolants_NodalStorage_2D( myNodal, lWest )
 ! S/R SetWesternInterpolants
 !
 !    
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS(NodalStorage_2D), INTENT(inout) :: myNodal
   REAL(prec), INTENT(in)                :: lWest(0:myNodal % nS)
    
      myNodal % lagWest = lWest

 END SUBROUTINE SetWesternInterpolants_NodalStorage_2D
!
!
!
 SUBROUTINE GetWesternInterpolants_NodalStorage_2D( myNodal, lWest )
 ! S/R GetWesternInterpolants
 !
 !    
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS(NodalStorage_2D), INTENT(in) :: myNodal
   REAL(prec), INTENT(out)            :: lWest(0:myNodal % nS)
    
      lWest = myNodal % lagWest

 END SUBROUTINE GetWesternInterpolants_NodalStorage_2D
!
!
!
 SUBROUTINE SetNorthernInterpolants_NodalStorage_2D( myNodal, lNorth )
 ! S/R SetNorthernInterpolants
 !
 !    
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS(NodalStorage_2D), INTENT(inout) :: myNodal
   REAL(prec), INTENT(in)                :: lNorth(0:myNodal % nP)
    
      myNodal % lagNorth = lNorth

 END SUBROUTINE SetNorthernInterpolants_NodalStorage_2D
!
!
!
 SUBROUTINE GetNorthernInterpolants_NodalStorage_2D( myNodal, lNorth )
 ! S/R GetNorthernInterpolants
 !
 !    
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS(NodalStorage_2D), INTENT(in) :: myNodal
   REAL(prec), INTENT(out)            :: lNorth(0:myNodal % nP)
    
      lNorth = myNodal % lagNorth

 END SUBROUTINE GetNorthernInterpolants_NodalStorage_2D
!
!
!
 SUBROUTINE SetSouthernInterpolants_NodalStorage_2D( myNodal, lSouth )
 ! S/R SetSouthernInterpolants
 !
 !    
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS(NodalStorage_2D), INTENT(inout) :: myNodal
   REAL(prec), INTENT(in)                :: lSouth(0:myNodal % nP)
    
      myNodal % lagSouth = lSouth

 END SUBROUTINE SetSouthernInterpolants_NodalStorage_2D
!
!
!
 SUBROUTINE GetSouthernInterpolants_NodalStorage_2D( myNodal, lSouth )
 ! S/R GetSouthernInterpolants
 !
 !    
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS(NodalStorage_2D), INTENT(in) :: myNodal
   REAL(prec), INTENT(out)            :: lSouth(0:myNodal % nP)
    
      lSouth = myNodal % lagSouth

 END SUBROUTINE GetSouthernInterpolants_NodalStorage_2D
!
!
!==================================================================================================!
!--------------------------- Type Specific Routines -----------------------------------------------!
!==================================================================================================!
!
!
 SUBROUTINE GenerateMappingMatrix_NodalStorage_2D( myNodal, nNew, xNew, plotMatS, plotMatP )
 ! S/R GenerateMappingMatrix_NodalStorage_2D
 !
 !    This subroutine generates the interpolation matrix to interpolate from the quadrature nodes
 !    to a grid with M+1 points. If the grid points are not given, a grid with uniform spacing is
 !    assumed
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS(NodalStorage_2D), INTENT(in) :: myNodal
   INTEGER, INTENT(in)                :: nNew
   REAL(prec), INTENT(in), OPTIONAL   :: xNew(0:nNew)
   REAL(prec), INTENT(out)            :: plotMatS(0:nNew,0:myNodal % nS)
   REAL(prec), INTENT(out)            :: plotMatP(0:nNew,0:myNodal % nP)
   ! Local
   INTEGER           :: N, M
   REAL(prec)        :: x(0:nNew)
   
   
      CALL myNodal % GetNumberOfNodes( N, M )
      
      IF( PRESENT(xNew) )THEN
         x = xNew
      ELSE
         x = UniformPoints( -ONE, ONE, nNew )
      ENDIF
      
      CALL myNodal % interp % CalculateInterpolationMatrix( nNew, nNew, x, x, plotMatS, plotMatP )
      
   
 END SUBROUTINE GenerateMappingMatrix_NodalStorage_2D
!
!
!
! SUBROUTINE Integrate_NodalStorage_2D( myNodal, f, intF )
 ! S/R Integrate_NodalStorage_2D
 !
 !    This subroutine integrates a function by multiplying its nodal values (f) by the quadrature 
 !    weights. The result is a real value stored in intF
 !
 ! =============================================================================================== !
 ! DECLARATIONS
!   IMPLICIT NONE
!   CLASS(NodalStorage_2D), INTENT(in) :: myNodal
!   REAL(prec), INTENT(in)             :: f(0:myNodal % nS)
!   REAL(prec), INTENT(out)            :: intF
   ! Local
!   INTEGER :: N
   
!      CALL myNodal % GetNumberOfNodes( N )
      
!      intF = DOT_PRODUCT( myNodal % qWeights(0:N), f(0:N) )
   
! END SUBROUTINE Integrate_NodalStorage_2D
!
!
!
 SUBROUTINE ApplyDerivativeMatrix_NodalStorage_2D( myNodal, f, dFds, dFdp )
 ! S/R ApplyDerivativeMatrix_NodalStorage_2D
 !
 !    This subroutine calculates the matrix-vector product between the derivative matrix and the 
 !    array of real values specified in "f". The result is an array of real values stored in "dF"
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS(NodalStorage_2D), INTENT(in) :: myNodal
   REAL(prec), INTENT(in)             :: f(0:myNodal % nS, 0:myNodal % nP)
   REAL(prec), INTENT(out)            :: dFds(0:myNodal % nS, 0:myNodal % nP)
   REAL(prec), INTENT(out)            :: dFdp(0:myNodal % nS, 0:myNodal % nP)
   ! Local
   INTEGER :: N, M, i
   
      CALL myNodal % GetNumberOfNodes( N, M )
      
      DO i = 0, M
         dFds(0:N,i) = MATMUL( myNodal % dMatS(0:N,0:N), f(0:N,i) )
      ENDDO
   
      DO i = 0, N
         dFdp(i,0:M) = MATMUL( myNodal % dMatP(0:M,0:M), f(i,0:M) )
      ENDDO
      
 END SUBROUTINE ApplyDerivativeMatrix_NodalStorage_2D
!
!
!
 SUBROUTINE CalculateAtBoundaries_NodalStorage_2D( myNodal, f, fWest, fEast, fSouth, fNorth )
 ! S/R CalculateAtBoundaries_NodalStorage_2D
 !
 !    This subroutine estimates a function at the s=-1 and s=1 whose nodal values are given in the 
 !    real array "f". fWest is the function value estimated at s=-1, and fEast is the function value
 !    estimated at s=1.
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS(NodalStorage_2D), INTENT(in) :: myNodal
   REAL(prec), INTENT(in)             :: f(0:myNodal % nS, 0:myNodal % nP)
   REAL(prec), INTENT(out)            :: fWest(0:myNodal % nP), fEast(0:myNodal % nP)
   REAL(prec), INTENT(out)            :: fSouth(0:myNodal % nS), fNorth(0:myNodal % nS)
   ! Local
   INTEGER :: N, M, i
   
      CALL myNodal % GetNumberOfNodes( N, M )
      
      DO i = 0,M
         fWest(i) = DOT_PRODUCT( myNodal % lagWest(0:N), f(0:N,i) )
         fEast(i) = DOT_PRODUCT( myNodal % lagEast(0:N), f(0:N,i) )
      ENDDO
      
      DO i = 0,N
         fSouth(i) = DOT_PRODUCT( myNodal % lagSouth(0:M), f(i,0:M) )
         fNorth(i) = DOT_PRODUCT( myNodal % lagNorth(0:M), f(i,0:M) )
      ENDDO
      
      
 END SUBROUTINE CalculateAtBoundaries_NodalStorage_2D
!
!
!
END MODULE NodalStorage_2D_Class
