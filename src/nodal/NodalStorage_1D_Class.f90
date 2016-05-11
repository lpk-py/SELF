! NodalStorage_1D_Class.f90
! 
! Copyright 2015 Joseph Schoonover <schoonover.numerics@gmail.com>
! 
! NodalStorage_1D_Class.f90 is part of the Spectral Element Libraries in Fortran (SELF).
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
 
 
MODULE NodalStorage_1D_Class
! ========================================= Logs ================================================= !
!2016-05-11  Joseph Schoonover  schoonover.numerics@gmail.com 
!
! //////////////////////////////////////////////////////////////////////////////////////////////// ! 
USE ModelPrecision
USE ConstantsDictionary
USE ModelFlags
USE Legendre
USE Lagrange_1D_Class

IMPLICIT NONE


    TYPE NodalStorage_1D
      INTEGER                          :: nS
      TYPE(Lagrange_1D)                :: interp      ! Interpolant
      REAL(prec), ALLOCATABLE, PRIVATE :: qWeight(:)  ! dimensioned (0:nS) : Quadrature weights in the "S" direction
      REAL(prec), ALLOCATABLE, PRIVATE :: dMat(:,:)   ! dimensioned (0:nS, 0:nS) : Contains the s-derivative matrix
      REAL(prec), ALLOCATABLE, PRIVATE :: lagWest(:)  ! x-lagrange polynomials at x=-1, dimensioned (0:M)
      REAL(prec), ALLOCATABLE, PRIVATE :: lagEast(:)  ! x-lagrange polynomials at x=1, dimensioned (0:M)

      CONTAINS

      ! Manual Constructors/Destructors
      PROCEDURE :: Build => BuildNodalStorage_1D
      PROCEDURE :: Trash => TrashNodalStorage_1D

      ! Accessors
      PROCEDURE :: SetNumberOfNodes => SetNumberOfNodes_NodalStorage_1D
      PROCEDURE :: GetNumberOfNodes => GetNumberOfNodes_NodalStorage_1D
      PROCEDURE :: SetNodes => SetNodes_NodalStorage_1D
      PROCEDURE :: GetNodes => GetNodes_NodalStorage_1D
      PROCEDURE :: SetQuadratureWeights => SetQuadratureWeights_NodalStorage_1D
      PROCEDURE :: GetQuadratureWeights => GetQuadratureWeights_NodalStorage_1D
      PROCEDURE :: GetDerivativeMatrix => GetDerivativeMatrix_NodalStorage_1D
      PROCEDURE :: SetEasternInterpolants => SetEasternInterpolants_NodalStorage_1D
      PROCEDURE :: GetEasternInterpolants => GetEasternInterpolants_NodalStorage_1D
      PROCEDURE :: SetWesternInterpolants => SetWesternInterpolants_NodalStorage_1D
      PROCEDURE :: GetWesternInterpolants => GetWesternInterpolants_NodalStorage_1D

      ! Type-Specific
      PROCEDURE :: GenerateMappingMatrix => GenerateMappingMatrix_NodalStorage_1D
      PROCEDURE :: ApplyDerivativeMatrix => ApplyDerivativeMatrix_NodalStorage_1D
      PROCEDURE :: CalculateAtBoundaries => CalculateAtBoundaries_NodalStorage_1D
      PROCEDURE :: Integrate => Integrate_NodalStorage_1D
      
      
    END TYPE NodalStorage_1D


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
 SUBROUTINE BuildNodalStorage_1D( myNodal, N, quadrature, approxForm  )
 ! S/R BuildNodalStorage_1D
 !  
 !   This constructor allocates space for the NodalStorage_1D class and fills in the quadrature and
 !   derivative matrix. The TYPE of quadrature can be set to GAUSS or GAUSS_LOBATTO; this determines
 !   the nodes that are generated and stored. The approximation form can be CG for Continuous 
 !   Galerkin or DG  for Discontinuous Galerkin. The number of nodes (minus 1) is N. These parameters
 !   are all optional. Therefore, to call this function with the optional parameters, use the 
 !   syntax :
 !
 !   CALL myNodal % Build( N          = numberOfPoints, &
 !                         quadrature = yourQuadrature, &
 !                         approxForm = yourForm          )
 !
 !   If none of the optional arguments are given, the Build defaults to a linear interpolant through
 !   s=-1 and s=1 and the Continuous Galerkin approximation is used giving the standard derivative
 !   matrix.
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS(NodalStorage_1D), INTENT(out) :: myNodal
   INTEGER, INTENT(in), optional       :: N, quadrature, approxForm
   !LOCAL
   INTEGER                 :: iP, jP, thisQuadrature, thisN, thisForm
   REAL(prec), ALLOCATABLE :: tempX(:), tempQ(:), tempDmat(:,:)


      ! Check and setup the optional arguments
      IF( PRESENT(N) )THEN
         thisN = N
      ELSE
         thisN = defaultN
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

      ALLOCATE( tempX(0:thisN), tempQ(0:thisN), tempDmat(0:thisN,0:thisN) )
      
      ALLOCATE( myNodal % lagWest(0:thisN), myNodal % lagEast(0:thisN) )

      ALLOCATE( myNodal % dMat(0:thisN,0:thisN) )

      ALLOCATE( myNodal % qWeight(0:thisN) )

      ! Generate the quadrature
      CALL GenerateLegendreQuadrature( thisN, tempX, tempQ, thisQuadrature )
   
      ! Store the quadrature weights
      CALL myNodal % SetQuadratureWeights( tempQ )

      ! Build and store the interpolant
      CALL myNodal % interp % Build( thisN, tempX )
   
      ! Calculate and store the interpolants evaluated at the endpoints
      myNodal % lagWest = myNodal % interp % EvaluateLagrangePolynomial( -ONE )

      myNodal % lagEast = myNodal % interp % EvaluateLagrangePolynomial( ONE )

      ! Build the derivative matrices
      CALL myNodal % interp % CalculateDerivativeMatrix( tempDmat )

      IF( thisForm == CG )then ! Continuous Galerkin, store the derivative matrices as is

         myNodal % dMat = tempDmat

      ELSEIF( thisForm == DG )then
      
         ! For Discontinuous Galerkin, the matrix is transposed and multiplied by a ratio of quadrature
         ! weights.
        DO jP = 0, thisN ! loop over the matrix rows
            DO iP = 0, thisN ! loop over the matrix columns

               myNodal % dMat(iP,jP) = -tempDmat(jP,iP)*&
                                       myNodal % qWeight(jP)/&
                                       myNodal % qWeight(iP)

            ENDDO
         ENDDO
        
      ELSE

         PRINT*,'Module NodalStorage_1D_Class.f90 : S/R BuildNodalStorage_1D : Invalid SEM form. Stopping'
         STOP

      ENDIF 

      DEALLOCATE( tempX, tempQ, tempDmat )

 END SUBROUTINE BuildNodalStorage_1D
!
!
!
 SUBROUTINE TrashNodalStorage_1D( myNodal)
 ! S/R TrashNodalStorage_1D
 ! 
 !  
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS(NodalStorage_1D), INTENT(inout) :: myNodal

   CALL myNodal % interp % TRASH( )

   DEALLOCATE( myNodal % lagWest, myNodal % lagEast )
 
   DEALLOCATE( myNodal % qWeight )

   DEALLOCATE( myNodal % dMat )


 END SUBROUTINE TrashNodalStorage_1D
!
!
!==================================================================================================!
!--------------------------------------- Accessors ------------------------------------------------!
!==================================================================================================!
!
!
 SUBROUTINE SetNumberOfNodes_NodalStorage_1D( myNodal, nS  )
 ! S/R SetNumberOfNodes
 !
 !    
 !
 ! =============================================================================================== !
 ! DECLARATIONS
    IMPLICIT NONE
    CLASS(NodalStorage_1D), INTENT(inout) :: myNodal
    INTEGER, INTENT(in)                   :: nS

       myNodal % nS = nS

 END SUBROUTINE SetNumberOfNodes_NodalStorage_1D
!
!
!
 SUBROUTINE GetNumberOfNodes_NodalStorage_1D( myNodal, nS  )
 ! S/R GetNumberOfNodes
 !
 !
 ! =============================================================================================== !
 ! DECLARATIONS
    IMPLICIT NONE
    CLASS(NodalStorage_1D), INTENT(in) :: myNodal
    INTEGER, INTENT(out)               :: nS

      nS = myNodal % nS

 END SUBROUTINE GetNumberOfNodes_NodalStorage_1D
!
!
!
 SUBROUTINE SetNodes_NodalStorage_1D( myNodal, nodes  )
 ! S/R SetQuadratureWeights_NodalStorage_1D
 !
 !    
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS(NodalStorage_1D), INTENT(inout) :: myNodal
   REAL(prec), INTENT(in)                :: nodes(0:myNodal % nS)
  
       CALL myNodal % interp % SetNodes( nodes )

 END SUBROUTINE SetNodes_NodalStorage_1D
!
!
!
 SUBROUTINE GetNodes_NodalStorage_1D( myNodal, nodes  )
 ! S/R GetNodes
 !
 !    
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS(NodalStorage_1D), INTENT(in)  :: myNodal
   REAL(prec), INTENT(out)             :: nodes(0:myNodal % nS)
  
       CALL myNodal % interp % GetNodes( nodes )

 END SUBROUTINE GetNodes_NodalStorage_1D
!
!
!
 SUBROUTINE SetQuadratureWeights_NodalStorage_1D( myNodal, quadWeight  )
 ! S/R SetNodes
 !
 !    
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS(NodalStorage_1D), INTENT(inout) :: myNodal
   REAL(prec), INTENT(in)                :: quadWeight(0:myNodal % nS)
  
       myNodal % qWeight = quadWeight

 END SUBROUTINE SetQuadratureWeights_NodalStorage_1D
!
!
!
 SUBROUTINE GetQuadratureWeights_NodalStorage_1D( myNodal, quadWeight  )
 ! S/R GetQuadratureWeights_NodalStorage_1D
 !
 !    
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS(NodalStorage_1D), INTENT(in)  :: myNodal
   REAL(prec), INTENT(out)             :: quadWeight(0:myNodal % nS)
  
       quadWeight = myNodal % qWeight

 END SUBROUTINE GetQuadratureWeights_NodalStorage_1D
!
!
!
 SUBROUTINE SetDerivativeMatrix_NodalStorage_1D( myNodal, derMat )
 ! S/R SetDerivativeMatrix_NodalStorage_1D
 !
 !    
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS(NodalStorage_1D), INTENT(inout) :: myNodal
   REAL(prec), INTENT(in)                :: derMat(0:myNodal % nS,0:myNodal % nS)
    
      myNodal % dMat = derMat

 END SUBROUTINE SetDerivativeMatrix_NodalStorage_1D
!
!
!
 SUBROUTINE GetDerivativeMatrix_NodalStorage_1D( myNodal, derMat )
 ! S/R GetDerivativeMatrix_NodalStorage_1D
 !
 !    
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS(NodalStorage_1D), INTENT(in) :: myNodal
   REAL(prec), INTENT(out)            :: derMat(0:myNodal % nS,0:myNodal % nS)
    
      derMat = myNodal % dMat

 END SUBROUTINE GetDerivativeMatrix_NodalStorage_1D
!
!
!
 SUBROUTINE SetEasternInterpolants_NodalStorage_1D( myNodal, lEast )
 ! S/R SetEasternInterpolants
 !
 !    
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS(NodalStorage_1D), INTENT(inout) :: myNodal
   REAL(prec), INTENT(in)                :: lEast(0:myNodal % nS)
    
      myNodal % lagEast = lEast

 END SUBROUTINE SetEasternInterpolants_NodalStorage_1D
!
!
!
 SUBROUTINE GetEasternInterpolants_NodalStorage_1D( myNodal, lEast )
 ! S/R GetEasternInterpolants
 !
 !    
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS(NodalStorage_1D), INTENT(in) :: myNodal
   REAL(prec), INTENT(out)            :: lEast(0:myNodal % nS)
    
      lEast = myNodal % lagEast

 END SUBROUTINE GetEasternInterpolants_NodalStorage_1D
!
!
!
 SUBROUTINE SetWesternInterpolants_NodalStorage_1D( myNodal, lWest )
 ! S/R SetWesternInterpolants
 !
 !    
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS(NodalStorage_1D), INTENT(inout) :: myNodal
   REAL(prec), INTENT(in)                :: lWest(0:myNodal % nS)
    
      myNodal % lagWest = lWest

 END SUBROUTINE SetWesternInterpolants_NodalStorage_1D
!
!
!
 SUBROUTINE GetWesternInterpolants_NodalStorage_1D( myNodal, lWest )
 ! S/R GetWesternInterpolants
 !
 !    
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS(NodalStorage_1D), INTENT(in) :: myNodal
   REAL(prec), INTENT(out)            :: lWest(0:myNodal % nS)
    
      lWest = myNodal % lagWest

 END SUBROUTINE GetWesternInterpolants_NodalStorage_1D
!
!
!==================================================================================================!
!--------------------------- Type Specific Routines -----------------------------------------------!
!==================================================================================================!
!
!
 SUBROUTINE GenerateMappingMatrix_NodalStorage_1D( myNodal, M, xNew, plotMat )
 ! S/R GenerateMappingMatrix_NodalStorage_1D
 !
 !    This subroutine generates the interpolation matrix to interpolate from the quadrature nodes
 !    to a grid with M+1 points. If the grid points are not given, a grid with uniform spacing is
 !    assumed
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS(NodalStorage_1D), INTENT(in) :: myNodal
   INTEGER, INTENT(in)                :: M
   REAL(prec), INTENT(in), OPTIONAL   :: xNew(0:M)
   REAL(prec), INTENT(out)            :: plotMat(0:M,0:myNodal % nS)
   ! Local
   INTEGER           :: N
   REAL(prec)        :: x(0:M)
   
   
      CALL myNodal % GetNumberOfNodes( N )
      
      IF( PRESENT(xNew) )THEN
         x = xNew
      ELSE
         x = UniformPoints( -ONE, ONE, M )
      ENDIF
      
      CALL myNodal % interp % CalculateInterpolationMatrix( M, x, plotMat )
      
   
 END SUBROUTINE GenerateMappingMatrix_NodalStorage_1D
!
!
!
 SUBROUTINE Integrate_NodalStorage_1D( myNodal, f, intF )
 ! S/R Integrate_NodalStorage_1D
 !
 !    This subroutine integrates a function by multiplying its nodal values (f) by the quadrature 
 !    weights. The result is a real value stored in intF
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS(NodalStorage_1D), INTENT(in) :: myNodal
   REAL(prec), INTENT(in)             :: f(0:myNodal % nS)
   REAL(prec), INTENT(out)            :: intF
   ! Local
   INTEGER :: N
   
      CALL myNodal % GetNumberOfNodes( N )
      
      intF = DOT_PRODUCT( myNodal % qWeight(0:N), f(0:N) )
   
 END SUBROUTINE Integrate_NodalStorage_1D
!
!
!
 SUBROUTINE ApplyDerivativeMatrix_NodalStorage_1D( myNodal, f, dF )
 ! S/R ApplyDerivativeMatrix_NodalStorage_1D
 !
 !    This subroutine calculates the matrix-vector product between the derivative matrix and the 
 !    array of real values specified in "f". The result is an array of real values stored in "dF"
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS(NodalStorage_1D), INTENT(in) :: myNodal
   REAL(prec), INTENT(in)             :: f(0:myNodal % nS)
   REAL(prec), INTENT(out)            :: dF(0:myNodal % nS)
   ! Local
   INTEGER :: N
   
      CALL myNodal % GetNumberOfNodes( N )
      
      dF(0:N) = MATMUL( myNodal % dMat(0:N,0:N), f(0:N) )
   
 END SUBROUTINE ApplyDerivativeMatrix_NodalStorage_1D
!
!
!
 SUBROUTINE CalculateAtBoundaries_NodalStorage_1D( myNodal, f, fWest, fEast )
 ! S/R CalculateAtBoundaries_NodalStorage_1D
 !
 !    This subroutine estimates a function at the s=-1 and s=1 whose nodal values are given in the 
 !    real array "f". fWest is the function value estimated at s=-1, and fEast is the function value
 !    estimated at s=1.
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS(NodalStorage_1D), INTENT(in) :: myNodal
   REAL(prec), INTENT(in)             :: f(0:myNodal % nS)
   REAL(prec), INTENT(out)            :: fWest, fEast
   ! Local
   INTEGER :: N
   
      CALL myNodal % GetNumberOfNodes( N )
      
      fWest = DOT_PRODUCT( myNodal % lagWest(0:N), f(0:N) )
      fEast = DOT_PRODUCT( myNodal % lagEast(0:N), f(0:N) )
   
 END SUBROUTINE CalculateAtBoundaries_NodalStorage_1D
!
!
!
END MODULE NodalStorage_1D_Class
