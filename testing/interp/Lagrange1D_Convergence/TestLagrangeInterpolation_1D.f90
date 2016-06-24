PROGRAM TestLagrangeInterpolation_1D
! 
! TestLagrangeInterpolation_1D.f90
!
! This program tests the functionality of the interpolation software. A gaussian in 1-, 2-, and 3- 
! Dimensions is interpolated and written to a file. Additionally, the derivative of the interpolant
! is calculated and written to file. We use the Legendre Gauss points to perform the interpolation 
! Legendre-Gauss Quadrature to perform the estimation of the L2-Norm.
!
! The L2-Norm is calculated by calculating the max error between the interpolant at a given 
! degreee (between nDeg and nDeg+nTests)
!

 USE ModelPrecision
 USE ConstantsDictionary
 USE CommonRoutines
 USE Legendre
 USE Lagrange_1D_Class

 IMPLICIT NONE
 
 INTEGER, PARAMETER :: nDeg = 2
 INTEGER, PARAMETER :: nDegRef = 40
 INTEGER, PARAMETER :: nTests = 15
 INTEGER, PARAMETER :: nPlot = 50
 CHARACTER(13), PARAMETER :: l2errFMT = '(I3,2x,E17.7)'
 TYPE(Lagrange_1D) :: thisInterp, plotInterp, refInterp

 
 REAL(prec), ALLOCATABLE :: testF(:), testDFDS(:)
 REAL(prec), ALLOCATABLE :: s(:), qWeights(:), dMat(:,:)
 REAL(prec), ALLOCATABLE :: thisToRef(:,:)
 REAL(prec), ALLOCATABLE :: thisToPlot(:,:)

 
 REAL(prec) :: L2Err(0:nTests), L2ErrDer(0:nTests)
 REAL(prec) :: sRef(0:nDegRef), qWeightRef(0:nDegRef)
 REAL(prec) :: actualF(0:nDegRef), actualDFDS(0:nDegRef) 
 REAL(prec) :: mappedF(0:nDegRef), mappedDFDS(0:nDegRef)
 REAL(prec) :: plotF(0:nPlot), plotDFDS(0:nPlot), sPlot(0:nPlot)

 INTEGER :: iTest, thisDeg, fUnit
 
 CHARACTER(5) :: degChar

 
    ! Build the reference quadrature grid for computing the L-2 norm
    CALL LegendreGaussQuad( nDegRef, sRef, qWeightRef )  
  
    ! Build the reference interpolant
    CALL refInterp % Build( nDegRef, sRef )
    
    ! Calculate the actual function on the reference grid
    CALL TestFunction( nDegRef, sRef, actualF, actualDFDS )
    
    
    ! Build the plotting grid
    CALL GenerateUniformGrid( nPlot, -ONE, ONE, sPlot )
    
    ! Build the plotting interpolant
    CALL plotInterp % Build( nPlot, sPlot )
 
    DO iTest = 0, nTests
    
       ! Calculate the polynomial degree for this test
       thisDeg = nDeg + iTest
       write(degChar,'(I5.5)') thisDeg
       
       
       ! Allocate space 
       ALLOCATE( testF(0:thisDeg), testDFDS(0:thisDeg), s(0:thisDeg), qWeights(0:thisDeg) )
       
       
       ! Generate the Legendre-points and quadrature weights for this interpolant
       CALL LegendreGaussQuad( thisDeg, s, qWeights ) 
       
       
       ! Build this interpolant
       CALL thisInterp % Build( thisDeg, s )
       
       CALL thisInterp % CalculateDerivativeMatrix( dMat ) ! Calculate the Derivative matrix
       
       ! Calculate the interpolation matrix from this grid to the reference grid
       CALL thisInterp % CalculateInterpolationMatrix( refInterp, thisToRef ) 
       
       ! Calculate the interpolation matrix from this grid to the plotting grid
       CALL thisInterp % CalculateInterpolationMatrix( plotInterp, thisToPlot ) 
       
       
       ! Generate the test function 
       CALL TestFunction( thisDeg, s, testF, testDFDS )
       
       ! Calculate dfds through matrix multiplication
       testDFDS(0:thisDeg) = MATMUL( dMat, testF )
       
       ! Interpolate to the reference grid through matrix multiplication
       CALL thisInterp % CoarseToFine( testF, thisToRef, nDegRef, mappedF )  
       CALL thisInterp % CoarseToFine( testDFDS, thisToRef, nDegRef, mappedDFDS )  
       
       ! Interpolate to the plotting grid through matrix multiplication
       CALL thisInterp % CoarseToFine( testF, thisToPlot, nPlot, plotF )
       CALL thisInterp % CoarseToFine( testDFDS, thisToPlot, nPlot, plotDFDS )
       
       
       ! Write the function and the derivative to curve files  
       CALL plotInterp % WriteCurve( plotF, 'f.'//degChar)
       CALL plotInterp % WriteCurve( plotDFDS, 'dfds.'//degChar)
       
       ! Calculate the error norms using the dot product of the square difference with the 
       ! Gauss quadrature weights
       mappedF    = ( mappedF - actualF )**2
       mappedDFDS = ( mappedDFDS - actualDFDS )**2
       
       L2Err(iTest)    = DOT_PRODUCT( mappedF, qWeightRef )
       L2ErrDer(iTest) = DOT_PRODUCT( mappedDFDS, qWeightRef ) 
       
       ! DeAllocate space to prepare for the next test 
       DEALLOCATE( testF, testDFDS, s, qWeights, thisToRef, thisToPlot, dMat )
       
       CALL thisInterp % Trash( )
    
    ENDDO
    
    OPEN( UNIT = NewUnit(fUnit), FILE = 'Lagrange1D_Test_L2Errf.curve')
    WRITE( fUnit, * ) '#L2-Error-f'
    
    DO iTest = 0, nTests    
       thisDeg = nDeg + iTest
       WRITE( fUnit, l2errFMT ) thisDeg, log10(L2Err(iTest))
    ENDDO
    CLOSE( fUnit )
    
    OPEN( UNIT = NewUnit(fUnit), FILE = 'Lagrange1D_Test_L2Errdfds.curve')
    WRITE( fUnit, * ) '#L2-Error-dfds'
    DO iTest = 0, nTests    
       thisDeg = nDeg + iTest
       WRITE( fUnit, l2errFMT ) thisDeg, log10(L2ErrDer(iTest))
    ENDDO

    CLOSE( fUnit )
 
 CONTAINS
 SUBROUTINE TestFunction( n, s , f, dfds )
 !
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   INTEGER, INTENT(in)     :: n
   REAL(prec), INTENT(in)  :: s(0:n)
   REAL(prec), INTENT(out) :: f(0:n), dfds(0:n)
   ! Local
   INTEGER :: j
   
      DO j = 0, n
   
         f(j)    = exp( -s(j)*s(j) )
         dfds(j) = -TWO*s(j)*f(j)
   
      ENDDO
      
 END SUBROUTINE TestFunction
 
 SUBROUTINE GenerateUniformGrid( n, a, b, s )
!
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   INTEGER, INTENT(in)      :: n
   REAL(prec), INTENT(in)   :: a, b
   REAL(prec), INTENT(out)  :: s(0:n)
   ! Local
   INTEGER    :: i
   REAL(prec) :: ds
   
      ds = (b-a)/REAL(n,prec)
      
      
      DO i = 0, n
      
         s(i) = a + ds*REAL(i,prec)
      
      ENDDO
   
 END SUBROUTINE GenerateUniformGrid
 
 
END PROGRAM TestLagrangeInterpolation_1D
