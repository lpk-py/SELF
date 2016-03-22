PROGRAM TestFilter_1D



USE ModelPrecision
USE ModelFlags

USE NodalStorage_1D_Class
USE AdaptiveFilterStorage1D_Class

IMPLICIT NONE

 TYPE(NodalStorage_1D)          :: dgStorage
 TYPE(AdaptiveFilterStorage_1D) :: filter
 REAL(prec), ALLOCATABLE :: f(:), fPlot(:), sUni(:), s(:)
 REAL(prec), ALLOCATABLE :: Tplot(:,:)
 REAL(prec) :: Lnorm, L, Ls
 INTEGER :: i, j, N, M, nPlot, fUnit

   OPEN( UNIT = NewUnit(fUnit), FILE = 'N2M.data' )
   READ( fUnit, * ) N, M, nPlot  
   CLOSE( fUnit )

   ALLOCATE( f(0:M), fPlot(0:nPlot) )
   ALLOCATE( Tplot(0:nPlot,0:M), sUni(0:nPlot), s(0:M) )

   CALL dgStorage % Build( M, GAUSS, DG )
   ! Generate the uniform points for plotting
   sUni = UniformPoints( -ONE, ONE, nPlot )

   ! Construct the mapping matrices
   CALL dgStorage % interp % CalculateInterpolationMatrix( nPlot, sUni, Tplot )

   ! Construct the filter
   CALL filter % Build( dgStorage, N )
   CALL dgStorage % GetNodes( s )
   
   f = ZERO
   DO i = 0, M
      IF( s(i) < ZERO )THEN
         f(i) = ONE
      ENDIF
   ENDDO

   fPlot = MATMUL( Tplot, f )
   OPEN( UNIT=NewUnit(fUnit), FILE='fPlot.curve')
   WRITE( fUnit, * ) '#f'
   DO i = 0, nPlot 
      WRITE(fUnit, *) sUni(i), fPlot(i)
   ENDDO
   
   f = filter % ApplyFilter( f )
   
   fPlot = MATMUL( Tplot, f )
   WRITE( fUnit, * ) '#f-filtered'
   DO i = 0, nPlot 
      WRITE(fUnit, *) sUni(i), fPlot(i)
   ENDDO

   CLOSE( fUnit )

   CALL dgStorage % Trash( )
   CALL filter % Trash( )


END PROGRAM TestFilter_1D
