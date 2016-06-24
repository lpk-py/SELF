
PROGRAM TestCoarseToFine_2D


 USE ModelPrecision
 USE ConstantsDictionary
 USE CommonRoutines
 USE Timing
 USE Legendre
 USE Lagrange_2D_Class

 IMPLICIT NONE
 
 INTEGER, PARAMETER :: nDeg = 10
 INTEGER, PARAMETER :: nPlot = 20
 INTEGER, PARAMETER :: nCycles = 100000

 TYPE(MultiTimers) :: timers
 TYPE(Lagrange_2D) :: thisInterp, plotInterp
 REAL(prec) :: s(0:nDeg), qWeight(0:nDeg)
 REAL(prec) :: F(0:nDeg,0:nDeg)
 REAL(prec) :: plotF(0:nPlot,0:nPlot), sPlot(0:nPlot)
 REAL(prec) :: thisToPlotS(0:nPlot,0:nDeg), thisToPlotP(0:nPlot,0:nDeg)
 INTEGER    :: i

   CALL timers % Build( )
   CALL timers % AddTimer( 'CoarseToFine_2D', 1 )
   
   CALL LegendreGauss( nDeg, s, qWeight )  
   sPlot = UniformPoints( -ONE, ONE, nPlot )

   ! Build this interpolant
   CALL thisInterp % Build( nDeg, nDeg, s )
   CALL plotInterp % Build( nPlot, nPlot, sPlot )
       
   ! Calculate the interpolation matrix from this grid to the reference grid
   CALL thisInterp % CalculateInterpolationMatrix( nPlot, nPlot, sPlot, sPlot, thisToPlotS, thisToPlotP ) 
      
   F = ZERO
   CALL timers % StartThisTimer( 1 )
   DO i = 1, nCycles
      ! Interpolate to the reference grid through matrix multiplication
      CALL thisInterp % CoarseToFine( F, thisToPlotS, thisToPlotP, nPlot, nPlot, plotF )   
   ENDDO
   CALL timers % StopThisTimer( 1 )

   CALL timers % Write_MultiTimers( )

   CALL timers % Trash( )
   CALL thisInterp % Trash( )
   CALL plotInterp % Trash( )
    


 
END PROGRAM TestCoarseToFine_2D
