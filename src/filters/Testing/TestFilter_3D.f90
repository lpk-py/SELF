PROGRAM TestFilter_3D



USE ModelPrecision
USE ModelFlags

USE NodalStorage_3D_Class
USE RollOffFilter3D_Class

IMPLICIT NONE

 TYPE(NodalStorage_3D)   :: dgStorage
 TYPE(RollOffFilter3D)   :: filter
 REAL(prec), ALLOCATABLE :: f(:,:,:), fPlot(:,:,:), fFiltered(:,:,:), sUni(:), s(:), p(:), q(:)
 REAL(prec), ALLOCATABLE :: T1plot(:,:), T2plot(:,:), T3plot(:,:) 
 REAL(prec) :: Lnorm, L, Ls
 INTEGER :: i, j, k, N, M, nPlot, fUnit
 CHARACTER(5) :: zoneID

   OPEN( UNIT = NewUnit(fUnit), FILE = 'N2M.data' )
   READ( fUnit, * ) N, M, nPlot  
   CLOSE( fUnit )

   ALLOCATE( f(0:M,0:M,0:M), fPlot(0:nPlot,0:nPlot,0:nPlot), fFiltered(0:nPlot,0:nPlot,0:nPlot) )
   ALLOCATE( sUni(0:nPlot), s(0:M), p(0:M), q(0:M) )
   ALLOCATE( T1plot(0:nPlot,0:M), T2plot(0:nPlot,0:M), T3plot(0:nPlot,0:M) )

   CALL dgStorage % Build( M, M, M, GAUSS, DG )
   ! Generate the uniform points for plotting
   sUni = UniformPoints( -ONE, ONE, nPlot )

   ! Construct the mapping matrices
   CALL dgStorage % interp % CalculateInterpolationMatrix( nPlot, nPlot, nPlot, & 
                                                           sUni, sUni, sUni, &
                                                           T1plot, T2plot, T3plot )

   ! Construct the filter
   CALL filter % Build( dgStorage, N )
   CALL dgStorage % GetNodes( s, p, q )
   
   f = ZERO
   DO k = 0, M
      DO j = 0, M
         DO i = 0, M
            IF( s(i) < ZERO .AND. p(j) < ZERO .AND. q(k) < ZERO )THEN
               f(i,j,k) = -ONE
            ELSE
               f(i,j,k) = ONE
            ENDIF
         ENDDO
      ENDDO
   ENDDO

   CALL dgStorage % interp % CoarseToFine( f, T1plot, T2plot, T3plot, nPlot, nPlot, nPlot, fPlot ) 
   f = filter % ApplyFilter( f )
   CALL dgStorage % interp % CoarseToFine( f, T1plot, T2plot, T3plot, nPlot, nPlot, nPlot, fFiltered ) 
 
   OPEN( UNIT=NEWUNIT(fUnit), &
         FILE= 'Filter3Dtest.tec', &
         FORM='formatted', &
         STATUS='replace')  
    
   WRITE(fUnit,*) 'VARIABLES = "X", "Y","Z", "F", "F-Filtered" '
 
        
   WRITE(zoneID,'(I5.5)') 0
   WRITE(fUnit,*) 'ZONE T="el'//trim(zoneID)//'", I=',nPlot+1,', J=', nPlot+1,', K=', nPlot+1,',F=POINT'

   DO k = 0, nPlot
      DO j = 0, nPlot
         DO i = 0, nPlot
            WRITE (fUnit,*) sUni(i), sUni(j), sUni(k), &
                            fPlot(i,j,k), fFiltered(i,j,k)
         ENDDO
      ENDDO
   ENDDO

   CLOSE(UNIT=fUnit)


   DEALLOCATE( f, fPlot, fFiltered )
   DEALLOCATE( sUni, s, p, q )
   DEALLOCATE( T1plot, T2plot, T3plot )

   CALL dgStorage % Trash( )
   CALL filter % Trash( )


END PROGRAM TestFilter_3D
