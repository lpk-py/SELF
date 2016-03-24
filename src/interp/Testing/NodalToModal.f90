PROGRAM NodalToModal



USE ModelPrecision
USE ModelFlags
USE CommonRoutines

USE Legendre
USE Lagrange_1D_Class

IMPLICIT NONE

 TYPE(Lagrange_1D)       :: interp
 REAL(prec), ALLOCATABLE :: fNodal(:),fModal(:), fFiltered(:), fPlot(:), fFiltPlot(:)
 REAL(prec), ALLOCATABLE :: s(:), sUni(:)
 REAL(prec), ALLOCATABLE :: q(:)
 REAL(prec), ALLOCATABLE :: P(:,:), V(:,:), Vinv(:,:) ! Filter, modal-projection matrix, and modal projection inverse matrix
 REAL(prec), ALLOCATABLE :: Tplot(:,:)
 REAL(prec) :: Lnorm, L, Ls
 INTEGER :: i, j, N, M, nPlot, fUnit

   OPEN( UNIT = NewUnit(fUnit), FILE = 'N2M.data' )
   READ( fUnit, * ) N, M, nPlot  
   CLOSE( fUnit )

   ALLOCATE( fNodal(0:M), fModal(0:M), fPlot(0:nPlot) )
   ALLOCATE( fFiltered(0:M), fFiltPlot(0:nPlot) )
   ALLOCATE( s(0:M), sUni(0:nPlot) )
   ALLOCATE( q(0:M) )
   ALLOCATE( P(0:M,0:M), V(0:M,0:M), Vinv(0:M,0:M), Tplot(0:nPlot,0:M) )

   ! Generate the Legendre Gauss Nodes
   CALL GenerateLegendreQuadrature( M, s, q, GAUSS )
   sUni = UniformPoints( -ONE, ONE, nPlot )

   ! Construct the mapping matrices
   CALL interp % Build( M, s )
   CALL interp % CalculateInterpolationMatrix( nPlot, sUni, Tplot )

   ! Construct the modal projection matrix
   P = ZERO
   DO i = 0, M
      Lnorm = ZERO
      
      ! Construct the filter matrix as a "box-cutoff matrix"
      IF( i <= N)THEN
         P(i,i) = ONE
      ENDIF 

      DO j = 0, M
         CALL LegendrePolynomial(i, s(j), L, Ls)
         Lnorm = Lnorm + L*L*q(j)
         V(i,j) = L*q(j)
         ! The inverse of the modal projection matrix is easy to build since the Legendre basis
         ! is an orthogonal basis
         Vinv(j,i) = L
      ENDDO
      V(i,0:M) = V(i,0:M)/Lnorm
   ENDDO

   ! Construct the modal coefficients
   fNodal = ZERO
   DO i = 0,M
      IF( s(i) < ZERO )THEN
         fNodal(i) = -ONE
      ELSE
         fNodal(i) = ONE
      ENDIF
     ! fNodal(i) = tanh( s(i)/0.2_prec)
   ENDDO

   fPlot = MATMUL( Tplot, fNodal )

   fModal = MATMUL( V, fNodal ) ! Obtain the modal coefficients 

   OPEN( UNIT=NewUnit(fUnit), FILE='fModal.curve')
   WRITE( fUnit, * ) '#f'
   DO i = 0, M 
      WRITE(fUnit, *) i, fModal(i)**2
   ENDDO


   fModal = MATMUL( P, fModal ) ! Apply the filter
   fFiltered = MATMUL( Vinv, fModal ) ! Apply the inverse of nodal-to-modal (wouldn't that be "modal-to-nodal" :P)

   fFiltPlot = MATMUL( Tplot, fFiltered ) 
   

   OPEN( UNIT=NewUnit(fUnit), FILE='fPlot.curve')
   WRITE( fUnit, * ) '#f'
   DO i = 0, nPlot 
      WRITE(fUnit, *) sUni(i), fPlot(i)
   ENDDO
   WRITE( fUnit, * ) '#f-filtered'
   DO i = 0, nPlot 
      WRITE(fUnit, *) sUni(i), fFiltPlot(i)
   ENDDO

   CLOSE( fUnit )


   CLOSE( fUnit )



   CALL interp % Trash( )

   DEALLOCATE( fNodal, fModal, fFiltered, fPlot, fFiltPlot )
   DEALLOCATE( s, sUni )
   DEALLOCATE( q )
   DEALLOCATE( P, V, Vinv, Tplot )




END PROGRAM NodalToModal
