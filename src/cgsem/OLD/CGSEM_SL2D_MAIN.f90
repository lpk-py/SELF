PROGRAM CGSEM_SL2D_MAIN

! 



! The model settings
USE MODEL_PRECISION
USE CONSTANTS_DICTIONARY
USE MODEL_FLAGS
!
 USE COMMONROUTINES
 USE RUN_PARAMS_CLASS
 USE MATRIXROUTINES

 USE LEGENDRE
 USE LAGRANGE_1D_CLASS
 USE LAGRANGE_2D_CLASS

 USE QUADMESH_CLASS
 USE QUADELEMENT_CLASS
 USE MAPPEDGEOM_2D_CLASS
 USE GEOMETRY_BASICS

 USE CGSEM2D_CLASS
 USE CGSEM_SL2D_CLASS
 
IMPLICIT NONE

  TYPE( CGSEM_SL2D )     :: thisCGSM
  TYPE( LAG_INTERP2D )   :: plotInterp, low
  

  real(prec), allocatable :: Tmat(:,:)
  real(prec), allocatable :: plotNodes(:)  
  real(prec)              :: thisX, thisY
  real(prec),allocatable  :: speeds(:,:)
  real(prec)              :: speedDiff, N, h, cActual

  integer :: nPlot, pDeg, iS, iP, iEl, iFree, iTest, nModes
  integer :: fUnit

  logical :: fileExists
  character(3) :: testChar, modeChar, pDegChar


 
   ! Setup


      ! Read in the model parameters from runtime.params
      CALL thisCGSM % BUILD( )

      nModes = thisCGSM % params % nModes

      nPlot = thisCGSM % params % nPlot

      thisCGSM % params % polyDeg = thisCGSM % params % polyDeg

      pDeg = thisCGSM % params % polyDeg

      write(pDegChar,'(I3.3)')pDeg

      ! Allocate the quadrature to plot mesh matrix
      ALLOCATE( Tmat(0:nPlot,0:pDeg) )
 
      ALLOCATE( plotNodes(0:nPlot) )


      ! Build the uniform plotting nodes
      do iS = 0, nPlot

         plotNodes(iS) = -ONE + (thisCGSM % params % dxPlot)*real(iS, prec)

      enddo


      CALL plotInterp % BUILD_WITHNODES( plotNodes, plotNodes )


   
   ! Build the data structure

      CALL thisCGSM % BUILD( )

 
      CALL thisCGSM % mesh % SCALE_THE_MESH( thisCGSM % params % xScale, & 
                                             thisCGSM % params % yScale, &
                                             thisCGSM % cgStorage )
       


      CALL PINTERP_MATRIX_1D(thisCGSM % cgStorage% interp % xInterp, plotInterp % xInterp, Tmat )  

      thisCGSM % source = ZERO 

     
     ! Calculate the eigenmodes
      print*, 'Calculating Eigenmodes...'

      CALL thisCGSM % CALC_EIGENMODES( )

      print*, 'Done!'

      CALL thisCGSM % WRITE_EVALS_TO_CURVE( )

      ! Write this mode to file

      CALL thisCGSM % WRITE_TECPLOT_MODES( nPlot, pDeg, plotInterp, Tmat, thisCGSM % params % nModes )

      print*, ' Writing the mesh to tecplot file.... '
      CALL thisCGSM % mesh % WRITE_TECPLOT( 'shelf'  )



   ! ------------------------------------------------------------------------ !
   ! Clean up memory

     DEALLOCATE( Tmat )
 
     DEALLOCATE( plotNodes  )


     CALL thisCGSM % TRASH( )

     CALL plotInterp % TRASH( )



 CONTAINS

 SUBROUTINE CALCULATE_RMS_DIFF( lowModes, highModes, lowToHigh, lowInt, highInt, polyDeg, nEl, error )

  integer, intent(in)              :: polyDeg, nEl
  real(prec), intent(in)           :: lowModes(0:polyDeg-1, 0:polyDeg-1, 1:nEl, 1:nModes)
  real(prec), intent(in)           :: highModes(0:polyDeg, 0:polyDeg, 1:nEl, 1:nModes)
  real(prec), intent(in)           :: lowToHigh(0:polyDeg, 0:polyDeg-1)
  TYPE( LAG_INTERP2D ), intent(in) :: lowInt, highInt
  real(prec), intent(inout)        :: error(1:nModes)
  ! LOCAL
  integer :: iS, iP, iEl, iMode
  real(prec) :: thisErr, corr, thisMag
  real(prec) :: thisMode(0:polyDeg, 0:polyDeg)

!$OMP PARALLEL
!$OMP DO
     do iMode = 1, nModes

        thisErr = ZERO
        corr = ZERO
        thisMag = ZERO


        !  Because the eigenvector/value software may produce the same structure with opposite sign
        !  we need to adjust the "difference" calculation by computing (thisMode - thatMode*sign(correlation) )
        !  where "sign(correlation)" is the sign of the correlation between the low Res and high Res modes.
        ! We first compute the correlation...
          
        do iEl = 1, nEl
           
           ! Interpolate the solutions onto a uniform mesh
           CALL COARSE_TO_FINE_2D( lowInt, lowModes(:,:,iEl,iMode), & 
                                   highInt, thisMode, lowToHigh, lowToHigh)

           do iP = 0, pDeg
              do iS = 0, pDeg

                  corr = corr + thisMode(iS,iP)*highModes(iS,iP,iEl,iMode) 
                  thisMag = thisMag + thisMode(iS,iP)*thisMode(iS,iP)

               enddo
            enddo

         enddo
       
        do iEl = 1, nEl
           
           ! Interpolate the solutions onto a uniform mesh
           CALL COARSE_TO_FINE_2D( lowInt, lowModes(:,:,iEl,iMode), & 
                                   highInt, thisMode, lowToHigh, lowToHigh)

           do iP = 0, pDeg
              do iS = 0, pDeg

                  thisErr = thisErr + (thisMode(iS,iP) - SIGN(ONE,corr)*highModes(iS,iP,iEl,iMode) )**2

               enddo
            enddo

         enddo

         error(iMode) = sqrt(thisErr)/sqrt(thisMag)

      enddo
!$OMP END DO
!$OMP END PARALLEL

 END SUBROUTINE CALCULATE_RMS_DIFF


END PROGRAM CGSEM_SL2D_MAIN
