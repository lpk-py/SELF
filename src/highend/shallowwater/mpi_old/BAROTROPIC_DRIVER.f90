PROGRAM BAROTROPIC_DRIVER
!
!
!


  
! "The basics"
USE COMMONDATA 
!
!
! Interpolation modules
USE LAGRANGE_1D_CLASS
USE LAGRANGE_2D_CLASS
!
!
! Storage Container modules
USE NODAL_STORAGE_2D_CLASS
USE DG2D_SOLUTION_STORAGE_CLASS
!
!
! Geometry Modules
USE MAPPEDGEOM_2D_CLASS
USE BARO_MAPPEDGEOM_CLASS
!
!
! Mesh Primitives
USE GEOMETRY_BASICS
USE BARO_QUADMESH_CLASS
!
!
! Barotropic Flow Data type
USE DGSEM_BAROTROPIC_CLASS

! Parallel libraries
USE OMP_LIB


  IMPLICIT NONE
  
  ! OPENMPI- declarations
  INCLUDE 'mpif.h'
  integer :: iError
  integer :: MPI_COMM
  integer :: MPI_RANK, MPI_SIZE


  ! Set the bathymetry gradient direction with these parameters
  TYPE( DGSEM_BAROTROPIC ) :: thisDGsem
  TYPE( LAG_INTERP2D )     :: plotInterp
 
  integer ::  iT, iS ! Loop Counters
  integer ::  nP, pDeg, iter0, nT, dFreq

  real(prec), allocatable :: sPlot(:)
  real(prec), allocatable :: Tx(:,:), Ty(:,:)
  
  real(prec) :: tn, dxP, deltaT

  character(len=10) :: iterChar
  character(len=5)  :: pIDchar
  integer :: numThreads
  real :: cpuStart, cpuEnd
  



     CALL MPI_INIT ( iError )
      
     CALL MPI_COMM_SIZE( MPI_COMM_WORLD, MPI_SIZE, iError )

     CALL MPI_COMM_RANK( MPI_COMM_WORLD, MPI_RANK, iError )
     
     ! Perform sanity check
     print *, "Greetings from Process :", MPI_RANK, "of", MPI_SIZE

     
     write(pIDchar,'(I5.5)')MPI_RANK  

     CALL MPI_BARRIER( MPI_COMM_WORLD, iError )
  

     CALL thisDGsem % BUILD(  3, .FALSE., MPI_COMM_WORLD, MPI_RANK )

     
     
   
!/////////////////////////////////////////////////////////////////!
!---------------- Set up the plot interpolants -------------------!
!/////////////////////////////////////////////////////////////////!

     nP = thisDGsem % params % nPlot
     pDeg = thisDGsem % params % polyDeg
     dxP = thisDGsem % params % dxPlot
 
     ALLOCATE( sPlot(0:nP), Tx(0:nP,0:pDeg), Ty(0:nP,0:pDeg) )

     !Use a uniform grid in the mapped (computational) domain 
     do iS = 0, nP

        sPlot(iS) = -1.0_prec + real(iS, prec)*dxP

     enddo


     CALL plotInterp % BUILD_WITHNODES( sPlot, sPlot )


      ! Build the interpolation matrices to interpolate onto uniform grid
      CALL PINTERP_MATRIX_1D( thisDGsem % dGStorage % interp % xInterp, &
                              plotInterp % xInterp, Tx )

      CALL PINTERP_MATRIX_1D( thisDGsem % dGStorage % interp % yInterp, &
                              plotInterp % yInterp, Ty )  


!/////////////////////////////////////////////////////////////////!
!---------------- Done with plot interpolant setup ---------------!
!/////////////////////////////////////////////////////////////////!

      iter0 = thisDGSEM % params % iterInit

      !  Write out the initial condition
  
       write(iterChar,'(I10.10)') iter0
       if( nP == pDeg )then
               
          CALL thisDGsem % WRITE_TECPLOT_ASIS( 'State.'//pIDchar//'.'//iterChar )

       else

          CALL thisDGsem % WRITE_TECPLOT(  Tx, Ty, pDeg, nP, plotInterp, &
                                                     'State.'//pIDchar//'.'//iterChar )

       endif


! ///////////////////// Time Stepping //////////////////////// !
       


      !CALL CPU_TIME( cpuStart )
      nT = thisDGSEM % params % nTimeSteps
      dFreq = thisDGSEM % params % dumpFreq
      deltaT = thisDGSEM % params % dt
      
      do iT = iter0, iter0+nT ! Loop over time-steps

          tn = real( iT, prec )*deltaT ! time at the beginning of this timestep

          CALL thisDGSEM % MAPPED_DGSEM_BAROTROPIC_RK3( tn, MPI_COMM_WORLD, MPI_RANK )


          if( mod( iT+1, dFreq) == 0 ) then

             write(iterChar,'(I10.10)') iT+1

             if( nP == pDeg )then
               
                CALL thisDGsem % WRITE_TECPLOT_ASIS( 'State.'//pIDchar//'.'//iterChar )

             else

                CALL thisDGsem % WRITE_TECPLOT(  Tx, Ty, pDeg, nP, plotInterp, &
                                               'State.'//pIDchar//'.'//iterChar )

             endif

                

          endif

       enddo
     
       CALL thisDGSEM % WRITE_PICKUP('State.'//pIDchar, iTer0+nT )

!//////////////////// DONE AND NOW CLEARING SPACE ///////////////////////// !
     
     CALL thisDGSEM % TRASH( )

     CALL plotInterp % TRASH( )

     DEALLOCATE( sPlot, Tx, Ty )
  


 close(unit=300)

 call MPI_FINALIZE ( iError )



END PROGRAM BAROTROPIC_DRIVER
