PROGRAM BAROTROPIC_INITIAL_CONDITION
! 

  
! "The basics"
USE COMMONDATA 
USE RUN_PARAMS_CLASS
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
USE DG2D_RELAXATION_FIELD_CLASS
!
!
! Geometry Modules
USE BARO_MAPPEDGEOM_2D_CLASS
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
 
  integer :: iEl! Loop Counters

  real(prec) :: tn

  character(len=5)  :: pIDchar

  
  



     CALL MPI_INIT ( iError )
      
     CALL MPI_COMM_SIZE( MPI_COMM_WORLD, MPI_SIZE, iError )

     CALL MPI_COMM_RANK( MPI_COMM_WORLD, MPI_RANK, iError )
     
     ! Perform sanity check
     print *, "Greetings from Process :", MPI_RANK, "of", MPI_SIZE

     
     write(pIDchar,'(I5.5)')MPI_RANK  

     CALL MPI_BARRIER( MPI_COMM_WORLD, iError )
  

     ! Build the DGSEM structure
     CALL thisDGsem % BUILD( 3,.TRUE., MPI_COMM_WORLD, MPI_RANK )

 
! ///////////////////// INITIAL CONDITION SETUP //////////////////////// !




!$OMP PARALLEL

!$OMP DO
      do iEl = 1, thisDGsem % mesh % nElems ! cycle over the elements
   
   
         CALL INIT( thisDGsem % sol(iEl), &
                    thisDGSEM % relax(iEl), &
                    thisDGsem % mesh % elements(iEl) % geometry, & 
                    thisDGsem % dgStorage )


      enddo ! iEl, cycle over the elements 
!$OMP END DO

!$OMP END PARALLEL

!//////////////////// WRITE THE INITIAL CONDITIONS TO PICKUP FILE ///////// !


     CALL thisDGSEM % WRITE_PICKUP('State.'//pIDchar, 0 )

!//////////////////// DONE AND NOW CLEARING SPACE ///////////////////////// !
     
     CALL thisDGSEM % TRASH( )
  


 close(unit=300)

 call MPI_FINALIZE ( iError )


CONTAINS

 SUBROUTINE INIT( swSol, rField, geometry, dgStorage )
    TYPE( DG2D_SOLUTION ), intent(inout)         :: swSol
    TYPE( DG2D_RELAXATION_FIELD ), intent(inout) :: rField
    TYPE( BARO_MAPPEDGEOM_2D ), intent(in)       :: geometry
    TYPE( NODAL_STORAGE_2D), intent(in)          :: dgStorage
    ! LOCAL
    integer :: iX, iY, iZ, nX, nY, nZ
    real(prec) :: x, y,z
    real(prec) :: tRamp, phase

       CALL dgStorage % GET_N_NODES( nX, nY )


          do iY = 0,nY
             do iX = 0,nX
    

                x = geometry % x(iX,iY)

                y = (geometry % y(iX,iY) - jetCenter)/hw


                
          
          
                swSol % solInterior(iX,iY,1) = -uMax*( tanh(y)*tanh(y) - 1.0_prec )

                swSol % solInterior(iX,iY,2) = 0.0_prec
                                                     

                swSol % solInterior(iX,iY,3) = -uMax*f0*hw*( tanh(y) + 1.0_prec ) + uMax*f0*hw

                ! Set the relaxation field
                rField % solField(iX,iY,1) = ZERO
                rField % solField(iX,iY,2) = ZERO
                rField % solField(iX,iY,3) = ZERO

                rField % rDrag(iX,iY) = 0.0_prec

             enddo
          enddo


 END SUBROUTINE INIT


END PROGRAM BAROTROPIC_INITIAL_CONDITION
