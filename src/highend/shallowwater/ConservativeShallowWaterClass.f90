MODULE ConservativeShallowWaterClass
! ShallowWaterClass.f90
!
! schoonover.numerics@gmail.com
! 
! o (ver 1.0) March 2014
! o (ver 2.1) Dec 2015
!
!
! 
!  
! ================================================================================================ !

! src/common/
USE ModelPrecision
USE ConstantsDictionary
USE ModelFlags
USE CommonRoutines
! src/interp/
USE Lagrange_2D_Class
! src/nodal/
USE NodalStorage_2D_Class
USE DGSEMSolutionStorageClass_2D
! src/geometry/
USE EdgeClass
USE QuadElementClass  
USE QuadMeshClass    
! src/highend/shallowwater/
USE CSWParamsClass
! Nocturnal Aviation classes and extensions
!USE FTTimerClass
!USE TIMING


! Parallel Libraries
USE OMP_LIB


IMPLICIT NONE
!
! Some parameters that are specific to this module

  INTEGER, parameter      :: nSWeq = 3  ! The number of prognostic variables for the shallow water equations
  INTEGER, parameter      :: nSWtimers = 5 ! the number of model timers
  REAL(prec), parameter   :: hDefault = ONE ! The default depth
     

! ================================================================================================ !
!     The ShallowWater is based on the DGSEM_2D framework with the addition of 
!     the AddOns data. Also, a few additional routines are added to facilitate
!     the calculation of "source" terms which include the coriolis force, relaxation terms, etc.
! ================================================================================================ !  

    TYPE ShallowWater
      INTEGER                            :: nEq, nPlot, nS, nP
      TYPE( QuadMesh )                   :: mesh
      TYPE( NodalStorage_2D )            :: dGStorage
      TYPE( DGSEMSolution_2D ), ALLOCATABLE :: sol(:)
      !TYPE( DGSEMSolution_2D ), ALLOCATABLE :: tracer(:)
      TYPE( DGSEMSolution_2D ), ALLOCATABLE :: relax(:)
      TYPE( DGSEMSolution_2D ), ALLOCATABLE :: bathymetry(:)
      TYPE( DGSEMSolution_2D ), ALLOCATABLE :: vorticity(:)
      TYPE( SWParams )                        :: params
      REAL(prec), ALLOCATABLE                :: plMatS(:,:), plMatP(:,:)
      REAL(prec)                              :: relativeFactor

      CONTAINS

      PROCEDURE :: Build => Build_ShallowWater
      PROCEDURE :: Trash => Trash_ShallowWater
      PROCEDURE :: BuildQuadMesh => BuildQuadMesh_ShallowWater
      
      ! DGSEMSolution_2DStorage_Class Wrapper Routines
      PROCEDURE :: GetSolution => GetSolution_ShallowWater
      PROCEDURE :: SetSolution => SetSolution_ShallowWater
      PROCEDURE :: GetSolutionAtNode => GetSolutionAtNode_ShallowWater
      PROCEDURE :: SetSolutionAtNode => SetSolutionAtNode_ShallowWater
      PROCEDURE :: GetSolutionAtNodeWithVarID => GetSolutionAtNodeWithVarID_ShallowWater
      PROCEDURE :: SetSolutionAtNodeWithVarID => SetSolutionAtNodeWithVarID_ShallowWater
      PROCEDURE :: GetTendency => GetTendency_ShallowWater
      PROCEDURE :: SetTendency => SetTendency_ShallowWater
      PROCEDURE :: GetTendencyAtNode => GetTendencyAtNode_ShallowWater
      PROCEDURE :: SetTendencyAtNode => SetTendencyAtNode_ShallowWater
      PROCEDURE :: GetTendencyWithVarID => GetTendencyWithVarID_ShallowWater
      PROCEDURE :: SetTendencyWithVarID => SetTendencyWithVarID_ShallowWater
      PROCEDURE :: GetTendencyAtNodeWithVarID => GetTendencyAtNodeWithVarID_ShallowWater
      PROCEDURE :: SetTendencyAtNodeWithVarID => SetTendencyAtNodeWithVarID_ShallowWater   
      PROCEDURE :: GetBoundarySolutionAtBoundary => GetBoundarySolutionAtBoundary_ShallowWater
      PROCEDURE :: SetBoundarySolutionAtBoundary => SetBoundarySolutionAtBoundary_ShallowWater
      PROCEDURE :: GetBoundaryFluxAtBoundary => GetBoundaryFluxAtBoundary_ShallowWater
      PROCEDURE :: SetBoundaryFluxAtBoundary => SetBoundaryFluxAtBoundary_ShallowWater
      PROCEDURE :: GetBoundaryFluxAtBoundaryNode => GetBoundaryFluxAtBoundaryNode_ShallowWater
      PROCEDURE :: SetBoundaryFluxAtBoundaryNode => SetBoundaryFluxAtBoundaryNode_ShallowWater
      ! Extensions of the DGSEMSolutionStorage_Class Wrapper Routines
      PROCEDURE :: GetVelocity => GetVelocity_ShallowWater
      PROCEDURE :: SetVelocity => SetVelocity_ShallowWater
      PROCEDURE :: GetVelocityAtBoundary => GetVelocityAtBoundary_ShallowWater
      PROCEDURE :: SetVelocityAtBoundary => SetVelocityAtBoundary_ShallowWater
      PROCEDURE :: GetRelaxationField => GetRelaxationField_ShallowWater
      PROCEDURE :: SetRelaxationField => SetRelaxationField_ShallowWater
      PROCEDURE :: GetRelaxationFieldAtNode => GetRelaxationFieldAtNode_ShallowWater
      PROCEDURE :: SetRelaxationFieldAtNode => SetRelaxationFieldAtNode_ShallowWater
      PROCEDURE :: GetRelaxationTimeScale => GetRelaxationTimeScale_ShallowWater
      PROCEDURE :: SetRelaxationTimeScale => SetRelaxationTimeScale_ShallowWater
      PROCEDURE :: GetRelaxationTimeScaleAtNode => GetRelaxationTimeScaleAtNode_ShallowWater
      PROCEDURE :: SetRelaxationTimeScaleAtNode => SetRelaxationTimeScaleAtNode_ShallowWater
      PROCEDURE :: GetBathymetry => GetBathymetry_ShallowWater
      PROCEDURE :: SetBathymetry => SetBathymetry_ShallowWater
      PROCEDURE :: GetBathymetryAtNode => GetBathymetryAtNode_ShallowWater
      PROCEDURE :: SetBathymetryAtNode => SetBathymetryAtNode_ShallowWater
      PROCEDURE :: GetBathymetryAtBoundary => GetBathymetryAtBoundary_ShallowWater
      PROCEDURE :: SetBathymetryAtBoundary => SetBathymetryAtBoundary_ShallowWater
      PROCEDURE :: GetPlanetaryVorticity => GetPlanetaryVorticity_ShallowWater
      PROCEDURE :: SetPlanetaryVorticity => SetPlanetaryVorticity_ShallowWater
      PROCEDURE :: GetPlanetaryVorticityAtNode => GetPlanetaryVorticityAtNode_ShallowWater
      PROCEDURE :: SetPlanetaryVorticityAtNode => SetPlanetaryVorticityAtNode_ShallowWater
      PROCEDURE :: GetPlanetaryVorticityAtBoundary => GetPlanetaryVorticityAtBoundary_ShallowWater
      PROCEDURE :: SetPlanetaryVorticityAtBoundary => SetPlanetaryVorticityAtBoundary_ShallowWater
      
       ! Type Specific Routines
      PROCEDURE :: GlobalTimeDerivative => GlobalTimeDerivative_ShallowWater
      PROCEDURE :: ForwardStepRK3 => ForwardStepRK3_ShallowWater
      PROCEDURE :: EdgeFlux => EdgeFlux_ShallowWater
      PROCEDURE :: MappedTimeDerivative => MappedTimeDerivative_ShallowWater
      PROCEDURE :: CalculateSolutionAtBoundaries => CalculateSolutionAtBoundaries_ShallowWater
      PROCEDURE :: CalculateBathymetryAtBoundaries => CalculateBathymetryAtBoundaries_ShallowWater
      PROCEDURE :: CalculateBathymetryGradient => CalculateBathymetryGradient_ShallowWater

      
      PROCEDURE :: CoarseToFine => CoarseToFine_ShallowWater
      PROCEDURE :: WriteTecplot => WriteTecplot_ShallowWater
      PROCEDURE :: WritePickup => WritePickup_ShallowWater
      PROCEDURE :: ReadPickup => ReadPickup_ShallowWater

    END TYPE ShallowWater



 CONTAINS
!
!
!==================================================================================================!
!------------------------------- Manual Constructors/Destructors ----------------------------------!
!==================================================================================================!
!
!
 SUBROUTINE Build_ShallowWater( myDGSEM, forceNoPickup )
 ! S/R Build
 ! 
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS(ShallowWater), INTENT(inout) :: myDGSEM
   INTEGER, INTENT(in), OPTIONAL      :: forceNoPickup
   !LOCAL
   INTEGER :: iEl, globElID, pID2, nS, nP, nPlot
   INTEGER :: rStartRel, rEndRel
   CHARACTER(10) :: pickupIterChar
   CHARACTER(40) :: meshFile
   REAL(prec), ALLOCATABLE :: h(:,:,:), sNew(:)

      CALL myDGSEM % params % Build( )
     
      myDGSEM % relativeFactor = ZERO
      IF( myDGSEM % params % ModelFormulation == SKEW_SYMMETRIC )THEN
         myDGSEM % relativeFactor = ONE
      ENDIF
      
      nS = myDGSEM % params % polyDeg
      nP = nS
      myDGSEM % nS      = nS
      myDGSEM % nP      = nP
      myDGSEM % nEq     = nSWEq
      myDGSEM % nPlot   = myDGSEM % params % nPlot
      
      nPlot = myDGSEM % params % nPlot
      !myDGSEM % nTracer = myDGSEM % params % nTracers
      
      ALLOCATE( h(0:nS,0:nP,3), sNew(0:nPlot) )
      h = hDefault
      h(:,:,2) = ZERO ! dhdx
      h(:,:,3) = ZERO ! dhdy
      CALL myDGSEM % dGStorage % Build( nS, nP, GAUSS, DG )
 
      CALL myDGSEM % BuildQuadMesh( )
      
      ! Set up the solution, relaxation fields, bathymetry, and vorticity
      ALLOCATE( myDGSEM % sol(1:myDGSEM % mesh % nElems) )
      ALLOCATE( myDGSEM % relax(1:myDGSEM % mesh % nElems) )
      ALLOCATE( myDGSEM % bathymetry(1:myDGSEM % mesh % nElems) )
      ALLOCATE( myDGSEM % vorticity(1:myDGSEM % mesh % nElems) )

      ! Build and initialize the solution and the relaxation fields to zero
      DO iEl = 1, myDGSEM % mesh % nElems
         CALL myDGSEM % sol(iEl) % Build( nS, nP, nSWEq )
         CALL myDGSEM % relax(iEl) % Build( nS, nP, nSWEq+1 ) 
         CALL myDGSEM % bathymetry(iEl) % Build( nS, nP, 3 )
         CALL myDGSEM % vorticity(iEl) % Build( nS, nP, 2 ) ! Contains both planetary and relative vorticity
         
         ! Initialize the bathymetry to the default value (to be safe)
         CALL myDGSEM % SetBathymetry( iEl, h )
      ENDDO

      
      ! In the event that this is a pickup run, we'll READ the pickup file here for the solution 
      ! and the relaxation fields. A negative initial iterate can be used to specify to start 
      ! from zeros.
      
      IF( .NOT. PRESENT(forceNoPickup) )then
         
         ! This call reads the solution, the addons and the relaxation-parameter
         CALL myDGSEM % ReadPickup( myDGSEM % params % iterInit )

         DO iEl = 1, myDGSEM % mesh % nElems
            CALL myDGSEM %  CalculateBathymetryGradient( iEl ) 
         ENDDO
      ENDIF
      
      DEALLOCATE(h)
      
      nPlot = myDGSEM % params % nPlot
      myDGSEM % nPlot = nPlot
      
      sNew = UniformPoints( -ONE, ONE, nPlot )
      ALLOCATE( myDGSEM % plMatS(0:nPlot,0:nS), myDGSEM % plMatP(0:nPlot,0:nP) )
      ! Build the plotting matrix
      CALL myDGSEM % dgStorage % interp % CalculateInterpolationMatrix( nPlot, nPlot, sNew, sNew, &
                                                                        myDGSEM % plMatS, &
                                                                        myDGSEM % plMatP )
                                                                        
      DEALLOCATE(sNew)
      
 END SUBROUTINE Build_ShallowWater
!
!
!
 SUBROUTINE Trash_ShallowWater( myDGSEM )
 ! S/R Trash
 ! 
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS(ShallowWater), INTENT(inout) :: myDGSEM
   ! LOCAL
   INTEGER :: iEl

     DO iEl = 1, myDGSEM % mesh % nElems
        CALL myDGSEM % sol(iEl) % Trash( )
        CALL myDGSEM % relax(iEl) % Trash( )
        CALL myDGSEM % bathymetry(iEl) % Trash( )
        CALL myDGSEM % vorticity(iEl) % Trash( )
     ENDDO

     CALL myDGSEM % dGStorage % Trash( )
     CALL myDGSEM % mesh % Trash( )

     DEALLOCATE( myDGSEM % sol ) 
     DEALLOCATE( myDGSEM % relax )
     DEALLOCATE( myDGSEM % bathymetry )
     DEALLOCATE( myDGSEM % vorticity )
     DEALLOCATE( myDGSEM % plMatS, myDGSEM % plMatP )

 END SUBROUTINE Trash_ShallowWater
!
!
!
 SUBROUTINE BuildQuadMesh_ShallowWater( myDGSEM )
 ! S/R BuildQuadMesh
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS( ShallowWater ), INTENT(inout) :: myDGSEM
   ! LOCAL
   
      PRINT*,'Module ShallowWaterClass.f90 : S/R BuildQuadMesh :'
      IF( TRIM( myDGSEM % params % SpecMeshFile ) == nada )THEN
         PRINT*,' Loading default mesh.'
         CALL myDGSEM % mesh % LoadDefaultMesh( myDGSEM % dgStorage % interp, &
                                                myDGSEM % params % nXelem, &
                                                myDGSEM % params % nYelem )
      ELSE
      ! Builds the lateral mesh
         PRINT*, 'Reading mesh from '//trim(myDGSEM % params % SpecMeshFile)//'.'
         CALL myDGSEM % mesh % ReadSpecMeshFile( myDGSEM % dgStorage % interp, &
                                                 myDGSEM % params % SpecMeshFile )
      ENDIF
      
      CALL myDGSEM % mesh % ScaleTheMesh( myDGSEM % dgStorage % interp, &
                                          myDGSEM % params % xScale, &
                                          myDGSEM % params % yScale )
 END SUBROUTINE BuildQuadMesh_ShallowWater
!
!
!==================================================================================================!
!--------------------------------------- Accessors ------------------------------------------------!
!==================================================================================================!
!
!
!
! ---------------------------------------- Solution ---------------------------------------------- !
! 
 SUBROUTINE GetSolution_ShallowWater( myDGSEM, iEl, theSolution  )
 ! S/R GetSolution
 !  
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS(ShallowWater), INTENT(in) :: myDGSEM
   INTEGER, INTENT(in)             :: iEl
   REAL(prec), INTENT(out)         :: theSolution(0:myDGSEM % nS, 0:myDGSEM % nP, 1:myDGSEM % nEq)

      CALL myDGSEM % sol(iEl) % GetSolution( theSolution )

 END SUBROUTINE GetSolution_ShallowWater
!
!
!
 SUBROUTINE SetSolution_ShallowWater( myDGSEM, iEl, theSolution  )
 ! S/R SetSolution
 !  
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS(ShallowWater), INTENT(inout) :: myDGSEM
   INTEGER, INTENT(in)                :: iEl
   REAL(prec), INTENT(in)             :: theSolution(0:myDGSEM % nS, 0:myDGSEM % nP, 1:nSWeq)

      CALL myDGSEM % sol(iEl) % SetSolution( theSolution )

 END SUBROUTINE SetSolution_ShallowWater
!
!
!
 SUBROUTINE GetSolutionAtNode_ShallowWater( myDGSEM, iEl, i, j, theSolution  )
 ! S/R GetSolution
 !  
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS(ShallowWater), INTENT(in) :: myDGSEM
   INTEGER, INTENT(in)             :: iEl
   INTEGER, INTENT(in)             :: i, j
   REAL(prec), INTENT(out)         :: theSolution(1:myDGSEM % nEq)

      CALL myDGSEM % sol(iEl) % GetSolutionAtNode( i, j, theSolution )

 END SUBROUTINE GetSolutionAtNode_ShallowWater
!
!
!
 SUBROUTINE SetSolutionAtNode_ShallowWater( myDGSEM, iEl, i, j, theSolution  )
 ! S/R SetSolution
 !  
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS(ShallowWater), INTENT(inout) :: myDGSEM
   INTEGER, INTENT(in)                :: iEl
   INTEGER, INTENT(in)                :: i, j
   REAL(prec), INTENT(in)             :: theSolution(1:myDGSEM % nEq)

      CALL myDGSEM % sol(iEl) % SetSolutionAtNode( i, j, theSolution )

 END SUBROUTINE SetSolutionAtNode_ShallowWater
!
!
!
 SUBROUTINE GetVelocity_ShallowWater( myDGSEM, iEl, u, v  )
 ! S/R GetVelocity
 !  
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS(ShallowWater), INTENT(in) :: myDGSEM
   INTEGER, INTENT(in)             :: iEl
   REAL(prec), INTENT(out)         :: u(0:myDGSEM % nS,0:myDGSEM % nP)
   REAL(prec), INTENT(out)         :: v(0:myDGSEM % nS,0:myDGSEM % nP)
   ! Local
   REAL(prec) :: h(0:myDGSEM % nS, 0:myDGSEM % nP,1:3)
   
      CALL myDGSEM % sol(iEl) % GetSolutionWithVarID( 1, u )
      CALL myDGSEM % sol(iEl) % GetSolutionWithVarID( 2, v )

      IF( myDGSEM % params % ModelFormulation == CONSERVATIVE )THEN
         CALL myDGSEM % bathymetry(iEl) % GetSolutionWithVarID( 1, h )
         u = u/h(:,:,1)
         v = v/h(:,:,1)
      ENDIF
      
 END SUBROUTINE GetVelocity_ShallowWater
!
!
!
 SUBROUTINE SetVelocity_ShallowWater( myDGSEM, iEl, u, v  )
 ! S/R SetVelocity
 !  
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS(ShallowWater), INTENT(inout) :: myDGSEM
   INTEGER, INTENT(in)                :: iEl
   REAL(prec), INTENT(in)             :: u(0:myDGSEM % nS,0:myDGSEM % nP)
   REAL(prec), INTENT(in)             :: v(0:myDGSEM % nS,0:myDGSEM % nP)

      CALL myDGSEM % sol(iEl) % SetSolutionWithVarID( 1, u )
      CALL myDGSEM % sol(iEL) % SetSolutionWithVarID( 2, v )

 END SUBROUTINE SetVelocity_ShallowWater
!
!
!
 SUBROUTINE GetVelocityAtBoundary_ShallowWater( myDGSEM, iEl, boundary, u, v  )
 ! S/R GetBoundarySolution
 !  
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS(ShallowWater), INTENT(in) :: myDGSEM
   INTEGER, INTENT(in)             :: iEl
   INTEGER, INTENT(in)             :: boundary
   REAL(prec), INTENT(out)         :: u(0:myDGSEM % nS)
   REAL(prec), INTENT(out)         :: v(0:myDGSEM % nS)
   ! Local
   REAL(prec) :: h(0:myDGSEM % nS,1:3)
   
      CALL myDGSEM % sol(iEl) % GetBoundarySolutionAtBoundaryWithVarID( boundary, 1, u )
      CALL myDGSEM % sol(iEl) % GetBoundarySolutionAtBoundaryWithVarID( boundary, 2, v )

      IF( myDGSEM % params % ModelFormulation == CONSERVATIVE )THEN
         CALL myDGSEM % bathymetry(iEl) % GetBoundarySolutionAtBoundaryWithVarID( boundary, 1, h )
         u = u/h(:,1)
         v = v/h(:,1)
      ENDIF
      
 END SUBROUTINE GetVelocityAtBoundary_ShallowWater
!
!
!
 SUBROUTINE SetVelocityAtBoundary_ShallowWater( myDGSEM, iEl, boundary, u, v  )
 ! S/R SetVelocityAtBoundary
 !  
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS(ShallowWater), INTENT(inout) :: myDGSEM
   INTEGER, INTENT(in)                :: iEl
   INTEGER, INTENT(in)                :: boundary
   REAL(prec), INTENT(in)             :: u(0:myDGSEM % nS)
   REAL(prec), INTENT(in)             :: v(0:myDGSEM % nS)
   
      CALL myDGSEM % sol(iEl) % SetBoundarySolutionAtBoundaryWithVarID( boundary, 1, u )
      CALL myDGSEM % sol(iEl) % SetBoundarySolutionAtBoundaryWithVarID( boundary, 2, v )
     
 END SUBROUTINE SetVelocityAtBoundary_ShallowWater
!
!
!
 SUBROUTINE GetSolutionAtNodeWithVarID_ShallowWater( myDGSEM, iEl, i, j, varID, theSolution  )
 ! S/R GetSolution
 !  
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS(ShallowWater), INTENT(in) :: myDGSEM
   INTEGER, INTENT(in)             :: iEl
   INTEGER, INTENT(in)             :: i, j, varID 
   REAL(prec), INTENT(out)         :: theSolution

      CALL myDGSEM % sol(iEl) % GetSolutionAtNodeWithVarID( i, j, varID, theSolution )

 END SUBROUTINE GetSolutionAtNodeWithVarID_ShallowWater
!
!
!
 SUBROUTINE SetSolutionAtNodeWithVarID_ShallowWater( myDGSEM, iEl, i, j, varID, theSolution  )
 ! S/R SetSolution
 !  
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS(ShallowWater), INTENT(inout) :: myDGSEM
   INTEGER, INTENT(in)                :: iEl
   INTEGER, INTENT(in)                :: i, j, varID 
   REAL(prec), INTENT(in)             :: theSolution

      CALL myDGSEM % sol(iEl) % SetSolutionAtNodeWithVarID( i, j, varID, theSolution )

 END SUBROUTINE SetSolutionAtNodeWithVarID_ShallowWater
!
! ----------------------------------------- Tendency --------------------------------------------- !
!
 SUBROUTINE GetTendency_ShallowWater( myDGSEM, iEl, theTend  )
 ! S/R GetTendency
 !  
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS(ShallowWater), INTENT(in) :: myDGSEM
   INTEGER, INTENT(in)             :: iEl
   REAL(prec), INTENT(out)         :: theTend(0:myDGSEM % nS,0:myDGSEM % nP, 1:myDGSEM % nEq)

      CALL myDGSEM % sol(iEl) % GetTendency( theTend )

 END SUBROUTINE GetTendency_ShallowWater
!
!
!
 SUBROUTINE SetTendency_ShallowWater( myDGSEM, iEl, theTend  )
 ! S/R SetTendency
 !  
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS(ShallowWater), INTENT(inout) :: myDGSEM
   INTEGER, INTENT(in)                :: iEl
   REAL(prec), INTENT(in)             :: theTend(0:myDGSEM % nS,0:myDGSEM % nP, 1:nSWeq)

      CALL myDGSEM % sol(iEl) % SetTendency( theTend )   

 END SUBROUTINE SetTendency_ShallowWater
!
!
!
 SUBROUTINE GetTendencyAtNode_ShallowWater( myDGSEM, iEl, i, j, theTend  )
 ! S/R GetTendencyAtNode
 !  
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS(ShallowWater), INTENT(in) :: myDGSEM
   INTEGER, INTENT(in)             :: iEl
   INTEGER, INTENT(in)             :: i, j 
   REAL(prec), INTENT(out)         :: theTend(1:myDGSEM % nEq)

      CALL myDGSEM % sol(iEl) % GetTendencyAtNode( i, j, theTend )

 END SUBROUTINE GetTendencyAtNode_ShallowWater
!
!
!
 SUBROUTINE SetTendencyAtNode_ShallowWater( myDGSEM, iEl, i, j, theTend  )
 ! S/R SetTendencyAtNode
 !  
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS(ShallowWater), INTENT(inout) :: myDGSEM
   INTEGER, INTENT(in)                :: iEl
   INTEGER, INTENT(in)                :: i, j 
   REAL(prec), INTENT(in)             :: theTend(1:myDGSEM % nEq)

      CALL myDGSEM % sol(iEl) % SetTendencyAtNode( i, j, theTend )

 END SUBROUTINE SetTendencyAtNode_ShallowWater
!
!
!
SUBROUTINE GetTendencyWithVarID_ShallowWater( myDGSEM, iEl, varID, theTend  )
 ! S/R GetTendency
 !  
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS(ShallowWater), INTENT(in) :: myDGSEM
   INTEGER, INTENT(in)             :: iEl
   INTEGER, INTENT(in)             :: varID 
   REAL(prec), INTENT(out)         :: theTend(0:myDGSEM % nS,0:myDGSEM % nP)

      CALL myDGSEM % sol(iEl) % GetTendencyWithVarID( varID, theTend )

 END SUBROUTINE GetTendencyWithVarID_ShallowWater
!
!
!
 SUBROUTINE SetTendencyWithVarID_ShallowWater( myDGSEM, iEl, varID, theTend  )
 ! S/R SetTendency
 !  
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS(ShallowWater), INTENT(inout) :: myDGSEM
   INTEGER, INTENT(in)                :: iEl
   INTEGER, INTENT(in)                :: varID 
   REAL(prec), INTENT(in)             :: theTend(0:myDGSEM % nS,0:myDGSEM % nP)

      CALL myDGSEM % sol(iEl) % SetTendencyWithVarID( varID, theTend )

 END SUBROUTINE SetTendencyWithVarID_ShallowWater
!
!
!
 SUBROUTINE GetTendencyAtNodeWithVarID_ShallowWater( myDGSEM, iEl, i, j, varID, theTend  )
 ! S/R GetTendency
 !  
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS(ShallowWater), INTENT(in) :: myDGSEM
   INTEGER, INTENT(in)             :: iEl
   INTEGER, INTENT(in)             :: i, j, varID 
   REAL(prec), INTENT(out)         :: theTend

      CALL myDGSEM % sol(iEl) % GetTendencyAtNodeWithVarID( i, j, varID, theTend )

 END SUBROUTINE GetTendencyAtNodeWithVarID_ShallowWater
!
!
!
 SUBROUTINE SetTendencyAtNodeWithVarID_ShallowWater( myDGSEM, iEl, i, j, varID, theTend  )
 ! S/R SetTendency
 !  
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS(ShallowWater), INTENT(inout) :: myDGSEM
   INTEGER, INTENT(in)                :: iEl
   INTEGER, INTENT(in)                :: i, j, varID 
   REAL(prec), INTENT(in)             :: theTend

      CALL myDGSEM % sol(iEl) % SetTendencyAtNodeWithVarID( i, j, varID, theTend )

 END SUBROUTINE SetTendencyAtNodeWithVarID_ShallowWater
!
! ------------------------------- Boundary Solution ---------------------------------------------- !
!
 SUBROUTINE GetBoundarySolutionAtBoundary_ShallowWater( myDGSEM, iEl, boundary, theSolution  )
 ! S/R GetBoundarySolution
 !  
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS(ShallowWater), INTENT(in) :: myDGSEM
   INTEGER, INTENT(in)             :: iEl
   INTEGER, INTENT(in)             :: boundary
   REAL(prec), INTENT(out)         :: theSolution(0:myDGSEM % nS,1:myDGSEM % nEq)

      CALL myDGSEM % sol(iEl) % GetBoundarySolutionAtBoundary( boundary, theSolution )

 END SUBROUTINE GetBoundarySolutionAtBoundary_ShallowWater
!
!
!
 SUBROUTINE SetBoundarySolutionAtBoundary_ShallowWater( myDGSEM, iEl, boundary, theSolution  )
 ! S/R SetBoundarySolution
 !  
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS(ShallowWater), INTENT(inout) :: myDGSEM
   INTEGER, INTENT(in)                :: iEl
   INTEGER, INTENT(in)                :: boundary
   REAL(prec), INTENT(in)             :: theSolution(0:myDGSEM % nS,1:myDGSEM % nEq)

      CALL myDGSEM % sol(iEl) % SetBoundarySolutionAtBoundary( boundary, theSolution )

 END SUBROUTINE SetBoundarySolutionAtBoundary_ShallowWater
!
! --------------------------------- Boundary Flux ------------------------------------------------ !
!
 SUBROUTINE GetBoundaryFluxAtBoundary_ShallowWater( myDGSEM, iEl, boundary, theFlux  )
 ! S/R GetBoundaryFlux
 !  
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS(ShallowWater), INTENT(in) :: myDGSEM
   INTEGER, INTENT(in)             :: iEl
   INTEGER, INTENT(in)             :: boundary
   REAL(prec), INTENT(out)         :: theFlux(0:myDGSEM % nS,1:myDGSEM % nEq)

      CALL myDGSEM % sol(iEl) % GetBoundaryFluxAtBoundary( boundary, theFlux )

 END SUBROUTINE GetBoundaryFluxAtBoundary_ShallowWater
!
!
!
 SUBROUTINE SetBoundaryFluxAtBoundary_ShallowWater( myDGSEM, iEl, boundary, theFlux  )
 ! S/R SetBoundaryFlux
 !  
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS(ShallowWater), INTENT(inout) :: myDGSEM
   INTEGER, INTENT(in)                :: iEl
   INTEGER, INTENT(in)                :: boundary
   REAL(prec), INTENT(in)             :: theFlux(0:myDGSEM % nS,1:myDGSEM % nEq)

      CALL myDGSEM % sol(iEl) % SetBoundaryFluxAtBoundary( boundary, theFlux )

 END SUBROUTINE SetBoundaryFluxAtBoundary_ShallowWater
!
!
!
 SUBROUTINE GetBoundaryFluxAtBoundaryNode_ShallowWater( myDGSEM, iEl, boundary, i, theFlux  )
 ! S/R GetBoundaryFlux
 !  
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS(ShallowWater), INTENT(in) :: myDGSEM
   INTEGER, INTENT(in)             :: iEl
   INTEGER, INTENT(in)             :: boundary, i
   REAL(prec), INTENT(out)         :: theFlux(1:myDGSEM % nEq)

      CALL myDGSEM % sol(iEl) % GetBoundaryFluxAtBoundaryNode( boundary, i, theFlux )

 END SUBROUTINE GetBoundaryFluxAtBoundaryNode_ShallowWater
!
!
!
 SUBROUTINE SetBoundaryFluxAtBoundaryNode_ShallowWater( myDGSEM, iEl, boundary, i, theFlux  )
 ! S/R SetBoundaryFlux
 !  
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS(ShallowWater), INTENT(inout) :: myDGSEM
   INTEGER, INTENT(in)                :: iEl
   INTEGER, INTENT(in)                :: boundary, i
   REAL(prec), INTENT(in)             :: theFlux(1:myDGSEM % nEq)

      CALL myDGSEM % sol(iEl) % SetBoundaryFluxAtBoundaryNode( boundary, i, theFlux )

 END SUBROUTINE SetBoundaryFluxAtBoundaryNode_ShallowWater
!
! --------------------------------- Relaxation Field --------------------------------------------- !
!
 SUBROUTINE GetRelaxationField_ShallowWater( myDGSEM, iEl, relaxField  )
 ! S/R GetRelaxationField
 !  
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS(ShallowWater), INTENT(in) :: myDGSEM
   INTEGER, INTENT(in)             :: iEl
   REAL(prec), INTENT(out)         :: relaxField(0:myDGSEM % nS, 0:myDGSEM % nP, 1:myDGSEM % nEq)
 
      CALL myDGSEM % relax(iEl) % GetSolutionWithVarID( 1, relaxField(:,:,1) )
      CALL myDGSEM % relax(iEl) % GetSolutionWithVarID( 2, relaxField(:,:,2) )
      CALL myDGSEM % relax(iEl) % GetSolutionWithVarID( 3, relaxField(:,:,3) )
      
 END SUBROUTINE GetRelaxationField_ShallowWater
!
!
!
 SUBROUTINE SetRelaxationField_ShallowWater( myDGSEM, iEl, relaxField  )
 ! S/R SetRelaxationField
 !  
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS(ShallowWater), INTENT(inout) :: myDGSEM
   INTEGER, INTENT(in)                :: iEl
   REAL(prec), INTENT(in)             :: relaxField(0:myDGSEM % nS, 0:myDGSEM % nP, 1:myDGSEM % nEq)

   
      CALL myDGSEM % relax(iEl) % SetSolutionWithVarID( 1, relaxField(:,:,1) )
      CALL myDGSEM % relax(iEl) % SetSolutionWithVarID( 2, relaxField(:,:,2) )
      CALL myDGSEM % relax(iEl) % SetSolutionWithVarID( 3, relaxField(:,:,3) )

 END SUBROUTINE SetRelaxationField_ShallowWater
!
!
!
 SUBROUTINE GetRelaxationFieldAtNode_ShallowWater( myDGSEM, iEl, i, j, relaxField  )
 ! S/R GetRelaxationField
 !  
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS(ShallowWater), INTENT(in) :: myDGSEM
   INTEGER, INTENT(in)             :: iEl, i, j
   REAL(prec), INTENT(out)         :: relaxField(1:myDGSEM % nEq)
 
      CALL myDGSEM % relax(iEl) % GetSolutionAtNodeWithVarID( i, j, 1, relaxField(1) )
      CALL myDGSEM % relax(iEl) % GetSolutionAtNodeWithVarID( i, j, 2, relaxField(2) )
      CALL myDGSEM % relax(iEl) % GetSolutionAtNodeWithVarID( i, j, 3, relaxField(3) )
      
 END SUBROUTINE GetRelaxationFieldAtNode_ShallowWater
!
!
!
 SUBROUTINE SetRelaxationFieldAtNode_ShallowWater( myDGSEM, iEl, i, j, relaxField  )
 ! S/R SetRelaxationFieldAtNode
 !  
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS(ShallowWater), INTENT(inout) :: myDGSEM
   INTEGER, INTENT(in)                :: iEl, i, j
   REAL(prec), INTENT(in)             :: relaxField(1:myDGSEM % nEq)

   
      CALL myDGSEM % relax(iEl) % SetSolutionAtNodeWithVarID( i, j, 1, relaxField(1) )
      CALL myDGSEM % relax(iEl) % SetSolutionAtNodeWithVarID( i, j, 2, relaxField(2) )
      CALL myDGSEM % relax(iEl) % SetSolutionAtNodeWithVarID( i, j, 3, relaxField(3) )

 END SUBROUTINE SetRelaxationFieldAtNode_ShallowWater
!
!
!
 SUBROUTINE GetRelaxationTimeScale_ShallowWater( myDGSEM, iEl, timeScale  )
 ! S/R GetRelaxationTimeScale
 !  
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS(ShallowWater), INTENT(in) :: myDGSEM
   INTEGER, INTENT(in)             :: iEl
   REAL(prec), INTENT(out)         :: timeScale(0:myDGSEM % nS, 0:myDGSEM % nP)

      CALL myDGSEM % relax(iEl) % GetSolutionWithVarID( 4, timeScale )
          
 END SUBROUTINE GetRelaxationTimeScale_ShallowWater
!
!
!
 SUBROUTINE SetRelaxationTimeScale_ShallowWater( myDGSEM, iEl, timeScale  )
 ! S/R SetRelaxationTimeScale
 !  
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS(ShallowWater), INTENT(inout) :: myDGSEM
   INTEGER, INTENT(in)                :: iEl
   REAL(prec), INTENT(in)             :: timeScale(0:myDGSEM % nS, 0:myDGSEM % nP)
   
      CALL myDGSEM % relax(iEl) % SetSolutionWithVarID( 4, timeScale )

 END SUBROUTINE SetRelaxationTimeScale_ShallowWater
!
!
!
 SUBROUTINE GetRelaxationTimeScaleAtNode_ShallowWater( myDGSEM, iEl, i, j, timescale  )
 ! S/R GetRelaxationTimeScale
 !  
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS(ShallowWater), INTENT(in) :: myDGSEM
   INTEGER, INTENT(in)             :: iEl, i, j
   REAL(prec), INTENT(out)         :: timescale
 
      CALL myDGSEM % relax(iEl) % GetSolutionAtNodeWithVarID( i, j, 4, timescale )
      
 END SUBROUTINE GetRelaxationTimeScaleAtNode_ShallowWater
!
!
!
 SUBROUTINE SetRelaxationTimeScaleAtNode_ShallowWater( myDGSEM, iEl, i, j, timescale  )
 ! S/R SetRelaxationTimeScaleAtNode
 !  
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS(ShallowWater), INTENT(inout) :: myDGSEM
   INTEGER, INTENT(in)                :: iEl, i, j
   REAL(prec), INTENT(in)             :: timescale

      CALL myDGSEM % relax(iEl) % SetSolutionAtNodeWithVarID( i, j, 4, timescale )

 END SUBROUTINE SetRelaxationTimeScaleAtNode_ShallowWater
!
! -------------------------------------- Bathymetry ---------------------------------------------- !
!
 SUBROUTINE GetBathymetry_ShallowWater( myDGSEM, iEl, theBathymetry  )
 ! S/R GetBathymetry
 !  
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS(ShallowWater), INTENT(in) :: myDGSEM
   INTEGER, INTENT(in)             :: iEl
   REAL(prec), INTENT(out)         :: theBathymetry(0:myDGSEM % nS, 0:myDGSEM % nP, 1:3)

      CALL myDGSEM % bathymetry(iEl) % GetSolution( theBathymetry )

 END SUBROUTINE GetBathymetry_ShallowWater
!
!
!
 SUBROUTINE SetBathymetry_ShallowWater( myDGSEM, iEl, theBathymetry  )
 ! S/R SetBathymetry
 !  
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS(ShallowWater), INTENT(inout) :: myDGSEM
   INTEGER, INTENT(in)                :: iEl
   REAL(prec), INTENT(in)             :: theBathymetry(0:myDGSEM % nS, 0:myDGSEM % nP, 1:3)

      CALL myDGSEM % bathymetry(iEl) % SetSolution( theBathymetry )

 END SUBROUTINE SetBathymetry_ShallowWater
!
!
!
 SUBROUTINE GetBathymetryAtNode_ShallowWater( myDGSEM, iEl, i, j, theBathymetry  )
 ! S/R GetBathymetry
 !  
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS(ShallowWater), INTENT(in) :: myDGSEM
   INTEGER, INTENT(in)             :: iEl
   INTEGER, INTENT(in)             :: i, j
   REAL(prec), INTENT(out)         :: theBathymetry(1:3)

      CALL myDGSEM % bathymetry(iEl) % GetSolutionAtNode( i, j, theBathymetry )

 END SUBROUTINE GetBathymetryAtNode_ShallowWater
!
!
!
 SUBROUTINE SetBathymetryAtNode_ShallowWater( myDGSEM, iEl, i, j, theBathymetry  )
 ! S/R SetBathymetry
 !  
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS(ShallowWater), INTENT(inout) :: myDGSEM
   INTEGER, INTENT(in)                :: iEl
   INTEGER, INTENT(in)                :: i, j
   REAL(prec), INTENT(in)             :: theBathymetry(1:3)

      CALL myDGSEM % bathymetry(iEl) % SetSolutionAtNode( i, j, theBathymetry )

 END SUBROUTINE SetBathymetryAtNode_ShallowWater
!
!
!
 SUBROUTINE GetBathymetryAtBoundary_ShallowWater( myDGSEM, iEl, boundary, theBathymetry  )
 ! S/R GetBathymetryAtBoundary
 !  
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS(ShallowWater), INTENT(in) :: myDGSEM
   INTEGER, INTENT(in)             :: iEl
   INTEGER, INTENT(in)             :: boundary
   REAL(prec), INTENT(out)         :: theBathymetry(0:myDGSEM % nS,1:3)

      CALL myDGSEM % bathymetry(iEl) % GetBoundarySolutionAtBoundary( boundary, theBathymetry )

 END SUBROUTINE GetBathymetryAtBoundary_ShallowWater
!
!
!
 SUBROUTINE SetBathymetryAtBoundary_ShallowWater( myDGSEM, iEl, boundary, theBathymetry  )
 ! S/R SetBathymetryAtBoundary
 !  
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS(ShallowWater), INTENT(inout) :: myDGSEM
   INTEGER, INTENT(in)                :: iEl
   INTEGER, INTENT(in)                :: boundary
   REAL(prec), INTENT(in)             :: theBathymetry(0:myDGSEM % nS,1:3)

      CALL myDGSEM % bathymetry(iEl) % SetBoundarySolutionAtBoundary( boundary, theBathymetry )

 END SUBROUTINE SetBathymetryAtBoundary_ShallowWater
!
! -------------------------------------- Vorticity ---------------------------------------------- !
!
 SUBROUTINE GetPlanetaryVorticity_ShallowWater( myDGSEM, iEl, theVorticity  )
 ! S/R GetVorticity
 !  
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS(ShallowWater), INTENT(in) :: myDGSEM
   INTEGER, INTENT(in)             :: iEl
   REAL(prec), INTENT(out)         :: theVorticity(0:myDGSEM % nS, 0:myDGSEM % nP)

      CALL myDGSEM % vorticity(iEl) % GetSolutionWithVarID( 1, theVorticity )

 END SUBROUTINE GetPlanetaryVorticity_ShallowWater
!
!
!
 SUBROUTINE SetPlanetaryVorticity_ShallowWater( myDGSEM, iEl, theVorticity  )
 ! S/R SetVorticity
 !  
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS(ShallowWater), INTENT(inout) :: myDGSEM
   INTEGER, INTENT(in)                :: iEl
   REAL(prec), INTENT(in)             :: theVorticity(0:myDGSEM % nS, 0:myDGSEM % nP)

      CALL myDGSEM % vorticity(iEl) % SetSolutionWithVarID( 1, theVorticity )

 END SUBROUTINE SetPlanetaryVorticity_ShallowWater
!
!
!
 SUBROUTINE GetPlanetaryVorticityAtNode_ShallowWater( myDGSEM, iEl, i, j, theVorticity  )
 ! S/R GetVorticity
 !  
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS(ShallowWater), INTENT(in) :: myDGSEM
   INTEGER, INTENT(in)             :: iEl
   INTEGER, INTENT(in)             :: i, j
   REAL(prec), INTENT(out)         :: theVorticity

      CALL myDGSEM % vorticity(iEl) % GetSolutionAtNodeWithVarID( i, j, 1, theVorticity )

 END SUBROUTINE GetPlanetaryVorticityAtNode_ShallowWater
!
!
!
 SUBROUTINE SetPlanetaryVorticityAtNode_ShallowWater( myDGSEM, iEl, i, j, theVorticity  )
 ! S/R SetVorticity
 !  
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS(ShallowWater), INTENT(inout) :: myDGSEM
   INTEGER, INTENT(in)                :: iEl
   INTEGER, INTENT(in)                :: i, j
   REAL(prec), INTENT(in)             :: theVorticity

      CALL myDGSEM % vorticity(iEl) % SetSolutionAtNodeWithVarID( i, j, 1, theVorticity )

 END SUBROUTINE SetPlanetaryVorticityAtNode_ShallowWater
!
!
!
 SUBROUTINE GetPlanetaryVorticityAtBoundary_ShallowWater( myDGSEM, iEl, boundary, theVorticity  )
 ! S/R GetVorticityAtBoundary
 !  
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS(ShallowWater), INTENT(in) :: myDGSEM
   INTEGER, INTENT(in)             :: iEl
   INTEGER, INTENT(in)             :: boundary
   REAL(prec), INTENT(out)         :: theVorticity(0:myDGSEM % nS)

      CALL myDGSEM % vorticity(iEl) % GetBoundarySolutionAtBoundaryWithVarID( boundary, 1, theVorticity )

 END SUBROUTINE GetPlanetaryVorticityAtBoundary_ShallowWater
!
!
!
 SUBROUTINE SetPlanetaryVorticityAtBoundary_ShallowWater( myDGSEM, iEl, boundary, theVorticity  )
 ! S/R SetVorticityAtBoundary
 !  
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS(ShallowWater), INTENT(inout) :: myDGSEM
   INTEGER, INTENT(in)                :: iEl
   INTEGER, INTENT(in)                :: boundary
   REAL(prec), INTENT(in)             :: theVorticity(0:myDGSEM % nS)

      CALL myDGSEM % vorticity(iEl) % SetBoundarySolutionAtBoundaryWithVarID( boundary, 1, theVorticity )

 END SUBROUTINE SetPlanetaryVorticityAtBoundary_ShallowWater
!
!
!==================================================================================================!
!--------------------------------- Type Specific Routines -----------------------------------------!
!==================================================================================================!
!
!
 SUBROUTINE CalculateSolutionAtBoundaries_ShallowWater( myDGSEM, iEl )
 ! S/R CalculateSolutionAtBoundaries
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS(ShallowWater), INTENT(inout) :: myDGSEM
   INTEGER, INTENT(in)                :: iEl
   
      CALL myDGSEM % sol(iEl) % CalculateSolutionAtBoundaries( myDGSEM % dgStorage )
      
 END SUBROUTINE CalculateSolutionAtBoundaries_ShallowWater
!
!
!
 SUBROUTINE CalculateVorticityAtBoundaries_ShallowWater( myDGSEM, iEl )
 ! S/R CalculateVorticityAtBoundaries
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS(ShallowWater), INTENT(inout) :: myDGSEM
   INTEGER, INTENT(in)                :: iEl
   
      CALL myDGSEM % vorticity(iEl) % CalculateSolutionAtBoundaries( myDGSEM % dgStorage )
      
 END SUBROUTINE CalculateVorticityAtBoundaries_ShallowWater
!
!
!
 SUBROUTINE CalculateBathymetryAtBoundaries_ShallowWater( myDGSEM, iEl )
 ! S/R CalculateBathymetryAtBoundaries
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS(ShallowWater), INTENT(inout) :: myDGSEM
   INTEGER, INTENT(in)                :: iEl
   
      CALL myDGSEM % bathymetry(iEl) % CalculateSolutionAtBoundaries( myDGSEM % dgStorage )
      
 END SUBROUTINE CalculateBathymetryAtBoundaries_ShallowWater
!
!
!
 SUBROUTINE ForwardStepRK3_ShallowWater( myDGSEM, tn )
 ! S/R ForwardStepRK3( 3rd order Runge-Kutta)
 ! 
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS(ShallowWater), INTENT(inout) :: myDGSEM
   REAL(prec), INTENT(in)                   :: tn
   ! LOCAL
   REAL(prec) :: t, dt
   REAL(prec) :: G2D(0:myDGSEM % nS,&
                     0:myDGSEM % nP,&
                     1:myDGSEM % nEq,&
                     1:myDGSEM % mesh % nElems) 
   REAL(prec) :: dSdt(0:myDGSEM % nS,&
                     0:myDGSEM % nP,&
                     1:myDGSEM % nEq )
   REAL(prec) :: dVDt(0:myDGSEM % nS,&
                      0:myDGSEM % nP,&
                      1:2 )
   INTEGER    :: m, iP, jP, nX, nY, nZ, nEq, iEqn, iEl

     nEq = myDGSEM % nEq 
     dt = myDGSEM % params % dt
     G2D = ZERO
    
     DO m = 1,3 ! Loop over RK3 steps

        t = tn + rk3_b(m)*dt
        ! Calculate the tendency
        CALL myDGSEM % GlobalTimeDerivative( t )
        
        DO iEl = 1, myDGSEM % mesh % nElems ! Loop over all of the elements

           CALL myDGSEM % GetTendency( iEl, dSdt )
           G2D(:,:,:,iEl) = rk3_a(m)*G2D(:,:,:, iEl) + dSdt

           myDGSEM % sol(iEl) % solution = myDGSEM % sol(iEl) % solution + rk3_g(m)*dt*G2D(:,:,:,iEl)
           

         ENDDO ! iEl, loop over all of the elements
         
         
      ENDDO ! m, loop over the RK3 steps
   
   
       
 END SUBROUTINE ForwardStepRK3_ShallowWater
!
!
!
 SUBROUTINE GlobalTimeDerivative_ShallowWater( myDGSEM, tn ) 
 ! S/R GlobalTimeDerivative_ShallowWater
 ! 
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS(ShallowWater), INTENT(inout) :: myDGSEM
   REAL(prec), INTENT(in)             :: tn
   ! Local
   INTEGER :: iEl, iEdge
   REAL(prec) :: cpuStart, cpuEnd


!$OMP PARALLEL
      ! Calculate the solution at the boundaries
!$OMP DO
      DO iEl = 1, myDGSEM % mesh % nElems
         CALL myDGSEM % CalculateSolutionAtBoundaries( iEl ) 
      ENDDO 
!$OMP END DO 
!$OMP FLUSH( myDGSEM )
      ! Calculate the flux along the element edges
!$OMP DO
      DO iEdge = 1, myDGSEM % mesh % nEdges
         CALL myDGSEM % EdgeFlux( iEdge, tn )
      ENDDO 
!$OMP END DO
!$OMP FLUSH( myDGSEM ) 


!$OMP DO
      DO iEl = 1, myDGSEM % mesh % nElems
         CALL myDGSEM % MappedTimeDerivative( iEl, tn )
      ENDDO
!$OMP END DO
!$OMP FLUSH( myDGSEM )
!$OMP END PARALLEL

 
 END SUBROUTINE GlobalTimeDerivative_ShallowWater
!
!
! 
 SUBROUTINE EdgeFlux_ShallowWater( myDGSEM, iEdge, tn )
 ! S/R EdgeFlux
 ! 
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS( ShallowWater ), INTENT(inout) :: myDGSEM
   INTEGER, INTENT(in)                  :: iEdge  
   REAL(prec), INTENT(in)               :: tn
   ! Local
   INTEGER :: k, iNode, iEqn
   INTEGER :: e1, s1, e2, s2, iError, start, inc, nS, nP
   REAL(prec) :: flux(1:nSWeq)
   REAL(prec) :: inState(0:myDGSEM % nS,1:nSWeq)
   REAL(prec) :: exState(0:myDGSEM % nS,1:nSWeq)
   REAL(prec) :: avgState(0:myDGSEM % nS,1:nSWeq)
   REAL(prec) :: avgState2(0:myDGSEM % nS,1:nSWeq)
   REAL(prec) :: inBathy(0:myDGSEM % nS,1:3)
   REAL(prec) :: exBathy(0:myDGSEM % nS,1:3)
   REAL(prec) :: wavespeed, h, x, y 
   REAL(preC) :: nHat(1:2), nHatLength 
    
      CALL myDGSEM % dgStorage % GetNumberOfNodes( nS, nP )

      CALL myDGSEM % mesh % GetEdgePrimaryElementID( iEdge, e1 )
      CALL myDGSEM % mesh % GetEdgePrimaryElementSide( iEdge, s1 )
      CALL myDGSEM % mesh % GetEdgeSecondaryElementID( iEdge, e2 )
      CALL myDGSEM % mesh % GetEdgeSecondaryElementSide( iEdge, s2 )
      
      s2 = ABS(s2)

      CALL myDGSEM % mesh % GetEdgeStart( iEdge, start )
      CALL myDGSEM % mesh % GetEdgeIncrement( iEdge, inc )
      
      CALL myDGSEM % GetBoundarySolutionAtBoundary( e1, s1, inState )
      CALL myDGSEM % GetBathymetryAtBoundary( e1, s1, inBathy )
      
      IF( e2 > 0 )then ! this is an interior edge
      
         k = start-inc
         CALL myDGSEM % GetBoundarySolutionAtBoundary( e2, s2, exState )
         CALL myDGSEM % GetBathymetryAtBoundary( e2, s2, exBathy )
         
         
         DO iNode = 0, nS ! Loop over the nodes
          
            h = HALF*(inBathy(iNode, 1) + exBathy(k, 1)) ! calculate the depth
            
            CALL myDGSEM % mesh % GetBoundaryNormalAtNode( e1, nHat, nHatLength, iNode, s1 ) ! Get nHat
           
            ! Calculate the RIEMANN flux
            flux = RiemannSolver( inState(iNode,:), exState(k,:), h, nHat, myDGSEM % params )*nHatLength

            ! Store the flux for the elements which share this edge
            CALL myDGSEM % SetBoundaryFluxAtBoundaryNode( e1, s1, iNode, flux )
            CALL myDGSEM % SetBoundaryFluxAtBoundaryNode( e2, s2, k, -flux )
             
            !avgState(iNode,:) = HALF*( inState(iNode,:) + exState(k,:) )
            !avgState2(k,:)    = HALF*( inState(iNode,:) + exState(k,:) )
            k = k + inc

         ENDDO ! iNode, loop over the nodes

         ! Set the boundary solutions to the average of the two solutions at the edge - this will 
         ! only affect the relative vorticity calculation. If you are using V2 of CalculateRelativeVorticity,
         ! the edge velocity is used in a weak form calculation of the relative vorticity.
        ! CALL myDGSEM % SetBoundarySolutionAtBoundary( e1, s1, avgState )
        ! CALL myDGSEM % SetBoundarySolutionAtBoundary( e2, s2, avgState2 )
         
      ELSE ! this edge is a boundary edge

         DO iNode = 0, nS ! loop over the nodes on this edge
       
            CALL myDGSEM % mesh % GetBoundaryNormalAtNode( e1, nHat, nHatLength, iNode, s1 ) ! Get nHat

            h = inBathy(iNode, 1) ! get the depth on the primary element boundary (s1)
            
            ! Get the boundary point locations
            CALL myDGSEM % mesh % GetBoundaryLocationAtNode( e1, x, y, iNode, s1 )

            ! Calculate the external state
            exState(iNode,1:nSWEq) = GetExternalState( nSWEq, nHat, x, y, h, tn, e2, &
                                                      inState(iNode,1:nSWeq), myDGSEM % params )

            ! Calculate the RIEMANN flux
            flux = RiemannSolver( inState(iNode,:), exState(iNode,:), h, nHat, myDGSEM % params )*nHatLength

            CALL myDGSEM % SetBoundaryFluxAtBoundaryNode( e1, s1, iNode, flux )

         ENDDO ! iNode, loop over the nodes on this edge

      ENDIF ! Choosing the edge type
 
 END SUBROUTINE EdgeFlux_ShallowWater
!
!
!
 SUBROUTINE MappedTimeDerivative_ShallowWater( myDGSEM, iEl, tn ) 
 ! S/R MappedTimeDerivative_ShallowWater
 ! 
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS(ShallowWater), INTENT(inout) :: myDGSEM
   INTEGER, INTENT(in)                :: iEl
   REAL(prec), INTENT(in)             :: tn
   ! LOCAL
   REAL(prec) :: dxds, dxdp, dyds, dydp
   REAL(prec) :: h, dhdx, dhdy, J, bathy(1:3), plv, rev, vorticity, x, y
   REAL(prec) :: xF(1:nSWEq), yF(1:nSWEq), sol(1:nSWEq), relaxSol(1:nSWeq), timescale
   REAL(prec) :: pContFlux(0:myDGSEM % nS,1:nSWEq)
   REAL(prec) :: sContFlux(0:myDGSEM % nS,1:nSWEq)
   REAL(prec) :: pContFluxDer(0:myDGSEM % nS,1:nSWEq)
   REAL(prec) :: sContFluxDer(0:myDGSEM % nS,1:nSWEq)
   REAL(prec) :: tend(0:myDGSEM % nS,0:myDGSEM % nP,1:nSWEq)
   REAL(prec) :: dMatX(0:myDGSEM % nS,0:myDGSEM % nS)
   REAL(prec) :: dMatY(0:myDGSEM % nP,0:myDGSEM % nP)
   REAL(prec) :: qWeightX(0:myDGSEM % nS)
   REAL(prec) :: qWeightY(0:myDGSEM % nP)
   REAL(prec) :: lsouth(0:myDGSEM % nP)
   REAL(prec) :: lnorth(0:myDGSEM % nP)
   REAL(prec) :: least(0:myDGSEM % nS)
   REAL(prec) :: lwest(0:myDGSEM % nS)
   REAL(prec) :: fL(1:nSWEq), fR(1:nSWEq)
   INTEGER    :: iS, iP, nS, nP, iEq

     CALL myDGSEM % dgStorage % GetNumberOfNodes( nS, nP )
     CALL myDGSEM % dgStorage % GetDerivativeMatrix( dMatX, dMatY )
     CALL myDGSEM % dgStorage % GetQuadratureWeights( qWeightX, qWeightY )
     CALL myDGSEM % dgStorage % GetSouthernInterpolants( lsouth )
     CALL myDGSEM % dgStorage % GetNorthernInterpolants( lnorth )
     CALL myDGSEM % dgStorage % GetEasternInterpolants( least )
     CALL myDGSEM % dgStorage % GetWesternInterpolants( lwest )
          
     DO iP = 0,nP ! Loop over the y-points

        DO iS = 0,nS ! Loop over the x-points

           ! Get the metric terms to calculate the contravariant flux
           CALL myDGSEM % mesh % GetCovariantMetricsAtNode( iEl, dxds, dxdp, dyds, dydp, iS, iP )

           ! Get the solution at x,y
           CALL myDGSEM % GetSolutionAtNode( iEl, iS, iP, sol )

           ! Get the depth - needed for wave speed calculation 
           CALL myDGSEM % GetBathymetryAtNode( iEl, iS, iP, bathy )
           h = bathy(1)

           ! Calculate the x and y fluxes
           xF = XFlux(  tn, h, myDGSEM % params, sol )
           yF = YFlux(  tn, h, myDGSEM % params, sol )

           !And now the contravariant flux
           sContFlux(iS,1:nSWEq) = dydp*xF(1:nSWEq) - dxdp*yF(1:nSWEq)

        ENDDO ! iS, Loop over the x-points

        ! Get the numerical fluxes at the boundaries (should be calculated before
        ! this routine is called )
        
        ! Get the "west" flux at iP
        CALL myDGSEM % GetBoundaryFluxAtBoundaryNode( iEl, 4, iP, fL )

        ! Get the "east" flux at iP
        CALL myDGSEM % GetBoundaryFluxAtBoundaryNode( iEl, 2, iP, fR )

        ! At this y-level, calculate the DG-advective derivative
        DO iEq = 1, nSWeq
           sContFluxDer(0:nS,iEq) = DGSystemDerivative(  nS, dMatX, qWeightX, fL(iEq), fR(iEq), &
                                                         sContFlux(:,iEq), lwest, least  )
        ENDDO
     
        DO iS = 0, nS ! Loop over the x-points
           DO iEq = 1, nSWEq  ! Loop over the number of equations
              tend(iS, iP, iEq) = -sContFluxDer(iS,iEq)
           ENDDO ! Loop over the number of equations
        ENDDO ! Loop over the x-points

      
     ENDDO ! iP  


     DO iS = 0,nS ! Loop over the x-points

        DO iP = 0,nP ! Loop over the y-points

           ! Get the metric terms to calculate the contravariant flux
           CALL myDGSEM % mesh % GetCovariantMetricsAtNode( iEl, dxds, dxdp, dyds, dydp, iS, iP )

           ! Get the solution at x,y
           CALL myDGSEM % GetSolutionAtNode( iEl, iS, iP, sol )

           ! Get the depth - needed for wave speed calculation 
           CALL myDGSEM % GetBathymetryAtNode( iEl, iS, iP, bathy )
           h = bathy(1)

           ! Calculate the x and y fluxes
           xF = XFlux(  tn, h, myDGSEM % params, sol )
           yF = YFlux(  tn, h, myDGSEM % params, sol )

           !And now the contravariant flux
           pContFlux(iP,1:nSWEq) = -dyds*xF(1:nSWEq) + dxds*yF(1:nSWEq)

        ENDDO ! iP, Loop over the y-points

        ! Get the numerical fluxes at the boundaries (should be calculated before
        ! this routine is called )
        
        ! Get the "south" flux at iS
        CALL  myDGSEM % GetBoundaryFluxAtBoundaryNode( iEl, 1, iS, fL )

        ! Get the "north" flux at iS
        CALL  myDGSEM % GetBoundaryFluxAtBoundaryNode( iEl, 3, iS, fR )


        ! At this x-level, calculate the y-DG-advective derivative
        DO iEq = 1, nSWeq
            pContFluxDer(0:nP,iEq) = DGSystemDerivative(  nP, dMatY, qWeightY, fL(iEq), fR(iEq), &
                                                          pContFlux(:,iEq), lsouth, lnorth  )
        ENDDO

         
        DO iP = 0, nP ! Loop over the y-points
           CALL myDGSEM % mesh % GetJacobianAtNode( iEl, J, iS, iP )
           DO iEq = 1, nSWEq  ! Loop over the number of equations
               tend(iS, iP, iEq) = ( tend(iS, iP, iEq)- pContFluxDer(iP,iEq) )/J
           ENDDO ! Loop over the number of equations
        ENDDO ! Loop over the x-points

      
     ENDDO ! iP
         
     ! COMPUTE THE SOURCE TERMS
     DO iS = 0, nS  ! Loop over the x-points
        DO iP = 0, nP ! Loop over the y-points

           CALL myDGSEM % GetSolutionAtNode( iEl, iS, iP, sol )
           CALL myDGSEM % GetRelaxationFieldAtNode( iEl, iS, iP, relaxSol )
           CALL myDGSEM % GetRelaxationTimeScaleAtNode( iEl, iS, iP, timescale )
           CALL myDGSEM % GetPlanetaryVorticityAtNode( iEl, iS, iP, plv )
           CALL myDGSEM % GetBathymetryAtNode( iEl, iS, iP, bathy )
           CALL myDGSEM % mesh % GetPositionAtNode( iEl, x, y, iS, iP )
           
           vorticity = plv
                  
           tend(iS,iP,1:nSWEq) = tend(iS,iP,1:nSWEq)+ Source( tn, sol, relaxSol, timescale, &
                                                              plv, bathy, &
                                                              x, y, myDGSEM % params  )

        ENDDO ! iP, Loop over the y-points
     ENDDO ! iS, Loop over the x-points

     CALL myDGSEM % SetTendency( iEl, tend )

 END SUBROUTINE MappedTimeDerivative_ShallowWater
!
!
!  
FUNCTION DGSystemDerivative(  nP, dMat, qWei, lFlux, rFlux, intFlux, lagLeft, lagRight  ) RESULT( tendency )
 ! FUNCTION DGSystemDerivative ( Discontinous Galerkin time Derivative)
 ! 
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   INTEGER    :: nP
   REAL(prec) :: dMat(0:nP,0:nP)
   REAL(prec) :: qWei(0:nP)
   REAL(prec) :: lFlux
   REAL(prec) :: rFlux
   REAL(prec) :: intFlux(0:nP)
   REAL(prec) :: lagLeft(0:nP), lagRight(0:nP)
   REAL(prec) :: tendency(0:nP)
   ! Local 
   INTEGER :: kP

      ! Apply the derivative operator
      tendency = MATMUL( dMat, intflux )

      DO kP = 0, nP ! Loop over the quadrature points
            tendency(kP) = tendency(kP) + ( rFlux*lagRight(kP) + lFlux*lagLeft(kP) )/qWei(kP)
      ENDDO ! kP, loop over the quadrature points
                                         
 END FUNCTION DGSystemDerivative
!
!
!
 FUNCTION RiemannSolver( inState, outState, h, nHat, locParams ) RESULT( numFlux )
 ! FUNCTION RiemannSolver 
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   REAL(prec) :: inState(1:nSWEq)
   REAL(prec) :: outState(1:nSWEq)
   REAL(prec) :: h
   REAL(prec) :: nHat(1:2) ! normal direction
   REAL(prec) :: numFlux(1:nSWEq)
   TYPE( SWParams ) :: locParams
   ! LOCAL
   REAL(prec) :: uL, uR, vL, vR, pL, pR, uOut, uIn
   REAL(prec) :: gamm1, gamm2, locG
   REAL(prec) :: jump(1:nSWEq), aS(1:nSWEq)
   REAL(prec) :: uNorm, fac


      jump = outState - inState
    
      locG = locParams % g
      uL = inState(1)
      vL = inState(2)
      pL = inState(3)

      uR = outState(1)
      vR = outState(2)
      pR = outState(3)

      uOut =  uR*nHat(1) + vR*nHat(2)

      uIn = uL*nHat(1) + vL*nHat(2)

      IF( locParams % ModelFormulation == LINEAR ) then
         
         fac = sqrt(locG*h)

         numFlux(1) = HALF*( locG*(pR + pL) + fac*(uIn-uOut) )*nHat(1)
         numFlux(2) = HALF*( locG*(pR + pL) + fac*(uIn-uOut) )*nHat(2)
         numFlux(3) = HALF*( h*(uIn + uOut) + fac*(pL - pR) )

      ELSE ! CONSERVATIVE

         fac = max( abs(HALF*( uOut + sqrt( uOut*uOut + 2.0_prec*locg*(pR*pR*pR) ) )/pR ), &
                    abs(HALF*( uOut - sqrt( uOut*uOut + 2.0_prec*locg*(pR*pR*pR) ) )/pR ), &
                    abs(HALF*( uIn + sqrt( uIn*uIn + 2.0_prec*locg*(pL*pL*pL) ) )/pL ), &
                    abs(HALF*( uIn - sqrt( uIn*uIn + 2.0_prec*locg*(pL*pL*pL) ) )/pL ) )

         aS(1) = (uOut/pR)*uR + HALF*locg*pR*pR*nHat(1) +&
                 (uIn/pL)*uL + HALF*locg*pL*pL*nHat(1)

         aS(2) = (uOut/pR)*vR + HALF*locg*pR*pR*nHat(2) +&
                 (uIn/pL)*vL + HALF*locg*pL*pL*nHat(2)

         aS(3) = uOut + uIn

         numFlux = HALF*( aS - fac*jump ) 

      ENDIF     

 END FUNCTION RiemannSolver
!
!
! 
 FUNCTION XFlux(  tn, h, params, solAtX ) RESULT( fx )
 ! FUNCTION XFlux
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   REAL(prec) :: tn, h
   REAL(prec) :: solAtX(1:nSWeq)
   REAL(prec) :: fx(1:nSWeq)
   TYPE(SWParams) :: params
   ! LOCAL
   REAL(prec) :: KE, g 

      g = params % g
      
      IF( params % modelFormulation == LINEAR ) then

         fx(1) = g*solAtX(3)

         fx(2) = ZERO

         fx(3) = h*solAtX(1)
         
      ELSE

         fx(1) = solAtX(1)*solAtX(1)/(solAtX(3)) + g*HALF*solAtX(3)*solAtX(3)

         fx(2) = solAtX(1)*solAtX(2)/(solAtX(3))

         fx(3) = solAtX(1)

      ENDIF                                                 

 END FUNCTION XFlux
!
!
!
 FUNCTION YFlux(  tn, h, params, solAtX ) RESULT( fy )
 ! FUNCTION YFlux
 ! 
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   REAL(prec) :: tn, h
   REAL(prec) :: solAtX(1:nSWeq)
   REAL(prec) :: fy(1:nSWeq)
   TYPE(SWParams) :: params
   ! LOCAL
   REAL(prec) :: KE, g 

      g = params % g

      IF( params % modelFormulation == LINEAR ) then

         fy(1) = ZERO

         fy(2) = g*solAtX(3)

         fy(3) = h*solAtX(2)
         
      ELSE

         fy(1) = solAtX(2)*solAtX(1)/(solAtX(3)) 

         fy(2) = solAtX(2)*solAtX(2)/(solAtX(3)) + g*HALF*solAtX(3)*solAtX(3)

         fy(3) = solAtX(2)

      ENDIF                                                

 END FUNCTION YFlux
!
!
!               
 FUNCTION Source( tn, sol, relax, timescale, vorticity, bathy, x, y, params ) RESULT( q )
 ! FUNCTION Source
 ! 
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   REAL(prec)     :: tn
   REAL(prec)     :: sol(1:nSWEq), relax(1:nSWEq), timescale
   REAL(prec)     :: vorticity, bathy(1:3)
   REAL(prec)     :: x, y
   TYPE(SWParams) :: params
   REAL(prec)     :: q(1:nSWeq)
   ! LOCAL
   REAL(prec) :: lDrag, h, dhdx, dhdy, g
  
     lDrag = params % LinearDrag
     g     = params % g
     h     = bathy(1)
     dhdx  = bathy(2)
     dhdy  = bathy(3)
     
     IF( params % ModelFormulation == LINEAR ) then

           ! u-momentum source = f*v, or (f+zeta)*v 
           q(1) =  vorticity*sol(2) - lDrag*sol(1) - (sol(1) - relax(1))*timescale                 
           ! v-momentum source = -f*u, or -(f+zeta)*v 
           q(2) = -vorticity*sol(1) - lDrag*sol(2) - (sol(2) - relax(2))*timescale              
           ! pressure source
           q(3) = -(sol(3) - relax(3))*timescale

     ELSE

           q(1) = vorticity*sol(2) + g*sol(3)*dhdx - lDrag*sol(1)/h -& 
                  (sol(1) - relax(1))*timescale                    !u-momentum source = f*v
    
           q(2) = -vorticity*sol(1) + g*sol(3)*dhdy - lDrag*sol(2)/h -&
                   (sol(2) - relax(2))*timescale                  !v-momentum source = -f*u
 
           q(3) = -(sol(3) - relax(3))*timescale   ! pressure source   

     ENDIF

 END FUNCTION Source
!
!
!                         
 FUNCTION GetExternalState( nEqn, nHat, x, y, h, t, bcFlag, intState, locParams) RESULT( extState )
 ! S/R GetExternalState
 ! 
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   INTEGER    :: nEqn
   REAL(prec) :: intState(1:nEqn)
   REAL(prec) :: x,y,t,h
   REAL(prec) :: nHat(1:2)
   INTEGER    :: bcFlag, formulation
   REAL(prec) :: extState(1:nEqn)
   TYPE( SWParams ) :: locParams

       IF( bcFlag == NO_NORMAL_FLOW ) then
          
          extState(1) = (nHat(2)*nHat(2) -nHat(1)*nHat(1))*intState(1) - 2.0_prec*nHat(1)*nHat(2)*intState(2) ! u velocity
          extState(2) = (nHat(1)*nHat(1) -nHat(2)*nHat(2))*intState(2) - 2.0_prec*nHat(1)*nHat(2)*intState(1) ! v velocity
          extState(3) = intState(3) ! barotropic pressure
             
       ELSEIF( bcFlag == RADIATION ) then
          extState = ZERO
          
       ELSE       
       
         PRINT*, 'FUNCTION : GetExternalState : Invalid bcflag =', bcflag
         extState = ZERO
         
       ENDIF

        
 END FUNCTION GetExternalState
!
!
!
 SUBROUTINE CalculateBathymetryGradient_ShallowWater( myDGSEM, iEl ) 
 ! S/R CalculateBathymetryGradient
 ! 
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS(ShallowWater), INTENT(inout) :: myDGSEM
   INTEGER, INTENT(in)                :: iEl
   ! LOCAL
   REAL(prec) :: dxds, dxdp, dyds, dydp, J, normLength, normVec(1:2)
   REAL(prec) :: h, bathy(1:3), plv, rev
   REAL(prec) :: xF(1:2), yF(1:2)
   REAL(prec) :: pContFlux(0:myDGSEM % nS,1:2)
   REAL(prec) :: sContFlux(0:myDGSEM % nS,1:2)
   REAL(prec) :: pContFluxDer(0:myDGSEM % nS,1:2)
   REAL(prec) :: sContFluxDer(0:myDGSEM % nS,1:2)
   REAL(prec) :: tend(0:myDGSEM % nS,0:myDGSEM % nP,1:2)
   REAL(prec) :: dMatX(0:myDGSEM % nS,0:myDGSEM % nS)
   REAL(prec) :: dMatY(0:myDGSEM % nP,0:myDGSEM % nP)
   REAL(prec) :: qWeightX(0:myDGSEM % nS)
   REAL(prec) :: qWeightY(0:myDGSEM % nP)
   REAL(prec) :: lsouth(0:myDGSEM % nP)
   REAL(prec) :: lnorth(0:myDGSEM % nP)
   REAL(prec) :: least(0:myDGSEM % nS)
   REAL(prec) :: lwest(0:myDGSEM % nS)
   REAL(prec) :: hsouth(0:myDGSEM % nP,1:3)
   REAL(prec) :: hnorth(0:myDGSEM % nP,1:3)
   REAL(prec) :: heast(0:myDGSEM % nS,1:3)
   REAL(prec) :: hwest(0:myDGSEM % nS,1:3)
   REAL(prec) :: fL(1:2), fR(1:2)
   INTEGER    :: iS, iP, nS, nP, iEq

     CALL myDGSEM % dgStorage % GetNumberOfNodes( nS, nP )
     CALL myDGSEM % dgStorage % GetDerivativeMatrix( dMatX, dMatY )
     CALL myDGSEM % dgStorage % GetQuadratureWeights( qWeightX, qWeightY )
     CALL myDGSEM % dgStorage % GetSouthernInterpolants( lsouth )
     CALL myDGSEM % dgStorage % GetNorthernInterpolants( lnorth )
     CALL myDGSEM % dgStorage % GetEasternInterpolants( least )
     CALL myDGSEM % dgStorage % GetWesternInterpolants( lwest )
         
     CALL myDGSEM % GetBathymetryAtBoundary( iEl, 1, hSouth )
     CALL myDGSEM % GetBathymetryAtBoundary( iEl, 2, hEast )
     CALL myDGSEM % GetBathymetryAtBoundary( iEl, 3, hNorth )
     CALL myDGSEM % GetBathymetryAtBoundary( iEl, 4, hWest ) 
     DO iP = 0,nP ! Loop over the y-points

        DO iS = 0,nS ! Loop over the x-points

           ! Get the metric terms to calculate the contravariant flux
           CALL myDGSEM % mesh % GetCovariantMetricsAtNode( iEl, dxds, dxdp, dyds, dydp, iS, iP )
           ! Get the depth - needed for wave speed calculation 
           CALL myDGSEM % GetBathymetryAtNode( iEl, iS, iP, bathy )
           h = bathy(1)

           ! Calculate the x and y fluxes
           xF = (/ h, ZERO /)
           yF = (/ ZERO, h /)

           !And now the contravariant flux
           sContFlux(iS,1:2) = dydp*xF - dxdp*yF

        ENDDO ! iS, Loop over the x-points

        ! East edge
         CALL myDGSEM % mesh % GetBoundaryNormalAtNode( iEl, normVec, normlength, iP, 2 ) 
         fR = ( hEast(iP,1)*normVec )*normLength
           
         ! West edge
         CALL myDGSEM % mesh % GetBoundaryNormalAtNode( iEl, normVec, normlength, iP, 4 ) 
         fL = ( hWest(iP,1)*normVec )*normLength
         
        ! At this y-level, calculate the DG-advective derivative
         DO iEq = 1, 2
            sContFluxDer(0:nS,iEq) = DGSystemDerivative(  nS, dMatX, qWeightX, fL(iEq), fR(iEq), &
                                                          sContFlux(:,iEq), lwest, least )
         ENDDO

         DO iS = 0, nS ! Loop over the x-points
            DO iEq = 1, 2 ! Loop over the number of equations
               tend(iS, iP, iEq) = sContFluxDer(iS,iEq)
            ENDDO ! Loop over the number of equations
         ENDDO ! Loop over the x-points

      
     ENDDO ! iP  


     DO iS = 0,nS ! Loop over the x-points

        DO iP = 0,nP ! Loop over the y-points

           ! Get the metric terms to calculate the contravariant flux
           CALL myDGSEM % mesh % GetCovariantMetricsAtNode( iEl, dxds, dxdp, dyds, dydp, iS, iP )

           ! Get the depth - needed for wave speed calculation 
           CALL myDGSEM % GetBathymetryAtNode( iEl, iS, iP, bathy )
           h = bathy(1)

           ! Calculate the x and y fluxes
           xF = (/ h, ZERO /)
           yF = (/ ZERO, h /)

           !And now the contravariant flux
           pContFlux(iP,:) = -dyds*xF + dxds*yF

        ENDDO ! iP, Loop over the y-points

        ! North edge
         CALL myDGSEM % mesh % GetBoundaryNormalAtNode( iEl, normVec, normlength, iS, 3 ) 
         fR = ( hNorth(iS,1)*normVec )*normLength
           
         ! South edge
         CALL myDGSEM % mesh % GetBoundaryNormalAtNode( iEl, normVec, normlength, iS, 1 ) 
         fL = ( hSouth(iS,1)*normVec )*normLength

        ! At this x-level, calculate the y-DG-advective derivative
         DO iEq = 1, 2
            pContFluxDer(0:nP,iEq) = DGSystemDerivative(  nP, dMatY, qWeightY, fL(iEq), fR(iEq), &
                                                          pContFlux(:,iEq), lsouth, lnorth  )
         ENDDO
         
         DO iP = 0, nP ! Loop over the y-points
            DO iEq = 1, 2  ! Loop over the number of equations
                CALL myDGSEM % mesh % GetJacobianAtNode( iEl, J, iS, iP )
                tend(iS, iP, iEq) = ( tend(iS, iP, iEq)+ pContFluxDer(iP,iEq) )/J
            ENDDO ! Loop over the number of equations
         ENDDO ! Loop over the x-points

      
     ENDDO ! iP
         
     CALL myDGSEM % bathymetry(iEl) % SetSolutionWithVarID( 2, tend(:,:,1) ) ! set dhdx
     CALL myDGSEM % bathymetry(iEl) % SetSolutionWithVarID( 3, tend(:,:,2) ) ! set dhdy

 END SUBROUTINE CalculateBathymetryGradient_ShallowWater                     
!
!
!==================================================================================================!
!-------------------------------------- FILE I/O ROUTINES -----------------------------------------!
!==================================================================================================!
!
!
 SUBROUTINE CoarseToFine_ShallowWater( myDGSEM, iEl, x, y, depth, u, v, eta, vorticity )
 ! CoarseToFine
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS( ShallowWater ), INTENT(in) :: myDGSEM
   INTEGER, INTENT(in)               :: iEl
   REAL(prec), INTENT(out)           :: x(0:myDGSEM % nPlot, 0:myDGSEM % nPlot)
   REAL(prec), INTENT(out)           :: y(0:myDGSEM % nPlot, 0:myDGSEM % nPlot)
   REAL(prec), INTENT(out)           :: depth(0:myDGSEM % nPlot, 0:myDGSEM % nPlot)
   REAL(prec), INTENT(out)           :: u(0:myDGSEM % nPlot, 0:myDGSEM % nPlot)
   REAL(prec), INTENT(out)           :: v(0:myDGSEM % nPlot, 0:myDGSEM % nPlot)
   REAL(prec), INTENT(out)           :: eta(0:myDGSEM % nPlot, 0:myDGSEM % nPlot)
   REAL(prec), INTENT(out)           :: vorticity(0:myDGSEM % nPlot, 0:myDGSEM % nPlot)
   ! Local
   REAL(prec) :: localX(0:myDGSEM % nS, 0:myDGSEM % nP)
   REAL(prec) :: localY(0:myDGSEM % nS, 0:myDGSEM % nP)
   REAL(prec) :: localDepth(0:myDGSEM % nS, 0:myDGSEM % nP,1:3)
   REAL(prec) :: localU(0:myDGSEM % nS, 0:myDGSEM % nP)
   REAL(prec) :: localv(0:myDGSEM % nS, 0:myDGSEM % nP)
   REAL(prec) :: sol(0:myDGSEM % nS, 0:myDGSEM % nP,1:nSWeq)
   REAL(prec) :: localVorticity(0:myDGSEM % nS, 0:myDGSEM % nP)
   
      CALL myDGSEM % GetSolution( iEl, sol )
      CALL myDGSEM % GetVelocity( iEl, localU, localV )
      CALL myDGSEM % GetBathymetry( iEl, localDepth )
      CALL myDGSEM % mesh % GetPositions( iEl, localX, localY )
      
      CALL myDGSEM % dgStorage % interp % CoarseToFine( localX, &
                                                        myDGSEM % plMatS, &
                                                        myDGSEM % plMatP, &
                                                        myDGSEM % nPlot, &
                                                        myDGSEM % nPlot, &
                                                        x )
      
      CALL myDGSEM % dgStorage % interp % CoarseToFine( localY, &
                                                        myDGSEM % plMatS, &
                                                        myDGSEM % plMatP, &
                                                        myDGSEM % nPlot, &
                                                        myDGSEM % nPlot, &
                                                        y )
                                                    
      CALL myDGSEM % dgStorage % interp % CoarseToFine( localDepth, &
                                                        myDGSEM % plMatS, &
                                                        myDGSEM % plMatP, &
                                                        myDGSEM % nPlot, &
                                                        myDGSEM % nPlot, &
                                                        depth )
    
      CALL myDGSEM % dgStorage % interp % CoarseToFine( localU, &
                                                        myDGSEM % plMatS, &
                                                        myDGSEM % plMatP, &
                                                        myDGSEM % nPlot, &
                                                        myDGSEM % nPlot, &
                                                        u )
                                                        
      CALL myDGSEM % dgStorage % interp % CoarseToFine( localV, &
                                                        myDGSEM % plMatS, &
                                                        myDGSEM % plMatP, &
                                                        myDGSEM % nPlot, &
                                                        myDGSEM % nPlot, &
                                                        v )
                                                        
      CALL myDGSEM % dgStorage % interp % CoarseToFine( sol(:,:,3), &
                                                        myDGSEM % plMatS, &
                                                        myDGSEM % plMatP, &
                                                        myDGSEM % nPlot, &
                                                        myDGSEM % nPlot, &
                                                        eta )                                 
      
      CALL myDGSEM % dgStorage % interp % CoarseToFine( localVorticity, &
                                                        myDGSEM % plMatS, &
                                                        myDGSEM % plMatP, &
                                                        myDGSEM % nPlot, &
                                                        myDGSEM % nPlot, &
                                                        vorticity )
 END SUBROUTINE CoarseToFine_ShallowWater
!
!
!
 SUBROUTINE WriteTecplot_ShallowWater( myDGSEM, filename )
 ! S/R WriteTecplot
 !  
 ! =============================================================================================== !
 ! DECLARATIONS
  IMPLICIT NONE
  CLASS( ShallowWater ), INTENT(in)  :: myDGsem
  CHARACTER(*), INTENT(in), OPTIONAL :: filename
  !LOCAL
  REAL(prec)  :: x(0:myDGSEM % nPlot,0:myDGSEM % nPlot)
  REAL(prec)  :: y(0:myDGSEM % nPlot,0:myDGSEM % nPlot)
  REAL(prec)  :: depth(0:myDGSEM % nPlot,0:myDGSEM % nPlot)
  REAL(prec)  :: u(0:myDGSEM % nPlot,0:myDGSEM % nPlot)
  REAL(prec)  :: v(0:myDGSEM % nPlot,0:myDGSEM % nPlot)
  REAL(prec)  :: eta(0:myDGSEM % nPlot,0:myDGSEM % nPlot)
  REAL(prec)  :: vorticity(0:myDGSEM % nPlot,0:myDGSEM % nPlot)
  INTEGER     :: iS, iP, iEl, fUnit, nPlot
  CHARACTER(len=5) :: zoneID

    nPlot = myDGSEM % nPlot
    
    IF( PRESENT(filename) )THEN
       OPEN( UNIT=NEWUNIT(fUnit), &
             FILE= TRIM(filename)//'.tec', &
             FORM='formatted', &
             STATUS='replace')
    ELSE
       OPEN( UNIT=NEWUNIT(fUnit), &
             FILE= 'ShallowWater.tec', &
             FORM='formatted', &
             STATUS='replace')  
    ENDIF
    
    WRITE(fUnit,*) 'VARIABLES = "X", "Y", "Depth" "U", "V", "Eta", "Vorticity" '
 
    DO iEl = 1, myDGsem % mesh % nElems

       CALL myDGSEM % CoarseToFine( iEl, x, y, depth, u, v, eta, vorticity )
        WRITE(zoneID,'(I5.5)') iEl
        WRITE(fUnit,*) 'ZONE T="el'//trim(zoneID)//'", I=',nPlot+1,', J=', nPlot+1,',F=POINT'

        DO iP = 0, nPlot
           DO iS = 0, nPlot
              WRITE (fUnit,*)  x( iS, iP ), y( iS, iP ), depth(iS,iP), &
                               u(iS,iP), v(iS,iP), eta(iS,iP), &
                               vorticity(iS,iP)
           ENDDO
        ENDDO
        
    ENDDO

    CLOSE(UNIT=fUnit)

 END SUBROUTINE WriteTecplot_ShallowWater
!
!
!
 SUBROUTINE WritePickup_ShallowWater( myDGSEM, iter )
 ! S/R WritePickup
 ! 
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS( ShallowWater ), INTENT(in) :: myDGSEM
   INTEGER, INTENT(in)               :: iter
  ! LOCAL
   REAL(prec)    :: sol(0:myDGSEM % nS, 0:myDGSEM % nP,1:nSWeq)
   REAL(prec)    :: bathy(0:myDGSEM % nS, 0:myDGSEM % nP, 1:3)
   REAL(prec)    :: plv(0:myDGSEM % nS, 0:myDGSEM % nP)
   CHARACTER(10) :: iterChar
   INTEGER       :: iEl
   INTEGER       :: thisRec, fUnit
   INTEGER       :: iS, iP, iEq, nS, nP

     nS = myDGSEM % nS
     nP = myDGSEM % nP
     
     WRITE(iterChar,'(I10.10)') iter

     OPEN( UNIT=NEWUNIT(fUnit), &
           FILE='ShallowWater.'//iterChar//'.pickup', &
           FORM='unformatted',&
           ACCESS='direct',&
           STATUS='replace',&
           ACTION='WRITE',&
           CONVERT='big_endian',&
           RECL=prec*(nS+1)*(nP+1) )

     thisRec = 1 
     DO iEl = 1, myDGSEM % mesh % nElems
        
        CALL myDGSEM % GetSolution( iEl, sol )
        DO iEq = 1, nSWeq
           WRITE( fUnit, REC=thisRec )sol(:,:,iEq) 
           thisRec = thisRec+1
        ENDDO
        
        CALL myDGSEM % GetBathymetry( iEl, bathy )
        WRITE( fUnit, REC=thisRec ) bathy(:,:,1)
        thisRec = thisRec+1
        
        CALL myDGSEM % GetPlanetaryVorticity( iEl, plv )
        WRITE( fUnit, REC=thisRec ) plv
        thisRec = thisRec+1

     ENDDO

     CLOSE(UNIT=fUnit)

 END SUBROUTINE WritePickup_ShallowWater
!
!
!
  SUBROUTINE ReadPickup_ShallowWater( myDGSEM, iter )
 ! S/R ReadPickup
 ! 
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS( ShallowWater ), INTENT(inout) :: myDGSEM
   INTEGER, INTENT(in)                  :: iter
  ! LOCAL
   REAL(prec)    :: sol(0:myDGSEM % nS, 0:myDGSEM % nP,1:nSWeq)
   REAL(prec)    :: bathy(0:myDGSEM % nS, 0:myDGSEM % nP, 1:3)
   REAL(prec)    :: plv(0:myDGSEM % nS, 0:myDGSEM % nP)
   CHARACTER(10) :: iterChar
   INTEGER       :: iEl
   INTEGER       :: thisRec, fUnit
   INTEGER       :: iS, iP, iEq, nS, nP

     nS = myDGSEM % nS
     nP = myDGSEM % nP
     bathy = ZERO
     
     WRITE(iterChar,'(I10.10)') iter

     OPEN( UNIT=NEWUNIT(fUnit), &
           FILE='ShallowWater.'//iterChar//'.pickup', &
           FORM='unformatted',&
           ACCESS='direct',&
           STATUS='old',&
           ACTION='READ',&
           CONVERT='big_endian',&
           RECL=prec*(nS+1)*(nP+1) )
     
     thisRec = 1
     DO iEl = 1, myDGSEM % mesh % nElems
        
        DO iEq = 1, nSWeq
           READ( fUnit, REC=thisRec )sol(:,:,iEq) 
           thisRec = thisRec+1
        ENDDO
        CALL myDGSEM % SetSolution( iEl, sol )
        
        READ( fUnit, REC=thisRec ) bathy(:,:,1)
        thisRec = thisRec+1
        CALL myDGSEM % SetBathymetry( iEl, bathy )
        
        READ( fUnit, REC=thisRec ) plv
        thisRec = thisRec+1
        CALL myDGSEM % SetPlanetaryVorticity( iEl, plv )
        
        CALL myDGSEM % CalculateBathymetryAtBoundaries( iEl )
        
     ENDDO

     CLOSE(UNIT=fUnit)

     
 END SUBROUTINE ReadPickup_ShallowWater
!
!
! 
 END MODULE ConservativeShallowWaterClass



