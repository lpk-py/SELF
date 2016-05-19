! LBSW_Class.f90
! 
! Copyright 2015 Joseph Schoonover <schoonover.numerics@gmail.com>
! 
! LBSW_Class.f90 is part of the Spectral Element Libraries in Fortran (SELF).
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
 
 
MODULE LBSWClass
! ========================================= Logs ================================================= !
!2016-05-19  Joseph Schoonover  schoonover.numerics@gmail.com 
!
! //////////////////////////////////////////////////////////////////////////////////////////////// ! 

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
! src/highend/LinearBarotropicShelf/
USE RunParamsClass
USE Vorticity_Class
! Nocturnal Aviation classes and extensions
!USE FTTimerClass
!USE TIMING


! Parallel Libraries
USE OMP_LIB


IMPLICIT NONE
!
! Some parameters that are specific to this module

 ! INTEGER, parameter      :: nSWtimers = 5 ! the number of model timers
  REAL(prec), parameter   :: hDefault = ONE ! The default depth
     

! ================================================================================================ !
!     The LBSW is based on the DGSEM_2D framework with the addition of 
!     the AddOns data. Also, a few additional routines are added to facilitate
!     the calculation of "source" terms which include the coriolis force, relaxation terms, etc.
! ================================================================================================ !  

    TYPE LBSW
      INTEGER                               :: nEq, nPlot, nS, nP
      TYPE( Vorticity )                     :: vortInverter
      TYPE( QuadMesh )                      :: mesh
      TYPE( NodalStorage_2D )               :: dGStorage
      TYPE( DGSEMSolution_2D ), ALLOCATABLE :: sol(:)
      TYPE( DGSEMSolution_2D ), ALLOCATABLE :: vBar(:)
      REAL(prec), ALLOCATABLE               :: u(:,:,:), v(:,:,:) 
      REAL(prec), ALLOCATABLE               :: source(:,:,:)
      REAL(prec), ALLOCATABLE               :: pvfac(:,:,:)
      TYPE( RunParams )                     :: params
      REAL(prec), ALLOCATABLE               :: plMatS(:,:), plMatP(:,:)

      CONTAINS

      PROCEDURE :: Build => Build_LBSW
      PROCEDURE :: Trash => Trash_LBSW
      PROCEDURE :: BuildQuadMesh => BuildQuadMesh_LBSW
      
      ! DGSEMSolution_2DStorage_Class Wrapper Routines
      PROCEDURE :: GetSolution => GetSolution_LBSW
      PROCEDURE :: SetSolution => SetSolution_LBSW
      PROCEDURE :: GetSolutionAtNode => GetSolutionAtNode_LBSW
      PROCEDURE :: SetSolutionAtNode => SetSolutionAtNode_LBSW
      PROCEDURE :: GetSolutionAtNodeWithVarID => GetSolutionAtNodeWithVarID_LBSW
      PROCEDURE :: SetSolutionAtNodeWithVarID => SetSolutionAtNodeWithVarID_LBSW
      PROCEDURE :: GetTendency => GetTendency_LBSW
      PROCEDURE :: SetTendency => SetTendency_LBSW
      PROCEDURE :: GetTendencyAtNode => GetTendencyAtNode_LBSW
      PROCEDURE :: SetTendencyAtNode => SetTendencyAtNode_LBSW
      PROCEDURE :: GetTendencyWithVarID => GetTendencyWithVarID_LBSW
      PROCEDURE :: SetTendencyWithVarID => SetTendencyWithVarID_LBSW
      PROCEDURE :: GetTendencyAtNodeWithVarID => GetTendencyAtNodeWithVarID_LBSW
      PROCEDURE :: SetTendencyAtNodeWithVarID => SetTendencyAtNodeWithVarID_LBSW 
      PROCEDURE :: GetBoundaryFluxAtBoundary => GetBoundaryFluxAtBoundary_LBSW
      PROCEDURE :: SetBoundaryFluxAtBoundary => SetBoundaryFluxAtBoundary_LBSW
      PROCEDURE :: GetBoundaryFluxAtBoundaryNode => GetBoundaryFluxAtBoundaryNode_LBSW
      PROCEDURE :: SetBoundaryFluxAtBoundaryNode => SetBoundaryFluxAtBoundaryNode_LBSW
     
      
       ! Type Specific Routines
      PROCEDURE :: GlobalTimeDerivative => GlobalTimeDerivative_LBSW
      PROCEDURE :: ForwardStepRK3 => ForwardStepRK3_LBSW
      PROCEDURE :: EdgeFlux => EdgeFlux_LBSW
      PROCEDURE :: MappedTimeDerivative => MappedTimeDerivative_LBSW 
      PROCEDURE :: CalculateSolutionAtBoundaries => CalculateSolutionAtBoundaries_LBSW
      PROCEDURE :: CalculateVelocity => CalculateVelocity_LBSW

      
      PROCEDURE :: CoarseToFine => CoarseToFine_LBSW
      PROCEDURE :: WriteTecplot => WriteTecplot_LBSW
      PROCEDURE :: WritePickup => WritePickup_LBSW
      PROCEDURE :: ReadPickup => ReadPickup_LBSW

    END TYPE LBSW



 CONTAINS
!
!
!==================================================================================================!
!------------------------------- Manual Constructors/Destructors ----------------------------------!
!==================================================================================================!
!
!
 SUBROUTINE Build_LBSW( myDGSEM, forceNoPickup )
 ! S/R Build
 ! 
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS(LBSW), INTENT(inout)    :: myDGSEM
   INTEGER, INTENT(in), OPTIONAL :: forceNoPickup
   !LOCAL
   INTEGER :: iEl, nS, nP, nPlot
   REAL(prec), ALLOCATABLE ::sNew(:)

      CALL myDGSEM % params % Build( )
      CALL myDGSEM % vortInverter % Build( myDGSEM % params )

      nS = myDGSEM % params % polyDeg
      nP = nS
      myDGSEM % nS      = nS
      myDGSEM % nP      = nP
      myDGSEM % nEq     = 1
      myDGSEM % nPlot   = myDGSEM % params % nPlot
      
      nPlot = myDGSEM % params % nPlot
      
      ALLOCATE( sNew(0:nPlot) )
      CALL myDGSEM % dGStorage % Build( nS, nP, GAUSS_LOBATTO, DG )
 
      CALL myDGSEM % BuildQuadMesh( )
      
      ! Set up the solution, relaxation fields, bathymetry, and vorticity
      ALLOCATE( myDGSEM % sol(1:myDGSEM % mesh % nElems) )
      ALLOCATE( myDGSEM % vBar(1:myDGSEM % mesh % nElems) )
      ALLOCATE( myDGSEM % u(0:nS, 0:nP, 1:myDGSEM % mesh % nElems) )
      ALLOCATE( myDGSEM % v(0:nS, 0:nP, 1:myDGSEM % mesh % nElems) )
      ALLOCATE( myDGSEM % source(0:nS, 0:nP, 1:myDGSEM % mesh % nElems) )
      ALLOCATE( myDGSEM % pvfac(0:nS, 0:nP, 1:myDGSEM % mesh % nElems) )
      
      myDGSEM % u = ZERO
      myDGSEM % v = ZERO
      myDGSEM % source = ZERO
      myDGSEM % pvfac = ZERO

      ! Build and initialize the solution and the relaxation fields to zero
      DO iEl = 1, myDGSEM % mesh % nElems
         CALL myDGSEM % sol(iEl) % Build( nS, nP, myDGSEM % nEq )
         CALL myDGSEM % vBar(iEl) % Build( nS, nP, myDGSEM % nEq )
      ENDDO

      
      ! In the event that this is a pickup run, we'll READ the pickup file here for the solution 
      ! and the relaxation fields. A negative initial iterate can be used to specify to start 
      ! from zeros.
      
      IF( .NOT. PRESENT(forceNoPickup) )then
         
         ! This call reads the solution, the addons and the relaxation-parameter
         CALL myDGSEM % ReadPickup( myDGSEM % params % iterInit )

      ENDIF
      
      
      nPlot = myDGSEM % params % nPlot
      myDGSEM % nPlot = nPlot
      
      sNew = UniformPoints( -ONE, ONE, nPlot )
      ALLOCATE( myDGSEM % plMatS(0:nPlot,0:nS), myDGSEM % plMatP(0:nPlot,0:nP) )
      ! Build the plotting matrix
      CALL myDGSEM % dgStorage % interp % CalculateInterpolationMatrix( nPlot, nPlot, sNew, sNew, &
                                                                        myDGSEM % plMatS, &
                                                                        myDGSEM % plMatP )
                                                                        
      DEALLOCATE(sNew)
      
 END SUBROUTINE Build_LBSW
!
!
!
 SUBROUTINE Trash_LBSW( myDGSEM )
 ! S/R Trash
 ! 
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS(LBSW), INTENT(inout) :: myDGSEM
   ! LOCAL
   INTEGER :: iEl

     DO iEl = 1, myDGSEM % mesh % nElems
        CALL myDGSEM % sol(iEl) % Trash( )
        CALL myDGSEM % vBar(iEl) % Trash( )
     ENDDO

     CALL myDGSEM % dGStorage % Trash( )
     CALL myDGSEM % mesh % Trash( )
     CALL myDGSEM % vortInverter % Trash( )

     DEALLOCATE( myDGSEM % sol ) 
     DEALLOCATE( myDGSEM % u ) 
     DEALLOCATE( myDGSEM % v ) 
     DEALLOCATE( myDGSEM % vBar )
     DEALLOCATE( myDGSEM % vorticity )
     DEALLOCATE( myDGSEM % source )
     DEALLOCATE( myDGSEM % pvFac )
     DEALLOCATE( myDGSEM % plMatS, myDGSEM % plMatP )

 END SUBROUTINE Trash_LBSW
!
!
!
 SUBROUTINE BuildQuadMesh_LBSW( myDGSEM )
 ! S/R BuildQuadMesh
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS( LBSW ), INTENT(inout) :: myDGSEM
   ! LOCAL
   INTEGER :: iEdge, s1, e2
   
      PRINT*,'Module LBSWClass.f90 : S/R BuildQuadMesh :'
      IF( TRIM( myDGSEM % params % SpecMeshFile ) == nada )THEN
         PRINT*,' Loading default mesh.'
         CALL myDGSEM % mesh % LoadDefaultMesh( myDGSEM % dgStorage % interp, &
                                                myDGSEM % params % nXelem, &
                                                myDGSEM % params % nYelem )

         DO iEdge = 1, myDGSEM % mesh % nEdges
            CALL myDGSEM % mesh % GetEdgeSecondaryElementID( iEdge, e2 )
            CALL myDGSEM % mesh % GetEdgePrimaryElementSide( iEdge, s1 )
            IF( e2 == NO_NORMAL_FLOW )THEN
               IF( s1 == west .OR. s1 == north )THEN
                  CALL myDGSEM % mesh % SetEdgeSecondaryElementID( iEdge, InflowOne )
               ELSEIF( s1 == east .OR. s1 == south)THEN
                  CALL myDGSEM % mesh % SetEdgeSecondaryElementID( iEdge, InflowTwo )
               ENDIF
            ENDIF
         ENDDO

      ELSE
      ! Builds the lateral mesh
         PRINT*, 'Reading mesh from '//trim(myDGSEM % params % SpecMeshFile)//'.'
         CALL myDGSEM % mesh % ReadSpecMeshFile( myDGSEM % dgStorage % interp, &
                                                 myDGSEM % params % SpecMeshFile )
      ENDIF
      
      CALL myDGSEM % mesh % ScaleTheMesh( myDGSEM % dgStorage % interp, &
                                          myDGSEM % params % xScale, &
                                          myDGSEM % params % yScale )

      

 END SUBROUTINE BuildQuadMesh_LBSW
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
 SUBROUTINE GetSolution_LBSW( myDGSEM, iEl, theSolution  )
 ! S/R GetSolution
 !  
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS(LBSW), INTENT(in) :: myDGSEM
   INTEGER, INTENT(in)          :: iEl
   REAL(prec), INTENT(out)      :: theSolution(0:myDGSEM % nS, 0:myDGSEM % nP, 1:myDGSEM % nEq)

      CALL myDGSEM % sol(iEl) % GetSolution( theSolution )

 END SUBROUTINE GetSolution_LBSW
!
!
!
 SUBROUTINE SetSolution_LBSW( myDGSEM, iEl, theSolution  )
 ! S/R SetSolution
 !  
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS(LBSW), INTENT(inout) :: myDGSEM
   INTEGER, INTENT(in)             :: iEl
   REAL(prec), INTENT(in)          :: theSolution(0:myDGSEM % nS, 0:myDGSEM % nP, 1:myDGSEM % nEq)

      CALL myDGSEM % sol(iEl) % SetSolution( theSolution )

 END SUBROUTINE SetSolution_LBSW
!
!
!
 SUBROUTINE GetSolutionAtNode_LBSW( myDGSEM, iEl, i, j, theSolution  )
 ! S/R GetSolution
 !  
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS(LBSW), INTENT(in) :: myDGSEM
   INTEGER, INTENT(in)          :: iEl
   INTEGER, INTENT(in)          :: i, j
   REAL(prec), INTENT(out)      :: theSolution(1:myDGSEM % nEq)

      CALL myDGSEM % sol(iEl) % GetSolutionAtNode( i, j, theSolution )

 END SUBROUTINE GetSolutionAtNode_LBSW
!
!
!
 SUBROUTINE SetSolutionAtNode_LBSW( myDGSEM, iEl, i, j, theSolution  )
 ! S/R SetSolution
 !  
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS(LBSW), INTENT(inout) :: myDGSEM
   INTEGER, INTENT(in)             :: iEl
   INTEGER, INTENT(in)             :: i, j
   REAL(prec), INTENT(in)          :: theSolution(1:myDGSEM % nEq)

      CALL myDGSEM % sol(iEl) % SetSolutionAtNode( i, j, theSolution )

 END SUBROUTINE SetSolutionAtNode_LBSW
!
!
!
 SUBROUTINE GetSolutionAtNodeWithVarID_LBSW( myDGSEM, iEl, i, j, varID, theSolution  )
 ! S/R GetSolution
 !  
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS(LBSW), INTENT(in) :: myDGSEM
   INTEGER, INTENT(in)          :: iEl
   INTEGER, INTENT(in)          :: i, j, varID 
   REAL(prec), INTENT(out)      :: theSolution

      CALL myDGSEM % sol(iEl) % GetSolutionAtNodeWithVarID( i, j, varID, theSolution )

 END SUBROUTINE GetSolutionAtNodeWithVarID_LBSW
!
!
!
 SUBROUTINE SetSolutionAtNodeWithVarID_LBSW( myDGSEM, iEl, i, j, varID, theSolution  )
 ! S/R SetSolution
 !  
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS(LBSW), INTENT(inout) :: myDGSEM
   INTEGER, INTENT(in)             :: iEl
   INTEGER, INTENT(in)             :: i, j, varID 
   REAL(prec), INTENT(in)          :: theSolution

      CALL myDGSEM % sol(iEl) % SetSolutionAtNodeWithVarID( i, j, varID, theSolution )

 END SUBROUTINE SetSolutionAtNodeWithVarID_LBSW
!
! ----------------------------------------- Tendency --------------------------------------------- !
!
 SUBROUTINE GetTendency_LBSW( myDGSEM, iEl, theTend  )
 ! S/R GetTendency
 !  
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS(LBSW), INTENT(in) :: myDGSEM
   INTEGER, INTENT(in)          :: iEl
   REAL(prec), INTENT(out)      :: theTend(0:myDGSEM % nS,0:myDGSEM % nP, 1:myDGSEM % nEq)

      CALL myDGSEM % sol(iEl) % GetTendency( theTend )

 END SUBROUTINE GetTendency_LBSW
!
!
!
 SUBROUTINE SetTendency_LBSW( myDGSEM, iEl, theTend  )
 ! S/R SetTendency
 !  
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS(LBSW), INTENT(inout) :: myDGSEM
   INTEGER, INTENT(in)             :: iEl
   REAL(prec), INTENT(in)          :: theTend(0:myDGSEM % nS,0:myDGSEM % nP, 1:myDGSEM % nEq)

      CALL myDGSEM % sol(iEl) % SetTendency( theTend )   

 END SUBROUTINE SetTendency_LBSW
!
!
!
 SUBROUTINE GetTendencyAtNode_LBSW( myDGSEM, iEl, i, j, theTend  )
 ! S/R GetTendencyAtNode
 !  
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS(LBSW), INTENT(in) :: myDGSEM
   INTEGER, INTENT(in)          :: iEl
   INTEGER, INTENT(in)          :: i, j 
   REAL(prec), INTENT(out)      :: theTend(1:myDGSEM % nEq)

      CALL myDGSEM % sol(iEl) % GetTendencyAtNode( i, j, theTend )

 END SUBROUTINE GetTendencyAtNode_LBSW
!
!
!
 SUBROUTINE SetTendencyAtNode_LBSW( myDGSEM, iEl, i, j, theTend  )
 ! S/R SetTendencyAtNode
 !  
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS(LBSW), INTENT(inout) :: myDGSEM
   INTEGER, INTENT(in)             :: iEl
   INTEGER, INTENT(in)             :: i, j 
   REAL(prec), INTENT(in)          :: theTend(1:myDGSEM % nEq)

      CALL myDGSEM % sol(iEl) % SetTendencyAtNode( i, j, theTend )

 END SUBROUTINE SetTendencyAtNode_LBSW
!
!
!
SUBROUTINE GetTendencyWithVarID_LBSW( myDGSEM, iEl, varID, theTend  )
 ! S/R GetTendency
 !  
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS(LBSW), INTENT(in) :: myDGSEM
   INTEGER, INTENT(in)          :: iEl
   INTEGER, INTENT(in)          :: varID 
   REAL(prec), INTENT(out)      :: theTend(0:myDGSEM % nS,0:myDGSEM % nP)

      CALL myDGSEM % sol(iEl) % GetTendencyWithVarID( varID, theTend )

 END SUBROUTINE GetTendencyWithVarID_LBSW
!
!
!
 SUBROUTINE SetTendencyWithVarID_LBSW( myDGSEM, iEl, varID, theTend  )
 ! S/R SetTendency
 !  
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS(LBSW), INTENT(inout) :: myDGSEM
   INTEGER, INTENT(in)             :: iEl
   INTEGER, INTENT(in)             :: varID 
   REAL(prec), INTENT(in)          :: theTend(0:myDGSEM % nS,0:myDGSEM % nP)

      CALL myDGSEM % sol(iEl) % SetTendencyWithVarID( varID, theTend )

 END SUBROUTINE SetTendencyWithVarID_LBSW
!
!
!
 SUBROUTINE GetTendencyAtNodeWithVarID_LBSW( myDGSEM, iEl, i, j, varID, theTend  )
 ! S/R GetTendency
 !  
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS(LBSW), INTENT(in) :: myDGSEM
   INTEGER, INTENT(in)          :: iEl
   INTEGER, INTENT(in)          :: i, j, varID 
   REAL(prec), INTENT(out)      :: theTend

      CALL myDGSEM % sol(iEl) % GetTendencyAtNodeWithVarID( i, j, varID, theTend )

 END SUBROUTINE GetTendencyAtNodeWithVarID_LBSW
!
!
!
 SUBROUTINE SetTendencyAtNodeWithVarID_LBSW( myDGSEM, iEl, i, j, varID, theTend  )
 ! S/R SetTendency
 !  
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS(LBSW), INTENT(inout) :: myDGSEM
   INTEGER, INTENT(in)             :: iEl
   INTEGER, INTENT(in)             :: i, j, varID 
   REAL(prec), INTENT(in)          :: theTend

      CALL myDGSEM % sol(iEl) % SetTendencyAtNodeWithVarID( i, j, varID, theTend )

 END SUBROUTINE SetTendencyAtNodeWithVarID_LBSW
!
! ------------------------------- Boundary Solution ---------------------------------------------- !
!
 SUBROUTINE GetBoundarySolutionAtBoundary_LBSW( myDGSEM, iEl, boundary, theSolution  )
 ! S/R GetBoundarySolution
 !  
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS(LBSW), INTENT(in) :: myDGSEM
   INTEGER, INTENT(in)          :: iEl
   INTEGER, INTENT(in)          :: boundary
   REAL(prec), INTENT(out)      :: theSolution(0:myDGSEM % nS,1:myDGSEM % nEq)
   ! LOCAL
   INTEGER :: k

      !CALL myDGSEM % sol(iEl) % GetBoundarySolutionAtBoundary( boundary, theSolution )
      k = myDGSEM % mesh % sideMap(boundary)
      IF( boundary == SOUTH .OR. boundary == NORTH )THEN
         theSolution = myDGSEM % sol(iEl) % boundarySolution(:,k,:)
      ELSE
         theSolution = myDGSEM % sol(iEl) % boundarySolution(k,:,:)
      ENDIF

 END SUBROUTINE GetBoundarySolutionAtBoundary_LBSW
!
!
!
 SUBROUTINE SetBoundarySolutionAtBoundary_LBSW( myDGSEM, iEl, boundary, theSolution  )
 ! S/R SetBoundarySolution
 !  
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS(LBSW), INTENT(inout) :: myDGSEM
   INTEGER, INTENT(in)             :: iEl
   INTEGER, INTENT(in)             :: boundary
   REAL(prec), INTENT(in)          :: theSolution(0:myDGSEM % nS,1:myDGSEM % nEq)

      CALL myDGSEM % sol(iEl) % SetBoundarySolutionAtBoundary( boundary, theSolution )

 END SUBROUTINE SetBoundarySolutionAtBoundary_LBSW
!
! --------------------------------- Boundary Flux ------------------------------------------------ !
!
 SUBROUTINE GetBoundaryFluxAtBoundary_LBSW( myDGSEM, iEl, boundary, theFlux  )
 ! S/R GetBoundaryFlux
 !  
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS(LBSW), INTENT(in) :: myDGSEM
   INTEGER, INTENT(in)          :: iEl
   INTEGER, INTENT(in)          :: boundary
   REAL(prec), INTENT(out)      :: theFlux(0:myDGSEM % nS,1:myDGSEM % nEq)

      CALL myDGSEM % sol(iEl) % GetBoundaryFluxAtBoundary( boundary, theFlux )

 END SUBROUTINE GetBoundaryFluxAtBoundary_LBSW
!
!
!
 SUBROUTINE SetBoundaryFluxAtBoundary_LBSW( myDGSEM, iEl, boundary, theFlux  )
 ! S/R SetBoundaryFlux
 !  
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS(LBSW), INTENT(inout) :: myDGSEM
   INTEGER, INTENT(in)             :: iEl
   INTEGER, INTENT(in)             :: boundary
   REAL(prec), INTENT(in)          :: theFlux(0:myDGSEM % nS,1:myDGSEM % nEq)

      CALL myDGSEM % sol(iEl) % SetBoundaryFluxAtBoundary( boundary, theFlux )

 END SUBROUTINE SetBoundaryFluxAtBoundary_LBSW
!
!
!
 SUBROUTINE GetBoundaryFluxAtBoundaryNode_LBSW( myDGSEM, iEl, boundary, i, theFlux  )
 ! S/R GetBoundaryFlux
 !  
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS(LBSW), INTENT(in) :: myDGSEM
   INTEGER, INTENT(in)          :: iEl
   INTEGER, INTENT(in)          :: boundary, i
   REAL(prec), INTENT(out)      :: theFlux(1:myDGSEM % nEq)

      CALL myDGSEM % sol(iEl) % GetBoundaryFluxAtBoundaryNode( boundary, i, theFlux )

 END SUBROUTINE GetBoundaryFluxAtBoundaryNode_LBSW
!
!
!
 SUBROUTINE SetBoundaryFluxAtBoundaryNode_LBSW( myDGSEM, iEl, boundary, i, theFlux  )
 ! S/R SetBoundaryFlux
 !  
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS(LBSW), INTENT(inout) :: myDGSEM
   INTEGER, INTENT(in)             :: iEl
   INTEGER, INTENT(in)             :: boundary, i
   REAL(prec), INTENT(in)          :: theFlux(1:myDGSEM % nEq)

      CALL myDGSEM % sol(iEl) % SetBoundaryFluxAtBoundaryNode( boundary, i, theFlux )

 END SUBROUTINE SetBoundaryFluxAtBoundaryNode_LBSW
!
!
!==================================================================================================!
!--------------------------------- Type Specific Routines -----------------------------------------!
!==================================================================================================!
!
!
 SUBROUTINE ForwardStepRK3_LBSW( myDGSEM, tn )
 ! S/R ForwardStepRK3( 3rd order Runge-Kutta)
 ! 
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS(LBSW), INTENT(inout) :: myDGSEM
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
   INTEGER    :: m, nEq, iEqn, iEl, ioerr

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

           myDGSEM % vortInverter % source(:,:,iEl) = myDGSEM % sol(iEl) % solution(:,:,:,1)

         ENDDO ! iEl, loop over all of the elements

         ! Invert the vorticity equation for the stream function
         CALL myDGSEM % vortInverter % Solve( ioerr ) 
         ! Use the stream function to update the velocity field
         CALL myDGSEM % CalculateVelocity( )
         
         
      ENDDO ! m, loop over the RK3 steps
   
   
       
 END SUBROUTINE ForwardStepRK3_LBSW
!
!
!
 SUBROUTINE CalculateSolutionAtBoundaries_LBSW( myDGSEM, iEl )
 ! S/R CalculateSolutionAtBoundaries
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS(LBSW), INTENT(inout) :: myDGSEM
   INTEGER, INTENT(in)        :: iEl
   
      CALL myDGSEM % sol(iEl) % CalculateSolutionAtBoundaries( myDGSEM % dgStorage )
      
 END SUBROUTINE CalculateSolutionAtBoundaries_LBSW
!
!
!
 SUBROUTINE CalculateVelocity_LBSW( myDGSEM )
 ! S/R CalculateVelocity
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS(LBSW), INTENT(inout) :: myDGSEM
   ! Local
   INTEGER    :: iEl
   REAL(prec) :: sx(0:myDGSEM % nS, 0:myDGSEM % nP)
   REAL(prec) :: sy(0:myDGSEM % nS, 0:myDGSEM % nP)
   

      DO iEl = 1, myDGSEM % mesh % nElems
      
         CALL myDGSEM % vortInverter % CalculateGradient( iEl, &
                                                          myDGSEM % vortInverter % solution(:,:,iEl), &
                                                          sx, sy )
         myDGSEM % u(:,:,iEl) = -sy*myDGSEM % vortInverter % fluxCoeff(:,:,iEl)
         myDGSEM % v(:,:,iEl) = sx*myDGSEM % vortInverter % fluxCoeff(:,:,iEl)
      
      ENDDO

 END SUBROUTINE CalculateVelocity_LBSW
!
!
!
 SUBROUTINE GlobalTimeDerivative_LBSW( myDGSEM, tn ) 
 ! S/R GlobalTimeDerivative_LBSW
 ! 
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS(LBSW), INTENT(inout) :: myDGSEM
   REAL(prec), INTENT(in)          :: tn
   ! Local
   INTEGER :: iEl, iEdge

      ! CALL myDGSEM % LoadVelocityField( tn )
!$OMP PARALLEL


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

 
 END SUBROUTINE GlobalTimeDerivative_LBSW
!
!
! 
 SUBROUTINE EdgeFlux_LBSW( myDGSEM, iEdge, tn )
 ! S/R EdgeFlux
 ! 
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS( LBSW ), INTENT(inout) :: myDGSEM
   INTEGER, INTENT(in)               :: iEdge  
   REAL(prec), INTENT(in)            :: tn
   ! Local
   INTEGER :: k, iNode, nEq
   INTEGER :: e1, s1, e2, s2, iError, start, inc, nS, nP
   REAL(prec) :: flux(1:myDGSEM % nEq)
   REAL(prec) :: inState(0:myDGSEM % nS,1:myDGSEM % nEq)
   REAL(prec) :: exState(0:myDGSEM % nS,1:myDGSEM % nEq)
   REAL(prec) :: x, y
   REAL(prec) :: u(0:myDGSEM % nS), v(0:myDGSEM % nS) 
   REAL(prec) :: nHat(1:2), nHatLength
    
      nEQ = myDGSEM % params % nTracers

      CALL myDGSEM % dgStorage % GetNumberOfNodes( nS, nP )

      CALL myDGSEM % mesh % GetEdgePrimaryElementID( iEdge, e1 )
      CALL myDGSEM % mesh % GetEdgePrimaryElementSide( iEdge, s1 )
      CALL myDGSEM % mesh % GetEdgeSecondaryElementID( iEdge, e2 )
      CALL myDGSEM % mesh % GetEdgeSecondaryElementSide( iEdge, s2 )
      
      k = myDGSEM % mesh % sideMap( s1 )
      s2 = ABS(s2)

      CALL myDGSEM % mesh % GetEdgeStart( iEdge, start )
      CALL myDGSEM % mesh % GetEdgeIncrement( iEdge, inc )
      CALL myDGSEM % GetBoundarySolutionAtBoundary( e1, s1, inState )
      !CALL myDGSEM % GetVelocityAtBoundary( e1, s1, u, v )
      IF( s1 == NORTH .OR. s1 == SOUTH )THEN
         u = ZERO
         v = myDGSEM % vBar(e1) % Solution(:,k,1)
      ELSE
         u = ZERO
         v = myDGSEM % vBar(e1) % Solution(k,:,1)
      ENDIF

      IF( e2 > 0 )then ! this is an interior edge
      
         k = start-inc
         CALL myDGSEM % GetBoundarySolutionAtBoundary( e2, s2, exState )
         
         DO iNode = 0, nS ! Loop over the nodes
            
            CALL myDGSEM % mesh % GetBoundaryNormalAtNode( e1, nHat, nHatLength, iNode, s1 ) ! Get nHat
           
            ! Calculate the RIEMANN flux

            flux = RiemannSolver( inState(iNode,:), exState(k,:), &
                                  u(iNode), v(iNode), nHat, nEq )*nHatLength

            ! Store the flux for the elements which share this edge
            CALL myDGSEM % SetBoundaryFluxAtBoundaryNode( e1, s1, iNode, flux )
            CALL myDGSEM % SetBoundaryFluxAtBoundaryNode( e2, s2, k, -flux )
             
            k = k + inc

         ENDDO ! iNode, loop over the nodes
         
      ELSE ! this edge is a boundary edge

         DO iNode = 0, nS ! loop over the nodes on this edge
       
            CALL myDGSEM % mesh % GetBoundaryNormalAtNode( e1, nHat, nHatLength, iNode, s1 ) ! Get nHat
            
            ! Get the boundary point locations
            CALL myDGSEM % mesh % GetBoundaryLocationAtNode( e1, x, y, iNode, s1 )

            ! Calculate the external state
            exState(iNode,1:nEq) = GetExternalState( nEq, x, y, tn, e2, inState(iNode,1:nEq))

            ! Calculate the RIEMANN flux
            flux = RiemannSolver( inState(iNode,1:nEq), exState(iNode,1:nEq), &
                                  u(iNode), v(iNode), nHat, nEq )*nHatLength
 
            CALL myDGSEM % SetBoundaryFluxAtBoundaryNode( e1, s1, iNode, flux )

         ENDDO ! iNode, loop over the nodes on this edge

      ENDIF ! Choosing the edge type
 
 END SUBROUTINE EdgeFlux_LBSW
!
!
!
 SUBROUTINE MappedTimeDerivative_LBSW( myDGSEM, iEl, tn ) 
 ! S/R MappedTimeDerivative_LBSW
 ! 
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS(LBSW), INTENT(inout) :: myDGSEM
   INTEGER, INTENT(in)             :: iEl
   REAL(prec), INTENT(in)          :: tn
   ! LOCAL
   REAL(prec) :: dxds, dxdp, dyds, dydp
   REAL(prec) :: J, x, y
   REAL(prec) :: xF(1:myDGSEM % nEq), yF(1:myDGSEM % nEq), sol(1:myDGSEM % nEq)
   REAL(prec) :: dcdx(1:myDGSEM % nEq), dcdy(1:myDGSEM % nEq)
   REAL(prec) :: relaxSol(1:myDGSEM % nEq), rFac(1:myDGSEM % nEq)
   REAL(prec) :: u(0:myDGSEM % nS,0:myDGSEM % nS)
   REAL(prec) :: v(0:myDGSEM % nP,0:myDGSEM % nP)
   REAL(prec) :: pContFlux(0:myDGSEM % nS,1:myDGSEM % nEq)
   REAL(prec) :: sContFlux(0:myDGSEM % nS,1:myDGSEM % nEq)
   REAL(prec) :: pContFluxDer(0:myDGSEM % nS,1:myDGSEM % nEq)
   REAL(prec) :: sContFluxDer(0:myDGSEM % nS,1:myDGSEM % nEq)
   REAL(prec) :: tend(0:myDGSEM % nS,0:myDGSEM % nP,1:myDGSEM % nEq)
   REAL(prec) :: dMatX(0:myDGSEM % nS,0:myDGSEM % nS)
   REAL(prec) :: dMatY(0:myDGSEM % nP,0:myDGSEM % nP)
   REAL(prec) :: qWeightX(0:myDGSEM % nS)
   REAL(prec) :: qWeightY(0:myDGSEM % nP)
   REAL(prec) :: lsouth(0:myDGSEM % nP)
   REAL(prec) :: lnorth(0:myDGSEM % nP)
   REAL(prec) :: least(0:myDGSEM % nS)
   REAL(prec) :: lwest(0:myDGSEM % nS)
   REAL(prec) :: fL(1:myDGSEM % nEq), fR(1:myDGSEM % nEq)
   INTEGER    :: iS, iP, nS, nP, iEq

     CALL myDGSEM % dgStorage % GetNumberOfNodes( nS, nP )
     CALL myDGSEM % dgStorage % GetDerivativeMatrix( dMatX, dMatY )
     CALL myDGSEM % dgStorage % GetQuadratureWeights( qWeightX, qWeightY )
     CALL myDGSEM % dgStorage % GetSouthernInterpolants( lsouth )
     CALL myDGSEM % dgStorage % GetNorthernInterpolants( lnorth )
     CALL myDGSEM % dgStorage % GetEasternInterpolants( least )
     CALL myDGSEM % dgStorage % GetWesternInterpolants( lwest )
     u = ZERO
     v = myDGSEM % vBar(iEl) % Solution(:,:,1)


     DO iP = 0,nP ! Loop over the y-points

        DO iS = 0,nS ! Loop over the x-points

           ! Get the metric terms to calculate the contravariant flux
           CALL myDGSEM % mesh % GetCovariantMetricsAtNode( iEl, dxds, dxdp, dyds, dydp, iS, iP )

           ! Get the solution at x,y
           CALL myDGSEM % GetSolutionAtNode( iEl, iS, iP, sol )
           ! Calculate the x and y fluxes
           xF = XFlux(  tn, myDGSEM % nEq, u(iS,iP), v(iS,iP), sol )
           yF = YFlux(  tn, myDGSEM % nEq, u(iS,iP), v(iS,iP), sol )

           !And now the contravariant flux
           sContFlux(iS,:) = dydp*xF - dxdp*yF

        ENDDO ! iS, Loop over the x-points

        ! Get the numerical fluxes at the boundaries (should be calculated before
        ! this routine is called )
        
        ! Get the "west" flux at iP
        CALL myDGSEM % GetBoundaryFluxAtBoundaryNode( iEl, 4, iP, fL )

        ! Get the "east" flux at iP
        CALL myDGSEM % GetBoundaryFluxAtBoundaryNode( iEl, 2, iP, fR )

        ! At this y-level, calculate the DG-advective derivative
        DO iEq = 1, myDGSEM % nEq
           sContFluxDer(0:nS,iEq) = DGSystemDerivative(  nS, dMatX, qWeightX, fL(iEq), fR(iEq), &
                                                         sContFlux(:,iEq), lwest, least  )
        ENDDO
     
        DO iS = 0, nS ! Loop over the x-points
           DO iEq = 1, myDGSEM % nEq  ! Loop over the number of equations
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


           ! Calculate the x and y fluxes
           xF = XFlux(  tn, myDGSEM % nEq, u(iS,iP), v(iS,iP), sol )
           yF = YFlux(  tn, myDGSEM % nEq, u(iS,iP), v(iS,iP), sol )

           !And now the contravariant flux
           pContFlux(iP,:) = -dyds*xF + dxds*yF

        ENDDO ! iP, Loop over the y-points

        ! Get the numerical fluxes at the boundaries (should be calculated before
        ! this routine is called )
        
        ! Get the "south" flux at iS
        CALL  myDGSEM % GetBoundaryFluxAtBoundaryNode( iEl, 1, iS, fL )

        ! Get the "north" flux at iS
        CALL  myDGSEM % GetBoundaryFluxAtBoundaryNode( iEl, 3, iS, fR )


        ! At this x-level, calculate the y-DG-advective derivative
        DO iEq = 1, myDGSEM % nEq
            pContFluxDer(0:nP,iEq) = DGSystemDerivative(  nP, dMatY, qWeightY, fL(iEq), fR(iEq), &
                                                          pContFlux(:,iEq), lsouth, lnorth  )
        ENDDO

         
        DO iP = 0, nP ! Loop over the y-points
           CALL myDGSEM % mesh % GetJacobianAtNode( iEl, J, iS, iP )
           DO iEq = 1, myDGSEM % nEq  ! Loop over the number of equations
               tend(iS, iP, iEq) = ( tend(iS, iP, iEq)- pContFluxDer(iP,iEq) )/J
           ENDDO ! Loop over the number of equations
        ENDDO ! Loop over the x-points

      
     ENDDO ! iP
         
     ! COMPUTE THE SOURCE TERMS
     DO iS = 0, nS  ! Loop over the x-points
        DO iP = 0, nP ! Loop over the y-points

           CALL myDGSEM % GetSolutionAtNode( iEl, iS, iP, sol )
         
                  
           tend(iS,iP,:) = tend(iS,iP,:) + myDGSEM % source(:,:,iEl) + &
                                           myDGSEM % pvFac(:,:,iEl)*myDGSEM % u(:,:,iEl)

        ENDDO ! iP, Loop over the y-points
     ENDDO ! iS, Loop over the x-points

     CALL myDGSEM % SetTendency( iEl, tend )

 END SUBROUTINE MappedTimeDerivative_LBSW
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
 FUNCTION RiemannSolver( inState, outState, u, v, nHat, nEq ) RESULT( numFlux )
 ! FUNCTION RiemannSolver 
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   INTEGER    :: nEq
   REAL(prec) :: inState(1:nEq)
   REAL(prec) :: outState(1:nEq)
   REAL(prec) :: u, v
   REAL(prec) :: nHat(1:2) ! normal direction
   REAL(prec) :: numFlux(1:nEq)
   
   ! LOCAL
   REAL(prec) :: uN, dcdn(1:nEq)
   REAL(prec) :: jump(1:nEq), aS(1:nEq)
   REAL(prec) :: uNorm, fac


      jump = outState - inState
    
      uN =  u*nHat(1) + v*nHat(2)

      numFlux = HALF*( uN*( outState + inState ) - abs(uN)*jump )


 END FUNCTION RiemannSolver
!
!
! 
 FUNCTION XFlux(  tn, nEq, u, v, solAtX ) RESULT( fx )
 ! FUNCTION XFlux
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE 
   INTEGER    :: nEq
   REAL(prec) :: tn
   REAL(prec) :: u, v
   REAL(prec) :: solAtX(1:nEq)
   REAL(prec) :: fx(1:nEq)

      fx = solAtX*u

 END FUNCTION XFlux
!
!
!
 FUNCTION YFlux(  tn, nEq, u, v, solAtX ) RESULT( fy )
 ! FUNCTION YFlux
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE 
   INTEGER    :: nEq
   REAL(prec) :: tn
   REAL(prec) :: u, v
   REAL(prec) :: solAtX(1:nEq)
   REAL(prec) :: fy(1:nEq)

      fy = solAtX*v

 END FUNCTION YFlux
!
!
!                        
 FUNCTION GetExternalState( nEq, x, y, t, bcFlag, intState ) RESULT( extState )
 ! S/R GetExternalState
 ! 
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   INTEGER    :: nEq
   REAL(prec) :: intState(1:nEq)
   REAL(prec) :: x,y,t
   REAL(prec) :: nHat(1:2)
   INTEGER    :: bcFlag, formulation
   REAL(prec) :: extState(1:nEq)

       extState = ZERO 

!       IF( bcFlag == NO_NORMAL_FLOW )THEN
      
!          extState = ZERO 

!       ELSEIF( bcFlag == INFLOW )THEN
         
         
  
!       ELSE      
       
!         PRINT*, 'FUNCTION : GetExternalState : Invalid bcflag =', bcflag
!         extState = ZERO
         
!       ENDIF

        
 END FUNCTION GetExternalState
!
!
!==================================================================================================!
!-------------------------------------- FILE I/O ROUTINES -----------------------------------------!
!==================================================================================================!
!
!
 SUBROUTINE CoarseToFine_LBSW( myDGSEM, iEl, x, y, u, v, s )
 ! CoarseToFine
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS( LBSW ), INTENT(in) :: myDGSEM
   INTEGER, INTENT(in)            :: iEl
   REAL(prec), INTENT(out)        :: x(0:myDGSEM % nPlot, 0:myDGSEM % nPlot)
   REAL(prec), INTENT(out)        :: y(0:myDGSEM % nPlot, 0:myDGSEM % nPlot)
   REAL(prec), INTENT(out)        :: u(0:myDGSEM % nPlot, 0:myDGSEM % nPlot)
   REAL(prec), INTENT(out)        :: v(0:myDGSEM % nPlot, 0:myDGSEM % nPlot)
   REAL(prec), INTENT(out)        :: s(0:myDGSEM % nPlot, 0:myDGSEM % nPlot,1:myDGSEM % nEq)
   ! Local
   REAL(prec) :: localX(0:myDGSEM % nS, 0:myDGSEM % nP)
   REAL(prec) :: localY(0:myDGSEM % nS, 0:myDGSEM % nP)
   REAL(prec) :: localU(0:myDGSEM % nS, 0:myDGSEM % nP)
   REAL(prec) :: localv(0:myDGSEM % nS, 0:myDGSEM % nP)
   REAL(prec) :: sol(0:myDGSEM % nS, 0:myDGSEM % nP,1:myDGSEM % nEq)
   INTEGER    :: iEq


      CALL myDGSEM % GetSolution( iEl, sol )
      CALL myDGSEM % GetVelocity( iEl, localU, localV )
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
               
      DO iEq = 1, myDGSEM % nEq                                         
         CALL myDGSEM % dgStorage % interp % CoarseToFine( sol(:,:,iEq), &
                                                           myDGSEM % plMatS, &
                                                           myDGSEM % plMatP, &
                                                           myDGSEM % nPlot, &
                                                           myDGSEM % nPlot, &
                                                           s(:,:,iEq) )  

      ENDDO                               
      
 END SUBROUTINE CoarseToFine_LBSW
!
!
!
 SUBROUTINE WriteTecplot_LBSW( myDGSEM, filename )
 ! S/R WriteTecplot
 !  
 ! =============================================================================================== !
 ! DECLARATIONS
  IMPLICIT NONE
  CLASS( LBSW ), INTENT(in)  :: myDGsem
  CHARACTER(*), INTENT(in), OPTIONAL :: filename
  !LOCAL
  REAL(prec)  :: x(0:myDGSEM % nPlot,0:myDGSEM % nPlot)
  REAL(prec)  :: y(0:myDGSEM % nPlot,0:myDGSEM % nPlot)
  REAL(prec)  :: u(0:myDGSEM % nPlot,0:myDGSEM % nPlot)
  REAL(prec)  :: v(0:myDGSEM % nPlot,0:myDGSEM % nPlot)
  REAL(prec)  :: s(0:myDGSEM % nPlot,0:myDGSEM % nPlot, 1:myDGSEM % nEq)
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
             FILE= 'LBSW.tec', &
             FORM='formatted', &
             STATUS='replace')  
    ENDIF
    
    WRITE(fUnit,*) 'VARIABLES = "X", "Y", "U", "V", "Vorticity"'
 
    DO iEl = 1, myDGsem % mesh % nElems

       CALL myDGSEM % CoarseToFine( iEl, x, y, u, v, s )
        WRITE(zoneID,'(I5.5)') iEl
        WRITE(fUnit,*) 'ZONE T="el'//trim(zoneID)//'", I=',nPlot+1,', J=', nPlot+1,',F=POINT'

        DO iP = 0, nPlot
           DO iS = 0, nPlot
              WRITE (fUnit,*)  x( iS, iP ), y( iS, iP ), &
                               u(iS,iP), v(iS,iP), &
                               s(iS,iP,1)
           ENDDO
        ENDDO
        
    ENDDO

    CLOSE(UNIT=fUnit)

 END SUBROUTINE WriteTecplot_LBSW
!
!
!
 SUBROUTINE WritePickup_LBSW( myDGSEM, iter )
 ! S/R WritePickup
 ! 
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS( LBSW ), INTENT(in) :: myDGSEM
   INTEGER, INTENT(in)               :: iter
  ! LOCAL
   REAL(prec)    :: sol(0:myDGSEM % nS, 0:myDGSEM % nP,1:myDGSEM % nEq)
   REAL(prec)    :: u(0:myDGSEM % nS, 0:myDGSEM % nP)
   REAL(prec)    :: v(0:myDGSEM % nS, 0:myDGSEM % nP)
   CHARACTER(10) :: iterChar
   INTEGER       :: iEl
   INTEGER       :: thisRec, fUnit
   INTEGER       :: iS, iP, iEq, nS, nP

     nS = myDGSEM % nS
     nP = myDGSEM % nP
     
     WRITE(iterChar,'(I10.10)') iter

     OPEN( UNIT=NEWUNIT(fUnit), &
           FILE='LBSW.'//iterChar//'.pickup', &
           FORM='unformatted',&
           ACCESS='direct',&
           STATUS='replace',&
           ACTION='WRITE',&
           CONVERT='big_endian',&
           RECL=prec*(nS+1)*(nP+1) )

     thisRec = 1 
     DO iEl = 1, myDGSEM % mesh % nElems
        
        CALL myDGSEM % GetSolution( iEl, sol )
        DO iEq = 1, myDGSEM % nEq
           WRITE( fUnit, REC=thisRec )sol(:,:,iEq) 
           thisRec = thisRec+1
        ENDDO

        WRITE( fUnit, REC=thisRec ) myDGSEM % u(:,:,iEl)
        thisRec = thisRec+1
        WRITE( fUnit, REC=thisRec ) myDGSEM % v(:,:,iEl) 
        thisRec = thisRec+1
        WRITE( fUnit, REC=thisRec ) myDGSEM % vBar(:,:,iEl) 
        thisRec = thisRec+1
        WRITE( fUnit, REC=thisRec ) myDGSEM % source(:,:,iEl) 
        thisRec = thisRec+1
        WRITE( fUnit, REC=thisRec ) myDGSEM % pvFac(:,:,iEl) 
        thisRec = thisRec+1
        WRITE( fUnit, REC=thisRec ) myDGSEM % pvInverter % h(:,:,iEl) 
        thisRec = thisRec+1

     ENDDO

     CLOSE(UNIT=fUnit)

 END SUBROUTINE WritePickup_LBSW
!
!
!
  SUBROUTINE ReadPickup_LBSW( myDGSEM, iter )
 ! S/R ReadPickup
 ! 
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS( LBSW ), INTENT(inout) :: myDGSEM
   INTEGER, INTENT(in)                  :: iter
  ! LOCAL
   REAL(prec)    :: sol(0:myDGSEM % nS, 0:myDGSEM % nP,1:myDGSEM % nEq)
   REAL(prec)    :: u(0:myDGSEM % nS, 0:myDGSEM % nP)
   REAL(prec)    :: v(0:myDGSEM % nS, 0:myDGSEM % nP)
   CHARACTER(10) :: iterChar
   INTEGER       :: iEl
   INTEGER       :: thisRec, fUnit
   INTEGER       :: iS, iP, iEq, nS, nP

     nS = myDGSEM % nS
     nP = myDGSEM % nP
     
     WRITE(iterChar,'(I10.10)') iter

     OPEN( UNIT=NEWUNIT(fUnit), &
           FILE='LBSW.'//iterChar//'.pickup', &
           FORM='unformatted',&
           ACCESS='direct',&
           STATUS='old',&
           ACTION='READ',&
           CONVERT='big_endian',&
           RECL=prec*(nS+1)*(nP+1) )
     
     thisRec = 1
     DO iEl = 1, myDGSEM % mesh % nElems
        
        DO iEq = 1, myDGSEM % nEq
           READ( fUnit, REC=thisRec )sol(:,:,iEq) 
           thisRec = thisRec+1
        ENDDO
        CALL myDGSEM % SetSolution( iEl, sol )

        READ( fUnit, REC=thisRec ) myDGSEM % u(:,:,iEl)
        thisRec = thisRec+1
        READ( fUnit, REC=thisRec ) myDGSEM % v(:,:,iEl) 
        thisRec = thisRec+1
        READ( fUnit, REC=thisRec ) myDGSEM % vBar(:,:,iEl) 
        thisRec = thisRec+1
        READ( fUnit, REC=thisRec ) myDGSEM % source(:,:,iEl) 
        thisRec = thisRec+1
        READ( fUnit, REC=thisRec ) myDGSEM % pvFac(:,:,iEl) 
        thisRec = thisRec+1
        READ( fUnit, REC=thisRec ) myDGSEM % pvInverter % h(:,:,iEl) 
        thisRec = thisRec+1
        
     ENDDO

     CLOSE(UNIT=fUnit)

     CALL myDGSEM % pvInverter % ResetFluxCoefficient( )
     
 END SUBROUTINE ReadPickup_LBSW
!
!
! 
 END MODULE LBSWClass



