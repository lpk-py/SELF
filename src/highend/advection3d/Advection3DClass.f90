! AdvectionClass.f90 ( new with v2.1 - 18 March 2016)
! 
! ====================================== LICENSE ================================================= !
!
!  This program is free software; you can redistribute it and/or modify
!  it under the terms of the GNU General Public License as published by
!  the Free Software Foundation; either version 2 of the License, or
!  (at your option) any later version.
!  
!  This program is distributed in the hope that it will be useful,
!  but WITHOUT ANY WARRANTY; without even the implied warranty of
!  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!  GNU General Public License for more details.
!  
!  You should have received a copy of the GNU General Public License
!  along with this program; if not, write to the Free Software
!  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
!  MA 02110-1301, USA.
! 
! ==================================== Module History ============================================ ! 
! 
! o  (ver 2.1) January 2016
!
! ========================================= Logs ================================================= !
! yyyy-mm-dd  Joe  <joe@clay>
!
!   
! 
! //////////////////////////////////////////////////////////////////////////////////////////////// !


MODULE Advection3DClass
! AdvectionClass.f90
!
! schoonover.numerics@gmail.com
! 
! o (ver 2.1) March 2016
!
! This module provides a basic template for DGSEM in 3-D. The test system provided here is the
! advection-diffusion equation. DG for advection diffusion is recommended for scenarios in which
! the Peclet number is greater than one, ie, advection dominates diffusion.
!
! The Timing module (common/Timing.f90) is used to estimate the workload balance in the global 
! time-derivative routine.
! 
!  
! ================================================================================================ !

! src/common/
USE ModelPrecision
USE ConstantsDictionary
USE ModelFlags
USE CommonRoutines
USE Timing
! src/interp/
USE Lagrange_3D_Class
! src/nodal/
USE NodalStorage_3D_Class
USE DGSEMSolutionStorageClass_3D
! src/filters/
USE RollOffFilter3D_Class
! src/geometry/
USE FaceClass
USE HexElementClass  
USE HexMeshClass    
! src/highend/advection2d/
USE AdvectionParamsClass
!


! Parallel Libraries
USE OMP_LIB


IMPLICIT NONE
!
! Some parameters that are specific to this module
     

! ================================================================================================ !
!     The Advection is based on the DGSEM_3D framework. An Adaptive Roll-off Filter (ARF) is applied
!     to prevent the buildup of aliasing errors that can lead to spurious numerical oscillations 
!     and possibly instability. The adaptive algorithm is similar to [ Flad, Beck, and Munz (2016) ] 
!     but has been modified to use a roll-off filter and is applied to a passive tracer.
! ================================================================================================ !  

    TYPE Advection
      INTEGER                               :: nEq, nPlot, nS, nP, nQ
      TYPE( HexMesh )                       :: mesh
      TYPE( NodalStorage_3D )               :: dGStorage
      TYPE( DGSEMSolution_3D ), ALLOCATABLE :: sol(:)
      TYPE( DGSEMSolution_3D ), ALLOCATABLE :: velocity(:)
      TYPE( DGSEMSolution_3D ), ALLOCATABLE :: relax(:)
      TYPE( DGSEMSolution_3D ), ALLOCATABLE :: relaxFactor(:)
      TYPE( RollOffFilter3D )               :: modalFilter
      REAL(prec), ALLOCATABLE               :: E1(:,:), E2(:,:), lim(:,:,:)
      REAL(prec), ALLOCATABLE               :: dE1(:,:), dE2(:,:)
      TYPE( AdvectionParams )               :: params
      TYPE( MultiTimers )                   :: clocks
      REAL(prec), ALLOCATABLE               :: plMatS(:,:), plMatP(:,:), plMatQ(:,:)

      CONTAINS

      PROCEDURE :: Build => Build_Advection
      PROCEDURE :: Trash => Trash_Advection
      PROCEDURE :: BuildHexMesh => BuildHexMesh_Advection
      
      ! DGSEMSolution_2DStorage_Class Wrapper Routines
      PROCEDURE :: GetSolution => GetSolution_Advection
      PROCEDURE :: SetSolution => SetSolution_Advection
      PROCEDURE :: GetSolutionWithVarID => GetSolutionWithVarID_Advection
      PROCEDURE :: SetSolutionWithVarID => SetSolutionWithVarID_Advection
      PROCEDURE :: GetSolutionAtNode => GetSolutionAtNode_Advection
      PROCEDURE :: SetSolutionAtNode => SetSolutionAtNode_Advection
      PROCEDURE :: GetSolutionAtNodeWithVarID => GetSolutionAtNodeWithVarID_Advection
      PROCEDURE :: SetSolutionAtNodeWithVarID => SetSolutionAtNodeWithVarID_Advection
      PROCEDURE :: GetTendency => GetTendency_Advection
      PROCEDURE :: SetTendency => SetTendency_Advection
      PROCEDURE :: GetTendencyAtNode => GetTendencyAtNode_Advection
      PROCEDURE :: SetTendencyAtNode => SetTendencyAtNode_Advection
      PROCEDURE :: GetTendencyWithVarID => GetTendencyWithVarID_Advection
      PROCEDURE :: SetTendencyWithVarID => SetTendencyWithVarID_Advection
      PROCEDURE :: GetTendencyAtNodeWithVarID => GetTendencyAtNodeWithVarID_Advection
      PROCEDURE :: SetTendencyAtNodeWithVarID => SetTendencyAtNodeWithVarID_Advection
      PROCEDURE :: GetBoundarySolutionAtBoundary => GetBoundarySolutionAtBoundary_Advection
      PROCEDURE :: SetBoundarySolutionAtBoundary => SetBoundarySolutionAtBoundary_Advection
      PROCEDURE :: GetBoundaryFluxAtBoundary => GetBoundaryFluxAtBoundary_Advection
      PROCEDURE :: SetBoundaryFluxAtBoundary => SetBoundaryFluxAtBoundary_Advection
      PROCEDURE :: GetBoundaryFluxAtBoundaryNode => GetBoundaryFluxAtBoundaryNode_Advection
      PROCEDURE :: SetBoundaryFluxAtBoundaryNode => SetBoundaryFluxAtBoundaryNode_Advection
      ! Extensions of the DGSEMSolutionStorage_Class Wrapper Routines
      PROCEDURE :: GetVelocity => GetVelocity_Advection
      PROCEDURE :: SetVelocity => SetVelocity_Advection
      PROCEDURE :: GetVelocityAtBoundary => GetVelocityAtBoundary_Advection
      PROCEDURE :: SetVelocityAtBoundary => SetVelocityAtBoundary_Advection
      PROCEDURE :: GetRelaxationField => GetRelaxationField_Advection
      PROCEDURE :: SetRelaxationField => SetRelaxationField_Advection
      PROCEDURE :: GetRelaxationFieldAtNode => GetRelaxationFieldAtNode_Advection
      PROCEDURE :: SetRelaxationFieldAtNode => SetRelaxationFieldAtNode_Advection
      PROCEDURE :: GetRelaxationFactor => GetRelaxationFactor_Advection
      PROCEDURE :: SetRelaxationFactor => SetRelaxationFactor_Advection
      PROCEDURE :: GetRelaxationFactorAtNode => GetRelaxationFactorAtNode_Advection
      PROCEDURE :: SetRelaxationFactorAtNode => SetRelaxationFactorAtNode_Advection
     
      
       ! Type Specific Routines
      PROCEDURE :: GlobalTimeDerivative => GlobalTimeDerivative_Advection
      PROCEDURE :: DoTheAdaptiveFiltering => DoTheAdaptiveFiltering_Advection
      PROCEDURE :: ForwardStepRK3 => ForwardStepRK3_Advection
      PROCEDURE :: EdgeFlux => EdgeFlux_Advection
      PROCEDURE :: MappedTimeDerivative => MappedTimeDerivative_Advection 
      PROCEDURE :: CalculateSolutionAtBoundaries => CalculateSolutionAtBoundaries_Advection
      PROCEDURE :: CalculateVelocityAtBoundaries => CalculateVelocityAtBoundaries_Advection

      
      PROCEDURE :: CoarseToFine => CoarseToFine_Advection
      PROCEDURE :: WriteTecplot => WriteTecplot_Advection
      PROCEDURE :: WritePickup => WritePickup_Advection
      PROCEDURE :: ReadPickup => ReadPickup_Advection

    END TYPE Advection


 INTEGER, PARAMETER, PRIVATE :: nDim = 3

 CONTAINS
!
!
!==================================================================================================!
!------------------------------- Manual Constructors/Destructors ----------------------------------!
!==================================================================================================!
!
!
 SUBROUTINE Build_Advection( myDGSEM, forceNoPickup )
 ! S/R Build
 ! 
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS(Advection), INTENT(inout) :: myDGSEM
   INTEGER, INTENT(in), OPTIONAL   :: forceNoPickup
   !LOCAL
   INTEGER :: iEl, nS, nP, nQ, nPlot, tID
   REAL(prec), ALLOCATABLE ::sNew(:)

      CALL myDGSEM % params % Build( )

      !  /////// Clock setup //////// !
      CALL myDGSEM % clocks % Build( )

      tID = 1
      CALL myDGSEM % clocks % AddTimer( 'DoTheAdaptiveFilter', tID )
      tID = tID + 1

      CALL myDGSEM % clocks % AddTimer( 'CalculateSolutionAtBoundaries', tID )
      tID = tID + 1

      CALL myDGSEM % clocks % AddTimer( 'EdgeFlux', tID )
      tID = tID + 1
      
      CALL myDGSEM % clocks % AddTimer( 'MappedTimeDerivative', tID )
      tID = tID + 1

      CALL myDGSEM % clocks % AddTimer( 'File-I/O', 0 )
      ! ///////////////////////////// !

      nS = myDGSEM % params % polyDeg
      nP = nS
      nQ = nS
      myDGSEM % nS      = nS
      myDGSEM % nP      = nP
      myDGSEM % nQ      = nQ
      myDGSEM % nEq     = myDGSEM % params % nTracers
      myDGSEM % nPlot   = myDGSEM % params % nPlot
      
      nPlot = myDGSEM % params % nPlot

      ALLOCATE( sNew(0:nPlot) )
      CALL myDGSEM % dGStorage % Build( nS, nP, nQ, GAUSS, DG )
      CALL myDGSEM % modalFilter % Build( myDGSEM % dgStorage, &
                                          myDGSEM % params % nCutoff ) 

      CALL myDGSEM % BuildHexMesh( )

      ! Set up the solution, relaxation fields, bathymetry, and vorticity
      ALLOCATE( myDGSEM % sol(1:myDGSEM % mesh % nElems) )
      ALLOCATE( myDGSEM % relax(1:myDGSEM % mesh % nElems) )
      ALLOCATE( myDGSEM % relaxFactor(1:myDGSEM % mesh % nElems) )
      ALLOCATE( myDGSEM % velocity(1:myDGSEM % mesh % nElems) )
      ALLOCATE( myDGSEM % E1(1:myDGSEM % mesh % nElems,1:myDGSEM % nEq) )
      ALLOCATE( myDGSEM % E2(1:myDGSEM % mesh % nElems,1:myDGSEM % nEq) )
      ALLOCATE( myDGSEM % dE1(1:myDGSEM % mesh % nElems,1:myDGSEM % nEq) )
      ALLOCATE( myDGSEM % dE2(1:myDGSEM % mesh % nElems,1:myDGSEM % nEq) )
      ALLOCATE( myDGSEM % lim(1:myDGSEM % nEq,1:myDGSEM % mesh % nElems,1:3) )
      myDGSEM % E1  = ZERO
      myDGSEM % E2  = ZERO
      myDGSEM % dE1 = ZERO
      myDGSEM % dE2 = ZERO 
      myDGSEM % lim = ZERO
 
      ! Build and initialize the solution and the relaxation fields to zero
      DO iEl = 1, myDGSEM % mesh % nElems
         CALL myDGSEM % sol(iEl) % Build( nS, nP, nQ, myDGSEM % nEq )
         CALL myDGSEM % relax(iEl) % Build( nS, nP, nQ, myDGSEM % nEq )
         CALL myDGSEM % relaxFactor(iEl) % Build( nS, nP, nQ, myDGSEM % nEq ) 
         CALL myDGSEM % velocity(iEl) % Build( nS, nP, nQ, nDim )
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
      ALLOCATE( myDGSEM % plMatS(0:nPlot,0:nS), myDGSEM % plMatP(0:nPlot,0:nP), myDGSEM % plMatQ(0:nPlot,0:nQ) )
      ! Build the plotting matrix
      CALL myDGSEM % dgStorage % interp % CalculateInterpolationMatrix( nPlot, nPlot, nPlot, &
                                                                        sNew, sNew, sNew, &
                                                                        myDGSEM % plMatS, &
                                                                        myDGSEM % plMatP, &
                                                                        myDGSEM % plMatQ )
                                                                        
      DEALLOCATE(sNew)

      
      
 END SUBROUTINE Build_Advection
!
!
!
 SUBROUTINE Trash_Advection( myDGSEM )
 ! S/R Trash
 ! 
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS(Advection), INTENT(inout) :: myDGSEM
   ! LOCAL
   INTEGER :: iEl

     CALL myDGSEM % clocks % Write_MultiTimers( )

     DO iEl = 1, myDGSEM % mesh % nElems
        CALL myDGSEM % sol(iEl) % Trash( )
        CALL myDGSEM % relax(iEl) % Trash( )
        CALL myDGSEM % relaxFactor(iEl) % Trash( )
        CALL myDGSEM % velocity(iEl) % Trash( )
     ENDDO
    
     CALL myDGSEM % dGStorage % Trash( )
     
     CALL myDGSEM % mesh % Trash( )
     
     CALL myDGSEM % clocks % Trash( )
     CALL myDGSEM % modalFilter % Trash( )
     
     DEALLOCATE( myDGSEM % sol ) 
     DEALLOCATE( myDGSEM % relax )
     DEALLOCATE( myDGSEM % relaxFactor )
     DEALLOCATE( myDGSEM % velocity )
     DEALLOCATE( myDGSEM % E1, myDGSEM % E2 )
     DEALLOCATE( myDGSEM % lim )
     DEALLOCATE( myDGSEM % dE1, myDGSEM % dE2 )
     DEALLOCATE( myDGSEM % plMatS, myDGSEM % plMatP, myDGSEM % plMatQ )

 END SUBROUTINE Trash_Advection
!
!
!
 SUBROUTINE BuildHexMesh_Advection( myDGSEM )
 ! S/R BuildHexMesh
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS( Advection ), INTENT(inout) :: myDGSEM
   ! LOCAL
   INTEGER :: iEdge, s1, e2
   
      PRINT*,'Module AdvectionClass.f90 : S/R BuildQuadMesh :'
     ! IF( TRIM( myDGSEM % params % SpecMeshFile ) == nada )THEN
         PRINT*,' Loading default mesh.'
         CALL myDGSEM % mesh % LoadDefaultMesh( myDGSEM % dgStorage % interp, &
                                                myDGSEM % params % nXelem, &
                                                myDGSEM % params % nYelem, &
                                                myDGSEM % params % nZelem )

      CALL myDGSEM % mesh % ScaleTheMesh( myDGSEM % dgStorage % interp, &
                                          myDGSEM % params % xScale, &
                                          myDGSEM % params % yScale, &
                                          myDGSEM % params % zScale )

      

 END SUBROUTINE BuildHexMesh_Advection
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
 SUBROUTINE GetSolution_Advection( myDGSEM, iEl, theSolution  )
 ! S/R GetSolution
 !  
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS(Advection), INTENT(in) :: myDGSEM
   INTEGER, INTENT(in)          :: iEl
   REAL(prec), INTENT(out)      :: theSolution(0:myDGSEM % nS, 0:myDGSEM % nP, 0:myDGSEM % nQ,  1:myDGSEM % nEq)

      CALL myDGSEM % sol(iEl) % GetSolution( theSolution )

 END SUBROUTINE GetSolution_Advection
!
!
!
 SUBROUTINE SetSolution_Advection( myDGSEM, iEl, theSolution  )
 ! S/R SetSolution
 !  
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS(Advection), INTENT(inout) :: myDGSEM
   INTEGER, INTENT(in)             :: iEl
   REAL(prec), INTENT(in)          :: theSolution(0:myDGSEM % nS, 0:myDGSEM % nP, 0:myDGSEM % nQ,  1:myDGSEM % nEq)

      CALL myDGSEM % sol(iEl) % SetSolution( theSolution )

 END SUBROUTINE SetSolution_Advection
!
!
!
 SUBROUTINE GetSolutionWithVarID_Advection( myDGSEM, iEl, varID, theSolution  )
 ! S/R GetSolution
 !  
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS(Advection), INTENT(in) :: myDGSEM
   INTEGER, INTENT(in)          :: iEl, varID
   REAL(prec), INTENT(out)      :: theSolution(0:myDGSEM % nS, 0:myDGSEM % nP, 0:myDGSEM % nQ)

      CALL myDGSEM % sol(iEl) % GetSolutionWithVarID( varID, theSolution )

 END SUBROUTINE GetSolutionWithVarID_Advection
!
!
!
 SUBROUTINE SetSolutionWithVarID_Advection( myDGSEM, iEl, varID, theSolution  )
 ! S/R SetSolution
 !  
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS(Advection), INTENT(inout) :: myDGSEM
   INTEGER, INTENT(in)             :: iEl, varID
   REAL(prec), INTENT(in)          :: theSolution(0:myDGSEM % nS, 0:myDGSEM % nP, 0:myDGSEM % nQ)

      CALL myDGSEM % sol(iEl) % SetSolutionWithVarID( varID, theSolution )

 END SUBROUTINE SetSolutionWithVarID_Advection
!
!
!
 SUBROUTINE GetSolutionAtNode_Advection( myDGSEM, iEl, i, j, k, theSolution  )
 ! S/R GetSolution
 !  
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS(Advection), INTENT(in) :: myDGSEM
   INTEGER, INTENT(in)          :: iEl
   INTEGER, INTENT(in)          :: i, j, k
   REAL(prec), INTENT(out)      :: theSolution(1:myDGSEM % nEq)

      CALL myDGSEM % sol(iEl) % GetSolutionAtNode( i, j, k, theSolution )

 END SUBROUTINE GetSolutionAtNode_Advection
!
!
!
 SUBROUTINE SetSolutionAtNode_Advection( myDGSEM, iEl, i, j, k, theSolution  )
 ! S/R SetSolution
 !  
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS(Advection), INTENT(inout) :: myDGSEM
   INTEGER, INTENT(in)             :: iEl
   INTEGER, INTENT(in)             :: i, j, k
   REAL(prec), INTENT(in)          :: theSolution(1:myDGSEM % nEq)

      CALL myDGSEM % sol(iEl) % SetSolutionAtNode( i, j, k, theSolution )

 END SUBROUTINE SetSolutionAtNode_Advection
!
!
!
 SUBROUTINE GetVelocity_Advection( myDGSEM, iEl, u, v, w  )
 ! S/R GetVelocity
 !  
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS(Advection), INTENT(in) :: myDGSEM
   INTEGER, INTENT(in)          :: iEl
   REAL(prec), INTENT(out)      :: u(0:myDGSEM % nS,0:myDGSEM % nP, 0:myDGSEM % nQ)
   REAL(prec), INTENT(out)      :: v(0:myDGSEM % nS,0:myDGSEM % nP, 0:myDGSEM % nQ)
   REAL(prec), INTENT(out)      :: w(0:myDGSEM % nS,0:myDGSEM % nP, 0:myDGSEM % nQ)
   
      CALL myDGSEM % velocity(iEl) % GetSolutionWithVarID( 1, u )
      CALL myDGSEM % velocity(iEl) % GetSolutionWithVarID( 2, v )
      CALL myDGSEM % velocity(iEl) % GetSolutionWithVarID( 3, w )
      
 END SUBROUTINE GetVelocity_Advection
!
!
!
 SUBROUTINE SetVelocity_Advection( myDGSEM, iEl, u, v, w  )
 ! S/R SetVelocity
 !  
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS(Advection), INTENT(inout) :: myDGSEM
   INTEGER, INTENT(in)             :: iEl
   REAL(prec), INTENT(in)          :: u(0:myDGSEM % nS,0:myDGSEM % nP, 0:myDGSEM % nQ)
   REAL(prec), INTENT(in)          :: v(0:myDGSEM % nS,0:myDGSEM % nP, 0:myDGSEM % nQ)
   REAL(prec), INTENT(in)          :: w(0:myDGSEM % nS,0:myDGSEM % nP, 0:myDGSEM % nQ)

      CALL myDGSEM % velocity(iEl) % SetSolutionWithVarID( 1, u )
      CALL myDGSEM % velocity(iEl) % SetSolutionWithVarID( 2, v )
      CALL myDGSEM % velocity(iEl) % SetSolutionWithVarID( 3, w )

 END SUBROUTINE SetVelocity_Advection
!
!
!
 SUBROUTINE GetVelocityAtBoundary_Advection( myDGSEM, iEl, boundary, u, v, w  )
 ! S/R GetBoundarySolution
 !  
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS(Advection), INTENT(in) :: myDGSEM
   INTEGER, INTENT(in)          :: iEl
   INTEGER, INTENT(in)          :: boundary
   REAL(prec), INTENT(out)      :: u(0:myDGSEM % nS, 0:myDGSEM % nP)
   REAL(prec), INTENT(out)      :: v(0:myDGSEM % nS, 0:myDGSEM % nP)
   REAL(prec), INTENT(out)      :: w(0:myDGSEM % nS, 0:myDGSEM % nP)

   
      CALL myDGSEM % velocity(iEl) % GetBoundarySolutionAtBoundaryWithVarID( boundary, 1, u )
      CALL myDGSEM % velocity(iEl) % GetBoundarySolutionAtBoundaryWithVarID( boundary, 2, v )
      CALL myDGSEM % velocity(iEl) % GetBoundarySolutionAtBoundaryWithVarID( boundary, 3, w )

      
 END SUBROUTINE GetVelocityAtBoundary_Advection
!
!
!
 SUBROUTINE SetVelocityAtBoundary_Advection( myDGSEM, iEl, boundary, u, v, w  )
 ! S/R SetVelocityAtBoundary
 !  
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS(Advection), INTENT(inout) :: myDGSEM
   INTEGER, INTENT(in)             :: iEl
   INTEGER, INTENT(in)             :: boundary
   REAL(prec), INTENT(in)          :: u(0:myDGSEM % nS, 0:myDGSEM % nP)
   REAL(prec), INTENT(in)          :: v(0:myDGSEM % nS, 0:myDGSEM % nP)
   REAL(prec), INTENT(in)          :: w(0:myDGSEM % nS, 0:myDGSEM % nP)
   
      CALL myDGSEM % velocity(iEl) % SetBoundarySolutionAtBoundaryWithVarID( boundary, 1, u )
      CALL myDGSEM % velocity(iEl) % SetBoundarySolutionAtBoundaryWithVarID( boundary, 2, v )
      CALL myDGSEM % velocity(iEl) % SetBoundarySolutionAtBoundaryWithVarID( boundary, 3, w )
     
 END SUBROUTINE SetVelocityAtBoundary_Advection
!
!
!
 SUBROUTINE GetSolutionAtNodeWithVarID_Advection( myDGSEM, iEl, i, j, k, varID, theSolution  )
 ! S/R GetSolution
 !  
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS(Advection), INTENT(in) :: myDGSEM
   INTEGER, INTENT(in)          :: iEl
   INTEGER, INTENT(in)          :: i, j, k, varID 
   REAL(prec), INTENT(out)      :: theSolution

      CALL myDGSEM % sol(iEl) % GetSolutionAtNodeWithVarID( i, j, k, varID, theSolution )

 END SUBROUTINE GetSolutionAtNodeWithVarID_Advection
!
!
!
 SUBROUTINE SetSolutionAtNodeWithVarID_Advection( myDGSEM, iEl, i, j, k, varID, theSolution  )
 ! S/R SetSolution
 !  
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS(Advection), INTENT(inout) :: myDGSEM
   INTEGER, INTENT(in)             :: iEl
   INTEGER, INTENT(in)             :: i, j, k, varID 
   REAL(prec), INTENT(in)          :: theSolution

      CALL myDGSEM % sol(iEl) % SetSolutionAtNodeWithVarID( i, j, k, varID, theSolution )

 END SUBROUTINE SetSolutionAtNodeWithVarID_Advection
!
! ----------------------------------------- Tendency --------------------------------------------- !
!
 SUBROUTINE GetTendency_Advection( myDGSEM, iEl, theTend  )
 ! S/R GetTendency
 !  
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS(Advection), INTENT(in) :: myDGSEM
   INTEGER, INTENT(in)          :: iEl
   REAL(prec), INTENT(out)      :: theTend(0:myDGSEM % nS, 0:myDGSEM % nP, 0:myDGSEM % nQ, 1:myDGSEM % nEq)

      CALL myDGSEM % sol(iEl) % GetTendency( theTend )

 END SUBROUTINE GetTendency_Advection
!
!
!
 SUBROUTINE SetTendency_Advection( myDGSEM, iEl, theTend  )
 ! S/R SetTendency
 !  
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS(Advection), INTENT(inout) :: myDGSEM
   INTEGER, INTENT(in)             :: iEl
   REAL(prec), INTENT(in)          :: theTend(0:myDGSEM % nS, 0:myDGSEM % nP, 0:myDGSEM % nQ, 1:myDGSEM % nEq)

      CALL myDGSEM % sol(iEl) % SetTendency( theTend )   

 END SUBROUTINE SetTendency_Advection
!
!
!
 SUBROUTINE GetTendencyAtNode_Advection( myDGSEM, iEl, i, j, k, theTend  )
 ! S/R GetTendencyAtNode
 !  
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS(Advection), INTENT(in) :: myDGSEM
   INTEGER, INTENT(in)          :: iEl
   INTEGER, INTENT(in)          :: i, j, k 
   REAL(prec), INTENT(out)      :: theTend(1:myDGSEM % nEq)

      CALL myDGSEM % sol(iEl) % GetTendencyAtNode( i, j, k, theTend )

 END SUBROUTINE GetTendencyAtNode_Advection
!
!
!
 SUBROUTINE SetTendencyAtNode_Advection( myDGSEM, iEl, i, j, k, theTend  )
 ! S/R SetTendencyAtNode
 !  
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS(Advection), INTENT(inout) :: myDGSEM
   INTEGER, INTENT(in)             :: iEl
   INTEGER, INTENT(in)             :: i, j, k 
   REAL(prec), INTENT(in)          :: theTend(1:myDGSEM % nEq)

      CALL myDGSEM % sol(iEl) % SetTendencyAtNode( i, j, k, theTend )

 END SUBROUTINE SetTendencyAtNode_Advection
!
!
!
SUBROUTINE GetTendencyWithVarID_Advection( myDGSEM, iEl, varID, theTend  )
 ! S/R GetTendency
 !  
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS(Advection), INTENT(in) :: myDGSEM
   INTEGER, INTENT(in)          :: iEl
   INTEGER, INTENT(in)          :: varID 
   REAL(prec), INTENT(out)      :: theTend(0:myDGSEM % nS,0:myDGSEM % nP, 0:myDGSEM % nQ)

      CALL myDGSEM % sol(iEl) % GetTendencyWithVarID( varID, theTend )

 END SUBROUTINE GetTendencyWithVarID_Advection
!
!
!
 SUBROUTINE SetTendencyWithVarID_Advection( myDGSEM, iEl, varID, theTend  )
 ! S/R SetTendency
 !  
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS(Advection), INTENT(inout) :: myDGSEM
   INTEGER, INTENT(in)             :: iEl
   INTEGER, INTENT(in)             :: varID 
   REAL(prec), INTENT(in)          :: theTend(0:myDGSEM % nS,0:myDGSEM % nP, 0:myDGSEM % nQ)

      CALL myDGSEM % sol(iEl) % SetTendencyWithVarID( varID, theTend )

 END SUBROUTINE SetTendencyWithVarID_Advection
!
!
!
 SUBROUTINE GetTendencyAtNodeWithVarID_Advection( myDGSEM, iEl, i, j, k, varID, theTend  )
 ! S/R GetTendency
 !  
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS(Advection), INTENT(in) :: myDGSEM
   INTEGER, INTENT(in)          :: iEl
   INTEGER, INTENT(in)          :: i, j, k, varID 
   REAL(prec), INTENT(out)      :: theTend

      CALL myDGSEM % sol(iEl) % GetTendencyAtNodeWithVarID( i, j, k, varID, theTend )

 END SUBROUTINE GetTendencyAtNodeWithVarID_Advection
!
!
!
 SUBROUTINE SetTendencyAtNodeWithVarID_Advection( myDGSEM, iEl, i, j, k, varID, theTend  )
 ! S/R SetTendency
 !  
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS(Advection), INTENT(inout) :: myDGSEM
   INTEGER, INTENT(in)             :: iEl
   INTEGER, INTENT(in)             :: i, j, k, varID 
   REAL(prec), INTENT(in)          :: theTend

      CALL myDGSEM % sol(iEl) % SetTendencyAtNodeWithVarID( i, j, k, varID, theTend )

 END SUBROUTINE SetTendencyAtNodeWithVarID_Advection
!
! ------------------------------- Boundary Solution ---------------------------------------------- !
!
 SUBROUTINE GetBoundarySolutionAtBoundary_Advection( myDGSEM, iEl, boundary, theSolution  )
 ! S/R GetBoundarySolution
 !  
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS(Advection), INTENT(in) :: myDGSEM
   INTEGER, INTENT(in)          :: iEl
   INTEGER, INTENT(in)          :: boundary
   REAL(prec), INTENT(out)      :: theSolution(0:myDGSEM % nS, 0:myDGSEM % nP, 1:myDGSEM % nEq)

      CALL myDGSEM % sol(iEl) % GetBoundarySolutionAtBoundary( boundary, theSolution )

 END SUBROUTINE GetBoundarySolutionAtBoundary_Advection
!
!
!
 SUBROUTINE SetBoundarySolutionAtBoundary_Advection( myDGSEM, iEl, boundary, theSolution  )
 ! S/R SetBoundarySolution
 !  
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS(Advection), INTENT(inout) :: myDGSEM
   INTEGER, INTENT(in)             :: iEl
   INTEGER, INTENT(in)             :: boundary
   REAL(prec), INTENT(in)          :: theSolution(0:myDGSEM % nS, 0:myDGSEM % nP, 1:myDGSEM % nEq)

      CALL myDGSEM % sol(iEl) % SetBoundarySolutionAtBoundary( boundary, theSolution )

 END SUBROUTINE SetBoundarySolutionAtBoundary_Advection
!
! --------------------------------- Boundary Flux ------------------------------------------------ !
!
 SUBROUTINE GetBoundaryFluxAtBoundary_Advection( myDGSEM, iEl, boundary, theFlux  )
 ! S/R GetBoundaryFlux
 !  
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS(Advection), INTENT(in) :: myDGSEM
   INTEGER, INTENT(in)          :: iEl
   INTEGER, INTENT(in)          :: boundary
   REAL(prec), INTENT(out)      :: theFlux(0:myDGSEM % nS, 0:myDGSEM % nP, 1:myDGSEM % nEq)

      CALL myDGSEM % sol(iEl) % GetBoundaryFluxAtBoundary( boundary, theFlux )

 END SUBROUTINE GetBoundaryFluxAtBoundary_Advection
!
!
!
 SUBROUTINE SetBoundaryFluxAtBoundary_Advection( myDGSEM, iEl, boundary, theFlux  )
 ! S/R SetBoundaryFlux
 !  
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS(Advection), INTENT(inout) :: myDGSEM
   INTEGER, INTENT(in)             :: iEl
   INTEGER, INTENT(in)             :: boundary
   REAL(prec), INTENT(in)          :: theFlux(0:myDGSEM % nS, 0:myDGSEM % nP, 1:myDGSEM % nEq)

      CALL myDGSEM % sol(iEl) % SetBoundaryFluxAtBoundary( boundary, theFlux )

 END SUBROUTINE SetBoundaryFluxAtBoundary_Advection
!
!
!
 SUBROUTINE GetBoundaryFluxAtBoundaryNode_Advection( myDGSEM, iEl, boundary, i, j, theFlux  )
 ! S/R GetBoundaryFlux
 !  
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS(Advection), INTENT(in) :: myDGSEM
   INTEGER, INTENT(in)          :: iEl
   INTEGER, INTENT(in)          :: boundary, i, j
   REAL(prec), INTENT(out)      :: theFlux(1:myDGSEM % nEq)

      CALL myDGSEM % sol(iEl) % GetBoundaryFluxAtBoundaryNode( boundary, i, j, theFlux )

 END SUBROUTINE GetBoundaryFluxAtBoundaryNode_Advection
!
!
!
 SUBROUTINE SetBoundaryFluxAtBoundaryNode_Advection( myDGSEM, iEl, boundary, i, j, theFlux  )
 ! S/R SetBoundaryFlux
 !  
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS(Advection), INTENT(inout) :: myDGSEM
   INTEGER, INTENT(in)             :: iEl
   INTEGER, INTENT(in)             :: boundary, i, j
   REAL(prec), INTENT(in)          :: theFlux(1:myDGSEM % nEq)

      CALL myDGSEM % sol(iEl) % SetBoundaryFluxAtBoundaryNode( boundary, i, j, theFlux )

 END SUBROUTINE SetBoundaryFluxAtBoundaryNode_Advection
!
! --------------------------------- Relaxation Field --------------------------------------------- !
!
 SUBROUTINE GetRelaxationField_Advection( myDGSEM, iEl, relaxField  )
 ! S/R GetRelaxationField
 !  
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS(Advection), INTENT(in) :: myDGSEM
   INTEGER, INTENT(in)          :: iEl
   REAL(prec), INTENT(out)      :: relaxField(0:myDGSEM % nS, 0:myDGSEM % nP, 0:myDGSEM % nQ,  1:myDGSEM % nEq)
 
      CALL myDGSEM % relax(iEl) % GetSolution( relaxField )
      
 END SUBROUTINE GetRelaxationField_Advection
!
!
!
 SUBROUTINE SetRelaxationField_Advection( myDGSEM, iEl, relaxField  )
 ! S/R SetRelaxationField
 !  
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS(Advection), INTENT(inout) :: myDGSEM
   INTEGER, INTENT(in)             :: iEl
   REAL(prec), INTENT(in)          :: relaxField(0:myDGSEM % nS, 0:myDGSEM % nP, 0:myDGSEM % nQ,  1:myDGSEM % nEq)

      CALL myDGSEM % relax(iEl) % SetSolution( relaxField )

 END SUBROUTINE SetRelaxationField_Advection
!
!
!
 SUBROUTINE GetRelaxationFieldAtNode_Advection( myDGSEM, iEl, i, j, k, relaxField  )
 ! S/R GetRelaxationField
 !  
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS(Advection), INTENT(in) :: myDGSEM
   INTEGER, INTENT(in)          :: iEl, i, j, k
   REAL(prec), INTENT(out)      :: relaxField(1:myDGSEM % nEq)
 
      CALL myDGSEM % relax(iEl) % GetSolutionAtNode( i, j, k, relaxField )
      
 END SUBROUTINE GetRelaxationFieldAtNode_Advection
!
!
!
 SUBROUTINE SetRelaxationFieldAtNode_Advection( myDGSEM, iEl, i, j, k, relaxField  )
 ! S/R SetRelaxationFieldAtNode
 !  
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS(Advection), INTENT(inout) :: myDGSEM
   INTEGER, INTENT(in)             :: iEl, i, j, k
   REAL(prec), INTENT(in)          :: relaxField(1:myDGSEM % nEq)

      CALL myDGSEM % relax(iEl) % SetSolutionAtNode( i, j, k, relaxField )

 END SUBROUTINE SetRelaxationFieldAtNode_Advection
!
!
!
 SUBROUTINE GetRelaxationFactor_Advection( myDGSEM, iEl, rFac  )
 ! S/R GetRelaxationFactor
 !  
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS(Advection), INTENT(in) :: myDGSEM
   INTEGER, INTENT(in)          :: iEl
   REAL(prec), INTENT(out)      :: rFac(0:myDGSEM % nS, 0:myDGSEM % nP, 0:myDGSEM % nQ,  1:myDGSEM % nEq)

      CALL myDGSEM % relax(iEl) % GetSolution( rFac )
          
 END SUBROUTINE GetRelaxationFactor_Advection
!
!
!
 SUBROUTINE SetRelaxationFactor_Advection( myDGSEM, iEl, rFac  )
 ! S/R SetRelaxationFactor
 !  
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS(Advection), INTENT(inout) :: myDGSEM
   INTEGER, INTENT(in)             :: iEl
   REAL(prec), INTENT(in)          :: rFac(0:myDGSEM % nS, 0:myDGSEM % nP, 0:myDGSEM % nQ,  1:myDGSEM % nEq)
   
      CALL myDGSEM % relax(iEl) % SetSolution( rFac )

 END SUBROUTINE SetRelaxationFactor_Advection
!
!
!
 SUBROUTINE GetRelaxationFactorAtNode_Advection( myDGSEM, iEl, i, j, k, rFac  )
 ! S/R GetRelaxationFactor
 !  
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS(Advection), INTENT(in) :: myDGSEM
   INTEGER, INTENT(in)          :: iEl, i, j, k
   REAL(prec), INTENT(out)      :: rFac(1:myDGSEM % nEq)
 
      CALL myDGSEM % relax(iEl) % GetSolutionAtNode( i, j, k, rFac )
      
 END SUBROUTINE GetRelaxationFactorAtNode_Advection
!
!
!
 SUBROUTINE SetRelaxationFactorAtNode_Advection( myDGSEM, iEl, i, j, k, rFac  )
 ! S/R SetRelaxationFactorAtNode
 !  
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS(Advection), INTENT(inout) :: myDGSEM
   INTEGER, INTENT(in)             :: iEl, i, j, k
   REAL(prec), INTENT(in)          :: rFac(1:myDGSEM % nEq)

      CALL myDGSEM % relax(iEl) % SetSolutionAtNode( i, j, k, rFac )

 END SUBROUTINE SetRelaxationFactorAtNode_Advection
!
!
!==================================================================================================!
!--------------------------------- Type Specific Routines -----------------------------------------!
!==================================================================================================!
!
!
 SUBROUTINE CalculateSolutionAtBoundaries_Advection( myDGSEM, iEl )
 ! S/R CalculateSolutionAtBoundaries
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS(Advection), INTENT(inout) :: myDGSEM
   INTEGER, INTENT(in)             :: iEl
   
      CALL myDGSEM % sol(iEl) % CalculateSolutionAtBoundaries( myDGSEM % dgStorage )
      
 END SUBROUTINE CalculateSolutionAtBoundaries_Advection
!
!
!
 SUBROUTINE CalculateVelocityAtBoundaries_Advection( myDGSEM, iEl )
 ! S/R CalculateVelocityAtBoundaries
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS(Advection), INTENT(inout) :: myDGSEM
   INTEGER, INTENT(in)             :: iEl
   
      CALL myDGSEM % velocity(iEl) % CalculateSolutionAtBoundaries( myDGSEM % dgStorage )
      
 END SUBROUTINE CalculateVelocityAtBoundaries_Advection
!
!
!
 SUBROUTINE ForwardStepRK3_Advection( myDGSEM, tn )
 ! S/R ForwardStepRK3( 3rd order Runge-Kutta)
 ! 
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS(Advection), INTENT(inout) :: myDGSEM
   REAL(prec), INTENT(in)                   :: tn
   ! LOCAL
   REAL(prec) :: t, dt
   REAL(prec) :: G2D(0:myDGSEM % nS,&
                     0:myDGSEM % nP,&
                     0:myDGSEM % nQ,&
                     1:myDGSEM % nEq,&
                     1:myDGSEM % mesh % nElems) 
   REAL(prec) :: dSdt(0:myDGSEM % nS,&
                      0:myDGSEM % nP,&
                      0:myDGSEM % nQ,&
                      1:myDGSEM % nEq )
   INTEGER    :: m, nEq, iEqn, iEl

     nEq = myDGSEM % nEq 
     dt = myDGSEM % params % dt
     G2D = ZERO
    
     DO m = 1,3 ! Loop over RK3 steps

        t = tn + rk3_b(m)*dt
        ! Calculate the tendency
        CALL myDGSEM % GlobalTimeDerivative( t, m )
        
        DO iEl = 1, myDGSEM % mesh % nElems ! Loop over all of the elements

           CALL myDGSEM % GetTendency( iEl, dSdt )
           G2D(:,:,:,:,iEl) = rk3_a(m)*G2D(:,:,:,:,iEl) + dSdt

           myDGSEM % sol(iEl) % solution = myDGSEM % sol(iEl) % solution + rk3_g(m)*dt*G2D(:,:,:,:,iEl)

         ENDDO ! iEl, loop over all of the elements
         
      ENDDO ! m, loop over the RK3 steps
   
   
       
 END SUBROUTINE ForwardStepRK3_Advection
!
!
!
 SUBROUTINE GlobalTimeDerivative_Advection( myDGSEM, tn, m ) 
 ! S/R GlobalTimeDerivative_Advection
 ! 
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS(Advection), INTENT(inout) :: myDGSEM
   REAL(prec), INTENT(in)          :: tn
   INTEGER, INTENT(in)             :: m
   ! Local
   INTEGER :: iEl, iFace, tID

      ! CALL myDGSEM % LoadVelocityField( tn )

!$OMP PARALLEL

!$OMP DO
      tID = 1
      CALL myDGSEM % clocks % StartThisTimer( tID )
      DO iEl = 1, myDGSEM % mesh % nElems
         CALL myDGSEM % DoTheAdaptiveFiltering( iEl, m ) 
      ENDDO 
      CALL myDGSEM % clocks % StopThisTimer( tID )
      CALL myDGSEM % clocks % AccumulateTimings( )
      tID = tID + 1
!$OMP END DO 
!$OMP FLUSH( myDGSEM )


      ! Calculate the solution at the boundaries
!$OMP DO
      CALL myDGSEM % clocks % StartThisTimer( tID )
      DO iEl = 1, myDGSEM % mesh % nElems

         CALL myDGSEM % CalculateSolutionAtBoundaries( iEl ) 
         CALL myDGSEM % CalculateVelocityAtBoundaries( iEl )

      ENDDO 
      CALL myDGSEM % clocks % StopThisTimer( tID )
      CALL myDGSEM % clocks % AccumulateTimings( )
      tID = tID + 1

!$OMP END DO 
!$OMP FLUSH( myDGSEM )

!$OMP DO

      CALL myDGSEM % clocks % StartThisTimer( tID )
      DO iFace = 1, myDGSEM % mesh % nFaces
         CALL myDGSEM % EdgeFlux( iFace, tn )
      ENDDO
      CALL myDGSEM % clocks % StopThisTimer( tID )
      CALL myDGSEM % clocks % AccumulateTimings( )
      tID = tID + 1
 
!$OMP END DO
!$OMP FLUSH( myDGSEM ) 
 

!$OMP DO
      CALL myDGSEM % clocks % StartThisTimer( tID )
      DO iEl = 1, myDGSEM % mesh % nElems
         CALL myDGSEM % MappedTimeDerivative( iEl, tn )
      ENDDO
      CALL myDGSEM % clocks % StopThisTimer( tID )
      CALL myDGSEM % clocks % AccumulateTimings( )
!$OMP END DO
!$OMP FLUSH( myDGSEM )
!$OMP END PARALLEL

 
 END SUBROUTINE GlobalTimeDerivative_Advection
!
!
! 
SUBROUTINE DoTheAdaptiveFiltering_Advection( myDGSEM, iEl, m )
 !
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS(Advection), INTENT(inout) :: myDGSEM
   INTEGER, INTENT(in)             :: iEl, m
   ! Local
   REAL(prec) :: sol(0:myDGSEM % nS, &
                     0:myDGSEM % nP, &
                     0:myDGSEM % nQ, &
                     1:myDGSEM % nEq)
   REAL(prec) :: solf(0:myDGSEM % nS, &
                      0:myDGSEM % nP, &
                      0:myDGSEM % nQ)
   REAL(prec) :: c1(0:myDGSEM % nS, &
                    0:myDGSEM % nP, &
                    0:myDGSEM % nQ)
   REAL(prec) :: c2(0:myDGSEM % nS, &
                    0:myDGSEM % nP, &
                    0:myDGSEM % nQ)
   REAL(prec) :: J(0:myDGSEM % nS, &
                   0:myDGSEM % nP, &
                   0:myDGSEM % nQ )
   REAL(prec) :: temp(0:myDGSEM % nS)
   REAL(prec) :: wS(0:myDGSEM % nS), wP(0:myDGSEM % nP), wQ(0:myDGSEM % nQ)
   REAL(prec) :: E1prior, E2prior, E1, E2, dE2, dE1, xi
   INTEGER    :: iS, iP, iQ, iEq


      CALL myDGSEM % GetSolution( iEl, sol )
      CALL myDGSEM % dgStorage % GetQuadratureWeights( wS, wP, wQ )
      CALL myDGSEM % mesh % GetJacobian( iEl, J )

      ! The first step in the process is to apply the filter to obtain the well resolved solution
      DO iEq = 1, myDGSEM % nEq
         solf = myDGSEM % modalFilter % ApplyFilter( sol(:,:,:,iEq) )
      
      ! This filtered solution is subtracted from the full solution to obtain the marginally-resolved
      ! portion of the solution
         c2 = (sol(:,:,:,iEq) - solf)**2 ! c2 is the energy in the marginally resolved modes
         c1 = solf**2         ! c1 is the energy in the well resolved modes

      
      ! Now that we have the two distinct components of the Legendre spectra, we want to calculate
      ! the energy in each component, and the change in the energy of each component from the 
      ! previous model state. 
      ! If the small scale (marginally resolved) exhibits a growth in energy, this should be balanced
      ! by a decay in the small scale energy. Aliasing errors may cause unphysical growth in the 
      ! energy associated with the marginally resolved. In this case, the solution is assigned to the
      ! filtered solution, effectively implying dissipation.
     
         E1prior = myDGSEM % E1(iEl,iEq)
         E2prior = myDGSEM % E2(iEl,iEq)

         E1 = ZERO
         E2 = ZERO

         ! For quicker volume integration the energies are pre-multiplied by the Jacobian.
         c1 = c1*J !
         c2 = c2*J
      ! Volume integration of the energy of the resolved and marginally resolved fields is done here
      
         DO iQ = 0, myDGSEM % nQ
            DO iP = 0, myDGSEM % nP
               temp = c1(:,iP,iQ)
               E1 = E1 + DOT_PRODUCT( temp, wS )*wP(iP)*wQ(iQ)
               temp = c2(:,iP,iQ)
               E2 = E2 + DOT_PRODUCT( temp, wS )*wP(iP)*wQ(iQ)
            ENDDO
         ENDDO

         myDGSEM % E1(iEl,iEq) = E1 
         myDGSEM % E2(iEl,iEq) = E2
      
         dE1 = E1-E1prior
         dE2 = E2-E2prior

         myDGSEM % dE1(iEl,iEq) = dE1
         myDGSEM % dE2(iEl,iEq) = dE2

         xi = (E2/E1)
  
         IF( dE2 > ZERO .AND. abs(dE1)/dE2 > ONE )THEN ! The energy in the small scales is growing faster than the large scale is giving it up
                            
            myDGSEM % lim(iEq,iEl,m) = xi
   
         ELSEIF( dE2 < ZERO )THEN
  
            myDGSEM % lim(iEq,iEl,m) = xi

         ELSE
         
            IF( m > 1)THEN
               myDGSEM % lim(iEq,iEl,m) = myDGSEM % lim(iEq,iEl,m-1)
            ENDIF

         ENDIF


         IF( xi > 1.05_prec*myDGSEM % lim(iEq,iEl,m) )THEN
            CALL myDGSEM % SetSolutionWithVarID( iEl, iEq, solf )
         !   PRINT*, 'Filtered!'
         ENDIF

      ENDDO

 END SUBROUTINE DoTheAdaptiveFiltering_Advection
!
!
!
 SUBROUTINE EdgeFlux_Advection( myDGSEM, iFace, tn )
 ! S/R EdgeFlux
 ! 
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS( Advection ), INTENT(inout) :: myDGSEM
   INTEGER, INTENT(in)               :: iFace  
   REAL(prec), INTENT(in)            :: tn
   ! Local
   INTEGER :: nEq
   INTEGER :: e1, s1, e2, s2, iError, nS, nP, nQ
   INTEGER :: iStart, jStart, iInc, jInc, swapDim
   INTEGER :: i, j, k, l
   REAL(prec) :: flux(1:myDGSEM % nEq)
   REAL(prec) :: inState(0:myDGSEM % nS,0:myDGSEM % nP, 1:myDGSEM % nEq)
   REAL(prec) :: exState(0:myDGSEM % nS,0:myDGSEM % nP, 1:myDGSEM % nEq)
   REAL(prec) :: inS(1:myDGSEM % nEq), outS(1:myDGSEM % nEq)
   REAL(prec) :: x, y, z
   REAL(prec) :: u(0:myDGSEM % nS, 0:myDGSEM % nP)
   REAL(prec) :: v(0:myDGSEM % nS, 0:myDGSEM % nP) 
   REAL(prec) :: w(0:myDGSEM % nS, 0:myDGSEM % nP) 
   REAL(prec) :: nHat(1:3), nHatLength
    
      nEQ = myDGSEM % params % nTracers

      CALL myDGSEM % dgStorage % GetNumberOfNodes( nS, nP, nQ )
      
      CALL myDGSEM % mesh % GetFacePrimaryElementID( iFace, e1 )
      CALL myDGSEM % mesh % GetFacePrimaryElementSide( iFace, s1 )
      CALL myDGSEM % mesh % GetFaceSecondaryElementID( iFace, e2 )
      CALL myDGSEM % mesh % GetFaceSecondaryElementSide( iFace, s2 )
      
      s2 = ABS(s2)

      CALL myDGSEM % mesh % GetFaceStart( iFace, iStart, jStart )
      CALL myDGSEM % mesh % GetFaceIncrement( iFace, iInc, jInc )
      CALL myDGSEM % mesh % GetSwapDimensions( iFace, swapDim )

      CALL myDGSEM % GetBoundarySolutionAtBoundary( e1, s1, inState )
      CALL myDGSEM % GetVelocityAtBoundary( e1, s1, u, v, w )

      IF( e2 > 0 )then ! this is an interior edge
      
         
         CALL myDGSEM % GetBoundarySolutionAtBoundary( e2, s2, exState )
        
     
         k = (iStart-iInc)*(1-swapDim) + (jStart-jInc)*(swapDim)
         l = (iStart-iInc)*(swapDim) + (jStart-jInc)*(1-swapDim)
         !print*, k,l
         DO j = 0, nS ! Loop over the nodes
         
            k = (k + (jInc)*(swapDim))*(swapDim) + (iStart-iInc)*(1-swapDim)
            l = (l + (jInc)*(1-swapDim))*(1-swapDim) + (iStart-iInc)*(swapDim)
            
            DO i = 0, nP
               k = k + (iInc)*(1-swapDim)
               l = l + (iInc)*(swapDim)
              
               CALL myDGSEM % mesh % GetBoundaryNormalAtNode( e1, nHat, nHatLength, i, j, s1 ) ! Get nHat

               inS  = inState(i,j,1:nEq)
               outS = exState(k,l,1:nEq)
               flux = RiemannSolver( inS, outS, &
                                     u(i,j), v(i,j), w(i,j), &
                                     nHat, nEq )*nHatLength

               ! Store the flux for the elements which share this edge
               CALL myDGSEM % SetBoundaryFluxAtBoundaryNode( e1, s1, i, j, flux )
               CALL myDGSEM % SetBoundaryFluxAtBoundaryNode( e2, s2, k, l, -flux )

               
            ENDDO ! i, 
         ENDDO ! j, loop over the nodes
         
      ELSE ! this edge is a boundary edge

         DO j = 0, nS ! loop over the nodes
            DO i = 0, nP
               
               CALL myDGSEM % mesh % GetBoundaryNormalAtNode( e1, nHat, nHatLength, i, j, s1 ) ! Get nHat
            
               ! Get the boundary point locations
               CALL myDGSEM % mesh % GetBoundaryLocationAtNode( e1, x, y, z, i, j, s1 )

               ! Calculate the external state
               inS = inState(i,j,1:nEq)
               exState(i,j,1:nEq) = GetExternalState( nEq, x, y, z, tn, e2, inS)
 
               ! Calculate the RIEMANN flux
               outS = exState(i,j,1:nEq)
               flux = RiemannSolver( inS, outS, &
                                     u(i,j), v(i,j), w(i,j), &
                                     nHat, nEq )*nHatLength
 
               CALL myDGSEM % SetBoundaryFluxAtBoundaryNode( e1, s1, i, j, flux )

            ENDDO 
         ENDDO ! j, loop over the nodes on this edge

      ENDIF ! Choosing the edge type
 
 END SUBROUTINE EdgeFlux_Advection
!
!
!
 SUBROUTINE MappedTimeDerivative_Advection( myDGSEM, iEl, tn ) 
 ! S/R MappedTimeDerivative_Advection
 ! 
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS(Advection), INTENT(inout) :: myDGSEM
   INTEGER, INTENT(in)             :: iEl
   REAL(prec), INTENT(in)          :: tn
   ! LOCAL
   REAL(prec) :: dxds, dxdp, dxdq, dyds, dydp, dydq, dzds, dzdp, dzdq
   REAL(prec) :: J, x, y, z
   REAL(prec) :: xF(1:myDGSEM % nEq), yF(1:myDGSEM % nEq), zF(1:myDGSEM % nEq), sol(1:myDGSEM % nEq)
   REAL(prec) :: relaxSol(1:myDGSEM % nEq), rFac(1:myDGSEM % nEq)
   REAL(prec) :: u(0:myDGSEM % nS,0:myDGSEM % nS,0:myDGSEM % nQ)
   REAL(prec) :: v(0:myDGSEM % nP,0:myDGSEM % nP,0:myDGSEM % nQ)
   REAL(prec) :: w(0:myDGSEM % nP,0:myDGSEM % nP,0:myDGSEM % nQ)
   REAL(prec) :: contFlux(0:myDGSEM % nS,1:myDGSEM % nEq)
   REAL(prec) :: contFluxDer(0:myDGSEM % nS,1:myDGSEM % nEq)
   REAL(prec) :: tend(0:myDGSEM % nS,0:myDGSEM % nP,0:myDGSEM % nQ,1:myDGSEM % nEq)
   REAL(prec) :: dMatX(0:myDGSEM % nS,0:myDGSEM % nS)
   REAL(prec) :: dMatY(0:myDGSEM % nP,0:myDGSEM % nP)
   REAL(prec) :: dMatZ(0:myDGSEM % nQ,0:myDGSEM % nQ)
   REAL(prec) :: qWeightX(0:myDGSEM % nS)
   REAL(prec) :: qWeightY(0:myDGSEM % nP)
   REAL(prec) :: qWeightZ(0:myDGSEM % nQ)
   REAL(prec) :: lsouth(0:myDGSEM % nP)
   REAL(prec) :: lnorth(0:myDGSEM % nP)
   REAL(prec) :: least(0:myDGSEM % nS)
   REAL(prec) :: lwest(0:myDGSEM % nS)
   REAL(prec) :: ltop(0:myDGSEM % nQ)
   REAL(prec) :: lbottom(0:myDGSEM % nQ)
   REAL(prec) :: fL(1:myDGSEM % nEq), fR(1:myDGSEM % nEq)
   REAL(prec) :: kappa
   INTEGER    :: iS, iP, iQ, nS, nP, nQ, iEq

     CALL myDGSEM % dgStorage % GetNumberOfNodes( nS, nP, nQ )
     CALL myDGSEM % dgStorage % GetDerivativeMatrix( dMatX, dMatY, dMatZ )
     CALL myDGSEM % dgStorage % GetQuadratureWeights( qWeightX, qWeightY, qWeightZ )
     CALL myDGSEM % dgStorage % GetSouthernInterpolants( lsouth )
     CALL myDGSEM % dgStorage % GetNorthernInterpolants( lnorth )
     CALL myDGSEM % dgStorage % GetEasternInterpolants( least )
     CALL myDGSEM % dgStorage % GetWesternInterpolants( lwest )
     CALL myDGSEM % dgStorage % GetTopInterpolants( ltop )
     CALL myDGSEM % dgStorage % GetBottomInterpolants( lbottom )
     CALL myDGSEM % GetVelocity( iEl, u, v, w )

      DO iQ = 0, nQ
         DO iP = 0,nP ! Loop over the y-points
            DO iS = 0,nS ! Loop over the x-points

               ! Get the metric terms to calculate the contravariant flux
               CALL myDGSEM % mesh % GetCovariantMetricsAtNode( iEl, &
                                                                dxds, dxdp, dxdq, &
                                                                dyds, dydp, dydq, &
                                                                dzds, dzdp, dzdq, &
                                                                iS, iP, iQ )

               ! Get the solution at x,y
               CALL myDGSEM % GetSolutionAtNode( iEl, iS, iP, iQ, sol )
               ! Calculate the x and y fluxes
               xF = XFlux(  tn, myDGSEM % nEq, u(iS,iP,iQ), v(iS,iP,iQ), w(iS,iP,iQ), sol )
               yF = YFlux(  tn, myDGSEM % nEq, u(iS,iP,iQ), v(iS,iP,iQ), w(iS,iP,iQ), sol )
               zF = ZFlux(  tn, myDGSEM % nEq, u(iS,iP,iQ), v(iS,iP,iQ), w(iS,iP,iQ), sol )

               !And now the contravariant flux
               contFlux(iS,:) = (dydp*dzdq - dzdp*dydq)*xF + &
                                 (dxdq*dzdp - dxdp*dzdq)*yF + &
                                 (dxdp*dydq - dxdq*dydp)*zF

            ENDDO ! iS, Loop over the x-points

            ! Get the numerical fluxes at the boundaries (should be calculated before
            ! this routine is called )
        
            ! Get the "west" flux at iP
            CALL myDGSEM % GetBoundaryFluxAtBoundaryNode( iEl, west, iP, iQ, fL )

            ! Get the "east" flux at iP
            CALL myDGSEM % GetBoundaryFluxAtBoundaryNode( iEl, east, iP, iQ, fR )

            ! At this y-level, calculate the DG-advective derivative
            DO iEq = 1, myDGSEM % nEq
               contFluxDer(0:nS,iEq) = DGSystemDerivative( nS, dMatX, qWeightX, fL(iEq), fR(iEq), &
                                                           contFlux(0:nS,iEq), lwest, least  )
            ENDDO
     
            DO iS = 0, nS ! Loop over the x-points
               DO iEq = 1, myDGSEM % nEq  ! Loop over the number of equations
                  tend(iS, iP, iQ, iEq) = -contFluxDer(iS,iEq)
               ENDDO ! Loop over the number of equations
            ENDDO ! Loop over the x-points

         ENDDO ! iP  


         DO iS = 0,nS ! Loop over the x-points
            DO iP = 0,nP ! Loop over the y-points

               ! Get the metric terms to calculate the contravariant flux
               CALL myDGSEM % mesh % GetCovariantMetricsAtNode( iEl, &
                                                                dxds, dxdp, dxdq, &
                                                                dyds, dydp, dydq, &
                                                                dzds, dzdp, dzdq, &
                                                                iS, iP, iQ )

               ! Get the solution at x,y
               CALL myDGSEM % GetSolutionAtNode( iEl, iS, iP, iQ, sol )
               ! Calculate the x and y fluxes
               xF = XFlux(  tn, myDGSEM % nEq, u(iS,iP,iQ), v(iS,iP,iQ), w(iS,iP,iQ), sol )
               yF = YFlux(  tn, myDGSEM % nEq, u(iS,iP,iQ), v(iS,iP,iQ), w(iS,iP,iQ), sol )
               zF = ZFlux(  tn, myDGSEM % nEq, u(iS,iP,iQ), v(iS,iP,iQ), w(iS,iP,iQ), sol )

               !And now the contravariant flux
               contFlux(iP,:) = (dydq*dzds - dzdq*dyds)*xF + &
                                 (dxds*dzdq - dxdq*dzds)*yF + &
                                 (dxdq*dyds - dxds*dydq)*zF

            ENDDO ! iP, Loop over the y-points

            ! Get the numerical fluxes at the boundaries (should be calculated before
            ! this routine is called )
        
            ! Get the "south" flux at iS
            CALL  myDGSEM % GetBoundaryFluxAtBoundaryNode( iEl, south, iS, iQ, fL )

            ! Get the "north" flux at iS
            CALL  myDGSEM % GetBoundaryFluxAtBoundaryNode( iEl, north, iS, iQ, fR )


            ! At this x-level, calculate the y-DG-advective derivative
            DO iEq = 1, myDGSEM % nEq
               contFluxDer(0:nP,iEq) = DGSystemDerivative(  nP, dMatY, qWeightY, fL(iEq), fR(iEq), &
                                                            contFlux(0:nP,iEq), lsouth, lnorth  )
            ENDDO

         
            DO iP = 0, nP 
               DO iEq = 1, myDGSEM % nEq  
                  tend(iS, iP, iQ, iEq) = tend(iS, iP, iQ, iEq)- contFluxDer(iP,iEq)
               ENDDO 
            ENDDO 

         ENDDO ! iS
      ENDDO ! iQ

      DO iP = 0 ,nP
         DO iS = 0, nS
            DO iQ = 0, nQ
               ! Get the metric terms to calculate the contravariant flux
               CALL myDGSEM % mesh % GetCovariantMetricsAtNode( iEl, &
                                                                dxds, dxdp, dxdq, &
                                                                dyds, dydp, dydq, &
                                                                dzds, dzdp, dzdq, &
                                                                iS, iP, iQ )

               ! Get the solution at x,y
               CALL myDGSEM % GetSolutionAtNode( iEl, iS, iP, iQ, sol )
               ! Calculate the x and y fluxes
               xF = XFlux(  tn, myDGSEM % nEq, u(iS,iP,iQ), v(iS,iP,iQ), w(iS,iP,iQ), sol )
               yF = YFlux(  tn, myDGSEM % nEq, u(iS,iP,iQ), v(iS,iP,iQ), w(iS,iP,iQ), sol )
               zF = ZFlux(  tn, myDGSEM % nEq, u(iS,iP,iQ), v(iS,iP,iQ), w(iS,iP,iQ), sol )

               !And now the contravariant flux
               contFlux(iQ,:) = (dyds*dzdp - dzds*dydp)*xF + &
                                 (dxdp*dzds - dxds*dzdp)*yF + &
                                 (dxds*dydp - dxdp*dyds)*zF

            ENDDO

            ! Get the numerical fluxes at the boundaries (should be calculated before
            ! this routine is called )
        
            ! Get the bottom flux
            CALL  myDGSEM % GetBoundaryFluxAtBoundaryNode( iEl, bottom, iS, iP, fL )

            ! Get the top flux
            CALL  myDGSEM % GetBoundaryFluxAtBoundaryNode( iEl, top, iS, iP, fR )


            ! At this x-level, calculate the y-DG-advective derivative
            DO iEq = 1, myDGSEM % nEq
               contFluxDer(0:nQ,iEq) = DGSystemDerivative(  nQ, dMatZ, qWeightY, fL(iEq), fR(iEq), &
                                                            contFlux(0:nP,iEq), lbottom, ltop  )
            ENDDO

         
            DO iQ = 0, nQ 
               CALL myDGSEM % mesh % GetJacobianAtNode( iEl, J, iS, iP, iQ )
               DO iEq = 1, myDGSEM % nEq  
                  tend(iS, iP, iQ, iEq) = ( tend(iS, iP, iQ, iEq)- contFluxDer(iQ,iEq) )/J
               ENDDO 
            ENDDO 

         ENDDO
      ENDDO
         
      ! COMPUTE THE SOURCE TERMS
      DO iQ = 0, nQ  
         DO iP = 0, nP
            DO iS = 0, nS

               CALL myDGSEM % GetSolutionAtNode( iEl, iS, iP, iQ, sol )
               CALL myDGSEM % GetRelaxationFieldAtNode( iEl, iS, iP, iQ, relaxSol )
               CALL myDGSEM % GetRelaxationFactorAtNode( iEl, iS, iP, iQ, rFac )
               CALL myDGSEM % mesh % GetPositionAtNode( iEl, x, y, z, iS, iP, iQ )
         
                  
               tend(iS,iP,iQ,:) = tend(iS,iP,iQ,:) + Source( tn, sol, relaxSol, rFac, &
                                                             x, y, z, myDGSEM % nEq  )
            ENDDO
         ENDDO 
      ENDDO 

     CALL myDGSEM % SetTendency( iEl, tend )

 END SUBROUTINE MappedTimeDerivative_Advection
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
 FUNCTION RiemannSolver( inState, outState, u, v, w, nHat, nEq ) RESULT( numFlux )
 ! FUNCTION RiemannSolver 
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   INTEGER    :: nEq
   REAL(prec) :: inState(1:nEq)
   REAL(prec) :: outState(1:nEq)
   REAL(prec) :: u, v, w
   REAL(prec) :: nHat(1:nDim) ! normal direction
   REAL(prec) :: numFlux(1:nEq)
   
   ! LOCAL
   REAL(prec) :: uN
   REAL(prec) :: jump(1:nEq), aS(1:nEq)


      jump = outState - inState
      uN =  u*nHat(1) + v*nHat(2) + w*nHat(3)

      numFlux = HALF*( uN*( outState + inState ) - abs(uN)*jump )



 END FUNCTION RiemannSolver
!
!
! 
 FUNCTION XFlux(  tn, nEq, u, v, w, solAtX ) RESULT( fx )
 ! FUNCTION XFlux
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE 
   INTEGER    :: nEq
   REAL(prec) :: tn
   REAL(prec) :: u, v, w
   REAL(prec) :: solAtX(1:nEq)
   REAL(prec) :: fx(1:nEq)

      fx = solAtX*u 

 END FUNCTION XFlux
!
!
!
 FUNCTION YFlux(  tn, nEq, u, v, w, solAtX ) RESULT( fy )
 ! FUNCTION YFlux
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE 
   INTEGER    :: nEq
   REAL(prec) :: tn
   REAL(prec) :: u, v, w
   REAL(prec) :: solAtX(1:nEq)
   REAL(prec) :: fy(1:nEq)

      fy = solAtX*v

 END FUNCTION YFlux
!
!
!
FUNCTION ZFlux(  tn, nEq, u, v, w, solAtX ) RESULT( fz )
 ! FUNCTION ZFlux
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE 
   INTEGER    :: nEq
   REAL(prec) :: tn
   REAL(prec) :: u, v, w
   REAL(prec) :: solAtX(1:nEq)
   REAL(prec) :: fz(1:nEq)

      fz = solAtX*w

 END FUNCTION ZFlux
!
!
!           
 FUNCTION Source( tn, sol, relax, rfac, x, y, z, nEq ) RESULT( q )
 ! FUNCTION Source
 ! 
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   INTEGER        :: nEq
   REAL(prec)     :: tn
   REAL(prec)     :: sol(1:nEq), relax(1:nEq), rFac(1:nEq)
   REAL(prec)     :: x, y, z
   REAL(prec)     :: q(1:nEq)
   ! LOCAL
   INTEGER :: iEq

      DO iEq = 1, nEq
         q(iEq) = (relax(iEq) - sol(iEq))*rFac(iEq)
      ENDDO

 END FUNCTION Source
!
!
!                         
 FUNCTION GetExternalState( nEq, x, y, z, t, bcFlag, intState ) RESULT( extState )
 ! S/R GetExternalState
 ! 
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   INTEGER    :: nEq
   REAL(prec) :: intState(1:nEq)
   REAL(prec) :: x, y, z,t
   REAL(prec) :: nHat(1:2)
   INTEGER    :: bcFlag, formulation
   REAL(prec) :: extState(1:nEq)

       extState = ZERO 

       IF( bcFlag == InflowOne  )THEN
          
          extState(1) = tanh(t) 
             
       ELSEIF( bcFlag == InflowTwo )THEN

          extState(2) =  tanh(t)
          
       ELSEIF( bcFlag == NO_NORMAL_FLOW )THEN
      
          extState = exp( -( (x-t)**2 + (y-t)**2 + (z-t)**2)/(0.2_prec**2) ) 
  
       ELSEIF( bcFlag == RADIATION )THEN

          extState = ZERO

       ELSE      
       
         PRINT*, 'FUNCTION : GetExternalState : Invalid bcflag =', bcflag
         extState = ZERO
         
       ENDIF

        
 END FUNCTION GetExternalState
!
!
!==================================================================================================!
!-------------------------------------- FILE I/O ROUTINES -----------------------------------------!
!==================================================================================================!
!
!
 SUBROUTINE CoarseToFine_Advection( myDGSEM, iEl, x, y, z, u, v, w, c  )
 ! CoarseToFine
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS( Advection ), INTENT(in) :: myDGSEM
   INTEGER, INTENT(in)            :: iEl
   REAL(prec), INTENT(out)        :: x(0:myDGSEM % nPlot, 0:myDGSEM % nPlot, 0:myDGSEM % nPlot)
   REAL(prec), INTENT(out)        :: y(0:myDGSEM % nPlot, 0:myDGSEM % nPlot, 0:myDGSEM % nPlot)
   REAL(prec), INTENT(out)        :: z(0:myDGSEM % nPlot, 0:myDGSEM % nPlot, 0:myDGSEM % nPlot)
   REAL(prec), INTENT(out)        :: u(0:myDGSEM % nPlot, 0:myDGSEM % nPlot, 0:myDGSEM % nPlot)
   REAL(prec), INTENT(out)        :: v(0:myDGSEM % nPlot, 0:myDGSEM % nPlot, 0:myDGSEM % nPlot)
   REAL(prec), INTENT(out)        :: w(0:myDGSEM % nPlot, 0:myDGSEM % nPlot, 0:myDGSEM % nPlot)
   REAL(prec), INTENT(out)        :: c(0:myDGSEM % nPlot, 0:myDGSEM % nPlot, 0:myDGSEM % nPlot,1:myDGSEM % nEq)
   ! Local
   REAL(prec) :: localX(0:myDGSEM % nS, 0:myDGSEM % nP, 0:myDGSEM % nQ)
   REAL(prec) :: localY(0:myDGSEM % nS, 0:myDGSEM % nP, 0:myDGSEM % nQ)
   REAL(prec) :: localZ(0:myDGSEM % nS, 0:myDGSEM % nP, 0:myDGSEM % nQ)
   REAL(prec) :: localU(0:myDGSEM % nS, 0:myDGSEM % nP, 0:myDGSEM % nQ)
   REAL(prec) :: localV(0:myDGSEM % nS, 0:myDGSEM % nP, 0:myDGSEM % nQ)
   REAL(prec) :: localW(0:myDGSEM % nS, 0:myDGSEM % nP, 0:myDGSEM % nQ)
   REAL(prec) :: sol(0:myDGSEM % nS, 0:myDGSEM % nP, 0:myDGSEM % nQ, 1:myDGSEM % nEq)
   INTEGER    :: iEq


      CALL myDGSEM % GetSolution( iEl, sol )
      CALL myDGSEM % GetVelocity( iEl, localU, localV, localW )
      CALL myDGSEM % mesh % GetPositions( iEl, localX, localY, localZ )
      
      CALL myDGSEM % dgStorage % interp % CoarseToFine( localX, &
                                                        myDGSEM % plMatS, &
                                                        myDGSEM % plMatP, &
                                                        myDGSEM % plMatQ, &
                                                        myDGSEM % nPlot, &
                                                        myDGSEM % nPlot, &
                                                        myDGSEM % nPlot, &
                                                        x )
      
      CALL myDGSEM % dgStorage % interp % CoarseToFine( localY, &
                                                        myDGSEM % plMatS, &
                                                        myDGSEM % plMatP, &
                                                        myDGSEM % plMatQ, &
                                                        myDGSEM % nPlot, &
                                                        myDGSEM % nPlot, &
                                                        myDGSEM % nPlot, &
                                                        y )
                                                    
      CALL myDGSEM % dgStorage % interp % CoarseToFine( localZ, &
                                                        myDGSEM % plMatS, &
                                                        myDGSEM % plMatP, &
                                                        myDGSEM % plMatQ, &
                                                        myDGSEM % nPlot, &
                                                        myDGSEM % nPlot, &
                                                        myDGSEM % nPlot, &
                                                        z )

      CALL myDGSEM % dgStorage % interp % CoarseToFine( localU, &
                                                        myDGSEM % plMatS, &
                                                        myDGSEM % plMatP, &
                                                        myDGSEM % plMatQ, &
                                                        myDGSEM % nPlot, &
                                                        myDGSEM % nPlot, &
                                                        myDGSEM % nPlot, &
                                                        u )
                                                        
      CALL myDGSEM % dgStorage % interp % CoarseToFine( localV, &
                                                        myDGSEM % plMatS, &
                                                        myDGSEM % plMatP, &
                                                        myDGSEM % plMatQ, &
                                                        myDGSEM % nPlot, &
                                                        myDGSEM % nPlot, &
                                                        myDGSEM % nPlot, &
                                                        v )

      CALL myDGSEM % dgStorage % interp % CoarseToFine( localW, &
                                                        myDGSEM % plMatS, &
                                                        myDGSEM % plMatP, &
                                                        myDGSEM % plMatQ, &
                                                        myDGSEM % nPlot, &
                                                        myDGSEM % nPlot, &
                                                        myDGSEM % nPlot, &
                                                        w )
               
      DO iEq = 1, myDGSEM % nEq                                         
         CALL myDGSEM % dgStorage % interp % CoarseToFine( sol(:,:,:,iEq), &
                                                           myDGSEM % plMatS, &
                                                           myDGSEM % plMatP, &
                                                           myDGSEM % plMatQ, &
                                                           myDGSEM % nPlot, &
                                                           myDGSEM % nPlot, &
                                                           myDGSEM % nPlot, &
                                                           c(:,:,:,iEq) )  
      ENDDO                               
      
 END SUBROUTINE CoarseToFine_Advection
!
!
!
 SUBROUTINE WriteTecplot_Advection( myDGSEM, filename )
 ! S/R WriteTecplot
 !  
 ! =============================================================================================== !
 ! DECLARATIONS
  IMPLICIT NONE
  CLASS( Advection ), INTENT(in)  :: myDGsem
  CHARACTER(*), INTENT(in), OPTIONAL :: filename
  !LOCAL
  REAL(prec)  :: x(0:myDGSEM % nPlot,0:myDGSEM % nPlot,0:myDGSEM % nPlot)
  REAL(prec)  :: y(0:myDGSEM % nPlot,0:myDGSEM % nPlot,0:myDGSEM % nPlot)
  REAL(prec)  :: z(0:myDGSEM % nPlot,0:myDGSEM % nPlot,0:myDGSEM % nPlot)
  REAL(prec)  :: u(0:myDGSEM % nPlot,0:myDGSEM % nPlot,0:myDGSEM % nPlot)
  REAL(prec)  :: v(0:myDGSEM % nPlot,0:myDGSEM % nPlot,0:myDGSEM % nPlot)
  REAL(prec)  :: w(0:myDGSEM % nPlot,0:myDGSEM % nPlot,0:myDGSEM % nPlot)
  REAL(prec)  :: c(0:myDGSEM % nPlot,0:myDGSEM % nPlot,0:myDGSEM % nPlot,1:myDGSEM % nEq)
  INTEGER     :: iS, iP, iQ, iEl, fUnit, nPlot
  CHARACTER(len=5) :: zoneID

    nPlot = myDGSEM % nPlot
    
    IF( PRESENT(filename) )THEN
       OPEN( UNIT=NEWUNIT(fUnit), &
             FILE= TRIM(filename)//'.tec', &
             FORM='formatted', &
             STATUS='replace')
    ELSE
       OPEN( UNIT=NEWUNIT(fUnit), &
             FILE= 'Advection.tec', &
             FORM='formatted', &
             STATUS='replace')  
    ENDIF
    
    WRITE(fUnit,*) 'VARIABLES = "X", "Y","Z", "U", "V", "W", "C1", "C2" '
 
    DO iEl = 1, myDGsem % mesh % nElems

        CALL myDGSEM % CoarseToFine( iEl, x, y, z, u, v, w, c )
        WRITE(zoneID,'(I5.5)') iEl
        WRITE(fUnit,*) 'ZONE T="el'//trim(zoneID)//'", I=',nPlot+1,', J=', nPlot+1,', K=', nPlot+1,',F=POINT'

        DO iQ = 0, nPlot
           DO iP = 0, nPlot
              DO iS = 0, nPlot
                 WRITE (fUnit,*)  x(iS,iP,iQ), y(iS,iP,iQ), z(iS,iP,iQ), &
                                  u(iS,iP,iQ), v(iS,iP,iQ), w(iS,iP,iQ), &
                                  c(iS,iP,iQ,1), c(iS,iP,iQ,2)
              ENDDO
           ENDDO
        ENDDO
        
    ENDDO

    CLOSE(UNIT=fUnit)

 END SUBROUTINE WriteTecplot_Advection
!
!
!
 SUBROUTINE WritePickup_Advection( myDGSEM, iter )
 ! S/R WritePickup
 ! 
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS( Advection ), INTENT(in) :: myDGSEM
   INTEGER, INTENT(in)               :: iter
  ! LOCAL
   REAL(prec)    :: sol(0:myDGSEM % nS, 0:myDGSEM % nP, 0:myDGSEM % nQ, 1:myDGSEM % nEq)
   REAL(prec)    :: u(0:myDGSEM % nS, 0:myDGSEM % nP, 0:myDGSEM % nQ)
   REAL(prec)    :: v(0:myDGSEM % nS, 0:myDGSEM % nP, 0:myDGSEM % nQ)
   REAL(prec)    :: w(0:myDGSEM % nS, 0:myDGSEM % nP, 0:myDGSEM % nQ)
   CHARACTER(10) :: iterChar
   INTEGER       :: iEl
   INTEGER       :: thisRec, fUnit
   INTEGER       :: iS, iP, iEq, nS, nP, nQ

     nS = myDGSEM % nS
     nP = myDGSEM % nP
     nQ = myDGSEM % nQ
     
     WRITE(iterChar,'(I10.10)') iter

     OPEN( UNIT=NEWUNIT(fUnit), &
           FILE='Advection.'//iterChar//'.pickup', &
           FORM='unformatted',&
           ACCESS='direct',&
           STATUS='replace',&
           ACTION='WRITE',&
           CONVERT='big_endian',&
           RECL=prec*(nS+1)*(nP+1)*(nQ+1) )

     thisRec = 1 
     DO iEl = 1, myDGSEM % mesh % nElems
        
        CALL myDGSEM % GetSolution( iEl, sol )
        DO iEq = 1, myDGSEM % nEq
           WRITE( fUnit, REC=thisRec )sol(:,:,:,iEq) 
           thisRec = thisRec+1
        ENDDO

        CALL myDGSEM % GetRelaxationField( iEl, sol )
        DO iEq = 1, myDGSEM % nEq
           WRITE( fUnit, REC=thisRec )sol(:,:,:,iEq) 
           thisRec = thisRec+1
        ENDDO

        CALL myDGSEM % GetRelaxationFactor( iEl, sol )
        DO iEq = 1, myDGSEM % nEq
           WRITE( fUnit, REC=thisRec )sol(:,:,:,iEq) 
           thisRec = thisRec+1
        ENDDO

        CALL myDGSEM % GetVelocity( iEl, u, v, w )
        WRITE( fUnit, REC=thisRec )u 
        thisRec = thisRec+1
        WRITE( fUnit, REC=thisRec )v 
        thisRec = thisRec+1
        WRITE( fUnit, REC=thisRec )w 
        thisRec = thisRec+1


     ENDDO

     CLOSE(UNIT=fUnit)

 END SUBROUTINE WritePickup_Advection
!
!
!
  SUBROUTINE ReadPickup_Advection( myDGSEM, iter )
 ! S/R ReadPickup
 ! 
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS( Advection ), INTENT(inout) :: myDGSEM
   INTEGER, INTENT(in)                  :: iter
  ! LOCAL
   REAL(prec)    :: sol(0:myDGSEM % nS, 0:myDGSEM % nP, 0:myDGSEM % nQ, 1:myDGSEM % nEq)
   REAL(prec)    :: u(0:myDGSEM % nS, 0:myDGSEM % nP, 0:myDGSEM % nQ)
   REAL(prec)    :: v(0:myDGSEM % nS, 0:myDGSEM % nP, 0:myDGSEM % nQ)
   REAL(prec)    :: w(0:myDGSEM % nS, 0:myDGSEM % nP, 0:myDGSEM % nQ)
   CHARACTER(10) :: iterChar
   INTEGER       :: iEl
   INTEGER       :: thisRec, fUnit
   INTEGER       :: iS, iP, iEq, nS, nP, nQ

     nS = myDGSEM % nS
     nP = myDGSEM % nP
     nQ = myDGSEM % nQ
     
     WRITE(iterChar,'(I10.10)') iter

     OPEN( UNIT=NEWUNIT(fUnit), &
           FILE='Advection.'//iterChar//'.pickup', &
           FORM='unformatted',&
           ACCESS='direct',&
           STATUS='old',&
           ACTION='READ',&
           CONVERT='big_endian',&
           RECL=prec*(nS+1)*(nP+1)*(nQ+1) )
     
     thisRec = 1
     DO iEl = 1, myDGSEM % mesh % nElems
        
        DO iEq = 1, myDGSEM % nEq
           READ( fUnit, REC=thisRec )sol(:,:,:,iEq) 
           thisRec = thisRec+1
        ENDDO
        CALL myDGSEM % SetSolution( iEl, sol )
        
        DO iEq = 1, myDGSEM % nEq
           READ( fUnit, REC=thisRec )sol(:,:,:,iEq) 
           thisRec = thisRec+1
        ENDDO
        CALL myDGSEM % SetRelaxationField( iEl, sol )

        
        DO iEq = 1, myDGSEM % nEq
           READ( fUnit, REC=thisRec )sol(:,:,:,iEq) 
           thisRec = thisRec+1
        ENDDO
        CALL myDGSEM % SetRelaxationFactor( iEl, sol )

        READ( fUnit, REC=thisRec )u 
        thisRec = thisRec+1
        READ( fUnit, REC=thisRec )v 
        thisRec = thisRec+1
        READ( fUnit, REC=thisRec )w 
        thisRec = thisRec+1
        CALL myDGSEM % SetVelocity( iEl, u, v, w )
        
     ENDDO

     CLOSE(UNIT=fUnit)

     
 END SUBROUTINE ReadPickup_Advection
!
!
! 
 END MODULE Advection3DClass



