! AdvectionClass.f90 ( new with v2.1 - 27 Jan. 2016)
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


MODULE AdvectionClass
! AdvectionClass.f90
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
! src/highend/advection2d/
USE AdvectionParamsClass
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
!     The Advection is based on the DGSEM_2D framework with the addition of 
!     the AddOns data. Also, a few additional routines are added to facilitate
!     the calculation of "source" terms which include the coriolis force, relaxation terms, etc.
! ================================================================================================ !  

    TYPE Advection
      INTEGER                               :: nEq, nPlot, nS, nP
      TYPE( QuadMesh )                      :: mesh
      TYPE( NodalStorage_2D )               :: dGStorage
      TYPE( NodalStorage_2D )               :: cGStorage
      TYPE( DGSEMSolution_2D ), ALLOCATABLE :: sol(:)
      TYPE( DGSEMSolution_2D ), ALLOCATABLE :: dCdx(:), dCdy(:) !, kappa(:)
      TYPE( DGSEMSolution_2D ), ALLOCATABLE :: velocity(:)
      TYPE( DGSEMSolution_2D ), ALLOCATABLE :: relax(:)
      TYPE( DGSEMSolution_2D ), ALLOCATABLE :: relaxFactor(:)
      TYPE( AdvectionParams )               :: params
      REAL(prec), ALLOCATABLE               :: plMatS(:,:), plMatP(:,:)

      CONTAINS

      PROCEDURE :: Build => Build_Advection
      PROCEDURE :: Trash => Trash_Advection
      PROCEDURE :: BuildQuadMesh => BuildQuadMesh_Advection
      
      ! DGSEMSolution_2DStorage_Class Wrapper Routines
      PROCEDURE :: GetSolution => GetSolution_Advection
      PROCEDURE :: SetSolution => SetSolution_Advection
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
      PROCEDURE :: GetGradient => GetGradient_Advection
      PROCEDURE :: SetGradient => SetGradient_Advection
      PROCEDURE :: GetGradientAtNode => GetGradientAtNode_Advection
      PROCEDURE :: SetGradientAtNode => SetGradientAtNode_Advection
      PROCEDURE :: GetGradientWithVarID => GetGradientWithVarID_Advection
      PROCEDURE :: SetGradientWithVarID => SetGradientWithVarID_Advection  
      PROCEDURE :: GetBoundarySolutionAtBoundary => GetBoundarySolutionAtBoundary_Advection
      PROCEDURE :: SetBoundarySolutionAtBoundary => SetBoundarySolutionAtBoundary_Advection
      PROCEDURE :: GetBoundaryGradientAtBoundary => GetBoundaryGradientAtBoundary_Advection
      PROCEDURE :: SetBoundaryGradientAtBoundary => SetBoundaryGradientAtBoundary_Advection
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
      PROCEDURE :: ForwardStepRK3 => ForwardStepRK3_Advection
      PROCEDURE :: EdgeFlux => EdgeFlux_Advection
      PROCEDURE :: EdgeGradientFlux => EdgeGradientFlux_Advection
      PROCEDURE :: MappedTimeDerivative => MappedTimeDerivative_Advection 
      PROCEDURE :: CalculateSolutionAtBoundaries => CalculateSolutionAtBoundaries_Advection
      PROCEDURE :: CalculateGradientAtBoundaries => CalculateGradientAtBoundaries_Advection
      PROCEDURE :: CalculateVelocityAtBoundaries => CalculateVelocityAtBoundaries_Advection
      PROCEDURE :: CalculateConcentrationGradient => CalculateConcentrationGradient_Advection

      
      PROCEDURE :: CoarseToFine => CoarseToFine_Advection
      PROCEDURE :: WriteTecplot => WriteTecplot_Advection
      PROCEDURE :: WritePickup => WritePickup_Advection
      PROCEDURE :: ReadPickup => ReadPickup_Advection

    END TYPE Advection



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
   INTEGER :: iEl, nS, nP, nPlot
   REAL(prec), ALLOCATABLE ::sNew(:)

      CALL myDGSEM % params % Build( )
     
      
      nS = myDGSEM % params % polyDeg
      nP = nS
      myDGSEM % nS      = nS
      myDGSEM % nP      = nP
      myDGSEM % nEq     = myDGSEM % params % nTracers
      myDGSEM % nPlot   = myDGSEM % params % nPlot
      
      nPlot = myDGSEM % params % nPlot
      
      ALLOCATE( sNew(0:nPlot) )
      CALL myDGSEM % dGStorage % Build( nS, nP, GAUSS, DG )
      CALL myDGSEM % cGStorage % Build( nS, nP, GAUSS, CG )
 
      CALL myDGSEM % BuildQuadMesh( )
      
      ! Set up the solution, relaxation fields, bathymetry, and vorticity
      ALLOCATE( myDGSEM % sol(1:myDGSEM % mesh % nElems) )
      ALLOCATE( myDGSEM % dCdx(1:myDGSEM % mesh % nElems) )
      ALLOCATE( myDGSEM % dCdy(1:myDGSEM % mesh % nElems) )
     ! ALLOCATE( myDGSEM % kappa(1:myDGSEM % mesh % nElems) )
      ALLOCATE( myDGSEM % relax(1:myDGSEM % mesh % nElems) )
      ALLOCATE( myDGSEM % relaxFactor(1:myDGSEM % mesh % nElems) )
      ALLOCATE( myDGSEM % velocity(1:myDGSEM % mesh % nElems) )

      ! Build and initialize the solution and the relaxation fields to zero
      DO iEl = 1, myDGSEM % mesh % nElems
         CALL myDGSEM % sol(iEl) % Build( nS, nP, myDGSEM % nEq )
         CALL myDGSEM % dCdx(iEl) % Build( nS, nP, myDGSEM % nEq )
         CALL myDGSEM % dCdy(iEl) % Build( nS, nP, myDGSEM % nEq )
         !CALL myDGSEM % kappa(iEl) % Build( nS, nP, myDGSEM % nEq )
         CALL myDGSEM % relax(iEl) % Build( nS, nP, myDGSEM % nEq )
         CALL myDGSEM % relaxFactor(iEl) % Build( nS, nP, myDGSEM % nEq ) 
         CALL myDGSEM % velocity(iEl) % Build( nS, nP, myDGSEM % nEq )
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

     DO iEl = 1, myDGSEM % mesh % nElems
        CALL myDGSEM % sol(iEl) % Trash( )
        CALL myDGSEM % dCdx(iEl) % Trash( )
        CALL myDGSEM % dCdy(iEl) % Trash( )
        CALL myDGSEM % relax(iEl) % Trash( )
        CALL myDGSEM % relaxFactor(iEl) % Trash( )
        CALL myDGSEM % velocity(iEl) % Trash( )
     ENDDO

     CALL myDGSEM % dGStorage % Trash( )
     CALL myDGSEM % cGStorage % Trash( )
     CALL myDGSEM % mesh % Trash( )

     DEALLOCATE( myDGSEM % sol ) 
     DEALLOCATE( myDGSEM % dCdx ) 
     DEALLOCATE( myDGSEM % dCdy ) 
     DEALLOCATE( myDGSEM % relax )
     DEALLOCATE( myDGSEM % relaxFactor )
     DEALLOCATE( myDGSEM % velocity )
     DEALLOCATE( myDGSEM % plMatS, myDGSEM % plMatP )

 END SUBROUTINE Trash_Advection
!
!
!
 SUBROUTINE BuildQuadMesh_Advection( myDGSEM )
 ! S/R BuildQuadMesh
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS( Advection ), INTENT(inout) :: myDGSEM
   ! LOCAL
   INTEGER :: iEdge, s1, e2
   
      PRINT*,'Module AdvectionClass.f90 : S/R BuildQuadMesh :'
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

      

 END SUBROUTINE BuildQuadMesh_Advection
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
   REAL(prec), INTENT(out)      :: theSolution(0:myDGSEM % nS, 0:myDGSEM % nP, 1:myDGSEM % nEq)

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
   REAL(prec), INTENT(in)          :: theSolution(0:myDGSEM % nS, 0:myDGSEM % nP, 1:myDGSEM % nEq)

      CALL myDGSEM % sol(iEl) % SetSolution( theSolution )

 END SUBROUTINE SetSolution_Advection
!
!
!
 SUBROUTINE GetSolutionAtNode_Advection( myDGSEM, iEl, i, j, theSolution  )
 ! S/R GetSolution
 !  
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS(Advection), INTENT(in) :: myDGSEM
   INTEGER, INTENT(in)          :: iEl
   INTEGER, INTENT(in)          :: i, j
   REAL(prec), INTENT(out)      :: theSolution(1:myDGSEM % nEq)

      CALL myDGSEM % sol(iEl) % GetSolutionAtNode( i, j, theSolution )

 END SUBROUTINE GetSolutionAtNode_Advection
!
!
!
 SUBROUTINE SetSolutionAtNode_Advection( myDGSEM, iEl, i, j, theSolution  )
 ! S/R SetSolution
 !  
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS(Advection), INTENT(inout) :: myDGSEM
   INTEGER, INTENT(in)             :: iEl
   INTEGER, INTENT(in)             :: i, j
   REAL(prec), INTENT(in)          :: theSolution(1:myDGSEM % nEq)

      CALL myDGSEM % sol(iEl) % SetSolutionAtNode( i, j, theSolution )

 END SUBROUTINE SetSolutionAtNode_Advection
!
!
!
 SUBROUTINE GetVelocity_Advection( myDGSEM, iEl, u, v  )
 ! S/R GetVelocity
 !  
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS(Advection), INTENT(in) :: myDGSEM
   INTEGER, INTENT(in)          :: iEl
   REAL(prec), INTENT(out)      :: u(0:myDGSEM % nS,0:myDGSEM % nP)
   REAL(prec), INTENT(out)      :: v(0:myDGSEM % nS,0:myDGSEM % nP)
   
      CALL myDGSEM % velocity(iEl) % GetSolutionWithVarID( 1, u )
      CALL myDGSEM % velocity(iEl) % GetSolutionWithVarID( 2, v )
      
 END SUBROUTINE GetVelocity_Advection
!
!
!
 SUBROUTINE SetVelocity_Advection( myDGSEM, iEl, u, v  )
 ! S/R SetVelocity
 !  
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS(Advection), INTENT(inout) :: myDGSEM
   INTEGER, INTENT(in)             :: iEl
   REAL(prec), INTENT(in)          :: u(0:myDGSEM % nS,0:myDGSEM % nP)
   REAL(prec), INTENT(in)          :: v(0:myDGSEM % nS,0:myDGSEM % nP)

      CALL myDGSEM % velocity(iEl) % SetSolutionWithVarID( 1, u )
      CALL myDGSEM % velocity(iEl) % SetSolutionWithVarID( 2, v )

 END SUBROUTINE SetVelocity_Advection
!
!
!
 SUBROUTINE GetVelocityAtBoundary_Advection( myDGSEM, iEl, boundary, u, v  )
 ! S/R GetBoundarySolution
 !  
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS(Advection), INTENT(in) :: myDGSEM
   INTEGER, INTENT(in)          :: iEl
   INTEGER, INTENT(in)          :: boundary
   REAL(prec), INTENT(out)      :: u(0:myDGSEM % nS)
   REAL(prec), INTENT(out)      :: v(0:myDGSEM % nS)

   
      CALL myDGSEM % velocity(iEl) % GetBoundarySolutionAtBoundaryWithVarID( boundary, 1, u )
      CALL myDGSEM % velocity(iEl) % GetBoundarySolutionAtBoundaryWithVarID( boundary, 2, v )

      
 END SUBROUTINE GetVelocityAtBoundary_Advection
!
!
!
 SUBROUTINE SetVelocityAtBoundary_Advection( myDGSEM, iEl, boundary, u, v  )
 ! S/R SetVelocityAtBoundary
 !  
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS(Advection), INTENT(inout) :: myDGSEM
   INTEGER, INTENT(in)             :: iEl
   INTEGER, INTENT(in)             :: boundary
   REAL(prec), INTENT(in)          :: u(0:myDGSEM % nS)
   REAL(prec), INTENT(in)          :: v(0:myDGSEM % nS)
   
      CALL myDGSEM % velocity(iEl) % SetBoundarySolutionAtBoundaryWithVarID( boundary, 1, u )
      CALL myDGSEM % velocity(iEl) % SetBoundarySolutionAtBoundaryWithVarID( boundary, 2, v )
     
 END SUBROUTINE SetVelocityAtBoundary_Advection
!
!
!
 SUBROUTINE GetSolutionAtNodeWithVarID_Advection( myDGSEM, iEl, i, j, varID, theSolution  )
 ! S/R GetSolution
 !  
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS(Advection), INTENT(in) :: myDGSEM
   INTEGER, INTENT(in)          :: iEl
   INTEGER, INTENT(in)          :: i, j, varID 
   REAL(prec), INTENT(out)      :: theSolution

      CALL myDGSEM % sol(iEl) % GetSolutionAtNodeWithVarID( i, j, varID, theSolution )

 END SUBROUTINE GetSolutionAtNodeWithVarID_Advection
!
!
!
 SUBROUTINE SetSolutionAtNodeWithVarID_Advection( myDGSEM, iEl, i, j, varID, theSolution  )
 ! S/R SetSolution
 !  
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS(Advection), INTENT(inout) :: myDGSEM
   INTEGER, INTENT(in)             :: iEl
   INTEGER, INTENT(in)             :: i, j, varID 
   REAL(prec), INTENT(in)          :: theSolution

      CALL myDGSEM % sol(iEl) % SetSolutionAtNodeWithVarID( i, j, varID, theSolution )

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
   REAL(prec), INTENT(out)      :: theTend(0:myDGSEM % nS,0:myDGSEM % nP, 1:myDGSEM % nEq)

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
   REAL(prec), INTENT(in)          :: theTend(0:myDGSEM % nS,0:myDGSEM % nP, 1:myDGSEM % nEq)

      CALL myDGSEM % sol(iEl) % SetTendency( theTend )   

 END SUBROUTINE SetTendency_Advection
!
!
!
 SUBROUTINE GetTendencyAtNode_Advection( myDGSEM, iEl, i, j, theTend  )
 ! S/R GetTendencyAtNode
 !  
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS(Advection), INTENT(in) :: myDGSEM
   INTEGER, INTENT(in)          :: iEl
   INTEGER, INTENT(in)          :: i, j 
   REAL(prec), INTENT(out)      :: theTend(1:myDGSEM % nEq)

      CALL myDGSEM % sol(iEl) % GetTendencyAtNode( i, j, theTend )

 END SUBROUTINE GetTendencyAtNode_Advection
!
!
!
 SUBROUTINE SetTendencyAtNode_Advection( myDGSEM, iEl, i, j, theTend  )
 ! S/R SetTendencyAtNode
 !  
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS(Advection), INTENT(inout) :: myDGSEM
   INTEGER, INTENT(in)             :: iEl
   INTEGER, INTENT(in)             :: i, j 
   REAL(prec), INTENT(in)          :: theTend(1:myDGSEM % nEq)

      CALL myDGSEM % sol(iEl) % SetTendencyAtNode( i, j, theTend )

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
   REAL(prec), INTENT(out)      :: theTend(0:myDGSEM % nS,0:myDGSEM % nP)

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
   REAL(prec), INTENT(in)          :: theTend(0:myDGSEM % nS,0:myDGSEM % nP)

      CALL myDGSEM % sol(iEl) % SetTendencyWithVarID( varID, theTend )

 END SUBROUTINE SetTendencyWithVarID_Advection
!
!
!
 SUBROUTINE GetTendencyAtNodeWithVarID_Advection( myDGSEM, iEl, i, j, varID, theTend  )
 ! S/R GetTendency
 !  
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS(Advection), INTENT(in) :: myDGSEM
   INTEGER, INTENT(in)          :: iEl
   INTEGER, INTENT(in)          :: i, j, varID 
   REAL(prec), INTENT(out)      :: theTend

      CALL myDGSEM % sol(iEl) % GetTendencyAtNodeWithVarID( i, j, varID, theTend )

 END SUBROUTINE GetTendencyAtNodeWithVarID_Advection
!
!
!
 SUBROUTINE SetTendencyAtNodeWithVarID_Advection( myDGSEM, iEl, i, j, varID, theTend  )
 ! S/R SetTendency
 !  
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS(Advection), INTENT(inout) :: myDGSEM
   INTEGER, INTENT(in)             :: iEl
   INTEGER, INTENT(in)             :: i, j, varID 
   REAL(prec), INTENT(in)          :: theTend

      CALL myDGSEM % sol(iEl) % SetTendencyAtNodeWithVarID( i, j, varID, theTend )

 END SUBROUTINE SetTendencyAtNodeWithVarID_Advection
!
! ----------------------------------------- Gradient --------------------------------------------- !
!
 SUBROUTINE GetGradient_Advection( myDGSEM, iEl, dcdx, dcdy  )
 ! S/R GetGradient
 !  
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS(Advection), INTENT(in) :: myDGSEM
   INTEGER, INTENT(in)          :: iEl
   REAL(prec), INTENT(out)      :: dcdx(0:myDGSEM % nS,0:myDGSEM % nP, 1:myDGSEM % nEq)
   REAL(prec), INTENT(out)      :: dcdy(0:myDGSEM % nS,0:myDGSEM % nP, 1:myDGSEM % nEq)

      CALL myDGSEM % dCdx(iEl) % GetSolution( dcdx )
      CALL myDGSEM % dCdy(iEl) % GetSolution( dcdy )

 END SUBROUTINE GetGradient_Advection
!
!
!
 SUBROUTINE SetGradient_Advection( myDGSEM, iEl, dcdx, dcdy  )
 ! S/R SetGradient
 !  
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS(Advection), INTENT(inout) :: myDGSEM
   INTEGER, INTENT(in)             :: iEl
   REAL(prec), INTENT(in)          :: dcdx(0:myDGSEM % nS,0:myDGSEM % nP, 1:myDGSEM % nEq)
   REAL(prec), INTENT(in)          :: dcdy(0:myDGSEM % nS,0:myDGSEM % nP, 1:myDGSEM % nEq)

      CALL myDGSEM % dCdx(iEl) % SetSolution( dcdx )
      CALL myDGSEM % dCdy(iEl) % SetSolution( dcdy )   

 END SUBROUTINE SetGradient_Advection
!
!
!
 SUBROUTINE GetGradientAtNode_Advection( myDGSEM, iEl, i, j, dcdx, dcdy  )
 ! S/R GetGradientAtNode
 !  
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS(Advection), INTENT(in) :: myDGSEM
   INTEGER, INTENT(in)          :: iEl
   INTEGER, INTENT(in)          :: i, j 
   REAL(prec), INTENT(out)      :: dcdx(1:myDGSEM % nEq), dcdy(1:myDGSEM % nEq)

      CALL myDGSEM % dCdx(iEl) % GetSolutionAtNode( i, j, dcdx )
      CALL myDGSEM % dCdy(iEl) % GetSolutionAtNode( i, j, dcdy )

 END SUBROUTINE GetGradientAtNode_Advection
!
!
!
 SUBROUTINE SetGradientAtNode_Advection( myDGSEM, iEl, i, j, dcdx, dcdy  )
 ! S/R SetGradientAtNode
 !  
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS(Advection), INTENT(inout) :: myDGSEM
   INTEGER, INTENT(in)             :: iEl
   INTEGER, INTENT(in)             :: i, j 
   REAL(prec), INTENT(in)          :: dcdx(1:myDGSEM % nEq), dcdy(1:myDGSEM % nEq)

      CALL myDGSEM % dCdx(iEl) % SetSolutionAtNode( i, j, dcdx )
      CALL myDGSEM % dCdy(iEl) % SetSolutionAtNode( i, j, dcdy )

 END SUBROUTINE SetGradientAtNode_Advection
!
!
!
SUBROUTINE GetGradientWithVarID_Advection( myDGSEM, iEl, varID, dcdx, dcdy  )
 ! S/R GetGradient
 !  
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS(Advection), INTENT(in) :: myDGSEM
   INTEGER, INTENT(in)          :: iEl
   INTEGER, INTENT(in)          :: varID 
   REAL(prec), INTENT(out)      :: dcdx(0:myDGSEM % nS,0:myDGSEM % nP)
   REAL(prec), INTENT(out)      :: dcdy(0:myDGSEM % nS,0:myDGSEM % nP)

      CALL myDGSEM % dCdx(iEl) % GetSolutionWithVarID( varID, dcdx )
      CALL myDGSEM % dCdy(iEl) % GetSolutionWithVarID( varID, dcdy )

 END SUBROUTINE GetGradientWithVarID_Advection
!
!
!
 SUBROUTINE SetGradientWithVarID_Advection( myDGSEM, iEl, varID, dcdx, dcdy  )
 ! S/R SetGradient
 !  
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS(Advection), INTENT(inout) :: myDGSEM
   INTEGER, INTENT(in)             :: iEl
   INTEGER, INTENT(in)             :: varID 
   REAL(prec), INTENT(in)          :: dcdx(0:myDGSEM % nS,0:myDGSEM % nP)
   REAL(prec), INTENT(in)          :: dcdy(0:myDGSEM % nS,0:myDGSEM % nP)

      CALL myDGSEM % dCdx(iEl) % SetSolutionWithVarID( varID, dcdx )
      CALL myDGSEM % dCdy(iEl) % SetSolutionWithVarID( varID, dcdy )

 END SUBROUTINE SetGradientWithVarID_Advection
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
   REAL(prec), INTENT(out)      :: theSolution(0:myDGSEM % nS,1:myDGSEM % nEq)

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
   REAL(prec), INTENT(in)          :: theSolution(0:myDGSEM % nS,1:myDGSEM % nEq)

      CALL myDGSEM % sol(iEl) % SetBoundarySolutionAtBoundary( boundary, theSolution )

 END SUBROUTINE SetBoundarySolutionAtBoundary_Advection
!
! ------------------------------- Boundary Gradient ---------------------------------------------- !
!
 SUBROUTINE GetBoundaryGradientAtBoundary_Advection( myDGSEM, iEl, boundary, dcdx, dcdy  )
 ! S/R GetBoundaryGradient
 !  
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS(Advection), INTENT(in) :: myDGSEM
   INTEGER, INTENT(in)          :: iEl
   INTEGER, INTENT(in)          :: boundary
   REAL(prec), INTENT(out)      :: dcdx(0:myDGSEM % nS,1:myDGSEM % nEq)
   REAL(prec), INTENT(out)      :: dcdy(0:myDGSEM % nS,1:myDGSEM % nEq)

      CALL myDGSEM % dcdx(iEl) % GetBoundarySolutionAtBoundary( boundary, dcdx )
      CALL myDGSEM % dcdy(iEl) % GetBoundarySolutionAtBoundary( boundary, dcdy )

 END SUBROUTINE GetBoundaryGradientAtBoundary_Advection
!
!
!
 SUBROUTINE SetBoundaryGradientAtBoundary_Advection( myDGSEM, iEl, boundary, dcdx, dcdy  )
 ! S/R SetBoundaryGradient
 !  
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS(Advection), INTENT(inout) :: myDGSEM
   INTEGER, INTENT(in)             :: iEl
   INTEGER, INTENT(in)             :: boundary
   REAL(prec), INTENT(in)          :: dcdx(0:myDGSEM % nS,1:myDGSEM % nEq)
   REAL(prec), INTENT(in)          :: dcdy(0:myDGSEM % nS,1:myDGSEM % nEq)

      CALL myDGSEM % dcdx(iEl) % SetBoundarySolutionAtBoundary( boundary, dcdx )
      CALL myDGSEM % dcdy(iEl) % SetBoundarySolutionAtBoundary( boundary, dcdy )

 END SUBROUTINE SetBoundaryGradientAtBoundary_Advection
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
   REAL(prec), INTENT(out)      :: theFlux(0:myDGSEM % nS,1:myDGSEM % nEq)

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
   REAL(prec), INTENT(in)          :: theFlux(0:myDGSEM % nS,1:myDGSEM % nEq)

      CALL myDGSEM % sol(iEl) % SetBoundaryFluxAtBoundary( boundary, theFlux )

 END SUBROUTINE SetBoundaryFluxAtBoundary_Advection
!
!
!
 SUBROUTINE GetBoundaryFluxAtBoundaryNode_Advection( myDGSEM, iEl, boundary, i, theFlux  )
 ! S/R GetBoundaryFlux
 !  
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS(Advection), INTENT(in) :: myDGSEM
   INTEGER, INTENT(in)          :: iEl
   INTEGER, INTENT(in)          :: boundary, i
   REAL(prec), INTENT(out)      :: theFlux(1:myDGSEM % nEq)

      CALL myDGSEM % sol(iEl) % GetBoundaryFluxAtBoundaryNode( boundary, i, theFlux )

 END SUBROUTINE GetBoundaryFluxAtBoundaryNode_Advection
!
!
!
 SUBROUTINE SetBoundaryFluxAtBoundaryNode_Advection( myDGSEM, iEl, boundary, i, theFlux  )
 ! S/R SetBoundaryFlux
 !  
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS(Advection), INTENT(inout) :: myDGSEM
   INTEGER, INTENT(in)             :: iEl
   INTEGER, INTENT(in)             :: boundary, i
   REAL(prec), INTENT(in)          :: theFlux(1:myDGSEM % nEq)

      CALL myDGSEM % sol(iEl) % SetBoundaryFluxAtBoundaryNode( boundary, i, theFlux )

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
   REAL(prec), INTENT(out)      :: relaxField(0:myDGSEM % nS, 0:myDGSEM % nP, 1:myDGSEM % nEq)
 
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
   REAL(prec), INTENT(in)          :: relaxField(0:myDGSEM % nS, 0:myDGSEM % nP, 1:myDGSEM % nEq)

      CALL myDGSEM % relax(iEl) % SetSolution( relaxField )

 END SUBROUTINE SetRelaxationField_Advection
!
!
!
 SUBROUTINE GetRelaxationFieldAtNode_Advection( myDGSEM, iEl, i, j, relaxField  )
 ! S/R GetRelaxationField
 !  
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS(Advection), INTENT(in) :: myDGSEM
   INTEGER, INTENT(in)          :: iEl, i, j
   REAL(prec), INTENT(out)      :: relaxField(1:myDGSEM % nEq)
 
      CALL myDGSEM % relax(iEl) % GetSolutionAtNode( i, j, relaxField )
      
 END SUBROUTINE GetRelaxationFieldAtNode_Advection
!
!
!
 SUBROUTINE SetRelaxationFieldAtNode_Advection( myDGSEM, iEl, i, j, relaxField  )
 ! S/R SetRelaxationFieldAtNode
 !  
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS(Advection), INTENT(inout) :: myDGSEM
   INTEGER, INTENT(in)                :: iEl, i, j
   REAL(prec), INTENT(in)             :: relaxField(1:myDGSEM % nEq)

      CALL myDGSEM % relax(iEl) % SetSolutionAtNode( i, j, relaxField )

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
   REAL(prec), INTENT(out)      :: rFac(0:myDGSEM % nS, 0:myDGSEM % nP, 1:myDGSEM % nEq)

      CALL myDGSEM % relaxFactor(iEl) % GetSolution( rFac )
          
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
   REAL(prec), INTENT(in)          :: rFac(0:myDGSEM % nS, 0:myDGSEM % nP, 1:myDGSEM % nEq)
   
      CALL myDGSEM % relaxFactor(iEl) % SetSolution( rFac )

 END SUBROUTINE SetRelaxationFactor_Advection
!
!
!
 SUBROUTINE GetRelaxationFactorAtNode_Advection( myDGSEM, iEl, i, j, rFac  )
 ! S/R GetRelaxationFactor
 !  
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS(Advection), INTENT(in) :: myDGSEM
   INTEGER, INTENT(in)          :: iEl, i, j
   REAL(prec), INTENT(out)      :: rFac(1:myDGSEM % nEq)
 
      CALL myDGSEM % relaxFactor(iEl) % GetSolutionAtNode( i, j, rFac )
      
 END SUBROUTINE GetRelaxationFactorAtNode_Advection
!
!
!
 SUBROUTINE SetRelaxationFactorAtNode_Advection( myDGSEM, iEl, i, j, rFac  )
 ! S/R SetRelaxationFactorAtNode
 !  
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS(Advection), INTENT(inout) :: myDGSEM
   INTEGER, INTENT(in)             :: iEl, i, j
   REAL(prec), INTENT(in)          :: rFac(1:myDGSEM % nEq)

      CALL myDGSEM % relaxFactor(iEl) % SetSolutionAtNode( i, j, rFac )

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
   INTEGER, INTENT(in)                :: iEl
   
      CALL myDGSEM % sol(iEl) % CalculateSolutionAtBoundaries( myDGSEM % dgStorage )
      
 END SUBROUTINE CalculateSolutionAtBoundaries_Advection
!
!
!
 SUBROUTINE CalculateGradientAtBoundaries_Advection( myDGSEM, iEl )
 ! S/R CalculateGradientAtBoundaries
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS(Advection), INTENT(inout) :: myDGSEM
   INTEGER, INTENT(in)             :: iEl
   
      CALL myDGSEM % dCdx(iEl) % CalculateSolutionAtBoundaries( myDGSEM % dgStorage )
      CALL myDGSEM % dCdy(iEl) % CalculateSolutionAtBoundaries( myDGSEM % dgStorage )
      
 END SUBROUTINE CalculateGradientAtBoundaries_Advection
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
 SUBROUTINE CalculateConcentrationGradient_Advection( myDGSEM, iEl ) 
 ! S/R CalculateConcentrationGradient
 ! 
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS(Advection), INTENT(inout) :: myDGSEM
   INTEGER, INTENT(in)                :: iEl
   ! LOCAL
   REAL(prec) :: dxds, dxdp, dyds, dydp, J, normLength, normVec(1:2)
   REAL(prec) :: h, sol, plv, rev
   REAL(prec) :: xF(1:2), yF(1:2)
   REAL(prec) :: pContFlux(0:myDGSEM % nS,1:2)
   REAL(prec) :: sContFlux(0:myDGSEM % nS,1:2)
   REAL(prec) :: pContFluxDer(0:myDGSEM % nS,1:2)
   REAL(prec) :: sContFluxDer(0:myDGSEM % nS,1:2)
   REAL(prec) :: dCdx(0:myDGSEM % nS,0:myDGSEM % nP)
   REAL(prec) :: dCdy(0:myDGSEM % nS,0:myDGSEM % nP)
   REAL(prec) :: dMatX(0:myDGSEM % nS,0:myDGSEM % nS)
   REAL(prec) :: dMatY(0:myDGSEM % nP,0:myDGSEM % nP)
   REAL(prec) :: qWeightX(0:myDGSEM % nS)
   REAL(prec) :: qWeightY(0:myDGSEM % nP)
   REAL(prec) :: lsouth(0:myDGSEM % nP)
   REAL(prec) :: lnorth(0:myDGSEM % nP)
   REAL(prec) :: least(0:myDGSEM % nS)
   REAL(prec) :: lwest(0:myDGSEM % nS)
   REAL(prec) :: csouth(0:myDGSEM % nP,1:myDGSEM % nEq)
   REAL(prec) :: cnorth(0:myDGSEM % nP,1:myDGSEM % nEq)
   REAL(prec) :: ceast(0:myDGSEM % nS,1:myDGSEM % nEq)
   REAL(prec) :: cwest(0:myDGSEM % nS,1:myDGSEM % nEq)
   REAL(prec) :: fL(1:2), fR(1:2), tt(1:myDGSEM % nEq)
   INTEGER    :: iS, iP, nS, nP, iEq, iDir, nEq

     nEq = myDGSEM % nEq
     CALL myDGSEM % dgStorage % GetNumberOfNodes( nS, nP )
     CALL myDGSEM % dgStorage % GetDerivativeMatrix( dMatX, dMatY )
     CALL myDGSEM % dgStorage % GetQuadratureWeights( qWeightX, qWeightY )
     CALL myDGSEM % dgStorage % GetSouthernInterpolants( lsouth )
     CALL myDGSEM % dgStorage % GetNorthernInterpolants( lnorth )
     CALL myDGSEM % dgStorage % GetEasternInterpolants( least )
     CALL myDGSEM % dgStorage % GetWesternInterpolants( lwest )
         
     CALL myDGSEM % GetBoundarySolutionAtBoundary( iEl, south, cSouth )
     CALL myDGSEM % GetBoundarySolutionAtBoundary( iEl, east, cEast )
     CALL myDGSEM % GetBoundarySolutionAtBoundary( iEl, north, cNorth )
     CALL myDGSEM % GetBoundarySolutionAtBoundary( iEl, west, cWest ) 

      DO iEq = 1, myDGSEM % nEq

         DO iP = 0,nP ! Loop over the y-points
            DO iS = 0,nS ! Loop over the x-points

               ! Get the metric terms to calculate the contravariant flux
               CALL myDGSEM % mesh % GetCovariantMetricsAtNode( iEl, dxds, dxdp, dyds, dydp, iS, iP )
               ! Get the depth - needed for wave speed calculation 
               CALL myDGSEM % GetSolutionAtNodeWithVarID( iEl, iS, iP, iEq, sol )

              ! Calculate the x and y fluxes
               xF(1:2) = (/ sol, ZERO /)
               yF(1:2) = (/ ZERO, sol /)

               !And now the contravariant flux
               sContFlux(iS,1:2) = dydp*xF - dxdp*yF

            ENDDO

            ! East edge
            CALL myDGSEM % dCdx(iEl) % GetBoundaryFluxAtBoundaryNode(  east, iP, tt(1:nEq) )
            fR(1) = tt(iEq) !( cEast(iP,iEq)*normVec )*normLength
            CALL myDGSEM % dCdy(iEl) % GetBoundaryFluxAtBoundaryNode(  east, iP, tt(1:nEq) )
            fR(2) = tt(iEq)

            ! West edge
            CALL myDGSEM % dCdx(iEl) % GetBoundaryFluxAtBoundaryNode(  west, iP, tt(1:nEq) )
            fL(1) = tt(iEq) !( cEast(iP,iEq)*normVec )*normLength
            CALL myDGSEM % dCdy(iEl) % GetBoundaryFluxAtBoundaryNode(  west, iP, tt(1:nEq) )
            fL(2) = tt(iEq)
         
            ! At this y-level, calculate the DG-advective derivative
            DO iDir = 1, 2
            
               sContFluxDer(0:nS,iDir) = DGSystemDerivative(  nS, dMatX, qWeightX, fL(iDir), fR(iDir), &
                                                              sContFlux(:,iDir), lwest, least )
            ENDDO

            DO iS = 0, nS ! Loop over the x-points
               dCdx(iS, iP) = sContFluxDer(iS,1)
               dCdy(iS, iP) = sContFluxDer(iS,2)
            ENDDO ! Loop over the number of equations

         ENDDO ! iP  


         DO iS = 0,nS ! Loop over the x-points
           DO iP = 0,nP ! Loop over the y-points

              ! Get the metric terms to calculate the contravariant flux
               CALL myDGSEM % mesh % GetCovariantMetricsAtNode( iEl, dxds, dxdp, dyds, dydp, iS, iP )

               CALL myDGSEM % GetSolutionAtNodeWithVarID( iEl, iS, iP, iEq, sol )

              ! Calculate the x and y fluxes
               xF(1:2) = (/ sol, ZERO /)
               yF(1:2) = (/ ZERO, sol /)

              !And now the contravariant flux
              pContFlux(iP,:) = -dyds*xF + dxds*yF

            ENDDO ! iP, Loop over the y-points

            ! North
            CALL myDGSEM % dCdx(iEl) % GetBoundaryFluxAtBoundaryNode(  north, iS, tt(1:nEq) )
            fR(1) = tt(iEq) !( cEast(iP,iEq)*normVec )*normLength
            CALL myDGSEM % dCdy(iEl) % GetBoundaryFluxAtBoundaryNode(  north, iS, tt(1:nEq) )
            fR(2) = tt(iEq)

            ! South
            CALL myDGSEM % dCdx(iEl) % GetBoundaryFluxAtBoundaryNode(  south, iS, tt(1:nEq) )
            fL(1) = tt(iEq) !( cEast(iP,iEq)*normVec )*normLength
            CALL myDGSEM % dCdy(iEl) % GetBoundaryFluxAtBoundaryNode(  south, iS, tt(1:nEq) )
            fL(2) = tt(iEq)

            ! At this x-level, calculate the y-DG-advective derivative
            DO iDir = 1, 2
               pContFluxDer(0:nP,iDir) = DGSystemDerivative( nP, dMatY, qWeightY, fL(iDir), fR(iDir), &
                                                             pContFlux(:,iDir), lsouth, lnorth  )
            ENDDO
         
            DO iP = 0, nP ! Loop over the y-points
                CALL myDGSEM % mesh % GetJacobianAtNode( iEl, J, iS, iP )
                dCdX(iS, iP) = ( dCdX(iS, iP)+ pContFluxDer(iP,1) )/J
                dCdY(iS, iP) = ( dCdY(iS, iP)+ pContFluxDer(iP,2) )/J
            ENDDO ! Loop over the number of equations

      
         ENDDO 
          
         CALL myDGSEM % SetGradientWithVarID( iEl, iEq, dcdx, dcdy  )
      ENDDO

 END SUBROUTINE CalculateConcentrationGradient_Advection
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
                     1:myDGSEM % nEq,&
                     1:myDGSEM % mesh % nElems) 
   REAL(prec) :: dSdt(0:myDGSEM % nS,&
                     0:myDGSEM % nP,&
                     1:myDGSEM % nEq )
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
   
   
       
 END SUBROUTINE ForwardStepRK3_Advection
!
!
!
 SUBROUTINE GlobalTimeDerivative_Advection( myDGSEM, tn ) 
 ! S/R GlobalTimeDerivative_Advection
 ! 
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS(Advection), INTENT(inout) :: myDGSEM
   REAL(prec), INTENT(in)          :: tn
   ! Local
   INTEGER :: iEl, iEdge

      ! CALL myDGSEM % LoadVelocityField( tn )
!$OMP PARALLEL
      ! Calculate the solution at the boundaries
!$OMP DO
      DO iEl = 1, myDGSEM % mesh % nElems

         CALL myDGSEM % CalculateSolutionAtBoundaries( iEl ) 
         CALL myDGSEM % CalculateVelocityAtBoundaries( iEl )

      ENDDO 
!$OMP END DO 
!$OMP FLUSH( myDGSEM )

      ! Calculate the flux along the element edges
!$OMP DO
      DO iEdge = 1, myDGSEM % mesh % nEdges
         CALL myDGSEM % EdgeGradientFlux( iEdge, tn )
      ENDDO 
!$OMP END DO
!$OMP FLUSH( myDGSEM ) 


!$OMP DO
      DO iEl = 1, myDGSEM % mesh % nElems

         CALL myDGSEM % CalculateConcentrationGradient( iEl )
         CALL myDGSEM % CalculateGradientAtBoundaries( iEl )

      ENDDO 
!$OMP END DO 
!$OMP FLUSH( myDGSEM )


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

 
 END SUBROUTINE GlobalTimeDerivative_Advection
!
!
! 
 SUBROUTINE EdgeFlux_Advection( myDGSEM, iEdge, tn )
 ! S/R EdgeFlux
 ! 
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS( Advection ), INTENT(inout) :: myDGSEM
   INTEGER, INTENT(in)               :: iEdge  
   REAL(prec), INTENT(in)            :: tn
   ! Local
   INTEGER :: k, iNode, nEq
   INTEGER :: e1, s1, e2, s2, iError, start, inc, nS, nP
   REAL(prec) :: flux(1:myDGSEM % nEq)
   REAL(prec) :: inState(0:myDGSEM % nS,1:myDGSEM % nEq)
   REAL(prec) :: exState(0:myDGSEM % nS,1:myDGSEM % nEq)
   REAL(prec) :: DcDx(0:myDGSEM % nS,1:myDGSEM % nEq)
   REAL(prec) :: exDcDx(0:myDGSEM % nS,1:myDGSEM % nEq)
   REAL(prec) :: DcDy(0:myDGSEM % nS,1:myDGSEM % nEq)
   REAL(prec) :: exDcDy(0:myDGSEM % nS,1:myDGSEM % nEq)
   REAL(prec) :: x, y
   REAL(prec) :: u(0:myDGSEM % nS), v(0:myDGSEM % nS) 
   REAL(prec) :: nHat(1:2), nHatLength , kappa
    
      nEQ = myDGSEM % params % nTracers

      CALL myDGSEM % dgStorage % GetNumberOfNodes( nS, nP )

      CALL myDGSEM % mesh % GetEdgePrimaryElementID( iEdge, e1 )
      CALL myDGSEM % mesh % GetEdgePrimaryElementSide( iEdge, s1 )
      CALL myDGSEM % mesh % GetEdgeSecondaryElementID( iEdge, e2 )
      CALL myDGSEM % mesh % GetEdgeSecondaryElementSide( iEdge, s2 )
      
      s2 = ABS(s2)

      CALL myDGSEM % mesh % GetEdgeStart( iEdge, start )
      CALL myDGSEM % mesh % GetEdgeIncrement( iEdge, inc )
      CALL myDGSEM % GetBoundarySolutionAtBoundary( e1, s1, inState )
      CALL myDGSEM % GetBoundaryGradientAtBoundary( e1, s1, dCdx, dCdy )
      CALL myDGSEM % GetVelocityAtBoundary( e1, s1, u, v )
      kappa = myDGSEM % params % kappa

      IF( e2 > 0 )then ! this is an interior edge
      
         k = start-inc
         CALL myDGSEM % GetBoundarySolutionAtBoundary( e2, s2, exState )
         CALL myDGSEM % GetBoundaryGradientAtBoundary( e2, s2, exDcDx, exDcDy )
         
         DO iNode = 0, nS ! Loop over the nodes
            
            CALL myDGSEM % mesh % GetBoundaryNormalAtNode( e1, nHat, nHatLength, iNode, s1 ) ! Get nHat
           
            ! Calculate the RIEMANN flux
            dCdx(iNode,:) = ( dCdx(iNode,:) + exDcDx(k,:) )*HALF
            dCdy(iNode,:) = ( dCdy(iNode,:) + exDcDy(k,:) )*HALF

            flux = RiemannSolver( inState(iNode,:), exState(k,:), &
                                  u(iNode), v(iNode), kappa, &
                                  dCdx(iNode,:), dCdy(iNode,:), nHat, nEq )*nHatLength

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

           ! IF( e2 == NO_NORMAL_FLOW )THEN
           !    flux = ZERO
           ! ELSE
               ! Calculate the RIEMANN flux
               flux = RiemannSolver( inState(iNode,1:nEq), exState(iNode,1:nEq), &
                                     u(iNode), v(iNode), kappa, &
                                     dCdx(iNode,:), dCdy(iNode,:), &
                                     nHat, nEq )*nHatLength
            !ENDIF
 
            CALL myDGSEM % SetBoundaryFluxAtBoundaryNode( e1, s1, iNode, flux )

         ENDDO ! iNode, loop over the nodes on this edge

      ENDIF ! Choosing the edge type
 
 END SUBROUTINE EdgeFlux_Advection
!
!
!
SUBROUTINE EdgeGradientFlux_Advection( myDGSEM, iEdge, tn )
 ! S/R EdgeGradientFlux
 ! 
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS( Advection ), INTENT(inout) :: myDGSEM
   INTEGER, INTENT(in)               :: iEdge  
   REAL(prec), INTENT(in)            :: tn
   ! Local
   INTEGER :: k, iNode, nEq
   INTEGER :: e1, s1, e2, s2, iError, start, inc, nS, nP
   REAL(prec) :: flux(1:2,1:myDGSEM % nEq)
   REAL(prec) :: inState(0:myDGSEM % nS,1:myDGSEM % nEq)
   REAL(prec) :: exState(0:myDGSEM % nS,1:myDGSEM % nEq)
   REAL(prec) :: x, y
   REAL(prec) :: nHat(1:2), nHatLength
    
      nEQ = myDGSEM % params % nTracers

      CALL myDGSEM % dgStorage % GetNumberOfNodes( nS, nP )

      CALL myDGSEM % mesh % GetEdgePrimaryElementID( iEdge, e1 )
      CALL myDGSEM % mesh % GetEdgePrimaryElementSide( iEdge, s1 )
      CALL myDGSEM % mesh % GetEdgeSecondaryElementID( iEdge, e2 )
      CALL myDGSEM % mesh % GetEdgeSecondaryElementSide( iEdge, s2 )
      
      s2 = ABS(s2)

      CALL myDGSEM % mesh % GetEdgeStart( iEdge, start )
      CALL myDGSEM % mesh % GetEdgeIncrement( iEdge, inc )
      
      CALL myDGSEM % GetBoundarySolutionAtBoundary( e1, s1, inState )

      IF( e2 > 0 )then ! this is an interior edge
      
         k = start-inc
         CALL myDGSEM % GetBoundarySolutionAtBoundary( e2, s2, exState )
         
         DO iNode = 0, nS ! Loop over the nodes
            
            CALL myDGSEM % mesh % GetBoundaryNormalAtNode( e1, nHat, nHatLength, iNode, s1 ) ! Get nHat
           
            
            flux(1,:) = HALF*(inState(iNode,:) + exState(k,:) )*nHat(1)*nHatLength
            flux(2,:) = HALF*(inState(iNode,:) + exState(k,:) )*nHat(2)*nHatLength

            ! Store the flux for the elements which share this edge
            CALL myDGSEM % dCdx(e1) % SetBoundaryFluxAtBoundaryNode(  s1, iNode, flux(1,:) )
            CALL myDGSEM % dCdx(e2) % SetBoundaryFluxAtBoundaryNode(  s2, k, -flux(1,:) )
            CALL myDGSEM % dCdy(e1) % SetBoundaryFluxAtBoundaryNode(  s1, iNode, flux(2,:) )
            CALL myDGSEM % dCdy(e2) % SetBoundaryFluxAtBoundaryNode(  s2, k, -flux(2,:) )
             
            k = k + inc

         ENDDO ! iNode, loop over the nodes
         
      ELSE ! this edge is a boundary edge

         DO iNode = 0, nS ! loop over the nodes on this edge
       
            CALL myDGSEM % mesh % GetBoundaryNormalAtNode( e1, nHat, nHatLength, iNode, s1 ) ! Get nHat
            
            ! Get the boundary point locations
            CALL myDGSEM % mesh % GetBoundaryLocationAtNode( e1, x, y, iNode, s1 )

            ! Calculate the external state
            exState(iNode,1:nEq) = GetExternalState( nEq, x, y, tn, e2, inState(iNode,1:nEq))
               
            flux(1,:) = HALF*(inState(iNode,:) + exState(iNode,:) )*nHat(1)*nHatLength
            flux(2,:) = HALF*(inState(iNode,:) + exState(iNode,:) )*nHat(2)*nHatLength

            ! Store the flux for the elements which share this edge
            CALL myDGSEM % dCdx(e1) % SetBoundaryFluxAtBoundaryNode(  s1, iNode, flux(1,:) )
            CALL myDGSEM % dCdy(e1) % SetBoundaryFluxAtBoundaryNode(  s1, iNode, flux(2,:) )

         ENDDO ! iNode, loop over the nodes on this edge

      ENDIF ! Choosing the edge type
 
 END SUBROUTINE EdgeGradientFlux_Advection
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
   REAL(prec) :: kappa
   INTEGER    :: iS, iP, nS, nP, iEq

     CALL myDGSEM % dgStorage % GetNumberOfNodes( nS, nP )
     CALL myDGSEM % dgStorage % GetDerivativeMatrix( dMatX, dMatY )
     CALL myDGSEM % dgStorage % GetQuadratureWeights( qWeightX, qWeightY )
     CALL myDGSEM % dgStorage % GetSouthernInterpolants( lsouth )
     CALL myDGSEM % dgStorage % GetNorthernInterpolants( lnorth )
     CALL myDGSEM % dgStorage % GetEasternInterpolants( least )
     CALL myDGSEM % dgStorage % GetWesternInterpolants( lwest )
     CALL myDGSEM % GetVelocity( iEl, u, v )

     kappa = myDGSEM % params % kappa

     DO iP = 0,nP ! Loop over the y-points

        DO iS = 0,nS ! Loop over the x-points

           ! Get the metric terms to calculate the contravariant flux
           CALL myDGSEM % mesh % GetCovariantMetricsAtNode( iEl, dxds, dxdp, dyds, dydp, iS, iP )

           ! Get the solution at x,y
           CALL myDGSEM % GetSolutionAtNode( iEl, iS, iP, sol )
           CALL myDGSEM % GetGradientAtNode( iEl, iS, iP, dcdx, dcdy )
           ! Calculate the x and y fluxes
           xF = XFlux(  tn, myDGSEM % nEq, u(iS,iP), v(iS,iP), kappa, dcdx, dcdy, sol )
           yF = YFlux(  tn, myDGSEM % nEq, u(iS,iP), v(iS,iP), kappa, dcdx, dcdy, sol )

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
           CALL myDGSEM % GetGradientAtNode( iEl, iS, iP, dcdx, dcdy )


           ! Calculate the x and y fluxes
           xF = XFlux(  tn, myDGSEM % nEq, u(iS,iP), v(iS,iP), kappa, dcdx, dcdy, sol )
           yF = YFlux(  tn, myDGSEM % nEq, u(iS,iP), v(iS,iP), kappa, dcdx, dcdy, sol )

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
           CALL myDGSEM % GetRelaxationFieldAtNode( iEl, iS, iP, relaxSol )
           CALL myDGSEM % GetRelaxationFactorAtNode( iEl, iS, iP, rFac )
           CALL myDGSEM % mesh % GetPositionAtNode( iEl, x, y, iS, iP )
         
                  
           tend(iS,iP,:) = tend(iS,iP,:)+ Source( tn, sol, relaxSol, rFac, &
                                                  x, y, myDGSEM % nEq  )

        ENDDO ! iP, Loop over the y-points
     ENDDO ! iS, Loop over the x-points

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
 FUNCTION RiemannSolver( inState, outState, u, v, kappa, dcdx, dcdy, nHat, nEq ) RESULT( numFlux )
 ! FUNCTION RiemannSolver 
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   INTEGER    :: nEq
   REAL(prec) :: inState(1:nEq)
   REAL(prec) :: outState(1:nEq)
   REAL(prec) :: dcdx(1:nEq), dcdy(1:nEq)
   REAL(prec) :: u, v, kappa
   REAL(prec) :: nHat(1:2) ! normal direction
   REAL(prec) :: numFlux(1:nEq)
   
   ! LOCAL
   REAL(prec) :: uN, dcdn(1:nEq)
   REAL(prec) :: jump(1:nEq), aS(1:nEq)
   REAL(prec) :: uNorm, fac


      jump = outState - inState
    
      uN =  u*nHat(1) + v*nHat(2)

      numFlux = HALF*( uN*( outState + inState ) - abs(uN)*jump )

      ! Diffusive part
      dcdn = dcdx*nHat(1) + dcdy*nHat(2)
      numFlux = numFlux - kappa*dcdn


 END FUNCTION RiemannSolver
!
!
! 
 FUNCTION XFlux(  tn, nEq, u, v, kappa,dcdx, dcdy,  solAtX ) RESULT( fx )
 ! FUNCTION XFlux
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE 
   INTEGER    :: nEq
   REAL(prec) :: tn
   REAL(prec) :: u, v, kappa
   REAL(prec) :: solAtX(1:nEq)
   REAL(prec) :: dcdx(1:nEq), dcdy(1:nEq)
   REAL(prec) :: fx(1:nEq)

      fx = solAtX*u  - kappa*dcdx

 END FUNCTION XFlux
!
!
!
 FUNCTION YFlux(  tn, nEq, u, v, kappa, dcdx, dcdy, solAtX ) RESULT( fy )
 ! FUNCTION YFlux
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE 
   INTEGER    :: nEq
   REAL(prec) :: tn
   REAL(prec) :: u, v, kappa
   REAL(prec) :: solAtX(1:nEq)
   REAL(prec) :: dcdx(1:nEq), dcdy(1:nEq)
   REAL(prec) :: fy(1:nEq)

      fy = solAtX*v - kappa*dcdy

 END FUNCTION YFlux
!
!
!            
 FUNCTION Source( tn, sol, relax, rfac, x, y, nEq ) RESULT( q )
 ! FUNCTION Source
 ! 
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   INTEGER        :: nEq
   REAL(prec)     :: tn
   REAL(prec)     :: sol(1:nEq), relax(1:nEq), rFac(1:nEq)
   REAL(prec)     :: x, y
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

       IF( bcFlag == NO_NORMAL_FLOW )THEN
      
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
 SUBROUTINE CoarseToFine_Advection( myDGSEM, iEl, x, y, u, v, c, dcdx, dcdy )
 ! CoarseToFine
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS( Advection ), INTENT(in) :: myDGSEM
   INTEGER, INTENT(in)            :: iEl
   REAL(prec), INTENT(out)        :: x(0:myDGSEM % nPlot, 0:myDGSEM % nPlot)
   REAL(prec), INTENT(out)        :: y(0:myDGSEM % nPlot, 0:myDGSEM % nPlot)
   REAL(prec), INTENT(out)        :: u(0:myDGSEM % nPlot, 0:myDGSEM % nPlot)
   REAL(prec), INTENT(out)        :: v(0:myDGSEM % nPlot, 0:myDGSEM % nPlot)
   REAL(prec), INTENT(out)        :: c(0:myDGSEM % nPlot, 0:myDGSEM % nPlot,1:myDGSEM % nEq)
   REAL(prec), INTENT(out)        :: dcdx(0:myDGSEM % nPlot, 0:myDGSEM % nPlot,1:myDGSEM % nEq)
   REAL(prec), INTENT(out)        :: dcdy(0:myDGSEM % nPlot, 0:myDGSEM % nPlot,1:myDGSEM % nEq)
   ! Local
   REAL(prec) :: localX(0:myDGSEM % nS, 0:myDGSEM % nP)
   REAL(prec) :: localY(0:myDGSEM % nS, 0:myDGSEM % nP)
   REAL(prec) :: localU(0:myDGSEM % nS, 0:myDGSEM % nP)
   REAL(prec) :: localv(0:myDGSEM % nS, 0:myDGSEM % nP)
   REAL(prec) :: sol(0:myDGSEM % nS, 0:myDGSEM % nP,1:myDGSEM % nEq)
   REAL(prec) :: dsdx(0:myDGSEM % nS, 0:myDGSEM % nP,1:myDGSEM % nEq)
   REAL(prec) :: dsdy(0:myDGSEM % nS, 0:myDGSEM % nP,1:myDGSEM % nEq)
   INTEGER    :: iEq


      CALL myDGSEM % GetSolution( iEl, sol )
      CALL myDGSEM % GetGradient( iEl, dsdx, dsdy )
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
                                                           c(:,:,iEq) )  

         CALL myDGSEM % dgStorage % interp % CoarseToFine( dsdx(:,:,iEq), &
                                                           myDGSEM % plMatS, &
                                                           myDGSEM % plMatP, &
                                                           myDGSEM % nPlot, &
                                                           myDGSEM % nPlot, &
                                                           dcdx(:,:,iEq) ) 

         CALL myDGSEM % dgStorage % interp % CoarseToFine( dsdy(:,:,iEq), &
                                                           myDGSEM % plMatS, &
                                                           myDGSEM % plMatP, &
                                                           myDGSEM % nPlot, &
                                                           myDGSEM % nPlot, &
                                                           dcdy(:,:,iEq) ) 
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
  REAL(prec)  :: x(0:myDGSEM % nPlot,0:myDGSEM % nPlot)
  REAL(prec)  :: y(0:myDGSEM % nPlot,0:myDGSEM % nPlot)
  REAL(prec)  :: u(0:myDGSEM % nPlot,0:myDGSEM % nPlot)
  REAL(prec)  :: v(0:myDGSEM % nPlot,0:myDGSEM % nPlot)
  REAL(prec)  :: c(0:myDGSEM % nPlot,0:myDGSEM % nPlot, 1:myDGSEM % nEq)
  REAL(prec)  :: dcdx(0:myDGSEM % nPlot,0:myDGSEM % nPlot, 1:myDGSEM % nEq)
  REAL(prec)  :: dcdy(0:myDGSEM % nPlot,0:myDGSEM % nPlot, 1:myDGSEM % nEq)
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
             FILE= 'Advection.tec', &
             FORM='formatted', &
             STATUS='replace')  
    ENDIF
    
    WRITE(fUnit,*) 'VARIABLES = "X", "Y", "U", "V", "C1", "C2", "C3", "C4" '
 
    DO iEl = 1, myDGsem % mesh % nElems

       CALL myDGSEM % CoarseToFine( iEl, x, y, u, v, c, dcdx, dcdy )
        WRITE(zoneID,'(I5.5)') iEl
        WRITE(fUnit,*) 'ZONE T="el'//trim(zoneID)//'", I=',nPlot+1,', J=', nPlot+1,',F=POINT'

        DO iP = 0, nPlot
           DO iS = 0, nPlot
              WRITE (fUnit,*)  x( iS, iP ), y( iS, iP ), &
                               u(iS,iP), v(iS,iP), &
                               c(iS,iP,1), c(iS,iP,2), &
                               c(iS,iP,3), c(iS,iP,4)
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
           FILE='Advection.'//iterChar//'.pickup', &
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

        CALL myDGSEM % GetRelaxationField( iEl, sol )
        DO iEq = 1, myDGSEM % nEq
           WRITE( fUnit, REC=thisRec )sol(:,:,iEq) 
           thisRec = thisRec+1
        ENDDO

        CALL myDGSEM % GetRelaxationFactor( iEl, sol )
        DO iEq = 1, myDGSEM % nEq
           WRITE( fUnit, REC=thisRec )sol(:,:,iEq) 
           thisRec = thisRec+1
        ENDDO

        CALL myDGSEM % GetVelocity( iEl, u, v )
        WRITE( fUnit, REC=thisRec )u 
        thisRec = thisRec+1
        WRITE( fUnit, REC=thisRec )v 
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
           FILE='Advection.'//iterChar//'.pickup', &
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
        
        DO iEq = 1, myDGSEM % nEq
           READ( fUnit, REC=thisRec )sol(:,:,iEq) 
           thisRec = thisRec+1
        ENDDO
        CALL myDGSEM % SetRelaxationField( iEl, sol )

        
        DO iEq = 1, myDGSEM % nEq
           READ( fUnit, REC=thisRec )sol(:,:,iEq) 
           thisRec = thisRec+1
        ENDDO
        CALL myDGSEM % SetRelaxationFactor( iEl, sol )

        READ( fUnit, REC=thisRec )u 
        thisRec = thisRec+1
        READ( fUnit, REC=thisRec )v 
        thisRec = thisRec+1
        CALL myDGSEM % SetVelocity( iEl, u, v )
        
     ENDDO

     CLOSE(UNIT=fUnit)

     
 END SUBROUTINE ReadPickup_Advection
!
!
! 
 END MODULE AdvectionClass



