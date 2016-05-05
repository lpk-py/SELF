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
   REAL(prec) :: h(0:myDGSEM % nS, 0:myDGSEM % nP)
   
      CALL myDGSEM % sol(iEl) % GetSolutionWithVarID( 1, u )
      CALL myDGSEM % sol(iEl) % GetSolutionWithVarID( 2, v )
!      CALL myDGSEM % sol(iEl) % GetSolutionWithVarID( 3, h )
!      u = u/h
!      v = v/h
      
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
      CALL myDGSEM % sol(iEl) % SetSolutionWithVarID( 2, v )

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
   REAL(prec) :: h(0:myDGSEM % nS)
   
      CALL myDGSEM % sol(iEl) % GetBoundarySolutionAtBoundaryWithVarID( boundary, 1, u )
      CALL myDGSEM % sol(iEl) % GetBoundarySolutionAtBoundaryWithVarID( boundary, 2, v )

      CALL myDGSEM % sol(iEl) % GetBoundarySolutionAtBoundaryWithVarID( boundary, 3, h )
      u = u/h
      v = v/h
      
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
