MODULE DGSEMSolutionStorageClass_2D
! DGSEMSolutionStorageClass_2D.f90
! 
! Module History 
! 
! o  (v1.0) March 2014
!    
! o  (v2.1) December 10, 2015
!
!  
! =======================================================================================
! src/common/
USE ConstantsDictionary
USE ModelPrecision
! src/nodal
USE NodalStorage_2D_Class

IMPLICIT NONE


    TYPE DGSEMSolution_2D
      INTEGER                 :: nEq, nS, nP, nMax
      REAL(prec), allocatable :: Solution(:,:,:)
      REAL(prec), allocatable :: tendency(:,:,:)
      REAL(prec), allocatable :: boundarySolution(:,:,:) ! Indexed over the boundaries. 
      REAL(prec), allocatable :: boundaryFlux(:,:,:)     ! Indexed over the boundaries

      !** Side ID's : 1 is south, 2 is east, 3 is north, 4 is west

      CONTAINS

      ! CONSTRUCTORS/DESTRUCTORS
      PROCEDURE :: Build => Build_DGSEMSolution_2D
      PROCEDURE :: Trash => Trash_DGSEMSolution_2D

      ! ACCESSORS
      PROCEDURE :: GetNumberOfEquations => GetNumberOfEquations_DGSEMSolution_2D
      PROCEDURE :: SetNumberOfEquations => SetNumberOfEquations_DGSEMSolution_2D
      
      PROCEDURE :: GetNumberOfNodes => GetNumberOfNodes_DGSEMSolution_2D
      PROCEDURE :: SetNumberOfNodes => SetNumberOfNodes_DGSEMSolution_2D

      PROCEDURE :: GetSolution => GetSolution_DGSEMSolution_2D
      PROCEDURE :: SetSolution => SetSolution_DGSEMSolution_2D
      PROCEDURE :: GetSolutionAtNode => GetSolutionAtNode_DGSEMSolution_2D
      PROCEDURE :: SetSolutionAtNode => SetSolutionAtNode_DGSEMSolution_2D
      PROCEDURE :: GetSolutionWithVarID => GetSolutionWithVarID_DGSEMSolution_2D
      PROCEDURE :: SetSolutionWithVarID => SetSolutionWithVarID_DGSEMSolution_2D
      PROCEDURE :: GetSolutionAtNodeWithVarID => GetSolutionAtNodeWithVarID_DGSEMSolution_2D
      PROCEDURE :: SetSolutionAtNodeWithVarID => SetSolutionAtNodeWithVarID_DGSEMSolution_2D
      
      PROCEDURE :: GetTendency => GetTendency_DGSEMSolution_2D
      PROCEDURE :: SetTendency => SetTendency_DGSEMSolution_2D
      PROCEDURE :: GetTendencyAtNode => GetTendencyAtNode_DGSEMSolution_2D
      PROCEDURE :: SetTendencyAtNode => SetTendencyAtNode_DGSEMSolution_2D
      PROCEDURE :: GetTendencyWithVarID => GetTendencyWithVarID_DGSEMSolution_2D
      PROCEDURE :: SetTendencyWithVarID => SetTendencyWithVarID_DGSEMSolution_2D
      PROCEDURE :: GetTendencyAtNodeWithVarID => GetTendencyAtNodeWithVarID_DGSEMSolution_2D
      PROCEDURE :: SetTendencyAtNodeWithVarID => SetTendencyAtNodeWithVarID_DGSEMSolution_2D
      
      PROCEDURE :: GetBoundarySolution => GetBoundarySolution_DGSEMSolution_2D
      PROCEDURE :: SetBoundarySolution => SetBoundarySolution_DGSEMSolution_2D
      PROCEDURE :: GetBoundarySolutionAtBoundary => GetBoundarySolutionAtBoundary_DGSEMSolution_2D
      PROCEDURE :: SetBoundarySolutionAtBoundary => SetBoundarySolutionAtBoundary_DGSEMSolution_2D
      PROCEDURE :: GetBoundarySolutionAtBoundaryWithVarID => GetBoundarySolutionAtBoundaryWithVarID_DGSEMSolution_2D
      PROCEDURE :: SetBoundarySolutionAtBoundaryWithVarID => SetBoundarySolutionAtBoundaryWithVarID_DGSEMSolution_2D

      PROCEDURE :: GetBoundaryFlux => GetBoundaryFlux_DGSEMSolution_2D
      PROCEDURE :: SetBoundaryFlux => SetBoundaryFlux_DGSEMSolution_2D
      PROCEDURE :: GetBoundaryFluxAtBoundary => GetBoundaryFluxAtBoundary_DGSEMSolution_2D
      PROCEDURE :: SetBoundaryFluxAtBoundary => SetBoundaryFluxAtBoundary_DGSEMSolution_2D
      PROCEDURE :: GetBoundaryFluxAtBoundaryNode => GetBoundaryFluxAtBoundaryNode_DGSEMSolution_2D
      PROCEDURE :: SetBoundaryFluxAtBoundaryNode => SetBoundaryFluxAtBoundaryNode_DGSEMSolution_2D
      PROCEDURE :: GetBoundaryFluxAtBoundaryWithVarID => GetBoundaryFluxAtBoundaryWithVarID_DGSEMSolution_2D
      PROCEDURE :: SetBoundaryFluxAtBoundaryWithVarID => SetBoundaryFluxAtBoundaryWithVarID_DGSEMSolution_2D 
      
      PROCEDURE :: CalculateSolutionAtBoundaries => CalculateSolutionAtBoundaries_DGSEMSolution_2D
    END TYPE DGSEMSolution_2D

 
 CONTAINS
!
!
!==================================================================================================!
!------------------------------- Manual Constructors/Destructors ----------------------------------!
!==================================================================================================!
!
!
 SUBROUTINE Build_DGSEMSolution_2D( myDGS, nS, nP, nEq )
 ! S/R Build_DGSEMSolution_2D
 !  
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS(DGSEMSolution_2D), INTENT(inout) :: myDGS
   INTEGER, INTENT(in)                    :: nS, nP, nEq

      CALL myDGS % SetNumberOfEquations( nEq )
      CALL myDGS % SetNumberOfNodes( nS, nP )
      
      ALLOCATE( myDGS % solution(0:nS,0:nP,1:nEq) )
      ALLOCATE( myDGS % tendency(0:nS,0:nP,1:nEq) )
      ALLOCATE( myDGS % boundarySolution(0:max(nS,nP),1:nEq,1:4) ) 
      ALLOCATE( myDGS % boundaryFlux(0:max(nS,nP),1:nEq,1:4) ) 

      
      myDGS % solution = ZERO

      myDGS % tendency = ZERO
 
      myDGS % boundarySolution = ZERO

      myDGS % boundaryFlux = ZERO


 END SUBROUTINE Build_DGSEMSolution_2D
!
!
!
 SUBROUTINE Trash_DGSEMSolution_2D( myDGS )
 ! S/R Trash_DGSEMSolution_2D
 !  
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS(DGSEMSolution_2D), INTENT(inout) :: myDGS

      DEALLOCATE( myDGS % Solution )
      DEALLOCATE( myDGS % tendency )
      DEALLOCATE( myDGS % boundarySolution ) 
      DEALLOCATE( myDGS % boundaryFlux ) 

 END SUBROUTINE Trash_DGSEMSolution_2D
!
!
!==================================================================================================!
!--------------------------------------- Accessors ------------------------------------------------!
!==================================================================================================!
!
!
! ------------------------------------ Number of Equations --------------------------------------- !
!
 SUBROUTINE GetNumberOfEquations_DGSEMSolution_2D( myDGS, nEq  )
 ! S/R GetNumberOfEquations_DGSEMSolution_2D
 !  
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS(DGSEMSolution_2D), INTENT(in) :: myDGS
   INTEGER, INTENT(out)                :: nEq 

      nEq = myDGS % nEq

 END SUBROUTINE GetNumberOfEquations_DGSEMSolution_2D
!
!
!
 SUBROUTINE SetNumberOfEquations_DGSEMSolution_2D( myDGS, nEq  )
 ! S/R SetNumberOfEquations_DGSEMSolution_2D
 !  
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS(DGSEMSolution_2D), INTENT(inout) :: myDGS
   INTEGER, INTENT(in)                    :: nEq ! 

      myDGS % nEq = nEq

 END SUBROUTINE SetNumberOfEquations_DGSEMSolution_2D
!
! ------------------------------------- Number of Nodes ------------------------------------------ !
!
 SUBROUTINE GetNumberOfNodes_DGSEMSolution_2D( myDGS, nS, nP  )
 ! S/R GetNumberOfNodes
 !  
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS(DGSEMSolution_2D), INTENT(in) :: myDGS
   INTEGER, INTENT(out)                :: nS, nP

      nS = myDGS % nS
      nP = myDGS % nP

 END SUBROUTINE GetNumberOfNodes_DGSEMSolution_2D
!
!
!
 SUBROUTINE SetNumberOfNodes_DGSEMSolution_2D( myDGS, nS, nP )
 ! S/R SetNumberOfNodes
 !  
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS(DGSEMSolution_2D), INTENT(inout) :: myDGS
   INTEGER, INTENT(in)                    :: nS, nP

      myDGS % nS = nS
      myDGS % nP = nP
      myDGS % nMax = max(nS,nP)
 END SUBROUTINE SetNumberOfNodes_DGSEMSolution_2D
!
! ---------------------------------------- Solution ---------------------------------------------- !
!
 SUBROUTINE GetSolution_DGSEMSolution_2D( myDGS, theSolution  )
 ! S/R GetSolution
 !  
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS(DGSEMSolution_2D), INTENT(in) :: myDGS
   REAL(prec), INTENT(out)             :: theSolution(0:myDGS % nS, 0:myDGS % nP, 1:myDGS % nEq)

      theSolution = myDGS % Solution

 END SUBROUTINE GetSolution_DGSEMSolution_2D
!
!
!
 SUBROUTINE SetSolution_DGSEMSolution_2D( myDGS, theSolution  )
 ! S/R SetSolution
 !  
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS(DGSEMSolution_2D), INTENT(inout) :: myDGS
   REAL(prec), INTENT(in)                 :: theSolution(0:myDGS % nS, 0:myDGS % nP, 1:myDGS % nEq)
   ! LOCAL

      myDGS % Solution = theSolution

 END SUBROUTINE SetSolution_DGSEMSolution_2D
!
!
!
 SUBROUTINE GetSolutionAtNode_DGSEMSolution_2D( myDGS, i, j, theSolution  )
 ! S/R GetSolution
 !  
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS(DGSEMSolution_2D), INTENT(in) :: myDGS
   INTEGER, INTENT(in)                 :: i, j
   REAL(prec), INTENT(out)             :: theSolution(1:myDGS % nEq)

      theSolution = myDGS % Solution(i,j,:)

 END SUBROUTINE GetSolutionAtNode_DGSEMSolution_2D
!
!
!
 SUBROUTINE SetSolutionAtNode_DGSEMSolution_2D( myDGS, i, j, theSolution  )
 ! S/R SetSolution
 !  
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS(DGSEMSolution_2D), INTENT(inout) :: myDGS
   INTEGER, INTENT(in)                    :: i, j
   REAL(prec), INTENT(in)                 :: theSolution(1:myDGS % nEq)
   ! LOCAL

      myDGS % Solution(i,j,:) = theSolution

 END SUBROUTINE SetSolutionAtNode_DGSEMSolution_2D
!
!
!
SUBROUTINE GetSolutionWithVarID_DGSEMSolution_2D( myDGS, varID, theSolution  )
 ! S/R GetSolution
 !  
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS(DGSEMSolution_2D), INTENT(in) :: myDGS
   INTEGER, INTENT(in)                 :: varID 
   REAL(prec), INTENT(out)             :: theSolution(0:myDGS % nS,0:myDGS % nP)

      theSolution = myDGS % Solution(:,:,varID)

 END SUBROUTINE GetSolutionWithVarID_DGSEMSolution_2D
!
!
!
 SUBROUTINE SetSolutionWithVarID_DGSEMSolution_2D( myDGS, varID, theSolution  )
 ! S/R SetSolution
 !  
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS(DGSEMSolution_2D), INTENT(inout) :: myDGS
   INTEGER, INTENT(in)                    :: varID 
   REAL(prec), INTENT(in)                 :: theSolution(0:myDGS % nS,0:myDGS % nP)

      myDGS % Solution(:,:,varID) = theSolution

 END SUBROUTINE SetSolutionWithVarID_DGSEMSolution_2D
!
!
!
 SUBROUTINE GetSolutionAtNodeWithVarID_DGSEMSolution_2D( myDGS, i, j, varID, theSolution  )
 ! S/R GetSolution
 !  
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS(DGSEMSolution_2D), INTENT(in) :: myDGS
   INTEGER, INTENT(in)                 :: i, j, varID 
   REAL(prec), INTENT(out)             :: theSolution

      theSolution = myDGS % Solution(i,j,varID)

 END SUBROUTINE GetSolutionAtNodeWithVarID_DGSEMSolution_2D
!
!
!
 SUBROUTINE SetSolutionAtNodeWithVarID_DGSEMSolution_2D( myDGS, i, j, varID, theSolution  )
 ! S/R SetSolution
 !  
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS(DGSEMSolution_2D), INTENT(inout) :: myDGS
   INTEGER, INTENT(in)                    :: i, j, varID 
   REAL(prec), INTENT(in)                 :: theSolution
   ! LOCAL

      myDGS % Solution(i,j,varID) = theSolution

 END SUBROUTINE SetSolutionAtNodeWithVarID_DGSEMSolution_2D
!
! ----------------------------------------- Tendency --------------------------------------------- !
!
 SUBROUTINE GetTendency_DGSEMSolution_2D( myDGS, theTend  )
 ! S/R GetTendency
 !  
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS(DGSEMSolution_2D), INTENT(in) :: myDGS
   REAL(prec), INTENT(out)             :: theTend(0:myDGS % nS,0:myDGS % nP, 1:myDGS % nEq)

      theTend = myDGS % tendency

 END SUBROUTINE GetTendency_DGSEMSolution_2D
!
!
!
 SUBROUTINE SetTendency_DGSEMSolution_2D( myDGS, theTend  )
 ! S/R SetTendency
 !  
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS(DGSEMSolution_2D), INTENT(inout) :: myDGS
   REAL(prec), INTENT(in)                 :: theTend(0:myDGS % nS,0:myDGS % nP, 1:myDGS % nEq)
   ! LOCAL

      myDGS % tendency = theTend

 END SUBROUTINE SetTendency_DGSEMSolution_2D
!
!
!
 SUBROUTINE GetTendencyAtNode_DGSEMSolution_2D( myDGS, i, j, theTend  )
 ! S/R GetTendency
 !  
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS(DGSEMSolution_2D), INTENT(in) :: myDGS
   INTEGER, INTENT(in)                 :: i, j 
   REAL(prec), INTENT(out)             :: theTend(1:myDGS % nEq)

      theTend = myDGS % tendency(i,j,:)

 END SUBROUTINE GetTendencyAtNode_DGSEMSolution_2D
!
!
!
 SUBROUTINE SetTendencyAtNode_DGSEMSolution_2D( myDGS, i, j, theTend  )
 ! S/R SetTendency
 !  
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS(DGSEMSolution_2D), INTENT(inout) :: myDGS
   INTEGER, INTENT(in)                    :: i, j 
   REAL(prec), INTENT(in)                 :: theTend(1:myDGS % nEq)
   ! LOCAL

      myDGS % tendency(i,j,:) = theTend

 END SUBROUTINE SetTendencyAtNode_DGSEMSolution_2D
!
!
!
SUBROUTINE GetTendencyWithVarID_DGSEMSolution_2D( myDGS, varID, theTend  )
 ! S/R GetTendency
 !  
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS(DGSEMSolution_2D), INTENT(in) :: myDGS
   INTEGER, INTENT(in)                 :: varID 
   REAL(prec), INTENT(out)             :: theTend(0:myDGS % nS,0:myDGS % nP)

      theTend = myDGS % tendency(:,:,varID)

 END SUBROUTINE GetTendencyWithVarID_DGSEMSolution_2D
!
!
!
 SUBROUTINE SetTendencyWithVarID_DGSEMSolution_2D( myDGS, varID, theTend  )
 ! S/R SetTendency
 !  
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS(DGSEMSolution_2D), INTENT(inout) :: myDGS
   INTEGER, INTENT(in)                    :: varID 
   REAL(prec), INTENT(in)                 :: theTend(0:myDGS % nS,0:myDGS % nP)

      myDGS % tendency(:,:,varID) = theTend

 END SUBROUTINE SetTendencyWithVarID_DGSEMSolution_2D
!
!
!
 SUBROUTINE GetTendencyAtNodeWithVarID_DGSEMSolution_2D( myDGS, i, j, varID, theTend  )
 ! S/R GetTendency
 !  
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS(DGSEMSolution_2D), INTENT(in) :: myDGS
   INTEGER, INTENT(in)                 :: i, j, varID 
   REAL(prec), INTENT(out)             :: theTend

      theTend = myDGS % tendency(i,j,varID)

 END SUBROUTINE GetTendencyAtNodeWithVarID_DGSEMSolution_2D
!
!
!
 SUBROUTINE SetTendencyAtNodeWithVarID_DGSEMSolution_2D( myDGS, i, j, varID, theTend  )
 ! S/R SetTendency
 !  
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS(DGSEMSolution_2D), INTENT(inout) :: myDGS
   INTEGER, INTENT(in)                    :: i, j, varID 
   REAL(prec), INTENT(in)                 :: theTend
   ! LOCAL

      myDGS % tendency(i,j,varID) = theTend

 END SUBROUTINE SetTendencyAtNodeWithVarID_DGSEMSolution_2D
!
! ------------------------------- Boundary Solution ---------------------------------------------- !
!
 SUBROUTINE GetBoundarySolution_DGSEMSolution_2D( myDGS, theSolution  )
 ! S/R GetBoundarySolution
 !  
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS(DGSEMSolution_2D), INTENT(in) :: myDGS
   REAL(prec), INTENT(out)             :: theSolution(0:myDGS % nMax,1:myDGS % nEq,1:4)

      theSolution = myDGS % boundarySolution

 END SUBROUTINE GetBoundarySolution_DGSEMSolution_2D
!
!
!
 SUBROUTINE SetBoundarySolution_DGSEMSolution_2D( myDGS, theSolution  )
 ! S/R SetBoundarySolution
 !  
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS(DGSEMSolution_2D), INTENT(inout) :: myDGS
   REAL(prec), INTENT(in)                 :: theSolution(0:myDGS % nMax,1:myDGS % nEq,1:4)
   ! LOCAL

      myDGS % boundarySolution = theSolution

 END SUBROUTINE SetBoundarySolution_DGSEMSolution_2D
!
!
!
 SUBROUTINE GetBoundarySolutionAtBoundary_DGSEMSolution_2D( myDGS, boundary, theSolution  )
 ! S/R GetBoundarySolution
 !  
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS(DGSEMSolution_2D), INTENT(in) :: myDGS
   INTEGER, INTENT(in)                 :: boundary
   REAL(prec), INTENT(out)             :: theSolution(0:myDGS % nMax,1:myDGS % nEq)

      theSolution = myDGS % boundarySolution(0:myDGS % nMax,1:myDGS % nEq,boundary)

 END SUBROUTINE GetBoundarySolutionAtBoundary_DGSEMSolution_2D
!
!
!
 SUBROUTINE SetBoundarySolutionAtBoundary_DGSEMSolution_2D( myDGS, boundary, theSolution  )
 ! S/R SetBoundarySolution
 !  
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS(DGSEMSolution_2D), INTENT(inout) :: myDGS
   INTEGER, INTENT(in)                    :: boundary
   REAL(prec), INTENT(in)                 :: theSolution(0:myDGS % nMax,1:myDGS % nEq)
   ! LOCAL

      myDGS % boundarySolution(0:myDGS % nMax,1:myDGS % nEq,boundary) = theSolution

 END SUBROUTINE SetBoundarySolutionAtBoundary_DGSEMSolution_2D
!
!
!
 SUBROUTINE GetBoundarySolutionAtBoundaryWithVarID_DGSEMSolution_2D( myDGS, boundary, varID, theSolution  )
 ! S/R GetBoundarySolution
 !  
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS(DGSEMSolution_2D), INTENT(in) :: myDGS
   INTEGER, INTENT(in)                 :: boundary, varID 
   REAL(prec), INTENT(out)             :: theSolution(0:myDGS % nMax)

      theSolution = myDGS % boundarySolution(0:myDGS % nMax,varID,boundary)

 END SUBROUTINE GetBoundarySolutionAtBoundaryWithVarID_DGSEMSolution_2D
!
!
!
 SUBROUTINE SetBoundarySolutionAtBoundaryWithVarID_DGSEMSolution_2D( myDGS, boundary, varID, theSolution  )
 ! S/R SetBoundarySolution
 !  
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS(DGSEMSolution_2D), INTENT(inout) :: myDGS
   INTEGER, INTENT(in)                    :: boundary, varID 
   REAL(prec), INTENT(in)                 :: theSolution(0:myDGS % nMax)
   ! LOCAL

      myDGS % boundarySolution(0:myDGS % nMax,varID, boundary) = theSolution

 END SUBROUTINE SetBoundarySolutionAtBoundaryWithVarID_DGSEMSolution_2D
!
! --------------------------------- Boundary Flux ------------------------------------------------ !
!
SUBROUTINE GetBoundaryFlux_DGSEMSolution_2D( myDGS, theFlux  )
 ! S/R GetBoundaryFlux
 !  
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS(DGSEMSolution_2D), INTENT(in) :: myDGS
   REAL(prec), INTENT(out)             :: theFlux(0:myDGS % nMax,1:myDGS % nEq,1:4)

      theFlux = myDGS % boundaryFlux

 END SUBROUTINE GetBoundaryFlux_DGSEMSolution_2D
!
!
!
 SUBROUTINE SetBoundaryFlux_DGSEMSolution_2D( myDGS, theFlux  )
 ! S/R SetBoundaryFlux
 !  
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS(DGSEMSolution_2D), INTENT(inout) :: myDGS
   REAL(prec), INTENT(in)                 :: theFlux(0:myDGS % nMax,1:myDGS % nEq,1:4)
   ! LOCAL

      myDGS % boundaryFlux = theFlux

 END SUBROUTINE SetBoundaryFlux_DGSEMSolution_2D
!
!
!
 SUBROUTINE GetBoundaryFluxAtBoundary_DGSEMSolution_2D( myDGS, boundary, theFlux  )
 ! S/R GetBoundaryFlux
 !  
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS(DGSEMSolution_2D), INTENT(in) :: myDGS
   INTEGER, INTENT(in)                 :: boundary
   REAL(prec), INTENT(out)             :: theFlux(0:myDGS % nMax,1:myDGS % nEq)

      theFlux = myDGS % boundaryFlux(0:myDGS % nMax,1:myDGS % nEq,boundary)

 END SUBROUTINE GetBoundaryFluxAtBoundary_DGSEMSolution_2D
!
!
!
 SUBROUTINE SetBoundaryFluxAtBoundary_DGSEMSolution_2D( myDGS, boundary, theFlux  )
 ! S/R SetBoundaryFlux
 !  
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS(DGSEMSolution_2D), INTENT(inout) :: myDGS
   INTEGER, INTENT(in)                    :: boundary
   REAL(prec), INTENT(in)                 :: theFlux(0:myDGS % nMax,1:myDGS % nEq)
   ! LOCAL

      myDGS % boundaryFlux(0:myDGS % nMax,1:myDGS % nEq,boundary) = theFlux

 END SUBROUTINE SetBoundaryFluxAtBoundary_DGSEMSolution_2D
!
!
!
  SUBROUTINE GetBoundaryFluxAtBoundaryNode_DGSEMSolution_2D( myDGS, boundary, i, theFlux  )
 ! S/R GetBoundaryFluxAtBoundaryNode
 !  
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS(DGSEMSolution_2D), INTENT(in) :: myDGS
   INTEGER, INTENT(in)                 :: boundary, i
   REAL(prec), INTENT(out)             :: theFlux(1:myDGS % nEq)

      theFlux = myDGS % boundaryFlux(i,1:myDGS % nEq,boundary)

 END SUBROUTINE GetBoundaryFluxAtBoundaryNode_DGSEMSolution_2D
!
!
!
 SUBROUTINE SetBoundaryFluxAtBoundaryNode_DGSEMSolution_2D( myDGS, boundary, i, theFlux  )
 ! S/R SetBoundaryFluxAtBoundaryNode
 !  
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS(DGSEMSolution_2D), INTENT(inout) :: myDGS
   INTEGER, INTENT(in)                    :: boundary, i
   REAL(prec), INTENT(in)                 :: theFlux(1:myDGS % nEq)
   ! LOCAL

      myDGS % boundaryFlux(i,1:myDGS % nEq,boundary) = theFlux

 END SUBROUTINE SetBoundaryFluxAtBoundaryNode_DGSEMSolution_2D
!
!
!
 SUBROUTINE GetBoundaryFluxAtBoundaryWithVarID_DGSEMSolution_2D( myDGS, boundary, varID, theFlux  )
 ! S/R GetBoundaryFlux
 !  
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS(DGSEMSolution_2D), INTENT(in) :: myDGS
   INTEGER, INTENT(in)                 :: boundary, varID 
   REAL(prec), INTENT(out)             :: theFlux(0:myDGS % nMax)

      theFlux = myDGS % boundaryFlux(0:myDGS % nMax,varID,boundary)

 END SUBROUTINE GetBoundaryFluxAtBoundaryWithVarID_DGSEMSolution_2D
!
!
!
 SUBROUTINE SetBoundaryFluxAtBoundaryWithVarID_DGSEMSolution_2D( myDGS, boundary, varID, theFlux  )
 ! S/R SetBoundaryFlux
 !  
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS(DGSEMSolution_2D), INTENT(inout) :: myDGS
   INTEGER, INTENT(in)                    :: boundary, varID 
   REAL(prec), INTENT(in)                 :: theFlux(0:myDGS % nMax)
   ! LOCAL

      myDGS % boundaryFlux(0:myDGS % nMax,varID, boundary) = theFlux

 END SUBROUTINE SetBoundaryFluxAtBoundaryWithVarID_DGSEMSolution_2D
!
!
!==================================================================================================!
!--------------------------------- Type Specific Routines -----------------------------------------!
!==================================================================================================!
!
!
 SUBROUTINE CalculateSolutionAtBoundaries_DGSEMSolution_2D( myDGS, dgStorage )
 ! S/R CalculateSolutionAtBoundaries
 !  
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS(DGSEMSolution_2D), INTENT(inout) :: myDGS
   TYPE(NodalStorage_2D), INTENT(in)      :: dgStorage
   ! LOCAL
   INTEGER :: iS, iP, iEq, nS, nP, nEq
   REAL(prec) :: lagSouth(0:myDGS % nP)
   REAL(prec) :: lagNorth(0:myDGS % nP)
   REAL(prec) :: lagWest(0:myDGS % nS)
   REAL(prec) :: lagEast(0:myDGS % nS)
   
       CALL myDGS % GetNumberOfEquations( nEq )
       CALL myDGS % GetNumberOfNodes( nS, nP )
       
       CALL dgStorage % GetSouthernInterpolants( lagSouth )
       CALL dgStorage % GetEasternInterpolants( lagEast )
       CALL dgStorage % GetNorthernInterpolants( lagNorth )
       CALL dgStorage % GetWesternInterpolants( lagWest )
       
       DO iP = 0, nP ! Loop over p-poinys
          DO iEq = 1, nEq ! Loop over the number of equations
          
             ! Setting the "west boundary"
              myDGS % boundarySolution(iP,iEq,4) = DOT_PRODUCT( lagWest, &
                                                        myDGS % solution(:,iP,iEq) )

            ! Setting the "east boundary"
              myDGS % boundarySolution(iP,iEq,2) = DOT_PRODUCT( lagEast, &
                                                        myDGS % solution(:,iP,iEq) )

          ENDDO ! iEq, loop over the number of equations
       ENDDO ! iP, loop over p-points
         

       DO iS = 0, nS ! Loop over s - points
          DO iEq = 1, nEq ! Loop over the number of equations

             ! Setting the "south boundary"
             myDGS % boundarySolution(iS,iEq,1)  = DOT_PRODUCT( lagSouth, &
                                                        myDGS % solution(iS,:,iEq) )

            ! Setting the "north boundary"
             myDGS % boundarySolution(iS,iEq,3)  = DOT_PRODUCT( lagNorth, &
                                                        myDGS % solution(iS,:,iEq) )

          ENDDO ! iEq, loop over the number of equations
       ENDDO ! iS, loop over s-points
   
 END SUBROUTINE CalculateSolutionAtBoundaries_DGSEMSolution_2D
!
!
!
END MODULE DGSEMSolutionStorageClass_2D
