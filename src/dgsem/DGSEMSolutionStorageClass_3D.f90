MODULE DGSEMSolutionStorageClass_3D
! DGSEMSolutionStorageClass_3D.f90
! 
! Module History 
! 
! o  (v1.0) March 2014
!    
! o  (v2.1) December 10, 2015
!
!  
! =======================================================================================
USE ConstantsDictionary
USE ModelPrecision


IMPLICIT NONE


    TYPE DGSEMSolution_3D
      INTEGER                 :: nEq, nS, nP, nQ, nMax
      REAL(prec), allocatable :: Solution(:,:,:,:)
      REAL(prec), allocatable :: tendency(:,:,:,:)
      REAL(prec), allocatable :: boundarySolution(:,:,:,:) ! Indexed over the boundaries. 
      REAL(prec), allocatable :: boundaryFlux(:,:,:,:)     ! Indexed over the boundaries

      !** Side ID's : 1 is south, 2 is east, 3 is north, 4 is west

      CONTAINS

      ! CONSTRUCTORS/DESTRUCTORS
      PROCEDURE :: Build => Build_DGSEMSolution_3D
      PROCEDURE :: Trash => Trash_DGSEMSolution_3D

      ! ACCESSORS
      PROCEDURE :: GetNumberOfEquations => GetNumberOfEquations_DGSEM_3D
      PROCEDURE :: SetNumberOfEquations => SetNumberOfEquations_DGSEM_3D
      
      PROCEDURE :: GetNumberOfNodes => GetNumberOfNodes_DGSEM_3D
      PROCEDURE :: SetNumberOfNodes => SetNumberOfNodes_DGSEM_3D

      PROCEDURE :: GetSolution => GetSolution_DGSEM_3D
      PROCEDURE :: SetSolution => SetSolution_DGSEM_3D
      PROCEDURE :: GetSolutionAtNode => GetSolutionAtNode_DGSEM_3D
      PROCEDURE :: SetSolutionAtNode => SetSolutionAtNode_DGSEM_3D
      PROCEDURE :: GetSolutionWithVarID => GetSolutionWithVarID_DGSEM_3D
      PROCEDURE :: SetSolutionWithVarID => SetSolutionWithVarID_DGSEM_3D
      PROCEDURE :: GetSolutionAtNodeWithVarID => GetSolutionAtNodeWithVarID_DGSEM_3D
      PROCEDURE :: SetSolutionAtNodeWithVarID => SetSolutionAtNodeWithVarID_DGSEM_3D
      
      PROCEDURE :: GetTendency => GetTendency_DGSEM_3D
      PROCEDURE :: SetTendency => SetTendency_DGSEM_3D
      PROCEDURE :: GetTendencyAtNode => GetTendencyAtNode_DGSEM_3D
      PROCEDURE :: SetTendencyAtNode => SetTendencyAtNode_DGSEM_3D
      PROCEDURE :: GetTendencyWithVarID => GetTendencyWithVarID_DGSEM_3D
      PROCEDURE :: SetTendencyWithVarID => SetTendencyWithVarID_DGSEM_3D
      PROCEDURE :: GetTendencyAtNodeWithVarID => GetTendencyAtNodeWithVarID_DGSEM_3D
      PROCEDURE :: SetTendencyAtNodeWithVarID => SetTendencyAtNodeWithVarID_DGSEM_3D
      
      PROCEDURE :: GetBoundarySolution => GetBoundarySolution_DGSEM_3D
      PROCEDURE :: SetBoundarySolution => SetBoundarySolution_DGSEM_3D
      PROCEDURE :: GetBoundarySolutionAtBoundary => GetBoundarySolutionAtBoundary_DGSEM_3D
      PROCEDURE :: SetBoundarySolutionAtBoundary => SetBoundarySolutionAtBoundary_DGSEM_3D
      PROCEDURE :: GetBoundarySolutionAtBoundaryWithVarID => GetBoundarySolutionAtBoundaryWithVarID_DGSEM_3D
      PROCEDURE :: SetBoundarySolutionAtBoundaryWithVarID => SetBoundarySolutionAtBoundaryWithVarID_DGSEM_3D

      PROCEDURE :: GetBoundaryFlux => GetBoundaryFlux_DGSEM_3D
      PROCEDURE :: SetBoundaryFlux => SetBoundaryFlux_DGSEM_3D
      PROCEDURE :: GetBoundaryFluxAtBoundary => GetBoundaryFluxAtBoundary_DGSEM_3D
      PROCEDURE :: SetBoundaryFluxAtBoundary => SetBoundaryFluxAtBoundary_DGSEM_3D
      PROCEDURE :: GetBoundaryFluxAtBoundaryWithVarID => GetBoundaryFluxAtBoundaryWithVarID_DGSEM_3D
      PROCEDURE :: SetBoundaryFluxAtBoundaryWithVarID => SetBoundaryFluxAtBoundaryWithVarID_DGSEM_3D 
      
    END TYPE DGSEMSolution_3D

 
 CONTAINS
!
!
!==================================================================================================!
!------------------------------- Manual Constructors/Destructors ----------------------------------!
!==================================================================================================!
!
!
 SUBROUTINE Build_DGSEMSolution_3D( myDGS, nS, nP, nQ, nEq )
 ! S/R Build_DGSEMSolution_3D
 !  
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS(DGSEMSolution_3D), INTENT(inout) :: myDGS
   INTEGER, INTENT(in)                    :: nS, nP, nQ, nEq
      
      ALLOCATE( myDGS % solution(0:nS,0:nP,0:nQ,1:nEq) )
      ALLOCATE( myDGS % tendency(0:nS,0:nP,0:nQ,1:nEq) )
      ALLOCATE( myDGS % boundarySolution(0:max(nS,nP,nQ),0:max(nS,nP,nQ),1:nEq,1:4) ) 
      ALLOCATE( myDGS % boundaryFlux(0:max(nS,nP,nQ),0:max(nS,nP,nQ),1:nEq,1:4) ) 

      
      myDGS % solution = ZERO

      myDGS % tendency = ZERO
 
      myDGS % boundarySolution = ZERO

      myDGS % boundaryFlux = ZERO

      CALL myDGS % SetNumberOfEquations( nEq )
      CALL myDGS % SetNumberOfNodes( nS, nP, nQ )

 END SUBROUTINE Build_DGSEMSolution_3D
!
!
!
 SUBROUTINE Trash_DGSEMSolution_3D( myDGS )
 ! S/R Trash_DGSEMSolution_3D
 !  
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS(DGSEMSolution_3D), INTENT(inout) :: myDGS

      DEALLOCATE( myDGS % Solution )
      DEALLOCATE( myDGS % tendency )
      DEALLOCATE( myDGS % boundarySolution ) 
      DEALLOCATE( myDGS % boundaryFlux ) 

 END SUBROUTINE Trash_DGSEMSolution_3D
!
!
!==================================================================================================!
!--------------------------------------- Accessors ------------------------------------------------!
!==================================================================================================!
!
!
! ------------------------------------ Number of Equations --------------------------------------- !
!
 SUBROUTINE GetNumberOfEquations_DGSEM_3D( myDGS, nEq  )
 ! S/R GetNumberOfEquations_DGSEM_3D
 !  
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS(DGSEMSolution_3D), INTENT(in) :: myDGS
   INTEGER, INTENT(out)                :: nEq 

      nEq = myDGS % nEq

 END SUBROUTINE GetNumberOfEquations_DGSEM_3D
!
!
!
 SUBROUTINE SetNumberOfEquations_DGSEM_3D( myDGS, nEq  )
 ! S/R SetNumberOfEquations_DGSEM_3D
 !  
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS(DGSEMSolution_3D), INTENT(inout) :: myDGS
   INTEGER, INTENT(in)                    :: nEq ! 

      myDGS % nEq = nEq

 END SUBROUTINE SetNumberOfEquations_DGSEM_3D
!
! ------------------------------------- Number of Nodes ------------------------------------------ !
!
 SUBROUTINE GetNumberOfNodes_DGSEM_3D( myDGS, nS, nP, nQ  )
 ! S/R GetNumberOfNodes
 !  
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS(DGSEMSolution_3D), INTENT(in) :: myDGS
   INTEGER, INTENT(out)                :: nS, nP, nQ

      nS = myDGS % nS
      nP = myDGS % nP
      nQ = myDGS % nQ
      
 END SUBROUTINE GetNumberOfNodes_DGSEM_3D
!
!
!
 SUBROUTINE SetNumberOfNodes_DGSEM_3D( myDGS, nS, nP, nQ )
 ! S/R SetNumberOfNodes
 !  
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS(DGSEMSolution_3D), INTENT(inout) :: myDGS
   INTEGER, INTENT(in)                    :: nS, nP, nQ

      myDGS % nS = nS
      myDGS % nP = nP
      myDGS % nQ = nQ

      myDGS % nMax = max( nS, nP, nQ )

 END SUBROUTINE SetNumberOfNodes_DGSEM_3D
!
! ---------------------------------------- Solution ---------------------------------------------- !
!
 SUBROUTINE GetSolution_DGSEM_3D( myDGS, theSolution  )
 ! S/R GetSolution
 !  
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS(DGSEMSolution_3D), INTENT(in) :: myDGS
   REAL(prec), INTENT(out)             :: theSolution(0:myDGS % nS, 0:myDGS % nP, 0:myDGS % nQ, 1:myDGS % nEq)

      theSolution = myDGS % Solution

 END SUBROUTINE GetSolution_DGSEM_3D
!
!
!
 SUBROUTINE SetSolution_DGSEM_3D( myDGS, theSolution  )
 ! S/R SetSolution
 !  
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS(DGSEMSolution_3D), INTENT(inout) :: myDGS
   REAL(prec), INTENT(in)                 :: theSolution(0:myDGS % nS, 0:myDGS % nP, 0:myDGS % nQ, 1:myDGS % nEq)
   ! LOCAL

      myDGS % Solution = theSolution

 END SUBROUTINE SetSolution_DGSEM_3D
!
!
!
 SUBROUTINE GetSolutionAtNode_DGSEM_3D( myDGS, i, j, k, theSolution  )
 ! S/R GetSolution
 !  
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS(DGSEMSolution_3D), INTENT(in) :: myDGS
   INTEGER, INTENT(in)                 :: i, j, k
   REAL(prec), INTENT(out)             :: theSolution(1:myDGS % nEq)

      theSolution = myDGS % Solution(i,j,k,:)

 END SUBROUTINE GetSolutionAtNode_DGSEM_3D
!
!
!
 SUBROUTINE SetSolutionAtNode_DGSEM_3D( myDGS, i, j, k, theSolution  )
 ! S/R SetSolution
 !  
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS(DGSEMSolution_3D), INTENT(inout) :: myDGS
   INTEGER, INTENT(in)                    :: i, j, k
   REAL(prec), INTENT(in)                 :: theSolution(1:myDGS % nEq)
   ! LOCAL

      myDGS % Solution(i,j,k,:) = theSolution

 END SUBROUTINE SetSolutionAtNode_DGSEM_3D
!
!
!
SUBROUTINE GetSolutionWithVarID_DGSEM_3D( myDGS, varID, theSolution  )
 ! S/R GetSolution
 !  
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS(DGSEMSolution_3D), INTENT(in) :: myDGS
   INTEGER, INTENT(in)                 :: varID 
   REAL(prec), INTENT(out)             :: theSolution(0:myDGS % nS,0:myDGS % nP,0:myDGS % nQ)

      theSolution = myDGS % Solution(:,:,:,varID)

 END SUBROUTINE GetSolutionWithVarID_DGSEM_3D
!
!
!
 SUBROUTINE SetSolutionWithVarID_DGSEM_3D( myDGS, varID, theSolution  )
 ! S/R SetSolution
 !  
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS(DGSEMSolution_3D), INTENT(inout) :: myDGS
   INTEGER, INTENT(in)                    :: varID 
   REAL(prec), INTENT(in)                 :: theSolution(0:myDGS % nS,0:myDGS % nP,0:myDGS % nQ)

      myDGS % Solution(:,:,:,varID) = theSolution

 END SUBROUTINE SetSolutionWithVarID_DGSEM_3D
!
!
!
 SUBROUTINE GetSolutionAtNodeWithVarID_DGSEM_3D( myDGS, i, j, k, varID, theSolution  )
 ! S/R GetSolution
 !  
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS(DGSEMSolution_3D), INTENT(in) :: myDGS
   INTEGER, INTENT(in)                 :: i, j, k, varID 
   REAL(prec), INTENT(out)             :: theSolution

      theSolution = myDGS % Solution(i,j,k,varID)

 END SUBROUTINE GetSolutionAtNodeWithVarID_DGSEM_3D
!
!
!
 SUBROUTINE SetSolutionAtNodeWithVarID_DGSEM_3D( myDGS, i, j, k, varID, theSolution  )
 ! S/R SetSolution
 !  
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS(DGSEMSolution_3D), INTENT(inout) :: myDGS
   INTEGER, INTENT(in)                    :: i, j, k, varID 
   REAL(prec), INTENT(in)                 :: theSolution
   ! LOCAL

      myDGS % Solution(i,j,k,varID) = theSolution

 END SUBROUTINE SetSolutionAtNodeWithVarID_DGSEM_3D
!
! ----------------------------------------- Tendency --------------------------------------------- !
!
 SUBROUTINE GetTendency_DGSEM_3D( myDGS, theTend  )
 ! S/R GetTendency
 !  
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS(DGSEMSolution_3D), INTENT(in) :: myDGS
   REAL(prec), INTENT(out)             :: theTend(0:myDGS % nS, 0:myDGS % nP, 0:myDGS % nQ, 1:myDGS % nEq)

      theTend = myDGS % tendency

 END SUBROUTINE GetTendency_DGSEM_3D
!
!
!
 SUBROUTINE SetTendency_DGSEM_3D( myDGS, theTend  )
 ! S/R SetTendency
 !  
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS(DGSEMSolution_3D), INTENT(inout) :: myDGS
   REAL(prec), INTENT(in)                 :: theTend(0:myDGS % nS, 0:myDGS % nP, 0:myDGS % nQ, 1:myDGS % nEq)
   ! LOCAL

      myDGS % tendency = theTend

 END SUBROUTINE SetTendency_DGSEM_3D
!
!
!
 SUBROUTINE GetTendencyAtNode_DGSEM_3D( myDGS, i, j, k, theTend  )
 ! S/R GetTendency
 !  
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS(DGSEMSolution_3D), INTENT(in) :: myDGS
   INTEGER, INTENT(in)                 :: i, j, k 
   REAL(prec), INTENT(out)             :: theTend(1:myDGS % nEq)

      theTend = myDGS % tendency(i,j,k,:)

 END SUBROUTINE GetTendencyAtNode_DGSEM_3D
!
!
!
 SUBROUTINE SetTendencyAtNode_DGSEM_3D( myDGS, i, j, k, theTend  )
 ! S/R SetTendency
 !  
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS(DGSEMSolution_3D), INTENT(inout) :: myDGS
   INTEGER, INTENT(in)                    :: i, j, k 
   REAL(prec), INTENT(in)                 :: theTend(1:myDGS % nEq)
   ! LOCAL

      myDGS % tendency(i,j,k,:) = theTend

 END SUBROUTINE SetTendencyAtNode_DGSEM_3D
!
!
!
SUBROUTINE GetTendencyWithVarID_DGSEM_3D( myDGS, varID, theTend  )
 ! S/R GetTendency
 !  
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS(DGSEMSolution_3D), INTENT(in) :: myDGS
   INTEGER, INTENT(in)                 :: varID 
   REAL(prec), INTENT(out)             :: theTend(0:myDGS % nS,0:myDGS % nP,0:myDGS % nQ)

      theTend = myDGS % tendency(:,:,:,varID)

 END SUBROUTINE GetTendencyWithVarID_DGSEM_3D
!
!
!
 SUBROUTINE SetTendencyWithVarID_DGSEM_3D( myDGS, varID, theTend  )
 ! S/R SetTendency
 !  
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS(DGSEMSolution_3D), INTENT(inout) :: myDGS
   INTEGER, INTENT(in)                    :: varID 
   REAL(prec), INTENT(in)                 :: theTend(0:myDGS % nS,0:myDGS % nP,0:myDGS % nQ)

      myDGS % tendency(:,:,:,varID) = theTend

 END SUBROUTINE SetTendencyWithVarID_DGSEM_3D
!
!
!
 SUBROUTINE GetTendencyAtNodeWithVarID_DGSEM_3D( myDGS, i, j, k, varID, theTend  )
 ! S/R GetTendency
 !  
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS(DGSEMSolution_3D), INTENT(in) :: myDGS
   INTEGER, INTENT(in)                 :: i, j, k, varID 
   REAL(prec), INTENT(out)             :: theTend

      theTend = myDGS % tendency(i,j,k,varID)

 END SUBROUTINE GetTendencyAtNodeWithVarID_DGSEM_3D
!
!
!
 SUBROUTINE SetTendencyAtNodeWithVarID_DGSEM_3D( myDGS, i, j, k, varID, theTend  )
 ! S/R SetTendency
 !  
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS(DGSEMSolution_3D), INTENT(inout) :: myDGS
   INTEGER, INTENT(in)                    :: i, j, k, varID 
   REAL(prec), INTENT(in)                 :: theTend
   ! LOCAL

      myDGS % tendency(i,j,k,varID) = theTend

 END SUBROUTINE SetTendencyAtNodeWithVarID_DGSEM_3D
!
! ------------------------------- Boundary Solution ---------------------------------------------- !
!
 SUBROUTINE GetBoundarySolution_DGSEM_3D( myDGS, theSolution  )
 ! S/R GetBoundarySolution
 !  
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS(DGSEMSolution_3D), INTENT(in) :: myDGS
   REAL(prec), INTENT(out)             :: theSolution(0:myDGS % nMax,0:myDGS % nMax,1:myDGS % nEq,1:6)

      theSolution = myDGS % boundarySolution

 END SUBROUTINE GetBoundarySolution_DGSEM_3D
!
!
!
 SUBROUTINE SetBoundarySolution_DGSEM_3D( myDGS, theSolution  )
 ! S/R SetBoundarySolution
 !  
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS(DGSEMSolution_3D), INTENT(inout) :: myDGS
   REAL(prec), INTENT(in)                 :: theSolution(0:myDGS % nMax,0:myDGS % nMax,1:myDGS % nEq,1:6)
   ! LOCAL

      myDGS % boundarySolution = theSolution

 END SUBROUTINE SetBoundarySolution_DGSEM_3D
!
!
!
 SUBROUTINE GetBoundarySolutionAtBoundary_DGSEM_3D( myDGS, boundary, theSolution  )
 ! S/R GetBoundarySolution
 !  
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS(DGSEMSolution_3D), INTENT(in) :: myDGS
   INTEGER, INTENT(in)                 :: boundary
   REAL(prec), INTENT(out)             :: theSolution(0:myDGS % nMax,0:myDGS % nMax,1:myDGS % nEq)

      theSolution = myDGS % boundarySolution(0:myDGS % nMax,0:myDGS % nMax,1:myDGS % nEq,boundary)

 END SUBROUTINE GetBoundarySolutionAtBoundary_DGSEM_3D
!
!
!
 SUBROUTINE SetBoundarySolutionAtBoundary_DGSEM_3D( myDGS, boundary, theSolution  )
 ! S/R SetBoundarySolution
 !  
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS(DGSEMSolution_3D), INTENT(inout) :: myDGS
   INTEGER, INTENT(in)                    :: boundary
   REAL(prec), INTENT(in)                 :: theSolution(0:myDGS % nMax,0:myDGS % nMax,1:myDGS % nEq)
   ! LOCAL

      myDGS % boundarySolution(0:myDGS % nMax,0:myDGS % nMax,1:myDGS % nEq,boundary) = theSolution

 END SUBROUTINE SetBoundarySolutionAtBoundary_DGSEM_3D
!
!
!
 SUBROUTINE GetBoundarySolutionAtBoundaryWithVarID_DGSEM_3D( myDGS, boundary, varID, theSolution  )
 ! S/R GetBoundarySolution
 !  
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS(DGSEMSolution_3D), INTENT(in) :: myDGS
   INTEGER, INTENT(in)                 :: boundary, varID 
   REAL(prec), INTENT(out)             :: theSolution(0:myDGS % nMax,0:myDGS % nMax)

      theSolution = myDGS % boundarySolution(0:myDGS % nMax, 0:myDGS % nMax, varID, boundary)

 END SUBROUTINE GetBoundarySolutionAtBoundaryWithVarID_DGSEM_3D
!
!
!
 SUBROUTINE SetBoundarySolutionAtBoundaryWithVarID_DGSEM_3D( myDGS, boundary, varID, theSolution  )
 ! S/R SetBoundarySolution
 !  
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS(DGSEMSolution_3D), INTENT(inout) :: myDGS
   INTEGER, INTENT(in)                    :: boundary, varID 
   REAL(prec), INTENT(in)                 :: theSolution(0:myDGS % nMax,0:myDGS % nMax)
   ! LOCAL

      myDGS % boundarySolution(0:myDGS % nMax, 0:myDGS % nMax, varID, boundary) = theSolution

 END SUBROUTINE SetBoundarySolutionAtBoundaryWithVarID_DGSEM_3D
!
! --------------------------------- Boundary Flux ------------------------------------------------ !
!
SUBROUTINE GetBoundaryFlux_DGSEM_3D( myDGS, theFlux  )
 ! S/R GetBoundaryFlux
 !  
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS(DGSEMSolution_3D), INTENT(in) :: myDGS
   REAL(prec), INTENT(out)             :: theFlux(0:myDGS % nMax,0:myDGS % nMax,1:myDGS % nEq,1:6)

      theFlux = myDGS % boundaryFlux

 END SUBROUTINE GetBoundaryFlux_DGSEM_3D
!
!
!
 SUBROUTINE SetBoundaryFlux_DGSEM_3D( myDGS, theFlux  )
 ! S/R SetBoundaryFlux
 !  
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS(DGSEMSolution_3D), INTENT(inout) :: myDGS
   REAL(prec), INTENT(in)                 :: theFlux(0:myDGS % nMax,0:myDGS % nMax,1:myDGS % nEq,1:6)
   ! LOCAL

      myDGS % boundaryFlux = theFlux

 END SUBROUTINE SetBoundaryFlux_DGSEM_3D
!
!
!
 SUBROUTINE GetBoundaryFluxAtBoundary_DGSEM_3D( myDGS, boundary, theFlux  )
 ! S/R GetBoundaryFlux
 !  
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS(DGSEMSolution_3D), INTENT(in) :: myDGS
   INTEGER, INTENT(in)                 :: boundary
   REAL(prec), INTENT(out)             :: theFlux(0:myDGS % nMax,0:myDGS % nMax,1:myDGS % nEq)

      theFlux = myDGS % boundaryFlux(0:myDGS % nMax, 0:myDGS % nMax, 1:myDGS % nEq, boundary)

 END SUBROUTINE GetBoundaryFluxAtBoundary_DGSEM_3D
!
!
!
 SUBROUTINE SetBoundaryFluxAtBoundary_DGSEM_3D( myDGS, boundary, theFlux  )
 ! S/R SetBoundaryFlux
 !  
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS(DGSEMSolution_3D), INTENT(inout) :: myDGS
   INTEGER, INTENT(in)                    :: boundary
   REAL(prec), INTENT(in)                 :: theFlux(0:myDGS % nMax,0:myDGS % nMax,1:myDGS % nEq)
   ! LOCAL

      myDGS % boundaryFlux(0:myDGS % nMax, 0:myDGS % nMax, 1:myDGS % nEq, boundary) = theFlux

 END SUBROUTINE SetBoundaryFluxAtBoundary_DGSEM_3D
!
!
!
 SUBROUTINE GetBoundaryFluxAtBoundaryWithVarID_DGSEM_3D( myDGS, boundary, varID, theFlux  )
 ! S/R GetBoundaryFlux
 !  
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS(DGSEMSolution_3D), INTENT(in) :: myDGS
   INTEGER, INTENT(in)                 :: boundary, varID 
   REAL(prec), INTENT(out)             :: theFlux(0:myDGS % nMax,0:myDGS % nMax)

      theFlux = myDGS % boundaryFlux(0:myDGS % nMax,0:myDGS % nMax, varID, boundary)

 END SUBROUTINE GetBoundaryFluxAtBoundaryWithVarID_DGSEM_3D
!
!
!
 SUBROUTINE SetBoundaryFluxAtBoundaryWithVarID_DGSEM_3D( myDGS, boundary, varID, theFlux  )
 ! S/R SetBoundaryFlux
 !  
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS(DGSEMSolution_3D), INTENT(inout) :: myDGS
   INTEGER, INTENT(in)                    :: boundary, varID 
   REAL(prec), INTENT(in)                 :: theFlux(0:myDGS % nMax,0:myDGS % nMax)
   ! LOCAL

      myDGS % boundaryFlux(0:myDGS % nMax, 0:myDGS % nMax, varID, boundary) = theFlux

 END SUBROUTINE SetBoundaryFluxAtBoundaryWithVarID_DGSEM_3D
!
!
!

END MODULE DGSEMSolutionStorageClass_3D
