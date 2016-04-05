MODULE CGSEMSolutionStorageClass_2D
! CGSEMSolutionStorageClass_2D.f90
! 
! Module History 
! 
! o  (v1.0) March 2014
!    
! o  (v2.1) December 10, 2015
!
!  
!  The CGSEM Solution Storage class is different from the CGSEM class in that it stores the solution
!  and the RHS (source) globally within the data-structure. Additionally, it includes a global
!  boundary flux that is indexed over the number of boundary edges (obtained from the BoundaryEdges
!  structure). The boundary flux is used to implement inhomogeneous Neumann boundary conditions.
!
! ================================================================================================ !
! src/common/
USE ConstantsDictionary
USE ModelPrecision
! src/nodal
USE NodalStorage_2D_Class

IMPLICIT NONE


    TYPE CGSEMSolution_2D
      INTEGER                 :: nBoundEdges, nEl, nS
      REAL(prec), allocatable :: solution(:,:,:)
      REAL(prec), allocatable :: source(:,:,:)
      REAL(prec), allocatable :: boundaryFlux(:,:)  ! Global boundary flux for inhomogeneous 
                                                    ! neumann boundary conditions. (0:nS,1:nBoundEdges)

 
      CONTAINS

      ! CONSTRUCTORS/DESTRUCTORS
      PROCEDURE :: Build => Build_CGSEMSolution_2D
      PROCEDURE :: Trash => Trash_CGSEMSolution_2D

      ! ACCESSORS
      PROCEDURE :: GetNumberOfNodes => GetNumberOfNodes_CGSEMSolution_2D
      PROCEDURE :: SetNumberOfNodes => SetNumberOfNodes_CGSEMSolution_2D

      PROCEDURE :: GetSolution => GetSolution_CGSEMSolution_2D
      PROCEDURE :: SetSolution => SetSolution_CGSEMSolution_2D
      PROCEDURE :: GetSolutionAtNode => GetSolutionAtNode_CGSEMSolution_2D
      PROCEDURE :: SetSolutionAtNode => SetSolutionAtNode_CGSEMSolution_2D
      
      PROCEDURE :: GetBoundaryFlux => GetBoundaryFlux_CGSEMSolution_2D
      PROCEDURE :: SetBoundaryFlux => SetBoundaryFlux_CGSEMSolution_2D
      PROCEDURE :: GetBoundaryFluxAtBoundary => GetBoundaryFluxAtBoundary_CGSEMSolution_2D
      PROCEDURE :: SetBoundaryFluxAtBoundary => SetBoundaryFluxAtBoundary_CGSEMSolution_2D
      PROCEDURE :: GetBoundaryFluxAtBoundaryNode => GetBoundaryFluxAtBoundaryNode_CGSEMSolution_2D
      PROCEDURE :: SetBoundaryFluxAtBoundaryNode => SetBoundaryFluxAtBoundaryNode_CGSEMSolution_2D
      PROCEDURE :: GetBoundaryFluxAtBoundaryWithVarID => GetBoundaryFluxAtBoundaryWithVarID_CGSEMSolution_2D
      PROCEDURE :: SetBoundaryFluxAtBoundaryWithVarID => SetBoundaryFluxAtBoundaryWithVarID_CGSEMSolution_2D 
      
    END TYPE CGSEMSolution_2D

 
 CONTAINS
!
!
!==================================================================================================!
!------------------------------- Manual Constructors/Destructors ----------------------------------!
!==================================================================================================!
!
!
 SUBROUTINE Build_CGSEMSolution_2D( myCGS, nS, nEl, nBoundEdges )
 ! S/R Build_CGSEMSolution_2D
 !  
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS(CGSEMSolution_2D), INTENT(inout) :: myCGS
   INTEGER, INTENT(in)                    :: nS, nEl, nBoundEdges
   
      CALL myCGS % SetNumberOfNodes( nS )
      CALL myCGS % SetNumberOfElements( nEl )
      CALL myCGS % SetNumberOfBoundaryEdges( nBoundEdges )
      
      ALLOCATE( myCGS % solution(0:nS,0:nS,1:nEl) )
      ALLOCATE( myCGS % source(0:nS,0:nS,1:nEl) )
      ALLOCATE( myCGS % boundaryFlux(0:nS,1:nBoundEdges) ) 

      
      myCGS % solution = ZERO
      myCGS % source = ZERO
      myCGS % boundaryFlux = ZERO


 END SUBROUTINE Build_CGSEMSolution_2D
!
!
!
 SUBROUTINE Trash_CGSEMSolution_2D( myCGS )
 ! S/R Trash_CGSEMSolution_2D
 !  
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS(CGSEMSolution_2D), INTENT(inout) :: myCGS

      DEALLOCATE( myCGS % solution )
      DEALLOCATE( myCGS % source ) 
      DEALLOCATE( myCGS % boundaryFlux ) 

 END SUBROUTINE Trash_CGSEMSolution_2D
!
!
!==================================================================================================!
!--------------------------------------- Accessors ------------------------------------------------!
!==================================================================================================!
!
!
! ------------------------------------ Number of Equations --------------------------------------- !
!
 SUBROUTINE GetNumberOfEquations_CGSEMSolution_2D( myCGS, nEq  )
 ! S/R GetNumberOfEquations_CGSEMSolution_2D
 !  
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS(CGSEMSolution_2D), INTENT(in) :: myCGS
   INTEGER, INTENT(out)                :: nEq 

      nEq = myCGS % nEq

 END SUBROUTINE GetNumberOfEquations_CGSEMSolution_2D
!
!
!
 SUBROUTINE SetNumberOfEquations_CGSEMSolution_2D( myCGS, nEq  )
 ! S/R SetNumberOfEquations_CGSEMSolution_2D
 !  
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS(CGSEMSolution_2D), INTENT(inout) :: myCGS
   INTEGER, INTENT(in)                    :: nEq ! 

      myCGS % nEq = nEq

 END SUBROUTINE SetNumberOfEquations_CGSEMSolution_2D
!
! ------------------------------------- Number of Nodes ------------------------------------------ !
!
 SUBROUTINE GetNumberOfNodes_CGSEMSolution_2D( myCGS, nS, nP  )
 ! S/R GetNumberOfNodes
 !  
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS(CGSEMSolution_2D), INTENT(in) :: myCGS
   INTEGER, INTENT(out)                :: nS, nP

      nS = myCGS % nS
      nP = myCGS % nP

 END SUBROUTINE GetNumberOfNodes_CGSEMSolution_2D
!
!
!
 SUBROUTINE SetNumberOfNodes_CGSEMSolution_2D( myCGS, nS, nP )
 ! S/R SetNumberOfNodes
 !  
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS(CGSEMSolution_2D), INTENT(inout) :: myCGS
   INTEGER, INTENT(in)                    :: nS, nP

      myCGS % nS = nS
      myCGS % nP = nP
      myCGS % nMax = max(nS,nP)
 END SUBROUTINE SetNumberOfNodes_CGSEMSolution_2D
!
! ---------------------------------------- Solution ---------------------------------------------- !
!
 SUBROUTINE GetSolution_CGSEMSolution_2D( myCGS, theSolution  )
 ! S/R GetSolution
 !  
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS(CGSEMSolution_2D), INTENT(in) :: myCGS
   REAL(prec), INTENT(out)             :: theSolution(0:myCGS % nS, 0:myCGS % nP, 1:myCGS % nEq)

      theSolution = myCGS % Solution

 END SUBROUTINE GetSolution_CGSEMSolution_2D
!
!
!
 SUBROUTINE SetSolution_CGSEMSolution_2D( myCGS, theSolution  )
 ! S/R SetSolution
 !  
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS(CGSEMSolution_2D), INTENT(inout) :: myCGS
   REAL(prec), INTENT(in)                 :: theSolution(0:myCGS % nS, 0:myCGS % nP, 1:myCGS % nEq)
   ! LOCAL

      myCGS % Solution = theSolution

 END SUBROUTINE SetSolution_CGSEMSolution_2D
!
!
!
 SUBROUTINE GetSolutionAtNode_CGSEMSolution_2D( myCGS, i, j, theSolution  )
 ! S/R GetSolution
 !  
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS(CGSEMSolution_2D), INTENT(in) :: myCGS
   INTEGER, INTENT(in)                 :: i, j
   REAL(prec), INTENT(out)             :: theSolution(1:myCGS % nEq)

      theSolution = myCGS % Solution(i,j,:)

 END SUBROUTINE GetSolutionAtNode_CGSEMSolution_2D
!
!
!
 SUBROUTINE SetSolutionAtNode_CGSEMSolution_2D( myCGS, i, j, theSolution  )
 ! S/R SetSolution
 !  
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS(CGSEMSolution_2D), INTENT(inout) :: myCGS
   INTEGER, INTENT(in)                    :: i, j
   REAL(prec), INTENT(in)                 :: theSolution(1:myCGS % nEq)
   ! LOCAL

      myCGS % Solution(i,j,:) = theSolution

 END SUBROUTINE SetSolutionAtNode_CGSEMSolution_2D
!
!
!
SUBROUTINE GetSolutionWithVarID_CGSEMSolution_2D( myCGS, varID, theSolution  )
 ! S/R GetSolution
 !  
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS(CGSEMSolution_2D), INTENT(in) :: myCGS
   INTEGER, INTENT(in)                 :: varID 
   REAL(prec), INTENT(out)             :: theSolution(0:myCGS % nS,0:myCGS % nP)

      theSolution = myCGS % Solution(:,:,varID)

 END SUBROUTINE GetSolutionWithVarID_CGSEMSolution_2D
!
!
!
 SUBROUTINE SetSolutionWithVarID_CGSEMSolution_2D( myCGS, varID, theSolution  )
 ! S/R SetSolution
 !  
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS(CGSEMSolution_2D), INTENT(inout) :: myCGS
   INTEGER, INTENT(in)                    :: varID 
   REAL(prec), INTENT(in)                 :: theSolution(0:myCGS % nS,0:myCGS % nP)

      myCGS % Solution(:,:,varID) = theSolution

 END SUBROUTINE SetSolutionWithVarID_CGSEMSolution_2D
!
!
!
 SUBROUTINE GetSolutionAtNodeWithVarID_CGSEMSolution_2D( myCGS, i, j, varID, theSolution  )
 ! S/R GetSolution
 !  
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS(CGSEMSolution_2D), INTENT(in) :: myCGS
   INTEGER, INTENT(in)                 :: i, j, varID 
   REAL(prec), INTENT(out)             :: theSolution

      theSolution = myCGS % Solution(i,j,varID)

 END SUBROUTINE GetSolutionAtNodeWithVarID_CGSEMSolution_2D
!
!
!
 SUBROUTINE SetSolutionAtNodeWithVarID_CGSEMSolution_2D( myCGS, i, j, varID, theSolution  )
 ! S/R SetSolution
 !  
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS(CGSEMSolution_2D), INTENT(inout) :: myCGS
   INTEGER, INTENT(in)                    :: i, j, varID 
   REAL(prec), INTENT(in)                 :: theSolution
   ! LOCAL

      myCGS % Solution(i,j,varID) = theSolution

 END SUBROUTINE SetSolutionAtNodeWithVarID_CGSEMSolution_2D
!
! ----------------------------------------- Tendency --------------------------------------------- !
!
 SUBROUTINE GetTendency_CGSEMSolution_2D( myCGS, theTend  )
 ! S/R GetTendency
 !  
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS(CGSEMSolution_2D), INTENT(in) :: myCGS
   REAL(prec), INTENT(out)             :: theTend(0:myCGS % nS,0:myCGS % nP, 1:myCGS % nEq)

      theTend = myCGS % tendency

 END SUBROUTINE GetTendency_CGSEMSolution_2D
!
!
!
 SUBROUTINE SetTendency_CGSEMSolution_2D( myCGS, theTend  )
 ! S/R SetTendency
 !  
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS(CGSEMSolution_2D), INTENT(inout) :: myCGS
   REAL(prec), INTENT(in)                 :: theTend(0:myCGS % nS,0:myCGS % nP, 1:myCGS % nEq)
   ! LOCAL

      myCGS % tendency = theTend

 END SUBROUTINE SetTendency_CGSEMSolution_2D
!
!
!
 SUBROUTINE GetTendencyAtNode_CGSEMSolution_2D( myCGS, i, j, theTend  )
 ! S/R GetTendency
 !  
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS(CGSEMSolution_2D), INTENT(in) :: myCGS
   INTEGER, INTENT(in)                 :: i, j 
   REAL(prec), INTENT(out)             :: theTend(1:myCGS % nEq)

      theTend = myCGS % tendency(i,j,:)

 END SUBROUTINE GetTendencyAtNode_CGSEMSolution_2D
!
!
!
 SUBROUTINE SetTendencyAtNode_CGSEMSolution_2D( myCGS, i, j, theTend  )
 ! S/R SetTendency
 !  
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS(CGSEMSolution_2D), INTENT(inout) :: myCGS
   INTEGER, INTENT(in)                    :: i, j 
   REAL(prec), INTENT(in)                 :: theTend(1:myCGS % nEq)
   ! LOCAL

      myCGS % tendency(i,j,:) = theTend

 END SUBROUTINE SetTendencyAtNode_CGSEMSolution_2D
!
!
!
SUBROUTINE GetTendencyWithVarID_CGSEMSolution_2D( myCGS, varID, theTend  )
 ! S/R GetTendency
 !  
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS(CGSEMSolution_2D), INTENT(in) :: myCGS
   INTEGER, INTENT(in)                 :: varID 
   REAL(prec), INTENT(out)             :: theTend(0:myCGS % nS,0:myCGS % nP)

      theTend = myCGS % tendency(:,:,varID)

 END SUBROUTINE GetTendencyWithVarID_CGSEMSolution_2D
!
!
!
 SUBROUTINE SetTendencyWithVarID_CGSEMSolution_2D( myCGS, varID, theTend  )
 ! S/R SetTendency
 !  
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS(CGSEMSolution_2D), INTENT(inout) :: myCGS
   INTEGER, INTENT(in)                    :: varID 
   REAL(prec), INTENT(in)                 :: theTend(0:myCGS % nS,0:myCGS % nP)

      myCGS % tendency(:,:,varID) = theTend

 END SUBROUTINE SetTendencyWithVarID_CGSEMSolution_2D
!
!
!
 SUBROUTINE GetTendencyAtNodeWithVarID_CGSEMSolution_2D( myCGS, i, j, varID, theTend  )
 ! S/R GetTendency
 !  
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS(CGSEMSolution_2D), INTENT(in) :: myCGS
   INTEGER, INTENT(in)                 :: i, j, varID 
   REAL(prec), INTENT(out)             :: theTend

      theTend = myCGS % tendency(i,j,varID)

 END SUBROUTINE GetTendencyAtNodeWithVarID_CGSEMSolution_2D
!
!
!
 SUBROUTINE SetTendencyAtNodeWithVarID_CGSEMSolution_2D( myCGS, i, j, varID, theTend  )
 ! S/R SetTendency
 !  
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS(CGSEMSolution_2D), INTENT(inout) :: myCGS
   INTEGER, INTENT(in)                    :: i, j, varID 
   REAL(prec), INTENT(in)                 :: theTend
   ! LOCAL

      myCGS % tendency(i,j,varID) = theTend

 END SUBROUTINE SetTendencyAtNodeWithVarID_CGSEMSolution_2D
!
! ------------------------------- Boundary Solution ---------------------------------------------- !
!
 SUBROUTINE GetBoundarySolution_CGSEMSolution_2D( myCGS, theSolution  )
 ! S/R GetBoundarySolution
 !  
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS(CGSEMSolution_2D), INTENT(in) :: myCGS
   REAL(prec), INTENT(out)             :: theSolution(0:myCGS % nMax,1:myCGS % nEq,1:4)

      theSolution = myCGS % boundarySolution

 END SUBROUTINE GetBoundarySolution_CGSEMSolution_2D
!
!
!
 SUBROUTINE SetBoundarySolution_CGSEMSolution_2D( myCGS, theSolution  )
 ! S/R SetBoundarySolution
 !  
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS(CGSEMSolution_2D), INTENT(inout) :: myCGS
   REAL(prec), INTENT(in)                 :: theSolution(0:myCGS % nMax,1:myCGS % nEq,1:4)
   ! LOCAL

      myCGS % boundarySolution = theSolution

 END SUBROUTINE SetBoundarySolution_CGSEMSolution_2D
!
!
!
 SUBROUTINE GetBoundarySolutionAtBoundary_CGSEMSolution_2D( myCGS, boundary, theSolution  )
 ! S/R GetBoundarySolution
 !  
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS(CGSEMSolution_2D), INTENT(in) :: myCGS
   INTEGER, INTENT(in)                 :: boundary
   REAL(prec), INTENT(out)             :: theSolution(0:myCGS % nMax,1:myCGS % nEq)

      theSolution = myCGS % boundarySolution(0:myCGS % nMax,1:myCGS % nEq,boundary)

 END SUBROUTINE GetBoundarySolutionAtBoundary_CGSEMSolution_2D
!
!
!
 SUBROUTINE SetBoundarySolutionAtBoundary_CGSEMSolution_2D( myCGS, boundary, theSolution  )
 ! S/R SetBoundarySolution
 !  
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS(CGSEMSolution_2D), INTENT(inout) :: myCGS
   INTEGER, INTENT(in)                    :: boundary
   REAL(prec), INTENT(in)                 :: theSolution(0:myCGS % nMax,1:myCGS % nEq)
   ! LOCAL

      myCGS % boundarySolution(0:myCGS % nMax,1:myCGS % nEq,boundary) = theSolution

 END SUBROUTINE SetBoundarySolutionAtBoundary_CGSEMSolution_2D
!
!
!
 SUBROUTINE GetBoundarySolutionAtBoundaryWithVarID_CGSEMSolution_2D( myCGS, boundary, varID, theSolution  )
 ! S/R GetBoundarySolution
 !  
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS(CGSEMSolution_2D), INTENT(in) :: myCGS
   INTEGER, INTENT(in)                 :: boundary, varID 
   REAL(prec), INTENT(out)             :: theSolution(0:myCGS % nMax)

      theSolution = myCGS % boundarySolution(0:myCGS % nMax,varID,boundary)

 END SUBROUTINE GetBoundarySolutionAtBoundaryWithVarID_CGSEMSolution_2D
!
!
!
 SUBROUTINE SetBoundarySolutionAtBoundaryWithVarID_CGSEMSolution_2D( myCGS, boundary, varID, theSolution  )
 ! S/R SetBoundarySolution
 !  
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS(CGSEMSolution_2D), INTENT(inout) :: myCGS
   INTEGER, INTENT(in)                    :: boundary, varID 
   REAL(prec), INTENT(in)                 :: theSolution(0:myCGS % nMax)
   ! LOCAL

      myCGS % boundarySolution(0:myCGS % nMax,varID, boundary) = theSolution

 END SUBROUTINE SetBoundarySolutionAtBoundaryWithVarID_CGSEMSolution_2D
!
! --------------------------------- Boundary Flux ------------------------------------------------ !
!
SUBROUTINE GetBoundaryFlux_CGSEMSolution_2D( myCGS, theFlux  )
 ! S/R GetBoundaryFlux
 !  
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS(CGSEMSolution_2D), INTENT(in) :: myCGS
   REAL(prec), INTENT(out)             :: theFlux(0:myCGS % nMax,1:myCGS % nEq,1:4)

      theFlux = myCGS % boundaryFlux

 END SUBROUTINE GetBoundaryFlux_CGSEMSolution_2D
!
!
!
 SUBROUTINE SetBoundaryFlux_CGSEMSolution_2D( myCGS, theFlux  )
 ! S/R SetBoundaryFlux
 !  
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS(CGSEMSolution_2D), INTENT(inout) :: myCGS
   REAL(prec), INTENT(in)                 :: theFlux(0:myCGS % nMax,1:myCGS % nEq,1:4)
   ! LOCAL

      myCGS % boundaryFlux = theFlux

 END SUBROUTINE SetBoundaryFlux_CGSEMSolution_2D
!
!
!
 SUBROUTINE GetBoundaryFluxAtBoundary_CGSEMSolution_2D( myCGS, boundary, theFlux  )
 ! S/R GetBoundaryFlux
 !  
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS(CGSEMSolution_2D), INTENT(in) :: myCGS
   INTEGER, INTENT(in)                 :: boundary
   REAL(prec), INTENT(out)             :: theFlux(0:myCGS % nMax,1:myCGS % nEq)

      theFlux = myCGS % boundaryFlux(0:myCGS % nMax,1:myCGS % nEq,boundary)

 END SUBROUTINE GetBoundaryFluxAtBoundary_CGSEMSolution_2D
!
!
!
 SUBROUTINE SetBoundaryFluxAtBoundary_CGSEMSolution_2D( myCGS, boundary, theFlux  )
 ! S/R SetBoundaryFlux
 !  
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS(CGSEMSolution_2D), INTENT(inout) :: myCGS
   INTEGER, INTENT(in)                    :: boundary
   REAL(prec), INTENT(in)                 :: theFlux(0:myCGS % nMax,1:myCGS % nEq)
   ! LOCAL

      myCGS % boundaryFlux(0:myCGS % nMax,1:myCGS % nEq,boundary) = theFlux

 END SUBROUTINE SetBoundaryFluxAtBoundary_CGSEMSolution_2D
!
!
!
  SUBROUTINE GetBoundaryFluxAtBoundaryNode_CGSEMSolution_2D( myCGS, boundary, i, theFlux  )
 ! S/R GetBoundaryFluxAtBoundaryNode
 !  
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS(CGSEMSolution_2D), INTENT(in) :: myCGS
   INTEGER, INTENT(in)                 :: boundary, i
   REAL(prec), INTENT(out)             :: theFlux(1:myCGS % nEq)

      theFlux = myCGS % boundaryFlux(i,1:myCGS % nEq,boundary)

 END SUBROUTINE GetBoundaryFluxAtBoundaryNode_CGSEMSolution_2D
!
!
!
 SUBROUTINE SetBoundaryFluxAtBoundaryNode_CGSEMSolution_2D( myCGS, boundary, i, theFlux  )
 ! S/R SetBoundaryFluxAtBoundaryNode
 !  
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS(CGSEMSolution_2D), INTENT(inout) :: myCGS
   INTEGER, INTENT(in)                    :: boundary, i
   REAL(prec), INTENT(in)                 :: theFlux(1:myCGS % nEq)
   ! LOCAL

      myCGS % boundaryFlux(i,1:myCGS % nEq,boundary) = theFlux

 END SUBROUTINE SetBoundaryFluxAtBoundaryNode_CGSEMSolution_2D
!
!
!
 SUBROUTINE GetBoundaryFluxAtBoundaryWithVarID_CGSEMSolution_2D( myCGS, boundary, varID, theFlux  )
 ! S/R GetBoundaryFlux
 !  
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS(CGSEMSolution_2D), INTENT(in) :: myCGS
   INTEGER, INTENT(in)                 :: boundary, varID 
   REAL(prec), INTENT(out)             :: theFlux(0:myCGS % nMax)

      theFlux = myCGS % boundaryFlux(0:myCGS % nMax,varID,boundary)

 END SUBROUTINE GetBoundaryFluxAtBoundaryWithVarID_CGSEMSolution_2D
!
!
!
 SUBROUTINE SetBoundaryFluxAtBoundaryWithVarID_CGSEMSolution_2D( myCGS, boundary, varID, theFlux  )
 ! S/R SetBoundaryFlux
 !  
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS(CGSEMSolution_2D), INTENT(inout) :: myCGS
   INTEGER, INTENT(in)                    :: boundary, varID 
   REAL(prec), INTENT(in)                 :: theFlux(0:myCGS % nMax)
   ! LOCAL

      myCGS % boundaryFlux(0:myCGS % nMax,varID, boundary) = theFlux

 END SUBROUTINE SetBoundaryFluxAtBoundaryWithVarID_CGSEMSolution_2D
!
!
!==================================================================================================!
!--------------------------------- Type Specific Routines -----------------------------------------!
!==================================================================================================!
!
!
 SUBROUTINE CalculateSolutionAtBoundaries_CGSEMSolution_2D( myCGS, dgStorage )
 ! S/R CalculateSolutionAtBoundaries
 !  
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS(CGSEMSolution_2D), INTENT(inout) :: myCGS
   TYPE(NodalStorage_2D), INTENT(in)      :: dgStorage
   ! LOCAL
   INTEGER :: iS, iP, iEq, nS, nP, nEq
   REAL(prec) :: lagSouth(0:myCGS % nP)
   REAL(prec) :: lagNorth(0:myCGS % nP)
   REAL(prec) :: lagWest(0:myCGS % nS)
   REAL(prec) :: lagEast(0:myCGS % nS)
   
       CALL myCGS % GetNumberOfEquations( nEq )
       CALL myCGS % GetNumberOfNodes( nS, nP )
       
       CALL dgStorage % GetSouthernInterpolants( lagSouth )
       CALL dgStorage % GetEasternInterpolants( lagEast )
       CALL dgStorage % GetNorthernInterpolants( lagNorth )
       CALL dgStorage % GetWesternInterpolants( lagWest )
       
       DO iP = 0, nP ! Loop over p-poinys
          DO iEq = 1, nEq ! Loop over the number of equations
          
             ! Setting the "west boundary"
              myCGS % boundarySolution(iP,iEq,4) = DOT_PRODUCT( lagWest, &
                                                        myCGS % solution(:,iP,iEq) )

            ! Setting the "east boundary"
              myCGS % boundarySolution(iP,iEq,2) = DOT_PRODUCT( lagEast, &
                                                        myCGS % solution(:,iP,iEq) )

          ENDDO ! iEq, loop over the number of equations
       ENDDO ! iP, loop over p-points
         

       DO iS = 0, nS ! Loop over s - points
          DO iEq = 1, nEq ! Loop over the number of equations

             ! Setting the "south boundary"
             myCGS % boundarySolution(iS,iEq,1)  = DOT_PRODUCT( lagSouth, &
                                                        myCGS % solution(iS,:,iEq) )

            ! Setting the "north boundary"
             myCGS % boundarySolution(iS,iEq,3)  = DOT_PRODUCT( lagNorth, &
                                                        myCGS % solution(iS,:,iEq) )

          ENDDO ! iEq, loop over the number of equations
       ENDDO ! iS, loop over s-points
   
 END SUBROUTINE CalculateSolutionAtBoundaries_CGSEMSolution_2D
!
!
!
END MODULE CGSEMSolutionStorageClass_2D
