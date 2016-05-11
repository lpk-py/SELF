!  CGsemElliptic_2D_Class.f90
!  
!  Copyright 2015 Joseph Schoonover <schoonover.numerics@gmail.com>
!  
!  CGsemElliptic_2D_Class.f90 is part of the Spectral Element Libraries in Fortran (SELF).
!
!  Permission is hereby granted, free of charge, to any person obtaining a copy of this software
!  and associated documentation files (the "Software"), to deal in the Software without restriction,
!  including without limitation the rights to use, copy, modify, merge, publish, distribute, 
!  sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is 
!  furnished to do so, subject to the following conditions:
!
!  The above copyright notice and this permission notice shall be included in all copies or 
!  substantial portions of the Software.
!
!  THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING
!  BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
!  NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM,
!  DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, 
!  OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
!  
! //////////////////////////////////////////////////////////////////////////////////////////////// !
!
! ========================================= Logs ================================================= !
! yyyy-mm-dd  Joe  <joe@clay>
!    This module defines the data structure used for implementing the Continuous Galerkin
!    Spectral Element Method. 
!
!    The module  is set up so that we solve div( Flux ) = Source
!   
! 
! //////////////////////////////////////////////////////////////////////////////////////////////// !


MODULE CGsemElliptic_2D_Class

! src/common/
USE ModelPrecision
USE ConstantsDictionary
USE ModelFlags
USE CommonRoutines
! src/nodal/
USE NodalStorage_2D_Class
! src/geom/
USE QuadMeshClass
USE EdgeClass
USE NodeClass_2D
USE VectorClass
USE MappedGeometryClass_2D
! src/cgsem/
USE CGsemParams_Class
USE CGsemPreconditioner_2D_Class


IMPLICIT NONE

      TYPE CGsemElliptic_2D
         INTEGER                        :: maxIters, nPlot, nS, nP
         REAL(prec)                     :: tolerance
         TYPE( CGsemPreconditioner_2D ) :: preconditioner
         TYPE( NodalStorage_2D )        :: SpectralOps
         TYPE( QuadMesh )               :: mesh           
         REAL(prec), ALLOCATABLE        :: solution(:,:,:)
         REAL(prec), ALLOCATABLE        :: source(:,:,:)
         REAL(prec), ALLOCATABLE        :: fluxCoeff(:,:,:)
         REAL(prec), ALLOCATABLE        :: resi(:)
         REAL(prec), ALLOCATABLE        :: plMatS(:,:), plMatP(:,:)

         CONTAINS

         PROCEDURE :: Build         => Build_CGsemElliptic_2D
         PROCEDURE :: BuildQuadMesh => BuildQuadMesh_CGsemElliptic_2D
         PROCEDURE :: Trash         => Trash_CGsemElliptic_2D
       
         

         PROCEDURE :: Mask                          => Mask_CGsemElliptic_2D
         PROCEDURE :: UnMask                        => UnMask_CGsemElliptic_2D
         PROCEDURE :: GlobalSum                     => GlobalSum_CGsemElliptic_2D
         PROCEDURE :: FluxDivergence                => FluxDivergence_CGsemElliptic_2D
         PROCEDURE :: SetThisDirichletEdge          => SetThisDirichletEdge_CGsemElliptic_2D
         PROCEDURE :: SetDirichletBoundaryCondition => SetDirichletBoundaryCondition_CGsemElliptic_2D
         PROCEDURE :: CalculateGradient             => CalculateGradient_CGsemElliptic_2D
         PROCEDURE :: MatrixAction                  => MatrixAction_CGsemElliptic_2D
         PROCEDURE :: Residual                      => Residual_CGsemElliptic_2D
         PROCEDURE :: DotProduct                    => DotProduct_CGsemElliptic_2D
         PROCEDURE :: Solve                         => Solve_PCConjugateGradient
         
         PROCEDURE :: CoarseToFine  => CoarseToFine_CGsemElliptic_2D
         PROCEDURE :: WriteTecplot  => WriteTecplot_CGsemElliptic_2D
         PROCEDURE :: WriteResidual => WriteResidual_CGsemElliptic_2D
         PROCEDURE :: ConstructMatrix
      END TYPE CGsemElliptic_2D

      

CONTAINS

 SUBROUTINE Build_CGsemElliptic_2D( myCGSEM, params )
 ! S/R Build
 !
 ! 
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS( CGsemElliptic_2D ), INTENT(inout) :: myCGSEM
   TYPE( CGsemParams ), INTENT(in)          :: params
   ! LOCAL
   INTEGER                 :: nS, nP, nPlot
   REAL(prec), ALLOCATABLE :: sNew(:)

      nS = params % polyDeg
      nP = nS
      
      nPlot = params % nPlot
      myCGSEM % tolerance = params % tolerance
      myCGSEM % maxIters  = params % maximumIterates
      myCGSEM % nS        = nS
      myCGSEM % nP        = nP 
      myCGSEM % nPlot     = nPlot
     
      CALL myCGSEM % SpectralOps % Build( nS, nP, GAUSS_LOBATTO, CG )
      CALL myCGSEM % BuildQuadMesh( params )
      
      ! Allocate space for the solution, source term, and diffusivity coefficient
      ALLOCATE( myCGSEM % solution(0:nS, 0:nP, 1:myCGSEM % mesh % nElems), &
                myCGSEM % source(0:nS, 0:nP, 1:myCGSEM % mesh % nElems), &
                myCGSEM % fluxCoeff(0:nS, 0:nP, 1:myCGSEM % mesh % nElems) )

      ALLOCATE( myCGSEM % resi(0:myCGSEM % maxIters) )

      myCGSEM % solution  = ZERO
      myCGSEM % source    = ZERO
      myCGSEM % fluxCoeff = ONE

      ALLOCATE( sNew(0:nPlot) )
      sNew = UniformPoints( -ONE, ONE, nPlot )
      ALLOCATE( myCGSEM % plMatS(0:nPlot,0:nS), myCGSEM % plMatP(0:nPlot,0:nP) )
      ! Build the plotting matrix
      CALL myCGSEM % SpectralOps % interp % CalculateInterpolationMatrix( nPlot, nPlot, sNew, sNew, &
                                                                          myCGSEM % plMatS, &
                                                                          myCGSEM % plMatP )
      DEALLOCATE( sNew )
      

     CALL myCGSEM % preconditioner % Build( params % pcTolerance, &
                                            params % pcIterates, &
                                            myCGSEM % mesh )

     CALL myCGSEM % preconditioner % MapFromOrigToPC( myCGSEM % fluxCoeff, &
                                                      myCGSEM % preconditioner % fluxCoeff, &
                                                      nS, nP, myCGSEM % mesh % nElems )
      
 END SUBROUTINE Build_CGsemElliptic_2D
!
!
!
 SUBROUTINE BuildQuadMesh_CGsemElliptic_2D( myCGSEM, params )
 ! S/R BuildQuadMesh
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS( CGsemElliptic_2D ), INTENT(inout) :: myCGSEM
   TYPE( CGsemParams ), INTENT(in)          :: params
   ! LOCAL
   INTEGER :: iEdge, eID, iNode, nID

      IF( TRIM( params % SpecMeshFile ) == nada )THEN
         PRINT*,' Loading default 2-D mesh.'
         CALL  myCGSEM % mesh % LoadDefaultMesh( myCGSEM % SpectralOps % interp, &
                                                 params % nXelem, &
                                                 params % nYelem )
      ELSE
      ! Builds the lateral mesh
         PRINT*, 'Reading 2-D mesh from '//trim(params % SpecMeshFile)//'.'
         CALL myCGSEM % mesh % ReadSpecMeshFile( myCGSEM % SpectralOps % interp, &
                                                 params % SpecMeshFile )
      ENDIF


      DO iEdge = 1, myCGSEM % mesh % nEdges
         CALL myCGSEM % mesh % GetEdgeSecondaryElementID( iEdge, eID )
         IF( eID == NO_NORMAL_FLOW )THEN
            CALL myCGSEM % mesh % SetEdgeSecondaryElementID( iEdge, DIRICHLET )
         ENDIF
      ENDDO

      DO iNode = 1, myCGSEM % mesh % nNodes
         CALL myCGSEM % mesh % GetNodeType( iNode, nID )
         IF( nID == NO_NORMAL_FLOW )THEN
            CALL myCGSEM % mesh % SetNodeType( iNode, DIRICHLET )
         ENDIF
      ENDDO

 END SUBROUTINE BuildQuadMesh_CGsemElliptic_2D
!
!
!
 SUBROUTINE Trash_CGsemElliptic_2D( myCGSEM )
 ! S/R Trash
 !
 ! 
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS( CGsemElliptic_2D ), INTENT(inout) :: myCGSEM


      ! Trash the nodal storage structure
      CALL myCGSEM % SpectralOps % Trash( )

      ! Trash the geometry
      CALL myCGSEM % mesh % Trash( )

      ! Deallocate space for the solution, source term, and diffusivity coefficient
      DEALLOCATE( myCGSEM % fluxCoeff, &
                  myCGSEM % solution,&
                  myCGSEM % source, &
                  myCGSEM % resi )
                  
      DEALLOCATE( myCGSEM % plMatS, myCGSEM % plMatP )

      CALL myCGSEM % preconditioner % Trash( )

 END SUBROUTINE Trash_CGsemElliptic_2D
!
!
!==================================================================================================!
!-------------------------------------- Type Specific ---------------------------------------------!
!==================================================================================================!
!
!
 SUBROUTINE SetDirichletBoundaryCondition_CGsemElliptic_2D( myCGSEM, dirF )
 ! S/R SetDirichletBoundaryCondition
 !
 !
 ! ============================================================================= !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS( CGsemElliptic_2D ), INTENT(inout) :: myCGSEM
   ! dirF is a function which is passed to this subroutine
   REAL(prec), EXTERNAL :: dirF !(x, y)
   !LOCAL
   INTEGER :: nS, nP, iEdge, nEl
   INTEGER :: e1, e2, s1, nodeType, i, j, eID, nID, iNode
   REAL(prec) :: x, y 

      nS = myCGSEM % nS
      nP = myCGSEM % nP
      nEl = myCGSEM % mesh % nElems
   
      ! Mask out the solution over the edges -- excluding the corner nodes
      DO iEdge = 1, myCGSEM % mesh % nEdges ! Loop over the edges
 
         CALL myCGSEM % mesh % GetEdgeSecondaryElementID( iEdge, e2 )
         IF( e2  == DIRICHLET )then ! This is a dirichlet boundary edge
      
            CALL myCGSEM % mesh % GetEdgePrimaryElementID( iEdge, e1 )
            CALL myCGSEM % mesh % GetEdgePrimaryElementSide( iEdge, s1 )

            CALL myCGSEM % SetThisDirichletEdge( e1, s1, dirF )

         ENDIF

      ENDDO ! iEdge, loop over the edges     


      ! Mask out the corner-nodes appropriately ( Leave the primary element intact )
      DO iNode = 1, myCGSEM % mesh % nNodes ! Loop over the corner nodes

         CALL myCGSEM % mesh % GetNodeType( iNode, nodeType )

         IF( nodeType == DIRICHLET ) then ! this is a prescribed (Dirichlet) boundary

            ! Rewind to the head of the linked list
            CALL myCGSEM % mesh % nodes(iNode) % RewindElementList( ) 

            ! Loop over the elements in this corner-node's connectivity list
            DO WHILE( myCGSEM % mesh % nodes(iNode) % CurrentElemExists( ) ) ! while there are elements in the node-to-element list

               CALL myCGSEM % mesh % nodes(iNode) % GetCurrentElementID( eID ) ! Get the element ID and the local node ID (1->4)
               CALL myCGSEM % mesh % GetLocalNodeID( iNode, nID )

               i = myCGSEM % mesh % cornerMap(1, nID) ! Grab the index over the first computational direction associated with this corner node
               j = myCGSEM % mesh % cornerMap(2, nID) ! Grab the index over the second computational direction associated with this corner node

               CALL myCGSEM % mesh % GetPositionAtNode( eID, x, y, i, j )

               myCGSEM % solution(i,j,eID) = dirF( x, y )

               CALL myCGSEM % mesh % nodes(iNode) % AdvanceElementList( )

            ENDDO ! while we have elements in the node-to-element list

         ENDIF

      ENDDO


 END SUBROUTINE SetDirichletBoundaryCondition_CGsemElliptic_2D
!
!
!
 SUBROUTINE SetThisDirichletEdge_CGsemElliptic_2D( myCGSEM, eID, sID, dirF )
 ! S/R SetThisDirichletEdge
 !
 !  This subroutine sets the solutions on side "sID" of element "eID",
 !  excluding the corner nodes ( Assumes Gauss-Lobatto points are used )
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS( CGsemElliptic_2D ), INTENT(inout) :: myCGSEM
   INTEGER, INTENT(in)                      :: eID, sID
   REAL(prec), EXTERNAL :: dirF !(x, y)
   ! LOCAL
   INTEGER    :: iS, iP, nS, nP,nElems
   REAL(prec) :: x, y
   

      nS     = myCGSEM % nS
      nP     = nS
      nElems = myCGSEM % mesh % nElems

      IF( sID == east .OR. sID == west )then ! east or west sides

         iS = myCGSEM % mesh % sideMap(sID) 

         DO iP = 1, nP-1 ! Loop over the edge, excluding the corner nodes
            CALL myCGSEM % mesh % GetPositionAtNode( eID, x, y, iS, iP )
            myCGSEM % solution(iS,iP,eID) = dirF( x, y )
         ENDDO ! iP, loop over the edge

      ELSE ! south or north sides

         iP = myCGSEM % mesh % sideMap(sID)

         DO iS = 1, nS-1 ! Loop over the edge, excluding the corner nodes
            CALL myCGSEM % mesh % GetPositionAtNode( eID, x, y, iS, iP )
            myCGSEM % solution(iS,iP,eID) = dirF( x, y )
         ENDDO

      ENDIF

 END SUBROUTINE SetThisDirichletEdge_CGsemElliptic_2D
!
!
!
 SUBROUTINE Mask_CGsemElliptic_2D( myCGSEM, u )
 ! S/R Mask_CGsemElliptic_2D
 !
 !
 ! ============================================================================= !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS( CGsemElliptic_2D ), INTENT(inout) :: myCGSEM
   REAL(prec), INTENT(inout)                :: u(0:myCGSEM % nS, &
                                                 0:myCGSEM % nP, &
                                                 1:myCGSEM % mesh % nElems )
   !LOCAL
   INTEGER :: nS, nP, iEdge, eID, iNode, nID, i, j, nEl, nEdges, nNodes
   INTEGER :: e1, e2, s1, s2, nodeType

      nS     = myCGSEM % nS
      nP     = myCGSEM % nP
      nEl    = myCGSEM % mesh % nElems
      nEdges = myCGSEM % mesh % nEdges
      nNodes = myCGSEM % mesh % nNodes
   
      ! Mask out the solution over the edges -- excluding the corner nodes
      DO iEdge = 1, nEdges ! Loop over the edges
 
         CALL myCGSEM % mesh % GetEdgeSecondaryElementID( iEdge, e2 )

         IF( e2 == DIRICHLET )THEN ! check the secondary element for boundary condition flag

            ! This is a prescribed boundary condition (Dirichlet) so we mask out the solution along this edge   

            CALL myCGSEM % mesh % GetEdgePrimaryElementID( iEdge, e1 )
            CALL myCGSEM % mesh % GetEdgePrimaryElementSide( iEdge, s1 )

            CALL MaskSide( e1, s1, myCGSEM % mesh, u, nS, nP, nEl )

         ELSE ! then this is not a prescribed boundary, and we choose to mask out the secondary element

            IF( eID > 0 )THEN ! this is an internal edge, and we should mask out the secondary element
               CALL myCGSEM % mesh % GetEdgeSecondaryElementID( iEdge, e2 )
               CALL myCGSEM % mesh % GetEdgeSecondaryElementSide( iEdge, s2 )
               CALL MaskSide( e2, s2, myCGSEM % mesh, u, nS, nP, nEl )
            ENDIF

         ENDIF

      ENDDO ! iEdge, loop over the edges     

      ! At this point, the secondary internal edges and prescribed boundary condition edges have been masked out.
      ! Now we mask out all but the primary element in the corner-node lists.

      ! Mask out the corner-nodes appropriately ( Leave the primary element intact )
      DO iNode = 1, nNodes ! Loop over the corner nodes

         CALL myCGSEM % mesh % nodes(iNode) % RewindElementList( )
         CALL myCGSEM % mesh % GetNodeType( iNode, nodeType )

         IF( nodeType /= DIRICHLET ) then ! this is NOT a prescribed (Dirichlet) boundary
                                                                          ! The data in this non-prescribed node is preserved

            CALL myCGSEM % mesh % nodes(iNode) % AdvanceElementList( ) ! skip the first element in the list so it's data is not masked

         ENDIF

         DO WHILE( myCGSEM % mesh % nodes(iNode) % CurrentElemExists( ) ) ! while there are elements in the node-to-element list

            CALL myCGSEM % mesh % nodes(iNode) % GetCurrentElementID( eID ) ! Get the element ID and the local node ID (1->4)
            CALL myCGSEM % mesh % GetLocalNodeID( iNode, nID )

            i = myCGSEM % mesh % cornerMap(1, nID) ! Grab the index over the first computational direction associated with this corner node
            j = myCGSEM % mesh % cornerMap(2, nID) ! Grab the index over the second computational direction associated with this corner node
 
            u(i,j,eID) = ZERO ! Mask out the solution at this element at this corner node

            CALL myCGSEM % mesh % nodes(iNode) % AdvanceElementList( )

         ENDDO ! while we have elements in the node-to-element list

      ENDDO


 END SUBROUTINE Mask_CGsemElliptic_2D
!
!
!
 SUBROUTINE MaskSide( eID, sID, mesh, u, nS, nP, nElems )
 ! S/R MaskSide
 !
 !  This subroutine masks the solutions on side "sID" of element "eID",
 !  excluding the corner nodes ( Assumes Gauss-Lobatto points are used )
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   INTEGER, INTENT(in)        :: eID, sID, nS, nP,nElems
   REAL(prec), INTENT(inout)  :: u(0:nS,0:nP,1:nElems)
   TYPE(QuadMesh), INTENT(in) :: mesh
   ! LOCAL
   INTEGER :: iS, iP, si
   
      si = ABS(sID)
      IF( si == East .OR. si == West )then ! east or west sides

         iS = mesh % sideMap(si) 
         DO iP = 1, nP-1 ! Loop over the edge, excluding the corner nodes
            u(iS,iP,eID) = ZERO
         ENDDO ! iP, loop over the edge

      ELSE ! south or north sides

         iP = mesh % sideMap(si)
         DO iS = 1, nS-1 ! Loop over the edge, excluding the corner nodes
            u(iS,iP,eID) = ZERO
         ENDDO

      ENDIF

 END SUBROUTINE MaskSide
!
!
!
 SUBROUTINE UnMask_CGsemElliptic_2D( myCGSEM, u )
 ! S/R UnMask
 !
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS( CGsemElliptic_2D ), INTENT(inout) :: myCGSEM
   REAL(prec), INTENT(inout)                :: u(0:myCGSEM % nS, &
                                                 0:myCGSEM % nP, &
                                                 1:myCGSEM % mesh % nElems )
   !LOCAL
   INTEGER :: nS, nP, iEdge, eID, iNode, nID, i, j, eID0, i0, j0, nEl, nEdges, nNodes
   INTEGER :: nodeType, e2

      nS     = myCGSEM % nS
      nP     = myCGSEM % nP
      nEl    = myCGSEM % mesh % nElems
      nEdges = myCGSEM % mesh % nEdges
      nNodes = myCGSEM % mesh % nNodes
   
      ! Unmask the solution over the interior edges only -- corner nodes are excluded in this step
      ! and are taken care of in the following step
      DO iEdge = 1, nEdges ! Loop over the edges
 
         CALL myCGSEM % mesh % GetEdgeSecondaryElementID( iEdge, e2 )

         IF( e2 > 0 )THEN ! This is an interior shared edge
            CALL UnMaskSide( iEdge, myCGSEM % mesh, u, nS, nP, nEl ) ! unmask the solution along this edge
         ENDIF

      ENDDO ! iEdge, loop over the edges     


      ! Unmask the corner-nodes
      DO iNode = 1, nNodes ! Loop over the corner nodes
         
         CALL myCGSEM % mesh % GetNodeType( iNode, nodeType )

         IF( nodeType /= DIRICHLET ) then ! this is NOT a prescribed (Dirichlet) boundary

            CALL myCGSEM % mesh % nodes(iNode) % RewindElementList( )

            ! Gather the element and local node IDs of the corner which contains the unmasked solution values.
            CALL myCGSEM % mesh % nodes(iNode) % GetCurrentElementID( eID0 ) 
            CALL myCGSEM % mesh % GetLocalNodeID( iNode, nID ) ! Get the element ID and the local node ID (1->4)

            i0 = myCGSEM % mesh % cornerMap(1, nID) ! Grab the index over the first computational direction associated with this corner node
            j0 = myCGSEM % mesh % cornerMap(2, nID) ! Grab the index over the second computational direction associated with this corner node

            CALL myCGSEM % mesh % nodes(iNode) % AdvanceElementList( ) ! move on to the next element which shares this node


            ! Loop over the elements in this corner-node's connectivity list
            DO WHILE( myCGSEM % mesh % nodes(iNode) % CurrentElemExists( ) ) ! while there are elements in the node-to-element list

               CALL myCGSEM % mesh % nodes(iNode) % GetCurrentElementID( eID ) ! Get the element ID and the local node ID (1->4)
               CALL myCGSEM % mesh % GetLocalNodeID( iNode, nID )
 
               i = myCGSEM % mesh % cornerMap(1, nID) ! Grab the index over the first computational direction associated with this corner node
               j = myCGSEM % mesh % cornerMap(2, nID) ! Grab the index over the second computational direction associated with this corner node
  
               u(i,j,eID) = u(i0,j0,eID0) ! copy the solution from the unmasked solution

               CALL myCGSEM % mesh % nodes(iNode) % AdvanceElementList( )

            ENDDO ! while we have elements in the node-to-element list

         ENDIF

      ENDDO

 END SUBROUTINE UnMask_CGsemElliptic_2D
!
!
!
 SUBROUTINE UnMaskSide( edgeID, mesh, u, nS, nP, nElems )
 ! S/R UnMaskSide
 !
 !  This subroutine unmasks the solutions on side "sID" of element "eID",
 !  excluding the corner nodes ( Assumes Gauss-Lobatto points are used )
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   INTEGER, INTENT(in)        :: edgeID
   INTEGER, INTENT(in)        :: nS, nP,nElems
   REAL(prec), INTENT(inout)  :: u(0:nS,0:nP,1:nElems)
   TYPE(QuadMesh), INTENT(in) :: mesh
   ! LOCAL
   INTEGER :: iS, iP, eID, sID, jS, jP, inc
   REAL(prec) :: temp(0:nS) ! assume nS == nP

      ! Get the primary element and side IDs
      CALL mesh % GetEdgePrimaryElementID( edgeID, eID )
      CALL mesh % GetEdgePrimaryElementSide( edgeID, sID )
 
      IF( sID == East .OR. sID == West )THEN ! east or west sides

         iS = mesh % sideMap(sID) 
         DO iP = 1, nP-1 ! Loop over the edge, excluding the corner nodes
            temp(iP) = u(iS,iP,eID) ! Copy the solution in the primary element
         ENDDO ! iP, loop over the edge

      ELSE ! south or north sides

         iP = mesh % sideMap(sID)
         DO iS = 1, nS-1 ! Loop over the edge, excluding the corner nodes
            temp(iS) = u(iS,iP,eID) ! Copy the solution in the primary element
         ENDDO

      ENDIF

       ! Get the secondary element and side IDs
      CALL mesh % GetEdgeSecondaryElementID( edgeID, eID )
      CALL mesh % GetEdgeSecondaryElementSide( edgeID, sID )
      sID = abs(sID)
      IF( sID == East .OR. sID == West )THEN ! east or west sides

         iS = mesh % sideMap(sID) 
         CALL mesh % GetEdgeStart( edgeID, jP )
         CALL mesh % GetEdgeIncrement( edgeID, inc ) 

         DO iP = 1, nP-1 ! Loop over the edge, excluding the corner nodes
            u(iS,jP,eID) = temp(iP) ! Copy the solution to the secondary element
            jP = jP + inc ! increment jP
         ENDDO ! iP, loop over the edge

      ELSE ! south or north sides

         iP = mesh % sideMap(sID)
         CALL mesh % GetEdgeStart( edgeID, jS )
         CALL mesh % GetEdgeIncrement( edgeID, inc )

         DO iS = 1, nS-1 ! Loop over the edge, excluding the corner nodes

            u(jS,iP,eID) = temp(iS) ! Copy the solution in the primary element
            jS = jS + inc ! increment jS

         ENDDO

      ENDIF

 END SUBROUTINE UnMaskSide
!
!
!
 SUBROUTINE GlobalSum_CGsemElliptic_2D( myCGSEM, u )
 ! S/R GlobalSum
 !
 !  Adds together the shared edge and node contributions for the CGSEM method.
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS( CGsemElliptic_2D ), INTENT(inout) :: myCGSEM
   REAL(prec), INTENT(inout)                :: u(0:myCGSEM % nS, &
                                                 0:myCGSEM % nP, &
                                                 1:myCGSEM % mesh % nElems )
   !LOCAL
   INTEGER    :: iEdge, eID, iNode, nID, i, j, nEdges, nNodes, nS, nP, nEl
   INTEGER    :: e2, nodeType 
   REAL(prec) :: theSum

      nS     = myCGSEM % nS
      nP     = myCGSEM % nP
      nEl    = myCGSEM % mesh % nElems
      nEdges = myCGSEM % mesh % nEdges
      nNodes = myCGSEM % mesh % nNodes

      DO iEdge = 1, nEdges ! loop over the edges

         CALL myCGSEM % mesh % GetEdgeSecondaryElementID( iEdge, e2 )
         IF( e2 > 0 )then !/= DIRICHLET )then ! This is an interior shared edge (or a neumann boundary edge)

            CALL SumSide( iEdge, myCGSEM % mesh, u, nS, nP, nEl )  ! add the contributions from the elements which share this edge
                                                                   ! The edge sum excludes corner points
         ENDIF

      ENDDO ! iEdge, loop over the edges


      DO iNode = 1, nNodes ! loop over the corner nodes
      
         CALL myCGSEM % mesh % GetNodeType( iNode, nodeType )

         IF( nodeType /= DIRICHLET ) then ! this node does NOT lie on a prescribed (Dirichlet) boundary

             ! Rewind to the head of the linked list
            CALL myCGSEM % mesh % nodes(iNode) % RewindElementList( )

            ! *** First, compute the sum from the contributing elements *** !
            theSum = ZERO 

            ! Loop over the elements in this corner-node's connectivity list
            DO WHILE( myCGSEM % mesh % nodes(iNode) % CurrentElemExists( ) ) ! while there are elements in the node-to-element list

               CALL myCGSEM % mesh % nodes(iNode) % GetCurrentElementID( eID )
               CALL myCGSEM % mesh % GetLocalNodeID( iNode, nID )

               i = myCGSEM % mesh % cornerMap(1, nID) ! Grab the index over the first computational direction associated with this corner node
               j = myCGSEM % mesh % cornerMap(2, nID) ! Grab the index over the second computational direction associated with this corner node
  
               theSum = theSum + u(i,j,eID) 

               CALL myCGSEM % mesh % nodes(iNode) % AdvanceElementList( )

            ENDDO ! while we have elements in the node-to-element list

            ! *** Now, copy the sum to the array "u" for each element which contributed to the sum *** !

             ! Rewind to the head of the linked list
            CALL myCGSEM % mesh % nodes(iNode) % RewindElementList( )

            DO WHILE( myCGSEM % mesh % nodes(iNode) % CurrentElemExists( ) ) ! while there are elements in the node-to-element list

               CALL myCGSEM % mesh % nodes(iNode) % GetCurrentElementID( eID )
               CALL myCGSEM % mesh % GetLocalNodeID( iNode, nID ) ! Get the element ID and the local node ID (1->4)

               i = myCGSEM % mesh % cornerMap(1, nID) ! Grab the index over the first computational direction associated with this corner node
               j = myCGSEM % mesh % cornerMap(2, nID) ! Grab the index over the second computational direction associated with this corner node
  
               u(i,j,eID) = theSum

               CALL myCGSEM % mesh % nodes(iNode) % AdvanceElementList( )

            ENDDO ! while we have elements in the node-to-element list

         ENDIF

      ENDDO

 END SUBROUTINE GlobalSum_CGsemElliptic_2D
!
!
!
 SUBROUTINE SumSide( edgeID, mesh, u, nS, nP, nElems )
 ! S/R SumSide
 !
 ! 
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   INTEGER, INTENT(in)        :: edgeID
   INTEGER, INTENT(in)        :: nS, nP,nElems
   REAL(prec), INTENT(inout)  :: u(0:nS,0:nP,1:nElems)
   TYPE(QUADMESH), INTENT(in) :: mesh
   ! LOCAL
   INTEGER :: iS, iP, k, n, inc, e(1:2), s(1:2)
   REAL(prec) :: temp(0:nS,1:2) ! assume nS == nP
   REAL(prec) :: theSum



      ! Get the primary/secondary element and side IDs
      CALL mesh % GetEdgePrimaryElementID( edgeID, e(1) )
      CALL mesh % GetEdgePrimaryElementSide( edgeID, s(1) )
      CALL mesh % GetEdgeSecondaryElementID( edgeID, e(2) )
      CALL mesh % GetEdgeSecondaryElementSide( edgeID, s(2) )   
      s(2) = abs(s(2))

      DO k = 1, 2
         IF( s(k) == east .OR. s(k) == west )then ! east or west sides

            iS = mesh % sideMap(s(k)) 
            DO iP = 1, nP-1 ! Loop over the edge, excluding the corner nodes
               temp(iP,k) = u(iS,iP,e(k)) ! 
            ENDDO ! iP, loop over the edge

         ELSE ! south or north sides

            iP = mesh % sideMap(s(k))
            DO iS = 1, nS-1 ! Loop over the edge, excluding the corner nodes
               temp(iS,k) = u(iS,iP,e(k)) ! 
            ENDDO

         ENDIF

      ENDDO

      ! Add the two contributions together
      CALL mesh % GetEdgeStart( edgeID, n )
      CALL mesh % GetEdgeIncrement( edgeID, inc )
     
      DO iS = 1, nS-1

         theSum = temp(iS,1) + temp(n,2)
         temp(iS,1) = theSum
         temp(n,2) = theSum

         n = n + inc ! increment the secondary element array adress

      ENDDO

      ! copy the sum into the two contributing element edges

      DO k = 1, 2
   
         IF( s(k) == east .OR. s(k) == west )THEN ! east or west sides

            iS = mesh % sideMap(s(k)) 
            DO iP = 1, nP-1 ! Loop over the edge, excluding the corner nodes
                u(iS,iP,e(k)) = temp(iP,k) !
            ENDDO ! iP, loop over the edge

         ELSE ! south or north sides

            iP = mesh % sideMap(s(k))
            DO iS = 1, nS-1 ! Loop over the edge, excluding the corner nodes
               u(iS,iP,e(k)) = temp(iS,k) ! 
            ENDDO

         ENDIF

      ENDDO

 END SUBROUTINE SumSide
!
!
!
 SUBROUTINE CalculateGradient_CGsemElliptic_2D( myCGSEM, iEl, u, dudx, dudy )
 ! S/R CalculateGradient
 !
 !  Calculates the gradient of the solution in computational coordinates within 
 !  a single element
 !
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS( CGsemElliptic_2D ), INTENT(in) :: myCGSEM 
   INTEGER, INTENT(in)                   :: iEl
   REAL(prec), INTENT(in)                :: u(0:myCGSEM % nS, 0:myCGSEM % nP)
   REAL(prec), INTENT(out)               :: dudx(0:myCGSEM % nS, 0:myCGSEM % nP)
   REAL(prec), INTENT(out)               :: dudy(0:myCGSEM % nS, 0:myCGSEM % nP) 
   ! LOCAL
   INTEGER :: iS, iP, nS, nP
   REAL(prec) :: duds(0:myCGSEM % nS, 0:myCGSEM % nP)
   REAL(prec) :: dudp(0:myCGSEM % nS, 0:myCGSEM % nP)
   REAL(prec) :: dMatS(0:myCGSEM % nS, 0:myCGSEM % nS)
   REAL(prec) :: dMatP(0:myCGSEM % nP, 0:myCGSEM % nP)
   REAL(prec) :: sTemp(0:myCGSEM % nS)
   REAL(prec) :: pTemp(0:myCGSEM % nP)
   REAL(prec) :: J, dxds, dxdp, dyds, dydp

      nS = myCGSEM % nS
      nP = myCGSEM % nP
      CALL myCGSEM % SpectralOps % GetDerivativeMatrix( dMatS, dMatP ) 

      DO iP = 0, nP ! Loop over the second computational direction
 
         sTemp = u(0:nS,iP)
         ! Calculate the derivative in the first computational direction
         duds(0:nS,iP) = MATMUL( dMatS, sTemp )

      ENDDO ! iP, loop over the second computational direction

      DO iS = 0, nS ! Loop over the first computational direction

         pTemp = u(iS,0:nP)
         ! Calculate the derivative in the second computational direction
         dudp(iS,0:nP) = MATMUL( dMatP, pTemp )

      ENDDO ! iS, loop over the first computational direction


      ! Calculate the gradient in physical space
      DO iP = 0, nP
         DO iS = 0, nS

            CALL myCGSEM % mesh % GetJacobianAtNode( iEl, J, iS, iP )
            CALL myCGSEM % mesh % GetCovariantMetricsAtNode( iEl, dxds, dxdp, dyds, dydp, iS, iP )

            dudx(iS,iP) = ( dydp*duds(iS,iP) - dyds*dudp(iS,iP) )/J
            dudy(iS,iP) = ( dxds*dudp(iS,iP) - dxdp*duds(iS,iP) )/J        

         ENDDO
      ENDDO

 END SUBROUTINE CalculateGradient_CGsemElliptic_2D
!
!
!
 SUBROUTINE FluxDivergence_CGsemElliptic_2D( myCGSEM, iEl, u, Lu )
 ! S/R FluxDivergence
 !
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS( CGsemElliptic_2D ), INTENT(inout) :: myCGSEM
   INTEGER, INTENT(in)                      :: iEl
   REAL(prec), INTENT(in)                   :: u(0:myCGSEM % nS, 0:myCGSEM % nP)
   REAL(prec), INTENT(out)                  :: Lu(0:myCGSEM % nS, 0:myCGSEM % nP)
   ! LOCAL
   INTEGER    :: iS, iP, nS, nP
   REAL(prec) :: F1(0:myCGSEM % nS, 0:myCGSEM % nP)
   REAL(prec) :: F2(0:myCGSEM % nS, 0:myCGSEM % nP)
   REAL(prec) :: dF1ds(0:myCGSEM % nS, 0:myCGSEM % nP)
   REAL(prec) :: dF2dp(0:myCGSEM % nS, 0:myCGSEM % nP)
   REAL(prec) :: dudx(0:myCGSEM % nS, 0:myCGSEM % nP)
   REAL(prec) :: dudy(0:myCGSEM % nS, 0:myCGSEM % nP)
   REAL(prec) :: dMatS(0:myCGSEM % nS, 0:myCGSEM % nS)
   REAL(prec) :: dMatSt(0:myCGSEM % nS, 0:myCGSEM % nS)
   REAL(prec) :: dMatP(0:myCGSEM % nP, 0:myCGSEM % nP)
   REAL(prec) :: dMatPt(0:myCGSEM % nP, 0:myCGSEM % nP)
   REAL(prec) :: ws(0:myCGSEM % nS)
   REAL(prec) :: wp(0:myCGSEM % nP)
   REAL(prec) :: sTemp(0:myCGSEM % nS)
   REAL(prec) :: pTemp(0:myCGSEM % nP)
   REAL(prec) :: J, dxds, dxdp, dyds, dydp

      nS = myCGSEM % nS
      nP = myCGSEM % nP

      CALL myCGSEM % SpectralOps % GetQuadratureWeights( ws, wp )
      CALL myCGSEM % SpectralOps % GetDerivativeMatrix( dMatS, dMatP )
      dMatSt = TRANSPOSE( dMatS )
      dMatPt = TRANSPOSE( dMatP )

      CALL myCGSEM % CalculateGradient( iEl, u, dudx, dudy ) ! calculate the gradient in physical coordinates

      ! Now apply the metric terms
      DO iP = 0, nP ! Loop over the second computational direction
         DO iS = 0, nS ! Loop over the first computational direction

            CALL myCGSEM % mesh % GetCovariantMetricsAtNode( iEl, dxds, dxdp, dyds, dydp, iS, iP )

            ! Contravariant flux calculation
            F1(iS,iP) = myCGSEM % fluxCoeff(iS,iP,iEl)*( dudx(iS,iP)*dydp - dudy(iS,iP)*dxdp )*ws(iS)
            F2(iS,iP) = myCGSEM % fluxCoeff(iS,iP,iEl)*( dudy(iS,iP)*dxds - dudx(iS,iP)*dyds )*wp(iP)        

         ENDDO ! iS
      ENDDO ! iP

      ! Now calculate the divergence of the flux

      DO iP = 0, nP ! loop over the second computational direction
         sTemp = F1(0:nS,iP)
         dF1ds(0:nS,iP) = MATMUL( dMatSt, sTemp )  
      ENDDO ! iP

      DO iS = 0, nS ! loop over the second computational direction
         pTemp = F2(iS,0:nP)
         dF2dp(iS,0:nP) = MATMUL( dMatPt, pTemp )  
      ENDDO


      DO iP = 0, nP
         DO iS = 0, nS
            Lu(iS,iP) = -( dF1ds(iS,iP)*wp(iP) + dF2dp(iS,iP)*ws(iS) )
         ENDDO 
      ENDDO 
       

 END SUBROUTINE FluxDivergence_CGsemElliptic_2D 
!
!
!
 SUBROUTINE MatrixAction_CGsemElliptic_2D( myCGSEM, u, Au )
 ! S/R MatrixAction
 ! 
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS(CGsemElliptic_2D), INTENT(inout) :: myCGSEM
   REAL(prec), INTENT(inout)              :: u(0:myCGSEM % nS, &
                                               0:myCGSEM % nP, &
                                               1:myCGSEM % mesh % nElems )
   REAL(prec), INTENT(out)                :: Au(0:myCGSEM % nS, &
                                                0:myCGSEM % nP, &
                                                1:myCGSEM % mesh % nElems  )
   ! Local
   INTEGER :: iEl, nS, nP, nEl
   REAL(prec) :: temp(0:myCGSEM % nS, &
                      0:myCGSEM % nP )
   REAL(prec) :: Atemp(0:myCGSEM % nS, &
                       0:myCGSEM % nP )

      nS  = myCGSEM % nS
      nP  = myCGSEM % nP
      nEl = myCGSEM % mesh % nElems
     
      CALL myCGSEM % UnMask( u ) ! unmask the solution

!$OMP PARALLEL PRIVATE( temp, Atemp )
!$OMP DO 
      DO iEl = 1, nEl
         temp(0:nS,0:nP) = u(0:nS,0:nP,iEl)
         CALL myCGSEM % FluxDivergence( iEl, temp, Atemp )
         Au(0:nS,0:nP,iEl) = Atemp(0:nS,0:nP)
      ENDDO
!$OMP END DO
!$OMP  FLUSH ( Au )
!$OMP END PARALLEL

     ! Add the contributions from the shared edges and corners
     CALL myCGSEM % GlobalSum( Au )

     ! Mask the edges in preparation for residual calculation
     CALL myCGSEM % Mask( Au )

     ! Re-mask the solution
     CALL myCGSEM % Mask( u )

 END SUBROUTINE MatrixAction_CGsemElliptic_2D
!
!
!
 SUBROUTINE Residual_CGsemElliptic_2D( myCGSEM, dirF, r )
 ! S/R Residual
 ! 
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS( CGsemElliptic_2D ), INTENT(inout) :: myCGSEM
   REAL(prec), EXTERNAL                     :: dirF
   REAL(prec), INTENT(out)                  :: r(0:myCGSEM % nS, &
                                                 0:myCGSEM % nP, &
                                                 1:myCGSEM % mesh % nElems )
   ! LOCAL
   INTEGER    :: iS, iP, nS, nP, nEl, iEl, nEdges
   REAL(prec) :: J
   REAL(prec) :: ws(0:myCGSEM % nS)
   REAL(prec) :: wp(0:myCGSEM % nP)
   REAL(prec) :: Au(0:myCGSEM % nS, &
                    0:myCGSEM % nP, &
                    1:myCGSEM % mesh % nElems )



      nS = myCGSEM % nS 
      nP = myCGSEM % nP
      nEl = myCGSEM % mesh % nElems 
      nEdges = myCGSEM % mesh % nEdges

      CALL myCGSEM % SpectralOps % GetQuadratureWeights( ws, wp )
      CALL myCGSEM % MatrixAction( myCGSEM % solution, Au )

!$OMP PARALLEL PRIVATE( J, locR )
!$OMP DO      
      DO iEl = 1, nEl
         DO iP = 0, nP
            DO iS = 0, nS
               CALL myCGSEM % mesh % GetJacobianAtNode( iEl, J, iS, iP ) 
               r(iS,iP,iEl) = myCGSEM % source(iS,iP,iEl)*J*ws(iS)*wp(iP) - Au(iS,iP,iEl)! gather source terms
            ENDDO
         ENDDO
      ENDDO
!$OMP END DO
!$OMP  FLUSH ( r )
!$OMP END PARALLEL

  ! Add in the boundary flux
!      DO iEdge = 1, nEdges
!         CALL myCGSEM % mesh % GetEdgeSecondaryElementID( iEdge, e2 )
!         IF( e2 == NEUMANN_WALL )then 

!            ! Gather the local addressing information
!            CALL myCGSEM % mesh % GetEdgePrimaryElementID( iEdge, e1 )
!            CALL myCGSEM % mesh % GetEdgePrimaryElementSide( iEdge, s1 )
!            nB = myCGSEM % mesh % sideMap(s1) ! get the local quadrature node ID for this boundary

!            IF( s1 == South .OR. s1 == North )THEN

!               DO iS = 0, nS
!                  CALL myCGSEM % mesh % GetBoundaryNormalAtNode( e1, nhat, length, iS, s1 ) 
!                  rhs = myCGSEM % boundaryFlux(e1,s1,iS)*length!*bfac(s1) 
!                  r(iS, nB, e1) = r(iS, nB, e1) - rhs*ws(iS)
!               ENDDO

!            ELSEIF( s1 == East .OR. s1 == West )then ! east or west boundary (respectively)

!               DO iP = 0, nP
!                  CALL myCGSEM % mesh % GetBoundaryNormalAtNode( e1, nhat, length, iS, s1 )
!                  rhs = myCGSEM % boundaryFlux(e1,s1,iP)*length!*bfac(s1) 
!                  r(nB, iP, e1) = r(nB, iP, e1) - rhs*wp(iP)
!               ENDDO

!            ENDIF
 
!         ENDIF

!      ENDDO ! iEdge, loop over the edges

     ! Add the contributions from the shared edges and corners
     CALL myCGSEM % GlobalSum( r )

     ! Mask the edges in preparation for residual calculation
     CALL myCGSEM % Mask( r )

     ! Re-mask the solution
     CALL myCGSEM % Mask( myCGSEM % solution )

     ! re-set the dirichlet boundary conditions
     CALL myCGSEM % SetDirichletBoundaryCondition( dirF )

 END SUBROUTINE Residual_CGsemElliptic_2D
!
!
!
 FUNCTION DotProduct_CGsemElliptic_2D( myCGSEM, u, v ) RESULT( uDotv )
 ! FUNCTION DotProduct
 ! Computes the dot product of two "vectors" whose data are
 ! stored in the "(iS,iP)" format
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS( CGsemElliptic_2D ) :: myCGSEM
   REAL(prec)                :: u(0:myCGSEM % nS, &
                                  0:myCGSEM % nP, &
                                  1:myCGSEM % mesh % nElems )
   REAL(prec)                :: v(0:myCGSEM % nS, &
                                  0:myCGSEM % nP, &
                                  1:myCGSEM % mesh % nElems )
   REAL(prec)                :: uDotv
   ! Local
   INTEGER :: iP, nS, nP, iEl, nEl
   REAL(prec) :: uloc(0:myCGSEM % nS)
   REAL(prec) :: vloc(0:myCGSEM % nS)
   
      nS  = myCGSEM % nS 
      nP  = myCGSEM % nP
      nEl = myCGSEM % mesh % nElems

      ! mask the arrays
      CALL myCGSEM % Mask( u )
      CALL myCGSEM % Mask( v )
 
      uDotv = ZERO
      DO iEl = 1, nEl
         DO iP = 0, nP
            uloc  = u(0:nS,iP,iEl)
            vLoc  = v(0:nS,iP,iEl)
            uDotV = uDotV + DOT_PRODUCT( uloc, vloc )
         ENDDO
      ENDDO

 END FUNCTION DotProduct_CGsemElliptic_2D
!
!
!
 SUBROUTINE Solve_PCConjugateGradient( myCGSEM, dirF, ioerr )
 !  S/R Solve
 !
 !  This subroutine solves the system Ax = b using the preconditioned conjugate gradient method.
 !  The matrix action, residual, preconditioning routines are supplied by a non-abstracted 
 !  type-extension of PCConjugateGradient. These routines should return an array indexed from 1 to 
 !  nDOF. Thus, in addition to a MatrixAction and Residual, the user should map their data-structure
 !  to a 1-D array.  
 !
 !  On output ioerr is set to an error checking flag. 
 !  If ioerr ==  0, the method converged within the maximum number of iterations.
 !     ioerr == -1, the method did not converge within the maximum number of iterations.
 !     ioerr == -2, something that is not caught by the current construct happened.
 !  
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS( CGsemElliptic_2D ), INTENT(inout) :: myCGSEM
   REAL(prec), EXTERNAL                     :: dirF
   INTEGER, INTENT(out)                     :: ioerr
   ! LOCAL
   INTEGER    :: nS, nP, nEl,pcioerr
   INTEGER    :: iter
   INTEGER    :: nIt
   REAL(prec) :: TOL
   REAL(prec) :: u(0:myCGSEM % nS, 0:myCGSEM % nP, 1:myCGSEM % mesh % nElems)
   REAL(prec) :: r(0:myCGSEM % nS, 0:myCGSEM % nP, 1:myCGSEM % mesh % nElems)
   REAL(prec) :: v(0:myCGSEM % nS, 0:myCGSEM % nP, 1:myCGSEM % mesh % nElems)
   REAL(prec) :: z(0:myCGSEM % nS, 0:myCGSEM % nP, 1:myCGSEM % mesh % nElems)
   REAL(prec) :: w(0:myCGSEM % nS, 0:myCGSEM % nP, 1:myCGSEM % mesh % nElems)
   REAL(prec) :: rNorm
   REAL(prec) :: a, b, num, den, r0
   
      nS  = myCGSEM % nS
      nP  = myCGSEM % nP
      nEl = myCGSEM % mesh % nElems

      ioerr = -2
      nIt = myCGSEM % maxIters
      TOL = myCGSEM % tolerance

      u = myCGSEM % solution
      CALL myCGSEM % Residual( dirF, r )
      r0 = sqrt( myCGSEM % DotProduct( r, r ) ) 
      
      myCGSEM % resi = ZERO
      myCGSEM % resi(0) = r0
      rNorm = r0
      IF( r0 < TOL )THEN
         RETURN
      ENDIF
 
      DO iter = 1,nIt ! Loop over the CG iterates

         CALL myCGSEM % UnMask( r )
         CALL myCGSEM % preconditioner % Solve( w, r, nS, nP, nEl, pcioerr )
        ! w = r
   
         num = myCGSEM % DotProduct( r, w )

         IF( iter == 1) THEN
            v = w
         ELSE
            b = num/den
            v = w + b*v
         ENDIF

         CALL myCGSEM % MatrixAction( v, z )
         a = num/myCGSEM % DotProduct( v, z )
         u = u + a*v
         r = r - a*z
         
         den = num

         rNorm = sqrt( myCGSEM % DotProduct(r,r) )
         myCGSEM % resi(iter) = rNorm
         IF( rNorm/r0 < TOL ) then
           EXIT
           ioerr = 0
         ENDIF

      ENDDO ! iter, loop over the CG iterates 

      IF( rNorm/r0 > TOL )THEN
         PRINT*, 'MODULE IterativeSolvers : ConjugateGradient failed to converge '
         PRINT*, 'Last L-2 residual : ', rNorm
         ioerr=-1
      ENDIF   

      CALL myCGSEM % UnMask( u )
      myCGSEM % solution = u
      CALL myCGSEM % SetDirichletBoundaryCondition( dirF )

 END SUBROUTINE Solve_PCConjugateGradient 
!
!
!==================================================================================================!
!---------------------------------------- File I/O ------------------------------------------------!
!==================================================================================================!
!
!
SUBROUTINE CoarseToFine_CGsemElliptic_2D( myCGSEM, iEl, x, y, sol, source )
 ! CoarseToFine
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS( CGsemElliptic_2D ), INTENT(in) :: myCGSEM
   INTEGER, INTENT(in)                :: iEl
   REAL(prec), INTENT(out)            :: x(0:myCGSEM % nPlot, 0:myCGSEM % nPlot)
   REAL(prec), INTENT(out)            :: y(0:myCGSEM % nPlot, 0:myCGSEM % nPlot)
   REAL(prec), INTENT(out)            :: sol(0:myCGSEM % nPlot, 0:myCGSEM % nPlot)
   REAL(prec), INTENT(out)            :: source(0:myCGSEM % nPlot, 0:myCGSEM % nPlot)
   ! Local
   REAL(prec) :: localX(0:myCGSEM % nS, 0:myCGSEM % nP)
   REAL(prec) :: localY(0:myCGSEM % nS, 0:myCGSEM % nP)
   
      CALL myCGSEM % mesh % GetPositions( iEl, localX, localY )
      
      CALL myCGSEM % SpectralOps % interp % CoarseToFine( localX, &
                                                        myCGSEM % plMatS, &
                                                        myCGSEM % plMatP, &
                                                        myCGSEM % nPlot, &
                                                        myCGSEM % nPlot, &
                                                        x )
      
      CALL myCGSEM % SpectralOps % interp % CoarseToFine( localY, &
                                                        myCGSEM % plMatS, &
                                                        myCGSEM % plMatP, &
                                                        myCGSEM % nPlot, &
                                                        myCGSEM % nPlot, &
                                                        y )
                                                                                       
      
      CALL myCGSEM % SpectralOps % interp % CoarseToFine( myCGSEM % solution(:,:,iEl), &
                                                        myCGSEM % plMatS, &
                                                        myCGSEM % plMatP, &
                                                        myCGSEM % nPlot, &
                                                        myCGSEM % nPlot, &
                                                        sol )
                                                      
      CALL myCGSEM % SpectralOps % interp % CoarseToFine( myCGSEM % source(:,:,iEl), &
                                                        myCGSEM % plMatS, &
                                                        myCGSEM % plMatP, &
                                                        myCGSEM % nPlot, &
                                                        myCGSEM % nPlot, &
                                                        source )
 END SUBROUTINE CoarseToFine_CGsemElliptic_2D
!
!
!
 SUBROUTINE WriteTecplot_CGsemElliptic_2D( myCGSEM, filename )
 ! S/R WriteTecplot
 !
 !
 ! =============================================================================================== !
   IMPLICIT NONE
   CLASS( CGsemElliptic_2D ), INTENT(inout) :: myCGSEM
   CHARACTER(*), INTENT(in), OPTIONAL       :: filename
   !LOCAL
   INTEGER :: iX, iY, iEl, nS, nEl, fUnit
   CHARACTER(len=5) :: zoneID
   REAL(prec)       :: x(0:myCGSEM % nPlot, 0:myCGSEM % nPlot)
   REAL(prec)       :: y(0:myCGSEM % nPlot, 0:myCGSEM % nPlot)
   REAL(prec)       :: sol(0:myCGSEM % nPlot, 0:myCGSEM % nPlot)
   REAL(prec)       :: source(0:myCGSEM % nPlot, 0:myCGSEM % nPlot)
  
      nEl = myCGSEM % mesh % nElems
      nS  = myCGSEM % nPlot

      IF( PRESENT(filename) )THEN
         OPEN( UNIT=NewUnit(fUnit), &
               FILE=TRIM(filename)//'.tec', &
               FORM='FORMATTED')
      ELSE
         OPEN( UNIT=NewUnit(fUnit), &
               FILE='CGsemSol.tec', &
               FORM='FORMATTED')
      ENDIF
      
      WRITE(fUnit,*) 'VARIABLES = "X", "Y", "sol", "source"'
    
      DO iEl = 1, myCGSEM % mesh % nElems

         WRITE(zoneID,'(I5.5)') iEl
         WRITE(fUnit,*) 'ZONE T="el'//trim(zoneID)//'", I=',nS+1,', J=', nS+1,',F=POINT'
         CALL myCGSEM % CoarseToFine( iEl, x, y, sol, source )

         DO iY = 0, nS
            DO iX = 0, nS
               WRITE(fUnit,*)  x(iX, iY), y(iX, iY), sol(iX,iY), source(iX,iY) 
            ENDDO
         ENDDO

      ENDDO
      
      CLOSE( fUnit )

 END SUBROUTINE WriteTecplot_CGsemElliptic_2D
!
!
!
 SUBROUTINE WriteResidual_CGsemElliptic_2D( myCGSEM, rFile )
 ! S/R WriteResidual
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS( CGsemElliptic_2D ), INTENT(in) :: myCGSEM
   CHARACTER(*), INTENT(in), OPTIONAL    :: rFile
   ! LOCAL
   INTEGER        :: n, i, fUnit
   CHARACTER(100) :: localFile
    
    n = myCGSEM % maxIters
    
    IF( PRESENT(rFile) )THEN
       localFile = rFile
    ELSE
       localFile = 'residual.curve'
    ENDIF
    
    OPEN( UNIT = NewUnit(fUnit), &
          FILE = TRIM(localFile), &
          FORM = 'FORMATTED' )
    WRITE( fUnit, * ) '#residual'
    DO i = 0, n
       WRITE(fUnit,'(I4,1x,E17.8)') i, myCGSEM % resi(i)
    ENDDO

    CLOSE(UNIT = fUnit)   
   
 END SUBROUTINE WriteResidual_CGsemElliptic_2D
!
!
!
 SUBROUTINE ConstructMatrix( myPC )
 !
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS( CGsemElliptic_2D ), INTENT(inout) :: myPC
   ! LOCAL
   INTEGER                 :: iEl, iS, iP, jEl, jS, jP, nEl, nS, nP
   INTEGER                 :: row, col, nfree, fUnit
   REAL(prec)              :: ei(0:myPC % nS, 0:myPC % nP, 1:myPC % mesh % nElems)
   REAL(prec)              :: Aei(0:myPC % nS, 0:myPC % nP, 1:myPC % mesh % nElems)
   REAL(prec)              :: thismask(0:myPC % nS, 0:myPC % nP, 1:myPC % mesh % nElems)
   REAL(prec)              :: s, p
   REAL(prec), ALLOCATABLE :: A(:,:)
   CHARACTER(5)            :: zoneID 

      nS  = myPC % nS
      nP  = myPC % nP
      nEl = myPC % mesh % nElems

      nfree = 0 
      thismask = ONE
      CALL myPC % Mask( thismask )
      DO iEl = 1, nEl
         DO iP = 0, nP
            DO iS = 0, nS
               nfree = nfree + thismask(iS,iP,iEl)
            ENDDO
         ENDDO
      ENDDO

      ALLOCATE(A(1:nfree,1:nfree))

      row = 0
      DO iEl = 1, nEl
         DO iP = 0, nP
            DO iS = 0, nS

               IF( AlmostEqual(thismask(iS,iP,iEl),ONE) )THEN
                  ei = ZERO
                  ei(iS,iP,iEl) = ONE
                  row = row + 1
                  
                  CALL myPC % MatrixAction( ei, Aei )

                  col = 0
                  DO jEl = 1, nEl
                     DO jP = 0, nP
                        DO jS = 0, nS
                           
                           IF( AlmostEqual( thisMask(jS,jP,jEl),ONE ) )THEN
                              col = col + 1
                              A(row,col) = Aei(jS,jP,jEl)
                           ENDIF
                   
                        ENDDO
                     ENDDO
                  ENDDO
  
               ENDIF

            ENDDO
         ENDDO
      ENDDO


      OPEN( UNIT = NewUnit(fUnit), &
            FILE = 'Omat.txt', &
            FORM = 'FORMATTED' )

      DO row = 1, nfree
         WRITE(fUnit,*) A(row,:)
      ENDDO

      CLOSE(fUnit)


           OPEN( UNIT=NewUnit(fUnit), &
            FILE='omask.tec', &
            FORM='FORMATTED')

      
      WRITE(fUnit,*) 'VARIABLES = "X", "Y", "mask"'
    
      DO iEl = 1, myPC % mesh % nElems

         WRITE(zoneID,'(I5.5)') iEl
         WRITE(fUnit,*) 'ZONE T="el'//trim(zoneID)//'", I=',nS+1,', J=', nP+1,',F=POINT'

         DO iP = 0, nS
            DO iS = 0, nS
               CALL myPC % mesh % GetPositionAtNode( iEl, s, p, iS, iP )
               WRITE(fUnit,*)  s, p, thismask(iS,iP,iEl)
            ENDDO
         ENDDO

      ENDDO
      
      CLOSE( fUnit )

               

 END SUBROUTINE ConstructMatrix


END MODULE CGsemElliptic_2D_Class
