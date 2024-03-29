! Stommel_Class.f90
! 
! Copyright 2015 Joseph Schoonover <schoonover.numerics@gmail.com>
! 
! Stommel_Class.f90 is part of the Spectral Element Libraries in Fortran (SELF).
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
 
 

MODULE Stommel_Class
! ========================================= Logs ================================================= !
!2016-05-13  Joseph Schoonover  schoonover.numerics@gmail.com 
!
! //////////////////////////////////////////////////////////////////////////////////////////////// ! 

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
USE StommelParams_Class
USE CGsemPreconditioner_2D_Class


IMPLICIT NONE

      TYPE Stommel
         INTEGER                        :: maxIters, mInnerIters, nPlot, nS, nP
         REAL(prec)                     :: tolerance
         REAL(prec)                     :: Cd
      !   TYPE( CGsemPreconditioner_2D ) :: preconditioner
         TYPE( NodalStorage_2D )        :: SpectralOps
         TYPE( QuadMesh )               :: mesh     
         REAL(prec), ALLOCATABLE        :: solution(:,:,:)
         REAL(prec), ALLOCATABLE        :: source(:,:,:)
         REAL(prec), ALLOCATABLE        :: fluxCoeff(:,:,:)
         REAL(prec), ALLOCATABLE        :: fCori(:,:,:), h(:,:,:)
         REAL(prec), ALLOCATABLE        :: Q(:,:,:)
         REAL(prec), ALLOCATABLE        :: resi(:)
         REAL(prec), ALLOCATABLE        :: plMatS(:,:), plMatP(:,:)

         CONTAINS

         PROCEDURE :: Build                => Build_Stommel
         PROCEDURE :: BuildQuadMesh        => BuildQuadMesh_Stommel
         PROCEDURE :: Trash                => Trash_Stommel
         PROCEDURE :: ResetFluxCoefficient => ResetFluxCoefficient_Stommel
         

         PROCEDURE :: Mask                          => Mask_Stommel
         PROCEDURE :: UnMask                        => UnMask_Stommel
         PROCEDURE :: GlobalSum                     => GlobalSum_Stommel
         PROCEDURE :: FluxDivergence                => FluxDivergence_Stommel
         PROCEDURE :: SetThisDirichletEdge          => SetThisDirichletEdge_Stommel
         PROCEDURE :: SetDirichletBoundaryCondition => SetDirichletBoundaryCondition_Stommel
         PROCEDURE :: CalculateGradient             => CalculateGradient_Stommel
         PROCEDURE :: MatrixAction                  => MatrixAction_Stommel
         PROCEDURE :: Residual                      => Residual_Stommel
         PROCEDURE :: DotProduct                    => DotProduct_Stommel
         PROCEDURE :: Solve                         => Solve_GMRES
         
         PROCEDURE :: CoarseToFine  => CoarseToFine_Stommel
         PROCEDURE :: WriteTecplot  => WriteTecplot_Stommel
         PROCEDURE :: WriteResidual => WriteResidual_Stommel

      END TYPE Stommel

      

CONTAINS

 SUBROUTINE Build_Stommel( myCGSEM, params )
 ! S/R Build
 !
 ! 
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS( Stommel ), INTENT(inout)  :: myCGSEM
   TYPE( StommelParams ), INTENT(in) :: params
   ! LOCAL
   INTEGER                 :: nS, nP, nPlot
   REAL(prec), ALLOCATABLE :: sNew(:)

      nS = params % polyDeg
      nP = nS
      
      nPlot = params % nPlot
      myCGSEM % tolerance   = params % tolerance
      myCGSEM % maxIters    = params % maximumIterates
      myCGSEM % mInnerIters = params % mInnerIters
      myCGSEM % nS          = nS
      myCGSEM % nP          = nP 
      myCGSEM % nPlot       = nPlot
      myCGSEM % Cd          = params % cDrag

      CALL myCGSEM % SpectralOps % Build( nS, nP, GAUSS_LOBATTO, CG )
      CALL myCGSEM % BuildQuadMesh( params )
      
      ! Allocate space for the solution, source term, and diffusivity coefficient
      ALLOCATE( myCGSEM % solution(0:nS, 0:nP, 1:myCGSEM % mesh % nElems), &
                myCGSEM % source(0:nS, 0:nP, 1:myCGSEM % mesh % nElems), &
                myCGSEM % fluxCoeff(0:nS, 0:nP, 1:myCGSEM % mesh % nElems), &
                myCGSEM % fCori(0:nS, 0:nP, 1:myCGSEM % mesh % nElems), &
                myCGSEM % h(0:nS, 0:nP, 1:myCGSEM % mesh % nElems), &
                myCGSEM % Q(0:nS, 0:nP, 1:myCGSEM % mesh % nElems) )

      ALLOCATE( myCGSEM % resi(0:myCGSEM % maxIters*myCGSEM % mInnerIters) )

      myCGSEM % solution  = ZERO
      myCGSEM % source    = ZERO      ! f*wEk --> vortex tube stretching
      myCGSEM % fluxCoeff = ONE
      myCGSEM % h         = ONE
      myCGSEM % fCori     = ZERO
      myCGSEM % Q         = ZERO

      ALLOCATE( sNew(0:nPlot) )
      sNew = UniformPoints( -ONE, ONE, nPlot )
      ALLOCATE( myCGSEM % plMatS(0:nPlot,0:nS), myCGSEM % plMatP(0:nPlot,0:nP) )
      ! Build the plotting matrix
      CALL myCGSEM % SpectralOps % interp % CalculateInterpolationMatrix( nPlot, nPlot, sNew, sNew, &
                                                                          myCGSEM % plMatS, &
                                                                          myCGSEM % plMatP )
      DEALLOCATE( sNew )
      

     !CALL myCGSEM % preconditioner % Build( params % pcTolerance, &
     !                                       params % pcIterates, &
     !                                       myCGSEM % mesh )

     
      
 END SUBROUTINE Build_Stommel
!
!
!
 SUBROUTINE BuildQuadMesh_Stommel( myCGSEM, params )
 ! S/R BuildQuadMesh
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS( Stommel ), INTENT(inout)   :: myCGSEM
   TYPE( StommelParams ), INTENT(in) :: params
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

      CALL myCGSEM % mesh % ScaleTheMesh( myCGSEM % SpectralOps % interp, &
                                          params % xScale, params % yScale  )

 END SUBROUTINE BuildQuadMesh_Stommel
!
!
!
 SUBROUTINE Trash_Stommel( myCGSEM )
 ! S/R Trash
 !
 ! 
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS( Stommel ), INTENT(inout) :: myCGSEM


      ! Trash the nodal storage structure
      CALL myCGSEM % SpectralOps % Trash( )

      ! Trash the geometry
      CALL myCGSEM % mesh % Trash( )

      ! Deallocate space for the solution, source term, and diffusivity coefficient
      DEALLOCATE( myCGSEM % fluxCoeff, &
                  myCGSEM % solution,&
                  myCGSEM % source, &
                  myCGSEM % resi, &
                  myCGSEM % fCori, &
                  myCGSEM % h, &
                  myCGSEM % Q )
                  
      DEALLOCATE( myCGSEM % plMatS, myCGSEM % plMatP )

      !CALL myCGSEM % preconditioner % Trash( )

 END SUBROUTINE Trash_Stommel
!
!
!
 SUBROUTINE ResetFluxCoefficient_Stommel( myCGSEM ) 
 ! S/R ResetFluxCoefficient
 !
 ! =============================================================================================== !
   IMPLICIT NONE
   CLASS( Stommel ), INTENT( inout ) :: myCGSEM 

      myCGSEM % fluxCoeff = myCGSEM % Cd/(myCGSEM % h**2)
     ! CALL myCGSEM % preconditioner % MapFromOrigToPC( myCGSEM % fluxCoeff, &
     !                                                  myCGSEM % preconditioner % fluxCoeff, &
     !                                                  myCGSEM % nS, myCGSEM % nP, &
     !                                                  myCGSEM % mesh % nElems )
 END SUBROUTINE ResetFluxCoefficient_Stommel
!
!
!==================================================================================================!
!-------------------------------------- Type Specific ---------------------------------------------!
!==================================================================================================!
!
!
 SUBROUTINE SetDirichletBoundaryCondition_Stommel( myCGSEM, dirF )
 ! S/R SetDirichletBoundaryCondition
 !
 !
 ! ============================================================================= !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS( Stommel ), INTENT(inout) :: myCGSEM
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


 END SUBROUTINE SetDirichletBoundaryCondition_Stommel
!
!
!
 SUBROUTINE SetThisDirichletEdge_Stommel( myCGSEM, eID, sID, dirF )
 ! S/R SetThisDirichletEdge
 !
 !  This subroutine sets the solutions on side "sID" of element "eID",
 !  excluding the corner nodes ( Assumes Gauss-Lobatto points are used )
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS( Stommel ), INTENT(inout) :: myCGSEM
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

 END SUBROUTINE SetThisDirichletEdge_Stommel
!
!
!
 SUBROUTINE Mask_Stommel( myCGSEM, u )
 ! S/R Mask_Stommel
 !
 !
 ! ============================================================================= !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS( Stommel ), INTENT(inout) :: myCGSEM
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

            IF( e2 > 0 )THEN ! this is an internal edge, and we should mask out the secondary element
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


 END SUBROUTINE Mask_Stommel
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
 SUBROUTINE UnMask_Stommel( myCGSEM, u )
 ! S/R UnMask
 !
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS( Stommel ), INTENT(inout) :: myCGSEM
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

 END SUBROUTINE UnMask_Stommel
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
 SUBROUTINE GlobalSum_Stommel( myCGSEM, u )
 ! S/R GlobalSum
 !
 !  Adds together the shared edge and node contributions for the CGSEM method.
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS( Stommel ), INTENT(inout) :: myCGSEM
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

 END SUBROUTINE GlobalSum_Stommel
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
 SUBROUTINE CalculateGradient_Stommel( myCGSEM, iEl, u, dudx, dudy )
 ! S/R CalculateGradient
 !
 !  Calculates the gradient of the solution in computational coordinates within 
 !  a single element
 !
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS( Stommel ), INTENT(in) :: myCGSEM 
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

 END SUBROUTINE CalculateGradient_Stommel
!
!
!
 SUBROUTINE FluxDivergence_Stommel( myCGSEM, iEl, u, Lu )
 ! S/R FluxDivergence
 !
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS( Stommel ), INTENT(inout) :: myCGSEM
   INTEGER, INTENT(in)             :: iEl
   REAL(prec), INTENT(in)          :: u(0:myCGSEM % nS, 0:myCGSEM % nP)
   REAL(prec), INTENT(out)         :: Lu(0:myCGSEM % nS, 0:myCGSEM % nP)
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
   REAL(prec) :: dxds, dxdp, dyds, dydp, Fx, Fy

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
!            F1(iS,iP) = myCGSEM % fluxCoeff(iS,iP,iEl)*( dudx(iS,iP)*dydp - dudy(iS,iP)*dxdp )*ws(iS)
!            F2(iS,iP) = myCGSEM % fluxCoeff(iS,iP,iEl)*( dudy(iS,iP)*dxds - dudx(iS,iP)*dyds )*wp(iP)  
            Fx = myCGSEM % fluxCoeff(iS,iP,iEl)*dudx(iS,iP) - &
                 myCGSEM % Q(iS,iP,iEl)*dudy(iS,iP)
            
            Fy = myCGSEM % fluxCoeff(iS,iP,iEl)*dudy(iS,iP) + &
                 myCGSEM % Q(iS,iP,iEl)*dudx(iS,iP)
                 
            F1(iS,iP) = myCGSEM % fluxCoeff(iS,iP,iEl)*( Fx*dydp - Fy*dxdp )*ws(iS)
            F2(iS,iP) = myCGSEM % fluxCoeff(iS,iP,iEl)*( Fy*dxds - Fx*dyds )*wp(iP) 

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
       

 END SUBROUTINE FluxDivergence_Stommel 
!
!
!
 SUBROUTINE MatrixAction_Stommel( myCGSEM, u, Au )
 ! S/R MatrixAction
 ! 
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS(Stommel), INTENT(inout) :: myCGSEM
   REAL(prec), INTENT(inout)     :: u( 0:myCGSEM % nS, &
                                       0:myCGSEM % nP, &
                                       1:myCGSEM % mesh % nElems )
   REAL(prec), INTENT(out)       :: Au(0:myCGSEM % nS, &
                                       0:myCGSEM % nP, &
                                       1:myCGSEM % mesh % nElems  )
   ! Local
   INTEGER :: iEl, iS, iP, nS, nP, nEl
   REAL(prec) :: temp(0:myCGSEM % nS, &
                      0:myCGSEM % nP )
   REAL(prec) :: Atemp(0:myCGSEM % nS, &
                       0:myCGSEM % nP )
   REAL(prec) :: J, qX, qY
   REAL(prec) :: ws(0:myCGSEM % nS)
   REAL(prec) :: wp(0:myCGSEM % nP)
   REAL(prec) :: ux(0:myCGSEM % nS, 0:myCGSEM % nP)
   REAL(prec) :: uy(0:myCGSEM % nS, 0:myCGSEM % nP)


      nS  = myCGSEM % nS
      nP  = myCGSEM % nP
      nEl = myCGSEM % mesh % nElems
      CALL myCGSEM % UnMask( u ) ! unmask the solution
      CALL myCGSEM % SpectralOps % GetQuadratureWeights( ws, wp )

!$OMP PARALLEL PRIVATE( temp, Atemp,J, locR, bX, bY, ux, uy )
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

 END SUBROUTINE MatrixAction_Stommel
!
!
!
 SUBROUTINE Residual_Stommel( myCGSEM, dirF, r )
 ! S/R Residual
 ! 
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS( Stommel ), INTENT(inout) :: myCGSEM
   REAL(prec), EXTERNAL            :: dirF
   REAL(prec), INTENT(out)         :: r(0:myCGSEM % nS, &
                                        0:myCGSEM % nP, &
                                        1:myCGSEM % mesh % nElems )
   ! LOCAL
   INTEGER    :: iS, iP, nS, nP, nEl, iEl, nEdges
   REAL(prec) :: J!, qX, qY
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

!$OMP PARALLEL PRIVATE( J )
!$OMP DO      
      DO iEl = 1, nEl

         DO iP = 0, nP
            DO iS = 0, nS
               CALL myCGSEM % mesh % GetJacobianAtNode( iEl, J, iS, iP )

               r(iS,iP,iEl) = myCGSEM % source(iS,iP,iEl)*J*ws(iS)*wp(iP)

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
     r = r - Au

     ! Mask the edges in preparation for residual calculation
     CALL myCGSEM % Mask( r )

     ! Re-mask the solution
     CALL myCGSEM % Mask( myCGSEM % solution )

     ! re-set the dirichlet boundary conditions
     CALL myCGSEM % SetDirichletBoundaryCondition( dirF )

 END SUBROUTINE Residual_Stommel
!
!
!
 FUNCTION DotProduct_Stommel( myCGSEM, u, v ) RESULT( uDotv )
 ! FUNCTION DotProduct
 ! Computes the dot product of two "vectors" whose data are
 ! stored in the "(iS,iP)" format
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS( Stommel ) :: myCGSEM
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

 END FUNCTION DotProduct_Stommel
!
!
!
 SUBROUTINE Solve_GMRES( myCGSEM, dirF, ioerr )
 !  S/R Solve
 !
 !  This subroutine solves the system Ax = b using the un-preconditioned GMRES.
 !  The matrix action and residual routines are supplied by a non-abstracted type-extension of
 !  GMRES. These routines should return an array indexed from 1 to nDOF. Thus,
 !  in addition to a MatrixAction and Residual, the user should map their data-structure to a 1-D 
 !  array.  
 !
 !  On output ioerr is set to an error checking flag. 
 !  If ioerr ==  0, the method converged within the maximum number of iterations.
 !     ioerr == -1, the method did not converge within the maximum number of iterations.
 !     ioerr == -2, something that is not caught by the current construct happened.
 !  
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS( Stommel ), INTENT(inout) :: myCGSEM
   REAL(prec), EXTERNAL            :: dirF
   INTEGER, INTENT(out)            :: ioerr
   ! LOCAL
   INTEGER    :: i, j, k, l, iS, iP, iEl
   INTEGER    :: nIt, m, nr 
   REAL(prec) :: TOL
   REAL(prec) :: r(0:myCGSEM % nS, 0:myCGSEM % nP, 1:myCGSEM % mesh % nElems)
   REAL(prec) :: v(0:myCGSEM % nS, 0:myCGSEM % nP, 1:myCGSEM % mesh % nElems, 1:myCGSEM % mInnerIters+1)
   REAL(prec) :: w(0:myCGSEM % nS, 0:myCGSEM % nP, 1:myCGSEM % mesh % nElems)
   REAL(prec) :: rho(1:myCGSEM % mInnerIters, 1:myCGSEM % mInnerIters)
   REAL(prec) :: h(1:myCGSEM % mInnerIters+1, 1:myCGSEM % mInnerIters)
   REAL(prec) :: c(1:myCGSEM % mInnerIters), s(1:myCGSEM % mInnerIters), bhat(1:myCGSEM % mInnerIters+1)
   REAL(prec) :: y(1:myCGSEM % mInnerIters+1)
   REAL(prec) :: b, d, g, r0, rc
   
      ioerr = -2
      nIt = myCGSEM % maxIters
      m   = myCGSEM % mInnerIters
      TOL = myCGSEM % tolerance
      r = ZERO

      CALL myCGSEM % Residual( dirF, r )
      r0 = sqrt(myCGSEM % DotProduct( r, r ) )
     
      l = 0 
      myCGSEM % resi = ZERO
      myCGSEM % resi(l) = r0
      
      v    = ZERO  
      rho  = ZERO
      bhat = ZERO 
      s    = ZERO
      c    = ZERO
      y    = ZERO
   
     
      DO j = 1,nIt

         b       = sqrt( myCGSEM % DotProduct( r, r ) )
         v(:,:,:,1)  = r/b
         bhat(1) = b

         DO i = 1, m
            l = l+1
            nr = i
            ! The first step in GMRES is to build the orthogonal basis up to order "i"
            ! with the accompanying upper hessenburg matrix.
            CALL myCGSEM % MatrixAction( v(:,:,:,i), w )
            
            ! The new basis vector is obtained by multiplying the previous basis vector by the matrix
            ! and orthogonalizing wrt to all of the previous basis vectors using a Gram-Schmidt process.
            DO k = 1, i
               h(k,i) = myCGSEM % DotProduct( v(:,:,:,k), w )
               w      = w - h(k,i)*v(:,:,:,k)
            ENDDO

            h(i+1,i) = sqrt( myCGSEM % DotProduct(w,w) )

            IF( AlmostEqual( h(i+1,i), ZERO )  )THEN
               EXIT
            ENDIF

            v(:,:,:,i+1) = w/h(i+1,i)
            rho(1,i) = h(1,i)

            ! Givens rotations are applied to the upper hessenburg matrix and to the residual vectors
            ! that are formed from the orthogonalization process. Here, they are done "on-the-fly"
            ! as opposed to building the entire upper hessenburg matrix and orthonormal basis
            ! before performing the rotations. This way, we can also tell if we have found an exact
            ! solution ( if h(i+1,i) = 0 ) with a smaller subspace than size m.
            DO k = 2, i

               g          = c(k-1)*rho(k-1,i) + s(k-1)*h(k,i)
               rho(k,i)   = -s(k-1)*rho(k-1,i) + c(k-1)*h(k,i)
               rho(k-1,i) = g 

            ENDDO

            ! Here the coefficients of the Givens rotation matrix are computed
            d = sqrt( rho(i,i)**2 + h(i+1,i)**2 )
            c(i) = rho(i,i)/d
            s(i) = h(i+1,i)/d

            rho(i,i) = c(i)*rho(i,i) + s(i)*h(i+1,i)
            ! And applied to the residual vector
            bhat(i+1) = -s(i)*bhat(i)
            bhat(i)   = c(i)*bhat(i)
 
         
            rc = abs( bhat(i+1) )

            myCGSEM % resi(l) = rc
            
            IF( rc/r0 <= TOL )THEN
               EXIT
            ENDIF

         ENDDO
         
         IF( rc/r0 > TOL )THEN
            nr = m
         ENDIF

         ! Back-substitution of the tridiagonal matrix that resulted from the rotations
         y(nr) = bhat(nr)/rho(nr,nr)
         DO k = nr-1, 1, -1

            y(k) = bhat(k)

            DO i = k+1, nr
               y(k) = y(k) - rho(k,i)*y(i)

            ENDDO

            y(k) = y(k)/rho(k,k)

         ENDDO
         
 
         DO iEl = 1,myCGSEM % mesh % nElems
            DO iP = 0, myCGSEM % nP
               DO iS = 0, myCGSEM % nS
                  myCGSEM % solution(iS,iP,iEl) = myCGSEM % solution(iS,iP,iEl) + &
                                                 DOT_PRODUCT( v(iS,iP,iEl,1:nr), y(1:nr) )
               ENDDO
            ENDDO
         ENDDO

         IF( rc/r0 <= TOL )THEN
            ioerr = 0
            EXIT
         ENDIF

         CALL myCGSEM % Residual( dirF, r )

      ENDDO 

      IF( rc/r0 > TOL )THEN
         PRINT*, 'MODULE IterativeSolvers : GMRES failed to converge '
         PRINT*, 'Last L-2 residual : ', sqrt(myCGSEM % DotProduct(r,r))
         ioerr=-1
      ENDIF   

 END SUBROUTINE Solve_GMRES
!
!
!==================================================================================================!
!---------------------------------------- File I/O ------------------------------------------------!
!==================================================================================================!
!
!
SUBROUTINE CoarseToFine_Stommel( myCGSEM, iEl, x, y, sol, source, depth, coriolis, u, v )
 ! CoarseToFine
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS( Stommel ), INTENT(in) :: myCGSEM
   INTEGER, INTENT(in)          :: iEl
   REAL(prec), INTENT(out)      :: x(0:myCGSEM % nPlot, 0:myCGSEM % nPlot)
   REAL(prec), INTENT(out)      :: y(0:myCGSEM % nPlot, 0:myCGSEM % nPlot)
   REAL(prec), INTENT(out)      :: sol(0:myCGSEM % nPlot, 0:myCGSEM % nPlot)
   REAL(prec), INTENT(out)      :: source(0:myCGSEM % nPlot, 0:myCGSEM % nPlot)
   REAL(prec), INTENT(out)      :: depth(0:myCGSEM % nPlot, 0:myCGSEM % nPlot)
   REAL(prec), INTENT(out)      :: coriolis(0:myCGSEM % nPlot, 0:myCGSEM % nPlot)
   REAL(prec), INTENT(out)      :: u(0:myCGSEM % nPlot, 0:myCGSEM % nPlot)
   REAL(prec), INTENT(out)      :: v(0:myCGSEM % nPlot, 0:myCGSEM % nPlot)

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

      CALL myCGSEM % SpectralOps % interp % CoarseToFine( myCGSEM % h(:,:,iEl), &
                                                        myCGSEM % plMatS, &
                                                        myCGSEM % plMatP, &
                                                        myCGSEM % nPlot, &
                                                        myCGSEM % nPlot, &
                                                        depth )

      CALL myCGSEM % SpectralOps % interp % CoarseToFine( myCGSEM % fCori(:,:,iEl), &
                                                        myCGSEM % plMatS, &
                                                        myCGSEM % plMatP, &
                                                        myCGSEM % nPlot, &
                                                        myCGSEM % nPlot, &
                                                        coriolis )

      CALL myCGSEM % CalculateGradient( iEl, myCGSEM % solution(:,:,iEl) , localX, localY )

      CALL myCGSEM % SpectralOps % interp % CoarseToFine( -localY, &
                                                          myCGSEM % plMatS, &
                                                          myCGSEM % plMatP, &
                                                          myCGSEM % nPlot, &
                                                          myCGSEM % nPlot, &
                                                          u )

      CALL myCGSEM % SpectralOps % interp % CoarseToFine( localX, &
                                                          myCGSEM % plMatS, &
                                                          myCGSEM % plMatP, &
                                                          myCGSEM % nPlot, &
                                                          myCGSEM % nPlot, &
                                                          v )


 END SUBROUTINE CoarseToFine_Stommel
!
!
!
 SUBROUTINE WriteTecplot_Stommel( myCGSEM, filename )
 ! S/R WriteTecplot
 !
 !
 ! =============================================================================================== !
   IMPLICIT NONE
   CLASS( Stommel ), INTENT(inout) :: myCGSEM
   CHARACTER(*), INTENT(in), OPTIONAL       :: filename
   !LOCAL
   INTEGER :: iX, iY, iEl, nS, nEl, fUnit
   CHARACTER(len=5) :: zoneID
   REAL(prec)       :: x(0:myCGSEM % nPlot, 0:myCGSEM % nPlot)
   REAL(prec)       :: y(0:myCGSEM % nPlot, 0:myCGSEM % nPlot)
   REAL(prec)       :: sol(0:myCGSEM % nPlot, 0:myCGSEM % nPlot)
   REAL(prec)       :: source(0:myCGSEM % nPlot, 0:myCGSEM % nPlot)
   REAL(prec)       :: depth(0:myCGSEM % nPlot, 0:myCGSEM % nPlot)
   REAL(prec)       :: coriolis(0:myCGSEM % nPlot, 0:myCGSEM % nPlot)
   REAL(prec)       :: u(0:myCGSEM % nPlot, 0:myCGSEM % nPlot)
   REAL(prec)       :: v(0:myCGSEM % nPlot, 0:myCGSEM % nPlot)

      nEl = myCGSEM % mesh % nElems
      nS  = myCGSEM % nPlot

      CALL myCGSEM % UnMask( myCGSEM % solution )

      IF( PRESENT(filename) )THEN
         OPEN( UNIT=NewUnit(fUnit), &
               FILE=TRIM(filename)//'.tec', &
               FORM='FORMATTED')
      ELSE
         OPEN( UNIT=NewUnit(fUnit), &
               FILE='CGsemSol.tec', &
               FORM='FORMATTED')
      ENDIF
      
      WRITE(fUnit,*) 'VARIABLES = "X", "Y", "StreamFunction", "U", "V", "Depth", "Coriolis", "source"'
    
      DO iEl = 1, myCGSEM % mesh % nElems

         WRITE(zoneID,'(I5.5)') iEl
         WRITE(fUnit,*) 'ZONE T="el'//trim(zoneID)//'", I=',nS+1,', J=', nS+1,',F=POINT'
         CALL myCGSEM % CoarseToFine( iEl, x, y, sol, source, depth, coriolis, u, v )

         DO iY = 0, nS
            DO iX = 0, nS
               WRITE(fUnit,*)  x(iX, iY), y(iX, iY), sol(iX, iY), &
                               u(iX, iY), v(iX, iY), depth(iX, iY), &
                               coriolis(iX,iY), source(iX,iY)
            ENDDO
         ENDDO

      ENDDO
      
      CLOSE( fUnit )

 END SUBROUTINE WriteTecplot_Stommel
!
!
!
 SUBROUTINE WriteResidual_Stommel( myCGSEM, rFile )
 ! S/R WriteResidual
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS( Stommel ), INTENT(in) :: myCGSEM
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
   
 END SUBROUTINE WriteResidual_Stommel
!
!
!
END MODULE Stommel_Class
