! CGsemElliptic_2D_2D_Class.f90 ( new with v2.1 - 25 March 2016)
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
! o  (ver 1.0) April 2015
! o  (ver 2.1) March 2016
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
USE RunParamsClass
! src/nodal/
USE NodalStorage_2D_Class
! src/geom/
USE QuadMeshClass
USE EdgeClass
USE NodeClass_2D
USE VectorClass
USE MappedGeometryClass_2D



IMPLICIT NONE

      TYPE CGsemElliptic_2D
         TYPE( NodalStorage_2D )       :: cgStorage
         TYPE( QuadMesh )              :: mesh           
         REAL(prec), ALLOCATABLE       :: cgSol(:,:,:)
         REAL(prec), ALLOCATABLE       :: source(:,:,:)
         REAL(prec), ALLOCATABLE       :: fluxCoeff(:,:,:)


         CONTAINS

         PROCEDURE :: Build         => Build_CGsemElliptic_2D
         PROCEDURE :: BuildQuadMesh => BuildQuadMesh_CGsemElliptic_2D
         PROCEDURE :: Trash         => Trash_CGsemElliptic_2D
       
         PROCEDURE :: SetDirichletBoundaryCondition => SetDirichletBoundaryCondition_CGsemElliptic_2D
         PROCEDURE :: SetThisDirichletEdge          => SetThisDirichletEdge_CGsemElliptic_2D

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
   TYPE( RunParams ), INTENT(in)         :: params
   ! LOCAL
   INTEGER :: nS, nP, iS, iP, iEl

      nS = params % polyDeg
      nP = nS

      CALL myCGSEM % cgStorage % Build( nS, nP, GAUSS_LOBATTO, CG )

      ! Allocate space for the solution, source term, and diffusivity coefficient
      ALLOCATE( myCGSEM % solution(0:nS, 0:nP, 1:myCGSEM % mesh % nElems), &
                myCGSEM % source(0:nS, 0:nP, 1:myCGSEM % mesh % nElems), &
                myCGSEM % fluxCoeff(0:nS, 0:nP, 1:myCGSEM % mesh % nElems) )

      myCGSEM % solution   =  ZERO
      myCGSEM % source     = ZERO
      myCGSEM % fluxCoeffs = ZERO

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
   TYPE( RunParams ), INTENT(in)         :: params

      IF( TRIM( params % SpecMeshFile ) == nada )THEN
         PRINT*,' Loading default 2-D mesh.'
         CALL  myCGSEM % mesh % LoadDefaultMesh( myCGSEM % cgStorage % interp, &
                                         params % nXelem, &
                                         params % nYelem )
      ELSE
      ! Builds the lateral mesh
         PRINT*, 'Reading 2-D mesh from '//trim(params % SpecMeshFile)//'.'
         CALL myCGSEM % mesh % ReadSpecMeshFile( myCGSEM % cgStorage % interp, params % SpecMeshFile )
      ENDIF

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
      CALL myCGSEM % cgStorage % Trash( )

      ! Trash the geometry
      CALL myCGSEM % mesh % Trash( )

      ! Deallocate space for the solution, source term, and diffusivity coefficient
      DEALLOCATE( myCGSEM % fluxCoeffs, &
                  myCGSEM % solution,&
                  myCGSEM % source )

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
   CLASS( CGSEM2D ), INTENT(inout) :: myCGSEM
   ! dirF is a function which is passed to this subroutine
   REAL(prec), EXTERNAL :: dirF !(x, y)
   !LOCAL
   INTEGER :: nS, nP, iEdge
   INTEGER :: e1, e2, s1, s2, nodeType, i, j
   real(prec) :: x, y 

      nS = myCGSEM % cgStorage % nS
      nP = myCGSEM % cgStorage % nP
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
            myCGSEM % mesh % nodes(iNode) % nodeToElement % current => myCGSEM % mesh % nodes(iNode) % nodeToElement % head

            ! Loop over the elements in this corner-node's connectivity list
            DO WHILE( ASSOCIATED( myCGSEM % mesh % nodes(iNode) % nodeToElement % current ) ) ! while there are elements in the node-to-element list

               CALL myCGSEM % mesh % nodes(iNode) % nodeToElement % GET_CURRENT_DATA( eID, nID ) ! Get the element ID and the local node ID (1->4)

               i = myCGSEM % mesh % cornerMap(1, nID) ! Grab the index over the first computational direction associated with this corner node
               j = myCGSEM % mesh % cornerMap(2, nID) ! Grab the index over the second computational direction associated with this corner node

               x = myCGSEM % mesh % elements(eID) % geometry % x(i,j)
               y = myCGSEM % mesh % elements(eID) % geometry % y(i,j)

               myCGSEM % solution(i,j,eID) = dirF( x, y )

               CALL myCGSEM % mesh % nodes(iNode) % nodeToElement % MoveToNext( )

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
   INTEGER, INTENT(in)                   :: eID, sID
   REAL(prec), EXTERNAL :: dirF !(x, y)
   ! LOCAL
   INTEGER    :: iS, iP, nS, nP,nElems
   REAL(prec) :: x, y
   

      nS     = myCGSEM % cgStorage % nS
      nP     = nS
      nElems = myCGSEM % mesh % nElems

      IF( sID == east .OR. sID == west )then ! east or west sides

         iS = myCGSEM % mesh % sideMap(sID) 

         DO iP = 1, nP-1 ! Loop over the edge, excluding the corner nodes
            CALL myDGSEM % mesh % GetPositionAtNode( eID, x, y, iS, iP )
            myCGSEM % cgSol(iS,iP,eID) = dirF( x, y )
         ENDDO ! iP, loop over the edge

      ELSE ! south or north sides

         iP = myCGSEM % mesh % sideMap(sID)

         DO iS = 1, nS-1 ! Loop over the edge, excluding the corner nodes
            CALL myDGSEM % mesh % GetPositionAtNode( eID, x, y, iS, iP )
            myCGSEM % cgSol(iS,iP,eID) = dirF( x, y )
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
   REAL(prec), INTENT(inout)             :: u(0:myCGSEM % cgStorage % nS, &
                                              0:myCGSEM % cgStorage % nP, &
                                              1:myCGSEM % mesh % nElems )
   !LOCAL
   INTEGER :: nS, nP, iEdge, eID, sID, iNode, nID, i, j, nEl, nEdges, nNodes
   INTEGER :: e1, e2, s1 s2

      nS     = myCGSEM % cgStorage % nS
      nP     = myCGSEM % cgStorage % nP
      nEl    = myCGSEM % mesh % nElems
      nEdges = myCGSEM % mesh % nEdges
      nNodes = myDGSEM % mesh % nNodes
   
      ! Mask out the solution over the edges -- excluding the corner nodes
      DO iEdge = 1, nEdges ! Loop over the edges
 
         CALL myCGSEM % mesh % GetEdgeSecondaryElementID( iEdge, e2 )

         IF( e2 == DIRICHLET )THEN ! check the secondary element for boundary condition flag

            ! This is a prescribed boundary condition (Dirichlet) so we mask out the solution along this edge   

            CALL myCGSEM % mesh % GetEdgePrimaryElementID( iEdge, e1 )
            CALL myCGSEM % mesh % GetEdgePrimaryElementSide( iEdge, s1 )

            CALL MaskSide( e1, s1, myCGSEM % mesh, u, nS, nP, nElems )

         ELSE ! then this is not a prescribed boundary, and we choose to mask out the secondary element

            IF( eID > 0 )THEN ! this is an internal edge, and we should mask out the secondary element
               CALL myCGSEM % mesh % GetEdgeSecondaryElementID( iEdge, e2 )
               CALL myCGSEM % mesh % GetEdgeSecondaryElementSide( iEdge, s2 )
               CALL MaskSide( e2, s2, myCGSEM % mesh, u, nS, nP, nElems )
            ENDIF

         ENDIF

      ENDDO ! iEdge, loop over the edges     

      ! At this point, the secondary internal edges and prescribed boundary condition edges have been masked out.
      ! Now we mask out all but the primary element in the corner-node lists.

      ! Mask out the corner-nodes appropriately ( Leave the primary element intact )
      DO iNode = 1, nNodes ! Loop over the corner nodes

         ! Rewind to the head of the linked list
         CALL myCGSEM % mesh % nodes(iNode) % nodeToElement % MoveToHead( ) 

         IF( myCGSEM % mesh % nodes(iNode) % nodeType /= DIRICHLET ) then ! this is NOT a prescribed (Dirichlet) boundary
                                                                          ! The data in this non-prescribed node is preserved

            CALL myCGSEM % mesh % nodes(iNode) % nodeToElement % MoveToNext( ) ! skip the first element in the list so it's data is not masked

         ENDIF

         DO WHILE( ASSOCIATED( myCGSEM % mesh % nodes(iNode) % nodeToElement % current ) ) ! while there are elements in the node-to-element list

            CALL myCGSEM % mesh % nodes(iNode) % nodeToElement % GetCurrentData( eID, nID ) ! Get the element ID and the local node ID (1->4)

            i = myCGSEM % mesh % cornerMap(1, nID) ! Grab the index over the first computational direction associated with this corner node
            j = myCGSEM % mesh % cornerMap(2, nID) ! Grab the index over the second computational direction associated with this corner node
 
            u(i,j,eID) = ZERO ! Mask out the solution at this element at this corner node

            CALL myCGSEM % mesh % nodes(iNode) % nodeToElement % MoveToNext( )

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
   INTEGER :: iS, iP
   
      IF( sID == 2 .OR. sID == 4)then ! east or west sides

         iS = mesh % sideMap(sID) 
         DO iP = 1, nP-1 ! Loop over the edge, excluding the corner nodes
            u(iS,iP,eID) = ZERO
         ENDDO ! iP, loop over the edge

      ELSE ! south or north sides

         iP = mesh % sideMap(sID)
         DO iS = 1, nS-1 ! Loop over the edge, excluding the corner nodes
            u(iS,iP,eID) = ZERO
         ENDDO

      ENDIF

 END SUBROUTINE MaskSide
!
!
!
 SUBROUTINE Unmask_CGsemElliptic_2D( myCGSEM, u )
 ! S/R UnMask
 !
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS( CGsemElliptic_2D ), INTENT(inout) :: myCGSEM
   REAL(prec), INTENT(inout)             :: u(0:myCGSEM % cgStorage % nS, &
                                              0:myCGSEM % cgStorage % nP, &
                                              1:myCGSEM % mesh % nElems )
   !LOCAL
   INTEGER :: nS, nP, iEdge, eID, sID, iNode, nID, i, j, eID0, i0, j0, nEl, nEdges, nNodes
   INTEGER :: e1, s1, e2, s2

      nS     = myCGSEM % cgStorage % nS
      nP     = myCGSEM % cgStorage % nP
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


         IF( myCGSEM % mesh % nodes(iNode) % nodeType /= DIRICHLET ) then ! this is NOT a prescribed (Dirichlet) boundary

            ! Rewind to the head of the linked list
            CALL myCGSEM % mesh % nodes(iNode) % nodeToElement % MoveToHead( ) 

            ! Gather the element and local node IDs of the corner which contains the unmasked solution values.
            CALL myCGSEM % mesh % nodes(iNode) % nodeToElement % GetCurrentData( eID0, nID ) ! Get the element ID and the local node ID (1->4)

            i0 = myCGSEM % mesh % cornerMap(1, nID) ! Grab the index over the first computational direction associated with this corner node
            j0 = myCGSEM % mesh % cornerMap(2, nID) ! Grab the index over the second computational direction associated with this corner node

            CALL myCGSEM % mesh % nodes(iNode) % nodeToElement % MoveToNext( ) ! move on to the next element which shares this node


            ! Loop over the elements in this corner-node's connectivity list
            DO WHILE( ASSOCIATED( myCGSEM % mesh % nodes(iNode) % nodeToElement % current ) ) ! while there are elements in the node-to-element list

               CALL myCGSEM % mesh % nodes(iNode) % nodeToElement % GetCurrentData( eID, nID ) ! Get the element ID and the local node ID (1->4)

               i = myCGSEM % mesh % cornerMap(1, nID) ! Grab the index over the first computational direction associated with this corner node
               j = myCGSEM % mesh % cornerMap(2, nID) ! Grab the index over the second computational direction associated with this corner node
  
               u(i,j,eID) = u(i0,j0,eID0) ! copy the solution from the unmasked solution

               CALL myCGSEM % mesh % nodes(iNode) % nodeToElement % MoveToNext( )

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
   INTEGER :: iS, iP, eID, sID, jS, jP
   REAL(prec) :: temp(0:nS) ! assume nS == nP

      ! Get the primary element and side IDs
      eID = mesh % edges(edgeID) % elementIDs(1)
      sID = mesh % edges(edgeID) % elementSides(1)
 
      IF( sID == 2 .OR. sID == 4)THEN ! east or west sides

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
      eID = mesh % edges(edgeID) % elementIDs(2)
      sID = abs( mesh % edges(edgeID) % elementSides(2) )
        
      IF( sID == 2 .OR. sID == 4)THEN ! east or west sides

         iS = mesh % sideMap(sID) 
         jP = mesh % edges(edgeID) % start ! accounting for the possibility that the node ordering may be different in the secondary element

         DO iP = 1, nP-1 ! Loop over the edge, excluding the corner nodes
            u(iS,jP,eID) = temp(iP) ! Copy the solution to the secondary element
            jP = jP + mesh % edges(edgeID) % inc ! increment jP
         ENDDO ! iP, loop over the edge

      ELSE ! south or north sides

         iP = mesh % sideMap(sID)
         jS = mesh % edges(edgeID) % start ! accounting for the possibility that the node ordering may be different in the secondary element

         DO iS = 1, nS-1 ! Loop over the edge, excluding the corner nodes

            u(jS,iP,eID) = temp(iS) ! Copy the solution in the primary element
            jS = jS + mesh % edges(edgeID) % inc ! increment jS

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
   REAL(prec), INTENT(inout)             :: u(0:myCGSEM % cgStorage % nS, &
                                              0:myCGSEM % cgStorage % nP, &
                                              1:myCGSEM % mesh % nElems )
   !LOCAL
   INTEGER    :: iEdge, eID, sID, iNode, nID, i, j, eID0, i0, j0, nEdges, nNodes, nS, nP, nEl
   INTEGER    :: e1, s1, e2, s2 
   REAL(prec) :: theSum

      nS     = myCGSEM % cgStorage % nS
      nP     = myCGSEM % cgStorage % nP
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

         IF( myCGSEM % mesh % nodes(iNode) % nodeType /= DIRICHLET ) then ! this node does NOT lie on a prescribed (Dirichlet) boundary

             ! Rewind to the head of the linked list
            CALL myCGSEM % mesh % nodes(iNode) % nodeToElement % MoveToHead( )

            ! *** First, compute the sum from the contributing elements *** !
            theSum = ZERO 

            ! Loop over the elements in this corner-node's connectivity list
            DO WHILE( ASSOCIATED( myCGSEM % mesh % nodes(iNode) % nodeToElement % current ) ) ! while there are elements in the node-to-element list

               CALL myCGSEM % mesh % nodes(iNode) % nodeToElement % GetCurrentData( eID, nID ) ! Get the element ID and the local node ID (1->4)

               i = myCGSEM % mesh % cornerMap(1, nID) ! Grab the index over the first computational direction associated with this corner node
               j = myCGSEM % mesh % cornerMap(2, nID) ! Grab the index over the second computational direction associated with this corner node
  
               theSum = theSum + u(i,j,eID) 

               CALL myCGSEM % mesh % nodes(iNode) % nodeToElement % MoveToNext( )

            ENDDO ! while we have elements in the node-to-element list

            ! *** Now, copy the sum to the array "u" for each element which contributed to the sum *** !

             ! Rewind to the head of the linked list
            CALL myCGSEM % mesh % nodes(iNode) % nodeToElement % MoveToHead( )

            DO WHILE( ASSOCIATED( myCGSEM % mesh % nodes(iNode) % nodeToElement % current ) ) ! while there are elements in the node-to-element list

               CALL myCGSEM % mesh % nodes(iNode) % nodeToElement % GetCurrentData( eID, nID ) ! Get the element ID and the local node ID (1->4)

               i = myCGSEM % mesh % cornerMap(1, nID) ! Grab the index over the first computational direction associated with this corner node
               j = myCGSEM % mesh % cornerMap(2, nID) ! Grab the index over the second computational direction associated with this corner node
  
               u(i,j,eID) = theSum

               CALL myCGSEM % mesh % nodes(iNode) % nodeToElement % MoveToNext( )

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
   INTEGER :: iS, iP, eID, sID, jS, jP, k, n
   REAL(prec) :: temp(0:nS,1:2) ! assume nS == nP
   REAL(prec) :: theSum


      DO k = 1,2

         ! Get the primary/secondary element and side IDs
         eID = mesh % edges(edgeID) % elementIDs(k)
         sID = abs( mesh % edges(edgeID) % elementSides(k) )
   
 
         IF( sID == east .OR. sID == west )then ! east or west sides

            iS = mesh % sideMap(sID) 
            DO iP = 1, nP-1 ! Loop over the edge, excluding the corner nodes
               temp(iP,k) = u(iS,iP,eID) ! 
            ENDDO ! iP, loop over the edge

         ELSE ! south or north sides

            iP = mesh % sideMap(sID)
            DO iS = 1, nS-1 ! Loop over the edge, excluding the corner nodes
               temp(iS,k) = u(iS,iP,eID) ! 
            ENDDO

         ENDIF

      ENDDO

      ! Add the two contributions together
      n = mesh % edges(edgeID) % start
     
      DO iS = 1, nS-1

         theSum = temp(iS,1) + temp(n,2)

         temp(iS,1) = theSum
         temp(n,2) = theSum

         n = n + mesh % edges(edgeID) % inc ! increment the secondary element array adress

      ENDDO

      ! copy the sum into the two contributing element edges

      DO k = 1,2

         ! Get the primary/secondary element and side IDs
         eID = mesh % edges(edgeID) % elementIDs(k)
         sID = abs( mesh % edges(edgeID) % elementSides(k) )
   
 
         IF( sID == east .OR. sID == west )THEN ! east or west sides

            iS = mesh % sideMap(sID) 
            DO iP = 1, nP-1 ! Loop over the edge, excluding the corner nodes
                u(iS,iP,eID) = temp(iP,k) !
            ENDDO ! iP, loop over the edge

         ELSE ! south or north sides

            iP = mesh % sideMap(sID)
            DO iS = 1, nS-1 ! Loop over the edge, excluding the corner nodes
               u(iS,iP,eID) = temp(iS,k) ! 
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
   CLASS( CGsemElliptic_2D ), INTENT(in)   :: myCGSEM 
   INTEGER, INTENT(in)                  :: iEl
   TYPE( MAPPEDGEOM_2D ), INTENT(in)    :: geometry
   REAL(prec), INTENT(in)               :: u(0:myCGSEM % cgStorage % nS, 0:myCGSEM % cgStorage % nP)
   REAL(prec), INTENT(out)              :: dudx(0:myCGSEM % cgStorage % nS, 0:myCGSEM % cgStorage % nP)
   REAL(prec), INTENT(out)              :: dudy(0:myCGSEM % cgStorage % nS, 0:myCGSEM % cgStorage % nP) 
   ! LOCAL
   INTEGER :: iS, iP, nS, nP
   REAL(prec) :: duds(0:myCGSEM % cgStorage % nS, 0:myCGSEM % cgStorage % nP)
   REAL(prec) :: dudp(0:myCGSEM % cgStorage % nS, 0:myCGSEM % cgStorage % nP)
   REAL(prec) :: dMatS(0:myCGSEM % cgStorage % nS, 0:myCGSEM % cgStorage % nS)
   REAL(prec) :: dMatP(0:myCGSEM % cgStorage % nP, 0:myCGSEM % cgStorage % nP)
   REAL(prec) :: sTemp(0:myCGSEM % cgStorage % nS)
   REAL(prec) :: pTemp(0:myCGSEM % cgStorage % nP)
   REAL(prec) :: J, dxds, dxdp, dyds, dydp

      nS = myCGSEM % cgStorage % nS
      nP = myCGSEM % cgStorage % nP
      CALL myCGSEM % cgStorage % GetDerivativeMatrix( dMatS, dMatP ) 

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

            dudx(iS,iP) = ( dydp*duds - dyds*dudp )/J
            dudy(iS,iP) = ( dxds*dudp - dxdp*duds )/J        

         ENDDO
      ENDDO

 END SUBROUTINE CalculateGradient_CGsemElliptic_2D
!
!
!
 SUBROUTINE FluxDivergence_CGsemElliptic_2D( myCGSEM, iEl, u, Lu, nS, nP )
 ! S/R FluxDivergence
 !
 !
 ! ============================================================================= !
 ! DECLARATIONS
   IMPLICIT NONE
   integer, INTENT(in)                  :: nS, nP
   TYPE( MAPPEDGEOM_2D ), INTENT(in)    :: geometry
   TYPE( NODAL_STORAGE_2D ), INTENT(in) :: cgStorage
   TYPE( FLUX_COEFFICIENTS), INTENT(in) :: coeffs
   REAL(prec), INTENT(in)               :: u(0:myCGSEM % cgStorage % nS, 0:myCGSEM % cgStorage % nP)
   REAL(prec), INTENT(out)              :: Lu(0:myCGSEM % cgStorage % nS, 0:myCGSEM % cgStorage % nP)
   ! LOCAL
   INTEGER    :: iS, iP
   REAL(prec) :: F1(0:myCGSEM % cgStorage % nS, 0:myCGSEM % cgStorage % nP)
   REAL(prec) :: F2(0:myCGSEM % cgStorage % nS, 0:myCGSEM % cgStorage % nP)
   REAL(prec) :: dF1ds(0:myCGSEM % cgStorage % nS, 0:myCGSEM % cgStorage % nP)
   REAL(prec) :: dF2dp(0:myCGSEM % cgStorage % nS, 0:myCGSEM % cgStorage % nP)
   REAL(prec) :: dudx(0:myCGSEM % cgStorage % nS, 0:myCGSEM % cgStorage % nP)
   REAL(prec) :: dudy(0:myCGSEM % cgStorage % nS, 0:myCGSEM % cgStorage % nP)
   REAL(prec) :: dMatS(0:myCGSEM % cgStorage % nS, 0:myCGSEM % cgStorage % nS)
   REAL(prec) :: dMatSt(0:myCGSEM % cgStorage % nS, 0:myCGSEM % cgStorage % nS)
   REAL(prec) :: dMatP(0:myCGSEM % cgStorage % nP, 0:myCGSEM % cgStorage % nP)
   REAL(prec) :: dMatPt(0:myCGSEM % cgStorage % nP, 0:myCGSEM % cgStorage % nP)
   REAL(prec) :: ws(0:myCGSEM % cgStorage % nS)
   REAL(prec) :: wp(0:myCGSEM % cgStorage % nP)
   REAL(prec) :: sTemp(0:myCGSEM % cgStorage % nS)
   REAL(prec) :: pTemp(0:myCGSEM % cgStorage % nP)
   REAL(prec) :: J, dxds, dxdp, dyds, dydp

      nS = myCGSEM % cgStorage % nS
      nP = myCGSEM % cgStorage % nP

      CALL myCGSEM % cgStorage % GetQuadratureWeights( ws, wp )
      CALL myCGSEM % cgStorage % GetDerivativeMatrix( dMatS, dMatP )
      dMatSt = TRANSPOSE( dMatS )
      dMatPt = TRANSPOSE( dMatP )

      CALL myCGSEM % CalculateGradient( iEl, u, dudx, dudy ) ! calculate the gradient in physical coordinates

      ! Now apply the metric terms
      DO iP = 0, nP ! Loop over the second computational direction
         DO iS = 0, nS ! Loop over the first computational direction

            CALL myCGSEM % mesh % GetJacobianAtNode( iEl, J, iS, iP )
            CALL myCGSEM % mesh % GetCovariantMetricsAtNode( iEl, dxds, dxdp, dyds, dydp, iS, iP )

            ! Contravariant flux calculation
            F1(iS,iP) = ( dudx(iS,iP)*dydp - dudy(iS,iP)*dxdp )*ws(iS)
            F2(iS,iP) = ( dudy(iS,iP)*dxds - dudx(iS,iP)*dyds )*wp(iP)        

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
 !
 ! 
 ! ==========================================================!
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS( CGSEM2D), INTENT(inout) :: myCGSEM
   real(prec), INTENT(inout)       :: u(0:myCGSEM % cgStorage % nS, &
                                        0:myCGSEM % cgStorage % nP, &
                                        1:myCGSEM % mesh % nElems )
   real(prec), INTENT(out)         :: Au(0:myCGSEM % cgStorage % nS, &
                                         0:myCGSEM % cgStorage % nP, &
                                         1:myCGSEM % mesh % nElems  )
   ! Local
   integer :: iEl, nS, nP
   real(prec) :: temp(0:myCGSEM % cgStorage % nS, &
                       0:myCGSEM % cgStorage % nP )
   real(prec) :: Atemp(0:myCGSEM % cgStorage % nS, &
                       0:myCGSEM % cgStorage % nP )


     nS = myCGSEM % cgStorage % nS
     nP = myCGSEM % cgStorage % nP
     
     CALL myCGSEM % UNMASK( u ) ! unmask the solution

!$OMP PARALLEL PRIVATE( temp, Atemp)
!$OMP DO 
     do iEl = 1, myCGSEM  % mesh % nElems

        temp(0:nS,0:nP) = u(0:nS,0:nP,iEl)


        CALL FLUX_DIVERGENCE( myCGSEM % mesh % elements(iEl) % geometry, &
                              myCGSEM % cgStorage, &
                              myCGSEM % fluxCoeffs(iEl), &
                              temp, &
                              Atemp, nS, nP )

        Au(0:nS,0:nP,iEl) = Atemp(0:nS,0:nP)

     enddo
!$OMP END DO

!$OMP  FLUSH ( Au )

!$OMP END PARALLEL

     ! Add the contributions from the shared edges and corners
     CALL myCGSEM % GLOBAL_SUM( Au )

     ! Mask the edges in preparation for residual calculation

     CALL myCGSEM % MASK( Au )

     ! Re-mask the solution

     CALL myCGSEM % MASK( u )


 END SUBROUTINE MATRIX_ACTION
!
!
!
 SUBROUTINE RESIDUAL( myCGSEM, dirF, r )
 !
 ! 
 ! ==========================================================!
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS( CGSEM2D ), INTENT(inout) :: myCGSEM
   real(prec), external            :: dirF
   real(prec), INTENT(out)         :: r(0:myCGSEM % cgStorage % nS, &
                                        0:myCGSEM % cgStorage % nP, &
                                        1:myCGSEM % mesh % nElems )
   ! LOCAL
   integer    :: iS, iP, nS, nP, nEl, iEl
   integer    :: iEdge, e1, s1, nB, nEdges, wallcount
   real(prec) :: w1, w2, J, length, rhs
   real(prec) :: nhat(1:2), x, y
   real(prec) :: temp(0:myCGSEM % cgStorage % nS, &
                      0:myCGSEM % cgStorage % nP, &
                      1:myCGSEM % mesh % nElems )
   real(prec) :: locR(0:myCGSEM % cgStorage % nS, &
                      0:myCGSEM % cgStorage % nP)
   real(prec) :: bfac(1:4)

      bfac = (/ -ONE, ONE, ONE, -ONE /)

      nS = myCGSEM % cgStorage % nS 
      nP = myCGSEM % cgStorage % nP
      nEl = myCGSEM % mesh % nElems 
      nEdges = myCGSEM % mesh % nEdges

      CALL myCGSEM % UNMASK( myCGSEM % solution )

!$OMP PARALLEL PRIVATE( w1, w2, J, locR )
!$OMP DO      
      do iEl = 1, nEl


        CALL FLUX_DIVERGENCE( myCGSEM % mesh % elements(iEl) % geometry, &
                               myCGSEM % cgStorage, &
                               myCGSEM % fluxCoeffs(iEl), & 
                               myCGSEM % solution(:,:,iEl),  &
                               locR, nS, nP )


         do iP = 0, nP
            do iS = 0, nS

               w1 = myCGSEM % cgStorage % qWeightX(iS)
               w2 = myCGSEM % cgStorage % qWeightY(iP)

               J = myCGSEM % mesh % elements(iEl) % geometry % J(iS,iP)
 
               r(iS,iP,iEl) = myCGSEM % source(iS,iP,iEl)*J*w1*w2 - locR(iS,iP)! gather source terms
           

           enddo
        enddo

     enddo
!$OMP END DO

!$OMP  FLUSH ( r )

!$OMP END PARALLEL

  ! Add in the boundary flux
       do iEdge = 1, nEdges
 
           if( myCGSEM % mesh % edges(iEdge) % elementIDs(2) == NEUMANN_WALL )then 


              ! Gather the local addressing information
              e1 = myCGSEM % mesh % edges(iEdge) % elementIDs(1)
              s1 = myCGSEM % mesh % edges(iEdge) % elementSides(1)
  
              nB = myCGSEM % mesh % sideMap(s1) ! get the local quadrature node ID for this boundary

              if( s1 == 1 .OR. s1 == 3 )then ! south or north boundary (respectively)

                  do iS = 0, nS

                     ! Get nHat
                     CALL myCGSEM % mesh % elements(e1) % geometry % GET_NHAT_ATPOINT( s1, iS, nhat, length ) 
 
                     x = myCGSEM % mesh % elements(e1) % geometry % x(iS,nB)
                     y = myCGSEM % mesh % elements(e1) % geometry % y(iS,nB)
                     w1 = myCGSEM % cgStorage % qWeightX(iS)


                     rhs = myCGSEM % boundaryFlux(e1,s1,iS)*length*bfac(s1) 
                     ! From benchmark05 -- toporesponse
                     ! rhs = -sin(pi*x/(10.0_prec**6))*length*bfac(s1)

                     r(iS, nB, e1) = r(iS, nB, e1) - rhs*w1

                  enddo


              elseif( s1 == 2 .OR. s1 == 4 )then ! east or west boundary (respectively)

                 do iP = 0, nP

                     ! Get nHat
                     CALL myCGSEM % mesh % elements(e1) % geometry % GET_NHAT_ATPOINT( s1, iP, nhat, length ) 
                     x = myCGSEM % mesh % elements(e1) % geometry % x(nB,iP)
                     y = myCGSEM % mesh % elements(e1) % geometry % y(nB,iP)
                     w2 = myCGSEM % cgStorage % qWeightY(iP)

                     rhs = myCGSEM % boundaryFlux(e1,s1,iP)*length*bfac(s1) 
                     ! From benchmark05 -- toporesponse
                     !rhs = -sin(pi*x/(10.0_prec**6))*length*bfac(s1) 

                     r(nB, iP, e1) = r(nB, iP, e1) - rhs*w2

                  enddo

              endif
 
           endif

        enddo ! iEdge, loop over the edges

 
     ! Add the contributions from the shared edges and corners
     CALL myCGSEM % GLOBAL_SUM( r )

     ! Mask the edges in preparation for residual calculation

     CALL myCGSEM % MASK( r )

     ! Re-mask the solution
     CALL myCGSEM % MASK( myCGSEM % solution )

     ! re-set the dirichlet boundary conditions
     CALL myCGSEM % SET_DIRICHLET_BOUNDARY_CONDITION( dirF )


 END SUBROUTINE RESIDUAL
!
!
!
  SUBROUTINE VECTOR_PRODUCT( myCGSEM, u, v, uDotv )
 !
 ! Computes the dot product of two "vectors" whose data are
 ! stored in the "(iS,iP)" format
 !
 ! ==========================================================!
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS( CGSEM2D ), INTENT(inout)  :: myCGSEM
   real(prec), INTENT(inout)        :: u(0:myCGSEM % cgStorage % nS, &
                                         0:myCGSEM % cgStorage % nP, &
                                         1:myCGSEM % mesh % nElems )
   real(prec), INTENT(inout)        :: v(0:myCGSEM % cgStorage % nS, &
                                         0:myCGSEM % cgStorage % nP, &
                                         1:myCGSEM % mesh % nElems )
   real(prec), INTENT(out)          :: uDotv
   ! Local
   integer :: iS, iP, nS, nP, iEl
   
      nS = myCGSEM % cgStorage % nS 
      nP = myCGSEM % cgStorage % nP

     ! mask the arrays
     CALL myCGSEM % MASK( u )
     CALL myCGSEM % MASK( v )
 
     uDotv = ZERO

     ! Compute the dot product
     do iEl = 1, myCGSEM % mesh % nElems
        do iP = 0, nP

           uDotV = uDotV + DOT_PRODUCT( u(0:nS,iP,iEl), v(0:nS,iP,iEl) )

        enddo
     enddo ! iEl, loop over the elements

     !CALL myCGSEM % UNMASK( u )
     !CALL myCGSEM % UNMASK( v )

 END SUBROUTINE VECTOR_PRODUCT
!
!
!
 
!
!
!================================================================================!
!=======================       FILE I/O      ===========================!
!================================================================================!
!
!
!
 SUBROUTINE WRITE_TECPLOT_CGSEM2D( myCGSEM ) ! nPlot, nOld, plotInterp, Tmat, filename )
 ! WRITE_TECPLOT_CGSEM2D
 !  Description :
 ! 
 !
 !    Subroutine dependencies :
 !    (NONE)
 !    
 !
 !  Input :
 !    type(LAG_INTERP2D) :: myPoly
 !
 !    real(prec) :: fAtNodes(:,:) - the interpolation nodes function values.
 !
 !    character(*) :: filename - name of output file
 !
 !
 ! 
 !------------------------------------------------------------------------
  IMPLICIT NONE
  CLASS( CGSEM2D ), INTENT(in)           :: myCGSEM
!  integer, INTENT(in)                   :: nOld, nPlot
!  real(prec), INTENT(in)                :: Tmat(0:nPlot, 0:nOld)
!  TYPE( LAG_INTERP2D ), INTENT(in)      :: plotInterp
!  character(*), INTENT(in)              :: filename
  !LOCAL
  integer :: iX, iY, iZ, iEl, nS
  character(len=5) :: zoneID

    nS = myCGSEM % cgStorage % nS

    open( unit=2, file= 'CheckState.tec', form='formatted',status='replace')

    write(2,*) 'VARIABLES = "X", "Y", "Sol"'
    
    do iEL = 1, myCGSEM % mesh % nElems

    
      ! CALL COARSE_TO_FINE_2D(myCGSEM % cgStorage % interp,&
      !                                myCGSEM % mesh % elements(iEl) %  geometry % x,&
      !                                plotInterp, x, Tmat, Tmat)


      ! CALL COARSE_TO_FINE_2D(myCGSEM % cgStorage % interp,&
      !                                myCGSEM % mesh % elements(iEl) %  geometry % y,&
      !                                plotInterp, y, Tmat, Tmat)

      ! CALL COARSE_TO_FINE_2D(myCGSEM % cgStorage % interp,&
      !                                myCGSEM % solution(:,:,iEl),&
      !                                plotInterp, u, Tmat, Tmat)
 
      !  CALL COARSE_TO_FINE_2D(myCGSEM % cgStorage % interp,&
      !                                myCGSEM % mesh % elements(iEl) % geometry % J,&
      !                                plotInterp, J, Tmat, Tmat)

        write(zoneID,'(I5.5)') iEl


        write(2,*) 'ZONE T="el'//trim(zoneID)//'", I=',nS+1,', J=', nS+1,',F=POINT'

        do iY = 0, nS
           do iX = 0, nS

              write (2,*)  myCGSEM % mesh % elements(iEl) % geometry % x( iX, iY ), &
                           myCGSEM % mesh % elements(iEl) % geometry % y( iX, iY ), &
                           myCGSEM % solution(iX,iY,iEl) 

           enddo
        enddo

    enddo

    close(unit=2)

 END SUBROUTINE WRITE_TECPLOT_CGSEM2D

END MODULE CGsemElliptic_2D_Class
