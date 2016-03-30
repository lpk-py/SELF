! QuadMeshClass.f90 ( new with v2.1 - 14 Dec. 2015)
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
! o  (ver 1.0) February 2014
! o  (ver 2.1) December 2015
!
! ========================================= Logs ================================================= !
! yyyy-mm-dd  Joe  <joe@clay>
!
!   
! 
! //////////////////////////////////////////////////////////////////////////////////////////////// !


MODULE QuadMeshClass

! src/common/
USE ModelPrecision
USE ConstantsDictionary
USE ModelFlags
USE LinkedListClass
USE HashTableClass
! src/interp/
USE Chebyshev
USE Lagrange_2D_Class
! src/geom/
USE MappedGeometryClass_2D
USE QuadElementClass
USE EdgeClass
USE NodeClass_2D



IMPLICIT NONE


! All properties are left as public in version 2.1 in order to allow the programmer the ability to
! have direct access to the attributes of the mesh class. This way, it can later be determined if
! additional accessor routines would be necessary.
 
    TYPE QuadMesh 
       INTEGER                          :: nElems, nNodes, nEdges
       TYPE( QuadElement ), ALLOCATABLE :: elements(:)
       TYPE( Node_2D ), ALLOCATABLE     :: nodes(:)  
       TYPE( Edge ), ALLOCATABLE        :: edges(:)
       INTEGER                          :: cornerMap(1:2,1:4) 
       INTEGER                          :: sideMap(1:4) 
       INTEGER                          :: edgeMap(1:2,1:4) 

       CONTAINS

       PROCEDURE :: Build => Build_QuadMesh
       PROCEDURE :: Trash => Trash_QuadMesh
       
       PROCEDURE :: SetNumberOfElements => SetNumberOfElements_QuadMesh
       PROCEDURE :: GetNumberOfElements => GetNumberOfElements_QuadMesh
       PROCEDURE :: SetNumberOfNodes => SetNumberOfNodes_QuadMesh
       PROCEDURE :: GetNumberOfNodes => GetNumberOfNodes_QuadMesh
       PROCEDURE :: SetNumberOfEdge => SetNumberOfEdges_QuadMesh
       PROCEDURE :: GetNumberOfEdges => GetNumberOfEdges_QuadMesh
       
       !QuadElement Wrapper Routines
       PROCEDURE :: SetElementNodeIDs => SetElementNodeIDs_QuadMesh
       PROCEDURE :: GetElementNodeIDs => GetElementNodeIDs_QuadMesh
       PROCEDURE :: SetElementNodeID => SetElementNodeID_QuadMesh
       PROCEDURE :: GetElementNodeID => GetElementNodeID_QuadMesh
       PROCEDURE :: SetElementID => SetElementID_QuadMesh
       PROCEDURE :: GetElementID => GetElementID_QuadMesh
       PROCEDURE :: SetElementNeighbor => SetElementNeighbor_QuadMesh
       PROCEDURE :: GetElementNeighbor => GetElementNeighbor_QuadMesh
       PROCEDURE :: SetElementSouthernNeighbor => SetElementSouthernNeighbor_QuadMesh
       PROCEDURE :: GetElementSouthernNeighbor => GetElementSouthernNeighbor_QuadMesh
       PROCEDURE :: SetElementNorthernNeighbor => SetElementNorthernNeighbor_QuadMesh
       PROCEDURE :: GetElementNorthernNeighbor => GetElementNorthernNeighbor_QuadMesh
       PROCEDURE :: SetElementEasternNeighbor  => SetElementEasternNeighbor_QuadMesh
       PROCEDURE :: GetElementEasternNeighbor  => GetElementEasternNeighbor_QuadMesh
       PROCEDURE :: SetElementWesternNeighbor  => SetElementWesternNeighbor_QuadMesh
       PROCEDURE :: GetElementWesternNeighbor  => GetElementWesternNeighbor_QuadMesh
      
       ! MappedGeometry_2D Wrapper Routines
       PROCEDURE :: SetNumberOfInternalNodes => SetNumberOfInternalNodes_QuadMesh
       PROCEDURE :: GetNumberOfInternalNodes => GetNumberOfInternalNodes_QuadMesh
       PROCEDURE :: SetPositions => SetPositions_QuadMesh
       PROCEDURE :: GetPositions => GetPositions_QuadMesh
       PROCEDURE :: SetPositionAtNode => SetPositionAtNode_QuadMesh
       PROCEDURE :: GetPositionAtNode => GetPositionAtNode_QuadMesh
       PROCEDURE :: SetJacobian => SetJacobian_QuadMesh
       PROCEDURE :: GetJacobian => GetJacobian_QuadMesh
       PROCEDURE :: SetJacobianAtNode => SetJacobianAtNode_QuadMesh
       PROCEDURE :: GetJacobianAtNode => GetJacobianAtNode_QuadMesh
       PROCEDURE :: SetCovariantMetrics => SetCovariantMetrics_QuadMesh
       PROCEDURE :: GetCovariantMetrics => GetCovariantMetrics_QuadMesh
       PROCEDURE :: SetCovariantMetricsAtNode => SetCovariantMetricsAtNode_QuadMesh
       PROCEDURE :: GetCovariantMetricsAtNode => GetCovariantMetricsAtNode_QuadMesh
       PROCEDURE :: SetBoundaryLocation => SetBoundaryLocation_QuadMesh
       PROCEDURE :: GetBoundaryLocation => GetBoundaryLocation_QuadMesh
       PROCEDURE :: SetBoundaryLocationAtNode => SetBoundaryLocationAtNode_QuadMesh
       PROCEDURE :: GetBoundaryLocationAtNode => GetBoundaryLocationAtNode_QuadMesh
       PROCEDURE :: SetBoundaryNormalAtNode => SetBoundaryNormalAtNode_QuadMesh
       PROCEDURE :: GetBoundaryNormalAtNode => GetBoundaryNormalAtNode_QuadMesh
       
       ! NodeClass_2D Wrapper Routines
       PROCEDURE :: SetNodeData => SetNodeData_QuadMesh
       PROCEDURE :: GetNodeData => GetNodeData_QuadMesh
       PROCEDURE :: SetNodeKey => SetNodeKey_QuadMesh
       PROCEDURE :: GetNodeKey => GetNodeKey_QuadMesh
       PROCEDURE :: SetNodeType => SetNodeType_QuadMesh
       PROCEDURE :: GetNodeType => GetNodeType_QuadMesh
       PROCEDURE :: SetNodePosition => SetNodePosition_QuadMesh
       PROCEDURE :: GetNodePosition => GetNodePosition_QuadMesh

       ! EdgeClass Wrapper Routines 
       PROCEDURE :: SetEdgeKey => SetEdgeKey_QuadMesh
       PROCEDURE :: GetEdgeKey => GetEdgeKey_QuadMesh
       PROCEDURE :: SetEdgeData => SetEdgeData_QuadMesh
       PROCEDURE :: GetEdgeData => GetEdgeData_QuadMesh
       PROCEDURE :: SetEdgeNodeIDs => SetEdgeNodeIDs_QuadMesh
       PROCEDURE :: GetEdgeNodeIDs => GetEdgeNodeIDs_QuadMesh
       PROCEDURE :: SetEdgeElementIDs => SetEdgeElementIDs_QuadMesh
       PROCEDURE :: GetEdgeElementIDs => GetEdgeElementIDs_QuadMesh
       PROCEDURE :: SetEdgePrimaryElementID => SetEdgePrimaryElementID_QuadMesh
       PROCEDURE :: GetEdgePrimaryElementID => GetEdgePrimaryElementID_QuadMesh
       PROCEDURE :: SetEdgeSecondaryElementID => SetEdgeSecondaryElementID_QuadMesh
       PROCEDURE :: GetEdgeSecondaryElementID => GetEdgeSecondaryElementID_QuadMesh
       PROCEDURE :: SetEdgeElementSides => SetEdgeElementSides_QuadMesh
       PROCEDURE :: GetEdgeElementSides => GetEdgeElementSides_QuadMesh
       PROCEDURE :: SetEdgePrimaryElementSide => SetEdgePrimaryElementSide_QuadMesh
       PROCEDURE :: GetEdgePrimaryElementSide => GetEdgePrimaryElementSide_QuadMesh
       PROCEDURE :: SetEdgeSecondaryElementSide => SetEdgeSecondaryElementSide_QuadMesh
       PROCEDURE :: GetEdgeSecondaryElementSide => GetEdgeSecondaryElementSide_QuadMesh
       PROCEDURE :: SetEdgeStart => SetEdgeStart_QuadMesh
       PROCEDURE :: GetEdgeStart => GetEdgeStart_QuadMesh
       PROCEDURE :: SetEdgeIncrement => SetEdgeIncrement_QuadMesh
       PROCEDURE :: GetEdgeIncrement => GetEdgeIncrement_QuadMesh
       
       PROCEDURE :: ConstructEdges => ConstructEdges_QuadMesh
       PROCEDURE :: GetNodeToElementConnectivity => GetNodeToElementConnectivity_QuadMesh
       PROCEDURE :: ScaleTheMesh => ScaleTheMesh_QuadMesh
       PROCEDURE :: LoadDefaultMesh => LoadDefaultMesh_QuadMesh
       
       PROCEDURE :: ReadSpecMeshFile => ReadSpecMeshFile_QuadMesh
       PROCEDURE :: WriteTecplot => WriteTecplot_Quadmesh
       
    END TYPE QuadMesh

 INTEGER, PRIVATE, PARAMETER    :: nXElemDefault = 5
 INTEGER, PRIVATE, PARAMETER    :: nYElemDefault = 5
 INTEGER, PRIVATE, PARAMETER    :: BoundaryFlagDefault = NO_NORMAL_FLOW
 REAL(prec), PRIVATE, PARAMETER :: dXDefault = ONE/(nXElemDefault)
 REAL(prec), PRIVATE, PARAMETER :: dYDefault = ONE/(nYElemDefault)


 CONTAINS
!
!
!==================================================================================================!
!------------------------------- Manual Constructors/Destructors ----------------------------------!
!==================================================================================================!
!
!
 SUBROUTINE Build_QuadMesh( myQuadMesh, nNodes, nElems, nEdges, nS )
 ! S/R Build
 ! 
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS(QuadMesh), INTENT(out) :: myQuadMesh
   INTEGER, INTENT(in)          :: nNodes, nElems, nEdges, nS
   !LOCAL
   INTEGER :: iNode

      myQuadmesh % sideMap(1:4) = (/ 0, nS, nS, 0 /)
      myQuadmesh % cornerMap(1, 1:4) = (/ 0, nS, nS, 0 /)
      myQuadmesh % cornerMap(2, 1:4) = (/ 0, 0, nS, nS /)
      myQuadmesh % edgeMap(1, 1:4) = (/ 1, 2, 4, 1 /)
      myQuadmesh % edgeMap(2, 1:4) = (/ 2, 3, 3, 4 /)

      myQuadmesh % nNodes = nNodes
      myQuadmesh % nElems = nElems

      ALLOCATE( myQuadmesh % elements(1:nElems) )
      ALLOCATE( myQuadmesh % nodes(1:nNodes) )
      ALLOCATE( myQuadmesh % edges(1:nEdges) )
     

      ! Default nodes to origin
      DO iNode = 1, myQuadmesh % nNodes ! loop over the number of nodes
         CALL myQuadmesh % nodes(iNode) % Build( ZERO, ZERO )
      ENDDO ! iNode, loop over the number of nodes
     

 END SUBROUTINE Build_QuadMesh
!
!
!
 SUBROUTINE Trash_QuadMesh( myQuadMesh )
 ! S/R Trash
 ! 
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS(QuadMesh), INTENT(inout) :: myQuadMesh
  ! LOCAL
   INTEGER :: iNode, iEl

      DO iEl = 1, myQuadMesh % nElems
         CALL myQuadMesh % elements(iEl) % TRASH( )
      ENDDO

      DO iNode = 1, myQuadMesh % nNodes
         CALL myQuadMesh % nodes(iNode) % TRASH( )
      ENDDO

      DEALLOCATE( myQuadMesh % nodes, myQuadMesh % elements, myQuadMesh % edges )
      

 END SUBROUTINE Trash_QuadMesh
!
!
!==================================================================================================!
!--------------------------------------- Accessors ------------------------------------------------!
!==================================================================================================!
!
! 
 SUBROUTINE SetNumberOfElements_QuadMesh( myQuadMesh, nElems )
 ! S/R SetNumberOfElements
 ! 
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS(QuadMesh), INTENT(inout) :: myQuadMesh
   INTEGER, INTENT(in)            :: nElems
   
      myQuadMesh % nElems = nElems
      
 END SUBROUTINE SetNumberOfElements_QuadMesh
!
!
!
 SUBROUTINE GetNumberOfElements_QuadMesh( myQuadMesh, nElems )
 ! S/R GetNumberOfElements
 ! 
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS(QuadMesh), INTENT(in) :: myQuadMesh
   INTEGER, INTENT(out)        :: nElems
   
      nElems = myQuadMesh % nElems
      
 END SUBROUTINE GetNumberOfElements_QuadMesh
!
!
!
 SUBROUTINE SetNumberOfNodes_QuadMesh( myQuadMesh, nNodes )
 ! S/R SetNumberOfNodes
 ! 
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS(QuadMesh), INTENT(inout) :: myQuadMesh
   INTEGER, INTENT(in)            :: nNodes
   
      myQuadMesh % nNodes = nNodes
      
 END SUBROUTINE SetNumberOfNodes_QuadMesh
!
!
!
 SUBROUTINE GetNumberOfNodes_QuadMesh( myQuadMesh, nNodes )
 ! S/R GetNumberOfNodes
 ! 
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS(QuadMesh), INTENT(in) :: myQuadMesh
   INTEGER, INTENT(out)        :: nNodes
   
      nNodes = myQuadMesh % nNodes
      
 END SUBROUTINE GetNumberOfNodes_QuadMesh
!
!
!
 SUBROUTINE SetNumberOfEdges_QuadMesh( myQuadMesh, nEdges )
 ! S/R SetNumberOfEdges
 ! 
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS(QuadMesh), INTENT(inout) :: myQuadMesh
   INTEGER, INTENT(in)            :: nEdges
   
      myQuadMesh % nEdges = nEdges
      
 END SUBROUTINE SetNumberOfEdges_QuadMesh
!
!
!
 SUBROUTINE GetNumberOfEdges_QuadMesh( myQuadMesh, nEdges )
 ! S/R GetNumberOfEdges
 ! 
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS(QuadMesh), INTENT(in) :: myQuadMesh
   INTEGER, INTENT(out)        :: nEdges
   
      nEdges = myQuadMesh % nEdges
      
 END SUBROUTINE GetNumberOfEdges_QuadMesh
!
! ---------------------------- QuadElement Wrapper Routines -------------------------------------- !
!
SUBROUTINE SetElementNodeIDs_Quadmesh( myQuadMesh, iEl, nodeIDs )
 ! S/R SetElementNodeIDs
 !
 ! =============================================================================================== !
   IMPLICIT NONE
   CLASS( Quadmesh ), INTENT(inout) :: myQuadMesh
   INTEGER, INTENT(in)              :: iEl
   INTEGER, INTENT(in)              :: nodeIDs(1:4)
   
      CALL myQuadMesh % elements(iEl)  % SetNodeIDs( nodeIDs )

 END SUBROUTINE SetElementNodeIDs_Quadmesh
!
!
!
 SUBROUTINE GetElementNodeIDs_Quadmesh( myQuadMesh, iEl, nodeIDs )
 ! S/R GetElementNodeIDs
 !
 ! =============================================================================================== !
   IMPLICIT NONE
   CLASS( Quadmesh ), INTENT(in) :: myQuadMesh
   INTEGER, INTENT(in)           :: iEl
   INTEGER, INTENT(out)          :: nodeIDs(1:4)
   
      CALL myQuadMesh % elements(iEl)  % GetNodeIDs( nodeIDs )

 END SUBROUTINE GetElementNodeIDs_Quadmesh
!
!
!
 SUBROUTINE SetElementNodeID_Quadmesh( myQuadMesh, iEl, localID, nodeID )
 ! S/R SetElementNodeID
 !
 ! =============================================================================================== !
   IMPLICIT NONE
   CLASS( Quadmesh ), INTENT(inout) :: myQuadMesh
   INTEGER, INTENT(in)              :: iEl
   INTEGER, INTENT(in)              :: localID, nodeID
   
      CALL myQuadMesh % elements(iEl)  % SetNodeID( localID, nodeID )

 END SUBROUTINE SetElementNodeID_Quadmesh
!
!
!
 SUBROUTINE GetElementNodeID_Quadmesh( myQuadMesh, iEl, localID, nodeID )
 ! S/R GetElementNodeID
 !
 ! =============================================================================================== !
   IMPLICIT NONE
   CLASS( Quadmesh ), INTENT(in) :: myQuadMesh
   INTEGER, INTENT(in)           :: iEl
   INTEGER, INTENT(in)           :: localID
   INTEGER, INTENT(out)          :: nodeID
   
      CALL myQuadMesh % elements(iEl) % GetNodeID( localID, nodeID )

 END SUBROUTINE GetElementNodeID_Quadmesh
!
!
!
 SUBROUTINE SetElementID_Quadmesh( myQuadMesh, iEl, eID )
 ! S/R SetElementID
 !
 ! =============================================================================================== !
   IMPLICIT NONE
   CLASS( Quadmesh ), INTENT(inout) :: myQuadMesh
   INTEGER, INTENT(in)              :: iEl
   INTEGER, INTENT(in)              :: eID
   
      CALL myQuadMesh % elements(iEl)  % SetElementID( eID )

 END SUBROUTINE SetElementID_Quadmesh
!
!
!
 SUBROUTINE GetElementID_Quadmesh( myQuadMesh, iEl, eID )
 ! S/R GetElementID
 !
 ! =============================================================================================== !
   IMPLICIT NONE
   CLASS( Quadmesh ), INTENT(in) :: myQuadMesh
   INTEGER, INTENT(in)           :: iEl
   INTEGER, INTENT(out)          :: eID
   
      CALL myQuadMesh % elements(iEl)  % GetElementID( eID )

 END SUBROUTINE GetElementID_Quadmesh
!
!
!
 SUBROUTINE SetElementNeighbor_Quadmesh( myQuadMesh, iEl, sID, eID )
  ! S/R SetElementNeighbor
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS( Quadmesh ), INTENT(inout) :: myQuadMesh
   INTEGER, INTENT(in)              :: iEl
   INTEGER, INTENT(in)              :: sID
   INTEGER, INTENT(in)              :: eID
   
      CALL myQuadMesh % elements(iEl) % SetNeighbor( sID, eID )
      
 END SUBROUTINE SetElementNeighbor_Quadmesh
!
!
!
 SUBROUTINE GetElementNeighbor_Quadmesh( myQuadMesh, iEl, sID, eID )
 ! S/R GetElementNeighbor
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS( Quadmesh ), INTENT(in) :: myQuadMesh
   INTEGER, INTENT(in)           :: iEl
   INTEGER, INTENT(in)           :: sID
   INTEGER, INTENT(out)          :: eID
   
      CALL myQuadMesh % elements(iEl) % GetNeighbor( sID, eID )
      
 END SUBROUTINE GetElementNeighbor_Quadmesh
!
!
!
 SUBROUTINE SetElementSouthernNeighbor_Quadmesh( myQuadMesh, iEl, eID )
 ! S/R SetElementSouthernNeighbor
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS( Quadmesh ), INTENT(inout) :: myQuadMesh
   INTEGER, INTENT(in)              :: iEl
   INTEGER, INTENT(in)              :: eID
   
      CALL myQuadMesh % elements(iEl) % SetSouthernNeighbor( eID )
      
 END SUBROUTINE SetElementSouthernNeighbor_Quadmesh
!
!
!
 SUBROUTINE GetElementSouthernNeighbor_Quadmesh( myQuadMesh, iEl, eID )
 ! S/R GetElementSouthernNeighbor
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS( Quadmesh ), INTENT(in) :: myQuadMesh
   INTEGER, INTENT(in)           :: iEl
   INTEGER, INTENT(out)          :: eID
   
      CALL myQuadMesh % elements(iEl) % GetSouthernNeighbor( eID )
      
 END SUBROUTINE GetElementSouthernNeighbor_Quadmesh
!
!
!
 SUBROUTINE SetElementNorthernNeighbor_Quadmesh( myQuadMesh, iEl, eID )
 ! S/R SetElementNorthernNeighbor
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS( Quadmesh ), INTENT(inout) :: myQuadMesh
   INTEGER, INTENT(in)              :: iEl
   INTEGER, INTENT(in)              :: eID
   
      CALL myQuadMesh % elements(iEl) % SetNorthernNeighbor( eID )
      
 END SUBROUTINE SetElementNorthernNeighbor_Quadmesh
!
!
!
 SUBROUTINE GetElementNorthernNeighbor_Quadmesh( myQuadMesh, iEl, eID )
 ! S/R GetElementNorthernNeighbor
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS( Quadmesh ), INTENT(in) :: myQuadMesh
   INTEGER, INTENT(in)           :: iEl 
   INTEGER, INTENT(out)          :: eID
   
      CALL myQuadMesh % elements(iEl) % GetNorthernNeighbor( eID )
      
 END SUBROUTINE GetElementNorthernNeighbor_Quadmesh
!
!
!
 SUBROUTINE SetElementEasternNeighbor_Quadmesh( myQuadMesh, iEl, eID )
 ! S/R SetElementEasternNeighbor
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS( Quadmesh ), INTENT(inout) :: myQuadMesh
   INTEGER, INTENT(in)              :: iEl
   INTEGER, INTENT(in)              :: eID
   
      CALL myQuadMesh % elements(iEl) % SetEasternNeighbor( eID )
      
 END SUBROUTINE SetElementEasternNeighbor_Quadmesh
!
!
!
 SUBROUTINE GetElementEasternNeighbor_Quadmesh( myQuadMesh, iEl, eID )
 ! S/R GetElementEasternNeighbor
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS( Quadmesh ), INTENT(in) :: myQuadMesh
   INTEGER, INTENT(in)           :: iEl
   INTEGER, INTENT(out)          :: eID
   
      CALL myQuadMesh % elements(iEl) % GetEasternNeighbor( eID )
      
 END SUBROUTINE GetElementEasternNeighbor_Quadmesh
!
!
!
 SUBROUTINE SetElementWesternNeighbor_Quadmesh( myQuadMesh, iEl, eID )
 ! S/R SetElementWesternNeighbor
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS( Quadmesh ), INTENT(inout) :: myQuadMesh
   INTEGER, INTENT(in)              :: iEl
   INTEGER, INTENT(in)              :: eID
   
      CALL myQuadMesh % elements(iEl) % SetWesternNeighbor( eID )
      
 END SUBROUTINE SetElementWesternNeighbor_Quadmesh
!
!
!
 SUBROUTINE GetElementWesternNeighbor_Quadmesh( myQuadMesh, iEl, eID )
 ! S/R GetElementWesternNeighbor
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS( Quadmesh ), INTENT(in) :: myQuadMesh
   INTEGER, INTENT(in)           :: iEl
   INTEGER, INTENT(out)          :: eID
   
      CALL myQuadMesh % elements(iEl) % GetWesternNeighbor( eID )
      
 END SUBROUTINE GetElementWesternNeighbor_Quadmesh
!
! ---------------------------- MappedGeometryClass_2D Wrappers ----------------------------------- !
!
SUBROUTINE SetNumberOfInternalNodes_Quadmesh( myQuadMesh, iEl, nS, nP )
 ! S/R SetNumberOfInternalNodes
 !
 ! =============================================================================================== !
 ! DECLARATIONS
  IMPLICIT NONE
  CLASS(Quadmesh), INTENT(inout) :: myQuadMesh
  INTEGER, INTENT(in)            :: iEl
  INTEGER, INTENT(in)            :: nS, nP
  
     CALL myQuadMesh % elements(iEl) % SetNumberOfNodes( nS, nP )
     
 END SUBROUTINE SetNumberOfInternalNodes_Quadmesh
!
!
!
 SUBROUTINE GetNumberOfInternalNodes_Quadmesh( myQuadMesh, iEl, nS, nP )
 ! S/R GetNumberOfInternalNodes
 !
 ! =============================================================================================== !
 ! DECLARATIONS
  IMPLICIT NONE
  CLASS(Quadmesh), INTENT(in) :: myQuadMesh
  INTEGER, INTENT(in)         :: iEl
  INTEGER, INTENT(out)        :: nS, nP
  
     CALL myQuadMesh % elements(iEl) % GetNumberOfNodes( nS, nP )
     
 END SUBROUTINE GetNumberOfInternalNodes_Quadmesh
!
!
!
 SUBROUTINE SetPositions_Quadmesh( myQuadMesh, iEl, x, y )
 ! S/R SetPositions
 !
 ! =============================================================================================== !
 ! DECLARATIONS
  IMPLICIT NONE
  CLASS(Quadmesh), INTENT(inout) :: myQuadMesh
  INTEGER, INTENT(in)            :: iEl
  REAL(prec), INTENT(in)         :: x(0:myQuadMesh % elements(iEl) % nS, 0:myQuadMesh % elements(iEl) % nP)
  REAL(prec), INTENT(in)         :: y(0:myQuadMesh % elements(iEl) % nS, 0:myQuadMesh % elements(iEl) % nP)
  
     CALL myQuadMesh % elements(iEl) % SetPositions( x, y )
     
 END SUBROUTINE SetPositions_Quadmesh
!
!
!
 SUBROUTINE GetPositions_Quadmesh( myQuadMesh, iEl, x, y )
 ! S/R GetPositions
 !
 ! =============================================================================================== !
 ! DECLARATIONS
  IMPLICIT NONE
  CLASS(Quadmesh), INTENT(in) :: myQuadMesh
  INTEGER, INTENT(in)         :: iEl
  REAL(prec), INTENT(out)     :: x(0:myQuadMesh % elements(iEl) % nS, 0:myQuadMesh % elements(iEl) % nP)
  REAL(prec), INTENT(out)     :: y(0:myQuadMesh % elements(iEl) % nS, 0:myQuadMesh % elements(iEl) % nP)
  
     CALL myQuadMesh % elements(iEl) % GetPositions( x, y )
     
 END SUBROUTINE GetPositions_Quadmesh
!
!
!
 SUBROUTINE SetPositionAtNode_Quadmesh( myQuadMesh, iEl, x, y, i, j )
 ! S/R SetPositionAtNode
 !
 ! =============================================================================================== !
 ! DECLARATIONS
  IMPLICIT NONE
  CLASS(Quadmesh), INTENT(inout) :: myQuadMesh
  INTEGER, INTENT(in)            :: iEl
  REAL(prec), INTENT(in)         :: x, y
  INTEGER, INTENT(in)            :: i, j
  
     CALL myQuadMesh % elements(iEl) % SetPositionAtNode( x, y, i, j )
     
 END SUBROUTINE SetPositionAtNode_Quadmesh
!
!
!
 SUBROUTINE GetPositionAtNode_Quadmesh( myQuadMesh, iEl, x, y, i, j )
 ! S/R GetPositionAtNode
 !
 ! =============================================================================================== !
 ! DECLARATIONS
  IMPLICIT NONE
  CLASS(Quadmesh), INTENT(in) :: myQuadMesh
  INTEGER, INTENT(in)         :: iEl
  REAL(prec), INTENT(out)     :: x, y
  INTEGER, INTENT(in)         :: i, j
  
     CALL myQuadMesh % elements(iEl) % GetPositionAtNode( x, y, i, j )
     
 END SUBROUTINE GetPositionAtNode_Quadmesh
!
!
!
 SUBROUTINE SetJacobian_Quadmesh( myQuadMesh, iEl, J )
 ! S/R SetJacobian
 !
 ! =============================================================================================== !
 ! DECLARATIONS
  IMPLICIT NONE
  CLASS(Quadmesh), INTENT(inout) :: myQuadMesh
  INTEGER, INTENT(in)            :: iEl
  REAL(prec), INTENT(in)         :: J(0:myQuadMesh % elements(iEl) % nS, 0:myQuadMesh % elements(iEl) % nP)
  
     CALL myQuadMesh % elements(iEl) % SetJacobian( J )
     
 END SUBROUTINE SetJacobian_Quadmesh
!
!
!
 SUBROUTINE GetJacobian_Quadmesh( myQuadMesh, iEl, J )
 ! S/R GetJacobian
 !
 ! =============================================================================================== !
 ! DECLARATIONS
  IMPLICIT NONE
  CLASS(Quadmesh), INTENT(in) :: myQuadMesh
  INTEGER, INTENT(in)         :: iEl
  REAL(prec), INTENT(out)     :: J(0:myQuadMesh % elements(iEl) % nS, 0:myQuadMesh % elements(iEl) % nP)
  
     CALL myQuadMesh % elements(iEl) % GetJacobian( J )
     
 END SUBROUTINE GetJacobian_Quadmesh
!
!
!
 SUBROUTINE SetJacobianAtNode_Quadmesh( myQuadMesh, iEl, J, iS, iP )
 ! S/R SetJacobianAtNode
 !
 ! =============================================================================================== !
 ! DECLARATIONS
  IMPLICIT NONE
  CLASS(Quadmesh), INTENT(inout) :: myQuadMesh
  INTEGER, INTENT(in)            :: iEl
  REAL(prec), INTENT(in)         :: J
  INTEGER, INTENT(in)            :: iS, iP
  
     CALL myQuadMesh % elements(iEl) % SetJacobianAtNode( J, iS, iP )
     
 END SUBROUTINE SetJacobianAtNode_Quadmesh
!
!
!
 SUBROUTINE GetJacobianAtNode_Quadmesh( myQuadMesh, iEl, J, iS, iP )
 ! S/R GetJacobianAtNode
 !
 ! =============================================================================================== !
 ! DECLARATIONS
  IMPLICIT NONE
  CLASS(Quadmesh), INTENT(in) :: myQuadMesh
  INTEGER, INTENT(in)         :: iEl
  REAL(prec), INTENT(out)     :: J
  INTEGER, INTENT(in)         :: iS, iP
  
     CALL myQuadMesh % elements(iEl) % GetJacobianAtNode( J, iS, iP )
     
 END SUBROUTINE GetJacobianAtNode_Quadmesh
!
!
!
 SUBROUTINE SetCovariantMetrics_Quadmesh( myQuadMesh, iEl, dxds, dxdp, dyds, dydp )
 ! S/R SetCovariantMetrics
 !
 ! =============================================================================================== !
 ! DECLARATIONS
  IMPLICIT NONE
  CLASS(Quadmesh), INTENT(inout) :: myQuadMesh
  INTEGER, INTENT(in)            :: iEl
  REAL(prec), INTENT(in)         :: dxds(0:myQuadMesh % elements(iEl) % nS, 0:myQuadMesh % elements(iEl) % nP)
  REAL(prec), INTENT(in)         :: dxdp(0:myQuadMesh % elements(iEl) % nS, 0:myQuadMesh % elements(iEl) % nP)
  REAL(prec), INTENT(in)         :: dyds(0:myQuadMesh % elements(iEl) % nS, 0:myQuadMesh % elements(iEl) % nP)
  REAL(prec), INTENT(in)         :: dydp(0:myQuadMesh % elements(iEl) % nS, 0:myQuadMesh % elements(iEl) % nP)
  
     CALL myQuadMesh % elements(iEl) % SetCovariantMetrics( dxds, dxdp, dyds, dydp )
     
 END SUBROUTINE SetCovariantMetrics_Quadmesh
!
!
!
 SUBROUTINE GetCovariantMetrics_Quadmesh( myQuadMesh, iEl, dxds, dxdp, dyds, dydp )
 ! S/R GetCovariantMetrics
 !
 ! =============================================================================================== !
 ! DECLARATIONS
  IMPLICIT NONE
  CLASS(Quadmesh), INTENT(in) :: myQuadMesh
  INTEGER, INTENT(in)         :: iEl
  REAL(prec), INTENT(out)     :: dxds(0:myQuadMesh % elements(iEl) % nS, 0:myQuadMesh % elements(iEl) % nP)
  REAL(prec), INTENT(out)     :: dxdp(0:myQuadMesh % elements(iEl) % nS, 0:myQuadMesh % elements(iEl) % nP)
  REAL(prec), INTENT(out)     :: dyds(0:myQuadMesh % elements(iEl) % nS, 0:myQuadMesh % elements(iEl) % nP)
  REAL(prec), INTENT(out)     :: dydp(0:myQuadMesh % elements(iEl) % nS, 0:myQuadMesh % elements(iEl) % nP)
  
     CALL myQuadMesh % elements(iEl) % GetCovariantMetrics( dxds, dxdp, dyds, dydp )
     
 END SUBROUTINE GetCovariantMetrics_Quadmesh
!
!
!
 SUBROUTINE SetCovariantMetricsAtNode_Quadmesh( myQuadMesh, iEl, dxds, dxdp, dyds, dydp, iS, iP )
 ! S/R SetCovariantMetricsAtNode
 !
 ! =============================================================================================== !
 ! DECLARATIONS
  IMPLICIT NONE
  CLASS(Quadmesh), INTENT(inout) :: myQuadMesh
  INTEGER, INTENT(in)            :: iEl
  REAL(prec), INTENT(in)         :: dxds, dxdp, dyds, dydp
  INTEGER, INTENT(in)            :: iS, iP
  
     CALL myQuadMesh % elements(iEl) % SetCovariantMetricsAtNode( dxds, dxdp, dyds, dydp, iS, iP )
     
 END SUBROUTINE SetCovariantMetricsAtNode_Quadmesh
!
!
!
 SUBROUTINE GetCovariantMetricsAtNode_Quadmesh( myQuadMesh, iEl, dxds, dxdp, dyds, dydp, iS, iP )
 ! S/R GetCovariantMetricsAtNode
 !
 ! =============================================================================================== !
 ! DECLARATIONS
  IMPLICIT NONE
  CLASS(Quadmesh), INTENT(in) :: myQuadMesh
  INTEGER, INTENT(in)         :: iEl
  REAL(prec), INTENT(out)     :: dxds, dxdp, dyds, dydp
  INTEGER, INTENT(in)         :: iS, iP
  
     CALL myQuadMesh % elements(iEl) % GetCovariantMetricsAtNode( dxds, dxdp, dyds, dydp, iS, iP )
     
 END SUBROUTINE GetCovariantMetricsAtNode_Quadmesh
!
!
!
 SUBROUTINE SetBoundaryLocation_Quadmesh( myQuadMesh, iEl, x, y, iBound )
 ! S/R SetBoundaryLocation
 !
 ! =============================================================================================== !
 ! DECLARATIONS
  IMPLICIT NONE
  CLASS(Quadmesh), INTENT(inout) :: myQuadMesh
  INTEGER, INTENT(in)            :: iEl
  INTEGER, INTENT(in)            :: iBound
  REAL(prec), INTENT(in)         :: x(0:myQuadMesh % elements(iEl) % nMax)
  REAL(prec), INTENT(in)         :: y(0:myQuadMesh % elements(iEl) % nMax)
  
     CALL myQuadMesh % elements(iEl) % SetBoundaryLocation( x, y, iBound )
     
 END SUBROUTINE SetBoundaryLocation_Quadmesh
!
!
!
 SUBROUTINE GetBoundaryLocation_Quadmesh( myQuadMesh, iEl, x, y, iBound )
 ! S/R GetBoundaryLocation
 !
 ! =============================================================================================== !
 ! DECLARATIONS
  IMPLICIT NONE
  CLASS(Quadmesh), INTENT(in) :: myQuadMesh
  INTEGER, INTENT(in)         :: iEl
  INTEGER, INTENT(in)         :: iBound
  REAL(prec), INTENT(out)     :: x(0:myQuadMesh % elements(iEl) % nMax)
  REAL(prec), INTENT(out)     :: y(0:myQuadMesh % elements(iEl) % nMax)
  
     CALL myQuadMesh % elements(iEl) % GetBoundaryLocation( x, y, iBound )
     
 END SUBROUTINE GetBoundaryLocation_Quadmesh
!
!
!
 SUBROUTINE SetBoundaryLocationAtNode_Quadmesh( myQuadMesh, iEl, x, y, i, iBound )
 ! S/R SetBoundaryLocationAtNode
 !
 ! =============================================================================================== !
 ! DECLARATIONS
  IMPLICIT NONE
  CLASS(Quadmesh), INTENT(inout) :: myQuadMesh
  INTEGER, INTENT(in)            :: iEl
  REAL(prec), INTENT(in)         :: x, y
  INTEGER, INTENT(in)            :: i, iBound
  
     CALL myQuadMesh % elements(iEl) % SetBoundaryLocationAtNode( x, y, i, iBound )
     
 END SUBROUTINE SetBoundaryLocationAtNode_Quadmesh
!
!
!
 SUBROUTINE GetBoundaryLocationAtNode_Quadmesh( myQuadMesh, iEl, x, y, i, iBound )
 ! S/R GetBoundaryLocationAtNode
 !
 ! =============================================================================================== !
 ! DECLARATIONS
  IMPLICIT NONE
  CLASS(Quadmesh), INTENT(in) :: myQuadMesh
  INTEGER, INTENT(in)         :: iEl
  REAL(prec), INTENT(out)     :: x, y
  INTEGER, INTENT(in)         :: i, iBound
  
     CALL myQuadMesh % elements(iEl) % GetBoundaryLocationAtNode( x, y, i, iBound )
     
 END SUBROUTINE GetBoundaryLocationAtNode_Quadmesh
!
!
!
 SUBROUTINE SetBoundaryNormalAtNode_Quadmesh( myQuadMesh, iEl, dir, i, iBound )
 ! S/R SetBoundaryNormalAtNode
 !
 ! =============================================================================================== !
 ! DECLARATIONS
  IMPLICIT NONE
  CLASS(Quadmesh), INTENT(inout) :: myQuadMesh
  INTEGER, INTENT(in)            :: iEl
  REAL(prec), INTENT(in)         :: dir(1:2)
  INTEGER, INTENT(in)            :: i, iBound
  
     CALL myQuadMesh % elements(iEl) % SetBoundaryNormalAtNode( dir, i, iBound )
     
 END SUBROUTINE SetBoundaryNormalAtNode_Quadmesh
!
!
!
 SUBROUTINE GetBoundaryNormalAtNode_Quadmesh( myQuadMesh, iEl, dir, length, i, iBound )
 ! S/R GetBoundaryNormalAtNode
 !
 ! =============================================================================================== !
 ! DECLARATIONS
  IMPLICIT NONE
  CLASS(Quadmesh), INTENT(in) :: myQuadMesh
  INTEGER, INTENT(in)         :: iEl
  REAL(prec), INTENT(out)     :: dir(1:2), length
  INTEGER, INTENT(in)         :: i, iBound
  
     CALL myQuadMesh % elements(iEl) % GetBoundaryNormalAtNode( dir, length, i, iBound )
     
 END SUBROUTINE GetBoundaryNormalAtNode_Quadmesh
!
! ------------------------------- Node Wrapper Routines ------------------------------------------ !
!
 SUBROUTINE SetNodeData_QuadMesh( myQuadMesh, iNode, x, y, nodeType )
 ! S/R SetNodeData
 !  
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS( QuadMesh ), INTENT(inout) :: myQuadMesh
   INTEGER, INTENT(in)              :: iNode
   REAL(prec), INTENT(in)           :: x, y
   INTEGER, INTENT(in)              :: nodeType

      CALL myQuadMesh % nodes(iNode) % SetData( x, y, nodeType )

 END SUBROUTINE SetNodeData_QuadMesh
!
!
!
 SUBROUTINE GetNodeData_QuadMesh( myQuadMesh, iNode, x, y, nodeType )
 ! S/R GetNodeData
 !  
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS( QuadMesh ), INTENT(in) :: myQuadMesh
   INTEGER, INTENT(in)           :: iNode
   REAL(prec), INTENT(out)       :: x, y
   INTEGER, INTENT(out)          :: nodeType

      CALL myQuadMesh % nodes(iNode) % GetData( x, y, nodeType )

 END SUBROUTINE GetNodeData_QuadMesh
!
!
!
 SUBROUTINE SetNodeKey_QuadMesh( myQuadMesh, iNode, key )
 ! S/R SetNodeKey
 !  
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS( QuadMesh ), INTENT(inout) :: myQuadMesh
   INTEGER, INTENT(in)              :: iNode
   INTEGER, INTENT(in)              :: key
   
      CALL myQuadMesh % nodes(iNode) % SetKey( key )
   
 END SUBROUTINE SetNodeKey_QuadMesh
!
!
!
 SUBROUTINE GetNodeKey_QuadMesh( myQuadMesh, iNode, key )
 ! S/R GetNodeKey
 !  
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS( QuadMesh ), INTENT(in) :: myQuadMesh
   INTEGER, INTENT(in)           :: iNode
   INTEGER, INTENT(out)          :: key
   
      CALL myQuadMesh % nodes(iNode) % GetKey( key )
   
 END SUBROUTINE GetNodeKey_QuadMesh
!
!
!
 SUBROUTINE SetNodeType_QuadMesh( myQuadMesh, iNode, nodeType )
 ! S/R SetNodeType_QuadMesh
 !  
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS( QuadMesh ), INTENT(inout) :: myQuadMesh
   INTEGER, INTENT(in)              :: iNode
   INTEGER, INTENT(in)              :: nodeType
   
      CALL myQuadMesh % nodes(iNode) % SetType( nodeType )
   
 END SUBROUTINE SetNodeType_QuadMesh
!
!
!
 SUBROUTINE GetNodeType_QuadMesh( myQuadMesh, iNode, nodeType )
 ! S/R GetNodeType_QuadMesh
 !  
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS( QuadMesh ), INTENT(in) :: myQuadMesh
   INTEGER, INTENT(in)           :: iNode
   INTEGER, INTENT(out)          :: nodeType
   
      CALL myQuadMesh % nodes(iNode) % GetType( nodeType )
   
 END SUBROUTINE GetNodeType_QuadMesh
!
!
!
 SUBROUTINE SetNodePosition_QuadMesh( myQuadMesh, iNode, x, y )
 ! S/R SetNodePosition_QuadMesh
 !  
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS( QuadMesh ), INTENT(inout) :: myQuadMesh
   INTEGER, INTENT(in)              :: iNode
   REAL(prec), INTENT(in)           :: x, y
   
      CALL myQuadMesh % nodes(iNode) % SetPosition( x, y )
   
 END SUBROUTINE SetNodePosition_QuadMesh
!
!
!
 SUBROUTINE GetNodePosition_QuadMesh( myQuadMesh, iNode, x, y )
 ! S/R GetNodePosition_QuadMesh
 !  
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS( QuadMesh ), INTENT(in) :: myQuadMesh
   INTEGER, INTENT(in)           :: iNode
   REAL(prec), INTENT(out)       :: x, y
   
      CALL myQuadMesh % nodes(iNode) % GetPosition( x, y )
   
 END SUBROUTINE GetNodePosition_QuadMesh
!
! --------------------------------- Edge Wrapper Routines ---------------------------------------- !
!
 SUBROUTINE SetEdgeData_QuadMesh( myQuadMesh, iEdge, nodeIDs, elementIDs, elementSides, key, start, inc )
 ! S/R SetEdgeData
 !  
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS( QuadMesh ), INTENT(inout) :: myQuadMesh
   INTEGER, INTENT(in)              :: iEdge
   INTEGER, INTENT(in)              :: nodeIDs(1:2), elementIDs(1:2), elementSides(1:2)
   INTEGER, INTENT(in)              :: key, start, inc

      CALL myQuadMesh % edges(iEdge) % SetData( nodeIDs, elementIDs, elementSides, key, start, inc )

 END SUBROUTINE SetEdgeData_QuadMesh
!
!
!
 SUBROUTINE GetEdgeData_QuadMesh( myQuadMesh, iEdge, nodeIDs, elementIDs, elementSides, key, start, inc )
 ! S/R GetEdgeData
 !  
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS( QuadMesh ), INTENT(in) :: myQuadMesh
   INTEGER, INTENT(in)           :: iEdge
   INTEGER, INTENT(out)          :: nodeIDs(1:2), elementIDs(1:2), elementSides(1:2)
   INTEGER, INTENT(out)          :: key, start, inc

      CALL myQuadMesh % edges(iEdge) % GetData( nodeIDs, elementIDs, elementSides, key, start, inc )
      
 END SUBROUTINE GetEdgeData_QuadMesh
!
!
!
SUBROUTINE SetEdgeKey_QuadMesh( myQuadMesh, iEdge, key )
 ! S/R SetEdgeKey
 !  
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS( QuadMesh ), INTENT(inout) :: myQuadMesh
   INTEGER, INTENT(in)              :: iEdge
   INTEGER, INTENT(in)              :: key

      CALL myQuadMesh % edges(iEdge) % SetKey( key )    

 END SUBROUTINE SetEdgeKey_QuadMesh
!
!
!
 SUBROUTINE GetEdgeKey_QuadMesh( myQuadMesh, iEdge, key )
 ! S/R GetEdgeKey
 !  
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS( QuadMesh ), INTENT(in) :: myQuadMesh
   INTEGER, INTENT(in)           :: iEdge
   INTEGER, INTENT(out)          :: key

      CALL myQuadMesh % edges(iEdge) % GetKey( key )

 END SUBROUTINE GetEdgeKey_QuadMesh
!
!
!
 SUBROUTINE SetEdgeNodeIDs_QuadMesh( myQuadMesh, iEdge, nodeIDs )
 ! S/R SetEdgeNodeIDs
 !  
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS( QuadMesh ), INTENT(inout) :: myQuadMesh
   INTEGER, INTENT(in)              :: iEdge
   INTEGER, INTENT(in)              :: nodeIDs(1:2)

      CALL myQuadMesh % edges(iEdge) % SetNodeIDs( nodeIds )

 END SUBROUTINE SetEdgeNodeIDs_QuadMesh
!
!
!
 SUBROUTINE GetEdgeNodeIDs_QuadMesh( myQuadMesh, iEdge, nodeIDs )
 ! S/R GetEdgeNodeIDs
 !  
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS( QuadMesh ), INTENT(in) :: myQuadMesh
   INTEGER, INTENT(in)           :: iEdge
   INTEGER, INTENT(out)          :: nodeIDs(1:2)

      CALL myQuadMesh % edges(iEdge) % GetNodeIDs( nodeIDs )

 END SUBROUTINE GetEdgeNodeIDs_QuadMesh
!
!
!
 SUBROUTINE SetEdgeElementIDs_QuadMesh( myQuadMesh, iEdge, elementIDs )
 ! S/R SetEdgeElementIDs
 !  
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS( QuadMesh ), INTENT(inout) :: myQuadMesh
   INTEGER, INTENT(in)              :: iEdge
   INTEGER, INTENT(in)              :: elementIDs(1:2)

      CALL myQuadMesh % edges(iEdge) % SetElementIDs( elementIDs )

 END SUBROUTINE SetEdgeElementIDs_QuadMesh
!
!
!
 SUBROUTINE GetEdgeElementIDs_QuadMesh( myQuadMesh, iEdge, elementIDs )
 ! S/R GetEdgeElementIDs
 !  
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS( QuadMesh ), INTENT(in) :: myQuadMesh
   INTEGER, INTENT(in)           :: iEdge
   INTEGER, INTENT(out)          :: elementIDs(1:2)

      CALL myQuadMesh % edges(iEdge) % GetElementIDs( elementIDs )

 END SUBROUTINE GetEdgeElementIDs_QuadMesh
!
!
!
 SUBROUTINE SetEdgePrimaryElementID_QuadMesh( myQuadMesh, iEdge, elementID )
 ! S/R SetEdgePrimaryElementID
 !  
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS( QuadMesh ), INTENT(inout) :: myQuadMesh
   INTEGER, INTENT(in)              :: iEdge
   INTEGER, INTENT(in)              :: elementID

      CALL myQuadMesh % edges(iEdge) % SetPrimaryElementID( elementID )
      
 END SUBROUTINE SetEdgePrimaryElementID_QuadMesh
!
!
!
 SUBROUTINE GetEdgePrimaryElementID_QuadMesh( myQuadMesh, iEdge, elementID )
 ! S/R GetEdgePrimaryElementID
 !  
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS( QuadMesh ), INTENT(in) :: myQuadMesh
   INTEGER, INTENT(in)           :: iEdge
   INTEGER, INTENT(out)          :: elementID

      CALL myQuadMesh % edges(iEdge) % GetPrimaryElementID( elementID )

 END SUBROUTINE GetEdgePrimaryElementID_QuadMesh
!
!
!
 SUBROUTINE SetEdgeSecondaryElementID_QuadMesh( myQuadMesh, iEdge, elementID )
 ! S/R SetEdgeSecondaryElementID
 !  
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS( QuadMesh ), INTENT(inout) :: myQuadMesh
   INTEGER, INTENT(in)              :: iEdge
   INTEGER, INTENT(in)              :: elementID

      CALL myQuadMesh % edges(iEdge) % SetSecondaryElementID( elementID )
      
 END SUBROUTINE SetEdgeSecondaryElementID_QuadMesh
!
!
!
 SUBROUTINE GetEdgeSecondaryElementID_QuadMesh( myQuadMesh, iEdge, elementID )
 ! S/R GetEdgeSecondaryElementID
 !  
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS( QuadMesh ), INTENT(in) :: myQuadMesh
   INTEGER, INTENT(in)           :: iEdge
   INTEGER, INTENT(out)          :: elementID

      CALL myQuadMesh % edges(iEdge) % GetSecondaryElementID( elementID )

 END SUBROUTINE GetEdgeSecondaryElementID_QuadMesh
!
!
!
 SUBROUTINE SetEdgeElementSides_QuadMesh( myQuadMesh, iEdge, elementSides )
 ! S/R SetEdgeElementSides
 !  
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS( QuadMesh ), INTENT(inout) :: myQuadMesh
   INTEGER, INTENT(in)              :: iEdge
   INTEGER, INTENT(in)              :: elementSides(1:2)

      CALL myQuadMesh % edges(iEdge) % SetElementSides( elementSides )

 END SUBROUTINE SetEdgeElementSides_QuadMesh
!
!
!
 SUBROUTINE GetEdgeElementSides_QuadMesh( myQuadMesh, iEdge, elementSides )
 ! S/R GetEdgeElementSides
 !  
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS( QuadMesh ), INTENT(in) :: myQuadMesh
   INTEGER, INTENT(in)           :: iEdge
   INTEGER, INTENT(out)          :: elementSides(1:2)

      CALL myQuadMesh % edges(iEdge) % GetElementSides( elementSides )

 END SUBROUTINE GetEdgeElementSides_QuadMesh
!
!
!
 SUBROUTINE SetEdgePrimaryElementSide_QuadMesh( myQuadMesh, iEdge, elementSide )
 ! S/R SetEdgePrimaryElementSide
 !  
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS( QuadMesh ), INTENT(inout) :: myQuadMesh
   INTEGER, INTENT(in)              :: iEdge
   INTEGER, INTENT(in)              :: elementSide

      CALL myQuadMesh % edges(iEdge) % SetPrimaryElementSide( elementSide )
      
 END SUBROUTINE SetEdgePrimaryElementSide_QuadMesh
!
!
!
 SUBROUTINE GetEdgePrimaryElementSide_QuadMesh( myQuadMesh, iEdge, elementSide )
 ! S/R GetEdgePrimaryElementSide
 !  
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS( QuadMesh ), INTENT(in) :: myQuadMesh
   INTEGER, INTENT(in)           :: iEdge
   INTEGER, INTENT(out)          :: elementSide

      CALL myQuadMesh % edges(iEdge) % GetPrimaryElementSide( elementSide )

 END SUBROUTINE GetEdgePrimaryElementSide_QuadMesh
!
!
!
 SUBROUTINE SetEdgeSecondaryElementSide_QuadMesh( myQuadMesh, iEdge, elementSide )
 ! S/R SetEdgeSecondaryElementSide
 !  
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS( QuadMesh ), INTENT(inout) :: myQuadMesh
   INTEGER, INTENT(in)              :: iEdge
   INTEGER, INTENT(in)              :: elementSide

      CALL myQuadMesh % edges(iEdge) % SetSecondaryElementSide( elementSide )
      
 END SUBROUTINE SetEdgeSecondaryElementSide_QuadMesh
!
!
!
 SUBROUTINE GetEdgeSecondaryElementSide_QuadMesh( myQuadMesh, iEdge, elementSide )
 ! S/R GetEdgeSecondaryElementSide
 !  
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS( QuadMesh ), INTENT(in) :: myQuadMesh
   INTEGER, INTENT(in)           :: iEdge
   INTEGER, INTENT(out)          :: elementSide

      CALL myQuadMesh % edges(iEdge) % GetSecondaryElementSide( elementSide )

 END SUBROUTINE GetEdgeSecondaryElementSide_QuadMesh
!
!
!
 SUBROUTINE SetEdgeStart_QuadMesh( myQuadMesh, iEdge, start )
 ! S/R SetEdgeStart
 !  
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS( QuadMesh ), INTENT(inout) :: myQuadMesh
   INTEGER, INTENT(in)              :: iEdge
   INTEGER, INTENT(in)              :: start
    
      CALL myQuadMesh % edges(iEdge) % SetStart( start )
      
 END SUBROUTINE SetEdgeStart_QuadMesh
!
!
!
 SUBROUTINE GetEdgeStart_QuadMesh( myQuadMesh, iEdge, start )
 ! S/R GetEdgeStart
 !  
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS( QuadMesh ), INTENT(in) :: myQuadMesh
   INTEGER, INTENT(in)           :: iEdge
   INTEGER, INTENT(out)          :: start

      CALL myQuadMesh % edges(iEdge) % GetStart( start )

 END SUBROUTINE GetEdgeStart_QuadMesh
!
!
!
 SUBROUTINE SetEdgeIncrement_QuadMesh( myQuadMesh, iEdge, inc )
 ! S/R SetEdgeIncrement
 !  
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS( QuadMesh ), INTENT(inout) :: myQuadMesh
   INTEGER, INTENT(in)              :: iEdge
   INTEGER, INTENT(in)              :: inc

      CALL myQuadMesh % edges(iEdge) % SetIncrement( inc )

 END SUBROUTINE SetEdgeIncrement_QuadMesh
!
!
!
 SUBROUTINE GetEdgeIncrement_QuadMesh( myQuadMesh, iEdge, inc )
 ! S/R GetEdgeIncrement
 !  
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS( QuadMesh ), INTENT(in) :: myQuadMesh
   INTEGER, INTENT(in)           :: iEdge
   INTEGER, INTENT(out)          :: inc

      CALL myQuadMesh % edges(iEdge) % GetIncrement( inc )

 END SUBROUTINE GetEdgeIncrement_QuadMesh
!
!
!
!==================================================================================================!
!--------------------------------- Type Specific Routines -----------------------------------------!
!==================================================================================================!
!
!
SUBROUTINE ConstructEdges_QuadMesh( myQuadMesh )
 ! S/R ConstructEdges
 !
 !    Takes in the mesh of quadrilaterals which has not filled in the 
 !    edge information, and finds all of the unique edges in the mesh.
 !    The space for the edges is REALlocated with the correct number of edges
 ! 
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS( QuadMesh ), INTENT(inout) :: myQuadMesh
   ! LOCAL
   TYPE( HashTable ) :: edgeTable
   INTEGER :: nEls, nNodes, iEl, nEdges, k  
   INTEGER :: l1, l2, startID, endID, key1, key2
   INTEGER :: e1, s1, edgeID, n1

      CALL myQuadMesh % GetNumberOfNodes( nNodes )   
      CALL myQuadMesh % GetNumberOfElements( nEls )
      nEdges = 0

      ! First, just count the number of edges
      CALL edgeTable % Build( nNodes )

      do iEl = 1, nEls ! Loop over the elements in the mesh

         do k = 1, 4 ! Loop over the sides of each element

            l1 = myQuadMesh % edgeMap(1,k) ! starting local node for this edge
            l2 = myQuadMesh % edgeMap(2,k) ! ending local node for this edge
            
            CALL myQuadMesh % GetElementNodeID( iEl, l1, startID )
            CALL myQuadMesh % GetElementNodeID( iEl, l2, endID )
            
            key1 = min( startID, endID )
            key2 = max( startID, endID )

            ! Add element to corner-node connectivity list
            IF( edgeTable % ContainsKeys( key1, key2 ) .EQV. .FALSE. )then ! this is a new edge
               
               ! Add the edge to the list
               nEdges = nEdges + 1
               CALL edgeTable % AddDataForKeys( nEdges, key1, key2 )
               
            ENDIF

         enddo ! k, Loop over the sides of each element
        
      enddo ! iEl, Loop over the elements in the mesh
 
      CALL edgeTable % Trash( ) ! TRASH the edgetable

      ! And rebuild it
      CALL edgeTable % Build( nNodes )
      
      ! Re-allocate space for the mesh edges

      DEALLOCATE( myQuadMesh % edges )

      ALLOCATE( myQuadMesh % edges( 1:nEdges ) )

      nEdges = 0 ! restart the edge counting

      do iEl = 1, nEls ! Loop over the elements in the mesh

         do k = 1, 4 ! Loop over the sides of each element

            l1 = myQuadMesh % edgeMap(1,k) ! starting local node for this edge
            l2 = myQuadMesh % edgeMap(2,k) ! ending local node for this edge

            CALL myQuadMesh % GetElementNodeID( iEl, l1, startID )
            CALL myQuadMesh % GetElementNodeID( iEl, l2, endID )

            key1 = min( startID, endID )
            key2 = max( startID, endID )

            IF( edgeTable % ContainsKeys( key1, key2 )  )then ! this edge already exists
               
               !Get the edgeID
               CALL edgeTable % GetDataForKeys( edgeID, key1, key2 )
               ! Find the primary element and the starting node for this element's edge
               ! This is compared with the secondary element's starting node to infer
               ! the relative orientation of the two elements.
               CALL myQuadMesh % GetEdgePrimaryElementID( edgeID, e1 )
               CALL myQuadMesh % GetEdgePrimaryElementSide( edgeID, s1 )
               
               l1 = myQuadMesh % edgeMap(1,s1)
               CALL myQuadMesh % GetElementNodeID( e1, l1, n1 )

               ! Set the secondary element information
               CALL myQuadMesh % SetEdgeSecondaryElementID( edgeID, iEl )
               
               !PRINT*, 'edgeID, primary, secondary:', edgeID, e1, iEl
               IF( startID == n1 ) then ! the elements are oriented the same direction
                  CALL myQuadMesh % SetEdgeSecondaryElementSide( edgeID, k )
                 

               ELSE ! the elements are oriented in the opposite direction

                  ! For these edges, we mark the side ID as negative
                  CALL myQuadMesh % SetEdgeSecondaryElementSide( edgeID, -k )
                 
               ENDIF

            ELSE ! this is a new edge

               ! Add the edge to the list
               nEdges = nEdges + 1
               
               edgeID = nEdges
               CALL myQuadMesh % edges(edgeID) % Build()
               CALL myQuadMesh % SetEdgePrimaryElementID( edgeID, iEl )
               CALL myQuadMesh % SetEdgePrimaryElementSide( edgeID, k )
               CALL myQuadMesh % SetEdgeNodeIDs( edgeID, (/startID, endID /) )

               ! Default the secondary information
               CALL myQuadMesh % SetEdgeSecondaryElementID( edgeID, BoundaryFlagDefault )
               
               CALL edgeTable % AddDataForKeys( edgeID, key1, key2 )
               
            ENDIF
         enddo ! k, Loop over the sides of each element
        
      enddo ! iEl, Loop over the elements in the mesh

      CALL edgeTable % Trash( )
      myQuadMesh % nEdges = nEdges

     ! do edgeID = 1, nEdges
      
     !    CALL myQuadMesh % GetEdgePrimaryElementID( edgeID, e1 )
     !    CALL myQuadMesh % GetEdgePrimaryElementSide( edgeID, s1 )
     !    CALL myQuadMesh % GetEdgeSecondaryElementID( edgeID, e2 )
     !    CALL myQuadMesh % GetEdgeSecondaryElementSide( edgeID, s2 )
         
     !    PRINT*, edgeID, e1, s1, e2, s2
     ! enddo
 END SUBROUTINE ConstructEdges_QuadMesh
!
!
!
 SUBROUTINE GetNodeToElementConnectivity_QuadMesh( myQuadMesh )
 ! S/R GetNodeToElementConnectivity
 ! 
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS( QuadMesh ), INTENT(inout) :: myQuadMesh
   ! LOCAL
   INTEGER :: nEls, nNodes, iEl, k  
   INTEGER :: nID

      CALL myQuadMesh % GetNumberOfNodes( nNodes )   
      CALL myQuadMesh % GetNumberOfElements( nEls )

      DO iEl = 1, nEls ! Loop over the elements in the mesh

         DO k = 1, 4 ! Loop over the nodes

            ! Add element to corner-node connectivity list
            CALL myQuadMesh % GetElementNodeID( iEl, k, nID )
            CALL myQuadMesh % nodes(nID) % nodeToElement % AddToList( iEl, k ) ! the key is  the local node ID, the data is the element ID

         ENDDO ! k, Loop over the nodes
        
      ENDDO ! iEl, Loop over the elements in the mesh

     
 END SUBROUTINE GetNodeToElementConnectivity_QuadMesh
!
!
!
 SUBROUTINE ConstructElementNeighbors_QuadMesh( myQuadMesh )
 ! S/R ConstructElementNeighbors
 ! 
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS( QuadMesh ), INTENT(inout) :: myQuadMesh
   ! LOCAL
   INTEGER :: e(1:2), s(1:2), nEdges, iEdge
  
      CALL myQuadMesh % GetNumberOfEdges( nEdges )
 
      DO iEdge = 1, nEdges ! Loop over the edges

         CALL myQuadMesh % GetEdgeElementIDs( iEdge, e )
         CALL myQuadMesh % GetEdgeElementSides( iEdge, s )

         CALL myQuadMesh % SetElementNeighbor( e(1), s(1), e(2) )
         CALL myQuadMesh % SetElementNeighbor( e(2), ABS(s(2)), e(1) )
         
      ENDDO ! iEdge, Loop over the edges

     
 END SUBROUTINE ConstructElementNeighbors_QuadMesh
!
!
!
 SUBROUTINE ScaleTheMesh_QuadMesh( myQuadMesh, interp, xScale, yScale  )
 ! S/R ScaleTheMesh
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS( QuadMesh ), INTENT(inout) :: myQuadMesh
   TYPE( Lagrange_2D ), INTENT(in)  :: interp
   REAL(prec), INTENT(in)           :: xScale, yScale
   ! LOCAL
   INTEGER :: nElems, nNodes, iEl, iNode
   
      CALL myQuadMesh % GetNumberOfElements( nElems )
      CALL myQuadMesh % GetNumberOfNodes( nNodes )

      DO iEl = 1, nElems
         CALL myQuadMesh % elements(iEl) % ScaleGeometry( interp, xScale, yScale )
      ENDDO
      
      DO iNode = 1, nNodes
         CALL myQuadMesh % nodes(iNode) % ScaleNodePosition( xScale, yScale )
      ENDDO

 END SUBROUTINE ScaleTheMesh_QuadMesh
!
!
!
 SUBROUTINE LoadDefaultMesh_QuadMesh( myQuadMesh, interp, nXelem, nYelem  )
 ! S/R LoadDefaultMesh
 !  
 !     Builds a square mesh with 5x5 elements. The side lengths are 1 unit.
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS( QuadMesh ), INTENT(inout)  :: myQuadMesh
   TYPE( Lagrange_2D ), INTENT(in)   :: interp
   INTEGER, INTENT(in)               :: nXelem, nYelem
   ! LOCAL
   TYPE( Curve_2D ) :: elBoundCurves(1:4)

   REAL(prec) :: x, y, dxElem, dyElem
   REAL(prec) :: x1, x2, y1, y2
   REAL(prec), ALLOCATABLE :: xc(:), yc(:), s(:)

   INTEGER :: nNodes, nElems, nEdges, gPolyDeg
   INTEGER :: nodes(1:4)
   INTEGER :: s2, n1, n2
   INTEGER :: iEdge, iNode, iEl, iSide, iX, iY
      
      dxElem = ONE/nXElem
      dyElem = ONE/nYElem
      
      ! ** "Hard-wired" values for a structured mesh with no holes ** !
      nNodes   = (nXElem+1)*(nYElem+1)
      nElems   = (nXElem)*(nYElem)
      nEdges   = (nXElem)*(nYElem+1) + (nXElem+1)*(nYElem)
      gPolyDeg = 1
      ! ************************************************************************* !

      PRINT*, 'nNodes    : ', nNodes
      PRINT*, 'nElems    : ', nElems
      PRINT*, 'gPolyDeg  : ', gPolyDeg

      ! Generate the chebyshev points of order gPolyDeg
      ! These are the points used to define the parametric
      ! curves for the element boundaries

      ALLOCATE( s(0:gPolyDeg), xc(0:gPolyDeg), yc(0:gPolyDeg) )
      CALL ChebyshevGaussLobatto( gPolyDeg, s, xc ) ! ** x is a placeholder here only. 

      ! ---- Build the quadrature mesh (empty) ---- !
  
      CALL myQuadMesh % Build( nNodes, nElems, nEdges, interp % nS ) 
      
      
      ! ---- Read in the corner nodes ---- !
    !  a = 2.5_prec
    !  b = -3.0_prec*HALF*a
    !  c = ONE-(a+b)
      DO iY = 1, nYElem+1
      
         y = dYElem*(REAL(iY-1,prec))
    !     y = a*y**3 + b*y**2 + c*y
         DO iX = 1, nXElem+1
            
            iNode = iX + (iY-1)*(nXElem+1)
            x = dXElem*(REAL(iX-1,prec))
            CALL myQuadMesh % SetNodePosition( iNode, x, y )
            
         ENDDO
         
      ENDDO
  
      ! Do the element information
 
      xc = ZERO
      yc = ZERO
      ! Do the initial build for the parametric curves
      DO iSide = 1, 4
         CALL elBoundCurves(iSide) % Build( xc, yc, s ) 
      ENDDO
   
     
      DO iY = 1, nYElem
         DO iX = 1, nXElem

            iEl = iX + (iY-1)*(nXElem)
            ! Calculate the global node IDs for this element.
            nodes(1) = iX + (iY-1)*(nXElem+1)    ! Southwest
            nodes(2) = iX + 1 + (iY-1)*(nXElem+1)! SouthEast
            nodes(3) = iX + 1 + (iY)*(nXElem+1)  ! NorthEast
            nodes(4) = iX + (iY)*(nXElem+1)      ! NorthWest
         
            DO iSide = 1, 4 ! Loop over the sides of the quads

               ! Build a straight curve for this side
               n1 = nodes( myQuadMesh % edgeMap(1,iSide) )
               n2 = nodes( myQuadMesh % edgeMap(2,iSide) )

               CALL myQuadMesh % GetNodePosition( n1, x1, y1 )
               CALL myQuadMesh % GetNodePosition( n2, x2, y2 )
               
               DO iNode = 0, gPolyDeg
                  xc(iNode) = x1 + (x2-x1)*HALF*( s(iNode) + ONE )
                  yc(iNode) = y1 + (y2-y1)*HALF*( s(iNode) + ONE )
               ENDDO
               CALL elBoundCurves(iSide) % SetNodes( xc, yc ) 

            ENDDO
            CALL myQuadMesh % elements(iEl) % Build( nodes, iEl, elBoundCurves, interp )
         ENDDO
      ENDDO ! iEl, cycle over the elements

      CALL myQuadMesh % ConstructEdges( )
      nEdges = myQuadMesh % nEdges
      PRINT*, 'nEdges    : ', nEdges
      
      ! Set the start and increment for the secondary element 
      DO iEdge = 1, nEdges
            CALL myQuadMesh % GetEdgeSecondaryElementSide(iEdge, s2)
            IF(s2 < 0)THEN
               CALL myQuadMesh % SetEdgeStart( iEdge, interp % nS-1 )
               CALL myQuadMesh % SetEdgeIncrement( iEdge, -1)
            ELSE
               CALL myQuadMesh % SetEdgeStart( iEdge, 1 )
               CALL myQuadMesh % SetEdgeIncrement( iEdge, 1 )
            ENDIF
            
      ENDDO


      CALL myQuadMesh % GetNodeToElementConnectivity( )

      ! Clear up memory
      DEALLOCATE( s, xc, yc )

      DO iSide = 1, 4
         CALL elBoundCurves(iSide) % Trash( )
      ENDDO
  
 END SUBROUTINE LoadDefaultMesh_QuadMesh
!
!
!==================================================================================================!
!--------------------------------- Mesh File I/O Routines -----------------------------------------!
!==================================================================================================!
!
!
 SUBROUTINE ReadSpecMeshFile_QuadMesh( myQuadMesh, interp, filename )
 ! S/R ReadSpecMeshFile
 !  Desription:
 !    Reads in the ISM-v2 version of the SpecMesh output file and builds the mesh with it's geometry 
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS( QuadMesh ), INTENT(out)  :: myQuadMesh
   TYPE( Lagrange_2D ), INTENT(in) :: interp
   CHARACTER(*), INTENT(in)        :: filename
   ! LOCAL
   TYPE( Curve_2D ) :: elBoundCurves(1:4)

   REAL(prec) :: x, y, z
   REAL(prec) :: x1, x2, y1, y2
   REAL(prec), ALLOCATABLE :: xc(:), yc(:), s(:)

   INTEGER :: nNodes, nElems, nEdges, gPolyDeg
   INTEGER :: bFlags(1:4), nodes(1:4)
   INTEGER :: e1, s1, s2, n1, n2, n(1:2), e(1:2), si(1:2)
   INTEGER :: iEdge, iNode, iEl, iSide
   INTEGER :: fUnit, iC, jC
   INTEGER, ALLOCATABLE :: sideFlags(:,:)

   CHARACTER(20) :: ISMversion
   CHARACTER(40) :: edgeNames
   CHARACTER(40) :: thisEdge

      PRINT*, 'Mesh File : '//TRIM( filename )
      
      ! Get a new file unit
      OPEN( UNIT=NEWUNIT(fUnit), &
            FILE= TRIM( filename ), &
            FORM='FORMATTED',&
            STATUS='OLD' )

      
      ! ---- Read in the file version ----- !

      READ( fUnit, '(A20)' ) ISMversion

      PRINT*, 'Version   : '//TRIM( ISMversion )


      ! ---- Gather the number of nodes, number of elements, and number of edges ---- !
      
      READ( fUnit, * ) nNodes, nEdges, nElems, gPolyDeg

      PRINT*, 'nNodes    : ', nNodes
      PRINT*, 'nElems    : ', nElems
      PRINT*, 'nEdges    : ', nEdges
      PRINT*, 'gPolyDeg  : ', gPolyDeg

      ! Generate the chebyshev points of order gPolyDeg
      ! These are the points used to define the parametric
      ! curves for the element boundaries

      ALLOCATE( s(0:gPolyDeg), xc(0:gPolyDeg), yc(0:gPolyDeg) )
      ALLOCATE( sideFlags(1:nElems,1:4) )
      CALL ChebyshevGaussLobatto( gPolyDeg, s, xc ) ! ** x is a placeholder here only.

      xc = ZERO
      yc = ZERO  

      ! ---- Build the quadrature mesh (empty) ---- !
  
      CALL myQuadMesh % Build( nNodes, nElems, nEdges, interp % nS ) 
      
      ! ---- Read in the corner nodes ---- !

      DO iNode = 1, nNodes  ! Loop over the nodes in the file
         READ( fUnit, * ) x, y, z
         CALL myQuadMesh % SetNodePosition( iNode, x, y )
      ENDDO


      ! ---- Read in the edge information ---- !

      DO iEdge = 1, nEdges

         READ( fUnit, * ) n, e, si 
         
         !CALL myQuadMesh % SetEdgeNodeIDs( iEdge, n )
         !CALL myQuadMesh % SetEdgeElementIDs( iEdge, e )
         !CALL myQuadMesh % SetEdgeElementSides( iEdge, si )
       
      ENDDO
      

      ! ---- Read in the element information ---- !
 
      xc = ZERO
      yc = ZERO
      ! Do the initial build for the parametric curves
      DO iSide = 1, 4
         CALL elBoundCurves(iSide) % Build( xc, yc, s ) 
      ENDDO
   
     sideFlags = INTERIOR
      DO iEl = 1, nElems

         READ( fUnit, * ) nodes(1:4)
         READ( fUnit, * ) bFlags(1:4)
         
         DO iSide = 1, 4 ! Loop over the sides of the quads

            
            IF( bFlags(iSide) == 0 )then ! this is an interior edge

               ! Build a straight curve for this side
               n1 = nodes( myQuadMesh % edgeMap(1,iSide) )
               n2 = nodes( myQuadMesh % edgeMap(2,iSide) )

               CALL myQuadMesh % GetNodePosition( n1, x1, y1 )
               CALL myQuadMesh % GetNodePosition( n2, x2, y2 )
               
               DO iNode = 0, gPolyDeg

                  xc(iNode) = x1 + (x2-x1)*HALF*( s(iNode) + ONE )
                  yc(iNode) = y1 + (y2-y1)*HALF*( s(iNode) + ONE )
  
               ENDDO

               CALL elBoundCurves(iSide) % SetNodes( xc, yc ) 

            ELSEIF( bFlags(iSide) == 1 )then ! this is a boundary edge

                ! Read in the parametric curve
                DO iNode = 0, gPolyDeg

                   READ( fUnit, * ) x, y, z
                   xc(iNode) = x
                   yc(iNode) = y

                ENDDO

                CALL elBoundCurves(iSide) % SetNodes( xc, yc ) 

             ELSE
 
                PRINT*,' S/R ReadSpecMeshFile : Impossible element boundary flag '

             ENDIF

         ENDDO

         ! Build this element's geometry
         CALL myQuadMesh % elements(iEl) % Build( nodes, iEl, elBoundCurves, interp )

         ! Read in and parse the edge names
         READ( fUnit, '(1x, A40)' )  edgeNames

         ! Parse the edge names into four edge names
         iSide = 1
         iC = 1
         jC = 1
         thisEdge = ' ' 
         DO while( iSide <= 4 )

            IF( edgeNames(iC:iC) == ' ' )then ! we have reached a blank space
     
               n1 = nodes( myQuadMesh % edgeMap(1,iSide) )
               n2 = nodes( myQuadMesh % edgeMap(2,iSide) )
               
               IF( TRIM(thisEdge) == 'DIRICHLET' )then
                  
                  sideFlags(iEl,iSide) = DIRICHLET

               ELSEIF( TRIM(thisEdge) == 'ROBIN' )then

                  sideFlags(iEl,iSide) = ROBIN

               ELSEIF( TRIM(thisEdge) == 'ROBIN_FORCED' )then

                  sideFlags(iEl,iSide) = ROBIN_FORCED

               ELSEIF( TRIM(thisEdge) == 'HOMOGENEOUS_NEUMANN' )then

                  sideFlags(iEl,iSide) = HOMOGENEOUS_NEUMANN

               ELSEIF( TRIM(thisEdge) == 'NEUMANN_WALL' )then

                  sideFlags(iEl,iSide) = NEUMANN_WALL

               ELSEIF( TRIM(thisEdge) == 'NEUMANN' )then

                  sideFlags(iEl,iSide) = NEUMANN
 
               ELSEIF( TRIM(thisEdge) == 'NO_NORMAL_FLOW' )then
                  
                  sideFlags(iEl,iSide) = NO_NORMAL_FLOW

               ELSEIF( TRIM(thisEdge) == 'PRESCRIBED' )then
               
                  sideFlags(iEl,iSide) = PRESCRIBED

               ELSEIF( TRIM(thisEdge) == 'RADIATION' )then
                  
                  sideFlags(iEl,iSide) = RADIATION

               ELSEIF( TRIM(thisEdge) == 'InflowOne' )then
                  
                  sideFlags(iEl,iSide) = InflowOne

               ELSEIF( TRIM(thisEdge) == 'InflowTwo' )then
                  
                  sideFlags(iEl,iSide) = InflowTwo

               ENDIF

               ! Reset thisEdge
               thisEdge = ' '
               jC = 1 
               iC = iC + 1
               ! Increment the side
               iSide = iSide + 1

            ELSE
            
               thisEdge(jC:jC) = edgeNames(iC:iC)
               iC = iC + 1
               jC = jC + 1

            ENDIF

         ENDDO ! while (iSide <= 4 )


      ENDDO ! iEl, cycle over the elements

      close( fUnit )
            CALL myQuadMesh % ConstructEdges( )
      ! Set the secondary element on boundary edges to the boundary flag
      DO iEdge = 1, nEdges

         CALL myQuadMesh % GetEdgePrimaryElementID( iEdge, e1 )
         CALL myQuadMesh % GetEdgePrimaryElementSide( iEdge, s1 )
      
         IF( sideFlags(e1,s1) == DIRICHLET )then

            CALL myQuadMesh % SetEdgeSecondaryElementID( iEdge, DIRICHLET )
            CALL myQuadMesh % GetEdgeNodeIDs( iEdge, n )
            CALL myQuadMesh % SetNodeType( n(1), DIRICHLET )
            CALL myQuadMesh % SetNodeType( n(2), DIRICHLET )

         ELSEIF( sideFlags(e1,s1) == ROBIN )then

            CALL myQuadMesh % SetEdgeSecondaryElementID( iEdge, ROBIN )
            CALL myQuadMesh % GetEdgeNodeIDs( iEdge, n )
            CALL myQuadMesh % SetNodeType( n(1), ROBIN )
            CALL myQuadMesh % SetNodeType( n(2), ROBIN )

         ELSEIF( sideFlags(e1,s1) == ROBIN_FORCED )then

            CALL myQuadMesh % SetEdgeSecondaryElementID( iEdge, ROBIN_FORCED )
            CALL myQuadMesh % GetEdgeNodeIDs( iEdge, n )
            CALL myQuadMesh % SetNodeType( n(1), ROBIN_FORCED )
            CALL myQuadMesh % SetNodeType( n(2), ROBIN_FORCED )

         ELSEIF( sideFlags(e1,s1) == HOMOGENEOUS_NEUMANN )then

            CALL myQuadMesh % SetEdgeSecondaryElementID( iEdge, HOMOGENEOUS_NEUMANN )
            CALL myQuadMesh % GetEdgeNodeIDs( iEdge, n )
            CALL myQuadMesh % SetNodeType( n(1), HOMOGENEOUS_NEUMANN )
            CALL myQuadMesh % SetNodeType( n(2), HOMOGENEOUS_NEUMANN )

         ELSEIF( sideFlags(e1,s1) == NEUMANN_WALL )then

            CALL myQuadMesh % SetEdgeSecondaryElementID( iEdge, NEUMANN_WALL )
            CALL myQuadMesh % GetEdgeNodeIDs( iEdge, n )
            CALL myQuadMesh % SetNodeType( n(1), NEUMANN_WALL )
            CALL myQuadMesh % SetNodeType( n(2), NEUMANN_WALL )

         ELSEIF( sideFlags(e1,s1) == NEUMANN )then

            CALL myQuadMesh % SetEdgeSecondaryElementID( iEdge, NEUMANN )
            CALL myQuadMesh % GetEdgeNodeIDs( iEdge, n )
            CALL myQuadMesh % SetNodeType( n(1), NEUMANN )
            CALL myQuadMesh % SetNodeType( n(2), NEUMANN )

         ELSEIF( sideFlags(e1,s1) == NO_NORMAL_FLOW )then

            CALL myQuadMesh % SetEdgeSecondaryElementID( iEdge, NO_NORMAL_FLOW )
            CALL myQuadMesh % GetEdgeNodeIDs( iEdge, n )
            CALL myQuadMesh % SetNodeType( n(1), NO_NORMAL_FLOW )
            CALL myQuadMesh % SetNodeType( n(2), NO_NORMAL_FLOW )

         ELSEIF( sideFlags(e1,s1) == PRESCRIBED )then

            CALL myQuadMesh % SetEdgeSecondaryElementID( iEdge, PRESCRIBED )
            CALL myQuadMesh % GetEdgeNodeIDs( iEdge, n )
            CALL myQuadMesh % SetNodeType( n(1), PRESCRIBED )
            CALL myQuadMesh % SetNodeType( n(2), PRESCRIBED )

         ELSEIF( sideFlags(e1,s1) == RADIATION )then

            CALL myQuadMesh % SetEdgeSecondaryElementID( iEdge, RADIATION )
            CALL myQuadMesh % GetEdgeNodeIDs( iEdge, n )
            CALL myQuadMesh % SetNodeType( n(1), RADIATION )
            CALL myQuadMesh % SetNodeType( n(2), RADIATION )

         ELSEIF( sideFlags(e1,s1) == InflowOne )then

            CALL myQuadMesh % SetEdgeSecondaryElementID( iEdge, InflowOne )
            CALL myQuadMesh % GetEdgeNodeIDs( iEdge, n )
            CALL myQuadMesh % SetNodeType( n(1), InflowOne )
            CALL myQuadMesh % SetNodeType( n(2), InflowOne )

         ELSEIF( sideFlags(e1,s1) == InflowTwo )then

            CALL myQuadMesh % SetEdgeSecondaryElementID( iEdge, InflowTwo )
            CALL myQuadMesh % GetEdgeNodeIDs( iEdge, n )
            CALL myQuadMesh % SetNodeType( n(1), InflowTwo )
            CALL myQuadMesh % SetNodeType( n(2), InflowTwo )
         
         ENDIF
 
      ENDDO

      DO iEdge = 1, nEdges
            CALL myQuadMesh % GetEdgeSecondaryElementSide(iEdge, s2)
            IF(s2 < 0)THEN
               CALL myQuadMesh % SetEdgeStart( iEdge, interp % nS-1 )
               CALL myQuadMesh % SetEdgeIncrement( iEdge, -1)
            ELSE
               CALL myQuadMesh % SetEdgeStart( iEdge, 1 )
               CALL myQuadMesh % SetEdgeIncrement( iEdge, 1 )
            ENDIF
            
      ENDDO
 ! If a node is at the corner of any boundary condition and a dirichlet boundary condition,
 ! we opt to keep it as a dirichlet boundary condition
    !  DO iEdge = 1, nEdges

    !     CALL myQuadMesh % edges(iEdge) % GetPrimaryElementID( e1 )
    !     CALL myQuadMesh % edges(iEdge) % GetPrimaryElementSide( s1 )

    !     IF( sideFlags(e1,s1) == DIRICHLET )then

    !        CALL myQuadMesh % edges(iEdge) % SetSecondaryElementID( DIRICHLET )
            
    !        CALL myQuadMesh % edges(iEdge) % GetNodeIDs( n )
    !        CALL myQuadMesh % nodes(n(1)) % SetType( DIRICHLET )
    !        CALL myQuadMesh % nodes(n(2)) % SetType( DIRICHLET )

    !     ENDIF
 
    !  ENDDO
      ! Get the node to element connectivity
      CALL myQuadMesh % GetNodeToElementConnectivity( )

      ! Clear up memory
      DEALLOCATE( s, xc, yc, sideFlags )

      DO iSide = 1, 4
         CALL elBoundCurves(iSide) % Trash( )
      ENDDO
  
 END SUBROUTINE ReadSpecMeshFile_QuadMesh
!
!
!
 SUBROUTINE WriteTecplot_Quadmesh( myQuadMesh, filename )
 ! S/R WriteTecplot
 !
 ! =============================================================================================== !
 ! DECLARATIONS
  IMPLICIT NONE
  CLASS(QuadMesh), INTENT(inout)     :: myQuadMesh
  CHARACTER(*), INTENT(in), OPTIONAL :: filename  
  ! Local
  INTEGER :: iS, iP, nS, nP, iEl,fUnit, eID
  REAL(prec) :: x, y, J, dxds, dxdp, dyds, dydp
  CHARACTER(7) :: zoneID

    CALL myQuadMesh % elements(1) % GetNumberOfNodes( nS, nP )

    IF( PRESENT(filename) )THEN
    
       OPEN( UNIT=NEWUNIT(fUnit), &
             FILE= TRIM(filename)//'.tec', &
             FORM='formatted')
    ELSE
    
       OPEN( UNIT=NEWUNIT(fUnit), &
             FILE= 'mesh.tec', &
             FORM='formatted')

    ENDIF
    
    WRITE(fUnit,*) 'VARIABLES = "X", "Y", "Jacobian", "dxds", "dxdp", "dyds", "dydp" '


    DO iEl = 1, myQuadMesh % nElems

       CALL myQuadMesh % elements(iEl) % GetElementID( eID )
       WRITE(zoneID,'(I7.7)') eID
       WRITE(fUnit,*) 'ZONE T="el'//trim(zoneID)//'", I=',nS+1,', J=', nP+1,',F=POINT'

       DO iP = 0, nP
          DO iS = 0,nS
             CALL myQuadMesh % GetPositionAtNode( iEl, x, y, iS, iP )
             CALL myQuadMesh % GetJacobianAtNode( iEl, J, iS, iP )
             CALL myQuadMesh % GetCovariantMetricsAtNode( iEl, dxds, dxdp, dyds, dydp, iS, iP )
             WRITE(fUnit,*)  x, y, J, dxds, dxdp, dyds, dydp
          ENDDO
      ENDDO

    ENDDO
    
    CLOSE(UNIT=fUnit)

    RETURN

 END SUBROUTINE WriteTecplot_Quadmesh
!
!
!
 

END MODULE QuadMeshClass
