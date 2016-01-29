! HexMeshClass.f90 ( new with v2.1 - 27 Jan. 2016)
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



MODULE HexMeshClass


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
USE HexElementClass
USE EdgeClass
USE NodeClass_2D



IMPLICIT NONE


! All properties are left as public in version 2.1 in order to allow the programmer the ability to
! have direct access to the attributes of the mesh class. This way, it can later be determined if
! additional accessor routines would be necessary.
 
    TYPE HexMesh 
       INTEGER                          :: nElems, nNodes, nEdges
       TYPE( HexElement ), ALLOCATABLE  :: elements(:)
       TYPE( Node ), ALLOCATABLE        :: nodes(:)  
       TYPE( Edge ), ALLOCATABLE        :: edges(:)
       INTEGER                          :: cornerMap(1:2,1:nHexNodes) 
       INTEGER                          :: sideMap(1:nHexFaces) 
       INTEGER                          :: faceMap(1:2,1:nHexFaces) 

       CONTAINS

       PROCEDURE :: Build => Build_HexMesh
       PROCEDURE :: Trash => Trash_HexMesh
       
       PROCEDURE :: SetNumberOfElements => SetNumberOfElements_HexMesh
       PROCEDURE :: GetNumberOfElements => GetNumberOfElements_HexMesh
       PROCEDURE :: SetNumberOfNodes => SetNumberOfNodes_HexMesh
       PROCEDURE :: GetNumberOfNodes => GetNumberOfNodes_HexMesh
       PROCEDURE :: SetNumberOfEdge => SetNumberOfEdges_HexMesh
       PROCEDURE :: GetNumberOfEdges => GetNumberOfEdges_HexMesh
       
       !HexElement Wrapper Routines
       PROCEDURE :: SetElementNodeIDs => SetElementNodeIDs_HexMesh
       PROCEDURE :: GetElementNodeIDs => GetElementNodeIDs_HexMesh
       PROCEDURE :: SetElementNodeID => SetElementNodeID_HexMesh
       PROCEDURE :: GetElementNodeID => GetElementNodeID_HexMesh
       PROCEDURE :: SetElementID => SetElementID_HexMesh
       PROCEDURE :: GetElementID => GetElementID_HexMesh
       PROCEDURE :: SetElementNeighbor => SetElementNeighbor_HexMesh
       PROCEDURE :: GetElementNeighbor => GetElementNeighbor_HexMesh
       PROCEDURE :: SetElementSouthernNeighbor => SetElementSouthernNeighbor_HexMesh
       PROCEDURE :: GetElementSouthernNeighbor => GetElementSouthernNeighbor_HexMesh
       PROCEDURE :: SetElementNorthernNeighbor => SetElementNorthernNeighbor_HexMesh
       PROCEDURE :: GetElementNorthernNeighbor => GetElementNorthernNeighbor_HexMesh
       PROCEDURE :: SetElementEasternNeighbor  => SetElementEasternNeighbor_HexMesh
       PROCEDURE :: GetElementEasternNeighbor  => GetElementEasternNeighbor_HexMesh
       PROCEDURE :: SetElementWesternNeighbor  => SetElementWesternNeighbor_HexMesh
       PROCEDURE :: GetElementWesternNeighbor  => GetElementWesternNeighbor_HexMesh
      
       ! MappedGeometry_2D Wrapper Routines
       PROCEDURE :: SetNumberOfInternalNodes => SetNumberOfInternalNodes_HexMesh
       PROCEDURE :: GetNumberOfInternalNodes => GetNumberOfInternalNodes_HexMesh
       PROCEDURE :: SetPositions => SetPositions_HexMesh
       PROCEDURE :: GetPositions => GetPositions_HexMesh
       PROCEDURE :: SetPositionAtNode => SetPositionAtNode_HexMesh
       PROCEDURE :: GetPositionAtNode => GetPositionAtNode_HexMesh
       PROCEDURE :: SetJacobian => SetJacobian_HexMesh
       PROCEDURE :: GetJacobian => GetJacobian_HexMesh
       PROCEDURE :: SetJacobianAtNode => SetJacobianAtNode_HexMesh
       PROCEDURE :: GetJacobianAtNode => GetJacobianAtNode_HexMesh
       PROCEDURE :: SetCovariantMetrics => SetCovariantMetrics_HexMesh
       PROCEDURE :: GetCovariantMetrics => GetCovariantMetrics_HexMesh
       PROCEDURE :: SetCovariantMetricsAtNode => SetCovariantMetricsAtNode_HexMesh
       PROCEDURE :: GetCovariantMetricsAtNode => GetCovariantMetricsAtNode_HexMesh
       PROCEDURE :: SetBoundaryLocation => SetBoundaryLocation_HexMesh
       PROCEDURE :: GetBoundaryLocation => GetBoundaryLocation_HexMesh
       PROCEDURE :: SetBoundaryLocationAtNode => SetBoundaryLocationAtNode_HexMesh
       PROCEDURE :: GetBoundaryLocationAtNode => GetBoundaryLocationAtNode_HexMesh
       PROCEDURE :: SetBoundaryNormalAtNode => SetBoundaryNormalAtNode_HexMesh
       PROCEDURE :: GetBoundaryNormalAtNode => GetBoundaryNormalAtNode_HexMesh
       
       ! NodeClass_2D Wrapper Routines
       PROCEDURE :: SetNodeData => SetNodeData_HexMesh
       PROCEDURE :: GetNodeData => GetNodeData_HexMesh
       PROCEDURE :: SetNodeKey => SetNodeKey_HexMesh
       PROCEDURE :: GetNodeKey => GetNodeKey_HexMesh
       PROCEDURE :: SetNodeType => SetNodeType_HexMesh
       PROCEDURE :: GetNodeType => GetNodeType_HexMesh
       PROCEDURE :: SetNodePosition => SetNodePosition_HexMesh
       PROCEDURE :: GetNodePosition => GetNodePosition_HexMesh

       ! EdgeClass Wrapper Routines 
       PROCEDURE :: SetEdgeKey => SetEdgeKey_HexMesh
       PROCEDURE :: GetEdgeKey => GetEdgeKey_HexMesh
       PROCEDURE :: SetEdgeData => SetEdgeData_HexMesh
       PROCEDURE :: GetEdgeData => GetEdgeData_HexMesh
       PROCEDURE :: SetEdgeNodeIDs => SetEdgeNodeIDs_HexMesh
       PROCEDURE :: GetEdgeNodeIDs => GetEdgeNodeIDs_HexMesh
       PROCEDURE :: SetEdgeElementIDs => SetEdgeElementIDs_HexMesh
       PROCEDURE :: GetEdgeElementIDs => GetEdgeElementIDs_HexMesh
       PROCEDURE :: SetEdgePrimaryElementID => SetEdgePrimaryElementID_HexMesh
       PROCEDURE :: GetEdgePrimaryElementID => GetEdgePrimaryElementID_HexMesh
       PROCEDURE :: SetEdgeSecondaryElementID => SetEdgeSecondaryElementID_HexMesh
       PROCEDURE :: GetEdgeSecondaryElementID => GetEdgeSecondaryElementID_HexMesh
       PROCEDURE :: SetEdgeElementSides => SetEdgeElementSides_HexMesh
       PROCEDURE :: GetEdgeElementSides => GetEdgeElementSides_HexMesh
       PROCEDURE :: SetEdgePrimaryElementSide => SetEdgePrimaryElementSide_HexMesh
       PROCEDURE :: GetEdgePrimaryElementSide => GetEdgePrimaryElementSide_HexMesh
       PROCEDURE :: SetEdgeSecondaryElementSide => SetEdgeSecondaryElementSide_HexMesh
       PROCEDURE :: GetEdgeSecondaryElementSide => GetEdgeSecondaryElementSide_HexMesh
       PROCEDURE :: SetEdgeStart => SetEdgeStart_HexMesh
       PROCEDURE :: GetEdgeStart => GetEdgeStart_HexMesh
       PROCEDURE :: SetEdgeIncrement => SetEdgeIncrement_HexMesh
       PROCEDURE :: GetEdgeIncrement => GetEdgeIncrement_HexMesh
       
       PROCEDURE :: ConstructEdges => ConstructEdges_HexMesh
       PROCEDURE :: GetNodeToElementConnectivity => GetNodeToElementConnectivity_HexMesh
       PROCEDURE :: ScaleTheMesh => ScaleTheMesh_HexMesh
       PROCEDURE :: LoadDefaultMesh => LoadDefaultMesh_HexMesh
       
       PROCEDURE :: ReadSpecMeshFile => ReadSpecMeshFile_HexMesh
       PROCEDURE :: WriteTecplot => WriteTecplot_Hexmesh
       
    END TYPE HexMesh

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
 SUBROUTINE Build_HexMesh( myHexMesh, nNodes, nElems, nEdges, nS )
 ! S/R Build
 ! 
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS(HexMesh), INTENT(out) :: myHexMesh
   INTEGER, INTENT(in)          :: nNodes, nElems, nEdges, nS
   !LOCAL
   INTEGER :: iNode

      myHexmesh % sideMap(1:4) = (/ 0, nS, nS, 0 /)
      myHexmesh % cornerMap(1, 1:4) = (/ 0, nS, nS, 0 /)
      myHexmesh % cornerMap(2, 1:4) = (/ 0, 0, nS, nS /)
      myHexmesh % faceMap(1, 1:4) = (/ 1, 2, 4, 1 /)
      myHexmesh % faceMap(2, 1:4) = (/ 2, 3, 3, 4 /)

      myHexmesh % nNodes = nNodes
      myHexmesh % nElems = nElems

      ALLOCATE( myHexmesh % elements(1:nElems) )
      ALLOCATE( myHexmesh % nodes(1:nNodes) )
      ALLOCATE( myHexmesh % edges(1:nEdges) )
     

      ! Default nodes to origin
      DO iNode = 1, myHexmesh % nNodes ! loop over the number of nodes
         CALL myHexmesh % nodes(iNode) % Build( ZERO, ZERO )
      ENDDO ! iNode, loop over the number of nodes
     

 END SUBROUTINE Build_HexMesh
!
!
!
 SUBROUTINE Trash_HexMesh( myHexMesh )
 ! S/R Trash
 ! 
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS(HexMesh), INTENT(inout) :: myHexMesh
  ! LOCAL
   INTEGER :: iNode, iEl

      DO iEl = 1, myHexMesh % nElems
         CALL myHexMesh % elements(iEl) % TRASH( )
      ENDDO

      DO iNode = 1, myHexMesh % nNodes
         CALL myHexMesh % nodes(iNode) % TRASH( )
      ENDDO

      DEALLOCATE( myHexMesh % nodes, myHexMesh % elements, myHexMesh % edges )
      

 END SUBROUTINE Trash_HexMesh
!
!
!==================================================================================================!
!--------------------------------------- Accessors ------------------------------------------------!
!==================================================================================================!
!
! 
 SUBROUTINE SetNumberOfElements_HexMesh( myHexMesh, nElems )
 ! S/R SetNumberOfElements
 ! 
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS(HexMesh), INTENT(inout) :: myHexMesh
   INTEGER, INTENT(in)            :: nElems
   
      myHexMesh % nElems = nElems
      
 END SUBROUTINE SetNumberOfElements_HexMesh
!
!
!
 SUBROUTINE GetNumberOfElements_HexMesh( myHexMesh, nElems )
 ! S/R GetNumberOfElements
 ! 
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS(HexMesh), INTENT(in) :: myHexMesh
   INTEGER, INTENT(out)        :: nElems
   
      nElems = myHexMesh % nElems
      
 END SUBROUTINE GetNumberOfElements_HexMesh
!
!
!
 SUBROUTINE SetNumberOfNodes_HexMesh( myHexMesh, nNodes )
 ! S/R SetNumberOfNodes
 ! 
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS(HexMesh), INTENT(inout) :: myHexMesh
   INTEGER, INTENT(in)            :: nNodes
   
      myHexMesh % nNodes = nNodes
      
 END SUBROUTINE SetNumberOfNodes_HexMesh
!
!
!
 SUBROUTINE GetNumberOfNodes_HexMesh( myHexMesh, nNodes )
 ! S/R GetNumberOfNodes
 ! 
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS(HexMesh), INTENT(in) :: myHexMesh
   INTEGER, INTENT(out)        :: nNodes
   
      nNodes = myHexMesh % nNodes
      
 END SUBROUTINE GetNumberOfNodes_HexMesh
!
!
!
 SUBROUTINE SetNumberOfEdges_HexMesh( myHexMesh, nEdges )
 ! S/R SetNumberOfEdges
 ! 
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS(HexMesh), INTENT(inout) :: myHexMesh
   INTEGER, INTENT(in)            :: nEdges
   
      myHexMesh % nEdges = nEdges
      
 END SUBROUTINE SetNumberOfEdges_HexMesh
!
!
!
 SUBROUTINE GetNumberOfEdges_HexMesh( myHexMesh, nEdges )
 ! S/R GetNumberOfEdges
 ! 
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS(HexMesh), INTENT(in) :: myHexMesh
   INTEGER, INTENT(out)        :: nEdges
   
      nEdges = myHexMesh % nEdges
      
 END SUBROUTINE GetNumberOfEdges_HexMesh
!
! ---------------------------- HexElement Wrapper Routines -------------------------------------- !
!
SUBROUTINE SetElementNodeIDs_Hexmesh( myHexMesh, iEl, nodeIDs )
 ! S/R SetElementNodeIDs
 !
 ! =============================================================================================== !
   IMPLICIT NONE
   CLASS( Hexmesh ), INTENT(inout) :: myHexMesh
   INTEGER, INTENT(in)              :: iEl
   INTEGER, INTENT(in)              :: nodeIDs(1:4)
   
      CALL myHexMesh % elements(iEl)  % SetNodeIDs( nodeIDs )

 END SUBROUTINE SetElementNodeIDs_Hexmesh
!
!
!
 SUBROUTINE GetElementNodeIDs_Hexmesh( myHexMesh, iEl, nodeIDs )
 ! S/R GetElementNodeIDs
 !
 ! =============================================================================================== !
   IMPLICIT NONE
   CLASS( Hexmesh ), INTENT(in) :: myHexMesh
   INTEGER, INTENT(in)           :: iEl
   INTEGER, INTENT(out)          :: nodeIDs(1:4)
   
      CALL myHexMesh % elements(iEl)  % GetNodeIDs( nodeIDs )

 END SUBROUTINE GetElementNodeIDs_Hexmesh
!
!
!
 SUBROUTINE SetElementNodeID_Hexmesh( myHexMesh, iEl, localID, nodeID )
 ! S/R SetElementNodeID
 !
 ! =============================================================================================== !
   IMPLICIT NONE
   CLASS( Hexmesh ), INTENT(inout) :: myHexMesh
   INTEGER, INTENT(in)              :: iEl
   INTEGER, INTENT(in)              :: localID, nodeID
   
      CALL myHexMesh % elements(iEl)  % SetNodeID( localID, nodeID )

 END SUBROUTINE SetElementNodeID_Hexmesh
!
!
!
 SUBROUTINE GetElementNodeID_Hexmesh( myHexMesh, iEl, localID, nodeID )
 ! S/R GetElementNodeID
 !
 ! =============================================================================================== !
   IMPLICIT NONE
   CLASS( Hexmesh ), INTENT(in) :: myHexMesh
   INTEGER, INTENT(in)           :: iEl
   INTEGER, INTENT(in)           :: localID
   INTEGER, INTENT(out)          :: nodeID
   
      CALL myHexMesh % elements(iEl) % GetNodeID( localID, nodeID )

 END SUBROUTINE GetElementNodeID_Hexmesh
!
!
!
 SUBROUTINE SetElementID_Hexmesh( myHexMesh, iEl, eID )
 ! S/R SetElementID
 !
 ! =============================================================================================== !
   IMPLICIT NONE
   CLASS( Hexmesh ), INTENT(inout) :: myHexMesh
   INTEGER, INTENT(in)              :: iEl
   INTEGER, INTENT(in)              :: eID
   
      CALL myHexMesh % elements(iEl)  % SetElementID( eID )

 END SUBROUTINE SetElementID_Hexmesh
!
!
!
 SUBROUTINE GetElementID_Hexmesh( myHexMesh, iEl, eID )
 ! S/R GetElementID
 !
 ! =============================================================================================== !
   IMPLICIT NONE
   CLASS( Hexmesh ), INTENT(in) :: myHexMesh
   INTEGER, INTENT(in)           :: iEl
   INTEGER, INTENT(out)          :: eID
   
      CALL myHexMesh % elements(iEl)  % GetElementID( eID )

 END SUBROUTINE GetElementID_Hexmesh
!
!
!
 SUBROUTINE SetElementNeighbor_Hexmesh( myHexMesh, iEl, sID, eID )
  ! S/R SetElementNeighbor
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS( Hexmesh ), INTENT(inout) :: myHexMesh
   INTEGER, INTENT(in)              :: iEl
   INTEGER, INTENT(in)              :: sID
   INTEGER, INTENT(in)              :: eID
   
      CALL myHexMesh % elements(iEl) % SetNeighbor( sID, eID )
      
 END SUBROUTINE SetElementNeighbor_Hexmesh
!
!
!
 SUBROUTINE GetElementNeighbor_Hexmesh( myHexMesh, iEl, sID, eID )
 ! S/R GetElementNeighbor
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS( Hexmesh ), INTENT(in) :: myHexMesh
   INTEGER, INTENT(in)           :: iEl
   INTEGER, INTENT(in)           :: sID
   INTEGER, INTENT(out)          :: eID
   
      CALL myHexMesh % elements(iEl) % GetNeighbor( sID, eID )
      
 END SUBROUTINE GetElementNeighbor_Hexmesh
!
!
!
 SUBROUTINE SetElementSouthernNeighbor_Hexmesh( myHexMesh, iEl, eID )
 ! S/R SetElementSouthernNeighbor
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS( Hexmesh ), INTENT(inout) :: myHexMesh
   INTEGER, INTENT(in)              :: iEl
   INTEGER, INTENT(in)              :: eID
   
      CALL myHexMesh % elements(iEl) % SetSouthernNeighbor( eID )
      
 END SUBROUTINE SetElementSouthernNeighbor_Hexmesh
!
!
!
 SUBROUTINE GetElementSouthernNeighbor_Hexmesh( myHexMesh, iEl, eID )
 ! S/R GetElementSouthernNeighbor
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS( Hexmesh ), INTENT(in) :: myHexMesh
   INTEGER, INTENT(in)           :: iEl
   INTEGER, INTENT(out)          :: eID
   
      CALL myHexMesh % elements(iEl) % GetSouthernNeighbor( eID )
      
 END SUBROUTINE GetElementSouthernNeighbor_Hexmesh
!
!
!
 SUBROUTINE SetElementNorthernNeighbor_Hexmesh( myHexMesh, iEl, eID )
 ! S/R SetElementNorthernNeighbor
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS( Hexmesh ), INTENT(inout) :: myHexMesh
   INTEGER, INTENT(in)              :: iEl
   INTEGER, INTENT(in)              :: eID
   
      CALL myHexMesh % elements(iEl) % SetNorthernNeighbor( eID )
      
 END SUBROUTINE SetElementNorthernNeighbor_Hexmesh
!
!
!
 SUBROUTINE GetElementNorthernNeighbor_Hexmesh( myHexMesh, iEl, eID )
 ! S/R GetElementNorthernNeighbor
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS( Hexmesh ), INTENT(in) :: myHexMesh
   INTEGER, INTENT(in)           :: iEl 
   INTEGER, INTENT(out)          :: eID
   
      CALL myHexMesh % elements(iEl) % GetNorthernNeighbor( eID )
      
 END SUBROUTINE GetElementNorthernNeighbor_Hexmesh
!
!
!
 SUBROUTINE SetElementEasternNeighbor_Hexmesh( myHexMesh, iEl, eID )
 ! S/R SetElementEasternNeighbor
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS( Hexmesh ), INTENT(inout) :: myHexMesh
   INTEGER, INTENT(in)              :: iEl
   INTEGER, INTENT(in)              :: eID
   
      CALL myHexMesh % elements(iEl) % SetEasternNeighbor( eID )
      
 END SUBROUTINE SetElementEasternNeighbor_Hexmesh
!
!
!
 SUBROUTINE GetElementEasternNeighbor_Hexmesh( myHexMesh, iEl, eID )
 ! S/R GetElementEasternNeighbor
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS( Hexmesh ), INTENT(in) :: myHexMesh
   INTEGER, INTENT(in)           :: iEl
   INTEGER, INTENT(out)          :: eID
   
      CALL myHexMesh % elements(iEl) % GetEasternNeighbor( eID )
      
 END SUBROUTINE GetElementEasternNeighbor_Hexmesh
!
!
!
 SUBROUTINE SetElementWesternNeighbor_Hexmesh( myHexMesh, iEl, eID )
 ! S/R SetElementWesternNeighbor
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS( Hexmesh ), INTENT(inout) :: myHexMesh
   INTEGER, INTENT(in)              :: iEl
   INTEGER, INTENT(in)              :: eID
   
      CALL myHexMesh % elements(iEl) % SetWesternNeighbor( eID )
      
 END SUBROUTINE SetElementWesternNeighbor_Hexmesh
!
!
!
 SUBROUTINE GetElementWesternNeighbor_Hexmesh( myHexMesh, iEl, eID )
 ! S/R GetElementWesternNeighbor
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS( Hexmesh ), INTENT(in) :: myHexMesh
   INTEGER, INTENT(in)           :: iEl
   INTEGER, INTENT(out)          :: eID
   
      CALL myHexMesh % elements(iEl) % GetWesternNeighbor( eID )
      
 END SUBROUTINE GetElementWesternNeighbor_Hexmesh
!
! ---------------------------- MappedGeometryClass_2D Wrappers ----------------------------------- !
!
SUBROUTINE SetNumberOfInternalNodes_Hexmesh( myHexMesh, iEl, nS, nP )
 ! S/R SetNumberOfInternalNodes
 !
 ! =============================================================================================== !
 ! DECLARATIONS
  IMPLICIT NONE
  CLASS(Hexmesh), INTENT(inout) :: myHexMesh
  INTEGER, INTENT(in)            :: iEl
  INTEGER, INTENT(in)            :: nS, nP
  
     CALL myHexMesh % elements(iEl) % SetNumberOfNodes( nS, nP )
     
 END SUBROUTINE SetNumberOfInternalNodes_Hexmesh
!
!
!
 SUBROUTINE GetNumberOfInternalNodes_Hexmesh( myHexMesh, iEl, nS, nP )
 ! S/R GetNumberOfInternalNodes
 !
 ! =============================================================================================== !
 ! DECLARATIONS
  IMPLICIT NONE
  CLASS(Hexmesh), INTENT(in) :: myHexMesh
  INTEGER, INTENT(in)         :: iEl
  INTEGER, INTENT(out)        :: nS, nP
  
     CALL myHexMesh % elements(iEl) % GetNumberOfNodes( nS, nP )
     
 END SUBROUTINE GetNumberOfInternalNodes_Hexmesh
!
!
!
 SUBROUTINE SetPositions_Hexmesh( myHexMesh, iEl, x, y )
 ! S/R SetPositions
 !
 ! =============================================================================================== !
 ! DECLARATIONS
  IMPLICIT NONE
  CLASS(Hexmesh), INTENT(inout) :: myHexMesh
  INTEGER, INTENT(in)            :: iEl
  REAL(prec), INTENT(in)         :: x(0:myHexMesh % elements(iEl) % nS, 0:myHexMesh % elements(iEl) % nP)
  REAL(prec), INTENT(in)         :: y(0:myHexMesh % elements(iEl) % nS, 0:myHexMesh % elements(iEl) % nP)
  
     CALL myHexMesh % elements(iEl) % SetPositions( x, y )
     
 END SUBROUTINE SetPositions_Hexmesh
!
!
!
 SUBROUTINE GetPositions_Hexmesh( myHexMesh, iEl, x, y )
 ! S/R GetPositions
 !
 ! =============================================================================================== !
 ! DECLARATIONS
  IMPLICIT NONE
  CLASS(Hexmesh), INTENT(in) :: myHexMesh
  INTEGER, INTENT(in)         :: iEl
  REAL(prec), INTENT(out)     :: x(0:myHexMesh % elements(iEl) % nS, 0:myHexMesh % elements(iEl) % nP)
  REAL(prec), INTENT(out)     :: y(0:myHexMesh % elements(iEl) % nS, 0:myHexMesh % elements(iEl) % nP)
  
     CALL myHexMesh % elements(iEl) % GetPositions( x, y )
     
 END SUBROUTINE GetPositions_Hexmesh
!
!
!
 SUBROUTINE SetPositionAtNode_Hexmesh( myHexMesh, iEl, x, y, i, j )
 ! S/R SetPositionAtNode
 !
 ! =============================================================================================== !
 ! DECLARATIONS
  IMPLICIT NONE
  CLASS(Hexmesh), INTENT(inout) :: myHexMesh
  INTEGER, INTENT(in)            :: iEl
  REAL(prec), INTENT(in)         :: x, y
  INTEGER, INTENT(in)            :: i, j
  
     CALL myHexMesh % elements(iEl) % SetPositionAtNode( x, y, i, j )
     
 END SUBROUTINE SetPositionAtNode_Hexmesh
!
!
!
 SUBROUTINE GetPositionAtNode_Hexmesh( myHexMesh, iEl, x, y, i, j )
 ! S/R GetPositionAtNode
 !
 ! =============================================================================================== !
 ! DECLARATIONS
  IMPLICIT NONE
  CLASS(Hexmesh), INTENT(in) :: myHexMesh
  INTEGER, INTENT(in)         :: iEl
  REAL(prec), INTENT(out)     :: x, y
  INTEGER, INTENT(in)         :: i, j
  
     CALL myHexMesh % elements(iEl) % GetPositionAtNode( x, y, i, j )
     
 END SUBROUTINE GetPositionAtNode_Hexmesh
!
!
!
 SUBROUTINE SetJacobian_Hexmesh( myHexMesh, iEl, J )
 ! S/R SetJacobian
 !
 ! =============================================================================================== !
 ! DECLARATIONS
  IMPLICIT NONE
  CLASS(Hexmesh), INTENT(inout) :: myHexMesh
  INTEGER, INTENT(in)            :: iEl
  REAL(prec), INTENT(in)         :: J(0:myHexMesh % elements(iEl) % nS, 0:myHexMesh % elements(iEl) % nP)
  
     CALL myHexMesh % elements(iEl) % SetJacobian( J )
     
 END SUBROUTINE SetJacobian_Hexmesh
!
!
!
 SUBROUTINE GetJacobian_Hexmesh( myHexMesh, iEl, J )
 ! S/R GetJacobian
 !
 ! =============================================================================================== !
 ! DECLARATIONS
  IMPLICIT NONE
  CLASS(Hexmesh), INTENT(in) :: myHexMesh
  INTEGER, INTENT(in)         :: iEl
  REAL(prec), INTENT(out)     :: J(0:myHexMesh % elements(iEl) % nS, 0:myHexMesh % elements(iEl) % nP)
  
     CALL myHexMesh % elements(iEl) % GetJacobian( J )
     
 END SUBROUTINE GetJacobian_Hexmesh
!
!
!
 SUBROUTINE SetJacobianAtNode_Hexmesh( myHexMesh, iEl, J, iS, iP )
 ! S/R SetJacobianAtNode
 !
 ! =============================================================================================== !
 ! DECLARATIONS
  IMPLICIT NONE
  CLASS(Hexmesh), INTENT(inout) :: myHexMesh
  INTEGER, INTENT(in)            :: iEl
  REAL(prec), INTENT(in)         :: J
  INTEGER, INTENT(in)            :: iS, iP
  
     CALL myHexMesh % elements(iEl) % SetJacobianAtNode( J, iS, iP )
     
 END SUBROUTINE SetJacobianAtNode_Hexmesh
!
!
!
 SUBROUTINE GetJacobianAtNode_Hexmesh( myHexMesh, iEl, J, iS, iP )
 ! S/R GetJacobianAtNode
 !
 ! =============================================================================================== !
 ! DECLARATIONS
  IMPLICIT NONE
  CLASS(Hexmesh), INTENT(in) :: myHexMesh
  INTEGER, INTENT(in)         :: iEl
  REAL(prec), INTENT(out)     :: J
  INTEGER, INTENT(in)         :: iS, iP
  
     CALL myHexMesh % elements(iEl) % GetJacobianAtNode( J, iS, iP )
     
 END SUBROUTINE GetJacobianAtNode_Hexmesh
!
!
!
 SUBROUTINE SetCovariantMetrics_Hexmesh( myHexMesh, iEl, dxds, dxdp, dyds, dydp )
 ! S/R SetCovariantMetrics
 !
 ! =============================================================================================== !
 ! DECLARATIONS
  IMPLICIT NONE
  CLASS(Hexmesh), INTENT(inout) :: myHexMesh
  INTEGER, INTENT(in)            :: iEl
  REAL(prec), INTENT(in)         :: dxds(0:myHexMesh % elements(iEl) % nS, 0:myHexMesh % elements(iEl) % nP)
  REAL(prec), INTENT(in)         :: dxdp(0:myHexMesh % elements(iEl) % nS, 0:myHexMesh % elements(iEl) % nP)
  REAL(prec), INTENT(in)         :: dyds(0:myHexMesh % elements(iEl) % nS, 0:myHexMesh % elements(iEl) % nP)
  REAL(prec), INTENT(in)         :: dydp(0:myHexMesh % elements(iEl) % nS, 0:myHexMesh % elements(iEl) % nP)
  
     CALL myHexMesh % elements(iEl) % SetCovariantMetrics( dxds, dxdp, dyds, dydp )
     
 END SUBROUTINE SetCovariantMetrics_Hexmesh
!
!
!
 SUBROUTINE GetCovariantMetrics_Hexmesh( myHexMesh, iEl, dxds, dxdp, dyds, dydp )
 ! S/R GetCovariantMetrics
 !
 ! =============================================================================================== !
 ! DECLARATIONS
  IMPLICIT NONE
  CLASS(Hexmesh), INTENT(in) :: myHexMesh
  INTEGER, INTENT(in)         :: iEl
  REAL(prec), INTENT(out)     :: dxds(0:myHexMesh % elements(iEl) % nS, 0:myHexMesh % elements(iEl) % nP)
  REAL(prec), INTENT(out)     :: dxdp(0:myHexMesh % elements(iEl) % nS, 0:myHexMesh % elements(iEl) % nP)
  REAL(prec), INTENT(out)     :: dyds(0:myHexMesh % elements(iEl) % nS, 0:myHexMesh % elements(iEl) % nP)
  REAL(prec), INTENT(out)     :: dydp(0:myHexMesh % elements(iEl) % nS, 0:myHexMesh % elements(iEl) % nP)
  
     CALL myHexMesh % elements(iEl) % GetCovariantMetrics( dxds, dxdp, dyds, dydp )
     
 END SUBROUTINE GetCovariantMetrics_Hexmesh
!
!
!
 SUBROUTINE SetCovariantMetricsAtNode_Hexmesh( myHexMesh, iEl, dxds, dxdp, dyds, dydp, iS, iP )
 ! S/R SetCovariantMetricsAtNode
 !
 ! =============================================================================================== !
 ! DECLARATIONS
  IMPLICIT NONE
  CLASS(Hexmesh), INTENT(inout) :: myHexMesh
  INTEGER, INTENT(in)            :: iEl
  REAL(prec), INTENT(in)         :: dxds, dxdp, dyds, dydp
  INTEGER, INTENT(in)            :: iS, iP
  
     CALL myHexMesh % elements(iEl) % SetCovariantMetricsAtNode( dxds, dxdp, dyds, dydp, iS, iP )
     
 END SUBROUTINE SetCovariantMetricsAtNode_Hexmesh
!
!
!
 SUBROUTINE GetCovariantMetricsAtNode_Hexmesh( myHexMesh, iEl, dxds, dxdp, dyds, dydp, iS, iP )
 ! S/R GetCovariantMetricsAtNode
 !
 ! =============================================================================================== !
 ! DECLARATIONS
  IMPLICIT NONE
  CLASS(Hexmesh), INTENT(in) :: myHexMesh
  INTEGER, INTENT(in)         :: iEl
  REAL(prec), INTENT(out)     :: dxds, dxdp, dyds, dydp
  INTEGER, INTENT(in)         :: iS, iP
  
     CALL myHexMesh % elements(iEl) % GetCovariantMetricsAtNode( dxds, dxdp, dyds, dydp, iS, iP )
     
 END SUBROUTINE GetCovariantMetricsAtNode_Hexmesh
!
!
!
 SUBROUTINE SetBoundaryLocation_Hexmesh( myHexMesh, iEl, x, y, iBound )
 ! S/R SetBoundaryLocation
 !
 ! =============================================================================================== !
 ! DECLARATIONS
  IMPLICIT NONE
  CLASS(Hexmesh), INTENT(inout) :: myHexMesh
  INTEGER, INTENT(in)            :: iEl
  INTEGER, INTENT(in)            :: iBound
  REAL(prec), INTENT(in)         :: x(0:myHexMesh % elements(iEl) % nMax)
  REAL(prec), INTENT(in)         :: y(0:myHexMesh % elements(iEl) % nMax)
  
     CALL myHexMesh % elements(iEl) % SetBoundaryLocation( x, y, iBound )
     
 END SUBROUTINE SetBoundaryLocation_Hexmesh
!
!
!
 SUBROUTINE GetBoundaryLocation_Hexmesh( myHexMesh, iEl, x, y, iBound )
 ! S/R GetBoundaryLocation
 !
 ! =============================================================================================== !
 ! DECLARATIONS
  IMPLICIT NONE
  CLASS(Hexmesh), INTENT(in) :: myHexMesh
  INTEGER, INTENT(in)         :: iEl
  INTEGER, INTENT(in)         :: iBound
  REAL(prec), INTENT(out)     :: x(0:myHexMesh % elements(iEl) % nMax)
  REAL(prec), INTENT(out)     :: y(0:myHexMesh % elements(iEl) % nMax)
  
     CALL myHexMesh % elements(iEl) % GetBoundaryLocation( x, y, iBound )
     
 END SUBROUTINE GetBoundaryLocation_Hexmesh
!
!
!
 SUBROUTINE SetBoundaryLocationAtNode_Hexmesh( myHexMesh, iEl, x, y, i, iBound )
 ! S/R SetBoundaryLocationAtNode
 !
 ! =============================================================================================== !
 ! DECLARATIONS
  IMPLICIT NONE
  CLASS(Hexmesh), INTENT(inout) :: myHexMesh
  INTEGER, INTENT(in)            :: iEl
  REAL(prec), INTENT(in)         :: x, y
  INTEGER, INTENT(in)            :: i, iBound
  
     CALL myHexMesh % elements(iEl) % SetBoundaryLocationAtNode( x, y, i, iBound )
     
 END SUBROUTINE SetBoundaryLocationAtNode_Hexmesh
!
!
!
 SUBROUTINE GetBoundaryLocationAtNode_Hexmesh( myHexMesh, iEl, x, y, i, iBound )
 ! S/R GetBoundaryLocationAtNode
 !
 ! =============================================================================================== !
 ! DECLARATIONS
  IMPLICIT NONE
  CLASS(Hexmesh), INTENT(in) :: myHexMesh
  INTEGER, INTENT(in)         :: iEl
  REAL(prec), INTENT(out)     :: x, y
  INTEGER, INTENT(in)         :: i, iBound
  
     CALL myHexMesh % elements(iEl) % GetBoundaryLocationAtNode( x, y, i, iBound )
     
 END SUBROUTINE GetBoundaryLocationAtNode_Hexmesh
!
!
!
 SUBROUTINE SetBoundaryNormalAtNode_Hexmesh( myHexMesh, iEl, dir, i, iBound )
 ! S/R SetBoundaryNormalAtNode
 !
 ! =============================================================================================== !
 ! DECLARATIONS
  IMPLICIT NONE
  CLASS(Hexmesh), INTENT(inout) :: myHexMesh
  INTEGER, INTENT(in)            :: iEl
  REAL(prec), INTENT(in)         :: dir(1:2)
  INTEGER, INTENT(in)            :: i, iBound
  
     CALL myHexMesh % elements(iEl) % SetBoundaryNormalAtNode( dir, i, iBound )
     
 END SUBROUTINE SetBoundaryNormalAtNode_Hexmesh
!
!
!
 SUBROUTINE GetBoundaryNormalAtNode_Hexmesh( myHexMesh, iEl, dir, length, i, iBound )
 ! S/R GetBoundaryNormalAtNode
 !
 ! =============================================================================================== !
 ! DECLARATIONS
  IMPLICIT NONE
  CLASS(Hexmesh), INTENT(in) :: myHexMesh
  INTEGER, INTENT(in)         :: iEl
  REAL(prec), INTENT(out)     :: dir(1:2), length
  INTEGER, INTENT(in)         :: i, iBound
  
     CALL myHexMesh % elements(iEl) % GetBoundaryNormalAtNode( dir, length, i, iBound )
     
 END SUBROUTINE GetBoundaryNormalAtNode_Hexmesh
!
! ------------------------------- Node Wrapper Routines ------------------------------------------ !
!
 SUBROUTINE SetNodeData_HexMesh( myHexMesh, iNode, x, y, nodeType )
 ! S/R SetNodeData
 !  
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS( HexMesh ), INTENT(inout) :: myHexMesh
   INTEGER, INTENT(in)              :: iNode
   REAL(prec), INTENT(in)           :: x, y
   INTEGER, INTENT(in)              :: nodeType

      CALL myHexMesh % nodes(iNode) % SetData( x, y, nodeType )

 END SUBROUTINE SetNodeData_HexMesh
!
!
!
 SUBROUTINE GetNodeData_HexMesh( myHexMesh, iNode, x, y, nodeType )
 ! S/R GetNodeData
 !  
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS( HexMesh ), INTENT(in) :: myHexMesh
   INTEGER, INTENT(in)           :: iNode
   REAL(prec), INTENT(out)       :: x, y
   INTEGER, INTENT(out)          :: nodeType

      CALL myHexMesh % nodes(iNode) % GetData( x, y, nodeType )

 END SUBROUTINE GetNodeData_HexMesh
!
!
!
 SUBROUTINE SetNodeKey_HexMesh( myHexMesh, iNode, key )
 ! S/R SetNodeKey
 !  
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS( HexMesh ), INTENT(inout) :: myHexMesh
   INTEGER, INTENT(in)              :: iNode
   INTEGER, INTENT(in)              :: key
   
      CALL myHexMesh % nodes(iNode) % SetKey( key )
   
 END SUBROUTINE SetNodeKey_HexMesh
!
!
!
 SUBROUTINE GetNodeKey_HexMesh( myHexMesh, iNode, key )
 ! S/R GetNodeKey
 !  
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS( HexMesh ), INTENT(in) :: myHexMesh
   INTEGER, INTENT(in)           :: iNode
   INTEGER, INTENT(out)          :: key
   
      CALL myHexMesh % nodes(iNode) % GetKey( key )
   
 END SUBROUTINE GetNodeKey_HexMesh
!
!
!
 SUBROUTINE SetNodeType_HexMesh( myHexMesh, iNode, nodeType )
 ! S/R SetNodeType_HexMesh
 !  
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS( HexMesh ), INTENT(inout) :: myHexMesh
   INTEGER, INTENT(in)              :: iNode
   INTEGER, INTENT(in)              :: nodeType
   
      CALL myHexMesh % nodes(iNode) % SetType( nodeType )
   
 END SUBROUTINE SetNodeType_HexMesh
!
!
!
 SUBROUTINE GetNodeType_HexMesh( myHexMesh, iNode, nodeType )
 ! S/R GetNodeType_HexMesh
 !  
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS( HexMesh ), INTENT(in) :: myHexMesh
   INTEGER, INTENT(in)           :: iNode
   INTEGER, INTENT(out)          :: nodeType
   
      CALL myHexMesh % nodes(iNode) % GetType( nodeType )
   
 END SUBROUTINE GetNodeType_HexMesh
!
!
!
 SUBROUTINE SetNodePosition_HexMesh( myHexMesh, iNode, x, y )
 ! S/R SetNodePosition_HexMesh
 !  
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS( HexMesh ), INTENT(inout) :: myHexMesh
   INTEGER, INTENT(in)              :: iNode
   REAL(prec), INTENT(in)           :: x, y
   
      CALL myHexMesh % nodes(iNode) % SetPosition( x, y )
   
 END SUBROUTINE SetNodePosition_HexMesh
!
!
!
 SUBROUTINE GetNodePosition_HexMesh( myHexMesh, iNode, x, y )
 ! S/R GetNodePosition_HexMesh
 !  
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS( HexMesh ), INTENT(in) :: myHexMesh
   INTEGER, INTENT(in)           :: iNode
   REAL(prec), INTENT(out)       :: x, y
   
      CALL myHexMesh % nodes(iNode) % GetPosition( x, y )
   
 END SUBROUTINE GetNodePosition_HexMesh
!
! --------------------------------- Edge Wrapper Routines ---------------------------------------- !
!
 SUBROUTINE SetEdgeData_HexMesh( myHexMesh, iEdge, nodeIDs, elementIDs, elementSides, key, start, inc )
 ! S/R SetEdgeData
 !  
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS( HexMesh ), INTENT(inout) :: myHexMesh
   INTEGER, INTENT(in)              :: iEdge
   INTEGER, INTENT(in)              :: nodeIDs(1:2), elementIDs(1:2), elementSides(1:2)
   INTEGER, INTENT(in)              :: key, start, inc

      CALL myHexMesh % edges(iEdge) % SetData( nodeIDs, elementIDs, elementSides, key, start, inc )

 END SUBROUTINE SetEdgeData_HexMesh
!
!
!
 SUBROUTINE GetEdgeData_HexMesh( myHexMesh, iEdge, nodeIDs, elementIDs, elementSides, key, start, inc )
 ! S/R GetEdgeData
 !  
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS( HexMesh ), INTENT(in) :: myHexMesh
   INTEGER, INTENT(in)           :: iEdge
   INTEGER, INTENT(out)          :: nodeIDs(1:2), elementIDs(1:2), elementSides(1:2)
   INTEGER, INTENT(out)          :: key, start, inc

      CALL myHexMesh % edges(iEdge) % GetData( nodeIDs, elementIDs, elementSides, key, start, inc )
      
 END SUBROUTINE GetEdgeData_HexMesh
!
!
!
SUBROUTINE SetEdgeKey_HexMesh( myHexMesh, iEdge, key )
 ! S/R SetEdgeKey
 !  
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS( HexMesh ), INTENT(inout) :: myHexMesh
   INTEGER, INTENT(in)              :: iEdge
   INTEGER, INTENT(in)              :: key

      CALL myHexMesh % edges(iEdge) % SetKey( key )    

 END SUBROUTINE SetEdgeKey_HexMesh
!
!
!
 SUBROUTINE GetEdgeKey_HexMesh( myHexMesh, iEdge, key )
 ! S/R GetEdgeKey
 !  
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS( HexMesh ), INTENT(in) :: myHexMesh
   INTEGER, INTENT(in)           :: iEdge
   INTEGER, INTENT(out)          :: key

      CALL myHexMesh % edges(iEdge) % GetKey( key )

 END SUBROUTINE GetEdgeKey_HexMesh
!
!
!
 SUBROUTINE SetEdgeNodeIDs_HexMesh( myHexMesh, iEdge, nodeIDs )
 ! S/R SetEdgeNodeIDs
 !  
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS( HexMesh ), INTENT(inout) :: myHexMesh
   INTEGER, INTENT(in)              :: iEdge
   INTEGER, INTENT(in)              :: nodeIDs(1:2)

      CALL myHexMesh % edges(iEdge) % SetNodeIDs( nodeIds )

 END SUBROUTINE SetEdgeNodeIDs_HexMesh
!
!
!
 SUBROUTINE GetEdgeNodeIDs_HexMesh( myHexMesh, iEdge, nodeIDs )
 ! S/R GetEdgeNodeIDs
 !  
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS( HexMesh ), INTENT(in) :: myHexMesh
   INTEGER, INTENT(in)           :: iEdge
   INTEGER, INTENT(out)          :: nodeIDs(1:2)

      CALL myHexMesh % edges(iEdge) % GetNodeIDs( nodeIDs )

 END SUBROUTINE GetEdgeNodeIDs_HexMesh
!
!
!
 SUBROUTINE SetEdgeElementIDs_HexMesh( myHexMesh, iEdge, elementIDs )
 ! S/R SetEdgeElementIDs
 !  
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS( HexMesh ), INTENT(inout) :: myHexMesh
   INTEGER, INTENT(in)              :: iEdge
   INTEGER, INTENT(in)              :: elementIDs(1:2)

      CALL myHexMesh % edges(iEdge) % SetElementIDs( elementIDs )

 END SUBROUTINE SetEdgeElementIDs_HexMesh
!
!
!
 SUBROUTINE GetEdgeElementIDs_HexMesh( myHexMesh, iEdge, elementIDs )
 ! S/R GetEdgeElementIDs
 !  
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS( HexMesh ), INTENT(in) :: myHexMesh
   INTEGER, INTENT(in)           :: iEdge
   INTEGER, INTENT(out)          :: elementIDs(1:2)

      CALL myHexMesh % edges(iEdge) % GetElementIDs( elementIDs )

 END SUBROUTINE GetEdgeElementIDs_HexMesh
!
!
!
 SUBROUTINE SetEdgePrimaryElementID_HexMesh( myHexMesh, iEdge, elementID )
 ! S/R SetEdgePrimaryElementID
 !  
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS( HexMesh ), INTENT(inout) :: myHexMesh
   INTEGER, INTENT(in)              :: iEdge
   INTEGER, INTENT(in)              :: elementID

      CALL myHexMesh % edges(iEdge) % SetPrimaryElementID( elementID )
      
 END SUBROUTINE SetEdgePrimaryElementID_HexMesh
!
!
!
 SUBROUTINE GetEdgePrimaryElementID_HexMesh( myHexMesh, iEdge, elementID )
 ! S/R GetEdgePrimaryElementID
 !  
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS( HexMesh ), INTENT(in) :: myHexMesh
   INTEGER, INTENT(in)           :: iEdge
   INTEGER, INTENT(out)          :: elementID

      CALL myHexMesh % edges(iEdge) % GetPrimaryElementID( elementID )

 END SUBROUTINE GetEdgePrimaryElementID_HexMesh
!
!
!
 SUBROUTINE SetEdgeSecondaryElementID_HexMesh( myHexMesh, iEdge, elementID )
 ! S/R SetEdgeSecondaryElementID
 !  
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS( HexMesh ), INTENT(inout) :: myHexMesh
   INTEGER, INTENT(in)              :: iEdge
   INTEGER, INTENT(in)              :: elementID

      CALL myHexMesh % edges(iEdge) % SetSecondaryElementID( elementID )
      
 END SUBROUTINE SetEdgeSecondaryElementID_HexMesh
!
!
!
 SUBROUTINE GetEdgeSecondaryElementID_HexMesh( myHexMesh, iEdge, elementID )
 ! S/R GetEdgeSecondaryElementID
 !  
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS( HexMesh ), INTENT(in) :: myHexMesh
   INTEGER, INTENT(in)           :: iEdge
   INTEGER, INTENT(out)          :: elementID

      CALL myHexMesh % edges(iEdge) % GetSecondaryElementID( elementID )

 END SUBROUTINE GetEdgeSecondaryElementID_HexMesh
!
!
!
 SUBROUTINE SetEdgeElementSides_HexMesh( myHexMesh, iEdge, elementSides )
 ! S/R SetEdgeElementSides
 !  
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS( HexMesh ), INTENT(inout) :: myHexMesh
   INTEGER, INTENT(in)              :: iEdge
   INTEGER, INTENT(in)              :: elementSides(1:2)

      CALL myHexMesh % edges(iEdge) % SetElementSides( elementSides )

 END SUBROUTINE SetEdgeElementSides_HexMesh
!
!
!
 SUBROUTINE GetEdgeElementSides_HexMesh( myHexMesh, iEdge, elementSides )
 ! S/R GetEdgeElementSides
 !  
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS( HexMesh ), INTENT(in) :: myHexMesh
   INTEGER, INTENT(in)           :: iEdge
   INTEGER, INTENT(out)          :: elementSides(1:2)

      CALL myHexMesh % edges(iEdge) % GetElementSides( elementSides )

 END SUBROUTINE GetEdgeElementSides_HexMesh
!
!
!
 SUBROUTINE SetEdgePrimaryElementSide_HexMesh( myHexMesh, iEdge, elementSide )
 ! S/R SetEdgePrimaryElementSide
 !  
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS( HexMesh ), INTENT(inout) :: myHexMesh
   INTEGER, INTENT(in)              :: iEdge
   INTEGER, INTENT(in)              :: elementSide

      CALL myHexMesh % edges(iEdge) % SetPrimaryElementSide( elementSide )
      
 END SUBROUTINE SetEdgePrimaryElementSide_HexMesh
!
!
!
 SUBROUTINE GetEdgePrimaryElementSide_HexMesh( myHexMesh, iEdge, elementSide )
 ! S/R GetEdgePrimaryElementSide
 !  
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS( HexMesh ), INTENT(in) :: myHexMesh
   INTEGER, INTENT(in)           :: iEdge
   INTEGER, INTENT(out)          :: elementSide

      CALL myHexMesh % edges(iEdge) % GetPrimaryElementSide( elementSide )

 END SUBROUTINE GetEdgePrimaryElementSide_HexMesh
!
!
!
 SUBROUTINE SetEdgeSecondaryElementSide_HexMesh( myHexMesh, iEdge, elementSide )
 ! S/R SetEdgeSecondaryElementSide
 !  
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS( HexMesh ), INTENT(inout) :: myHexMesh
   INTEGER, INTENT(in)              :: iEdge
   INTEGER, INTENT(in)              :: elementSide

      CALL myHexMesh % edges(iEdge) % SetSecondaryElementSide( elementSide )
      
 END SUBROUTINE SetEdgeSecondaryElementSide_HexMesh
!
!
!
 SUBROUTINE GetEdgeSecondaryElementSide_HexMesh( myHexMesh, iEdge, elementSide )
 ! S/R GetEdgeSecondaryElementSide
 !  
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS( HexMesh ), INTENT(in) :: myHexMesh
   INTEGER, INTENT(in)           :: iEdge
   INTEGER, INTENT(out)          :: elementSide

      CALL myHexMesh % edges(iEdge) % GetSecondaryElementSide( elementSide )

 END SUBROUTINE GetEdgeSecondaryElementSide_HexMesh
!
!
!
 SUBROUTINE SetEdgeStart_HexMesh( myHexMesh, iEdge, start )
 ! S/R SetEdgeStart
 !  
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS( HexMesh ), INTENT(inout) :: myHexMesh
   INTEGER, INTENT(in)              :: iEdge
   INTEGER, INTENT(in)              :: start
    
      CALL myHexMesh % edges(iEdge) % SetStart( start )
      
 END SUBROUTINE SetEdgeStart_HexMesh
!
!
!
 SUBROUTINE GetEdgeStart_HexMesh( myHexMesh, iEdge, start )
 ! S/R GetEdgeStart
 !  
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS( HexMesh ), INTENT(in) :: myHexMesh
   INTEGER, INTENT(in)           :: iEdge
   INTEGER, INTENT(out)          :: start

      CALL myHexMesh % edges(iEdge) % GetStart( start )

 END SUBROUTINE GetEdgeStart_HexMesh
!
!
!
 SUBROUTINE SetEdgeIncrement_HexMesh( myHexMesh, iEdge, inc )
 ! S/R SetEdgeIncrement
 !  
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS( HexMesh ), INTENT(inout) :: myHexMesh
   INTEGER, INTENT(in)              :: iEdge
   INTEGER, INTENT(in)              :: inc

      CALL myHexMesh % edges(iEdge) % SetIncrement( inc )

 END SUBROUTINE SetEdgeIncrement_HexMesh
!
!
!
 SUBROUTINE GetEdgeIncrement_HexMesh( myHexMesh, iEdge, inc )
 ! S/R GetEdgeIncrement
 !  
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS( HexMesh ), INTENT(in) :: myHexMesh
   INTEGER, INTENT(in)           :: iEdge
   INTEGER, INTENT(out)          :: inc

      CALL myHexMesh % edges(iEdge) % GetIncrement( inc )

 END SUBROUTINE GetEdgeIncrement_HexMesh
!
!
!
!==================================================================================================!
!--------------------------------- Type Specific Routines -----------------------------------------!
!==================================================================================================!
!
!
SUBROUTINE ConstructEdges_HexMesh( myHexMesh )
 ! S/R ConstructEdges
 !
 !    Takes in the mesh of quadrilaterals which has not filled in the 
 !    edge information, and finds all of the unique edges in the mesh.
 !    The space for the edges is REALlocated with the correct number of edges
 ! 
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS( HexMesh ), INTENT(inout) :: myHexMesh
   ! LOCAL
   TYPE( HashTable ) :: edgeTable
   INTEGER :: nEls, nNodes, iEl, nEdges, k  
   INTEGER :: l1, l2, startID, endID, key1, key2
   INTEGER :: e1, e2, s1, s2, edgeID, n1, nID

      CALL myHexMesh % GetNumberOfNodes( nNodes )   
      CALL myHexMesh % GetNumberOfElements( nEls )
      nEdges = 0

      ! First, just count the number of edges
      CALL edgeTable % Build( nNodes )

      do iEl = 1, nEls ! Loop over the elements in the mesh

         do k = 1, 4 ! Loop over the sides of each element

            l1 = myHexMesh % faceMap(1,k) ! starting local node for this edge
            l2 = myHexMesh % faceMap(2,k) ! ending local node for this edge
            
            CALL myHexMesh % GetElementNodeID( iEl, l1, startID )
            CALL myHexMesh % GetElementNodeID( iEl, l2, endID )
            
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

      DEALLOCATE( myHexMesh % edges )

      ALLOCATE( myHexMesh % edges( 1:nEdges ) )

      nEdges = 0 ! restart the edge counting

      do iEl = 1, nEls ! Loop over the elements in the mesh

         do k = 1, 4 ! Loop over the sides of each element

            l1 = myHexMesh % faceMap(1,k) ! starting local node for this edge
            l2 = myHexMesh % faceMap(2,k) ! ending local node for this edge

            CALL myHexMesh % GetElementNodeID( iEl, l1, startID )
            CALL myHexMesh % GetElementNodeID( iEl, l2, endID )

            key1 = min( startID, endID )
            key2 = max( startID, endID )

            IF( edgeTable % ContainsKeys( key1, key2 )  )then ! this edge already exists
               
               !Get the edgeID
               CALL edgeTable % GetDataForKeys( edgeID, key1, key2 )
               ! Find the primary element and the starting node for this element's edge
               ! This is compared with the secondary element's starting node to infer
               ! the relative orientation of the two elements.
               CALL myHexMesh % GetEdgePrimaryElementID( edgeID, e1 )
               CALL myHexMesh % GetEdgePrimaryElementSide( edgeID, s1 )
               
               l1 = myHexMesh % faceMap(1,s1)
               CALL myHexMesh % GetElementNodeID( e1, l1, n1 )

               ! Set the secondary element information
               CALL myHexMesh % SetEdgeSecondaryElementID( edgeID, iEl )
               
               !PRINT*, 'edgeID, primary, secondary:', edgeID, e1, iEl
               IF( startID == n1 ) then ! the elements are oriented the same direction
                  CALL myHexMesh % SetEdgeSecondaryElementSide( edgeID, k )
                 

               ELSE ! the elements are oriented in the opposite direction

                  ! For these edges, we mark the side ID as negative
                  CALL myHexMesh % SetEdgeSecondaryElementSide( edgeID, -k )
                 
               ENDIF

            ELSE ! this is a new edge

               ! Add the edge to the list
               nEdges = nEdges + 1
               
               edgeID = nEdges
               CALL myHexMesh % edges(edgeID) % Build()
               CALL myHexMesh % SetEdgePrimaryElementID( edgeID, iEl )
               CALL myHexMesh % SetEdgePrimaryElementSide( edgeID, k )
               CALL myHexMesh % SetEdgeNodeIDs( edgeID, (/startID, endID /) )

               ! Default the secondary information
               CALL myHexMesh % SetEdgeSecondaryElementID( edgeID, BoundaryFlagDefault )
               
               CALL edgeTable % AddDataForKeys( edgeID, key1, key2 )
               
            ENDIF
         enddo ! k, Loop over the sides of each element
        
      enddo ! iEl, Loop over the elements in the mesh

      CALL edgeTable % Trash( )
      myHexMesh % nEdges = nEdges

     ! do edgeID = 1, nEdges
      
     !    CALL myHexMesh % GetEdgePrimaryElementID( edgeID, e1 )
     !    CALL myHexMesh % GetEdgePrimaryElementSide( edgeID, s1 )
     !    CALL myHexMesh % GetEdgeSecondaryElementID( edgeID, e2 )
     !    CALL myHexMesh % GetEdgeSecondaryElementSide( edgeID, s2 )
         
     !    PRINT*, edgeID, e1, s1, e2, s2
     ! enddo
 END SUBROUTINE ConstructEdges_HexMesh
!
!
!
 SUBROUTINE GetNodeToElementConnectivity_HexMesh( myHexMesh )
 ! S/R GetNodeToElementConnectivity
 ! 
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS( HexMesh ), INTENT(inout) :: myHexMesh
   ! LOCAL
   INTEGER :: nEls, nNodes, iEl, nEdges, k  
   INTEGER :: l1, l2, startID, endID, key1, key2
   INTEGER :: e1, e2, s1, s2, edgeID, n1, nID

      CALL myHexMesh % GetNumberOfNodes( nNodes )   
      CALL myHexMesh % GetNumberOfElements( nEls )
      nEdges = 0

      DO iEl = 1, nEls ! Loop over the elements in the mesh

         DO k = 1, 4 ! Loop over the nodes

            ! Add element to corner-node connectivity list
            CALL myHexMesh % GetElementNodeID( iEl, k, nID )
            CALL myHexMesh % nodes(nID) % nodeToElement % AddToList( iEl, k ) ! the key is  the local node ID, the data is the element ID

         ENDDO ! k, Loop over the nodes
        
      ENDDO ! iEl, Loop over the elements in the mesh

     
 END SUBROUTINE GetNodeToElementConnectivity_HexMesh
!
!
!
 SUBROUTINE ConstructElementNeighbors_HexMesh( myHexMesh )
 ! S/R ConstructElementNeighbors
 ! 
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS( HexMesh ), INTENT(inout) :: myHexMesh
   ! LOCAL
   INTEGER :: e(1:2), s(1:2), nEdges, iEdge
  
      CALL myHexMesh % GetNumberOfEdges( nEdges )
 
      DO iEdge = 1, nEdges ! Loop over the edges

         CALL myHexMesh % GetEdgeElementIDs( iEdge, e )
         CALL myHexMesh % GetEdgeElementSides( iEdge, s )

         CALL myHexMesh % SetElementNeighbor( e(1), s(1), e(2) )
         CALL myHexMesh % SetElementNeighbor( e(2), ABS(s(2)), e(1) )
         
      ENDDO ! iEdge, Loop over the edges

     
 END SUBROUTINE ConstructElementNeighbors_HexMesh
!
!
!
 SUBROUTINE ScaleTheMesh_HexMesh( myHexMesh, interp, xScale, yScale  )
 ! S/R ScaleTheMesh
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS( HexMesh ), INTENT(inout) :: myHexMesh
   TYPE( Lagrange_2D ), INTENT(in)  :: interp
   REAL(prec), INTENT(in)           :: xScale, yScale
   ! LOCAL
   INTEGER :: nElems, nNodes, iEl, iNode
   
      CALL myHexMesh % GetNumberOfElements( nElems )
      CALL myHexMesh % GetNumberOfNodes( nNodes )

      DO iEl = 1, nElems
         CALL myHexMesh % elements(iEl) % ScaleGeometry( interp, xScale, yScale )
      ENDDO
      
      DO iNode = 1, nNodes
         CALL myHexMesh % nodes(iNode) % ScaleNodePosition( xScale, yScale )
      ENDDO

 END SUBROUTINE ScaleTheMesh_HexMesh
!
!
!
 SUBROUTINE LoadDefaultMesh_HexMesh( myHexMesh, interp, nXelem, nYelem  )
 ! S/R LoadDefaultMesh
 !  
 !     Builds a square mesh with 5x5 elements. The side lengths are 1 unit.
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS( HexMesh ), INTENT(inout)  :: myHexMesh
   TYPE( Lagrange_2D ), INTENT(in)   :: interp
   INTEGER, INTENT(in)               :: nXelem, nYelem
   ! LOCAL
   TYPE( Curve_2D ) :: elBoundCurves(1:4)

   REAL(prec) :: x, y, z, dxElem, dyElem
   REAL(prec) :: x1, x2, y1, y2, a, b, c
   REAL(prec), ALLOCATABLE :: xc(:), yc(:), s(:)

   INTEGER :: nNodes, nElems, nEdges, gPolyDeg
   INTEGER :: bFlags(1:4), nodes(1:4)
   INTEGER :: e1, e2, s1, s2, n1, n2, n(1:2), e(1:2), si(1:2)
   INTEGER :: iEdge, iNode, iEl, iSide, iX, iY
   INTEGER :: fUnit, iC, jC

   CHARACTER(20) :: ISMversion
   CHARACTER(40) :: edgeNames
   CHARACTER(40) :: thisEdge
      
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
  
      CALL myHexMesh % Build( nNodes, nElems, nEdges, interp % nS ) 
      
      
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
            CALL myHexMesh % SetNodePosition( iNode, x, y )
            
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
               n1 = nodes( myHexMesh % faceMap(1,iSide) )
               n2 = nodes( myHexMesh % faceMap(2,iSide) )

               CALL myHexMesh % GetNodePosition( n1, x1, y1 )
               CALL myHexMesh % GetNodePosition( n2, x2, y2 )
               
               DO iNode = 0, gPolyDeg
                  xc(iNode) = x1 + (x2-x1)*HALF*( s(iNode) + ONE )
                  yc(iNode) = y1 + (y2-y1)*HALF*( s(iNode) + ONE )
               ENDDO
               CALL elBoundCurves(iSide) % SetNodes( xc, yc ) 

            ENDDO
            CALL myHexMesh % elements(iEl) % Build( nodes, iEl, elBoundCurves, interp )
         ENDDO
      ENDDO ! iEl, cycle over the elements

      CALL myHexMesh % ConstructEdges( )
      nEdges = myHexMesh % nEdges
      PRINT*, 'nEdges    : ', nEdges
      
      ! Set the start and increment for the secondary element 
      DO iEdge = 1, nEdges
            CALL myHexMesh % GetEdgeSecondaryElementSide(iEdge, s2)
            IF(s2 < 0)THEN
               CALL myHexMesh % SetEdgeStart( iEdge, interp % nS-1 )
               CALL myHexMesh % SetEdgeIncrement( iEdge, -1)
            ELSE
               CALL myHexMesh % SetEdgeStart( iEdge, 1 )
               CALL myHexMesh % SetEdgeIncrement( iEdge, 1 )
            ENDIF
            
      ENDDO


      CALL myHexMesh % GetNodeToElementConnectivity( )

      ! Clear up memory
      DEALLOCATE( s, xc, yc )

      DO iSide = 1, 4
         CALL elBoundCurves(iSide) % Trash( )
      ENDDO
  
 END SUBROUTINE LoadDefaultMesh_HexMesh
!
!
!==================================================================================================!
!--------------------------------- Mesh File I/O Routines -----------------------------------------!
!==================================================================================================!
!
!
 SUBROUTINE ReadSpecMeshFile_HexMesh( myHexMesh, interp, filename )
 ! S/R ReadSpecMeshFile
 !  Desription:
 !    Reads in the ISM-v2 version of the SpecMesh output file and builds the mesh with it's geometry 
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS( HexMesh ), INTENT(out)  :: myHexMesh
   TYPE( Lagrange_2D ), INTENT(in) :: interp
   CHARACTER(*), INTENT(in)        :: filename
   ! LOCAL
   TYPE( Curve_2D ) :: elBoundCurves(1:4)

   REAL(prec) :: x, y, z
   REAL(prec) :: x1, x2, y1, y2
   REAL(prec), ALLOCATABLE :: xc(:), yc(:), s(:)

   INTEGER :: nNodes, nElems, nEdges, gPolyDeg
   INTEGER :: bFlags(1:4), nodes(1:4)
   INTEGER :: e1, e2, s1, s2, n1, n2, n(1:2), e(1:2), si(1:2)
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
  
      CALL myHexMesh % Build( nNodes, nElems, nEdges, interp % nS ) 
      
      ! ---- Read in the corner nodes ---- !

      DO iNode = 1, nNodes  ! Loop over the nodes in the file
         READ( fUnit, * ) x, y, z
         CALL myHexMesh % SetNodePosition( iNode, x, y )
      ENDDO


      ! ---- Read in the edge information ---- !

      DO iEdge = 1, nEdges

         READ( fUnit, * ) n, e, si 
         
         !CALL myHexMesh % SetEdgeNodeIDs( iEdge, n )
         !CALL myHexMesh % SetEdgeElementIDs( iEdge, e )
         !CALL myHexMesh % SetEdgeElementSides( iEdge, si )
       
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
               n1 = nodes( myHexMesh % faceMap(1,iSide) )
               n2 = nodes( myHexMesh % faceMap(2,iSide) )

               CALL myHexMesh % GetNodePosition( n1, x1, y1 )
               CALL myHexMesh % GetNodePosition( n2, x2, y2 )
               
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
         CALL myHexMesh % elements(iEl) % Build( nodes, iEl, elBoundCurves, interp )

         ! Read in and parse the edge names
         READ( fUnit, '(1x, A40)' )  edgeNames

         ! Parse the edge names into four edge names
         iSide = 1
         iC = 1
         jC = 1

         DO while( iSide <= 4 )

            IF( edgeNames(iC:iC) == ' ' )then ! we have reached a blank space
     
               n1 = nodes( myHexMesh % faceMap(1,iSide) )
               n2 = nodes( myHexMesh % faceMap(2,iSide) )
               
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
            CALL myHexMesh % ConstructEdges( )
      ! Set the secondary element on boundary edges to the boundary flag
      DO iEdge = 1, nEdges

         CALL myHexMesh % GetEdgePrimaryElementID( iEdge, e1 )
         CALL myHexMesh % GetEdgePrimaryElementSide( iEdge, s1 )
      
         IF( sideFlags(e1,s1) == DIRICHLET )then

            CALL myHexMesh % SetEdgeSecondaryElementID( iEdge, DIRICHLET )
            CALL myHexMesh % GetEdgeNodeIDs( iEdge, n )
            CALL myHexMesh % SetNodeType( n(1), DIRICHLET )
            CALL myHexMesh % SetNodeType( n(2), DIRICHLET )

         ELSEIF( sideFlags(e1,s1) == ROBIN )then

            CALL myHexMesh % SetEdgeSecondaryElementID( iEdge, ROBIN )
            CALL myHexMesh % GetEdgeNodeIDs( iEdge, n )
            CALL myHexMesh % SetNodeType( n(1), ROBIN )
            CALL myHexMesh % SetNodeType( n(2), ROBIN )

         ELSEIF( sideFlags(e1,s1) == ROBIN_FORCED )then

            CALL myHexMesh % SetEdgeSecondaryElementID( iEdge, ROBIN_FORCED )
            CALL myHexMesh % GetEdgeNodeIDs( iEdge, n )
            CALL myHexMesh % SetNodeType( n(1), ROBIN_FORCED )
            CALL myHexMesh % SetNodeType( n(2), ROBIN_FORCED )

         ELSEIF( sideFlags(e1,s1) == HOMOGENEOUS_NEUMANN )then

            CALL myHexMesh % SetEdgeSecondaryElementID( iEdge, HOMOGENEOUS_NEUMANN )
            CALL myHexMesh % GetEdgeNodeIDs( iEdge, n )
            CALL myHexMesh % SetNodeType( n(1), HOMOGENEOUS_NEUMANN )
            CALL myHexMesh % SetNodeType( n(2), HOMOGENEOUS_NEUMANN )

         ELSEIF( sideFlags(e1,s1) == NEUMANN_WALL )then

            CALL myHexMesh % SetEdgeSecondaryElementID( iEdge, NEUMANN_WALL )
            CALL myHexMesh % GetEdgeNodeIDs( iEdge, n )
            CALL myHexMesh % SetNodeType( n(1), NEUMANN_WALL )
            CALL myHexMesh % SetNodeType( n(2), NEUMANN_WALL )

         ELSEIF( sideFlags(e1,s1) == NEUMANN )then

            CALL myHexMesh % SetEdgeSecondaryElementID( iEdge, NEUMANN )
            CALL myHexMesh % GetEdgeNodeIDs( iEdge, n )
            CALL myHexMesh % SetNodeType( n(1), NEUMANN )
            CALL myHexMesh % SetNodeType( n(2), NEUMANN )

         ELSEIF( sideFlags(e1,s1) == NO_NORMAL_FLOW )then

            CALL myHexMesh % SetEdgeSecondaryElementID( iEdge, NO_NORMAL_FLOW )
            CALL myHexMesh % GetEdgeNodeIDs( iEdge, n )
            CALL myHexMesh % SetNodeType( n(1), NO_NORMAL_FLOW )
            CALL myHexMesh % SetNodeType( n(2), NO_NORMAL_FLOW )

         ELSEIF( sideFlags(e1,s1) == PRESCRIBED )then

            CALL myHexMesh % SetEdgeSecondaryElementID( iEdge, PRESCRIBED )
            CALL myHexMesh % GetEdgeNodeIDs( iEdge, n )
            CALL myHexMesh % SetNodeType( n(1), PRESCRIBED )
            CALL myHexMesh % SetNodeType( n(2), PRESCRIBED )

         ELSEIF( sideFlags(e1,s1) == RADIATION )then

            CALL myHexMesh % SetEdgeSecondaryElementID( iEdge, RADIATION )
            CALL myHexMesh % GetEdgeNodeIDs( iEdge, n )
            CALL myHexMesh % SetNodeType( n(1), RADIATION )
            CALL myHexMesh % SetNodeType( n(2), RADIATION )
         
         ENDIF
 
      ENDDO

      DO iEdge = 1, nEdges
            CALL myHexMesh % GetEdgeSecondaryElementSide(iEdge, s2)
            IF(s2 < 0)THEN
               CALL myHexMesh % SetEdgeStart( iEdge, interp % nS-1 )
               CALL myHexMesh % SetEdgeIncrement( iEdge, -1)
            ELSE
               CALL myHexMesh % SetEdgeStart( iEdge, 1 )
               CALL myHexMesh % SetEdgeIncrement( iEdge, 1 )
            ENDIF
            
      ENDDO
 ! If a node is at the corner of any boundary condition and a dirichlet boundary condition,
 ! we opt to keep it as a dirichlet boundary condition
    !  DO iEdge = 1, nEdges

    !     CALL myHexMesh % edges(iEdge) % GetPrimaryElementID( e1 )
    !     CALL myHexMesh % edges(iEdge) % GetPrimaryElementSide( s1 )

    !     IF( sideFlags(e1,s1) == DIRICHLET )then

    !        CALL myHexMesh % edges(iEdge) % SetSecondaryElementID( DIRICHLET )
            
    !        CALL myHexMesh % edges(iEdge) % GetNodeIDs( n )
    !        CALL myHexMesh % nodes(n(1)) % SetType( DIRICHLET )
    !        CALL myHexMesh % nodes(n(2)) % SetType( DIRICHLET )

    !     ENDIF
 
    !  ENDDO
      ! Get the node to element connectivity
      CALL myHexMesh % GetNodeToElementConnectivity( )

      ! Clear up memory
      DEALLOCATE( s, xc, yc, sideFlags )

      DO iSide = 1, 4
         CALL elBoundCurves(iSide) % Trash( )
      ENDDO
  
 END SUBROUTINE ReadSpecMeshFile_HexMesh
!
!
!
 SUBROUTINE WriteTecplot_Hexmesh( myHexMesh, filename )
 ! S/R WriteTecplot
 !
 ! =============================================================================================== !
 ! DECLARATIONS
  IMPLICIT NONE
  CLASS(HexMesh), INTENT(inout)     :: myHexMesh
  CHARACTER(*), INTENT(in), OPTIONAL :: filename  
  ! Local
  INTEGER :: iS, iP, nS, nP, iEl,fUnit, eID
  REAL(prec) :: x, y, J, dxds, dxdp, dyds, dydp
  CHARACTER(7) :: zoneID

    CALL myHexMesh % elements(1) % GetNumberOfNodes( nS, nP )

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


    DO iEl = 1, myHexMesh % nElems

       CALL myHexMesh % elements(iEl) % GetElementID( eID )
       WRITE(zoneID,'(I7.7)') eID
       WRITE(fUnit,*) 'ZONE T="el'//trim(zoneID)//'", I=',nS+1,', J=', nP+1,',F=POINT'

       DO iP = 0, nP
          DO iS = 0,nS
             CALL myHexMesh % GetPositionAtNode( iEl, x, y, iS, iP )
             CALL myHexMesh % GetJacobianAtNode( iEl, J, iS, iP )
             CALL myHexMesh % GetCovariantMetricsAtNode( iEl, dxds, dxdp, dyds, dydp, iS, iP )
             WRITE(fUnit,*)  x, y, J, dxds, dxdp, dyds, dydp
          ENDDO
      ENDDO

    ENDDO
    
    CLOSE(UNIT=fUnit)

    RETURN

 END SUBROUTINE WriteTecplot_Hexmesh
!
!
!
 

END MODULE HexMeshClass
