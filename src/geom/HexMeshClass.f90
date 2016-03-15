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
USE Lagrange_3D_Class
! src/geom/
USE MappFaceometryClass_3D
USE HexElementClass
USE FaceClass
USE NodeClass_3D



IMPLICIT NONE


! All properties are left as public in version 2.1 in order to allow the programmer the ability to
! have direct access to the attributes of the mesh class. This way, it can later be determined if
! additional accessor routines would be necessary.
 
    TYPE HexMesh 
       INTEGER                          :: nElems, nNodes, nFaces
       TYPE( HexElement ), ALLOCATABLE  :: elements(:)
       TYPE( Node ), ALLOCATABLE        :: nodes(:)  
       TYPE( Face ), ALLOCATABLE        :: faces(:)
       INTEGER                          :: cornerMap(1:3,1:nHexNodes) 
       INTEGER                          :: sideMap(1:nHexFaces) 
       INTEGER                          :: faceMap(1:nQuadNodes,1:nHexFaces) 
       INTEGER                          :: edgeFaceMap(1:2,1:nQuadEdges)

       CONTAINS

       PROCEDURE :: Build => Build_HexMesh
       PROCEDURE :: Trash => Trash_HexMesh
       
       PROCEDURE :: SetNumberOfElements => SetNumberOfElements_HexMesh
       PROCEDURE :: GetNumberOfElements => GetNumberOfElements_HexMesh
       PROCEDURE :: SetNumberOfNodes => SetNumberOfNodes_HexMesh
       PROCEDURE :: GetNumberOfNodes => GetNumberOfNodes_HexMesh
       PROCEDURE :: SetNumberOfFace => SetNumberOfFaces_HexMesh
       PROCEDURE :: GetNumberOfFaces => GetNumberOfFaces_HexMesh
       
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
       PROCEDURE :: SetElementBottomNeighbor  => SetElementBottomNeighbor_HexMesh
       PROCEDURE :: GetElementBottomNeighbor  => GetElementBottomNeighbor_HexMesh
       PROCEDURE :: SetElementTopNeighbor  => SetElementTopNeighbor_HexMesh
       PROCEDURE :: GetElementTopNeighbor  => GetElementTopNeighbor_HexMesh

       ! MappedGeometry_3D Wrapper Routines
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
       
       ! NodeClass_3D Wrapper Routines
       PROCEDURE :: SetNodeData => SetNodeData_HexMesh
       PROCEDURE :: GetNodeData => GetNodeData_HexMesh
       PROCEDURE :: SetNodeKey => SetNodeKey_HexMesh
       PROCEDURE :: GetNodeKey => GetNodeKey_HexMesh
       PROCEDURE :: SetNodeType => SetNodeType_HexMesh
       PROCEDURE :: GetNodeType => GetNodeType_HexMesh
       PROCEDURE :: SetNodePosition => SetNodePosition_HexMesh
       PROCEDURE :: GetNodePosition => GetNodePosition_HexMesh

       ! FaceClass Wrapper Routines 
       PROCEDURE :: SetFaceKey => SetFaceKey_HexMesh
       PROCEDURE :: GetFaceKey => GetFaceKey_HexMesh
       PROCEDURE :: SetFaceData => SetFaceData_HexMesh
       PROCEDURE :: GetFaceData => GetFaceData_HexMesh
       PROCEDURE :: SetFaceNodeIDs => SetFaceNodeIDs_HexMesh
       PROCEDURE :: GetFaceNodeIDs => GetFaceNodeIDs_HexMesh
       PROCEDURE :: SetFaceElementIDs => SetFaceElementIDs_HexMesh
       PROCEDURE :: GetFaceElementIDs => GetFaceElementIDs_HexMesh
       PROCEDURE :: SetFacePrimaryElementID => SetFacePrimaryElementID_HexMesh
       PROCEDURE :: GetFacePrimaryElementID => GetFacePrimaryElementID_HexMesh
       PROCEDURE :: SetFaceSecondaryElementID => SetFaceSecondaryElementID_HexMesh
       PROCEDURE :: GetFaceSecondaryElementID => GetFaceSecondaryElementID_HexMesh
       PROCEDURE :: SetFaceElementSides => SetFaceElementSides_HexMesh
       PROCEDURE :: GetFaceElementSides => GetFaceElementSides_HexMesh
       PROCEDURE :: SetFacePrimaryElementSide => SetFacePrimaryElementSide_HexMesh
       PROCEDURE :: GetFacePrimaryElementSide => GetFacePrimaryElementSide_HexMesh
       PROCEDURE :: SetFaceSecondaryElementSide => SetFaceSecondaryElementSide_HexMesh
       PROCEDURE :: GetFaceSecondaryElementSide => GetFaceSecondaryElementSide_HexMesh
       PROCEDURE :: SetFaceStart => SetFaceStart_HexMesh
       PROCEDURE :: GetFaceStart => GetFaceStart_HexMesh
       PROCEDURE :: SetFaceIncrement => SetFaceIncrement_HexMesh
       PROCEDURE :: GetFaceIncrement => GetFaceIncrement_HexMesh
       PROCEDURE :: SetSwapDimensions => SetSwapDimensions_HexMesh 
       PROCEDURE :: GetSwapDimensions => GetSwapDimensions_HexMesh
       
       PROCEDURE :: ConstructFaces => ConstructFaces_HexMesh
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
 SUBROUTINE Build_HexMesh( myHexMesh, nNodes, nElems, nFaces, nS )
 ! S/R Build
 ! 
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS(HexMesh), INTENT(out) :: myHexMesh
   INTEGER, INTENT(in)         :: nNodes, nElems, nFaces, nS
   !LOCAL
   INTEGER :: iNode

    ! A hexahedron element (hex-element for short) has six faces. Each face has geometry that 
    ! requires the use of two computational coordinates. The third computational coordinate is 
    ! fixed. The sideMap gives the value of the remaining computational coordinate for each face
    ! of the hex-element.
    ! 
    ! The face ordering is  1 - South, 2 - East, 3 - North, 4 - West, 5 - Bottom, 6 - Top 
    myHexmesh % sideMap(1:nHexFaces) = (/ 0, nS, nS, 0, 0, nS /)
      
    ! The eight corner nodes in an approximation that uses Gauss-Lobatto points, typically CG-type
    ! methods, have fixed computational coordinates. The corner node numbering  starts in the 
    ! southwest corner (in the computational grid) of the bottom face, proceeds counter clockwise
    ! around the face, and then repeats for the top face. This gives the local ID for the corner
    ! nodes as
    ! 
    ! Bottom, SouthWest = 1
    ! Bottom, SouthEast = 2
    ! Bottom, NorthEast = 3
    ! Bottom, NorthWest = 4
    ! Top, SouthWest = 5
    ! Top, SouthEast = 6
    ! Top, NorthEast = 7
    ! Top, NorthWest = 8
    !
    ! The computational coordinates for the corner nodes is given assuming a Gauss-Lobatto 
    ! computational mesh is used. Note that for a Gauss mesh, the corner nodes are not included.
     
    myHexmesh % cornerMap(1, 1:nHexNodes) = (/ 0, nS, nS,  0,  0, nS, nS,  0 /)
    myHexmesh % cornerMap(2, 1:nHexNodes) = (/ 0,  0, nS, nS,  0,  0, nS, nS /)
    myHexmesh % cornerMap(3, 1:nHexNodes) = (/ 0,  0,  0,  0, nS, nS, nS, nS /)
      
    ! Mesh construction usually begins with the specification of elements and the corner nodes, in
    ! addition to the element geometry. From the element-to-node connectivity, we need to construct
    ! the unique faces in the mesh and specify the abutting elements and their relative orientation.
    ! This procedure is aided by a "convenience array" that lists the local corner node IDs in the
    ! local counter-clockwise direction beginning in the local southwest corner of the face. When 
    ! two elements share a face, the global node IDs for each element can be found using this 
    ! convenience array (called "faceMap") and the relative orientation of the neighboring elements
    ! can be determined. The first index cycles over the nodes which make up the face in a 
    ! counterclockwise direction. The second index cycles over the faces in the element.
      
    myHexMesh % faceMap(1:nQuadNodes, south)  = (/ 1, 2, 6, 5 /)
    myHexMesh % faceMap(1:nQuadNodes, east)   = (/ 2, 3, 7, 6 /)
    myHexMesh % faceMap(1:nQuadNodes, north)  = (/ 4, 3, 7, 8 /)
    myHexMesh % faceMap(1:nQuadNodes, west)   = (/ 1, 4, 8, 5 /)
    myHexMesh % faceMap(1:nQuadNodes, bottom) = (/ 1, 2, 3, 4 /)
    myHexMesh % faceMap(1:nQuadNodes, top)    = (/ 5, 6, 7, 8 /)
      
    ! Each of the faces can be identified by their four corner nodes. The geometry of the faces
    ! is described using two computational coordinates between [-1,1]X[-1,1]. This 2-D computational
    ! grid has its own "southwest", "southeast", "northeast", and "northwest" identifications. The
    ! The corner nodes of the faces are labeled in the order mentioned in the previous sentence
    ! ( counterclockwise starting from southwest ). For quick referencing when producing tri-linear
    ! elements, a book-keeping array is useful for ensuring that we reference each edge of the 
    ! face in the order of increasing computational coordinate. This is identical to what is done
    ! in the QuadMeshClass.f90 for the "edgeMap". Here, it is called the "edgeFaceMap".
    ! The first index references the starting(1) or ending(2) node. The second index references
    ! the edge of the face with the first being the southern edge and increasing the second index
    ! proceeds counter-clockwise.
      
    myHexmesh % edgeFaceMap(1, 1:nQuadEdges) = (/ 1, 2, 4, 1 /)
    myHexmesh % edgeFaceMap(2, 1:nQuadEdges) = (/ 2, 3, 3, 4 /)

    ! The number of nodes, the number of elements, and the number of faces are stored in this data
    ! structure for convenience. In another implementation (planned for the next version), the 
    ! number of elements, nodes, and faces is dynamic; that implementation 
    ! requires the use of dynamic storage, e.g. a linked-list like structure for elements, edges,
    ! and nodes.
    CALL myHexMesh % SetNumberOfNodes( nNodes )
    CALL myHexMesh % SetNumberOfElements( nElems )
    CALL myHexMesh % SetNumberOfFaces( nFaces )

    ALLOCATE( myHexmesh % elements(1:nElems) )
    ALLOCATE( myHexmesh % nodes(1:nNodes) )
    ALLOCATE( myHexmesh % Faces(1:nFaces) )
     
      ! Default nodes to origin
    DO iNode = 1, myHexmesh % nNodes ! loop over the number of nodes
       CALL myHexmesh % nodes(iNode) % Build( ZERO, ZERO, ZERO )
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

    DEALLOCATE( myHexMesh % nodes, myHexMesh % elements, myHexMesh % Faces )
      

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
 SUBROUTINE SetNumberOfFaces_HexMesh( myHexMesh, nFaces )
 ! S/R SetNumberOfFaces
 ! 
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS(HexMesh), INTENT(inout) :: myHexMesh
   INTEGER, INTENT(in)            :: nFaces
   
    myHexMesh % nFaces = nFaces
      
 END SUBROUTINE SetNumberOfFaces_HexMesh
!
!
!
 SUBROUTINE GetNumberOfFaces_HexMesh( myHexMesh, nFaces )
 ! S/R GetNumberOfFaces
 ! 
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS(HexMesh), INTENT(in) :: myHexMesh
   INTEGER, INTENT(out)        :: nFaces
   
    nFaces = myHexMesh % nFaces
      
 END SUBROUTINE GetNumberOfFaces_HexMesh
!
! ------------------------------------------------------------------------------------------------ !
! ---------------------------- HexElement Wrapper Routines --------------------------------------- !
! ------------------------------------------------------------------------------------------------ !
!
SUBROUTINE SetElementNodeIDs_Hexmesh( myHexMesh, iEl, nodeIDs )
 ! S/R SetElementNodeIDs
 !
 ! =============================================================================================== !
   IMPLICIT NONE
   CLASS( Hexmesh ), INTENT(inout) :: myHexMesh
   INTEGER, INTENT(in)              :: iEl
   INTEGER, INTENT(in)              :: nodeIDs(1:nHexNodes)
   
    CALL myHexMesh % elements(iEl)  % SetNodeIDs( nodeIDs )

 END SUBROUTINE SetElementNodeIDs_Hexmesh
!
!
!
 FUNCTION GetElementNodeIDs_Hexmesh( myHexMesh, iEl ) RESULT( nodeIDs )
 ! S/R GetElementNodeIDs
 !
 ! =============================================================================================== !
   IMPLICIT NONE
   CLASS( Hexmesh ) :: myHexMesh
   INTEGER          :: iEl
   INTEGER          :: nodeIDs(1:nHexNodes)
   
    nodeIDs =  myHexMesh % elements(iEl)  % GetNodeIDs( )

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
 FUNCTION GetElementNodeID_Hexmesh( myHexMesh, iEl, localID ) RESULT( nodeID )
 ! S/R GetElementNodeID
 !
 ! =============================================================================================== !
   IMPLICIT NONE
   CLASS( Hexmesh ) :: myHexMesh
   INTEGER          :: iEl
   INTEGER          :: localID
   INTEGER          :: nodeID
   
    nodeID = myHexMesh % elements(iEl) % GetNodeID( localID )

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
 FUNCTION GetElementID_Hexmesh( myHexMesh, iEl ) RESULT( eID )
 ! S/R GetElementID
 !
 ! =============================================================================================== !
   IMPLICIT NONE
   CLASS( Hexmesh ) :: myHexMesh
   INTEGER          :: iEl
   INTEGER          :: eID
   
    eID = myHexMesh % elements(iEl) % GetElementID( )

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
 FUNCTION GetElementNeighbor_Hexmesh( myHexMesh, iEl, sID ) RESULT( eID )
 ! S/R GetElementNeighbor
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS( Hexmesh ) :: myHexMesh
   INTEGER          :: iEl
   INTEGER          :: sID
   INTEGER          :: eID
   
    eID = myHexMesh % elements(iEl) % GetNeighbor( sID )
      
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
 FUNCTION GetElementSouthernNeighbor_Hexmesh( myHexMesh, iEl ) RESULT( eID )
 ! S/R GetElementSouthernNeighbor
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS( Hexmesh ) :: myHexMesh
   INTEGER          :: iEl
   INTEGER          :: eID
   
    eID = myHexMesh % elements(iEl) % GetSouthernNeighbor( )
      
 END FUNCTION GetElementSouthernNeighbor_Hexmesh
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
 FUNCTION GetElementNorthernNeighbor_Hexmesh( myHexMesh, iEl ) RESULT( eID )
 ! S/R GetElementNorthernNeighbor
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS( Hexmesh ) :: myHexMesh
   INTEGER          :: iEl 
   INTEGER          :: eID
   
    eID = myHexMesh % elements(iEl) % GetNorthernNeighbor( )
      
 END FUNCTION GetElementNorthernNeighbor_Hexmesh
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
 FUNCTION GetElementEasternNeighbor_Hexmesh( myHexMesh, iEl ) RESULT( eID )
 ! S/R GetElementEasternNeighbor
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS( Hexmesh ) :: myHexMesh
   INTEGER          :: iEl
   INTEGER          :: eID
   
    eID = myHexMesh % elements(iEl) % GetEasternNeighbor(  )
      
 END FUNCTION GetElementEasternNeighbor_Hexmesh
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
   INTEGER, INTENT(in)             :: iEl
   INTEGER, INTENT(in)             :: eID
   
    CALL myHexMesh % elements(iEl) % SetWesternNeighbor( eID )
      
 END SUBROUTINE SetElementWesternNeighbor_Hexmesh
!
!
!
 FUNCTION GetElementWesternNeighbor_Hexmesh( myHexMesh, iEl ) RESULT ( eID )
 ! S/R GetElementWesternNeighbor
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS( Hexmesh ) :: myHexMesh
   INTEGER          :: iEl
   INTEGER          :: eID
   
    eID = myHexMesh % elements(iEl) % GetWesternNeighbor(  )
      
 END FUNCTION GetElementWesternNeighbor_Hexmesh
!
!
!
 SUBROUTINE SetElementBottomNeighbor_Hexmesh( myHexMesh, iEl, eID )
 ! S/R SetElementBottomNeighbor
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS( Hexmesh ), INTENT(inout) :: myHexMesh
   INTEGER, INTENT(in)             :: iEl
   INTEGER, INTENT(in)             :: eID
   
    CALL myHexMesh % elements(iEl) % SetbottomNeighbor( eID )
      
 END SUBROUTINE SetElementBottomNeighbor_Hexmesh
!
!
!
 FUNCTION GetElementBottomNeighbor_Hexmesh( myHexMesh, iEl ) RESULT ( eID )
 ! S/R GetElementBottomNeighbor
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS( Hexmesh ) :: myHexMesh
   INTEGER          :: iEl
   INTEGER          :: eID
   
    eID = myHexMesh % elements(iEl) % GetbottomNeighbor(  )
      
 END FUNCTION GetElementBottomNeighbor_Hexmesh
!
!
!
 SUBROUTINE SetElementTopNeighbor_Hexmesh( myHexMesh, iEl, eID )
 ! S/R SetElementTopNeighbor
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS( Hexmesh ), INTENT(inout) :: myHexMesh
   INTEGER, INTENT(in)             :: iEl
   INTEGER, INTENT(in)             :: eID
   
    CALL myHexMesh % elements(iEl) % SetTopNeighbor( eID )
      
 END SUBROUTINE SetElementTopNeighbor_Hexmesh
!
!
!
 FUNCTION GetElementTopNeighbor_Hexmesh( myHexMesh, iEl ) RESULT ( eID )
 ! S/R GetElementTopNeighbor
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS( Hexmesh ) :: myHexMesh
   INTEGER          :: iEl
   INTEGER          :: eID
   
    eID = myHexMesh % elements(iEl) % GetTopNeighbor(  )
      
 END FUNCTION GetElementTopNeighbor_Hexmesh
!
! ------------------------------------------------------------------------------------------------ !
! ---------------------------- MappedGeometryClass_3D Wrappers ----------------------------------- !
! ------------------------------------------------------------------------------------------------ !
!
SUBROUTINE SetNumberOfInternalNodes_Hexmesh( myHexMesh, iEl, nS, nP, nQ )
 ! S/R SetNumberOfInternalNodes
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS(Hexmesh), INTENT(inout) :: myHexMesh
   INTEGER, INTENT(in)           :: iEl
   INTEGER, INTENT(in)           :: nS, nP, nQ
  
    CALL myHexMesh % elements(iEl) % SetNumberOfNodes( nS, nP, nQ )
     
 END SUBROUTINE SetNumberOfInternalNodes_Hexmesh
!
!
!
 SUBROUTINE GetNumberOfInternalNodes_Hexmesh( myHexMesh, iEl, nS, nP, nQ )
 ! S/R GetNumberOfInternalNodes
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS(Hexmesh), INTENT(in) :: myHexMesh
   INTEGER, INTENT(in)        :: iEl
   INTEGER, INTENT(out)       :: nS, nP, nQ
  
    CALL myHexMesh % elements(iEl) % GetNumberOfNodes( nS, nP, nQ )
     
 END SUBROUTINE GetNumberOfInternalNodes_Hexmesh
!
!
!
 SUBROUTINE SetPositions_Hexmesh( myHexMesh, iEl, x, y, z )
 ! S/R SetPositions
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS(Hexmesh), INTENT(inout) :: myHexMesh
   INTEGER, INTENT(in)           :: iEl
   REAL(prec), INTENT(in)        :: x(0:myHexMesh % elements(iEl) % nS, &
                                      0:myHexMesh % elements(iEl) % nP, & 
                                      0:myHexMesh % elements(iEl) % nQ )
   REAL(prec), INTENT(in)        :: y(0:myHexMesh % elements(iEl) % nS, &
                                      0:myHexMesh % elements(iEl) % nP, & 
                                      0:myHexMesh % elements(iEl) % nQ )
   REAL(prec), INTENT(in)        :: z(0:myHexMesh % elements(iEl) % nS, &
                                      0:myHexMesh % elements(iEl) % nP, & 
                                      0:myHexMesh % elements(iEl) % nQ )
  
    CALL myHexMesh % elements(iEl) % SetPositions( x, y, z )
     
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
   INTEGER, INTENT(in)        :: iEl
   REAL(prec), INTENT(out)    :: x(0:myHexMesh % elements(iEl) % nS, &
                                   0:myHexMesh % elements(iEl) % nP, & 
                                   0:myHexMesh % elements(iEl) % nQ )
   REAL(prec), INTENT(out)    :: y(0:myHexMesh % elements(iEl) % nS, &
                                   0:myHexMesh % elements(iEl) % nP, & 
                                   0:myHexMesh % elements(iEl) % nQ )
   REAL(prec), INTENT(out)    :: z(0:myHexMesh % elements(iEl) % nS, &
                                   0:myHexMesh % elements(iEl) % nP, & 
                                   0:myHexMesh % elements(iEl) % nQ )
   
    CALL myHexMesh % elements(iEl) % GetPositions( x, y, z )
     
 END SUBROUTINE GetPositions_Hexmesh
!
!
!
 SUBROUTINE SetPositionAtNode_Hexmesh( myHexMesh, iEl, x, y, z, i, j, k )
 ! S/R SetPositionAtNode
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS(Hexmesh), INTENT(inout) :: myHexMesh
   INTEGER, INTENT(in)            :: iEl
   REAL(prec), INTENT(in)         :: x, y, z
   INTEGER, INTENT(in)            :: i, j, k
  
    CALL myHexMesh % elements(iEl) % SetPositionAtNode( x, y, z, i, j, k )
     
 END SUBROUTINE SetPositionAtNode_Hexmesh
!
!
!
 SUBROUTINE GetPositionAtNode_Hexmesh( myHexMesh, iEl, x, y, z, i, j, k )
 ! S/R GetPositionAtNode
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS(Hexmesh), INTENT(in) :: myHexMesh
   INTEGER, INTENT(in)         :: iEl
   REAL(prec), INTENT(out)     :: x, y, z
   INTEGER, INTENT(in)         :: i, j, k
  
    CALL myHexMesh % elements(iEl) % GetPositionAtNode( x, y, z, i, j, k )
     
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
   INTEGER, INTENT(in)           :: iEl
   REAL(prec), INTENT(in)        :: J(0:myHexMesh % elements(iEl) % nS, &
                                      0:myHexMesh % elements(iEl) % nP, &
                                      0:myHexMesh % elements(iEl) % nQ )
  
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
   INTEGER, INTENT(in)        :: iEl
   REAL(prec), INTENT(out)    :: J(0:myHexMesh % elements(iEl) % nS, &
                                   0:myHexMesh % elements(iEl) % nP, &
                                   0:myHexMesh % elements(iEl) % nQ )
  
    CALL myHexMesh % elements(iEl) % GetJacobian( J )
     
 END SUBROUTINE GetJacobian_Hexmesh
!
!
!
 SUBROUTINE SetJacobianAtNode_Hexmesh( myHexMesh, iEl, J, iS, iP, iQ )
 ! S/R SetJacobianAtNode
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS(Hexmesh), INTENT(inout) :: myHexMesh
   INTEGER, INTENT(in)           :: iEl
   REAL(prec), INTENT(in)        :: J
   INTEGER, INTENT(in)           :: iS, iP, iQ
  
    CALL myHexMesh % elements(iEl) % SetJacobianAtNode( J, iS, iP, iQ )
     
 END SUBROUTINE SetJacobianAtNode_Hexmesh
!
!
!
 SUBROUTINE GetJacobianAtNode_Hexmesh( myHexMesh, iEl, J, iS, iP, iQ )
 ! S/R GetJacobianAtNode
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS(Hexmesh), INTENT(in) :: myHexMesh
   INTEGER, INTENT(in)        :: iEl
   REAL(prec), INTENT(out)    :: J
   INTEGER, INTENT(in)        :: iS, iP, iQ
  
    CALL myHexMesh % elements(iEl) % GetJacobianAtNode( J, iS, iP, iQ )
     
 END SUBROUTINE GetJacobianAtNode_Hexmesh
!
!
!
 SUBROUTINE SetCovariantMetrics_Hexmesh( myHexMesh, iEl, dxds, dxdp, dxdq, &
                                                         dyds, dydp, dydq, &
                                                         dzds, dzdp, dzdq )
 ! S/R SetCovariantMetrics
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS(Hexmesh), INTENT(inout) :: myHexMesh
   INTEGER, INTENT(in)           :: iEl
   REAL(prec), INTENT(in)        :: dxds(0:myHexMesh % elements(iEl) % nS, &
                                         0:myHexMesh % elements(iEl) % nP, &
                                         0:myHexMesh % elements(iEl) % nQ )
   REAL(prec), INTENT(in)        :: dxdp(0:myHexMesh % elements(iEl) % nS, &
                                         0:myHexMesh % elements(iEl) % nP, &
                                         0:myHexMesh % elements(iEl) % nQ )
   REAL(prec), INTENT(in)        :: dxdq(0:myHexMesh % elements(iEl) % nS, &
                                         0:myHexMesh % elements(iEl) % nP, &
                                         0:myHexMesh % elements(iEl) % nQ )
   REAL(prec), INTENT(in)        :: dyds(0:myHexMesh % elements(iEl) % nS, &
                                         0:myHexMesh % elements(iEl) % nP, &
                                         0:myHexMesh % elements(iEl) % nQ )
   REAL(prec), INTENT(in)        :: dydp(0:myHexMesh % elements(iEl) % nS, &
                                         0:myHexMesh % elements(iEl) % nP, &
                                         0:myHexMesh % elements(iEl) % nQ )
   REAL(prec), INTENT(in)        :: dydq(0:myHexMesh % elements(iEl) % nS, &
                                         0:myHexMesh % elements(iEl) % nP, &
                                         0:myHexMesh % elements(iEl) % nQ )
   REAL(prec), INTENT(in)        :: dzds(0:myHexMesh % elements(iEl) % nS, &
                                         0:myHexMesh % elements(iEl) % nP, &
                                         0:myHexMesh % elements(iEl) % nQ )
   REAL(prec), INTENT(in)        :: dzdp(0:myHexMesh % elements(iEl) % nS, &
                                         0:myHexMesh % elements(iEl) % nP, &
                                         0:myHexMesh % elements(iEl) % nQ )
   REAL(prec), INTENT(in)        :: dzdq(0:myHexMesh % elements(iEl) % nS, &
                                         0:myHexMesh % elements(iEl) % nP, &
                                         0:myHexMesh % elements(iEl) % nQ ) 

    CALL myHexMesh % elements(iEl) % SetCovariantMetrics( dxds, dxdp, dxdq, &
                                                          dyds, dydp, dydq, &
                                                          dzds, dzdp, dzdq )
     
 END SUBROUTINE SetCovariantMetrics_Hexmesh
!
!
!
 SUBROUTINE GetCovariantMetrics_Hexmesh( myHexMesh, iEl, dxds, dxdp, dxdq, &
                                                         dyds, dydp, dydq, &
                                                         dzds, dzdp, dzdq )
 ! S/R GetCovariantMetrics
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS(Hexmesh), INTENT(in) :: myHexMesh
   INTEGER, INTENT(in)        :: iEl
   REAL(prec), INTENT(out)    :: dxds(0:myHexMesh % elements(iEl) % nS, &
                                      0:myHexMesh % elements(iEl) % nP, &
                                      0:myHexMesh % elements(iEl) % nQ )
   REAL(prec), INTENT(out)    :: dxdp(0:myHexMesh % elements(iEl) % nS, &
                                      0:myHexMesh % elements(iEl) % nP, &
                                      0:myHexMesh % elements(iEl) % nQ )
   REAL(prec), INTENT(out)    :: dxdq(0:myHexMesh % elements(iEl) % nS, &
                                      0:myHexMesh % elements(iEl) % nP, &
                                      0:myHexMesh % elements(iEl) % nQ )
   REAL(prec), INTENT(out)    :: dyds(0:myHexMesh % elements(iEl) % nS, &
                                      0:myHexMesh % elements(iEl) % nP, &
                                      0:myHexMesh % elements(iEl) % nQ )
   REAL(prec), INTENT(out)    :: dydp(0:myHexMesh % elements(iEl) % nS, &
                                      0:myHexMesh % elements(iEl) % nP, &
                                      0:myHexMesh % elements(iEl) % nQ )
   REAL(prec), INTENT(out)    :: dydq(0:myHexMesh % elements(iEl) % nS, &
                                      0:myHexMesh % elements(iEl) % nP, &
                                      0:myHexMesh % elements(iEl) % nQ )
   REAL(prec), INTENT(out)    :: dzds(0:myHexMesh % elements(iEl) % nS, &
                                      0:myHexMesh % elements(iEl) % nP, &
                                      0:myHexMesh % elements(iEl) % nQ )
   REAL(prec), INTENT(out)    :: dzdp(0:myHexMesh % elements(iEl) % nS, &
                                      0:myHexMesh % elements(iEl) % nP, &
                                      0:myHexMesh % elements(iEl) % nQ )
   REAL(prec), INTENT(out)    :: dzdq(0:myHexMesh % elements(iEl) % nS, &
                                      0:myHexMesh % elements(iEl) % nP, &
                                      0:myHexMesh % elements(iEl) % nQ )
  
    CALL myHexMesh % elements(iEl) % GetCovariantMetrics( dxds, dxdp, dxdq, &
                                                          dyds, dydp, dydq, &
                                                          dzds, dzdp, dzdq )
    
 END SUBROUTINE GetCovariantMetrics_Hexmesh
!
!
!
 SUBROUTINE SetCovariantMetricsAtNode_Hexmesh( myHexMesh, iEl, dxds, dxdp, dxdq, &
                                                               dyds, dydp, dydq, &
                                                               dzds, dzdp, dzdq, &
                                                               iS, iP, iQ )
 ! S/R SetCovariantMetricsAtNode
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS(Hexmesh), INTENT(inout) :: myHexMesh
   INTEGER, INTENT(in)           :: iEl
   REAL(prec), INTENT(in)        :: dxds, dxdp, dxdq, dyds, dydp, dydq, dzds, dzdp, dzdq
   INTEGER, INTENT(in)           :: iS, iP, iQ
   
    CALL myHexMesh % elements(iEl) % SetCovariantMetricsAtNode( dxds, dxdp, dxdq, &
                                                                dyds, dydp, dydq, &
                                                                dzds, dzdp, dzdq, &
                                                                iS, iP, iQ )
     
 END SUBROUTINE SetCovariantMetricsAtNode_Hexmesh
!
!
!
 SUBROUTINE GetCovariantMetricsAtNode_Hexmesh( myHexMesh, iEl, dxds, dxdp, dxdq, &
                                                               dyds, dydp, dydq, &
                                                               dzds, dzdp, dzdq, &
                                                               iS, iP, iQ )
 ! S/R GetCovariantMetricsAtNode
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS(Hexmesh), INTENT(in) :: myHexMesh
   INTEGER, INTENT(in)        :: iEl
   REAL(prec), INTENT(out)    :: dxds, dxdp, dxdq, dyds, dydp, dydq, dzds, dzdp, dzdq
   INTEGER, INTENT(in)        :: iS, iP, iQ
   
    CALL myHexMesh % elements(iEl) % GetCovariantMetricsAtNode( dxds, dxdp, dxdq, &
                                                                dyds, dydp, dydq, &
                                                                dzds, dzdp, dzdq, &
                                                                iS, iP, iQ )
     
 END SUBROUTINE GetCovariantMetricsAtNode_Hexmesh
!
!
!
 SUBROUTINE SetBoundaryLocation_Hexmesh( myHexMesh, iEl, x, y, z, iBound )
 ! S/R SetBoundaryLocation
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS(Hexmesh), INTENT(inout) :: myHexMesh
   INTEGER, INTENT(in)           :: iEl
   INTEGER, INTENT(in)           :: iBound
   REAL(prec), INTENT(in)        :: x(0:myHexMesh % elements(iEl) % nMax, &
                                      0:myHexMesh % elements(iEl) % nMax)
   REAL(prec), INTENT(in)        :: y(0:myHexMesh % elements(iEl) % nMax, &
                                      0:myHexMesh % elements(iEl) % nMax)
   REAL(prec), INTENT(in)        :: z(0:myHexMesh % elements(iEl) % nMax, &
                                      0:myHexMesh % elements(iEl) % nMax)  

    CALL myHexMesh % elements(iEl) % SetBoundaryLocation( x, y, z, iBound )
     
 END SUBROUTINE SetBoundaryLocation_Hexmesh
!
!
!
 SUBROUTINE GetBoundaryLocation_Hexmesh( myHexMesh, iEl, x, y, z, iBound )
 ! S/R GetBoundaryLocation
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS(Hexmesh), INTENT(in) :: myHexMesh
   INTEGER, INTENT(in)        :: iEl
   INTEGER, INTENT(in)        :: iBound
   REAL(prec), INTENT(out)    :: x(0:myHexMesh % elements(iEl) % nMax, &
                                   0:myHexMesh % elements(iEl) % nMax)
   REAL(prec), INTENT(out)    :: y(0:myHexMesh % elements(iEl) % nMax, &
                                   0:myHexMesh % elements(iEl) % nMax)
   REAL(prec), INTENT(out)    :: z(0:myHexMesh % elements(iEl) % nMax, &
                                   0:myHexMesh % elements(iEl) % nMax)  
  
    CALL myHexMesh % elements(iEl) % GetBoundaryLocation( x, y, z, iBound )
     
 END SUBROUTINE GetBoundaryLocation_Hexmesh
!
!
!
 SUBROUTINE SetBoundaryLocationAtNode_Hexmesh( myHexMesh, iEl, x, y, z, i, j, iBound )
 ! S/R SetBoundaryLocationAtNode
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS(Hexmesh), INTENT(inout) :: myHexMesh
   INTEGER, INTENT(in)           :: iEl
   REAL(prec), INTENT(in)        :: x, y, z
   INTEGER, INTENT(in)           :: i, j, iBound
  
    CALL myHexMesh % elements(iEl) % SetBoundaryLocationAtNode( x, y, z, i, j, iBound )
     
 END SUBROUTINE SetBoundaryLocationAtNode_Hexmesh
!
!
!
 SUBROUTINE GetBoundaryLocationAtNode_Hexmesh( myHexMesh, iEl, x, y, z, i, j, iBound )
 ! S/R GetBoundaryLocationAtNode
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS(Hexmesh), INTENT(in) :: myHexMesh
   INTEGER, INTENT(in)        :: iEl
   REAL(prec), INTENT(out)    :: x, y, z
   INTEGER, INTENT(in)        :: i, j, iBound
  
    CALL myHexMesh % elements(iEl) % GetBoundaryLocationAtNode( x, y, z, i, j, iBound )
     
 END SUBROUTINE GetBoundaryLocationAtNode_Hexmesh
!
!
!
 SUBROUTINE SetBoundaryNormalAtNode_Hexmesh( myHexMesh, iEl, dir, i, j, iBound )
 ! S/R SetBoundaryNormalAtNode
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS(Hexmesh), INTENT(inout) :: myHexMesh
   INTEGER, INTENT(in)           :: iEl
   REAL(prec), INTENT(in)        :: dir(1:nDims)
   INTEGER, INTENT(in)           :: i, j, iBound
  
    CALL myHexMesh % elements(iEl) % SetBoundaryNormalAtNode( dir, i, j, iBound )
     
 END SUBROUTINE SetBoundaryNormalAtNode_Hexmesh
!
!
!
 SUBROUTINE GetBoundaryNormalAtNode_Hexmesh( myHexMesh, iEl, dir, length, i, j, iBound )
 ! S/R GetBoundaryNormalAtNode
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS(Hexmesh), INTENT(in) :: myHexMesh
   INTEGER, INTENT(in)        :: iEl
   REAL(prec), INTENT(out)    :: dir(1:nDims), length
   INTEGER, INTENT(in)        :: i, j, iBound
  
    CALL myHexMesh % elements(iEl) % GetBoundaryNormalAtNode( dir, length, i, j, iBound )
     
 END SUBROUTINE GetBoundaryNormalAtNode_Hexmesh
!
! ------------------------------------------------------------------------------------------------ !
! ------------------------------- Node Wrapper Routines ------------------------------------------ !
! ------------------------------------------------------------------------------------------------ !
!
 SUBROUTINE SetNodeData_HexMesh( myHexMesh, iNode, x, y, z, nodeType )
 ! S/R SetNodeData
 !  
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS( HexMesh ), INTENT(inout) :: myHexMesh
   INTEGER, INTENT(in)             :: iNode
   REAL(prec), INTENT(in)          :: x, y, z
   INTEGER, INTENT(in)             :: nodeType

    CALL myHexMesh % nodes(iNode) % SetData( x, y, z, nodeType )

 END SUBROUTINE SetNodeData_HexMesh
!
!
!
 SUBROUTINE GetNodeData_HexMesh( myHexMesh, iNode, x, y, z, nodeType )
 ! S/R GetNodeData
 !  
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS( HexMesh ), INTENT(in) :: myHexMesh
   INTEGER, INTENT(in)          :: iNode
   REAL(prec), INTENT(out)      :: x, y
   INTEGER, INTENT(out)         :: nodeType

    CALL myHexMesh % nodes(iNode) % GetData( x, y, z, nodeType )

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
   INTEGER, INTENT(in)             :: iNode
   INTEGER, INTENT(in)             :: key
   
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
   INTEGER, INTENT(in)          :: iNode
   INTEGER, INTENT(out)         :: key
   
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
   INTEGER, INTENT(in)             :: iNode
   INTEGER, INTENT(in)             :: nodeType
   
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
   INTEGER, INTENT(in)          :: iNode
   INTEGER, INTENT(out)         :: nodeType
   
    CALL myHexMesh % nodes(iNode) % GetType( nodeType )
   
 END SUBROUTINE GetNodeType_HexMesh
!
!
!
 SUBROUTINE SetNodePosition_HexMesh( myHexMesh, iNode, x, y, z )
 ! S/R SetNodePosition_HexMesh
 !  
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS( HexMesh ), INTENT(inout) :: myHexMesh
   INTEGER, INTENT(in)             :: iNode
   REAL(prec), INTENT(in)          :: x, y, z
   
    CALL myHexMesh % nodes(iNode) % SetPosition( x, y, z )
   
 END SUBROUTINE SetNodePosition_HexMesh
!
!
!
 SUBROUTINE GetNodePosition_HexMesh( myHexMesh, iNode, x, y, z )
 ! S/R GetNodePosition_HexMesh
 !  
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS( HexMesh ), INTENT(in) :: myHexMesh
   INTEGER, INTENT(in)          :: iNode
   REAL(prec), INTENT(out)      :: x, y, z
   
    CALL myHexMesh % nodes(iNode) % GetPosition( x, y, z )
   
 END SUBROUTINE GetNodePosition_HexMesh
!
! ------------------------------------------------------------------------------------------------ !
! --------------------------------- Face Wrapper Routines ---------------------------------------- !
! ------------------------------------------------------------------------------------------------ !
!
 SUBROUTINE SetFaceData_HexMesh( myHexMesh, iFace, nodeIDs, elementIDs, elementSides, &
                                 key, iStart, jStart, iInc, jInc )
 ! S/R SetFaceData
 !  
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS( HexMesh ), INTENT(inout) :: myHexMesh
   INTEGER, INTENT(in)             :: iFace
   INTEGER, INTENT(in)             :: nodeIDs(1:nQuadNodes), elementIDs(1:2), elementSides(1:2)
   INTEGER, INTENT(in)             :: key, iStart, jStart, iInc, jInc

    CALL myHexMesh % Faces(iFace) % SetData( nodeIDs, elementIDs, elementSides, &
                                             key, iStart, jStart, iInc, jInc )

 END SUBROUTINE SetFaceData_HexMesh
!
!
!
 SUBROUTINE GetFaceData_HexMesh( myHexMesh, iFace, nodeIDs, elementIDs, elementSides, &
                                 key, iStart, jStart, iInc, jInc )
 ! S/R GetFaceData
 !  
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS( HexMesh ), INTENT(in) :: myHexMesh
   INTEGER, INTENT(in)          :: iFace
   INTEGER, INTENT(out)         :: nodeIDs(1:nQuadNodes), elementIDs(1:2), elementSides(1:2)
   INTEGER, INTENT(out)         :: key, iStart, jStart, iInc, jInc

    CALL myHexMesh % Faces(iFace) % GetData( nodeIDs, elementIDs, elementSides, &
                                             key, iStart, jStart, iInc, jInc )
      
 END SUBROUTINE GetFaceData_HexMesh
!
!
!
SUBROUTINE SetFaceKey_HexMesh( myHexMesh, iFace, key )
 ! S/R SetFaceKey
 !  
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS( HexMesh ), INTENT(inout) :: myHexMesh
   INTEGER, INTENT(in)             :: iFace
   INTEGER, INTENT(in)             :: key

    CALL myHexMesh % Faces(iFace) % SetKey( key )    

 END SUBROUTINE SetFaceKey_HexMesh
!
!
!
 SUBROUTINE GetFaceKey_HexMesh( myHexMesh, iFace, key )
 ! S/R GetFaceKey
 !  
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS( HexMesh ), INTENT(in) :: myHexMesh
   INTEGER, INTENT(in)          :: iFace
   INTEGER, INTENT(out)         :: key

    CALL myHexMesh % Faces(iFace) % GetKey( key )

 END SUBROUTINE GetFaceKey_HexMesh
!
!
!
 SUBROUTINE SetFaceNodeIDs_HexMesh( myHexMesh, iFace, nodeIDs )
 ! S/R SetFaceNodeIDs
 !  
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS( HexMesh ), INTENT(inout) :: myHexMesh
   INTEGER, INTENT(in)             :: iFace
   INTEGER, INTENT(in)             :: nodeIDs(1:nQuadNodes)

    CALL myHexMesh % Faces(iFace) % SetNodeIDs( nodeIds )

 END SUBROUTINE SetFaceNodeIDs_HexMesh
!
!
!
 FUNCTION GetFaceNodeIDs_HexMesh( myHexMesh, iFace ) RESULT( nodeIDs )
 ! FUNCTION GetFaceNodeIDs
 !  
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS( HexMesh ) :: myHexMesh
   INTEGER          :: iFace
   INTEGER          :: nodeIDs(1:nQuadNodes)

    CALL myHexMesh % Faces(iFace) % GetNodeIDs( nodeIDs )

 END FUNCTION GetFaceNodeIDs_HexMesh
!
!
!
 SUBROUTINE SetFaceElementIDs_HexMesh( myHexMesh, iFace, elementIDs )
 ! S/R SetFaceElementIDs
 !  
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS( HexMesh ), INTENT(inout) :: myHexMesh
   INTEGER, INTENT(in)             :: iFace
   INTEGER, INTENT(in)             :: elementIDs(1:2)

    CALL myHexMesh % Faces(iFace) % SetElementIDs( elementIDs )

 END SUBROUTINE SetFaceElementIDs_HexMesh
!
!
!
 SUBROUTINE GetFaceElementIDs_HexMesh( myHexMesh, iFace, elementIDs )
 ! S/R GetFaceElementIDs
 !  
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS( HexMesh ), INTENT(in) :: myHexMesh
   INTEGER, INTENT(in)          :: iFace
   INTEGER, INTENT(out)         :: elementIDs(1:2)

    CALL myHexMesh % Faces(iFace) % GetElementIDs( elementIDs )

 END SUBROUTINE GetFaceElementIDs_HexMesh
!
!
!
 SUBROUTINE SetFacePrimaryElementID_HexMesh( myHexMesh, iFace, elementID )
 ! S/R SetFacePrimaryElementID
 !  
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS( HexMesh ), INTENT(inout) :: myHexMesh
   INTEGER, INTENT(in)             :: iFace
   INTEGER, INTENT(in)             :: elementID

    CALL myHexMesh % Faces(iFace) % SetPrimaryElementID( elementID )
      
 END SUBROUTINE SetFacePrimaryElementID_HexMesh
!
!
!
 SUBROUTINE GetFacePrimaryElementID_HexMesh( myHexMesh, iFace, elementID )
 ! S/R GetFacePrimaryElementID
 !  
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS( HexMesh ), INTENT(in) :: myHexMesh
   INTEGER, INTENT(in)           :: iFace
   INTEGER, INTENT(out)          :: elementID

    CALL myHexMesh % Faces(iFace) % GetPrimaryElementID( elementID )

 END SUBROUTINE GetFacePrimaryElementID_HexMesh
!
!
!
 SUBROUTINE SetFaceSecondaryElementID_HexMesh( myHexMesh, iFace, elementID )
 ! S/R SetFaceSecondaryElementID
 !  
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS( HexMesh ), INTENT(inout) :: myHexMesh
   INTEGER, INTENT(in)             :: iFace
   INTEGER, INTENT(in)             :: elementID

    CALL myHexMesh % Faces(iFace) % SetSecondaryElementID( elementID )
      
 END SUBROUTINE SetFaceSecondaryElementID_HexMesh
!
!
!
 SUBROUTINE GetFaceSecondaryElementID_HexMesh( myHexMesh, iFace, elementID )
 ! S/R GetFaceSecondaryElementID
 !  
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS( HexMesh ), INTENT(in) :: myHexMesh
   INTEGER, INTENT(in)          :: iFace
   INTEGER, INTENT(out)         :: elementID

    CALL myHexMesh % Faces(iFace) % GetSecondaryElementID( elementID )

 END SUBROUTINE GetFaceSecondaryElementID_HexMesh
!
!
!
 SUBROUTINE SetFaceElementSides_HexMesh( myHexMesh, iFace, elementSides )
 ! S/R SetFaceElementSides
 !  
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS( HexMesh ), INTENT(inout) :: myHexMesh
   INTEGER, INTENT(in)             :: iFace
   INTEGER, INTENT(in)             :: elementSides(1:2)

    CALL myHexMesh % Faces(iFace) % SetElementSides( elementSides )

 END SUBROUTINE SetFaceElementSides_HexMesh
!
!
!
 SUBROUTINE GetFaceElementSides_HexMesh( myHexMesh, iFace, elementSides )
 ! S/R GetFaceElementSides
 !  
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS( HexMesh ), INTENT(in) :: myHexMesh
   INTEGER, INTENT(in)          :: iFace
   INTEGER, INTENT(out)         :: elementSides(1:2)

    CALL myHexMesh % Faces(iFace) % GetElementSides( elementSides )

 END SUBROUTINE GetFaceElementSides_HexMesh
!
!
!
 SUBROUTINE SetFacePrimaryElementSide_HexMesh( myHexMesh, iFace, elementSide )
 ! S/R SetFacePrimaryElementSide
 !  
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS( HexMesh ), INTENT(inout) :: myHexMesh
   INTEGER, INTENT(in)             :: iFace
   INTEGER, INTENT(in)             :: elementSide

    CALL myHexMesh % Faces(iFace) % SetPrimaryElementSide( elementSide )
      
 END SUBROUTINE SetFacePrimaryElementSide_HexMesh
!
!
!
 SUBROUTINE GetFacePrimaryElementSide_HexMesh( myHexMesh, iFace, elementSide )
 ! S/R GetFacePrimaryElementSide
 !  
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS( HexMesh ), INTENT(in) :: myHexMesh
   INTEGER, INTENT(in)          :: iFace
   INTEGER, INTENT(out)         :: elementSide

    CALL myHexMesh % Faces(iFace) % GetPrimaryElementSide( elementSide )

 END SUBROUTINE GetFacePrimaryElementSide_HexMesh
!
!
!
 SUBROUTINE SetFaceSecondaryElementSide_HexMesh( myHexMesh, iFace, elementSide )
 ! S/R SetFaceSecondaryElementSide
 !  
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS( HexMesh ), INTENT(inout) :: myHexMesh
   INTEGER, INTENT(in)             :: iFace
   INTEGER, INTENT(in)             :: elementSide

    CALL myHexMesh % Faces(iFace) % SetSecondaryElementSide( elementSide )
      
 END SUBROUTINE SetFaceSecondaryElementSide_HexMesh
!
!
!
 SUBROUTINE GetFaceSecondaryElementSide_HexMesh( myHexMesh, iFace, elementSide )
 ! S/R GetFaceSecondaryElementSide
 !  
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS( HexMesh ), INTENT(in) :: myHexMesh
   INTEGER, INTENT(in)          :: iFace
   INTEGER, INTENT(out)         :: elementSide

    CALL myHexMesh % Faces(iFace) % GetSecondaryElementSide( elementSide )

 END SUBROUTINE GetFaceSecondaryElementSide_HexMesh
!
!
!
 SUBROUTINE SetFaceStart_HexMesh( myHexMesh, iFace, iStart, jStart )
 ! S/R SetFaceStart
 !  
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS( HexMesh ), INTENT(inout) :: myHexMesh
   INTEGER, INTENT(in)             :: iFace
   INTEGER, INTENT(in)             :: iStart, jStart
    
    CALL myHexMesh % Faces(iFace) % SetStart( iStart, jStart )
      
 END SUBROUTINE SetFaceStart_HexMesh
!
!
!
 SUBROUTINE GetFaceStart_HexMesh( myHexMesh, iFace, iStart, jStart )
 ! S/R GetFaceStart
 !  
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS( HexMesh ), INTENT(in) :: myHexMesh
   INTEGER, INTENT(in)          :: iFace
   INTEGER, INTENT(out)         :: iStart, jStart

    CALL myHexMesh % Faces(iFace) % GetStart( iStart, jStart )

 END SUBROUTINE GetFaceStart_HexMesh
!
!
!
 SUBROUTINE SetFaceIncrement_HexMesh( myHexMesh, iFace, iInc, jInc )
 ! S/R SetFaceIncrement
 !  
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS( HexMesh ), INTENT(inout) :: myHexMesh
   INTEGER, INTENT(in)             :: iFace
   INTEGER, INTENT(in)             :: iInc, jInc

    CALL myHexMesh % Faces(iFace) % SetIncrement( iInc, jInc )

 END SUBROUTINE SetFaceIncrement_HexMesh
!
!
!
 SUBROUTINE GetFaceIncrement_HexMesh( myHexMesh, iFace, iInc, jInc )
 ! S/R GetFaceIncrement
 !  
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS( HexMesh ), INTENT(in) :: myHexMesh
   INTEGER, INTENT(in)          :: iFace
   INTEGER, INTENT(out)         :: iInc, jInc

    CALL myHexMesh % Faces(iFace) % GetIncrement( iInc, jInc )

 END SUBROUTINE GetFaceIncrement_HexMesh
!
!
!
 SUBROUTINE SetSwapDimensions_HexMesh( myHexMesh, iFace, swap )
 ! S/R SetSwapDimensions
 !  
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS( HexMesh ), INTENT(inout) :: myHexMesh
   INTEGER, INTENT(in)             :: iFace
   INTEGER, INTENT(in)             :: swap

    CALL myHexMesh % Faces(iFace) % SetSwapDimensions( swap )

 END SUBROUTINE SetSwapDimensions_HexMesh
!
!
!
 SUBROUTINE GetSwapDimensions_HexMesh( myHexMesh, iFace, swap )
 ! S/R GetSwapDimensions
 !  
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS( HexMesh ), INTENT(in) :: myHexMesh
   INTEGER, INTENT(in)          :: iFace
   INTEGER, INTENT(out)         :: swap

    CALL myHexMesh % Faces(iFace) % GetSwapDimensions( swap )

 END SUBROUTINE GetSwapDimensions_HexMesh
!
!
!==================================================================================================!
!--------------------------------- Type Specific Routines -----------------------------------------!
!==================================================================================================!
!
!
 SUBROUTINE ConstructFaces_HexMesh( myHexMesh )
 ! S/R ConstructFaces
 !
 !    Takes in the mesh of hexahedrons which has not filled in the 
 !    Face information, and finds all of the unique Faces in the mesh.
 !    The space for the Faces is REALlocated with the correct number of Faces
 ! 
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS( HexMesh ), INTENT(inout) :: myHexMesh
   ! LOCAL
   TYPE( KeyRing ) :: KeyCabinet(1:myHexMesh % nNodes)
   INTEGER :: nEls, nNodes, iEl, nFaces, k, j  
   INTEGER :: locNodeIDs(1:nQuadNodes), globNodeIDs(1:nQuadNodes)
   INTEGER :: keyRingID, globFaceID
   INTEGER :: e1, e2, s1, s2, FaceID, n1, nID
   LOGICAL :: keyExists

      CALL myHexMesh % GetNumberOfNodes( nNodes )   
      CALL myHexMesh % GetNumberOfElements( nEls )
      nFaces = 0

      DO k = 1, nNodes
         CALL keyCabinet(keyRingID) % Build( )
      ENDDO

      DO iEl = 1, nEls ! Loop over the elements in the mesh

         DO k = 1, nHexFaces ! Loop over the faces of each element

            ! In this first step, we want to identify the unique global node ID's for this face
            ! To do this, we start by gathering the local (to the element) node ID's.
            DO j = 1, nQuadNodes   
               locNodeIDs(j) = myHexMesh % faceMap(j,k) ! starting local node for this Face
            ENDDO
            
            ! Now, we extract, for this element, the global node ID's for each local node ID
            DO j = 1, nQuadNodes 
               globNodeIDs(j) = myHexMesh % GetElementNodeID( iEl, locNodeIDs(j) )
            ENDDO
            
            ! Our key cabinet has many key-rings with each key-ring containing a set of notched keys.
            ! Each notched key corresponds to a unique face in the mesh. To enable fast searching
            ! for a unique face, we address our key-rings according to the minimum global node ID
            ! that resides on the current face. If another face shares this minimum global node ID,
            ! then it is possible that the face has already been generated. If not, our search is 
            ! limited only to the key-ring with the same key-ring ID. Here, we grab the key-ring ID.
            keyRingID = MIN( globNodeIDs )
            
            ! Now, we check to see if a notched key already exists with the same set of global ID's
            ! This is done by making a call to "FindDataForNotches". This routine searches through
            ! the current key ring for the key which has the same global ID's (though not 
            ! necessarily in the same order). If the face has been built already, the face ID is 
            ! returned in globFaceID and the logical, "keyExists", is set to TRUE. Otherwise, 
            ! globFaceID is set to zero and "keyExsists" is set to FALSE.
            CALL KeyCabinet(keyRingID) % FindDataForNotches( globNodeIDS, nQuadNodes, &
                                                             globFaceID, keyExists )

            ! Here is where the conditional processing begins. 
            ! 
            ! If this face has already been found then we do nothing
            !
            ! If this is a new face, we increment the number of faces and add to the key-ring.

            ! Add element to corner-node connectivity list
            IF( .NOT.(keyExists) )then ! this is a new face

               nFaces = nFaces + 1
               CALL KeyCabinet(keyRingID) % AddToList( nFaces, globNodeIDs, nQuadNodes )
               
            ENDIF

         ENDDO ! k, Loop over the faces of each element
        
      ENDDO! iEl, Loop over the elements in the mesh
 

      DO k = 1, nNodes
         CALL keyCabinet(keyRingID) % Trash( ) ! Trash the Facetable
         CALL keyCabinet(keyRingID) % Build( ) ! and rebuild a blank cabinet 
      ENDDO
      
      ! Re-allocate space for the mesh Faces

      DEALLOCATE( myHexMesh % Faces )

      ALLOCATE( myHexMesh % Faces( 1:nFaces ) )

      nFaces = 0

      DO iEl = 1, nEls ! Loop over the elements in the mesh

         DO k = 1, nHexFaces ! Loop over the faces of each element

            ! In this first step, we want to identify the unique global node ID's for this face
            ! To do this, we start by gathering the local (to the element) node ID's.
            DO j = 1, nQuadNodes   
               locNodeIDs(j) = myHexMesh % faceMap(j,k) ! starting local node for this Face
            ENDDO
            
            ! Now, we extract, for this element, the global node ID's for each local node ID
            DO j = 1, nQuadNodes 
               globNodeIDs(j) = myHexMesh % GetElementNodeID( iEl, locNodeIDs(j) )
            ENDDO
            
            ! Our key cabinet has many key-rings with each key-ring containing a set of notched keys.
            ! Each notched key corresponds to a unique face in the mesh. To enable fast searching
            ! for a unique face, we address our key-rings according to the minimum global node ID
            ! that resides on the current face. If another face shares this minimum global node ID,
            ! then it is possible that the face has already been generated. If not, our search is 
            ! limited only to the key-ring with the same key-ring ID. Here, we grab the key-ring ID.
            keyRingID = MIN( globNodeIDs )
            
            ! Now, we check to see if a notched key already exists with the same set of global ID's
            ! This is done by making a call to "FindDataForNotches". This routine searches through
            ! the current key ring for the key which has the same global ID's (though not 
            ! necessarily in the same order). If the face has been built already, the face ID is 
            ! returned in globFaceID and the logical, "keyExists", is set to TRUE. Otherwise, 
            ! globFaceID is set to zero and "keyExsists" is set to FALSE.
            CALL KeyCabinet(keyRingID) % FindDataForNotches( globNodeIDS, nQuadNodes, &
                                                             globFaceID, keyExists )

            ! Here is where the conditional processing begins. 
            ! 
            ! If this face has already been found then we need to determine the face orientation
            ! of the secondary element relative to the primary element.
            !
            ! If this is a new face, we set the primary element information, and default the
            ! secondary element information.

            ! Add element to corner-node connectivity list
            IF( keyExists )then ! this face has already been found

               ! Since this face exists, we need to compare the relative face orientation of the 
               ! secondary element to the primary element. .
               !
               ! *** March 4, 2016 ***
               ! Joe : I still need to write this routine that determines the orientation of the 
               !       secondary element.
               CALL myHexMesh % DetermineOrientation( globFaceID, globNodeIDs )

               myHexMesh % faces( globFaceID ) % elementIDs(2)   = iEl
               myHexMesh % faces( globFaceID ) % elementSides(2) = k
               
            ELSE ! This is a new face

               ! First, we store the key-ring information
               nFaces = nFaces + 1
               CALL KeyCabinet(keyRingID) % AddToList( nFaces, globNodeIDs, nQuadNodes )

               ! Now, we set the primary element information
               myHexMesh % faces( nFaces ) % key             = nFaces
               myHexMesh % faces( nFaces ) % nodeIDs         = globNodeIDs
               myHexMesh % faces( nFaces ) % elementIDs(1)   = iEl
               myHexMesh % faces( nFaces ) % elementSides(1) = k
               ! Now we default the secondary element information and the swap flag
               myHexMesh % faces( nFaces ) % elementIDs(2)   = BoundaryFlagDefault
               myHexMesh % faces( nFaces ) % elementSides(2) = k
               myHexMesh % faces( nFaces ) % iStart          = 1
               myHexMesh % faces( nFaces ) % iInc            = 1
               myHexMesh % faces( nFaces ) % jStart          = 1
               myHexMesh % faces( nFaces ) % jInc            = 1
               myHexMesh % faces( nFaces ) % swapDimensions  = 0

            ENDIF

         ENDDO ! k, Loop over the faces of each element
        
      ENDDO! iEl, Loop over the elements in the mesh
      CALL FaceTable % Trash( )
      myHexMesh % nFaces = nFaces


 END SUBROUTINE ConstructFaces_HexMesh
!
!
!
 SUBROUTINE DetermineOrientation_HexMesh( myHexMesh, faceID, secondaryNodes ) 
 ! S/R DetermineOrientation
 !
 !  This support routine takes as input the mesh, a face ID number, and a list of node IDs.
 !  This routine assumes that the primary element information for the given face ID has been
 !  filled upon issuing a call to this routine. If the primary element shares this face with
 !  another element, this routine is called and the secondary element node ID's are passed in.
 !  The order of the secondary node ID's relative to the order of the primary node ID's determines
 !  the orientation of the secondary element relative to the primary element. We need to know
 !  if the roles of the computational coordinates are flipped, and (simultaneously) if the 
 !  incrementing of the computational coordinates are reversed (or not).
 !  This routine determines the orientation by determining how many times the ordered secondary
 !  node ID's need to be shifted in order to match the ordered primary node ID's.
 !  
 ! 
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS( HexMesh ), INTENT(inout) :: myHexMesh
   INTEGER, INTENT(in)             :: faceID
   INTEGER, INTENT(in)             :: secondaryNodes(1:nQuadNodes)
   ! Local
   INTEGER :: primaryNodes(1:nQuadNodes)
   INTEGER :: nShifts, shiftInc, i 
   
      primaryNodes = myHexMesh % GetFaceNodeIDs( faceID )

      nShifts = 0

      DO i = 1, nQuadNodes

         ! First, we compare the primary and secondary nodes. This routine returns a zero 
         ! if the arrays match, and a one if they do not match.
         shiftInc = CompareArray( primaryNodes, secondaryNodes, nQuadNodes )
         
         IF( shiftInc == 0 )THEN
            EXIT
         ELSE
            
      

         
      ENDDO
   

 END SUBROUTINE DetermineOrientation_HexMesh
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
   INTEGER :: nEls, nNodes, iEl, nFaces, k  
   INTEGER :: l1, l2, startID, endID, key1, key2
   INTEGER :: e1, e2, s1, s2, FaceID, n1, nID

      CALL myHexMesh % GetNumberOfNodes( nNodes )   
      CALL myHexMesh % GetNumberOfElements( nEls )
      nFaces = 0

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
   INTEGER :: e(1:2), s(1:2), nFaces, iFace
  
      CALL myHexMesh % GetNumberOfFaces( nFaces )
 
      DO iFace = 1, nFaces ! Loop over the Faces

         CALL myHexMesh % GetFaceElementIDs( iFace, e )
         CALL myHexMesh % GetFaceElementSides( iFace, s )

         CALL myHexMesh % SetElementNeighbor( e(1), s(1), e(2) )
         CALL myHexMesh % SetElementNeighbor( e(2), ABS(s(2)), e(1) )
         
      ENDDO ! iFace, Loop over the Faces

     
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

   INTEGER :: nNodes, nElems, nFaces, gPolyDeg
   INTEGER :: bFlags(1:4), nodes(1:4)
   INTEGER :: e1, e2, s1, s2, n1, n2, n(1:2), e(1:2), si(1:2)
   INTEGER :: iFace, iNode, iEl, iSide, iX, iY
   INTEGER :: fUnit, iC, jC

   CHARACTER(20) :: ISMversion
   CHARACTER(40) :: FaceNames
   CHARACTER(40) :: thisFace
      
      dxElem = ONE/nXElem
      dyElem = ONE/nYElem
      
      ! ** "Hard-wired" values for a structured mesh with no holes ** !
      nNodes   = (nXElem+1)*(nYElem+1)
      nElems   = (nXElem)*(nYElem)
      nFaces   = (nXElem)*(nYElem+1) + (nXElem+1)*(nYElem)
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
  
      CALL myHexMesh % Build( nNodes, nElems, nFaces, interp % nS ) 
      
      
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

      CALL myHexMesh % ConstructFaces( )
      nFaces = myHexMesh % nFaces
      PRINT*, 'nFaces    : ', nFaces
      
      ! Set the start and increment for the secondary element 
      DO iFace = 1, nFaces
            CALL myHexMesh % GetFaceSecondaryElementSide(iFace, s2)
            IF(s2 < 0)THEN
               CALL myHexMesh % SetFaceStart( iFace, interp % nS-1 )
               CALL myHexMesh % SetFaceIncrement( iFace, -1)
            ELSE
               CALL myHexMesh % SetFaceStart( iFace, 1 )
               CALL myHexMesh % SetFaceIncrement( iFace, 1 )
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

   INTEGER :: nNodes, nElems, nFaces, gPolyDeg
   INTEGER :: bFlags(1:4), nodes(1:4)
   INTEGER :: e1, e2, s1, s2, n1, n2, n(1:2), e(1:2), si(1:2)
   INTEGER :: iFace, iNode, iEl, iSide
   INTEGER :: fUnit, iC, jC
   INTEGER, ALLOCATABLE :: sideFlags(:,:)

   CHARACTER(20) :: ISMversion
   CHARACTER(40) :: FaceNames
   CHARACTER(40) :: thisFace

      PRINT*, 'Mesh File : '//TRIM( filename )
      
      ! Get a new file unit
      OPEN( UNIT=NEWUNIT(fUnit), &
            FILE= TRIM( filename ), &
            FORM='FORMATTED',&
            STATUS='OLD' )

      
      ! ---- Read in the file version ----- !

      READ( fUnit, '(A20)' ) ISMversion

      PRINT*, 'Version   : '//TRIM( ISMversion )


      ! ---- Gather the number of nodes, number of elements, and number of Faces ---- !
      
      READ( fUnit, * ) nNodes, nFaces, nElems, gPolyDeg

      PRINT*, 'nNodes    : ', nNodes
      PRINT*, 'nElems    : ', nElems
      PRINT*, 'nFaces    : ', nFaces
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
  
      CALL myHexMesh % Build( nNodes, nElems, nFaces, interp % nS ) 
      
      ! ---- Read in the corner nodes ---- !

      DO iNode = 1, nNodes  ! Loop over the nodes in the file
         READ( fUnit, * ) x, y, z
         CALL myHexMesh % SetNodePosition( iNode, x, y )
      ENDDO


      ! ---- Read in the Face information ---- !

      DO iFace = 1, nFaces

         READ( fUnit, * ) n, e, si 
         
         !CALL myHexMesh % SetFaceNodeIDs( iFace, n )
         !CALL myHexMesh % SetFaceElementIDs( iFace, e )
         !CALL myHexMesh % SetFaceElementSides( iFace, si )
       
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

            
            IF( bFlags(iSide) == 0 )then ! this is an interior Face

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

            ELSEIF( bFlags(iSide) == 1 )then ! this is a boundary Face

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

         ! Read in and parse the Face names
         READ( fUnit, '(1x, A40)' )  FaceNames

         ! Parse the Face names into four Face names
         iSide = 1
         iC = 1
         jC = 1

         DO while( iSide <= 4 )

            IF( FaceNames(iC:iC) == ' ' )then ! we have reached a blank space
     
               n1 = nodes( myHexMesh % faceMap(1,iSide) )
               n2 = nodes( myHexMesh % faceMap(2,iSide) )
               
               IF( TRIM(thisFace) == 'DIRICHLET' )then
                  
                  sideFlags(iEl,iSide) = DIRICHLET

               ELSEIF( TRIM(thisFace) == 'ROBIN' )then

                  sideFlags(iEl,iSide) = ROBIN

               ELSEIF( TRIM(thisFace) == 'ROBIN_FORCED' )then

                  sideFlags(iEl,iSide) = ROBIN_FORCED

               ELSEIF( TRIM(thisFace) == 'HOMOGENEOUS_NEUMANN' )then

                  sideFlags(iEl,iSide) = HOMOGENEOUS_NEUMANN

               ELSEIF( TRIM(thisFace) == 'NEUMANN_WALL' )then

                  sideFlags(iEl,iSide) = NEUMANN_WALL

               ELSEIF( TRIM(thisFace) == 'NEUMANN' )then

                  sideFlags(iEl,iSide) = NEUMANN
 
               ELSEIF( TRIM(thisFace) == 'NO_NORMAL_FLOW' )then
                  
                  sideFlags(iEl,iSide) = NO_NORMAL_FLOW

               ELSEIF( TRIM(thisFace) == 'PRESCRIBED' )then
               
                  sideFlags(iEl,iSide) = PRESCRIBED

               ELSEIF( TRIM(thisFace) == 'RADIATION' )then
                  
                  sideFlags(iEl,iSide) = RADIATION

               ENDIF

               ! Reset thisFace
               thisFace = ' '
               jC = 1 
               iC = iC + 1
               ! Increment the side
               iSide = iSide + 1

            ELSE
            
               thisFace(jC:jC) = FaceNames(iC:iC)
               iC = iC + 1
               jC = jC + 1

            ENDIF

         ENDDO ! while (iSide <= 4 )


      ENDDO ! iEl, cycle over the elements

      close( fUnit )
            CALL myHexMesh % ConstructFaces( )
      ! Set the secondary element on boundary Faces to the boundary flag
      DO iFace = 1, nFaces

         CALL myHexMesh % GetFacePrimaryElementID( iFace, e1 )
         CALL myHexMesh % GetFacePrimaryElementSide( iFace, s1 )
      
         IF( sideFlags(e1,s1) == DIRICHLET )then

            CALL myHexMesh % SetFaceSecondaryElementID( iFace, DIRICHLET )
            CALL myHexMesh % GetFaceNodeIDs( iFace, n )
            CALL myHexMesh % SetNodeType( n(1), DIRICHLET )
            CALL myHexMesh % SetNodeType( n(2), DIRICHLET )

         ELSEIF( sideFlags(e1,s1) == ROBIN )then

            CALL myHexMesh % SetFaceSecondaryElementID( iFace, ROBIN )
            CALL myHexMesh % GetFaceNodeIDs( iFace, n )
            CALL myHexMesh % SetNodeType( n(1), ROBIN )
            CALL myHexMesh % SetNodeType( n(2), ROBIN )

         ELSEIF( sideFlags(e1,s1) == ROBIN_FORCED )then

            CALL myHexMesh % SetFaceSecondaryElementID( iFace, ROBIN_FORCED )
            CALL myHexMesh % GetFaceNodeIDs( iFace, n )
            CALL myHexMesh % SetNodeType( n(1), ROBIN_FORCED )
            CALL myHexMesh % SetNodeType( n(2), ROBIN_FORCED )

         ELSEIF( sideFlags(e1,s1) == HOMOGENEOUS_NEUMANN )then

            CALL myHexMesh % SetFaceSecondaryElementID( iFace, HOMOGENEOUS_NEUMANN )
            CALL myHexMesh % GetFaceNodeIDs( iFace, n )
            CALL myHexMesh % SetNodeType( n(1), HOMOGENEOUS_NEUMANN )
            CALL myHexMesh % SetNodeType( n(2), HOMOGENEOUS_NEUMANN )

         ELSEIF( sideFlags(e1,s1) == NEUMANN_WALL )then

            CALL myHexMesh % SetFaceSecondaryElementID( iFace, NEUMANN_WALL )
            CALL myHexMesh % GetFaceNodeIDs( iFace, n )
            CALL myHexMesh % SetNodeType( n(1), NEUMANN_WALL )
            CALL myHexMesh % SetNodeType( n(2), NEUMANN_WALL )

         ELSEIF( sideFlags(e1,s1) == NEUMANN )then

            CALL myHexMesh % SetFaceSecondaryElementID( iFace, NEUMANN )
            CALL myHexMesh % GetFaceNodeIDs( iFace, n )
            CALL myHexMesh % SetNodeType( n(1), NEUMANN )
            CALL myHexMesh % SetNodeType( n(2), NEUMANN )

         ELSEIF( sideFlags(e1,s1) == NO_NORMAL_FLOW )then

            CALL myHexMesh % SetFaceSecondaryElementID( iFace, NO_NORMAL_FLOW )
            CALL myHexMesh % GetFaceNodeIDs( iFace, n )
            CALL myHexMesh % SetNodeType( n(1), NO_NORMAL_FLOW )
            CALL myHexMesh % SetNodeType( n(2), NO_NORMAL_FLOW )

         ELSEIF( sideFlags(e1,s1) == PRESCRIBED )then

            CALL myHexMesh % SetFaceSecondaryElementID( iFace, PRESCRIBED )
            CALL myHexMesh % GetFaceNodeIDs( iFace, n )
            CALL myHexMesh % SetNodeType( n(1), PRESCRIBED )
            CALL myHexMesh % SetNodeType( n(2), PRESCRIBED )

         ELSEIF( sideFlags(e1,s1) == RADIATION )then

            CALL myHexMesh % SetFaceSecondaryElementID( iFace, RADIATION )
            CALL myHexMesh % GetFaceNodeIDs( iFace, n )
            CALL myHexMesh % SetNodeType( n(1), RADIATION )
            CALL myHexMesh % SetNodeType( n(2), RADIATION )
         
         ENDIF
 
      ENDDO

      DO iFace = 1, nFaces
            CALL myHexMesh % GetFaceSecondaryElementSide(iFace, s2)
            IF(s2 < 0)THEN
               CALL myHexMesh % SetFaceStart( iFace, interp % nS-1 )
               CALL myHexMesh % SetFaceIncrement( iFace, -1)
            ELSE
               CALL myHexMesh % SetFaceStart( iFace, 1 )
               CALL myHexMesh % SetFaceIncrement( iFace, 1 )
            ENDIF
            
      ENDDO
 ! If a node is at the corner of any boundary condition and a dirichlet boundary condition,
 ! we opt to keep it as a dirichlet boundary condition
    !  DO iFace = 1, nFaces

    !     CALL myHexMesh % Faces(iFace) % GetPrimaryElementID( e1 )
    !     CALL myHexMesh % Faces(iFace) % GetPrimaryElementSide( s1 )

    !     IF( sideFlags(e1,s1) == DIRICHLET )then

    !        CALL myHexMesh % Faces(iFace) % SetSecondaryElementID( DIRICHLET )
            
    !        CALL myHexMesh % Faces(iFace) % GetNodeIDs( n )
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
