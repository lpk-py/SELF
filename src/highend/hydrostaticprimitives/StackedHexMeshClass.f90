! StackedHexMeshClass.f90 ( new with v2.1 - 29 March 2016)
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



MODULE StackedHexMeshClass


! src/common/
USE ModelPrecision
USE ConstantsDictionary
USE ModelFlags
! src/interp/
USE Lagrange_2D_Class
USE Lagrange_3D_Class
! src/geom/
USE QuadMeshClass
USE MappedGeometryClass_3D
USE SurfaceClass_3D
USE HexElementClass
USE FaceClass
USE NodeClass_3D
USE HexMeshClass



IMPLICIT NONE


! All properties are left as public in version 2.1 in order to allow the programmer the ability to
! have direct access to the attributes of the mesh class. This way, it can later be determined if
! additional accessor routines would be necessary.
 
    TYPE, EXTENDS( HexMesh ) :: StackedHexMesh 
   !    INTEGER                          :: nElems, nNodes, nFaces
   !    TYPE( HexElement ), ALLOCATABLE  :: elements(:)
   !    TYPE( Node ), ALLOCATABLE        :: nodes(:)  
   !    TYPE( Face ), ALLOCATABLE        :: faces(:)
   !    INTEGER                          :: cornerMap(1:3,1:nHexNodes) 
   !    INTEGER                          :: sideMap(1:nHexFaces) 
   !    INTEGER                          :: faceMap(1:nQuadNodes,1:nHexFaces) 
   !    INTEGER                          :: edgeFaceMap(1:2,1:nQuadEdges)
        INTEGER              :: nHElems, nLayers
        INTEGER, ALLOCATABLE :: stackMap(:,:)

       CONTAINS

       PROCEDURE :: Initialize => Initialize_StackedHexMesh
       PROCEDURE :: Trash => Trash_StackedHexMesh
       
       PROCEDURE :: TerrainFollowingMesh => TerrainFollowingMesh_StackedHexMesh
       
       
    END TYPE StackedHexMesh

 INTEGER, PRIVATE, PARAMETER    :: nDims = 3
 INTEGER, PRIVATE, PARAMETER    :: BoundaryFlagDefault = NO_NORMAL_FLOW



 CONTAINS
!
!
!==================================================================================================!
!------------------------------- Manual Constructors/Destructors ----------------------------------!
!==================================================================================================!
!
!
 SUBROUTINE Initialize_StackedHexMesh( myHexMesh, nLayers, nS, thisQuadMesh )
 ! S/R Initialize
 ! 
 !  A "StackedHexMesh" is a 3-D mesh of hexahedrons that is constructed by extruding a 2-D mesh.
 !  Because of this, a pre-constructed quadrilateral mesh in 2-D and the number of vertical layers
 !  is needed to allocate space and initialize the hexmesh.
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS(StackedHexMesh), INTENT(out) :: myHexMesh
   INTEGER, INTENT(in)                :: nLayers, nS
   TYPE(QuadMesh), INTENT(in)         :: thisQuadMesh
   !LOCAL
   INTEGER :: iNode, nHElems, nHnodes, nHEdges, nElems, nNodes, nFaces

    nHElems = thisQuadMesh % nElems
    nHNodes = thisQuadMesh % nNodes
    nHEdges = thisQuadMesh % nEdges

    myHexMesh % nHElems = nHElems
    myHexMesh % nLayers = nLayers

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
    
    nNodes = nHNodes*(nLayers+1)
    nElems = nHElems*(nLayers)
    nFaces = nHElems*(nLayers+1) + nHEdges*(nLayers)
 
    CALL myHexMesh % SetNumberOfNodes( nNodes )
    CALL myHexMesh % SetNumberOfElements( nElems )
    CALL myHexMesh % SetNumberOfFaces( nFaces )

    ALLOCATE( myHexmesh % elements(1:nElems) )
    ALLOCATE( myHexmesh % nodes(1:nNodes) )
    ALLOCATE( myHexmesh % Faces(1:nFaces) )

    ALLOCATE( myHexmesh % stackMap(1:nHElems,1:nLayers) )
     
      ! Default nodes to origin
    DO iNode = 1, myHexmesh % nNodes ! loop over the number of nodes
       CALL myHexmesh % nodes(iNode) % Build( ZERO, ZERO, ZERO )
    ENDDO ! iNode, loop over the number of nodes
     

 END SUBROUTINE Initialize_StackedHexMesh
!
!
!
 SUBROUTINE Trash_StackedHexMesh( myHexMesh )
 ! S/R Trash
 ! 
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS(StackedHexMesh), INTENT(inout) :: myHexMesh
  ! LOCAL
   INTEGER :: iNode, iEl

    DO iEl = 1, myHexMesh % nElems
       CALL myHexMesh % elements(iEl) % Trash( )
    ENDDO
     
    DO iNode = 1, myHexMesh % nNodes
       CALL myHexMesh % nodes(iNode) % Trash( )
    !   PRINT*, iNode
    ENDDO
    
    DEALLOCATE( myHexMesh % nodes, myHexMesh % elements, myHexMesh % Faces )
    DEALLOCATE( myHexMesh % stackMap )  

 END SUBROUTINE Trash_StackedHexMesh
!
!
!==================================================================================================!
!--------------------------------------- Accessors ------------------------------------------------!
!==================================================================================================!
!
! 
 
!
!
!==================================================================================================!
!--------------------------------- Type Specific Routines -----------------------------------------!
!==================================================================================================!
!
!
 SUBROUTINE TerrainFollowingMesh_StackedHexMesh( myHexMesh, interp, thisQuadMesh, quadInterp, nLayers, h  )
 ! S/R TerrainFollowingMesh
 !  
 !   This subroutine uses the quadrilateral mesh (in 2-D), and the field "h" (for the intended 
 !   purpose, h is the bathymetry field), to construct the 3-D mesh. The 3-D mesh is constructed
 !   by extruding the 2-D mesh from z=0 to z=-h(x,y) so that a "terrain-following" mesh is built.
 !   ! **BE AWARE** ! 
 !       The variable "thisQuadMesh" should be constructed at Gauss-Lobatto points so that the
 !       element edge and corner node locations are included in the mesh locations. Consequently,
 !       "quadInterp" should reflect this in the computational domain.
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS( StackedHexMesh ), INTENT(inout) :: myHexMesh
   TYPE( Lagrange_3D ), INTENT(in)        :: interp
   TYPE( QuadMesh ), INTENT(in)           :: thisQuadMesh
   TYPE( Lagrange_2D ), INTENT(in)        :: quadInterp
   INTEGER, INTENT(in)                    :: nLayers
   REAL(prec), INTENT(in)                 :: h(0:quadInterp % nS, &
                                               0:quadInterp % nS, &
                                               1:thisQuadMesh % nElems)
   ! LOCAL
   TYPE( Surface_3D ) :: boundSurfs(1:nHexFaces)

   REAL(prec) :: x, y, z
   REAL(prec) :: x1, x2, x3, x4, y1, y2, y3, y4, z1, z2, z3, z4
   REAL(prec) :: c1, c2
   REAL(prec), ALLOCATABLE :: xc(:,:), yc(:,:), zc(:,:), s(:), p(:)
   REAL(prec), ALLOCATABLE :: x2d(:,:), y2d(:,:)
   REAL(prec) :: zCorner(1:nHexNodes), locH(1:nQuadNodes)
   INTEGER :: nNodes, nElems, nFaces, gPolyDeg, nHElems, nHNodes
   INTEGER :: nodes(1:nHexNodes), quadNodeIDs(1:nQuadNodes)
   INTEGER :: nID, elID, n1, n2, n3, n4, nS
   INTEGER :: iLayer, iNode, iEl, iSide, i, j

      gPolyDeg = quadInterp % nS
      nS       = interp % nS
      nHElems  = thisQuadMesh % nElems
      nHNodes  = thisQuadMesh % nNodes


      ! ---- Initialize the quadrature mesh (empty) ---- !
      CALL myHexMesh % Initialize( nLayers, nS, thisQuadMesh )

      nNodes   = myHexMesh % nNodes
      nElems   = myHexMesh % nElems
      nFaces   = myHexMesh % nFaces

      ! ************************************************************************* !
      PRINT*, ' ======================================= '
      PRINT*, ' S/R TerrainFollowingMesh : Extruding 2-D mesh '
      PRINT*, 'nNodes    : ', nNodes
      PRINT*, 'nElems    : ', nElems
      PRINT*, 'nFaces    : ', nFaces
      PRINT*, 'gPolyDeg  : ', gPolyDeg

      
      ALLOCATE( s(0:gPolyDeg), p(0:gPolyDeg) ) 
      ALLOCATE( xc(0:gPolyDeg,0:gPolyDeg), yc(0:gPolyDeg,0:gPolyDeg), zc(0:gPolyDeg,0:gPolyDeg) )
      ALLOCATE( x2d(0:gPolyDeg,0:gPolyDeg), y2d(0:gPolyDeg,0:gPolyDeg) )
       
      CALL quadInterp % GetNodes( s, p )
      ! ---- Set the corner nodes ---- !
      !
      DO iLayer = 0, nLayers
         DO iNode = 1, nHNodes
            nID = iNode + iLayer*nHNodes
            CALL thisQuadMesh % GetNodePosition( iNode, x, y )
            ! For now, the vertical position is zet to zero. In the loop over the elements, the
            ! vertical position will be added
            z = ZERO
            CALL myHexMesh % SetNodePosition( nID, x, y, z )
         ENDDO
      ENDDO
  
      ! Do the element information
 
      xc = ZERO
      yc = ZERO
      zc = ZERO
      ! Do the initial build for the parametric surfaces
      DO iSide = 1, nHexFaces
         CALL boundSurfs(iSide) % Build( xc, yc, zc, s, p ) 
      ENDDO
   
      DO iLayer = 1, nLayers
         DO iEl = 1, nHElems
      
            elID = iEl + (iLayer-1)*(nHElems)

            ! The "stackMap" stores the "global" element IDs for elements that are stacked on top 
            ! each other -- these elements have the same lateral positions
            myHexmesh % stackMap(iEl,iLayer) = elID

            CALL thisQuadMesh % GetElementNodeIDs( iEl, quadNodeIDs )

            ! Calculate the global node IDs for this element.
            nodes(1) = quadNodeIDs(1) + (iLayer-1)*nHNodes   ! Southwest-"Bottom"
            nodes(2) = quadNodeIDs(2) + (iLayer-1)*nHNodes   ! SouthEast-"Bottom"
            nodes(3) = quadNodeIDs(3) + (iLayer-1)*nHNodes   ! NorthEast-"Bottom"
            nodes(4) = quadNodeIDs(4) + (iLayer-1)*nHNodes   ! NorthWest-"Bottom"
            nodes(5) = quadNodeIDs(1) + iLayer*nHNodes       ! Southwest-"Top"
            nodes(6) = quadNodeIDs(2) + iLayer*nHNodes       ! SouthEast-"Top"
            nodes(7) = quadNodeIDs(3) + iLayer*nHNodes       ! NorthEast-"Top"
            nodes(8) = quadNodeIDs(4) + iLayer*nHNodes       ! NorthWest-"Top"

            DO iNode = 1, nQuadNodes
               i = thisQuadMesh % cornerMap(1,iNode)
               j = thisQuadMesh % cornerMap(2,iNode)
               locH(iNode) = h(i,j,iEl)
            ENDDO
 
            ! Sigma-coordinate mapping for the vertical position of the corner nodes
            zCorner(1:4) = ( locH/REAL(nLayers,prec) )*REAL(iLayer-1,prec) - locH ! Bottom
            zCorner(5:8) = ( locH/REAL(nLayers,prec) )*REAL(iLayer,prec) - locH   ! Top

            DO iNode = 1, nHexNodes
               CALL myHexMesh % GetNodePosition( nodes(iNode), x, y, z )
               ! Reassign z
               z = zCorner(iNode)
               CALL myHexMesh % SetNodePosition( nodes(iNode), x, y, z )
            ENDDO

            CALL thisQuadMesh % GetPositions( iEl, x2d, y2d )

            DO iSide = 1, nHexFaces ! Loop over the sides of the quads

               ! To build the current face, we construct a plane that passes through
               ! the four corner nodes. Here, we grab the global node ID's for the four
               ! corner nodes.
               n1 = nodes( myHexMesh % faceMap(1,iSide) )
               n2 = nodes( myHexMesh % faceMap(2,iSide) )
               n3 = nodes( myHexMesh % faceMap(3,iSide) ) 
               n4 = nodes( myHexMesh % faceMap(4,iSide) ) 

               CALL myHexMesh % GetNodePosition( n1, x1, y1, z1 )
               CALL myHexMesh % GetNodePosition( n2, x2, y2, z2 )
               CALL myHexMesh % GetNodePosition( n3, x3, y3, z3 )
               CALL myHexMesh % GetNodePosition( n4, x4, y4, z4 )
               
               IF( iLayer == 1) THEN

                  IF( iSide == south .OR. iSide == north )THEN  

                     nID = thisQuadMesh % sideMap(iSide)
                     DO j = 0, gPolyDeg
                        DO i = 0, gPolyDeg
                  
                           xc(i,j) = x2d(i,nID) 
                           yc(i,j) = y2d(i,nID)

                           c1 = -h(i,nID,iEl) ! Bottom curve of this face is the bathymetry 
                           c2 = ( HALF*(z3-z4)*(ONE+s(i)) + z4 )
                           zc(i,j) = HALF*(c2-c1)*(ONE+p(j)) + c1

                        ENDDO
                     ENDDO

                  ELSEIF( iSide == west .OR. iSide == east )THEN

                     nID = thisQuadMesh % sideMap(iSide)
                     DO j = 0, gPolyDeg
                        DO i = 0, gPolyDeg
                  
                           xc(i,j) = x2d(nID,i) 
                           yc(i,j) = y2d(nID,i)

                           c1 = -h(nID,i,iEl) !! Bottom curve of this face is the bathymetry 
                           c2 = ( HALF*(z3-z4)*(ONE+s(i)) + z4 )
                           zc(i,j) = HALF*(c2-c1)*(ONE+p(j)) + c1

                        ENDDO
                     ENDDO

                  ELSEIF( iSide == bottom )THEN

                     xc = x2d
                     yc = y2d
                     zc = -h(:,:,iEl)

                  ELSE ! Top of the bottom element
                    
                     xc = x2d
                     yc = y2d
                     DO j = 0, gPolyDeg
                        DO i = 0, gPolyDeg

                           c1 = ( HALF*(z2-z1)*(ONE+s(i)) + z1 ) 
                           c2 = ( HALF*(z3-z4)*(ONE+s(i)) + z4 )
                           zc(i,j) = HALF*(c2-c1)*(ONE+p(j)) + c1

                        ENDDO
                     ENDDO

                  ENDIF


               ELSE 

                  IF( iSide == south .OR. iSide == north )THEN  

                     nID = thisQuadMesh % sideMap(iSide)
                     
                     DO j = 0, gPolyDeg
                        DO i = 0, gPolyDeg
                  
                           xc(i,j) = x2d(i,nID) 
                           yc(i,j) = y2d(i,nID)

                           c1 = ( HALF*(z2-z1)*(ONE+s(i)) + z1 ) 
                           c2 = ( HALF*(z3-z4)*(ONE+s(i)) + z4 )
                           zc(i,j) = HALF*(c2-c1)*(ONE+p(j)) + c1

                        ENDDO
                     ENDDO

                  ELSEIF( iSide == west .OR. iSide == east )THEN

                     nID = thisQuadMesh % sideMap(iSide)
                     DO j = 0, gPolyDeg
                        DO i = 0, gPolyDeg
                  
                           xc(i,j) = x2d(nID,i) 
                           yc(i,j) = y2d(nID,i)

                           c1 = ( HALF*(z2-z1)*(ONE+s(i)) + z1 ) 
                           c2 = ( HALF*(z3-z4)*(ONE+s(i)) + z4 )
                           zc(i,j) = HALF*(c2-c1)*(ONE+p(j)) + c1

                        ENDDO
                     ENDDO

                  ELSE

                     xc = x2d
                     yc = y2d

                     DO j = 0, gPolyDeg
                        DO i = 0, gPolyDeg
                  
                           c1 = ( HALF*(z2-z1)*(ONE+s(i)) + z1 ) 
                           c2 = ( HALF*(z3-z4)*(ONE+s(i)) + z4 )
                           zc(i,j) = HALF*(c2-c1)*(ONE+p(j)) + c1

                        ENDDO
                     ENDDO

                  ENDIF

               ENDIF

               CALL boundSurfs(iSide) % SetNodes( xc, yc, zc ) 
            
            ENDDO

            CALL myHexMesh % elements(elID) % Build( nodes, elID, boundSurfs, interp )

         ENDDO
      ENDDO 

      PRINT*, 'Constructing face connectivity...'
      CALL myHexMesh % ConstructFaces( )
      nFaces = myHexMesh % nFaces
      PRINT*, 'nFaces (verification)    : ', nFaces
      

      ! Clear up memory
      DEALLOCATE( s, p, xc, yc, zc )

      DO iSide = 1, nHexFaces
         CALL boundSurfs(iSide) % Trash( )
      ENDDO
  
 END SUBROUTINE TerrainFollowingMesh_StackedHexMesh
!
!
!==================================================================================================!
!--------------------------------- Mesh File I/O Routines -----------------------------------------!
!==================================================================================================!
!
!
 
!
!
!
END MODULE StackedHexMeshClass
