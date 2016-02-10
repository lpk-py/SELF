MODULE CGSEM2D_FEM_PRECONDITIONER
!    
! September 2015
!
!    
!
!
!
!
! ========================================================================================== !

! The model settings
USE MODEL_PRECISION
USE CONSTANTS_DICTIONARY
USE MODEL_FLAGS
!
USE COMMONROUTINES
USE MATRIXROUTINES
USE RUN_PARAMS_CLASS

USE NODAL_STORAGE_2D_CLASS

USE QUADMESH_CLASS
USE EDGE_CLASS
USE CORNERNODE_2D_CLASS
USE GEOMETRY_BASICS
USE VECTOR_CLASS
USE MAPPEDGEOM_2D_CLASS



IMPLICIT NONE

      TYPE FEM_FLUX_COEFFICIENTS
      ! This structure holds the coefficients for calculatin the flux divergence
      !  div( F ) = ( a*dsdx + b*dsdy )_x + ( c*dsdx + d*dsdy )_y
         real(prec), allocatable :: a(:,:), b(:,:), c(:,:), d(:,:)
      
          CONTAINS
          PROCEDURE :: BUILD => BUILD_FEM_FLUX_COEFFS
          PROCEDURE :: TRASH => TRASH_FEM_FLUX_COEFFS

      END TYPE FEM_FLUX_COEFFICIENTS

      TYPE FEM2D_PRECONDITIONER
         TYPE( QUADMESH )                         :: mesh           ! Boundary condition flags are stored in the secondary element ID's of the edges.
         real(prec), allocatable                  :: solution(:,:,:)
         real(prec), allocatable                  :: source(:,:,:)     

         type(FEM_FLUX_COEFFICIENTS), allocatable :: fluxcoeffs(:)

         CONTAINS

         PROCEDURE :: BUILD => BUILD_CGSEM2D
         PROCEDURE :: TRASH => TRASH_CGSEM2D
       
         PROCEDURE :: MASK => MASK_FEM2D
         PROCEDURE :: UNMASK => UNMASK_FEM2D
         PROCEDURE :: GLOBAL_SUM => GLOBALSUM_FEM2D

         PROCEDURE :: MATRIX_ACTION => MATRIX_ACTION_FEM2D
         PROCEDURE :: RESIDUAL

         PROCEDURE :: VECTOR_PRODUCT
         PROCEDURE :: CONJUGATE_GRADIENT => CONJUGATE_GRADIENT_CGSEM2D


      END TYPE FEM2D_PRECONDITIONER


      

CONTAINS

 SUBROUTINE BUILD_CGSEM2D( myCGSEM )
 ! S/R BUILD_CGSEM2D
 !
 ! 
 ! ================================================ !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS( CGSEM2D ), intent(inout) :: myCGSEM
   ! LOCAL
   integer :: nS, nP, iS, iP, iEl
   character(40) :: meshFile
   character(10) :: pickupIterChar

      CALL myCGSEM % params % BUILD( )
      nS = myCGSEM % params % polyDeg
      nP = nS

      ! Build the nodal storage structure
      CALL myCGSEM % cgStorage % BUILD( nS, nP, GAUSS_LOBATTO, CG )


     ! Builds the lateral mesh
      if( TRIM(myCGSEM % params % ISMmeshFile) == ' ' )then

         meshFile = 'mesh'
         CALL myCGSEM % mesh % READ_PICKUP( trim(meshfile), myCGSEM % cgStorage  ) 

      else

         CALL myCGSEM % mesh % READ_SPECMESH_FILE( myCGSEM % params, myCGSEM % cgStorage )

      endif

      CALL myCGSEM % mesh % SCALE_THE_MESH( myCGSEM % params % xScale, &
                                            myCGSEM % params % yScale, &
                                            myCGSEM % cgStorage )

      ! Allocate space for the solution, source term, and diffusivity coefficient
      ALLOCATE( myCGSEM % boundaryFlux(1:myCGSEM % mesh % nElems, 1:4, 0:nS ), &
                myCGSEM % solution(0:nS, 0:nP, 1:myCGSEM % mesh % nElems), &
                myCGSEM % source(0:nS, 0:nP, 1:myCGSEM % mesh % nElems), &
                myCGSEM % fluxcoeffs(1:myCGSEM % mesh % nElems) )

      myCGSEM % solution =  ZERO
      myCGSEM % source = ZERO
      myCGSEM % boundaryFlux = ZERO

     do iEl = 1, myCGSEM % mesh % nElems
   
        CALL myCGSEM % fluxCoeffs(iEl) % BUILD( nS, nP )
 
     enddo


    !  if( myCGSEM % params % iterInit > 0 )then 
    !     write(pickupIterChar,'(I10.10)')  myCGSEM % params % iterInit
    !     CALL myCGSEM  % READ_PICKUP( 'Solution.'//pickupIterchar//'.pickup' )
    !  endif
 
   
      CALL myCGSEM % SET_DIRICHLET_BOUNDARY_CONDITION( DIRICHLET_FUNCTION )

      CALL myCGSEM % diag % BUILD( nP, myCGSEM % mesh % nElems )


 END SUBROUTINE BUILD_CGSEM2D
!
!
!
 SUBROUTINE TRASH_CGSEM2D( myCGSEM )
 ! S/R TRASH_CGSEM2D
 !
 ! 
 ! ================================================ !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS( CGSEM2D ), intent(inout) :: myCGSEM


      ! Trash the nodal storage structure
      CALL myCGSEM % cgStorage % TRASH( )

      ! Trash the geometry
      CALL myCGSEM % mesh % TRASH( )

      ! Trash the diagonal preconditioner
      CALL myCGSEM % diag % TRASH( )

      ! Deallocate space for the solution, source term, and diffusivity coefficient
      DEALLOCATE( myCGSEM % boundaryFlux, &
                  myCGSEM % solution,&
                  myCGSEM % source )

 END SUBROUTINE TRASH_CGSEM2D
!
!
!
 SUBROUTINE BUILD_FLUX_COEFFS(myCoeffs, nS, nP)
 ! S/R BUILD_FLUX_COEFFS
 !
 !   Allocates space for the flux coefficients
 !
 ! ============================================================================= !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS( FLUX_COEFFICIENTS ), intent(out) :: myCoeffs
   integer, intent(in)                     :: nS, nP


      ALLOCATE( myCoeffs % a(0:nS,0:nP), &
                myCoeffs % b(0:nS,0:nP), &
                myCoeffs % c(0:nS,0:nP), &
                myCoeffs % d(0:nS,0:nP) )

      myCoeffs % a = ZERO
      myCoeffs % b = ZERO
      myCoeffs % c = ZERO
      myCoeffs % d = ZERO

  END SUBROUTINE BUILD_FLUX_COEFFS
!
!
!
  SUBROUTINE TRASH_FLUX_COEFFS(myCoeffs)
 ! S/R TRASH_FLUX_COEFFS
 !
 !   Deallocates space for the flux coefficients
 !
 ! ============================================================================= !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS( FLUX_COEFFICIENTS ), intent(inout) :: myCoeffs


      DEALLOCATE( myCoeffs % a, &
                  myCoeffs % b, &
                  myCoeffs % c, &
                  myCoeffs % d )


  END SUBROUTINE TRASH_FLUX_COEFFS
!
!
!
!
!================================================================================!
!=======================   PROBLEM SPECIFIC ROUTINES    =========================!
!================================================================================!
!
!
!
 SUBROUTINE SET_DIAG_COEFFS( myCGSEM, polyDeg )
 ! S/R SET_DIAG_COEFFS
 !
 !   Sets the coefficients for the diagonal preconditioner.
 !   An impulse response function is used to obtain a matrix column of the laplacian
 !   operator. 
 !
 ! ============================================================================= !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS( CGSEM2D ), intent(inout) :: myCGSEM
   integer, intent(in)             :: polyDeg
   ! LOCAL
   integer :: iS, iP, iEl
   real(prec) :: IRF(0:polyDeg, 0:polyDeg, 1:myCGSEM % mesh % nElems)
   real(prec) :: aIRF(0:polyDeg, 0:polyDeg, 1:myCGSEM % mesh % nElems)
 

      do iEl = 1, myCGSEM % mesh % nElems
         do iP = 0, polyDeg
            do iS = 0, polyDeg

               ! Set the impulse response function
               IRF = ZERO
               IRF(iS,iP, iEl) = ONE
            
               ! Calculate the Laplacian of the impulse response function
               CALL myCGSEM % MATRIX_ACTION(IRF, aIRF )
            

               ! The diagonal contribution is the matrix action at iS,iP
               myCGSEM % diag % coeffs(iS,iP,iEl) = aIRF(iS,iP,iEl)
      
            enddo
         enddo
      enddo
 

      CALL myCGSEM % UNMASK( myCGSEM % diag % coeffs )

 END SUBROUTINE SET_DIAG_COEFFS
!
!
!================================================================================!
!=======================       SUPPORT ROUTINES       ===========================!
!================================================================================!
!
!
!
 FUNCTION DIRICHLET_FUNCTION( x, y ) RESULT( u )
 ! DIRICHLET_FUNCTION
 !
 ! ============================================================================= !
  IMPLICIT NONE
  real(prec) :: x, y
  real(prec) :: u

     u = ZERO

 END FUNCTION DIRICHLET_FUNCTION
!
!
!
 SUBROUTINE SET_DIRICHLET_BOUNDARY_CONDITION( myCGSEM, dirF )
 ! S/R SET_DIRICHLET_BOUNDARY_CONDITION
 !
 !
 ! ============================================================================= !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS( CGSEM2D ), intent(inout) :: myCGSEM
   ! dirF is a function which is passed to this subroutine
   real(prec), external :: dirF !(x, y)
   !LOCAL
   integer :: nS, nP, iEdge, eID, sID, iNode, nID, i, j, eID0, i0, j0, nEl
   real(prec) :: x, y 

      nS = myCGSEM % cgStorage % nS
      nP = myCGSEM % cgStorage % nP
      nEl = myCGSEM % mesh % nElems
   
      ! Mask out the solution over the edges -- excluding the corner nodes
      do iEdge = 1, myCGSEM % mesh % nEdges ! Loop over the edges
 
         if( myCGSEM % mesh % edges(iEdge) % elementIDS(2)  == DIRICHLET )then ! This is a dirichlet boundary edge

            eID = myCGSEM % mesh % edges(iEdge) % elementIDS(1)
            sID = myCGSEM % mesh % edges(iEdge) % elementSides(1)

            CALL SET_DIRICHLET_SIDE( eID, sID, myCGSEM % mesh, myCGSEM % solution, nS, nP, nEl, dirF )

         endif

      enddo ! iEdge, loop over the edges     


      ! Mask out the corner-nodes appropriately ( Leave the primary element intact )
      do iNode = 1, myCGSEM % mesh % nNodes ! Loop over the corner nodes

         if( myCGSEM % mesh % nodes(iNode) % nodeType == DIRICHLET ) then ! this is a prescribed (Dirichlet) boundary

            ! Rewind to the head of the linked list
            myCGSEM % mesh % nodes(iNode) % nodeToElement % current => myCGSEM % mesh % nodes(iNode) % nodeToElement % head

            ! Loop over the elements in this corner-node's connectivity list
            do while( ASSOCIATED( myCGSEM % mesh % nodes(iNode) % nodeToElement % current ) ) ! while there are elements in the node-to-element list

               CALL myCGSEM % mesh % nodes(iNode) % nodeToElement % GET_CURRENT_DATA( eID, nID ) ! Get the element ID and the local node ID (1->4)

               i = myCGSEM % mesh % cornerMap(1, nID) ! Grab the index over the first computational direction associated with this corner node
               j = myCGSEM % mesh % cornerMap(2, nID) ! Grab the index over the second computational direction associated with this corner node

               x = myCGSEM % mesh % elements(eID) % geometry % x(i,j)
               y = myCGSEM % mesh % elements(eID) % geometry % y(i,j)

               myCGSEM % solution(i,j,eID) = dirF( x, y)

               CALL myCGSEM % mesh % nodes(iNode) % nodeToElement % MOVE_TO_NEXT( )

            enddo ! while we have elements in the node-to-element list



         endif


      enddo


 END SUBROUTINE SET_DIRICHLET_BOUNDARY_CONDITION
!
!
!
 SUBROUTINE SET_DIRICHLET_SIDE( eID, sID, mesh, u, nS, nP, nElems, dirF )
 ! S/R SET_DIRICHLET_SIDE
 !
 !  This subroutine sets the solutions on side "sID" of element "eID",
 !  excluding the corner nodes ( Assumes Gauss-Lobatto points are used )
 !
 ! ======================================================================
 ! DECLARATIONS
   IMPLICIT NONE
   integer, intent(in)        :: eID, sID, nS, nP,nElems
   real(prec), intent(inout)  :: u(0:nS,0:nP,1:nElems)
   TYPE(QUADMESH), intent(in) :: mesh
   ! dirF is a function which is passed to this subroutine
   real(prec), external :: dirF !(x, y)
   ! LOCAL
   integer    :: iS, iP
   real(prec) :: x, y
   
      if( sID == 2 .OR. sID == 4)then ! east or west sides

         iS = mesh % sideMap(sID) 

         do iP = 1, nP-1 ! Loop over the edge, excluding the corner nodes

            x = mesh % elements(eID) % geometry % x(iS,iP)
            y = mesh % elements(eID) % geometry % y(iS,iP)

            u(iS,iP,eID) = dirF( x, y)

         enddo ! iP, loop over the edge

      else ! south or north sides

         iP = mesh % sideMap(sID)

         do iS = 1, nS-1 ! Loop over the edge, excluding the corner nodes

            x = mesh % elements(eID) % geometry % x(iS,iP)
            y = mesh % elements(eID) % geometry % y(iS,iP)

            u(iS,iP,eID) = dirF( x, y)

         enddo

      endif

 END SUBROUTINE SET_DIRICHLET_SIDE
!
!
!
 SUBROUTINE MASK_CGSEM2D( myCGSEM, u )
 ! S/R MASK_CGSEM2D
 !
 !
 ! ============================================================================= !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS( CGSEM2D ), intent(inout) :: myCGSEM
   real(prec), intent(inout)       :: u(0:myCGSEM % cgStorage % nS, &
                                        0:myCGSEM % cgStorage % nP, &
                                        1:myCGSEM % mesh % nElems )
   !LOCAL
   integer :: nS, nP, iEdge, eID, sID, iNode, nID, i, j, nEl

      nS = myCGSEM % cgStorage % nS
      nP = myCGSEM % cgStorage % nP
      nEl = myCGSEM % mesh % nElems
   
      ! Mask out the solution over the edges -- excluding the corner nodes
      do iEdge = 1, myCGSEM % mesh % nEdges ! Loop over the edges
 
         if( myCGSEM % mesh % edges(iEdge) % elementIDS(2) == DIRICHLET )then ! check the secondary element for boundary condition flag

            ! This is a prescribed boundary condition (Dirichlet) so we mask out the solution along this edge   

            eID = myCGSEM % mesh % edges(iEdge) % elementIDS(1)   ! Find the primary element ID

            sID = myCGSEM % mesh % edges(iEdge) % elementSides(1) ! Get the side of the element which should be masked

            CALL MASK_SIDE( eID, sID, myCGSEM % mesh, u, nS, nP, nEl )

         else ! then this is not a prescribed boundary, and we choose to mask out the secondary element

            eID = myCGSEM % mesh % edges(iEdge) % elementIDS(2)  ! Get the secondary element ID

            sID = abs(myCGSEM % mesh % edges(iEdge) % elementSides(2)) ! Get the side of the secondary element which should be masked

            if( eID > 0 )then ! this is an internal edge, and we should mask out the secondary element
               CALL MASK_SIDE( eID, sID, myCGSEM % mesh, u, nS, nP, nEl )
            endif

         endif


      enddo ! iEdge, loop over the edges     

      ! At this point, the secondary internal edges and prescribed boundary condition edges have been masked out.
      ! Now we mask out all but the primary element in the corner-node lists.

      ! Mask out the corner-nodes appropriately ( Leave the primary element intact )
      do iNode = 1, myCGSEM % mesh % nNodes ! Loop over the corner nodes

         ! Rewind to the head of the linked list
         myCGSEM % mesh % nodes(iNode) % nodeToElement % current => myCGSEM % mesh % nodes(iNode) % nodeToElement % head

         if( myCGSEM % mesh % nodes(iNode) % nodeType /= DIRICHLET ) then ! this is NOT a prescribed (Dirichlet) boundary
                                                                          ! The data in this non-prescribed node is preserved

            CALL myCGSEM % mesh % nodes(iNode) % nodeToElement % MOVE_TO_NEXT( ) ! skip the first element in the list so it's data is not masked

         endif

         do while( ASSOCIATED( myCGSEM % mesh % nodes(iNode) % nodeToElement % current ) ) ! while there are elements in the node-to-element list

            CALL myCGSEM % mesh % nodes(iNode) % nodeToElement % GET_CURRENT_DATA( eID, nID ) ! Get the element ID and the local node ID (1->4)

            i = myCGSEM % mesh % cornerMap(1, nID) ! Grab the index over the first computational direction associated with this corner node
            j = myCGSEM % mesh % cornerMap(2, nID) ! Grab the index over the second computational direction associated with this corner node
 
            u(i,j,eID) = ZERO ! Mask out the solution at this element at this corner node

            CALL myCGSEM % mesh % nodes(iNode) % nodeToElement % MOVE_TO_NEXT( )

         enddo ! while we have elements in the node-to-element list

      enddo


 END SUBROUTINE MASK_CGSEM2D
!
!
!
 SUBROUTINE MASK_SIDE( eID, sID, mesh, u, nS, nP, nElems )
 ! S/R MASK_SIDE
 !
 !  This subroutine masks the solutions on side "sID" of element "eID",
 !  excluding the corner nodes ( Assumes Gauss-Lobatto points are used )
 !
 ! ======================================================================
 ! DECLARATIONS
   IMPLICIT NONE
   integer, intent(in)        :: eID, sID, nS, nP,nElems
   real(prec), intent(inout)  :: u(0:nS,0:nP,1:nElems)
   TYPE(QUADMESH), intent(in) :: mesh
   ! LOCAL
   integer :: iS, iP
   
      if( sID == 2 .OR. sID == 4)then ! east or west sides

         iS = mesh % sideMap(sID) 

         do iP = 1, nP-1 ! Loop over the edge, excluding the corner nodes

            u(iS,iP,eID) = ZERO

         enddo ! iP, loop over the edge

      else ! south or north sides

         iP = mesh % sideMap(sID)

         do iS = 1, nS-1 ! Loop over the edge, excluding the corner nodes

            u(iS,iP,eID) = ZERO

         enddo

      endif

 END SUBROUTINE MASK_SIDE
!
!
!
 SUBROUTINE UNMASK_CGSEM2D( myCGSEM, u )
 ! S/R UNMASK_CGSEM2D
 !
 !
 ! ============================================================================= !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS( CGSEM2D ), intent(inout) :: myCGSEM
   real(prec), intent(inout)       :: u(0:myCGSEM % cgStorage % nS, &
                                        0:myCGSEM % cgStorage % nP, &
                                        1:myCGSEM % mesh % nElems )
   !LOCAL
   integer :: nS, nP, iEdge, eID, sID, iNode, nID, i, j, eID0, i0, j0, nEl

      nS = myCGSEM % cgStorage % nS
      nP = myCGSEM % cgStorage % nP
      nEl = myCGSEM % mesh % nElems
   
      ! Unmask the solution over the interior edges only -- corner nodes are excluded in this step
      ! and are taken care of in the following step
      do iEdge = 1, myCGSEM % mesh % nEdges ! Loop over the edges
 
         eID = myCGSEM % mesh % edges(iEdge) % elementIDS(2)

         if( eID > 0 )then ! This is an interior shared edge

            CALL UNMASK_SIDE( myCGSEM % mesh % edges(iEdge), myCGSEM % mesh, u, nS, nP, nEl ) ! unmask the solution along this edge

         endif

      enddo ! iEdge, loop over the edges     


      ! Unmask the corner-nodes
      do iNode = 1, myCGSEM % mesh % nNodes ! Loop over the corner nodes


         if( myCGSEM % mesh % nodes(iNode) % nodeType /= DIRICHLET ) then ! this is NOT a prescribed (Dirichlet) boundary

            ! Rewind to the head of the linked list
            myCGSEM % mesh % nodes(iNode) % nodeToElement % current => myCGSEM % mesh % nodes(iNode) % nodeToElement % head

            ! Gather the element and local node IDs of the corner which contains the unmasked solution values.
            CALL myCGSEM % mesh % nodes(iNode) % nodeToElement % GET_CURRENT_DATA( eID0, nID ) ! Get the element ID and the local node ID (1->4)

            i0 = myCGSEM % mesh % cornerMap(1, nID) ! Grab the index over the first computational direction associated with this corner node
            j0 = myCGSEM % mesh % cornerMap(2, nID) ! Grab the index over the second computational direction associated with this corner node

            CALL myCGSEM % mesh % nodes(iNode) % nodeToElement % MOVE_TO_NEXT( ) ! move on to the next element which shares this node


            ! Loop over the elements in this corner-node's connectivity list
            do while( ASSOCIATED( myCGSEM % mesh % nodes(iNode) % nodeToElement % current ) ) ! while there are elements in the node-to-element list

               CALL myCGSEM % mesh % nodes(iNode) % nodeToElement % GET_CURRENT_DATA( eID, nID ) ! Get the element ID and the local node ID (1->4)

               i = myCGSEM % mesh % cornerMap(1, nID) ! Grab the index over the first computational direction associated with this corner node
               j = myCGSEM % mesh % cornerMap(2, nID) ! Grab the index over the second computational direction associated with this corner node
  
               u(i,j,eID) = u(i0,j0,eID0) ! copy the solution from the unmasked solution

               CALL myCGSEM % mesh % nodes(iNode) % nodeToElement % MOVE_TO_NEXT( )

            enddo ! while we have elements in the node-to-element list



         endif


      enddo


 END SUBROUTINE UNMASK_CGSEM2D
!
!
!
 SUBROUTINE UNMASK_SIDE( thisEdge, mesh, u, nS, nP, nElems )
 ! S/R UNMASK_SIDE
 !
 !  This subroutine unmasks the solutions on side "sID" of element "eID",
 !  excluding the corner nodes ( Assumes Gauss-Lobatto points are used )
 !
 ! ======================================================================
 ! DECLARATIONS
   IMPLICIT NONE
   TYPE(EDGE), intent(in)     :: thisEdge
   integer, intent(in)        :: nS, nP,nElems
   real(prec), intent(inout)  :: u(0:nS,0:nP,1:nElems)
   TYPE(QUADMESH), intent(in) :: mesh
   ! LOCAL
   integer :: iS, iP, eID, sID, jS, jP
   real(prec) :: temp(0:nS) ! assume nS == nP
   

      ! Get the primary element and side IDs
      eID = thisEdge % elementIDs(1)
      sID = thisEdge % elementSides(1)
   
 
      if( sID == 2 .OR. sID == 4)then ! east or west sides

         iS = mesh % sideMap(sID) 

         do iP = 1, nP-1 ! Loop over the edge, excluding the corner nodes

            temp(iP) = u(iS,iP,eID) ! Copy the solution in the primary element

         enddo ! iP, loop over the edge

      else ! south or north sides

         iP = mesh % sideMap(sID)

         do iS = 1, nS-1 ! Loop over the edge, excluding the corner nodes

            temp(iS) = u(iS,iP,eID) ! Copy the solution in the primary element

         enddo

      endif

       ! Get the secondary element and side IDs
      eID = thisEdge % elementIDs(2)
      sID = abs( thisEdge % elementSides(2) )
        
      if( sID == 2 .OR. sID == 4)then ! east or west sides

         iS = mesh % sideMap(sID) 

         jP = thisEdge % start ! accounting for the possibility that the node ordering may be different in the secondary element

         do iP = 1, nP-1 ! Loop over the edge, excluding the corner nodes

            u(iS,jP,eID) = temp(iP) ! Copy the solution to the secondary element

            jP = jP + thisEdge % inc ! increment jP

         enddo ! iP, loop over the edge

      else ! south or north sides

         iP = mesh % sideMap(sID)
       
         jS = thisEdge % start ! accounting for the possibility that the node ordering may be different in the secondary element

         do iS = 1, nS-1 ! Loop over the edge, excluding the corner nodes

            u(jS,iP,eID) = temp(iS) ! Copy the solution in the primary element

            jS = jS + thisEdge % inc ! increment jS

         enddo

      endif

 END SUBROUTINE UNMASK_SIDE
!
!
!
 SUBROUTINE GLOBALSUM_CGSEM2D( myCGSEM, u )
 ! S/R GLOBALSUM_CGSEM2D
 !
 !  Adds together the shared edge and node contributions for the CGSEM method.
 !
 ! ============================================================================= !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS( CGSEM2D ), intent(inout) :: myCGSEM
   real(prec), intent(inout)       :: u(0:myCGSEM % cgStorage % nS, &
                                        0:myCGSEM % cgStorage % nP, &
                                        1:myCGSEM % mesh % nElems )
   !LOCAL
   integer :: iEdge, eID, sID, iNode, nID, i, j, eID0, i0, j0
   real(prec) :: theSum


      do iEdge = 1, myCGSEM % mesh % nEdges ! loop over the edges

         if( myCGSEM % mesh % edges(iEdge)  % elementIDS(2) > 0 )then !/= DIRICHLET )then ! This is an interior shared edge (or a neumann boundary edge)

            CALL SUM_SIDE( myCGSEM % mesh % edges(iEdge), myCGSEM % mesh, u, &
                           myCGSEM % cgStorage % nS, myCGSEM % cgStorage % nP, &
                           myCGSEM % mesh % nElems )                             ! add the contributions from the elements which share this edge
                                                                                 ! The edge sum excludes corner points

         endif

      enddo ! iEdge, loop over the edges


      do iNode = 1, myCGSEM % mesh % nNodes ! loop over the corner nodes

         if( myCGSEM % mesh % nodes(iNode) % nodeType /= DIRICHLET ) then ! this node does NOT lie on a prescribed (Dirichlet) boundary

             ! Rewind to the head of the linked list
            myCGSEM % mesh % nodes(iNode) % nodeToElement % current => myCGSEM % mesh % nodes(iNode) % nodeToElement % head

            ! *** First, compute the sum from the contributing elements *** !

            theSum = ZERO 


            ! Loop over the elements in this corner-node's connectivity list
            do while( ASSOCIATED( myCGSEM % mesh % nodes(iNode) % nodeToElement % current ) ) ! while there are elements in the node-to-element list

               CALL myCGSEM % mesh % nodes(iNode) % nodeToElement % GET_CURRENT_DATA( eID, nID ) ! Get the element ID and the local node ID (1->4)

               i = myCGSEM % mesh % cornerMap(1, nID) ! Grab the index over the first computational direction associated with this corner node
               j = myCGSEM % mesh % cornerMap(2, nID) ! Grab the index over the second computational direction associated with this corner node
  
               theSum = theSum + u(i,j,eID) 

               CALL myCGSEM % mesh % nodes(iNode) % nodeToElement % MOVE_TO_NEXT( )

            enddo ! while we have elements in the node-to-element list


            ! *** Now, copy the sum to the array "u" for each element which contributed to the sum *** !

             ! Rewind to the head of the linked list
            myCGSEM % mesh % nodes(iNode) % nodeToElement % current => myCGSEM % mesh % nodes(iNode) % nodeToElement % head

            do while( ASSOCIATED( myCGSEM % mesh % nodes(iNode) % nodeToElement % current ) ) ! while there are elements in the node-to-element list

               CALL myCGSEM % mesh % nodes(iNode) % nodeToElement % GET_CURRENT_DATA( eID, nID ) ! Get the element ID and the local node ID (1->4)

               i = myCGSEM % mesh % cornerMap(1, nID) ! Grab the index over the first computational direction associated with this corner node
               j = myCGSEM % mesh % cornerMap(2, nID) ! Grab the index over the second computational direction associated with this corner node
  
               u(i,j,eID) = theSum

               CALL myCGSEM % mesh % nodes(iNode) % nodeToElement % MOVE_TO_NEXT( )

            enddo ! while we have elements in the node-to-element list



         endif

      enddo

 END SUBROUTINE GLOBALSUM_CGSEM2D
!
!
!
 SUBROUTINE SUM_SIDE( thisEdge, mesh, u, nS, nP, nElems )
 ! S/R SUM_SIDE
 !
 ! 
 !
 ! ======================================================================
 ! DECLARATIONS
   IMPLICIT NONE
   TYPE(EDGE), intent(in)     :: thisEdge
   integer, intent(in)        :: nS, nP,nElems
   real(prec), intent(inout)  :: u(0:nS,0:nP,1:nElems)
   TYPE(QUADMESH), intent(in) :: mesh
   ! LOCAL
   integer :: iS, iP, eID, sID, jS, jP, k, n
   real(prec) :: temp(0:nS,1:2) ! assume nS == nP
   real(prec) :: theSum


      do k = 1,2

         ! Get the primary/secondary element and side IDs
         eID = thisEdge % elementIDs(k)
         sID = abs( thisEdge % elementSides(k) )
   
 
         if( sID == 2 .OR. sID == 4)then ! east or west sides

            iS = mesh % sideMap(sID) 

            do iP = 1, nP-1 ! Loop over the edge, excluding the corner nodes

               temp(iP,k) = u(iS,iP,eID) ! 

            enddo ! iP, loop over the edge

         else ! south or north sides

            iP = mesh % sideMap(sID)

            do iS = 1, nS-1 ! Loop over the edge, excluding the corner nodes

               temp(iS,k) = u(iS,iP,eID) ! 

            enddo

         endif

      enddo

      ! Add the two contributions together
      n = thisEdge % start
     
      do iS = 1, nS-1

         theSum = temp(iS,1) + temp(n,2)

         temp(iS,1) = theSum

         temp(n,2) = theSum

         n = n + thisEdge % inc ! increment the secondary element array adress

      enddo

      ! copy the sum into the two contributing element edges

      do k = 1,2

         ! Get the primary/secondary element and side IDs
         eID = thisEdge % elementIDs(k)
         sID = abs( thisEdge % elementSides(k) )
   
 
         if( sID == 2 .OR. sID == 4)then ! east or west sides

            iS = mesh % sideMap(sID) 

            do iP = 1, nP-1 ! Loop over the edge, excluding the corner nodes

                u(iS,iP,eID) = temp(iP,k) !

            enddo ! iP, loop over the edge

         else ! south or north sides

            iP = mesh % sideMap(sID)

            do iS = 1, nS-1 ! Loop over the edge, excluding the corner nodes

               u(iS,iP,eID) = temp(iS,k) ! 

            enddo

         endif

      enddo


 END SUBROUTINE SUM_SIDE
!
!
!
 SUBROUTINE CALC_GRADIENT( geometry, cgStorage, u, dudx, dudy, nS, nP )
 ! S/R CALC_GRADIENT
 !
 !  Calculates the gradient of the solution in computational coordinates within 
 !  a single element
 !
 !
 ! ============================================================================= !
 ! DECLARATIONS
   IMPLICIT NONE
   integer, intent(in)                  :: nS, nP
   TYPE( MAPPEDGEOM_2D ), intent(in)    :: geometry
   TYPE( NODAL_STORAGE_2D ), intent(in) :: cgStorage
   real(prec), intent(in)               :: u(0:nS, 0:nP)
   real(prec), intent(out)              :: dudx(0:nS, 0:nP)
   real(prec), intent(out)              :: dudy(0:nS, 0:nP) 
   ! LOCAL
   integer :: iS, iP
   real(prec) :: duds(0:nS, 0:nP)
   real(prec) :: dudp(0:nS, 0:nP)
   real(prec) :: J, dxds, dxdp, dyds, dydp


      do iP = 0, nP ! Loop over the second computational direction

         ! Calculate the derivative in the first computational direction
         duds(0:nS,iP) = MATVECMUL( cgStorage % dMatX(0:nS,0:nS), &
                                    u(0:nS,iP), nS, nS )

      enddo ! iP, loop over the second computational direction



      do iS = 0, nS ! Loop over the first computational direction

         ! Calculate the derivative in the second computational direction
         dudp(iS,0:nP) = MATVECMUL( cgStorage % dMatY(0:nP,0:nP), &
                                    u(iS,0:nP), nP, nP )

      enddo ! iS, loop over the first computational direction


      ! Calculate the gradient in physical space
      do iP = 0, nP
         do iS = 0, nS

            J = geometry % J(iS,iP)

            dxds = geometry % dxds(iS,iP)
            dxdp = geometry % dxdp(iS,iP)
            dyds = geometry % dyds(iS,iP)
            dydp = geometry % dydp(iS,iP)

            dudx(iS,iP) = ( dydp*duds(iS,iP) - dyds*dudp(iS,iP) )/J

            dudy(iS,iP) = ( dxds*dudp(iS,iP) - dxdp*duds(iS,iP) )/J        

         enddo
      enddo


 END SUBROUTINE CALC_GRADIENT
!
!
!
 FUNCTION FLUX_X( u, dudx, dudy, a, b) RESULT( fx )
 ! FUNCTION FLUX_X
 !
 !   In the problem Div( F ) = s, the x-component of the flux is calculated here.
 !
 !
 ! ============================================================================= !
 ! DECLARATIONS
   IMPLICIT NONE
   real(prec) :: u, dudx, dudy, fx, a, b


      fx = a*dudx + b*dudy


 END FUNCTION FLUX_X
!
!
!
  FUNCTION FLUX_Y( u, dudx, dudy, c, d ) RESULT( fy )
 ! FUNCTION FLUX_Y
 !
 !   In the problem Div( F ) = s, the y-component of the flux is calculated here.
 !
 !
 ! ============================================================================= !
 ! DECLARATIONS
   IMPLICIT NONE
   real(prec) :: u, dudx, dudy, fy, c, d


      fy = c*dudx + d*dudy


 END FUNCTION FLUX_Y
!
!
!
 SUBROUTINE FLUX_DIVERGENCE( geometry, cgStorage, coeffs, u, Lu, nS, nP )
 ! S/R FLUX_DIVERGENCE
 !
 !
 ! ============================================================================= !
 ! DECLARATIONS
   IMPLICIT NONE
   integer, intent(in)                  :: nS, nP
   TYPE( MAPPEDGEOM_2D ), intent(in)    :: geometry
   TYPE( NODAL_STORAGE_2D ), intent(in) :: cgStorage
   TYPE( FLUX_COEFFICIENTS), intent(in) :: coeffs
   real(prec), intent(in)               :: u(0:nS, 0:nP)
   real(prec), intent(out)              :: Lu(0:nS, 0:nP)
   ! LOCAL
   integer    :: iS, iP
   real(prec) :: F1(0:nS,0:nP)
   real(prec) :: F2(0:nS,0:nP)
   real(prec) :: dF1ds(0:nS, 0:nP)
   real(prec) :: dF2dp(0:nS, 0:nP)
   real(prec) :: dudx(0:nS, 0:nP)
   real(prec) :: dudy(0:nS, 0:nP)
   real(prec) :: J, dxds, dxdp, dyds, dydp, Fx, Fy


      CALL CALC_GRADIENT( geometry, cgStorage, u, dudx, dudy, nS, nP ) ! calculate the gradient in physical coordinates

      ! Now apply the metric terms
      do iP = 0, nP ! Loop over the second computational direction
         do iS = 0, nS ! Loop over the first computational direction

            J = geometry % J(iS,iP)

            dxds = geometry % dxds(iS,iP)
            dxdp = geometry % dxdp(iS,iP)
            dyds = geometry % dyds(iS,iP)
            dydp = geometry % dydp(iS,iP)

            ! Calculate the conservative flux
            Fx = FLUX_X( u(iS,iP), dudx(iS,iP), dudy(iS,iP), coeffs % a(iS,iP), coeffs % b(iS,iP) )
            Fy = FLUX_Y( u(iS,iP), dudx(iS,iP), dudy(iS,iP), coeffs % c(iS,iP), coeffs % d(iS,iP) )

            ! Contravariant flux calculation
            F1(iS,iP) = ( Fx*dydp - Fy*dxdp )*cgStorage % qWeightX(iS)

            F2(iS,iP) = ( Fy*dxds - Fx*dyds )*cgStorage % qWeightY(iP)
        

         enddo ! iS
      enddo ! iP


   ! Now calculate the divergence of the flux

   do iP = 0, nP ! loop over the second computational direction

      dF1ds(0:nS,iP) = MATVECMUL_TRANSPOSE( cgStorage % dMatX(0:nS,0:nS), F1(0:nS,iP), nS )  

   enddo ! iP

   do iS = 0, nS ! loop over the second computational direction

      dF2dp(iS,0:nP) = MATVECMUL_TRANSPOSE( cgStorage % dMatY(0:nP,0:nP), F2(iS,0:nP), nP )  

   enddo ! iP

   ! Now calculate the Laplacian

   do iP = 0, nP
      do iS = 0, nS

         Lu(iS,iP) = -( dF1ds(iS,iP)*cgStorage % qWeightY(iP) +&
                        dF2dp(iS,iP)*cgStorage % qWeightX(iS) )

      enddo ! iS
   enddo ! iP
       

 END SUBROUTINE FLUX_DIVERGENCE 
!
!
!
 SUBROUTINE MATRIX_ACTION( myCGSEM, u, Au )
 !
 ! 
 ! ==========================================================!
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS( CGSEM2D), intent(inout) :: myCGSEM
   real(prec), intent(inout)       :: u(0:myCGSEM % cgStorage % nS, &
                                        0:myCGSEM % cgStorage % nP, &
                                        1:myCGSEM % mesh % nElems )
   real(prec), intent(out)         :: Au(0:myCGSEM % cgStorage % nS, &
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
   CLASS( CGSEM2D ), intent(inout) :: myCGSEM
   real(prec), external            :: dirF
   real(prec), intent(out)         :: r(0:myCGSEM % cgStorage % nS, &
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
   CLASS( CGSEM2D ), intent(inout)  :: myCGSEM
   real(prec), intent(inout)        :: u(0:myCGSEM % cgStorage % nS, &
                                         0:myCGSEM % cgStorage % nP, &
                                         1:myCGSEM % mesh % nElems )
   real(prec), intent(inout)        :: v(0:myCGSEM % cgStorage % nS, &
                                         0:myCGSEM % cgStorage % nP, &
                                         1:myCGSEM % mesh % nElems )
   real(prec), intent(out)          :: uDotv
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
 SUBROUTINE CONJUGATE_GRADIENT_CGSEM2D( myCGSEM, dirF, resi )
 !
 ! Inverts the system implied by the routine matrix action and residual 
 ! using the conjugate gradient algorithm.
 !
 ! =============================================================================!
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS( CGSEM2D ), intent(inout)    :: myCGSEM
   real(prec), external               :: dirF
   real(prec), intent(out), optional :: resi(0:myCGSEM % params % cgIterMax)
   ! LOCAL
   integer    :: nS, nP, iEl, nEl
   integer    :: iter
   integer    :: nIt
   real(prec) :: TOL
   real(prec) :: r(0:myCGSEM % cgStorage % nS, &
                   0:myCGSEM % cgStorage % nP, &
                   1:myCGSEM % mesh % nElems )
   real(prec) :: v(0:myCGSEM % cgStorage % nS, &
                   0:myCGSEM % cgStorage % nP, &
                   1:myCGSEM % mesh % nElems )
   real(prec) :: z(0:myCGSEM % cgStorage % nS, &
                   0:myCGSEM % cgStorage % nP, &
                   1:myCGSEM % mesh % nElems )
   real(prec) :: Ad(0:myCGSEM % cgStorage % nS, &
                    0:myCGSEM % cgStorage % nP, &
                    1:myCGSEM % mesh % nElems )
   real(prec) :: rNorm
   real(prec) :: a, b, num, den, r0
   
   

      nS = myCGSEM % cgStorage % nS 
      nP = myCGSEM % cgStorage % nP
      nEl = myCGSEM % mesh % nElems
      
      nIt = myCGSEM % params % cgIterMax
      TOL = myCGSEM % params % cgTol

      ! Calculate the residual
      CALL myCGSEM % RESIDUAL( dirF, r )
      
      CALL myCGSEM % VECTOR_PRODUCT( r, r, r0 )

      !print*, r0
      if( sqrt(r0) < TOL ) then
      !      print*, 'MODULE CGSEM2D_CLASS : CONJUGATE_GRADIENT : Initial residual small enough. '
            CALL myCGSEM % UNMASK( myCGSEM % solution )
            RETURN
      endif   
      if( PRESENT( resi ) )then
         resi = ZERO
         resi(0) = r0
      endif

      ! Apply the preconditioner
      CALL myCGSEM % diag % SOLVE( r, z )  ! 
     
      CALL myCGSEM % VECTOR_PRODUCT( r, z, num ) ! numerator   r (DOT) (H^(-1)r )

      v = z ! Copy 

      do iter = 1,nIt ! Loop over the PCG iterates
 
         ! Compute Ad matrix-vector productz
         CALL myCGSEM % MATRIX_ACTION( v, z ) ! z = A*v

         ! Compute the search-direction magnitude
         CALL myCGSEM % VECTOR_PRODUCT( v, z, den) ! denominator

         a = num/den

         ! Update the solution guess
         myCGSEM % solution(0:nS,0:nP,1:nEl) = myCGSEM % solution(0:nS,0:nP,1:nEl) + a*v(0:nS, 0:nP,1:nEl)

          
         ! Update the residual
         r(0:nS,0:nP,1:nEl) = r(0:nS,0:nP,1:nEl) - a*z(0:nS,0:nP,1:nEl)
         
         CALL myCGSEM % VECTOR_PRODUCT( r, r, rNorm ) 

         if( PRESENT( resi ) )then
            resi(iter) = sqrt(rNorm)
         endif

        ! print*, rNorm
         if( sqrt(rNorm)/r0 < TOL ) then
            EXIT
         endif   

         ! Apply the preconidionter to the residual
         CALL myCGSEM % diag % SOLVE( r, z )! solves H z = r, to get z = H^(-1) r

         den = num ! r(DOT)[ H^(-1)r ] = r(DOT)d

         ! Calculate the change in the search direction
         CALL myCGSEM % VECTOR_PRODUCT( r, z, num ) 
         
         ! Update the search direction
         b = num/den

         v(0:nS,0:nP,1:nEl) = z(0:nS,0:nP,1:nEl) + b*v(0:nS,0:nP,1:nEl)
         
      enddo ! iter, loop over the PCG iterates 


      if( sqrt(rNorm) > TOL ) then
         print*, 'MODULE CGSEM2D_CLASS : CONJUGATE_GRADIENT failed to converge '
         print*, 'Last L-2 residual : ', sqrt(rNorm)

      endif   


      CALL myCGSEM % UNMASK( myCGSEM % solution )

 END SUBROUTINE CONJUGATE_GRADIENT_CGSEM2D
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
  CLASS( CGSEM2D ), intent(in)           :: myCGSEM
!  integer, intent(in)                   :: nOld, nPlot
!  real(prec), intent(in)                :: Tmat(0:nPlot, 0:nOld)
!  TYPE( LAG_INTERP2D ), intent(in)      :: plotInterp
!  character(*), intent(in)              :: filename
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

END MODULE CGSEM2D_CLASS
