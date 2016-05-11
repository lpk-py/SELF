!  CGsemCoarsePreconditioner_Class.f90
!  
!  Copyright 2015 Joseph Schoonover <schoonover.numerics@gmail.com>
!  
!  CGsemCoarsePreconditioner_Class.f90 is part of the Spectral Element Libraries in Fortran (SELF).
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


! ========================================= Logs ================================================= !
! yyyy-mm-dd  Joe  <joe@clay>
!    This module defines the data structure used for implementing the Continuous Galerkin
!    Spectral Element Method. 
!
!    The module  is set up so that we solve div( Flux ) = Source
!   
! 
! //////////////////////////////////////////////////////////////////////////////////////////////// !


MODULE CGsemCoarsePreconditioner_Class

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
!/src/filters
USE ModalCutoffFilter2D_Class
! src/cgsem/
USE CGsemParams_Class


IMPLICIT NONE


      ! The Preconditioner class implements a 2nd order CGSEM where the preconditioner mesh is 
      ! generated from the higher order SE mesh. A mapping between the linear mesh and the 
      ! original mesh is generated. By reducing the polynomial order and increasing the number
      ! of elements, the CG elliptic operator becomes sparse and the condition number is reduced.
      TYPE CGsemPreconditioner_2D
         INTEGER                       :: maxIters, nS, nP
         REAL(prec)                    :: tolerance
         TYPE( NodalStorage_2D )       :: SpectralOps
         TYPE( QuadMesh )              :: mesh         
         REAL(prec), ALLOCATABLE       :: fluxCoeff(:,:,:)
         TYPE( ModalCutoffFilter2D )   :: filter
!         REAL(prec), ALLOCATABLE       :: lowToHighS(:,:), lowToHighP(:,:)
!         REAL(prec), ALLOCATABLE       :: highToLowS(:,:), highToLowP(:,:)


         CONTAINS

         PROCEDURE :: Build         => Build_CGsemPreconditioner_2D
         PROCEDURE :: BuildQuadMesh => BuildQuadMesh_CGsemPreconditioner_2D
         PROCEDURE :: Trash         => Trash_CGsemPreconditioner_2D
         PROCEDURE :: ConstructPCMatrix
         
         PROCEDURE :: MapFromOrigToPC => MapFromOrigToPC_Preconditioner_2D
         PROCEDURE :: MapFromPCToOrig => MapFromPCToOrig_Preconditioner_2D

         PROCEDURE :: Mask              => Mask_CGsemPreconditioner_2D
         PROCEDURE :: UnMask            => UnMask_CGsemPreconditioner_2D
         PROCEDURE :: GlobalSum         => GlobalSum_CGsemPreconditioner_2D
         PROCEDURE :: FluxDivergence    => FluxDivergence_CGsemPreconditioner_2D
         PROCEDURE :: CalculateGradient => CalculateGradient_CGsemPreconditioner_2D
         PROCEDURE :: MatrixAction      => MatrixAction_CGsemPreconditioner_2D
         PROCEDURE :: DotProduct        => DotProduct_CGsemPreconditioner_2D
         PROCEDURE :: Solve             => Solve_CGsemPreconditioner_2D 

      END TYPE CGsemPreconditioner_2D

      

CONTAINS

 SUBROUTINE Build_CGsemPreconditioner_2D( myPC, highOps, params )
 ! S/R Build
 !
 ! 
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS( CGsemPreconditioner_2D ), INTENT(inout) :: myPC
   TYPE( NodalStorage_2D ), INTENT(in)            :: highOps
   TYPE( CGsemParams ), INTENT(in)                :: params
   ! LOCAL
   INTEGER :: nS, nP, nID, iNode, highPdeg
 !  REAL(prec), ALLOCATABLE :: lowS(:), lowP(:), highS(:), highP(:)

      nS = params % polyDeg
      nP = nS

      highPDeg = params % polyDeg
      
      myPC % tolerance = params % pcTolerance
      myPC % maxIters  = params % pcIterates
      myPC % nS        = nS
      myPC % nP        = nP 
     
      CALL myPC % SpectralOps % Build( nS, nP, GAUSS_LOBATTO, CG )
      CALL myPC % BuildQuadMesh( params )

      DO iNode = 1, myPC % mesh % nNodes
         CALL myPC % mesh % GetNodeType( iNode, nID )
         IF( nID == NO_NORMAL_FLOW )THEN
            CALL myPC % mesh % SetNodeType( iNode, DIRICHLET )
         ENDIF
      ENDDO
      
      ! Allocate space for the solution, source term, and diffusivity coefficient
      ALLOCATE( myPC % fluxCoeff(0:nS, 0:nP, 1:myPC % mesh % nElems) )
      myPC % fluxCoeff = ONE

      CALL myPC % filter % Build( myPC % SpectralOps, params % pcDegree, params % pcDegree )
!      ALLOCATE( lowS(0:nS), lowP(0:nP), highS(0:highPDeg), highP(0:highPDeg) )

!      CALL highOps % interp % GetNodes( highS, highP )
!      CALL myPC % SpectralOps % interp % GetNodes( lowS, lowP )

!      ALLOCATE( myPC % lowToHighS(0:highPDeg,0:nS), myPC % lowToHighP(0:highPDeg,0:nP) )
!      ALLOCATE( myPC % highToLowS(0:nS,0:highPDeg), myPC % highToLowP(0:nP,0:highPDeg) )
      
!      CALL myPC % SpectralOps % interp % CalculateInterpolationMatrix( highPDeg, highPDeg, &
!                                                                       highS, highP, &
!                                                                       myPC % lowToHighS, &
!                                                                       myPC % lowToHighP )
!      CALL highOps % interp % CalculateInterpolationMatrix( nS, nP, &
!                                                            lowS, lowP, &
!                                                            myPC % highToLowS, &
!                                                            myPC % highToLowP ) 

!      DEALLOCATE( lowS, lowP, highS, highP )

 END SUBROUTINE Build_CGsemPreconditioner_2D
!
!
!
 SUBROUTINE BuildQuadMesh_CGsemPreconditioner_2D( myCGSEM, params )
 ! S/R BuildQuadMesh
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS( CGsemPreconditioner_2D ), INTENT(inout) :: myCGSEM
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

 END SUBROUTINE BuildQuadMesh_CGsemPreconditioner_2D
!
!
!
 SUBROUTINE Trash_CGsemPreconditioner_2D( myPC )
 ! S/R Trash
 !
 ! 
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS( CGsemPreconditioner_2D ), INTENT(inout) :: myPC


      ! Trash the nodal storage structure
      CALL myPC % SpectralOps % Trash( )

      ! Trash the geometry
      CALL myPC % mesh % Trash( )

      ! Deallocate space for the solution, source term, and diffusivity coefficient
      DEALLOCATE( myPC % fluxCoeff )!, &
!                  myPC % lowToHighS, &
!                  myPC % lowToHighP, &
!                  myPC % highToLowS, & 
!                  myPC % highToLowP )
                  
      CALL myPC % filter % Trash( )

 END SUBROUTINE Trash_CGsemPreconditioner_2D
!
!
!==================================================================================================!
!-------------------------------------- Type Specific ---------------------------------------------!
!==================================================================================================!
!
!
 SUBROUTINE MapFromOrigToPC_Preconditioner_2D( myPC, uHigh, uLow, nS, nP, nEl )
 ! S/R MapFromOrigToPC
 ! 
 !  This subroutine maps the array "uOrig" to the "uLin".
 !  "uLin" is an array indexed over the linear (preconditioner) mesh and "uOrig"
 !  is the array indexed over the higher order SE mesh.
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS( CGsemPreconditioner_2D ), INTENT(in) :: myPC
   INTEGER, INTENT(in)                         :: nS, nP, nEl
   REAL(prec), INTENT(in)                      :: uHigh(0:nS,0:nP,1:nEl)
   REAL(prec), INTENT(out)                     :: uLow(0:myPC % nS,0:myPC % nP,1:myPC % mesh % nElems)
   ! LOCAL
   INTEGER :: iEl, iS, iP
   INTEGER :: i, j, k
   REAL(prec) :: uInt(0:nS,0:myPC % nP)


      uLow = uHigh

 END SUBROUTINE MapFromOrigToPC_Preconditioner_2D
!
!
!
 SUBROUTINE MapFromPCToOrig_Preconditioner_2D( myPC, uHigh, uLow, nS, nP, nEl )
 ! S/R MapFromPCToOrig
 ! 
 !  This subroutine maps the array "uLin" to the array "uOrig".
 !  The "solution" is an array indexed over the linear (preconditioner) mesh and "uOrig"
 !  is the array indexed over the higher order SE mesh.
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS( CGsemPreconditioner_2D ), INTENT(in) :: myPC
   INTEGER, INTENT(in)                         :: nS, nP, nEl
   REAL(prec), INTENT(out)                     :: uHigh(0:nS,0:nP,1:nEl)
   REAL(prec), INTENT(in)                      :: uLow(0:myPC % nS,0:myPC % nP,1:myPC % mesh % nElems)
   ! LOCAL
   INTEGER :: iEl, iS, iP
   INTEGER :: i, j, k

      uHigh = uLow

 END SUBROUTINE MapFromPCToOrig_Preconditioner_2D
!
!
!
 SUBROUTINE Mask_CGsemPreconditioner_2D( myPC, u )
 ! S/R Mask_CGsemPreconditioner_2D
 !
 !
 ! ============================================================================= !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS( CGsemPreconditioner_2D ), INTENT(inout) :: myPC
   REAL(prec), INTENT(inout)                :: u(0:myPC % nS, &
                                                 0:myPC % nP, &
                                                 1:myPC % mesh % nElems )
   !LOCAL
   INTEGER :: nS, nP, eID, iNode, nID, i, j, nEl, nNodes
   INTEGER :: nodeType

      nS     = myPC % nS
      nP     = myPC % nP
      nEl    = myPC % mesh % nElems
      nNodes = myPC % mesh % nNodes
   
   

      ! Mask out the corner-nodes appropriately ( Leave the primary element intact )
      DO iNode = 1, nNodes ! Loop over the corner nodes

         CALL myPC % mesh % nodes(iNode) % RewindElementList( )
         CALL myPC % mesh % GetNodeType( iNode, nodeType )

         IF( nodeType /= DIRICHLET ) then ! this is NOT a prescribed (Dirichlet) boundary
                                                                          ! The data in this non-prescribed node is preserved

            CALL myPC % mesh % nodes(iNode) % AdvanceElementList( ) ! skip the first element in the list so it's data is not masked

         ENDIF

         DO WHILE( myPC % mesh % nodes(iNode) % CurrentElemExists( ) ) ! while there are elements in the node-to-element list

            CALL myPC % mesh % nodes(iNode) % GetCurrentElementID( eID ) ! Get the element ID and the local node ID (1->4)
            CALL myPC % mesh % GetLocalNodeID( iNode, nID )

            i = myPC % mesh % cornerMap(1, nID) ! Grab the index over the first computational direction associated with this corner node
            j = myPC % mesh % cornerMap(2, nID) ! Grab the index over the second computational direction associated with this corner node
 
            u(i,j,eID) = ZERO ! Mask out the solution at this element at this corner node

            CALL myPC % mesh % nodes(iNode) % AdvanceElementList( )

         ENDDO ! while we have elements in the node-to-element list

      ENDDO


 END SUBROUTINE Mask_CGsemPreconditioner_2D
!
!
!
 SUBROUTINE UnMask_CGsemPreconditioner_2D( myPC, u )
 ! S/R UnMask
 !
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS( CGsemPreconditioner_2D ), INTENT(inout) :: myPC
   REAL(prec), INTENT(inout)                :: u(0:myPC % nS, &
                                                 0:myPC % nP, &
                                                 1:myPC % mesh % nElems )
   !LOCAL
   INTEGER :: nS, nP, eID, iNode, nID, i, j, eID0, i0, j0, nEl, nNodes
   INTEGER :: nodeType

      nS     = myPC % nS
      nP     = myPC % nP
      nEl    = myPC % mesh % nElems
      nNodes = myPC % mesh % nNodes
   
      ! Unmask the corner-nodes
      DO iNode = 1, nNodes ! Loop over the corner nodes
         
         CALL myPC % mesh % GetNodeType( iNode, nodeType )

         IF( nodeType /= DIRICHLET ) then ! this is NOT a prescribed (Dirichlet) boundary

            CALL myPC % mesh % nodes(iNode) % RewindElementList( )

            ! Gather the element and local node IDs of the corner which contains the unmasked solution values.
            CALL myPC % mesh % nodes(iNode) % GetCurrentElementID( eID0 ) 
            CALL myPC % mesh % GetLocalNodeID( iNode, nID ) ! Get the element ID and the local node ID (1->4)

            i0 = myPC % mesh % cornerMap(1, nID) ! Grab the index over the first computational direction associated with this corner node
            j0 = myPC % mesh % cornerMap(2, nID) ! Grab the index over the second computational direction associated with this corner node

            CALL myPC % mesh % nodes(iNode) % AdvanceElementList( ) ! move on to the next element which shares this node


            ! Loop over the elements in this corner-node's connectivity list
            DO WHILE( myPC % mesh % nodes(iNode) % CurrentElemExists( ) ) ! while there are elements in the node-to-element list

               CALL myPC % mesh % nodes(iNode) % GetCurrentElementID( eID ) ! Get the element ID and the local node ID (1->4)
               CALL myPC % mesh % GetLocalNodeID( iNode, nID )
 
               i = myPC % mesh % cornerMap(1, nID) ! Grab the index over the first computational direction associated with this corner node
               j = myPC % mesh % cornerMap(2, nID) ! Grab the index over the second computational direction associated with this corner node
  
               u(i,j,eID) = u(i0,j0,eID0) ! copy the solution from the unmasked solution

               CALL myPC % mesh % nodes(iNode) % AdvanceElementList( )

            ENDDO ! while we have elements in the node-to-element list

         ENDIF

      ENDDO

 END SUBROUTINE UnMask_CGsemPreconditioner_2D
!
!
!
 SUBROUTINE GlobalSum_CGsemPreconditioner_2D( myPC, u )
 ! S/R GlobalSum
 !
 !  Adds together the shared edge and node contributions for the CGSEM method.
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS( CGsemPreconditioner_2D ), INTENT(inout) :: myPC
   REAL(prec), INTENT(inout)                :: u(0:myPC % nS, &
                                                 0:myPC % nP, &
                                                 1:myPC % mesh % nElems )
   !LOCAL
   INTEGER    :: eID, iNode, nID, i, j, nNodes, nS, nP, nEl
   INTEGER    :: nodeType 
   REAL(prec) :: theSum

      nS     = myPC % nS
      nP     = myPC % nP
      nEl    = myPC % mesh % nElems
      nNodes = myPC % mesh % nNodes


      DO iNode = 1, nNodes ! loop over the corner nodes
      
         CALL myPC % mesh % GetNodeType( iNode, nodeType )

         IF( nodeType /= DIRICHLET ) then ! this node does NOT lie on a prescribed (Dirichlet) boundary

             ! Rewind to the head of the linked list
            CALL myPC % mesh % nodes(iNode) % RewindElementList( )

            ! *** First, compute the sum from the contributing elements *** !
            theSum = ZERO 

            ! Loop over the elements in this corner-node's connectivity list
            DO WHILE( myPC % mesh % nodes(iNode) % CurrentElemExists( ) ) ! while there are elements in the node-to-element list

               CALL myPC % mesh % nodes(iNode) % GetCurrentElementID( eID )
               CALL myPC % mesh % GetLocalNodeID( iNode, nID )

               i = myPC % mesh % cornerMap(1, nID) ! Grab the index over the first computational direction associated with this corner node
               j = myPC % mesh % cornerMap(2, nID) ! Grab the index over the second computational direction associated with this corner node
  
               theSum = theSum + u(i,j,eID) 

               CALL myPC % mesh % nodes(iNode) % AdvanceElementList( )

            ENDDO ! while we have elements in the node-to-element list

            ! *** Now, copy the sum to the array "u" for each element which contributed to the sum *** !

             ! Rewind to the head of the linked list
            CALL myPC % mesh % nodes(iNode) % RewindElementList( )

            DO WHILE( myPC % mesh % nodes(iNode) % CurrentElemExists( ) ) ! while there are elements in the node-to-element list

               CALL myPC % mesh % nodes(iNode) % GetCurrentElementID( eID )
               CALL myPC % mesh % GetLocalNodeID( iNode, nID ) ! Get the element ID and the local node ID (1->4)

               i = myPC % mesh % cornerMap(1, nID) ! Grab the index over the first computational direction associated with this corner node
               j = myPC % mesh % cornerMap(2, nID) ! Grab the index over the second computational direction associated with this corner node
  
               u(i,j,eID) = theSum

               CALL myPC % mesh % nodes(iNode) % AdvanceElementList( )

            ENDDO ! while we have elements in the node-to-element list

         ENDIF

      ENDDO

 END SUBROUTINE GlobalSum_CGsemPreconditioner_2D
!
!
!
 SUBROUTINE CalculateGradient_CGsemPreconditioner_2D( myPC, iEl, u, dudx, dudy )
 ! S/R CalculateGradient
 !
 !  Calculates the gradient of the solution in computational coordinates within 
 !  a single element
 !
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS( CGsemPreconditioner_2D ), INTENT(in) :: myPC 
   INTEGER, INTENT(in)                   :: iEl
   REAL(prec), INTENT(in)                :: u(0:myPC % nS, 0:myPC % nP)
   REAL(prec), INTENT(out)               :: dudx(0:myPC % nS, 0:myPC % nP)
   REAL(prec), INTENT(out)               :: dudy(0:myPC % nS, 0:myPC % nP) 
   ! LOCAL
   INTEGER :: iS, iP, nS, nP
   REAL(prec) :: locU(0:myPC % nS, 0:myPC % nP)
   REAL(prec) :: duds(0:myPC % nS, 0:myPC % nP)
   REAL(prec) :: dudp(0:myPC % nS, 0:myPC % nP)
   REAL(prec) :: dMatS(0:myPC % nS, 0:myPC % nS)
   REAL(prec) :: dMatP(0:myPC % nP, 0:myPC % nP)
   REAL(prec) :: sTemp(0:myPC % nS)
   REAL(prec) :: pTemp(0:myPC % nP)
   REAL(prec) :: J, dxds, dxdp, dyds, dydp

      nS = myPC % nS
      nP = myPC % nP
      CALL myPC % SpectralOps % GetDerivativeMatrix( dMatS, dMatP ) 
      locU = myPC % filter % ApplyFilter( u )
   
      DO iP = 0, nP ! Loop over the second computational direction
 
         sTemp = locU(0:nS,iP)
         ! Calculate the derivative in the first computational direction
         duds(0:nS,iP) = MATMUL( dMatS, sTemp )

      ENDDO ! iP, loop over the second computational direction

      DO iS = 0, nS ! Loop over the first computational direction

         pTemp = locU(iS,0:nP)
         ! Calculate the derivative in the second computational direction
         dudp(iS,0:nP) = MATMUL( dMatP, pTemp )

      ENDDO ! iS, loop over the first computational direction


      ! Calculate the gradient in physical space
      DO iP = 0, nP
         DO iS = 0, nS

            CALL myPC % mesh % GetJacobianAtNode( iEl, J, iS, iP )
            CALL myPC % mesh % GetCovariantMetricsAtNode( iEl, dxds, dxdp, dyds, dydp, iS, iP )

            dudx(iS,iP) = ( dydp*duds(iS,iP) - dyds*dudp(iS,iP) )/J
            dudy(iS,iP) = ( dxds*dudp(iS,iP) - dxdp*duds(iS,iP) )/J        

         ENDDO
      ENDDO

 END SUBROUTINE CalculateGradient_CGsemPreconditioner_2D
!
!
!
 SUBROUTINE FluxDivergence_CGsemPreconditioner_2D( myPC, iEl, u, Lu )
 ! S/R FluxDivergence
 !
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS( CGsemPreconditioner_2D ), INTENT(inout) :: myPC
   INTEGER, INTENT(in)                      :: iEl
   REAL(prec), INTENT(in)                   :: u(0:myPC % nS, 0:myPC % nP)
   REAL(prec), INTENT(out)                  :: Lu(0:myPC % nS, 0:myPC % nP)
   ! LOCAL
   INTEGER    :: iS, iP, nS, nP
   REAL(prec) :: F1(0:myPC % nS, 0:myPC % nP)
   REAL(prec) :: F2(0:myPC % nS, 0:myPC % nP)
   REAL(prec) :: dF1ds(0:myPC % nS, 0:myPC % nP)
   REAL(prec) :: dF2dp(0:myPC % nS, 0:myPC % nP)
   REAL(prec) :: dudx(0:myPC % nS, 0:myPC % nP)
   REAL(prec) :: dudy(0:myPC % nS, 0:myPC % nP)
   REAL(prec) :: dMatS(0:myPC % nS, 0:myPC % nS)
   REAL(prec) :: dMatSt(0:myPC % nS, 0:myPC % nS)
   REAL(prec) :: dMatP(0:myPC % nP, 0:myPC % nP)
   REAL(prec) :: dMatPt(0:myPC % nP, 0:myPC % nP)
   REAL(prec) :: ws(0:myPC % nS)
   REAL(prec) :: wp(0:myPC % nP)
   REAL(prec) :: sTemp(0:myPC % nS)
   REAL(prec) :: pTemp(0:myPC % nP)
   REAL(prec) :: dxds, dxdp, dyds, dydp

      nS = myPC % nS
      nP = myPC % nP

      CALL myPC % SpectralOps % GetQuadratureWeights( ws, wp )
      CALL myPC % SpectralOps % GetDerivativeMatrix( dMatS, dMatP )
      dMatSt = TRANSPOSE( dMatS )
      dMatPt = TRANSPOSE( dMatP )

      CALL myPC % CalculateGradient( iEl, u, dudx, dudy ) ! calculate the gradient in physical coordinates

      ! Now apply the metric terms
      DO iP = 0, nP ! Loop over the second computational direction
         DO iS = 0, nS ! Loop over the first computational direction

            CALL myPC % mesh % GetCovariantMetricsAtNode( iEl, dxds, dxdp, dyds, dydp, iS, iP )

            ! Contravariant flux calculation
            F1(iS,iP) = myPC % fluxCoeff(iS,iP,iEl)*( dudx(iS,iP)*dydp - dudy(iS,iP)*dxdp )*ws(iS)
            F2(iS,iP) = myPC % fluxCoeff(iS,iP,iEl)*( dudy(iS,iP)*dxds - dudx(iS,iP)*dyds )*wp(iP)        

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
       

 END SUBROUTINE FluxDivergence_CGsemPreconditioner_2D 
!
!
!
 SUBROUTINE MatrixAction_CGsemPreconditioner_2D( myPC, u, Au )
 ! S/R MatrixAction
 ! 
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS(CGsemPreconditioner_2D), INTENT(inout) :: myPC
   REAL(prec), INTENT(inout)              :: u(0:myPC % nS, &
                                               0:myPC % nP, &
                                               1:myPC % mesh % nElems )
   REAL(prec), INTENT(out)                :: Au(0:myPC % nS, &
                                                0:myPC % nP, &
                                                1:myPC % mesh % nElems  )
   ! Local
   INTEGER :: iEl, nS, nP, nEl
   REAL(prec) :: temp(0:myPC % nS, &
                      0:myPC % nP )
   REAL(prec) :: Atemp(0:myPC % nS, &
                       0:myPC % nP )

      nS  = myPC % nS
      nP  = myPC % nP
      nEl = myPC % mesh % nElems
     
      CALL myPC % UnMask( u ) ! unmask the solution

!$OMP PARALLEL PRIVATE( temp, Atemp )
!$OMP DO 
      DO iEl = 1, nEl
         temp(0:nS,0:nP) = u(0:nS,0:nP,iEl)
         CALL myPC % FluxDivergence( iEl, temp, Atemp )
         Au(0:nS,0:nP,iEl) =  Atemp(0:nS,0:nP)
      ENDDO
!$OMP END DO
!$OMP  FLUSH ( Au )
!$OMP END PARALLEL

     ! Add the contributions from the shared edges and corners
     CALL myPC % GlobalSum( Au )

     ! Mask the edges in preparation for residual calculation
     CALL myPC % Mask( Au )

     ! Re-mask the solution
     CALL myPC % Mask( u )

 END SUBROUTINE MatrixAction_CGsemPreconditioner_2D
!
!
!
 FUNCTION DotProduct_CGsemPreconditioner_2D( myPC, u, v ) RESULT( uDotv )
 ! FUNCTION DotProduct
 ! Computes the dot product of two "vectors" whose data are
 ! stored in the "(iS,iP)" format
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS( CGsemPreconditioner_2D ) :: myPC
   REAL(prec)                      :: u(0:myPC % nS, &
                                        0:myPC % nP, &
                                        1:myPC % mesh % nElems )
   REAL(prec)                      :: v(0:myPC % nS, &
                                        0:myPC % nP, &
                                        1:myPC % mesh % nElems )
   REAL(prec)                      :: uDotv
   ! Local
   INTEGER :: iP, nS, nP, iEl, nEl
   REAL(prec) :: uu(0:myPC % nS, &
                    0:myPC % nP, &
                    1:myPC % mesh % nElems )
   REAL(prec) :: vv(0:myPC % nS, &
                    0:myPC % nP, &
                    1:myPC % mesh % nElems )
   REAL(prec) :: uloc(0:myPC % nS)
   REAL(prec) :: vloc(0:myPC % nS)
   
      nS  = myPC % nS 
      nP  = myPC % nP
      nEl = myPC % mesh % nElems

      uu = u
      vv = v
      ! mask the arrays
      CALL myPC % Mask( uu )
      CALL myPC % Mask( vv )
 
      uDotv = ZERO
      DO iEl = 1, nEl
         DO iP = 0, nP
            uloc  = uu(0:nS,iP,iEl)
            vLoc  = vv(0:nS,iP,iEl)
            uDotV = uDotV + DOT_PRODUCT( uloc, vloc )
         ENDDO
      ENDDO

 END FUNCTION DotProduct_CGsemPreconditioner_2D
!
!
!
 SUBROUTINE Solve_CGsemPreconditioner_2D( myPC, x, b, nS, nP, nEl, ioerr )
 !  S/R Solve_CGsemPreconditioner_2D
 !
 !  This subroutine solves the system Hx = b using the conjugate gradient method.
 !  The solver is used as a preconditioner, and does not need to be iterated to convergence
 !
 !  On output ioerr is set to an error checking flag. 
 !  If ioerr ==  0, the method converged within the maximum number of iterations.
 !     ioerr == -1, the method did not converge within the maximum number of iterations.
 !     ioerr == -2, something that is not caught by the current construct happened.
 !  
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS( CGsemPreconditioner_2D ), INTENT(inout) :: myPC
   INTEGER, INTENT(in)                            :: nS, nP, nEl
   REAL(prec), INTENT(out)                        :: x(0:nS,0:nP,1:nEl)
   REAL(prec), INTENT(in)                         :: b(0:nS,0:nP,1:nEl)
   INTEGER, INTENT(out)                           :: ioerr
   ! LOCAL
   INTEGER    :: iter
   INTEGER    :: nIt
   REAL(prec) :: TOL
   REAL(prec) :: u(0:myPC % nS, 0:myPC % nP, 1:myPC % mesh % nElems)
   REAL(prec) :: r(0:myPC % nS, 0:myPC % nP, 1:myPC % mesh % nElems)
   REAL(prec) :: v(0:myPC % nS, 0:myPC % nP, 1:myPC % mesh % nElems)
   REAL(prec) :: z(0:myPC % nS, 0:myPC % nP, 1:myPC % mesh % nElems)
   REAL(prec) :: w(0:myPC % nS, 0:myPC % nP, 1:myPC % mesh % nElems)
   REAL(prec) :: rNorm
   REAL(prec) :: a, c, num, den, r0
   
      ioerr = -2
      nIt = myPC % maxIters
      TOL = myPC % tolerance

      ! Map the right-hand-side to the preconditioner mesh
      CALL myPC % MapFromOrigToPC( b, r, nS, nP, nEl )
      
      u = ZERO
      r0 = sqrt( myPC % DotProduct( r, r ) ) 
!      print*, r0
      DO iter = 1,nIt ! Loop over the CG iterates
      
         w = r
         num = myPC % DotProduct( r, w )

         IF( iter == 1) THEN
            v = w
         ELSE
            c = num/den
            v = w + c*v
         ENDIF

         CALL myPC % MatrixAction( v, z )
         a = num/myPC % DotProduct( v, z )
         u = u + a*v
         r = r - a*z
         
         den = num

         rNorm = sqrt( myPC % DotProduct(r,r) )
!         print*, rNorm/r0, rNorm/r0 < TOL
         IF( rNorm/r0 < TOL ) then
           EXIT
           ioerr = 0
         ENDIF

      ENDDO ! iter, loop over the CG iterates 

      IF( SQRT(rNorm)/r0 > TOL )THEN
      !   PRINT*, 'MODULE IterativeSolvers : ConjugateGradient failed to converge '
      !   PRINT*, 'Last L-2 residual : ', sqrt(rNorm)
         ioerr=-1
      ENDIF   

      CALL myPC % UnMask( u )
      CALL myPC % MapFromPCtoOrig( x, u, nS, nP, nEl ) !

 END SUBROUTINE Solve_CGsemPreconditioner_2D
!
!
! 
 SUBROUTINE ConstructPCMatrix( myPC )
 !
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS( CGsemPreconditioner_2D ), INTENT(inout) :: myPC
   ! LOCAL
   INTEGER                 :: iEl, iS, iP, jEl, jS, jP, nEl, nS, nP
   INTEGER                 :: row, col, nfree, fUnit
   REAL(prec)              :: ei(0:myPC % nS, 0:myPC % nP, 1:myPC % mesh % nElems)
   REAL(prec)              :: Aei(0:myPC % nS, 0:myPC % nP, 1:myPC % mesh % nElems)
   REAL(prec)              :: thismask(0:myPC % nS, 0:myPC % nP, 1:myPC % mesh % nElems)
   REAL(prec), ALLOCATABLE :: A(:,:) 

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
            FILE = 'PCmat.txt', &
            FORM = 'FORMATTED' )

      DO row = 1, nfree
         WRITE(fUnit,*) A(row,:)
      ENDDO

      CLOSE(fUnit)

               

 END SUBROUTINE

END MODULE CGsemCoarsePreconditioner_Class
