! HexElementClass.f90
! 
! Copyright 2015 Joseph Schoonover <schoonover.numerics@gmail.com>
! 
! HexElementClass.f90 is part of the Spectral Element Libraries in Fortran (SELF).
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

MODULE HexElementClass
! ========================================= Logs ================================================= !
!2016-05-11  Joseph Schoonover  schoonover.numerics@gmail.com 
!
! //////////////////////////////////////////////////////////////////////////////////////////////// ! 

! src/common/
USE ModelPrecision
USE ConstantsDictionary
USE ModelFlags
! src/interp/
USE Lagrange_3D_Class
! src/geom/
USE SurfaceClass_3D
USE MappedGeometryClass_3D


IMPLICIT NONE

   TYPE HexElement
      INTEGER                         :: nS, nP, nQ, nmax
      INTEGER, PRIVATE                :: nodeIDs(1:nHexNodes)   ! Corner Node ID's
      INTEGER, PRIVATE                :: neighbors(1:nHexFaces) ! Elements IDs for the neighbors
      INTEGER, PRIVATE                :: globElID               ! Global Element ID 
      TYPE(MappedGeometry_3D),PRIVATE :: geometry               ! Contains the element metric terms, etc.


      CONTAINS

      PROCEDURE :: Build => Build_HexElement
      PROCEDURE :: Trash => Trash_HexElement
    
      PROCEDURE :: SetNodeIDs => SetNodeIDs_HexElement
      PROCEDURE :: GetNodeIDs => GetNodeIDs_HexElement
      PROCEDURE :: SetNodeID => SetNodeID_HexElement
      PROCEDURE :: GetNodeID => GetNodeID_HexElement
      PROCEDURE :: SetElementID => SetElementID_HexElement
      PROCEDURE :: GetElementID => GetElementID_HexElement
      PROCEDURE :: SetNeighbor => SetNeighbor_HexElement
      PROCEDURE :: GetNeighbor => GetNeighbor_HexElement
      PROCEDURE :: SetSouthernNeighbor => SetSouthernNeighbor_HexElement
      PROCEDURE :: GetSouthernNeighbor => GetSouthernNeighbor_HexElement
      PROCEDURE :: SetNorthernNeighbor => SetNorthernNeighbor_HexElement
      PROCEDURE :: GetNorthernNeighbor => GetNorthernNeighbor_HexElement
      PROCEDURE :: SetEasternNeighbor  => SetEasternNeighbor_HexElement
      PROCEDURE :: GetEasternNeighbor  => GetEasternNeighbor_HexElement
      PROCEDURE :: SetWesternNeighbor  => SetWesternNeighbor_HexElement
      PROCEDURE :: GetWesternNeighbor  => GetWesternNeighbor_HexElement
      PROCEDURE :: SetBottomNeighbor   => SetBottomNeighbor_HexElement
      PROCEDURE :: GetBottomNeighbor   => GetBottomNeighbor_HexElement
      PROCEDURE :: SetTopNeighbor      => SetTopNeighbor_HexElement
      PROCEDURE :: GetTopNeighbor      => GetTopNeighbor_HexElement      

      ! MappedGeometry_3D Wrapper Routines
      PROCEDURE :: SetNumberOfNodes => SetNumberOfNodes_HexElement
      PROCEDURE :: GetNumberOfNodes => GetNumberOfNodes_HexElement
      PROCEDURE :: SetPositions => SetPositions_HexElement
      PROCEDURE :: GetPositions => GetPositions_HexElement
      PROCEDURE :: SetPositionAtNode => SetPositionAtNode_HexElement
      PROCEDURE :: GetPositionAtNode => GetPositionAtNode_HexElement
      PROCEDURE :: SetJacobian => SetJacobian_HexElement
      PROCEDURE :: GetJacobian => GetJacobian_HexElement
      PROCEDURE :: SetJacobianAtNode => SetJacobianAtNode_HexElement
      PROCEDURE :: GetJacobianAtNode => GetJacobianAtNode_HexElement
      PROCEDURE :: SetCovariantMetrics => SetCovariantMetrics_HexElement
      PROCEDURE :: GetCovariantMetrics => GetCovariantMetrics_HexElement
      PROCEDURE :: SetCovariantMetricsAtNode => SetCovariantMetricsAtNode_HexElement
      PROCEDURE :: GetCovariantMetricsAtNode => GetCovariantMetricsAtNode_HexElement
      PROCEDURE :: SetBoundaryLocation => SetBoundaryLocation_HexElement
      PROCEDURE :: GetBoundaryLocation => GetBoundaryLocation_HexElement
      PROCEDURE :: SetBoundaryLocationAtNode => SetBoundaryLocationAtNode_HexElement
      PROCEDURE :: GetBoundaryLocationAtNode => GetBoundaryLocationAtNode_HexElement
      PROCEDURE :: SetBoundaryNormalAtNode => SetBoundaryNormalAtNode_HexElement
      PROCEDURE :: GetBoundaryNormalAtNode => GetBoundaryNormalAtNode_HexElement
      PROCEDURE :: ScaleGeometry => ScaleGeometry_HexElement
      
      PROCEDURE :: IsInside => IsInside_HexElement
      
   END TYPE HexElement

 INTEGER, PRIVATE :: nDims = 3

 CONTAINS
!
!
!==================================================================================================!
!------------------------------- Manual Constructors/Destructors ----------------------------------!
!==================================================================================================!
!
!
 SUBROUTINE Build_HexElement( myElement, nodeIDs, eID, boundSurfaces, interp  )
 ! S/R Build
 ! 
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS( HexElement ), INTENT(inout) :: myElement
   INTEGER, INTENT(in)                :: nodeIDs(1:nHexNodes), eID
   TYPE( Surface_3D ), INTENT(in)     :: boundSurfaces(1:nHexFaces)
   TYPE( Lagrange_3D ), INTENT(in)    :: interp
   
      CALL myElement % SetNodeIDs( nodeIDs )

      CALL myElement % SetElementID( eID )

      CALL myElement % geometry % Build( interp, boundSurfaces  )
      
      myElement % nS = interp % nS
      myElement % nP = interp % nP
      myElement % nQ = interp % nQ
      myElement % nMax = MAX( interp % nS, interp % nP, interp % nQ )

 END SUBROUTINE Build_HexElement
!
!
!
 SUBROUTINE Trash_HexElement( myElement )
 ! S/R Trash
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS( HexElement ), INTENT(inout) :: myElement

      CALL myElement % geometry % Trash( )

 END SUBROUTINE Trash_HexElement
!
!
!==================================================================================================!
!--------------------------------------- Accessors ------------------------------------------------!
!==================================================================================================!
!
!
 SUBROUTINE SetNodeIDs_HexElement( myElement, nodeIDs )
 ! S/R SetNodeIDs
 !
 ! =============================================================================================== !
   IMPLICIT NONE
   CLASS( HexElement ), INTENT(inout) :: myElement
   INTEGER, INTENT(in)                :: nodeIDs(1:nHexNodes)
   
      myElement % nodeIDs = nodeIDs

 END SUBROUTINE SetNodeIDs_HexElement
!
!
!
 FUNCTION GetNodeIDs_HexElement( myElement ) RESULT( nodeIDs )
 ! S/R GetNodeIDs
 !
 ! =============================================================================================== !
   IMPLICIT NONE
   CLASS( HexElement ) :: myElement
   INTEGER             :: nodeIDs(1:nHexNodes)
   
      nodeIDs = myElement % nodeIDs

 END FUNCTION GetNodeIDs_HexElement
!
!
!
 SUBROUTINE SetNodeID_HexElement( myElement, localID, nodeID )
 ! S/R SetNodeID
 !
 ! =============================================================================================== !
   IMPLICIT NONE
   CLASS( HexElement ), INTENT(inout) :: myElement
   INTEGER, INTENT(in)                 :: localID, nodeID
   
      myElement % nodeIDs(localID) = nodeID

 END SUBROUTINE SetNodeID_HexElement
!
!
!
 FUNCTION GetNodeID_HexElement( myElement, localID ) RESULT( nodeID )
 ! S/R GetNodeID
 !
 ! =============================================================================================== !
   IMPLICIT NONE
   CLASS( HexElement ) :: myElement
   INTEGER             :: localID
   INTEGER             :: nodeID
   
      nodeID = myElement % nodeIDs(localID)

 END FUNCTION GetNodeID_HexElement
!
!
!
 SUBROUTINE SetElementID_HexElement( myElement, eID )
 ! S/R SetElementID
 !
 ! =============================================================================================== !
   IMPLICIT NONE
   CLASS( HexElement ), INTENT(inout) :: myElement
   INTEGER, INTENT(in)                 :: eID
   
      myElement % globElID = eID

 END SUBROUTINE SetElementID_HexElement
!
!
!
 FUNCTION GetElementID_HexElement( myElement ) RESULT( eID )
 ! S/R GetElementID
 !
 ! =============================================================================================== !
   IMPLICIT NONE
   CLASS( HexElement ) :: myElement
   INTEGER             :: eID
   
      eID = myElement % globElID

 END FUNCTION GetElementID_HexElement
!
!
!
 SUBROUTINE SetNeighbor_HexElement( myElement, sID, eID )
  ! S/R SetNeighbor
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS( HexElement ), INTENT(inout) :: myElement
   INTEGER, INTENT(in)                 :: sID
   INTEGER, INTENT(in)                 :: eID
   
      myElement % neighbors(sID) = eID
      
 END SUBROUTINE SetNeighbor_HexElement
!
!
!
 FUNCTION GetNeighbor_HexElement( myElement, sID )RESULT( eID )
 ! S/R GetNeighbor
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS( HexElement ) :: myElement
   INTEGER             :: sID
   INTEGER             :: eID
   
      eID = myElement % neighbors(sID)
      
 END FUNCTION GetNeighbor_HexElement
!
!
!
 SUBROUTINE SetSouthernNeighbor_HexElement( myElement, eID )
 ! S/R SetSouthernNeighbor
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS( HexElement ), INTENT(inout) :: myElement
   INTEGER, INTENT(in)                :: eID
   
      myElement % neighbors(south) = eID
      
 END SUBROUTINE SetSouthernNeighbor_HexElement
!
!
!
 FUNCTION GetSouthernNeighbor_HexElement( myElement )RESULT( eID )
 ! S/R GetSouthernNeighbor
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS( HexElement ) :: myElement
   INTEGER             :: eID
   
      eID = myElement % neighbors(south)
      
 END FUNCTION GetSouthernNeighbor_HexElement
!
!
!
 SUBROUTINE SetNorthernNeighbor_HexElement( myElement, eID )
 ! S/R SetNorthernNeighbor
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS( HexElement ), INTENT(inout) :: myElement
   INTEGER, INTENT(in)                :: eID
   
      myElement % neighbors(north) = eID
      
 END SUBROUTINE SetNorthernNeighbor_HexElement
!
!
!
 FUNCTION GetNorthernNeighbor_HexElement( myElement )RESULT( eID )
 ! S/R GetNorthernNeighbor
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS( HexElement ) :: myElement
   INTEGER             :: eID
   
      eID = myElement % neighbors(north)
      
 END FUNCTION GetNorthernNeighbor_HexElement
!
!
!
 SUBROUTINE SetEasternNeighbor_HexElement( myElement, eID )
 ! S/R SetEasternNeighbor
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS( HexElement ), INTENT(inout) :: myElement
   INTEGER, INTENT(in)                :: eID
   
      myElement % neighbors(east) = eID
      
 END SUBROUTINE SetEasternNeighbor_HexElement
!
!
!
 FUNCTION GetEasternNeighbor_HexElement( myElement )RESULT( eID )
 ! S/R GetEasternNeighbor
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS( HexElement ) :: myElement
   INTEGER             :: eID
   
      eID = myElement % neighbors(east)
      
 END FUNCTION GetEasternNeighbor_HexElement
!
!
!
 SUBROUTINE SetWesternNeighbor_HexElement( myElement, eID )
 ! S/R SetWesternNeighbor
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS( HexElement ), INTENT(inout) :: myElement
   INTEGER, INTENT(in)                :: eID
   
      myElement % neighbors(west) = eID
      
 END SUBROUTINE SetWesternNeighbor_HexElement
!
!
!
 FUNCTION GetWesternNeighbor_HexElement( myElement )RESULT( eID )
 ! S/R GetWesternNeighbor
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS( HexElement ) :: myElement
   INTEGER             :: eID
   
      eID = myElement % neighbors(west)
      
 END FUNCTION GetWesternNeighbor_HexElement
!
!
!
 SUBROUTINE SetBottomNeighbor_HexElement( myElement, eID )
 ! S/R SetBottomNeighbor
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS( HexElement ), INTENT(inout) :: myElement
   INTEGER, INTENT(in)                :: eID
   
      myElement % neighbors(bottom) = eID
      
 END SUBROUTINE SetBottomNeighbor_HexElement
!
!
!
 FUNCTION GetBottomNeighbor_HexElement( myElement )RESULT( eID )
 ! S/R GetBottomNeighbor
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS( HexElement ) :: myElement
   INTEGER             :: eID
   
      eID = myElement % neighbors(bottom)
      
 END FUNCTION GetBottomNeighbor_HexElement
!
!
!
 SUBROUTINE SetTopNeighbor_HexElement( myElement, eID )
 ! S/R SetTopNeighbor
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS( HexElement ), INTENT(inout) :: myElement
   INTEGER, INTENT(in)                :: eID
   
      myElement % neighbors(top) = eID
      
 END SUBROUTINE SetTopNeighbor_HexElement
!
!
!
 FUNCTION GetTopNeighbor_HexElement( myElement )RESULT( eID )
 ! S/R GetTopNeighbor
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS( HexElement ) :: myElement
   INTEGER             :: eID
   
      eID = myElement % neighbors(top)
      
 END FUNCTION GetTopNeighbor_HexElement
!
! ---------------------------- MappedGeometryClass_3D Wrappers ----------------------------------- !
!
SUBROUTINE SetNumberOfNodes_HexElement( myElement, nS, nP, nQ )
 ! S/R SetNumberOfNodes
 !
 ! =============================================================================================== !
 ! DECLARATIONS
  IMPLICIT NONE
  CLASS(HexElement), INTENT(inout) :: myElement
  INTEGER, INTENT(in)              :: nS, nP, nQ
  
     CALL myElement % geometry % SetNumberOfNodes( nS, nP, nQ )
     
 END SUBROUTINE SetNumberOfNodes_HexElement
!
!
!
 SUBROUTINE GetNumberOfNodes_HexElement( myElement, nS, nP, nQ )
 ! S/R GetNumberOfNodes
 !
 ! =============================================================================================== !
 ! DECLARATIONS
  IMPLICIT NONE
  CLASS(HexElement), INTENT(in) :: myElement
  INTEGER, INTENT(out)          :: nS, nP, nQ
  
     CALL myElement % geometry % GetNumberOfNodes( nS, nP, nQ )
     
 END SUBROUTINE GetNumberOfNodes_HexElement
!
!
!
 SUBROUTINE SetPositions_HexElement( myElement, x, y, z )
 ! S/R SetPositions
 !
 ! =============================================================================================== !
 ! DECLARATIONS
  IMPLICIT NONE
  CLASS(HexElement), INTENT(inout) :: myElement
  REAL(prec), INTENT(in)           :: x(0:myElement % nS, 0:myElement % nP, 0:myElement % nQ)
  REAL(prec), INTENT(in)           :: y(0:myElement % nS, 0:myElement % nP, 0:myElement % nQ)
  REAL(prec), INTENT(in)           :: z(0:myElement % nS, 0:myElement % nP, 0:myElement % nQ)
  
     CALL myElement % geometry % SetPositions( x, y, z )
     
 END SUBROUTINE SetPositions_HexElement
!
!
!
 SUBROUTINE GetPositions_HexElement( myElement, x, y, z )
 ! S/R GetPositions
 !
 ! =============================================================================================== !
 ! DECLARATIONS
  IMPLICIT NONE
  CLASS(HexElement), INTENT(in) :: myElement
  REAL(prec), INTENT(out)       :: x(0:myElement % nS, 0:myElement % nP, 0:myElement % nQ)
  REAL(prec), INTENT(out)       :: y(0:myElement % nS, 0:myElement % nP, 0:myElement % nQ)
  REAL(prec), INTENT(out)       :: z(0:myElement % nS, 0:myElement % nP, 0:myElement % nQ)
  
     CALL myElement % geometry % GetPositions( x, y, z )
     
 END SUBROUTINE GetPositions_HexElement
!
!
!
 SUBROUTINE SetPositionAtNode_HexElement( myElement, x, y, z, i, j, k )
 ! S/R SetPositionAtNode
 !
 ! =============================================================================================== !
 ! DECLARATIONS
  IMPLICIT NONE
  CLASS(HexElement), INTENT(inout) :: myElement
  REAL(prec), INTENT(in)            :: x, y, z
  INTEGER, INTENT(in)               :: i, j, k
  
     CALL myElement % geometry % SetPositionAtNode( x, y, z, i, j, k )
     
 END SUBROUTINE SetPositionAtNode_HexElement
!
!
!
 SUBROUTINE GetPositionAtNode_HexElement( myElement, x, y, z, i, j, k )
 ! S/R GetPositionAtNode
 !
 ! =============================================================================================== !
 ! DECLARATIONS
  IMPLICIT NONE
  CLASS(HexElement), INTENT(in) :: myElement
  REAL(prec), INTENT(out)       :: x, y, z
  INTEGER, INTENT(in)           :: i, j, k
  
     CALL myElement % geometry % GetPositionAtNode( x, y, z, i, j, k )
     
 END SUBROUTINE GetPositionAtNode_HexElement
!
!
!
 SUBROUTINE SetJacobian_HexElement( myElement, J )
 ! S/R SetJacobian
 !
 ! =============================================================================================== !
 ! DECLARATIONS
  IMPLICIT NONE
  CLASS(HexElement), INTENT(inout) :: myElement
  REAL(prec), INTENT(in)           :: J(0:myElement % nS, 0:myElement % nP, 0:myElement % nQ)
  
     CALL myElement % geometry % SetJacobian( J )
     
 END SUBROUTINE SetJacobian_HexElement
!
!
!
 SUBROUTINE GetJacobian_HexElement( myElement, J )
 ! S/R GetJacobian
 !
 ! =============================================================================================== !
 ! DECLARATIONS
  IMPLICIT NONE
  CLASS(HexElement), INTENT(in) :: myElement
  REAL(prec), INTENT(out)       :: J(0:myElement % nS, 0:myElement % nP, 0:myElement % nQ)
  
     CALL myElement % geometry % GetJacobian( J )
     
 END SUBROUTINE GetJacobian_HexElement
!
!
!
 SUBROUTINE SetJacobianAtNode_HexElement( myElement, J, iS, iP, iQ )
 ! S/R SetJacobianAtNode
 !
 ! =============================================================================================== !
 ! DECLARATIONS
  IMPLICIT NONE
  CLASS(HexElement), INTENT(inout) :: myElement
  REAL(prec), INTENT(in)           :: J
  INTEGER, INTENT(in)              :: iS, iP, iQ
  
     CALL myElement % geometry % SetJacobianAtNode( J, iS, iP, iQ )
     
 END SUBROUTINE SetJacobianAtNode_HexElement
!
!
!
 SUBROUTINE GetJacobianAtNode_HexElement( myElement, J, iS, iP, iQ )
 ! S/R GetJacobianAtNode
 !
 ! =============================================================================================== !
 ! DECLARATIONS
  IMPLICIT NONE
  CLASS(HexElement), INTENT(in) :: myElement
  REAL(prec), INTENT(out)       :: J
  INTEGER, INTENT(in)           :: iS, iP, iQ
  
     CALL myElement % geometry % GetJacobianAtNode( J, iS, iP, iQ )
     
 END SUBROUTINE GetJacobianAtNode_HexElement
!
!
!
 SUBROUTINE SetCovariantMetrics_HexElement( myElement, dxds, dxdp, dxdq, &
                                                       dyds, dydp, dydq, &
                                                       dzds, dzdp, dzdq )
 ! S/R SetCovariantMetrics
 !
 ! =============================================================================================== !
 ! DECLARATIONS
  IMPLICIT NONE
  CLASS(HexElement), INTENT(inout) :: myElement
  REAL(prec), INTENT(in)           :: dxds(0:myElement % nS, 0:myElement % nP, 0:myElement % nQ)
  REAL(prec), INTENT(in)           :: dxdp(0:myElement % nS, 0:myElement % nP, 0:myElement % nQ)
  REAL(prec), INTENT(in)           :: dxdq(0:myElement % nS, 0:myElement % nP, 0:myElement % nQ)
  REAL(prec), INTENT(in)           :: dyds(0:myElement % nS, 0:myElement % nP, 0:myElement % nQ)
  REAL(prec), INTENT(in)           :: dydp(0:myElement % nS, 0:myElement % nP, 0:myElement % nQ)
  REAL(prec), INTENT(in)           :: dydq(0:myElement % nS, 0:myElement % nP, 0:myElement % nQ)
  REAL(prec), INTENT(in)           :: dzds(0:myElement % nS, 0:myElement % nP, 0:myElement % nQ)
  REAL(prec), INTENT(in)           :: dzdp(0:myElement % nS, 0:myElement % nP, 0:myElement % nQ)
  REAL(prec), INTENT(in)           :: dzdq(0:myElement % nS, 0:myElement % nP, 0:myElement % nQ)
  
     CALL myElement % geometry % SetCovariantMetrics( dxds, dxdp, dxdq, &
                                                      dyds, dydp, dydq, &
                                                      dzds, dzdp, dzdq )
     
 END SUBROUTINE SetCovariantMetrics_HexElement
!
!
!
 SUBROUTINE GetCovariantMetrics_HexElement( myElement, dxds, dxdp, dxdq, &
                                                       dyds, dydp, dydq, &
                                                       dzds, dzdp, dzdq )
 ! S/R GetCovariantMetrics
 !
 ! =============================================================================================== !
 ! DECLARATIONS
  IMPLICIT NONE
  CLASS(HexElement), INTENT(in) :: myElement
  REAL(prec), INTENT(out)       :: dxds(0:myElement % nS, 0:myElement % nP, 0:myElement % nQ)
  REAL(prec), INTENT(out)       :: dxdp(0:myElement % nS, 0:myElement % nP, 0:myElement % nQ)
  REAL(prec), INTENT(out)       :: dxdq(0:myElement % nS, 0:myElement % nP, 0:myElement % nQ)
  REAL(prec), INTENT(out)       :: dyds(0:myElement % nS, 0:myElement % nP, 0:myElement % nQ)
  REAL(prec), INTENT(out)       :: dydp(0:myElement % nS, 0:myElement % nP, 0:myElement % nQ)
  REAL(prec), INTENT(out)       :: dydq(0:myElement % nS, 0:myElement % nP, 0:myElement % nQ)
  REAL(prec), INTENT(out)       :: dzds(0:myElement % nS, 0:myElement % nP, 0:myElement % nQ)
  REAL(prec), INTENT(out)       :: dzdp(0:myElement % nS, 0:myElement % nP, 0:myElement % nQ)
  REAL(prec), INTENT(out)       :: dzdq(0:myElement % nS, 0:myElement % nP, 0:myElement % nQ)
  
     CALL myElement % geometry % GetCovariantMetrics( dxds, dxdp, dxdq, &
                                                      dyds, dydp, dydq, &
                                                      dzds, dzdp, dzdq )
     
 END SUBROUTINE GetCovariantMetrics_HexElement
!
!
!
 SUBROUTINE SetCovariantMetricsAtNode_HexElement( myElement, dxds, dxdp, dxdq, &
                                                             dyds, dydp, dydq, &
                                                             dzds, dzdp, dzdq, &
                                                             iS, iP, iQ )
 ! S/R SetCovariantMetricsAtNode
 !
 ! =============================================================================================== !
 ! DECLARATIONS
  IMPLICIT NONE
  CLASS(HexElement), INTENT(inout) :: myElement
  REAL(prec), INTENT(in)           :: dxds, dxdp, dxdq, dyds, dydp, dydq, dzds, dzdp, dzdq
  INTEGER, INTENT(in)              :: iS, iP, iQ
  
     CALL myElement % geometry % SetCovariantMetricsAtNode( dxds, dxdp, dxdq, &
                                                            dyds, dydp, dydq, &
                                                            dzds, dzdp, dzdq, &
                                                            iS, iP, iQ )
     
 END SUBROUTINE SetCovariantMetricsAtNode_HexElement
!
!
!
 SUBROUTINE GetCovariantMetricsAtNode_HexElement( myElement, dxds, dxdp, dxdq, &
                                                             dyds, dydp, dydq, &
                                                             dzds, dzdp, dzdq, &
                                                             iS, iP, iQ )
 ! S/R GetCovariantMetricsAtNode
 !
 ! =============================================================================================== !
 ! DECLARATIONS
  IMPLICIT NONE
  CLASS(HexElement), INTENT(in) :: myElement
  REAL(prec), INTENT(out)        :: dxds, dxdp, dxdq, dyds, dydp, dydq, dzds, dzdp, dzdq
  INTEGER, INTENT(in)            :: iS, iP, iQ
  
     CALL myElement % geometry % GetCovariantMetricsAtNode( dxds, dxdp, dxdq, &
                                                            dyds, dydp, dydq, &
                                                            dzds, dzdp, dzdq, &
                                                            iS, iP, iQ )
     
 END SUBROUTINE GetCovariantMetricsAtNode_HexElement
!
!
!
 SUBROUTINE SetBoundaryLocation_HexElement( myElement, x, y, z, iBound )
 ! S/R SetBoundaryLocation
 !
 ! =============================================================================================== !
 ! DECLARATIONS
  IMPLICIT NONE
  CLASS(HexElement), INTENT(inout) :: myElement
  INTEGER, INTENT(in)              :: iBound
  REAL(prec), INTENT(in)           :: x(0:myElement % geometry % nMax, 0:myElement % geometry % nMax)
  REAL(prec), INTENT(in)           :: y(0:myElement % geometry % nMax, 0:myElement % geometry % nMax)
  REAL(prec), INTENT(in)           :: z(0:myElement % geometry % nMax, 0:myElement % geometry % nMax)
  
     CALL myElement % geometry % SetBoundaryLocation( x, y, z, iBound )
     
 END SUBROUTINE SetBoundaryLocation_HexElement
!
!
!
 SUBROUTINE GetBoundaryLocation_HexElement( myElement, x, y, z, iBound )
 ! S/R GetBoundaryLocation
 !
 ! =============================================================================================== !
 ! DECLARATIONS
  IMPLICIT NONE
  CLASS(HexElement), INTENT(in) :: myElement
  INTEGER, INTENT(in)           :: iBound
  REAL(prec), INTENT(out)       :: x(0:myElement % geometry % nMax, 0:myElement % geometry % nMax)
  REAL(prec), INTENT(out)       :: y(0:myElement % geometry % nMax, 0:myElement % geometry % nMax)
  REAL(prec), INTENT(out)       :: z(0:myElement % geometry % nMax, 0:myElement % geometry % nMax)
  
     CALL myElement % geometry % GetBoundaryLocation( x, y, z, iBound )
     
 END SUBROUTINE GetBoundaryLocation_HexElement
!
!
!
 SUBROUTINE SetBoundaryLocationAtNode_HexElement( myElement, x, y, z, i, j, iBound )
 ! S/R SetBoundaryLocationAtNode
 !
 ! =============================================================================================== !
 ! DECLARATIONS
  IMPLICIT NONE
  CLASS(HexElement), INTENT(inout) :: myElement
  REAL(prec), INTENT(in)            :: x, y, z
  INTEGER, INTENT(in)               :: i, j, iBound
  
     CALL myElement % geometry % SetBoundaryLocationAtNode( x, y, z, i, j, iBound )
     
 END SUBROUTINE SetBoundaryLocationAtNode_HexElement
!
!
!
 SUBROUTINE GetBoundaryLocationAtNode_HexElement( myElement, x, y, z, i, j, iBound )
 ! S/R GetBoundaryLocationAtNode
 !
 ! =============================================================================================== !
 ! DECLARATIONS
  IMPLICIT NONE
  CLASS(HexElement), INTENT(in) :: myElement
  REAL(prec), INTENT(out)       :: x, y, z
  INTEGER, INTENT(in)           :: i, j, iBound
  
     CALL myElement % geometry % GetBoundaryLocationAtNode( x, y, z, i, j, iBound )
     
 END SUBROUTINE GetBoundaryLocationAtNode_HexElement
!
!
!
 SUBROUTINE SetBoundaryNormalAtNode_HexElement( myElement, dir, i, j, iBound )
 ! S/R SetBoundaryNormalAtNode
 !
 ! =============================================================================================== !
 ! DECLARATIONS
  IMPLICIT NONE
  CLASS(HexElement), INTENT(inout) :: myElement
  REAL(prec), INTENT(in)           :: dir(1:nDims)
  INTEGER, INTENT(in)              :: i, j, iBound
  
     CALL myElement % geometry % SetBoundaryNormalAtNode( dir, i, j, iBound )
     
 END SUBROUTINE SetBoundaryNormalAtNode_HexElement
!
!
!
 SUBROUTINE GetBoundaryNormalAtNode_HexElement( myElement, dir, length, i, j, iBound )
 ! S/R GetBoundaryNormalAtNode
 !
 ! =============================================================================================== !
 ! DECLARATIONS
  IMPLICIT NONE
  CLASS(HexElement), INTENT(in) :: myElement
  REAL(prec), INTENT(out)       :: dir(1:nDims), length
  INTEGER, INTENT(in)           :: i, j, iBound
  
     CALL myElement % geometry % GetBoundaryNormalAtNode( dir, length, i, j, iBound )
     
 END SUBROUTINE GetBoundaryNormalAtNode_HexElement
!
!
!
 SUBROUTINE ScaleGeometry_HexElement( myElement, interp, xScale, yScale, zScale )
 !
 !
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS( HexElement ), INTENT(inout) :: myElement
   TYPE( Lagrange_3D ), INTENT(in)    :: interp
   REAL(prec), INTENT(in)             :: xScale, yScale, zScale
   
         CALL myElement % geometry % ScaleGeometry( interp, xScale, yScale, zScale )
         
 END SUBROUTINE ScaleGeometry_HexElement
!
!
!==================================================================================================!
!--------------------------------- Type Specific Routines -----------------------------------------!
!==================================================================================================!
!
!
 SUBROUTINE IsInside_HexElement( myElement, interp, x, y, z, s, p, q, isInElement )
 ! S/R IsInside
 !
 !    This routine determines if the point (x,y) is contained within the element "myElement".
 !    A logical is returned which indicates whether or not the point is contained within the 
 !    element. Additionally, the computational coordinates are returned.
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS( HexElement ), INTENT(in) :: myElement
   TYPE( Lagrange_3D ), INTENT(in) :: interp
   REAL(prec), INTENT(in)          :: x, y, z
   REAL(prec),INTENT(out)          :: s, p, q
   LOGICAL, INTENT(out)            :: isInElement
   ! LOCAL
   REAL(prec) :: thisS, thisP, thisQ
   LOGICAL    :: successful
   
      CALL myElement % geometry % CalculateComputationalCoordinates( interp, x, y, z, &
                                                                     thisS, thisP, thisQ, &
                                                                     successful )
      s = thisS
      p = thisP
      q = thisQ
      
      IF( successful )THEN
         IF( ABS(thisS) <= ONE .AND. ABS(thisP) <= ONE .AND. ABS(thisQ) <= ONE )THEN
            isInElement = .TRUE.
         ELSE
            isInElement = .FALSE.
         ENDIF
      ELSE
         isInElement = .FALSE.
      ENDIF
      
      
 END SUBROUTINE IsInside_HexElement
!
!
!
END MODULE HexElementClass
