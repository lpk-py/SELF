! QuadElementClass.f90
! 
! Copyright 2015 Joseph Schoonover <schoonover.numerics@gmail.com>
! 
! QuadElementClass.f90 is part of the Spectral Element Libraries in Fortran (SELF).
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
 
 
MODULE QuadElementClass
! ========================================= Logs ================================================= !
!2016-05-11  Joseph Schoonover  schoonover.numerics@gmail.com 
!
! //////////////////////////////////////////////////////////////////////////////////////////////// ! 
! src/common/
USE ModelPrecision
USE ConstantsDictionary
USE ModelFlags
! src/interp/
USE Lagrange_2D_Class
! src/geom/
USE CurveClass_2D
USE GeometrySupportRoutines
USE MappedGeometryClass_2D


IMPLICIT NONE

   TYPE QuadElement
      INTEGER                         :: nS, nP, nmax
      INTEGER, PRIVATE                :: nodeIDs(1:4)   ! Corner Node ID's
      INTEGER, PRIVATE                :: neighbors(1:4) ! Elements IDs for the neighbors
      INTEGER, PRIVATE                :: globElID       ! Global Element ID 
      TYPE(MappedGeometry_2D)         :: geometry       ! Contains the element metric terms, etc.


      CONTAINS

      PROCEDURE :: Build => Build_QuadElement
      PROCEDURE :: Trash => Trash_QuadElement
    
      PROCEDURE :: SetNodeIDs => SetNodeIDs_QuadElement
      PROCEDURE :: GetNodeIDs => GetNodeIDs_QuadElement
      PROCEDURE :: SetNodeID => SetNodeID_QuadElement
      PROCEDURE :: GetNodeID => GetNodeID_QuadElement
      PROCEDURE :: SetElementID => SetElementID_QuadElement
      PROCEDURE :: GetElementID => GetElementID_QuadElement
      PROCEDURE :: SetNeighbor => SetNeighbor_QuadElement
      PROCEDURE :: GetNeighbor => GetNeighbor_QuadElement
      PROCEDURE :: SetSouthernNeighbor => SetSouthernNeighbor_QuadElement
      PROCEDURE :: GetSouthernNeighbor => GetSouthernNeighbor_QuadElement
      PROCEDURE :: SetNorthernNeighbor => SetNorthernNeighbor_QuadElement
      PROCEDURE :: GetNorthernNeighbor => GetNorthernNeighbor_QuadElement
      PROCEDURE :: SetEasternNeighbor  => SetEasternNeighbor_QuadElement
      PROCEDURE :: GetEasternNeighbor  => GetEasternNeighbor_QuadElement
      PROCEDURE :: SetWesternNeighbor  => SetWesternNeighbor_QuadElement
      PROCEDURE :: GetWesternNeighbor  => GetWesternNeighbor_QuadElement
      
      ! MappedGeometry_2D Wrapper Routines
      PROCEDURE :: SetNumberOfNodes => SetNumberOfNodes_QuadElement
      PROCEDURE :: GetNumberOfNodes => GetNumberOfNodes_QuadElement
      PROCEDURE :: SetPositions => SetPositions_QuadElement
      PROCEDURE :: GetPositions => GetPositions_QuadElement
      PROCEDURE :: SetPositionAtNode => SetPositionAtNode_QuadElement
      PROCEDURE :: GetPositionAtNode => GetPositionAtNode_QuadElement
      PROCEDURE :: SetJacobian => SetJacobian_QuadElement
      PROCEDURE :: GetJacobian => GetJacobian_QuadElement
      PROCEDURE :: SetJacobianAtNode => SetJacobianAtNode_QuadElement
      PROCEDURE :: GetJacobianAtNode => GetJacobianAtNode_QuadElement
      PROCEDURE :: SetCovariantMetrics => SetCovariantMetrics_QuadElement
      PROCEDURE :: GetCovariantMetrics => GetCovariantMetrics_QuadElement
      PROCEDURE :: SetCovariantMetricsAtNode => SetCovariantMetricsAtNode_QuadElement
      PROCEDURE :: GetCovariantMetricsAtNode => GetCovariantMetricsAtNode_QuadElement
      PROCEDURE :: SetBoundaryLocation => SetBoundaryLocation_QuadElement
      PROCEDURE :: GetBoundaryLocation => GetBoundaryLocation_QuadElement
      PROCEDURE :: SetBoundaryLocationAtNode => SetBoundaryLocationAtNode_QuadElement
      PROCEDURE :: GetBoundaryLocationAtNode => GetBoundaryLocationAtNode_QuadElement
      PROCEDURE :: SetBoundaryNormalAtNode => SetBoundaryNormalAtNode_QuadElement
      PROCEDURE :: GetBoundaryNormalAtNode => GetBoundaryNormalAtNode_QuadElement
      PROCEDURE :: ScaleGeometry => ScaleGeometry_QuadElement
      
      PROCEDURE :: IsInside => IsInside_QuadElement
      
   END TYPE QuadElement



 CONTAINS
!
!
!==================================================================================================!
!------------------------------- Manual Constructors/Destructors ----------------------------------!
!==================================================================================================!
!
!
 SUBROUTINE Build_QuadElement( myElement, nodeIDs, eID, boundCurves, interp  )
 ! S/R Build
 ! 
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS( QuadElement ), INTENT(inout) :: myElement
   INTEGER, INTENT(in)                 :: nodeIDs(1:4), eID
   TYPE( Curve_2D ), INTENT(in)        :: boundCurves(1:4)
   TYPE( Lagrange_2D ), INTENT(in)     :: interp
   
      CALL myElement % SetNodeIDs( nodeIDs )

      CALL myElement % SetElementID( eID )

      CALL myElement % geometry % Build( interp, boundCurves  )
      
      myElement % nS = interp % nS
      myElement % nP = interp % nP
      myElement % nMax = MAX( interp % nS, interp % nP )

 END SUBROUTINE Build_QuadElement
!
!
!
 SUBROUTINE Trash_QuadElement( myElement )
 ! S/R Trash
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS( QuadElement ), INTENT(inout) :: myElement

      CALL myElement % geometry % Trash( )

 END SUBROUTINE Trash_QuadElement
!
!
!==================================================================================================!
!--------------------------------------- Accessors ------------------------------------------------!
!==================================================================================================!
!
!
 SUBROUTINE SetNodeIDs_QuadElement( myElement, nodeIDs )
 ! S/R SetNodeIDs
 !
 ! =============================================================================================== !
   IMPLICIT NONE
   CLASS( QuadElement ), INTENT(inout) :: myElement
   INTEGER, INTENT(in)                 :: nodeIDs(1:4)
   
      myElement % nodeIDs = nodeIDs

 END SUBROUTINE SetNodeIDs_QuadElement
!
!
!
 SUBROUTINE GetNodeIDs_QuadElement( myElement, nodeIDs )
 ! S/R GetNodeIDs
 !
 ! =============================================================================================== !
   IMPLICIT NONE
   CLASS( QuadElement ), INTENT(in) :: myElement
   INTEGER, INTENT(out)             :: nodeIDs(1:4)
   
      nodeIDs = myElement % nodeIDs

 END SUBROUTINE GetNodeIDs_QuadElement
!
!
!
 SUBROUTINE SetNodeID_QuadElement( myElement, localID, nodeID )
 ! S/R SetNodeID
 !
 ! =============================================================================================== !
   IMPLICIT NONE
   CLASS( QuadElement ), INTENT(inout) :: myElement
   INTEGER, INTENT(in)                 :: localID, nodeID
   
      myElement % nodeIDs(localID) = nodeID

 END SUBROUTINE SetNodeID_QuadElement
!
!
!
 SUBROUTINE GetNodeID_QuadElement( myElement, localID, nodeID )
 ! S/R GetNodeID
 !
 ! =============================================================================================== !
   IMPLICIT NONE
   CLASS( QuadElement ), INTENT(in) :: myElement
   INTEGER, INTENT(in)              :: localID
   INTEGER, INTENT(out)             :: nodeID
   
      nodeID = myElement % nodeIDs(localID)

 END SUBROUTINE GetNodeID_QuadElement
!
!
!
 SUBROUTINE SetElementID_QuadElement( myElement, eID )
 ! S/R SetElementID
 !
 ! =============================================================================================== !
   IMPLICIT NONE
   CLASS( QuadElement ), INTENT(inout) :: myElement
   INTEGER, INTENT(in)                 :: eID
   
      myElement % globElID = eID

 END SUBROUTINE SetElementID_QuadElement
!
!
!
 SUBROUTINE GetElementID_QuadElement( myElement, eID )
 ! S/R GetElementID
 !
 ! =============================================================================================== !
   IMPLICIT NONE
   CLASS( QuadElement ), INTENT(in) :: myElement
   INTEGER, INTENT(out)             :: eID
   
      eID = myElement % globElID

 END SUBROUTINE GetElementID_QuadElement
!
!
!
 SUBROUTINE SetNeighbor_QuadElement( myElement, sID, eID )
  ! S/R SetNeighbor
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS( QuadElement ), INTENT(inout) :: myElement
   INTEGER, INTENT(in)                 :: sID
   INTEGER, INTENT(in)                 :: eID
   
      myElement % neighbors(sID) = eID
      
 END SUBROUTINE SetNeighbor_QuadElement
!
!
!
 SUBROUTINE GetNeighbor_QuadElement( myElement, sID, eID )
 ! S/R GetNeighbor
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS( QuadElement ), INTENT(in) :: myElement
   INTEGER, INTENT(in)              :: sID
   INTEGER, INTENT(out)             :: eID
   
      eID = myElement % neighbors(sID)
      
 END SUBROUTINE GetNeighbor_QuadElement
!
!
!
 SUBROUTINE SetSouthernNeighbor_QuadElement( myElement, eID )
 ! S/R SetSouthernNeighbor
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS( QuadElement ), INTENT(inout) :: myElement
   INTEGER, INTENT(in)                 :: eID
   
      myElement % neighbors(1) = eID
      
 END SUBROUTINE SetSouthernNeighbor_QuadElement
!
!
!
 SUBROUTINE GetSouthernNeighbor_QuadElement( myElement, eID )
 ! S/R GetSouthernNeighbor
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS( QuadElement ), INTENT(in) :: myElement
   INTEGER, INTENT(out)             :: eID
   
      eID = myElement % neighbors(1)
      
 END SUBROUTINE GetSouthernNeighbor_QuadElement
!
!
!
 SUBROUTINE SetNorthernNeighbor_QuadElement( myElement, eID )
 ! S/R SetNorthernNeighbor
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS( QuadElement ), INTENT(inout) :: myElement
   INTEGER, INTENT(in)                 :: eID
   
      myElement % neighbors(3) = eID
      
 END SUBROUTINE SetNorthernNeighbor_QuadElement
!
!
!
 SUBROUTINE GetNorthernNeighbor_QuadElement( myElement, eID )
 ! S/R GetNorthernNeighbor
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS( QuadElement ), INTENT(in) :: myElement
   INTEGER, INTENT(out)             :: eID
   
      eID = myElement % neighbors(3)
      
 END SUBROUTINE GetNorthernNeighbor_QuadElement
!
!
!
 SUBROUTINE SetEasternNeighbor_QuadElement( myElement, eID )
 ! S/R SetEasternNeighbor
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS( QuadElement ), INTENT(inout) :: myElement
   INTEGER, INTENT(in)                 :: eID
   
      myElement % neighbors(2) = eID
      
 END SUBROUTINE SetEasternNeighbor_QuadElement
!
!
!
 SUBROUTINE GetEasternNeighbor_QuadElement( myElement, eID )
 ! S/R GetEasternNeighbor
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS( QuadElement ), INTENT(in) :: myElement
   INTEGER, INTENT(out)             :: eID
   
      eID = myElement % neighbors(2)
      
 END SUBROUTINE GetEasternNeighbor_QuadElement
!
!
!
 SUBROUTINE SetWesternNeighbor_QuadElement( myElement, eID )
 ! S/R SetWesternNeighbor
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS( QuadElement ), INTENT(inout) :: myElement
   INTEGER, INTENT(in)                 :: eID
   
      myElement % neighbors(4) = eID
      
 END SUBROUTINE SetWesternNeighbor_QuadElement
!
!
!
 SUBROUTINE GetWesternNeighbor_QuadElement( myElement, eID )
 ! S/R GetWesternNeighbor
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS( QuadElement ), INTENT(in) :: myElement
   INTEGER, INTENT(out)             :: eID
   
      eID = myElement % neighbors(4)
      
 END SUBROUTINE GetWesternNeighbor_QuadElement
!
! ---------------------------- MappedGeometryClass_2D Wrappers ----------------------------------- !
!
SUBROUTINE SetNumberOfNodes_QuadElement( myElement, nS, nP )
 ! S/R SetNumberOfNodes
 !
 ! =============================================================================================== !
 ! DECLARATIONS
  IMPLICIT NONE
  CLASS(QuadElement), INTENT(inout) :: myElement
  INTEGER, INTENT(in)               :: nS, nP
  
     CALL myElement % geometry % SetNumberOfNodes( nS, nP )
     
 END SUBROUTINE SetNumberOfNodes_QuadElement
!
!
!
 SUBROUTINE GetNumberOfNodes_QuadElement( myElement, nS, nP )
 ! S/R GetNumberOfNodes
 !
 ! =============================================================================================== !
 ! DECLARATIONS
  IMPLICIT NONE
  CLASS(QuadElement), INTENT(in) :: myElement
  INTEGER, INTENT(out)           :: nS, nP
  
     CALL myElement % geometry % GetNumberOfNodes( nS, nP )
     
 END SUBROUTINE GetNumberOfNodes_QuadElement
!
!
!
 SUBROUTINE SetPositions_QuadElement( myElement, x, y )
 ! S/R SetPositions
 !
 ! =============================================================================================== !
 ! DECLARATIONS
  IMPLICIT NONE
  CLASS(QuadElement), INTENT(inout) :: myElement
  REAL(prec), INTENT(in)            :: x(0:myElement % geometry % nS, 0:myElement % geometry % nP)
  REAL(prec), INTENT(in)            :: y(0:myElement % geometry % nS, 0:myElement % geometry % nP)
  
     CALL myElement % geometry % SetPositions( x, y )
     
 END SUBROUTINE SetPositions_QuadElement
!
!
!
 SUBROUTINE GetPositions_QuadElement( myElement, x, y )
 ! S/R GetPositions
 !
 ! =============================================================================================== !
 ! DECLARATIONS
  IMPLICIT NONE
  CLASS(QuadElement), INTENT(in) :: myElement
  REAL(prec), INTENT(out)        :: x(0:myElement % geometry % nS, 0:myElement % geometry % nP)
  REAL(prec), INTENT(out)        :: y(0:myElement % geometry % nS, 0:myElement % geometry % nP)
  
     CALL myElement % geometry % GetPositions( x, y )
     
 END SUBROUTINE GetPositions_QuadElement
!
!
!
 SUBROUTINE SetPositionAtNode_QuadElement( myElement, x, y, i, j )
 ! S/R SetPositionAtNode
 !
 ! =============================================================================================== !
 ! DECLARATIONS
  IMPLICIT NONE
  CLASS(QuadElement), INTENT(inout) :: myElement
  REAL(prec), INTENT(in)            :: x, y
  INTEGER, INTENT(in)               :: i, j
  
     CALL myElement % geometry % SetPositionAtNode( x, y, i, j )
     
 END SUBROUTINE SetPositionAtNode_QuadElement
!
!
!
 SUBROUTINE GetPositionAtNode_QuadElement( myElement, x, y, i, j )
 ! S/R GetPositionAtNode
 !
 ! =============================================================================================== !
 ! DECLARATIONS
  IMPLICIT NONE
  CLASS(QuadElement), INTENT(in) :: myElement
  REAL(prec), INTENT(out)        :: x, y
  INTEGER, INTENT(in)            :: i, j
  
     CALL myElement % geometry % GetPositionAtNode( x, y, i, j )
     
 END SUBROUTINE GetPositionAtNode_QuadElement
!
!
!
 SUBROUTINE SetJacobian_QuadElement( myElement, J )
 ! S/R SetJacobian
 !
 ! =============================================================================================== !
 ! DECLARATIONS
  IMPLICIT NONE
  CLASS(QuadElement), INTENT(inout) :: myElement
  REAL(prec), INTENT(in)            :: J(0:myElement % geometry % nS, 0:myElement % geometry % nP)
  
     CALL myElement % geometry % SetJacobian( J )
     
 END SUBROUTINE SetJacobian_QuadElement
!
!
!
 SUBROUTINE GetJacobian_QuadElement( myElement, J )
 ! S/R GetJacobian
 !
 ! =============================================================================================== !
 ! DECLARATIONS
  IMPLICIT NONE
  CLASS(QuadElement), INTENT(in) :: myElement
  REAL(prec), INTENT(out)        :: J(0:myElement % geometry % nS, 0:myElement % geometry % nP)
  
     CALL myElement % geometry % GetJacobian( J )
     
 END SUBROUTINE GetJacobian_QuadElement
!
!
!
 SUBROUTINE SetJacobianAtNode_QuadElement( myElement, J, iS, iP )
 ! S/R SetJacobianAtNode
 !
 ! =============================================================================================== !
 ! DECLARATIONS
  IMPLICIT NONE
  CLASS(QuadElement), INTENT(inout) :: myElement
  REAL(prec), INTENT(in)            :: J
  INTEGER, INTENT(in)               :: iS, iP
  
     CALL myElement % geometry % SetJacobianAtNode( J, iS, iP )
     
 END SUBROUTINE SetJacobianAtNode_QuadElement
!
!
!
 SUBROUTINE GetJacobianAtNode_QuadElement( myElement, J, iS, iP )
 ! S/R GetJacobianAtNode
 !
 ! =============================================================================================== !
 ! DECLARATIONS
  IMPLICIT NONE
  CLASS(QuadElement), INTENT(in) :: myElement
  REAL(prec), INTENT(out)        :: J
  INTEGER, INTENT(in)            :: iS, iP
  
     CALL myElement % geometry % GetJacobianAtNode( J, iS, iP )
     
 END SUBROUTINE GetJacobianAtNode_QuadElement
!
!
!
 SUBROUTINE SetCovariantMetrics_QuadElement( myElement, dxds, dxdp, dyds, dydp )
 ! S/R SetCovariantMetrics
 !
 ! =============================================================================================== !
 ! DECLARATIONS
  IMPLICIT NONE
  CLASS(QuadElement), INTENT(inout) :: myElement
  REAL(prec), INTENT(in)            :: dxds(0:myElement % geometry % nS, 0:myElement % geometry % nP)
  REAL(prec), INTENT(in)            :: dxdp(0:myElement % geometry % nS, 0:myElement % geometry % nP)
  REAL(prec), INTENT(in)            :: dyds(0:myElement % geometry % nS, 0:myElement % geometry % nP)
  REAL(prec), INTENT(in)            :: dydp(0:myElement % geometry % nS, 0:myElement % geometry % nP)
  
     CALL myElement % geometry % SetCovariantMetrics( dxds, dxdp, dyds, dydp )
     
 END SUBROUTINE SetCovariantMetrics_QuadElement
!
!
!
 SUBROUTINE GetCovariantMetrics_QuadElement( myElement, dxds, dxdp, dyds, dydp )
 ! S/R GetCovariantMetrics
 !
 ! =============================================================================================== !
 ! DECLARATIONS
  IMPLICIT NONE
  CLASS(QuadElement), INTENT(in) :: myElement
  REAL(prec), INTENT(out)        :: dxds(0:myElement % geometry % nS, 0:myElement % geometry % nP)
  REAL(prec), INTENT(out)        :: dxdp(0:myElement % geometry % nS, 0:myElement % geometry % nP)
  REAL(prec), INTENT(out)        :: dyds(0:myElement % geometry % nS, 0:myElement % geometry % nP)
  REAL(prec), INTENT(out)        :: dydp(0:myElement % geometry % nS, 0:myElement % geometry % nP)
  
     CALL myElement % geometry % GetCovariantMetrics( dxds, dxdp, dyds, dydp )
     
 END SUBROUTINE GetCovariantMetrics_QuadElement
!
!
!
 SUBROUTINE SetCovariantMetricsAtNode_QuadElement( myElement, dxds, dxdp, dyds, dydp, iS, iP )
 ! S/R SetCovariantMetricsAtNode
 !
 ! =============================================================================================== !
 ! DECLARATIONS
  IMPLICIT NONE
  CLASS(QuadElement), INTENT(inout) :: myElement
  REAL(prec), INTENT(in)            :: dxds, dxdp, dyds, dydp
  INTEGER, INTENT(in)               :: iS, iP
  
     CALL myElement % geometry % SetCovariantMetricsAtNode( dxds, dxdp, dyds, dydp, iS, iP )
     
 END SUBROUTINE SetCovariantMetricsAtNode_QuadElement
!
!
!
 SUBROUTINE GetCovariantMetricsAtNode_QuadElement( myElement, dxds, dxdp, dyds, dydp, iS, iP )
 ! S/R GetCovariantMetricsAtNode
 !
 ! =============================================================================================== !
 ! DECLARATIONS
  IMPLICIT NONE
  CLASS(QuadElement), INTENT(in) :: myElement
  REAL(prec), INTENT(out)        :: dxds, dxdp, dyds, dydp
  INTEGER, INTENT(in)            :: iS, iP
  
     CALL myElement % geometry % GetCovariantMetricsAtNode( dxds, dxdp, dyds, dydp, iS, iP )
     
 END SUBROUTINE GetCovariantMetricsAtNode_QuadElement
!
!
!
 SUBROUTINE SetBoundaryLocation_QuadElement( myElement, x, y, iBound )
 ! S/R SetBoundaryLocation
 !
 ! =============================================================================================== !
 ! DECLARATIONS
  IMPLICIT NONE
  CLASS(QuadElement), INTENT(inout) :: myElement
  INTEGER, INTENT(in)               :: iBound
  REAL(prec), INTENT(in)            :: x(0:myElement % geometry % nMax)
  REAL(prec), INTENT(in)            :: y(0:myElement % geometry % nMax)
  
     CALL myElement % geometry % SetBoundaryLocation( x, y, iBound )
     
 END SUBROUTINE SetBoundaryLocation_QuadElement
!
!
!
 SUBROUTINE GetBoundaryLocation_QuadElement( myElement, x, y, iBound )
 ! S/R GetBoundaryLocation
 !
 ! =============================================================================================== !
 ! DECLARATIONS
  IMPLICIT NONE
  CLASS(QuadElement), INTENT(in) :: myElement
  INTEGER, INTENT(in)            :: iBound
  REAL(prec), INTENT(out)        :: x(0:myElement % geometry % nMax)
  REAL(prec), INTENT(out)        :: y(0:myElement % geometry % nMax)
  
     CALL myElement % geometry % GetBoundaryLocation( x, y, iBound )
     
 END SUBROUTINE GetBoundaryLocation_QuadElement
!
!
!
 SUBROUTINE SetBoundaryLocationAtNode_QuadElement( myElement, x, y, i, iBound )
 ! S/R SetBoundaryLocationAtNode
 !
 ! =============================================================================================== !
 ! DECLARATIONS
  IMPLICIT NONE
  CLASS(QuadElement), INTENT(inout) :: myElement
  REAL(prec), INTENT(in)            :: x, y
  INTEGER, INTENT(in)               :: i, iBound
  
     CALL myElement % geometry % SetBoundaryLocationAtNode( x, y, i, iBound )
     
 END SUBROUTINE SetBoundaryLocationAtNode_QuadElement
!
!
!
 SUBROUTINE GetBoundaryLocationAtNode_QuadElement( myElement, x, y, i, iBound )
 ! S/R GetBoundaryLocationAtNode
 !
 ! =============================================================================================== !
 ! DECLARATIONS
  IMPLICIT NONE
  CLASS(QuadElement), INTENT(in) :: myElement
  REAL(prec), INTENT(out)        :: x, y
  INTEGER, INTENT(in)            :: i, iBound
  
     CALL myElement % geometry % GetBoundaryLocationAtNode( x, y, i, iBound )
     
 END SUBROUTINE GetBoundaryLocationAtNode_QuadElement
!
!
!
 SUBROUTINE SetBoundaryNormalAtNode_QuadElement( myElement, dir, i, iBound )
 ! S/R SetBoundaryNormalAtNode
 !
 ! =============================================================================================== !
 ! DECLARATIONS
  IMPLICIT NONE
  CLASS(QuadElement), INTENT(inout) :: myElement
  REAL(prec), INTENT(in)            :: dir(1:2)
  INTEGER, INTENT(in)               :: i, iBound
  
     CALL myElement % geometry % SetBoundaryNormalAtNode( dir, i, iBound )
     
 END SUBROUTINE SetBoundaryNormalAtNode_QuadElement
!
!
!
 SUBROUTINE GetBoundaryNormalAtNode_QuadElement( myElement, dir, length, i, iBound )
 ! S/R GetBoundaryNormalAtNode
 !
 ! =============================================================================================== !
 ! DECLARATIONS
  IMPLICIT NONE
  CLASS(QuadElement), INTENT(in) :: myElement
  REAL(prec), INTENT(out)        :: dir(1:2), length
  INTEGER, INTENT(in)            :: i, iBound
  
     CALL myElement % geometry % GetBoundaryNormalAtNode( dir, length, i, iBound )
     
 END SUBROUTINE GetBoundaryNormalAtNode_QuadElement
!
!
!
 SUBROUTINE ScaleGeometry_QuadElement( myElement, interp, xScale, yScale )
 !
 !
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS( QuadElement ), INTENT(inout) :: myElement
   TYPE( Lagrange_2D ), INTENT(in)     :: interp
   REAL(prec), INTENT(in)              :: xScale, yScale
   
         CALL myElement % geometry % ScaleGeometry( interp, xScale, yScale )
         
 END SUBROUTINE ScaleGeometry_QuadElement
!
!
!==================================================================================================!
!--------------------------------- Type Specific Routines -----------------------------------------!
!==================================================================================================!
!
!
 SUBROUTINE IsInside_QuadElement( myElement, interp, x, y, s, p, isInElement )
 ! S/R IsInside
 !
 !    This routine determines if the point (x,y) is contained within the element "myElement".
 !    A logical is returned which indicates whether or not the point is contained within the 
 !    element. Additionally, the computational coordinates are returned.
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS( QuadElement ), INTENT(in) :: myElement
   TYPE( Lagrange_2D ), INTENT(in)  :: interp
   REAL(prec), INTENT(in)           :: x, y
   REAL(prec),INTENT(out)           :: s, p
   LOGICAL, INTENT(out)             :: isInElement
   ! LOCAL
   REAL(prec) :: thisS, thisP
   LOGICAL    :: successful
   
      CALL myElement % geometry % CalculateComputationalCoordinates( interp, x, y, thisS, thisP, successful )
      s = thisS
      p = thisP
      
      IF( successful )THEN
         IF( ABS(thisS) <= ONE .AND. ABS(thisP) <= ONE )THEN
            isInElement = .TRUE.
         ELSE
            isInElement = .FALSE.
         ENDIF
      ELSE
         isInElement = .FALSE.
      ENDIF
      
      
 END SUBROUTINE IsInside_QuadElement
!
!
!
END MODULE QuadElementClass
