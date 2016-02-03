! MappedGeometryClass_3D.f90 ( v2.1 - 14 Dec. 2015)
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
! o (ver 1.0) March 2014 
! o (ver 2.1) December 2015
!
! ========================================= Logs ================================================= !
! yyyy-mm-dd  Name  <user@computer> : e-mail address
!
! 
!   
! 
! //////////////////////////////////////////////////////////////////////////////////////////////// !

MODULE MappedGeometryClass_3D
! MappedGeometryClass_3D.f90 (v2.1 - 14 Dec. 2015) 
! 
! schoonover.numerics@gmail.com
! 
! o (ver 1.0) March 2014
! o (ver 2.1) Dec 2015
!
!
! 
!  
! ================================================================================================ !
! src/common/
USE ModelPrecision
USE ConstantsDictionary
USE ModelFlags
USE CommonRoutines
! src/interp/
USE Lagrange_3D_Class
! src/geom/
USE SurfaceClass_3D
USE VectorClass



IMPLICIT NONE



   TYPE MappedGeometry_3D
      INTEGER                             :: nS, nP, nQ, nMax
      TYPE(Vector), ALLOCATABLE, PRIVATE  :: nHat(:,:,:) 
      REAL(prec), ALLOCATABLE, PRIVATE    :: xBound(:,:,:)
      REAL(prec), ALLOCATABLE, PRIVATE    :: yBound(:,:,:) 
      REAL(prec), ALLOCATABLE, PRIVATE    :: zBound(:,:,:) 
      REAL(prec), ALLOCATABLE, PRIVATE    :: x(:,:,:), y(:,:,:), z(:,:,:)
      REAL(prec), ALLOCATABLE, PRIVATE    :: J(:,:,:)    
      REAL(prec), ALLOCATABLE, PRIVATE    :: dxds(:,:,:), dxdp(:,:,:), dxdq(:,:,:)
      REAL(prec), ALLOCATABLE, PRIVATE    :: dyds(:,:,:), dydp(:,:,:), dydq(:,:,:)
      REAL(prec), ALLOCATABLE, PRIVATE    :: dzds(:,:,:), dzdp(:,:,:), dzdq(:,:,:)

      CONTAINS

      PROCEDURE :: Build => Build_MappedGeometry_3D
      PROCEDURE :: Trash => Trash_MappedGeometry_3D
      
      PROCEDURE :: SetNumberOfNodes => SetNumberOfNodes_MappedGeometry_3D
      PROCEDURE :: GetNumberOfNodes => GetNumberOfNodes_MappedGeometry_3D
      PROCEDURE :: SetPositions => SetPositions_MappedGeometry_3D
      PROCEDURE :: GetPositions => GetPositions_MappedGeometry_3D
      PROCEDURE :: SetPositionAtNode => SetPositionAtNode_MappedGeometry_3D
      PROCEDURE :: GetPositionAtNode => GetPositionAtNode_MappedGeometry_3D
      PROCEDURE :: SetJacobian => SetJacobian_MappedGeometry_3D
      PROCEDURE :: GetJacobian => GetJacobian_MappedGeometry_3D
      PROCEDURE :: SetJacobianAtNode => SetJacobianAtNode_MappedGeometry_3D
      PROCEDURE :: GetJacobianAtNode => GetJacobianAtNode_MappedGeometry_3D
      PROCEDURE :: SetCovariantMetrics => SetCovariantMetrics_MappedGeometry_3D
      PROCEDURE :: GetCovariantMetrics => GetCovariantMetrics_MappedGeometry_3D
      PROCEDURE :: SetCovariantMetricsAtNode => SetCovariantMetricsAtNode_MappedGeometry_3D
      PROCEDURE :: GetCovariantMetricsAtNode => GetCovariantMetricsAtNode_MappedGeometry_3D
      PROCEDURE :: SetBoundaryLocation => SetBoundaryLocation_MappedGeometry_3D
      PROCEDURE :: GetBoundaryLocation => GetBoundaryLocation_MappedGeometry_3D
      PROCEDURE :: SetBoundaryLocationAtNode => SetBoundaryLocationAtNode_MappedGeometry_3D
      PROCEDURE :: GetBoundaryLocationAtNode => GetBoundaryLocationAtNode_MappedGeometry_3D
      PROCEDURE :: SetBoundaryNormalAtNode => SetBoundaryNormalAtNode_MappedGeometry_3D
      PROCEDURE :: GetBoundaryNormalAtNode => GetBoundaryNormalAtNode_MappedGeometry_3D

      PROCEDURE :: GenerateMesh => GenerateMesh_MappedGeometry_3D
      PROCEDURE :: GenerateBoundaryPositions => GenerateBoundaryPositions_MappedGeometry_3D 
      PROCEDURE :: GenerateMetrics => GenerateMetrics_MappedGeometry_3D
      PROCEDURE :: GenerateBoundaryMetrics => GenerateBoundaryMetrics_MappedGeometry_3D
      PROCEDURE :: ScaleGeometry => ScaleGeometry_MappedGeometry_3D
      PROCEDURE :: CalculateLocation => CalculateLocation_MappedGeometry_3D
      PROCEDURE :: CalculateMetrics => CalculateMetrics_MappedGeometry_3D
      PROCEDURE :: CalculateComputationalCoordinates => CalculateComputationalCoordinates_MappedGeometry_3D
      
      PROCEDURE :: WriteTecplot => WriteTecplot_MappedGeometry_3D
   END TYPE MappedGeometry_3D


 INTEGER, PRIVATE :: nDims = 3
 
 CONTAINS
!
!
!==================================================================================================!
!------------------------------- Manual Constructors/Destructors ----------------------------------!
!==================================================================================================!
!
!
 SUBROUTINE Build_MappedGeometry_3D( myGeom, interp, mySurfaces )
 ! S/R Build
 !
 ! =============================================================================================== !
 ! DECLARATIONS
  IMPLICIT NONE
  CLASS(MappedGeometry_3D), INTENT(out) :: myGeom
  TYPE(Lagrange_3D), INTENT(in)         :: interp
  TYPE(Surface_3D), INTENT(in)          :: mySurfaces(1:nHexFaces)
  !LOCAL
  INTEGER :: nS, nP, nQ, i, j
   

      CALL interp % GetNumberOfNodes( nS, nP, nQ )
      CALL myGeom % SetNumberOfNodes( nS, nP, nQ )
       
      ! Allocate space
      ALLOCATE( myGeom % dxds(0:nS,0:nP,0:nQ), myGeom % dxdp(0:nS,0:nP,0:nQ), myGeom % dxdq(0:nS,0:nP,0:nQ) )
      ALLOCATE( myGeom % dyds(0:nS,0:nP,0:nQ), myGeom % dydp(0:nS,0:nP,0:nQ), myGeom % dydq(0:nS,0:nP,0:nQ) )
      ALLOCATE( myGeom % dzds(0:nS,0:nP,0:nQ), myGeom % dzdp(0:nS,0:nP,0:nQ), myGeom % dzdq(0:nS,0:nP,0:nQ) )
      ALLOCATE( myGeom % J(0:nS,0:nP,0:nQ) )
      ALLOCATE( myGeom % x(0:nS,0:nP,0:nQ), myGeom % y(0:nS,0:nP,0:nQ), myGeom % z(0:nS,0:nP,0:nQ) )
      ALLOCATE( myGeom % xBound(0:max(nS,nP,nQ),0:max(nS,nP,nQ),1:nHexFaces) )
      ALLOCATE( myGeom % yBound(0:max(nS,nP,nQ),0:max(nS,nP,nQ),1:nHexFaces) )
      ALLOCATE( myGeom % zBound(0:max(nS,nP,nQ),0:max(nS,nP,nQ),1:nHexFaces) )
      ALLOCATE( myGeom % nHat(0:max(nS,nP,nQ),0:max(nS,nP,nQ),1:nHexFaces) )
     
      DO j = 0,max(nS,nP,nQ)
         DO i = 0,max(nS,nP,nQ)
            CALL myGeom % nHat(i,j,south)  % Build( nDims ) ! Allocate space for the south nHats
            CALL myGeom % nHat(i,j,east)   % Build( nDims ) ! Allocate space for the east nHats
            CALL myGeom % nHat(i,j,north)  % Build( nDims ) ! Allocate space for the north nHats
            CALL myGeom % nHat(i,j,west)   % Build( nDims ) ! Allocate space for the west nHats
            CALL myGeom % nHat(i,j,bottom) % Build( nDims ) ! Allocate space for the bottom nHats
            CALL myGeom % nHat(i,j,top)    % Build( nDims ) ! Allocate space for the top nHats
         ENDDO
      ENDDO
      
      ! Generate the mesh locations, and the mapping metrics

      CALL myGeom % GenerateMesh( interp, mySurfaces )
      CALL myGeom % GenerateMetrics( interp )
 
 END SUBROUTINE Build_MappedGeometry_3D
!
!
!
 SUBROUTINE Trash_MappedGeometry_3D( myGeom )
 ! S/R Trash
 !
 ! =============================================================================================== !
 ! DECLARATIONS
  IMPLICIT NONE
  CLASS(MappedGeometry_3D), INTENT(inout)  :: myGeom
  ! LOCAL
  INTEGER :: i, j

      DEALLOCATE( myGeom % dxds, myGeom % dxdp, myGeom % dxdq )
      DEALLOCATE( myGeom % dyds, myGeom % dydp, myGeom % dydq )
      DEALLOCATE( myGeom % dzds, myGeom % dzdp, myGeom % dzdq )
      DEALLOCATE( myGeom % J, myGeom % x, myGeom % y, myGeom % z )
      DEALLOCATE( myGeom % xBound, myGeom % yBound, myGeom % zBound )
     
      DO j = 0, myGeom % nMax
         DO i = 0, myGeom % nMax
            CALL myGeom % nHat(i,j,south)  % Trash( )
            CALL myGeom % nHat(i,j,east)   % Trash( )
            CALL myGeom % nHat(i,j,north)  % Trash( )
            CALL myGeom % nHat(i,j,west)   % Trash( )
            CALL myGeom % nHat(i,j,bottom) % Trash( )
            CALL myGeom % nHat(i,j,top)    % Trash( )
         ENDDO
      ENDDO

      DEALLOCATE( myGeom % nHat )

 
 END SUBROUTINE Trash_MappedGeometry_3D
!
!
!==================================================================================================!
!--------------------------------------- Accessors ------------------------------------------------!
!==================================================================================================!
!
!
 SUBROUTINE SetNumberOfNodes_MappedGeometry_3D( myGeom, nS, nP, nQ )
 ! S/R SetNumberOfNodes
 !
 ! =============================================================================================== !
 ! DECLARATIONS
  IMPLICIT NONE
  CLASS(MappedGeometry_3D), INTENT(inout) :: myGeom
  INTEGER, INTENT(in)                     :: nS, nP, nQ
  
     myGeom % nS   = nS
     myGeom % nP   = nP
     myGeom % nQ   = nQ
     myGeom % nMax = max(nS,nP,nQ)
     
 END SUBROUTINE SetNumberOfNodes_MappedGeometry_3D
!
!
!
 SUBROUTINE GetNumberOfNodes_MappedGeometry_3D( myGeom, nS, nP, nQ )
 ! S/R GetNumberOfNodes
 !
 ! =============================================================================================== !
 ! DECLARATIONS
  IMPLICIT NONE
  CLASS(MappedGeometry_3D), INTENT(in) :: myGeom
  INTEGER, INTENT(out)                 :: nS, nP, nQ
  
     nS = myGeom % nS
     nP = myGeom % nP
     nQ = myGeom % nQ
     
 END SUBROUTINE GetNumberOfNodes_MappedGeometry_3D
!
!
!
 SUBROUTINE SetPositions_MappedGeometry_3D( myGeom, x, y, z )
 ! S/R SetPositions
 !
 ! =============================================================================================== !
 ! DECLARATIONS
  IMPLICIT NONE
  CLASS(MappedGeometry_3D), INTENT(inout) :: myGeom
  REAL(prec), INTENT(in)                  :: x(0:myGeom % nS, 0:myGeom % nP, 0:myGeom % nQ)
  REAL(prec), INTENT(in)                  :: y(0:myGeom % nS, 0:myGeom % nP, 0:myGeom % nQ)
  REAL(prec), INTENT(in)                  :: z(0:myGeom % nS, 0:myGeom % nP, 0:myGeom % nQ)
  
     myGeom % x = x
     myGeom % y = y
     myGeom % z = z
     
 END SUBROUTINE SetPositions_MappedGeometry_3D
!
!
!
 SUBROUTINE GetPositions_MappedGeometry_3D( myGeom, x, y, z )
 ! S/R GetPositions
 !
 ! =============================================================================================== !
 ! DECLARATIONS
  IMPLICIT NONE
  CLASS(MappedGeometry_3D), INTENT(in) :: myGeom
  REAL(prec), INTENT(out)              :: x(0:myGeom % nS, 0:myGeom % nP, 0:myGeom % nQ)
  REAL(prec), INTENT(out)              :: y(0:myGeom % nS, 0:myGeom % nP, 0:myGeom % nQ)
  REAL(prec), INTENT(out)              :: z(0:myGeom % nS, 0:myGeom % nP, 0:myGeom % nQ)
  
     x = myGeom % x
     y = myGeom % y
     z = myGeom % z
     
 END SUBROUTINE GetPositions_MappedGeometry_3D
!
!
!
 SUBROUTINE SetPositionAtNode_MappedGeometry_3D( myGeom, x, y, z, i, j, k )
 ! S/R SetPositionAtNode
 !
 ! =============================================================================================== !
 ! DECLARATIONS
  IMPLICIT NONE
  CLASS(MappedGeometry_3D), INTENT(inout) :: myGeom
  REAL(prec), INTENT(in)                  :: x, y, z
  INTEGER, INTENT(in)                     :: i, j, k
  
     myGeom % x(i,j,k) = x
     myGeom % y(i,j,k) = y
     myGeom % z(i,j,k) = z
     
 END SUBROUTINE SetPositionAtNode_MappedGeometry_3D
!
!
!
 SUBROUTINE GetPositionAtNode_MappedGeometry_3D( myGeom, x, y, z, i, j, k )
 ! S/R GetPositionAtNode
 !
 ! =============================================================================================== !
 ! DECLARATIONS
  IMPLICIT NONE
  CLASS(MappedGeometry_3D), INTENT(in) :: myGeom
  REAL(prec), INTENT(out)              :: x, y, z
  INTEGER, INTENT(in)                  :: i, j, k
  
     x = myGeom % x(i,j,k)
     y = myGeom % y(i,j,k)
     z = myGeom % z(i,j,k)
     
 END SUBROUTINE GetPositionAtNode_MappedGeometry_3D
!
!
!
 SUBROUTINE SetJacobian_MappedGeometry_3D( myGeom, J )
 ! S/R SetJacobian
 !
 ! =============================================================================================== !
 ! DECLARATIONS
  IMPLICIT NONE
  CLASS(MappedGeometry_3D), INTENT(inout) :: myGeom
  REAL(prec), INTENT(in)                  :: J(0:myGeom % nS, 0:myGeom % nP, 0:myGeom % nQ)
  
     myGeom % J = J
     
 END SUBROUTINE SetJacobian_MappedGeometry_3D
!
!
!
 SUBROUTINE GetJacobian_MappedGeometry_3D( myGeom, J )
 ! S/R GetJacobian
 !
 ! =============================================================================================== !
 ! DECLARATIONS
  IMPLICIT NONE
  CLASS(MappedGeometry_3D), INTENT(in) :: myGeom
  REAL(prec), INTENT(out)              :: J(0:myGeom % nS, 0:myGeom % nP, 0:myGeom % nQ)
  
     J = myGeom % J
     
 END SUBROUTINE GetJacobian_MappedGeometry_3D
!
!
!
 SUBROUTINE SetJacobianAtNode_MappedGeometry_3D( myGeom, J, iS, iP, iQ )
 ! S/R SetJacobianAtNode
 !
 ! =============================================================================================== !
 ! DECLARATIONS
  IMPLICIT NONE
  CLASS(MappedGeometry_3D), INTENT(inout) :: myGeom
  REAL(prec), INTENT(in)                  :: J
  INTEGER, INTENT(in)                     :: iS, iP, iQ
  
     myGeom % J(iS,iP,iQ) = J
     
 END SUBROUTINE SetJacobianAtNode_MappedGeometry_3D
!
!
!
 SUBROUTINE GetJacobianAtNode_MappedGeometry_3D( myGeom, J, iS, iP, iQ )
 ! S/R GetJacobianAtNode
 !
 ! =============================================================================================== !
 ! DECLARATIONS
  IMPLICIT NONE
  CLASS(MappedGeometry_3D), INTENT(in) :: myGeom
  REAL(prec), INTENT(out)              :: J
  INTEGER, INTENT(in)                  :: iS, iP, iQ
  
     J = myGeom % J(iS,iP, iQ)
     
 END SUBROUTINE GetJacobianAtNode_MappedGeometry_3D
!
!
!
 SUBROUTINE SetCovariantMetrics_MappedGeometry_3D( myGeom, dxds, dxdp, dxdq, dyds, dydp, dydq, &
                                                           dzds, dzdp, dzdq )
 ! S/R SetCovariantMetrics
 !
 ! =============================================================================================== !
 ! DECLARATIONS
  IMPLICIT NONE
  CLASS(MappedGeometry_3D), INTENT(inout) :: myGeom
  REAL(prec), INTENT(in)                  :: dxds(0:myGeom % nS, 0:myGeom % nP, 0:myGeom % nQ)
  REAL(prec), INTENT(in)                  :: dxdp(0:myGeom % nS, 0:myGeom % nP, 0:myGeom % nQ)
  REAL(prec), INTENT(in)                  :: dxdq(0:myGeom % nS, 0:myGeom % nP, 0:myGeom % nQ)
  REAL(prec), INTENT(in)                  :: dyds(0:myGeom % nS, 0:myGeom % nP, 0:myGeom % nQ)
  REAL(prec), INTENT(in)                  :: dydp(0:myGeom % nS, 0:myGeom % nP, 0:myGeom % nQ)
  REAL(prec), INTENT(in)                  :: dydq(0:myGeom % nS, 0:myGeom % nP, 0:myGeom % nQ)
  REAL(prec), INTENT(in)                  :: dzds(0:myGeom % nS, 0:myGeom % nP, 0:myGeom % nQ)
  REAL(prec), INTENT(in)                  :: dzdp(0:myGeom % nS, 0:myGeom % nP, 0:myGeom % nQ)
  REAL(prec), INTENT(in)                  :: dzdq(0:myGeom % nS, 0:myGeom % nP, 0:myGeom % nQ)
  
     myGeom % dxds = dxds
     myGeom % dxdp = dxdp
     myGeom % dxdq = dxdq
     myGeom % dyds = dyds
     myGeom % dydp = dydp
     myGeom % dydq = dydq
     myGeom % dzds = dzds
     myGeom % dzdp = dzdp
     myGeom % dzdq = dzdq
     
 END SUBROUTINE SetCovariantMetrics_MappedGeometry_3D
!
!
!
 SUBROUTINE GetCovariantMetrics_MappedGeometry_3D( myGeom, dxds, dxdp, dxdq, dyds, dydp, dydq, &
                                                           dzds, dzdp, dzdq )
 ! S/R GetCovariantMetrics
 !
 ! =============================================================================================== !
 ! DECLARATIONS
  IMPLICIT NONE
  CLASS(MappedGeometry_3D), INTENT(in) :: myGeom
  REAL(prec), INTENT(out)              :: dxds(0:myGeom % nS, 0:myGeom % nP, 0:myGeom % nQ)
  REAL(prec), INTENT(out)              :: dxdp(0:myGeom % nS, 0:myGeom % nP, 0:myGeom % nQ)
  REAL(prec), INTENT(out)              :: dxdq(0:myGeom % nS, 0:myGeom % nP, 0:myGeom % nQ)
  REAL(prec), INTENT(out)              :: dyds(0:myGeom % nS, 0:myGeom % nP, 0:myGeom % nQ)
  REAL(prec), INTENT(out)              :: dydp(0:myGeom % nS, 0:myGeom % nP, 0:myGeom % nQ)
  REAL(prec), INTENT(out)              :: dydq(0:myGeom % nS, 0:myGeom % nP, 0:myGeom % nQ)
  REAL(prec), INTENT(out)              :: dzds(0:myGeom % nS, 0:myGeom % nP, 0:myGeom % nQ)
  REAL(prec), INTENT(out)              :: dzdp(0:myGeom % nS, 0:myGeom % nP, 0:myGeom % nQ)
  REAL(prec), INTENT(out)              :: dzdq(0:myGeom % nS, 0:myGeom % nP, 0:myGeom % nQ)
  
     dxds = myGeom % dxds
     dxdp = myGeom % dxdp
     dxdq = myGeom % dxdq
     dyds = myGeom % dyds
     dydp = myGeom % dydp
     dydq = myGeom % dydq
     dzds = myGeom % dzds
     dzdp = myGeom % dzdp
     dzdq = myGeom % dzdq
     
     
 END SUBROUTINE GetCovariantMetrics_MappedGeometry_3D
!
!
!
 SUBROUTINE SetCovariantMetricsAtNode_MappedGeometry_3D( myGeom, dxds, dxdp, dxdq, dyds, dydp, dydq, &
                                                                 dzds, dzdp, dzdq, iS, iP, iQ )
 ! S/R SetCovariantMetricsAtNode
 !
 ! =============================================================================================== !
 ! DECLARATIONS
  IMPLICIT NONE
  CLASS(MappedGeometry_3D), INTENT(inout) :: myGeom
  REAL(prec), INTENT(in)                  :: dxds, dxdp, dxdq, dyds, dydp, dydq, dzds, dzdp, dzdq
  INTEGER, INTENT(in)                     :: iS, iP, iQ
  
     myGeom % dxds(iS,iP,iQ) = dxds
     myGeom % dxdp(iS,iP,iQ) = dxdp
     myGeom % dxdq(iS,iP,iQ) = dxdq
     myGeom % dyds(iS,iP,iQ) = dyds
     myGeom % dydp(iS,iP,iQ) = dydp
     myGeom % dydq(iS,iP,iQ) = dydq
     myGeom % dzds(iS,iP,iQ) = dzds
     myGeom % dzdp(iS,iP,iQ) = dzdp
     myGeom % dzdq(iS,iP,iQ) = dzdq
     
 END SUBROUTINE SetCovariantMetricsAtNode_MappedGeometry_3D
!
!
!
 SUBROUTINE GetCovariantMetricsAtNode_MappedGeometry_3D( myGeom, dxds, dxdp, dxdq, dyds, dydp, dydq, &
                                                                 dzds, dzdp, dzdq, iS, iP, iQ )
 ! S/R GetCovariantMetricsAtNode
 !
 ! =============================================================================================== !
 ! DECLARATIONS
  IMPLICIT NONE
  CLASS(MappedGeometry_3D), INTENT(in) :: myGeom
  REAL(prec), INTENT(out)              :: dxds, dxdp, dxdq, dyds, dydp, dydq, dzds, dzdp, dzdq
  INTEGER, INTENT(in)                  :: iS, iP, iQ
  
     dxds = myGeom % dxds(iS,iP,iQ)
     dxdp = myGeom % dxdp(iS,iP,iQ)
     dxdq = myGeom % dxdq(iS,iP,iQ)
     dyds = myGeom % dyds(iS,iP,iQ)
     dydp = myGeom % dydp(iS,iP,iQ)
     dydq = myGeom % dydq(iS,iP,iQ)
     dzds = myGeom % dzds(iS,iP,iQ)
     dzdp = myGeom % dzdp(iS,iP,iQ)
     dzdq = myGeom % dzdq(iS,iP,iQ)
     
 END SUBROUTINE GetCovariantMetricsAtNode_MappedGeometry_3D
!
!
!
 SUBROUTINE SetBoundaryLocation_MappedGeometry_3D( myGeom, x, y, z, iBound )
 ! S/R SetBoundaryLocation
 !
 ! =============================================================================================== !
 ! DECLARATIONS
  IMPLICIT NONE
  CLASS(MappedGeometry_3D), INTENT(inout) :: myGeom
  INTEGER, INTENT(in)                     :: iBound
  REAL(prec), INTENT(in)                  :: x(0:myGeom % nMax,0:myGeom % nMax)
  REAL(prec), INTENT(in)                  :: y(0:myGeom % nMax,0:myGeom % nMax)
  REAL(prec), INTENT(in)                  :: z(0:myGeom % nMax,0:myGeom % nMax)
  
     myGeom % xBound(:,:,iBound) = x
     myGeom % yBound(:,:,iBound) = y
     myGeom % zBound(:,:,iBound) = z
     
 END SUBROUTINE SetBoundaryLocation_MappedGeometry_3D
!
!
!
 SUBROUTINE GetBoundaryLocation_MappedGeometry_3D( myGeom, x, y, z, iBound )
 ! S/R GetBoundaryLocation
 !
 ! =============================================================================================== !
 ! DECLARATIONS
  IMPLICIT NONE
  CLASS(MappedGeometry_3D), INTENT(in) :: myGeom
  INTEGER, INTENT(in)                  :: iBound
  REAL(prec), INTENT(out)              :: x(0:myGeom % nMax,0:myGeom % nMax)
  REAL(prec), INTENT(out)              :: y(0:myGeom % nMax,0:myGeom % nMax)
  REAL(prec), INTENT(out)              :: z(0:myGeom % nMax,0:myGeom % nMax)
  
     x = myGeom % xBound(:,:,iBound)
     y = myGeom % yBound(:,:,iBound)
     z = myGeom % zBound(:,:,iBound)
     
 END SUBROUTINE GetBoundaryLocation_MappedGeometry_3D
!
!
!
 SUBROUTINE SetBoundaryLocationAtNode_MappedGeometry_3D( myGeom, x, y, z, i, j, iBound )
 ! S/R SetBoundaryLocationAtNode
 !
 ! =============================================================================================== !
 ! DECLARATIONS
  IMPLICIT NONE
  CLASS(MappedGeometry_3D), INTENT(inout) :: myGeom
  REAL(prec), INTENT(in)                  :: x, y, z
  INTEGER, INTENT(in)                     :: i, j, iBound
  
     myGeom % xBound(i,j,iBound) = x
     myGeom % yBound(i,j,iBound) = y
     myGeom % zBound(i,j,iBound) = z
     
 END SUBROUTINE SetBoundaryLocationAtNode_MappedGeometry_3D
!
!
!
 SUBROUTINE GetBoundaryLocationAtNode_MappedGeometry_3D( myGeom, x, y, z, i, j, iBound )
 ! S/R GetBoundaryLocationAtNode
 !
 ! =============================================================================================== !
 ! DECLARATIONS
  IMPLICIT NONE
  CLASS(MappedGeometry_3D), INTENT(in) :: myGeom
  REAL(prec), INTENT(out)              :: x, y, z
  INTEGER, INTENT(in)                  :: i, j, iBound
  
     x = myGeom % xBound(i,j,iBound)
     y = myGeom % yBound(i,j,iBound)
     z = myGeom % zBound(i,j,iBound)
     
 END SUBROUTINE GetBoundaryLocationAtNode_MappedGeometry_3D
!
!
!
 SUBROUTINE SetBoundaryNormalAtNode_MappedGeometry_3D( myGeom, dir, i, j, iBound )
 ! S/R SetBoundaryNormalAtNode
 !
 ! =============================================================================================== !
 ! DECLARATIONS
  IMPLICIT NONE
  CLASS(MappedGeometry_3D), INTENT(inout) :: myGeom
  REAL(prec), INTENT(in)                  :: dir(1:nDims)
  INTEGER, INTENT(in)                     :: i, j, iBound
  
     CALL myGeom % nHat(i,j,iBound) % SetDirection( dir )

     CALL myGeom % nHat(i,j,iBound) % SetLength( SQRT(DOT_PRODUCT(dir,dir)) )
     
 END SUBROUTINE SetBoundaryNormalAtNode_MappedGeometry_3D
!
!
!
 SUBROUTINE GetBoundaryNormalAtNode_MappedGeometry_3D( myGeom, dir, length, i, j, iBound )
 ! S/R GetBoundaryNormalAtNode
 !
 ! =============================================================================================== !
 ! DECLARATIONS
  IMPLICIT NONE
  CLASS(MappedGeometry_3D), INTENT(in) :: myGeom
  REAL(prec), INTENT(out)              :: dir(1:nDims), length
  INTEGER, INTENT(in)                  :: i, j, iBound
  
     CALL myGeom % nHat(i,j,iBound) % GetDirection( dir )

     CALL myGeom % nHat(i,j,iBound) % GetLength( length )
     
 END SUBROUTINE GetBoundaryNormalAtNode_MappedGeometry_3D
!
!
!==================================================================================================!
!--------------------------------- Type Specific Routines -----------------------------------------!
!==================================================================================================!
!
!
 SUBROUTINE GenerateMesh_MappedGeometry_3D( myGeom, interp, theSurfaces )
 ! S/R GenerateMesh
 !
 ! =============================================================================================== !
 ! DECLARATIONS
  IMPLICIT NONE
  CLASS( MappedGeometry_3D ), INTENT(inout) :: myGeom
  TYPE( Lagrange_3D ), INTENT(in)           :: interp
  TYPE( Surface_3D ), INTENT(in)            :: theSurfaces(1:nHexFaces)
  ! Local
  INTEGER    :: iS, iP, iQ, nS, nP, nQ
  REAL(prec) :: s, p, q, x(1:ndims)
   
      CALL interp % GetNumberOfNodes( nS, nP, nQ )
      
      DO iQ = 0, nQ
         DO iP = 0,nP
            DO iS = 0, nS

               CALL interp % GetNode( iS, iP, iQ, s, p, q ) 
               ! Mesh location is computed via transfinite interpolation
               x = TransfiniteInterpolation( theSurfaces, s, p, q )
            
               CALL myGeom % SetPositionAtNode( x(1), x(2), x(3), iS, iP, iQ )
               
            ENDDO ! iS
         ENDDO ! iP
      ENDDO ! iQ
      
      CALL myGeom % GenerateBoundaryPositions( interp, theSurfaces )
      
 END SUBROUTINE GenerateMesh_MappedGeometry_3D
!
!
!
SUBROUTINE GenerateBoundaryPositions_MappedGeometry_3D( myGeom, interp, theSurfaces  )
 ! S/R GenerateBoundaryPositions
 !
 ! =============================================================================================== !
 ! DECLARATIONS
  IMPLICIT NONE
  CLASS( MappedGeometry_3D ), INTENT(inout) :: myGeom
  TYPE( Lagrange_3D ), INTENT(in)           :: interp
  TYPE( Surface_3D ), INTENT(in)            :: theSurfaces(1:nHexFaces)
   ! Local
  INTEGER    :: iS, iP, iQ, nS, nP, nQ
  REAL(prec) :: s, p, q, x(1:nDims)
   
      CALL interp % GetNumberOfNodes( nS, nP, nQ )
  
      
      ! Now the boundary metrics are calculated.
      DO iQ = 0, nQ
         DO iS = 0, nS
            CALL interp % GetNode( iS, 0, iQ, s, p, q )
            
            p = -ONE  ! south boundary
            x = TransfiniteInterpolation( theSurfaces, s, p, q )
            CALL myGeom % SetBoundaryLocationAtNode( x(1), x(2), x(3), iS, iQ, south )
            
            p = ONE  ! south boundary
            x = TransfiniteInterpolation( theSurfaces, s, p, q )
            CALL myGeom % SetBoundaryLocationAtNode( x(1), x(2), x(3), iS, iQ, north )
            
         ENDDO
      ENDDO
      
       DO iQ = 0, nQ
         DO iP = 0, nP
            CALL interp % GetNode( 0, iP, iQ, s, p, q )
            
            s = -ONE  ! west boundary
            x = TransfiniteInterpolation( theSurfaces, s, p, q )
            CALL myGeom % SetBoundaryLocationAtNode( x(1), x(2), x(3), iP, iQ, west )
            
            s = ONE  ! east boundary
            x = TransfiniteInterpolation( theSurfaces, s, p, q )
            CALL myGeom % SetBoundaryLocationAtNode( x(1), x(2), x(3), iP, iQ, east )
            
         ENDDO
      ENDDO
      
      DO iP = 0, nP
         DO iS = 0, nS
            CALL interp % GetNode( iS, iP, 0, s, p, q )
            
            q = -ONE  ! bottom boundary
            x = TransfiniteInterpolation( theSurfaces, s, p, q )
            CALL myGeom % SetBoundaryLocationAtNode( x(1), x(2), x(3), iS, iP, bottom )
            
            q = ONE  ! top boundary
            x = TransfiniteInterpolation( theSurfaces, s, p, q )
            CALL myGeom % SetBoundaryLocationAtNode( x(1), x(2), x(3), iS, iP, top )
            
         ENDDO
      ENDDO
      
 END SUBROUTINE GenerateBoundaryPositions_MappedGeometry_3D
!
!
!
 SUBROUTINE GenerateMetrics_MappedGeometry_3D( myGeom, interp )
 ! S/R GenerateMetrics
 !
 ! =============================================================================================== !
 ! DECLARATIONS
  IMPLICIT NONE
  CLASS( MappedGeometry_3D ), INTENT(inout) :: myGeom
  TYPE( Lagrange_3D ), INTENT(in)           :: interp
  ! Local
  INTEGER    :: iS, iP, iQ, nS, nP, nQ
  REAL(prec) :: s, p, q, covT(1:3,1:3), J
   
      CALL interp % GetNumberOfNodes( nS, nP, nQ )
      
      DO iQ = 0, nQ
         DO iP = 0, nP
            DO iS = 0, nS
      
               CALL interp % GetNode( iS, iP, iQ, s, p, q ) 
            
               ! The metric terms can be computed using a transfinite interpolation formula or can
               ! be computed using the derivative of the interpolant. Of course, the latter is only
               ! possible if the element positions are already given. It is assumed that this is the
               ! case, so that we can save on writing another subroutine. This way we can call the
               ! MappedGeometryClass_3D routine "CalculateMetrics" to obtain the covariant metric
               ! tensor.
               CALL myGeom % CalculateMetrics( interp, s, p, q, covT )
               CALL myGeom % SetCovariantMetricsAtNode( covT(1,1), & ! dxds
                                                        covT(1,2), & ! dxdp
                                                        covT(1,3), & ! dxdq
                                                        covT(2,1), & ! dyds
                                                        covT(2,2), & ! dydp
                                                        covT(2,3), & ! dydq
                                                        covT(3,1), & ! dzds
                                                        covT(3,2), & ! dzdp
                                                        covT(3,3), & ! dzdq 
                                                        iS, iP, iQ )
            
               J = Determinant( covT, 3 )
               CALL myGeom % SetJacobianAtNode( J, iS, iP, iQ )

            ENDDO ! iQ
         ENDDO ! iP
      ENDDO ! iS
      
      CALL myGeom % GenerateBoundaryMetrics( interp )
    
 END SUBROUTINE GenerateMetrics_MappedGeometry_3D
!
!
! 
 SUBROUTINE GenerateBoundaryMetrics_MappedGeometry_3D( myGeom, interp  )
 ! S/R GenerateBoundaryMetrics
 !
 ! =============================================================================================== !
 ! DECLARATIONS
  IMPLICIT NONE
  CLASS( MappedGeometry_3D ), INTENT(inout) :: myGeom
  TYPE( Lagrange_3D ), INTENT(in)           :: interp
   ! Local
  INTEGER    :: iS, iP, iQ, nS, nP, nQ
  REAL(prec) :: s, p, q, covT(1:3,1:3), J, signJ, nhat(1:nDims)
   
      CALL interp % GetNumberOfNodes( nS, nP, nQ )
  
      
      ! Now the boundary metrics are calculated.
      DO iQ = 0, nQ
         DO iS = 0, nS
            CALL interp % GetNode( iS, 0, iQ, s, p, q )
            p = -ONE  ! south boundary
         
            CALL myGeom % CalculateMetrics( interp, s, p, q, covT )
         
            J = Determinant( covT, 3 )

            signJ = abs(J)/J       
            !                      dydq*dzds      -      dzdq*dyds            
            nHat(1) = -signJ*(covT(2,3)*covT(3,1) - covT(3,3)*covT(2,1))
            !                      dzdq*dxds      -      dxdq*dzds
            nHat(2) = -signJ*(covT(3,3)*covT(1,1) - covT(1,3)*covT(3,1))
            !                      dxdq*dyds      -      dydq*dxds
            nHat(3) = -signJ*(covT(1,3)*covT(2,1) - covT(2,3)*covT(1,1))
            
            CALL myGeom % SetBoundaryNormalAtNode( nHat, iS, iQ, south )
            CALL myGeom % nHat(iS, iQ, south) % Normalize( )
          
            p = ONE ! north boundary
            CALL myGeom % CalculateMetrics( interp, s, p, q, covT )
         
            J = Determinant( covT, 3 )

            signJ = abs(J)/J       
            !                     dydq*dzds      -      dzdq*dyds            
            nHat(1) = signJ*(covT(2,3)*covT(3,1) - covT(3,3)*covT(2,1))
            !                     dzdq*dxds      -      dxdq*dzds
            nHat(2) = signJ*(covT(3,3)*covT(1,1) - covT(1,3)*covT(3,1))
            !                     dxdq*dyds      -      dydq*dxds
            nHat(3) = signJ*(covT(1,3)*covT(2,1) - covT(2,3)*covT(1,1))
            
            CALL myGeom % SetBoundaryNormalAtNode( nHat, iS, iQ, north )
            CALL myGeom % nHat(iS, iQ, north) % Normalize( )
         ENDDO
      ENDDO
      
       DO iQ = 0, nQ
         DO iP = 0, nP
            CALL interp % GetNode( 0, iP, iQ, s, p, q )
            s = -ONE  ! west boundary
         
            CALL myGeom % CalculateMetrics( interp, s, p, q, covT )
         
            J = Determinant( covT, 3 )

            signJ = abs(J)/J       
            !                      dydp*dzdq      -      dzdp*dydq            
            nHat(1) = -signJ*(covT(2,2)*covT(3,3) - covT(3,2)*covT(2,3))
            !                      dzdp*dxdq      -      dxdp*dzdq
            nHat(2) = -signJ*(covT(3,2)*covT(1,3) - covT(1,2)*covT(3,3))
            !                      dxdp*dydq      -      dydp*dxdq
            nHat(3) = -signJ*(covT(1,2)*covT(2,3) - covT(2,2)*covT(1,3))
            
            CALL myGeom % SetBoundaryNormalAtNode( nHat, iP, iQ, west )
            CALL myGeom % nHat(iP, iQ, west) % Normalize( )
          
            s = ONE ! east boundary
            CALL myGeom % CalculateMetrics( interp, s, p, q, covT )
         
            J = Determinant( covT, 3 )

            signJ = abs(J)/J       
            !                     dydp*dzdq      -      dzdp*dydq            
            nHat(1) = signJ*(covT(2,2)*covT(3,3) - covT(3,2)*covT(2,3))
            !                     dzdp*dxdq      -      dxdp*dzdq
            nHat(2) = signJ*(covT(3,2)*covT(1,3) - covT(1,2)*covT(3,3))
            !                     dxdp*dydq      -      dydp*dxdq
            nHat(3) = signJ*(covT(1,2)*covT(2,3) - covT(2,2)*covT(1,3))
            
            CALL myGeom % SetBoundaryNormalAtNode( nHat, iP, iQ, east )
            CALL myGeom % nHat(iP, iQ, east) % Normalize( )
         ENDDO
      ENDDO
      
      DO iP = 0, nP
         DO iS = 0, nS
            CALL interp % GetNode( iS, iP, 0, s, p, q )
            q = -ONE  ! bottom boundary
         
            CALL myGeom % CalculateMetrics( interp, s, p, q, covT )
         
            J = Determinant( covT, 3 )

            signJ = abs(J)/J       
            !                      dyds*dzdp      -      dzds*dydp            
            nHat(1) = -signJ*(covT(2,1)*covT(3,2) - covT(3,1)*covT(2,2))
            !                      dzds*dxdp      -      dxds*dzdp
            nHat(2) = -signJ*(covT(3,1)*covT(1,2) - covT(1,1)*covT(3,2))
            !                      dxds*dydp      -      dyds*dxdp
            nHat(3) = -signJ*(covT(1,1)*covT(2,2) - covT(2,1)*covT(1,2))
            
            CALL myGeom % SetBoundaryNormalAtNode( nHat, iS, iP, bottom )
            CALL myGeom % nHat(iS, iP, bottom) % Normalize( )
          
            q = ONE ! top boundary
            CALL myGeom % CalculateMetrics( interp, s, p, q, covT )
         
            J = Determinant( covT, 3 )

            signJ = abs(J)/J       
            !                      dyds*dzdp      -      dzds*dydp            
            nHat(1) = signJ*(covT(2,1)*covT(3,2) - covT(3,1)*covT(2,2))
            !                      dzds*dxdp      -      dxds*dzdp
            nHat(2) = signJ*(covT(3,1)*covT(1,2) - covT(1,1)*covT(3,2))
            !                      dxds*dydp      -      dyds*dxdp
            nHat(3) = signJ*(covT(1,1)*covT(2,2) - covT(2,1)*covT(1,2))
            
            CALL myGeom % SetBoundaryNormalAtNode( nHat, iS, iP, top )
            CALL myGeom % nHat(iS, iP, top) % Normalize( )
         ENDDO
      ENDDO
      
 END SUBROUTINE GenerateBoundaryMetrics_MappedGeometry_3D
!
!
!
 FUNCTION TransfiniteInterpolation( surfaces, a, b, c ) RESULT( P )
 ! TransfiniteInterpolation
 !  Takes in the six surfaces (south, east, north, west, bottom, top) and evaluates the 
 !  bidirectional mapping at xi^1 = a, xi^2 = b, xi^3 = c. The boundary of the computational 
 !  coordinate system is assumed to be at +/- 1 in each direction.
 !
 !   A volume in 3-D is the assumed geometrical object generated from this evaluation.
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   TYPE( Surface_3D ) :: surfaces(1:nHexFaces)
   REAL(prec)       :: a, b, c
   REAL(prec)       :: P(1:nDims)
   ! LOCAL
   REAL(prec)  :: P1(1:nDims), P2(1:nDims), P3(1:nDims)
   REAL(prec)  :: sSurf(1:nDims), nSurf(1:nDims), eSurf(1:nDims), wSurf(1:nDims), bSurf(1:nDims), tSurf(1:nDims)
   REAL(prec)  :: l1(1:2), l2(1:2), l3(1:2)
   REAL(prec)  :: ref(1:2)
   INTEGER     :: i, j
   
      ref = (/ -ONE, ONE /)

      ! Transfinite interpolation with linear blending uses linear lagrange interpolating polynomials
      ! to blend the bounding surfaces.
      ! The linear blending weights in the first computational direction are calculated.
      l1 = LinearBlend( a )
      
      ! The bounding surfaces need to be evaluated at the provided computational coordinates
      CALL surfaces(west) % Evaluate( b, c, wSurf(1), wSurf(2), wSurf(3) )
      CALL surfaces(east) % Evaluate( b, c, eSurf(1), eSurf(2), eSurf(3) )
      
      ! P1 contains the interpolation in the first computational coordinate
      ! The first computational coordinate is assumed to vary between the (computational) east and
      ! west boundaries.
      P1 = wSurf*l1(1) + eSurf*l1(2)     

      ! The linear blending weights in the second computational direction are calculated.
      l2 = LinearBlend( b )
      
      ! The bounding surfaces need to be evaluated at the provided computational coordinates
      CALL surfaces(south) % Evaluate( a, c, sSurf(1), sSurf(2), sSurf(3) )
      CALL surfaces(north) % Evaluate( a, c, nSurf(1), nSurf(2), nSurf(3) )
      
      ! P2 contains the interpolation in the second computational coordinate
      ! The second computational coordinate is assumed to vary between the (computational) south and
      ! north boundaries.
      P2 = sSurf*l2(1) + nSurf*l2(2)  

      ! The linear blending weights in the third computational direction are calculated.
      l3 = LinearBlend( c )
      
      ! The bounding surfaces need to be evaluated at the provided computational coordinates
      CALL surfaces(bottom) % Evaluate( a, b, bSurf(1), bSurf(2), bSurf(3) )
      CALL surfaces(top) % Evaluate( a, b, tSurf(1), tSurf(2), tSurf(3) )
      
      ! P3 contains the interpolation in the first computational coordinate
      ! The first computational coordinate is assumed to vary between the (computational) bottom and
      ! top boundaries.
      P3 = bSurf*l3(1) + tSurf*l3(2)  
     
      DO i = 1, 2
         ! Now we need to compute the tensor product of the first and second computational direction 
         ! interpolants and subtract from P1.
         CALL surfaces(west) % Evaluate( ref(i), c, wSurf(1), wSurf(2), wSurf(3) )
         CALL surfaces(east) % Evaluate( ref(i), c, eSurf(1), eSurf(2), eSurf(3) )
         P1 = P1 - l2(i)*( wSurf*l1(1) + eSurf*l1(2) )
         
         ! Now we need to compute the tensor product of the first and third computational direction 
         ! interpolants and subtract from P1.
         CALL surfaces(west) % Evaluate( b, ref(i), wSurf(1), wSurf(2), wSurf(3) )
         CALL surfaces(east) % Evaluate( b, ref(i), eSurf(1), eSurf(2), eSurf(3) )

         P1 = P1 - l3(i)*( wSurf*l1(1) + eSurf*l1(2) )
      
         ! Now we need to compute the tensor product of the second and third computational direction 
         ! interpolants and subtract from P2.
         CALL surfaces(south) % Evaluate( a, ref(i), sSurf(1), sSurf(2), sSurf(3) )
         CALL surfaces(north) % Evaluate( a, ref(i), nSurf(1), nSurf(2), nSurf(3) )

         P2 = P2 - l3(i)*( sSurf*l2(1) + nSurf*l2(2) )
      
      ENDDO

      ! Next, the compounded tensor product is computed and added to P3.
      DO j = 1,2
         DO i = 1,2
         
            CALL surfaces(west) % Evaluate( ref(i), ref(j), wSurf(1), wSurf(2), wSurf(3) )
            CALL surfaces(east) % Evaluate( ref(i), ref(j), eSurf(1), eSurf(2), eSurf(3) )
            P3 = P3 + l2(i)*l3(j)*( wSurf*l1(1) + eSurf*l1(2) )
         
         ENDDO
      ENDDO
      
      !Finally, the sum the interpolants is computed to yield the computational coordinate
      P = P1 + P2 + P3
      
 END FUNCTION TransfiniteInterpolation
!
!
!
 FUNCTION LinearBlend( a ) RESULT( weights )
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   REAL(prec) :: a
   REAL(prec) :: weights(1:2)

       weights(1) = HALF*(ONE - a)
       weights(2) = HALF*(ONE + a)
    
 END FUNCTION LinearBlend
!
!
!
 SUBROUTINE CalculateLocation_MappedGeometry_3D( myGeom, interp, s, p, q, x, y, z )
 ! S/R CalculateLocation
 !  Description :
 ! 
 !     This subroutine calculates the physical coordinates (x,y) given the computational coordinates
 !     (s,p).
 ! =============================================================================================== !
  IMPLICIT NONE
  CLASS( MappedGeometry_3D ), INTENT(in) :: myGeom
  TYPE( Lagrange_3D ), INTENT(in)        :: interp
  REAL(prec), INTENT(in)                 :: s, p, q
  REAL(prec), INTENT(out)                :: x, y, z
  
     x = interp % EvaluateInterpolant( s, p, q, myGeom % x )
     y = interp % EvaluateInterpolant( s, p, q, myGeom % y )
     z = interp % EvaluateInterpolant( s, p, q, myGeom % z )
  
 END SUBROUTINE CalculateLocation_MappedGeometry_3D
!
!
!
 SUBROUTINE CalculateMetrics_MappedGeometry_3D( myGeom, interp, s, p, q, covT )
 ! S/R CalculateMetrics
 !  Description :
 ! 
 !     This subroutine calculates the covariant metric terms given the computational coordinates
 !     (s,p).
 !     covT(1,1) -> dxds
 !     covT(1,2) -> dxdp
 !     covT(1,3) -> dxdq
 !     covT(2,1) -> dyds
 !     covT(2,2) -> dydp
 !     covT(2,3) -> dydq
 !     covT(3,1) -> dzds
 !     covT(3,2) -> dzdp
 !     covT(3,3) -> dzdq
 ! =============================================================================================== !
  IMPLICIT NONE
  CLASS( MappedGeometry_3D ), INTENT(in) :: myGeom
  TYPE( Lagrange_3D ), INTENT(in)        :: interp
  REAL(prec), INTENT(in)                 :: s, p, q
  REAL(prec), INTENT(out)                :: covT(1:3,1:3)

     covT(1,1:3) = interp % EvaluateDerivative( s, p, q, myGeom % x )
     covT(2,1:3) = interp % EvaluateDerivative( s, p, q, myGeom % y )
     covT(3,1:3) = interp % EvaluateDerivative( s, p, q, myGeom % z )
     
 END SUBROUTINE CalculateMetrics_MappedGeometry_3D
!
!
!
 SUBROUTINE CalculateComputationalCoordinates_MappedGeometry_3D( myGeom, interp, x, y, z, s, p, q, success )
 ! S/R CalculateComputationalCoordinates
 !  Description :
 ! 
 !     Given the physical coordinates (x*,y*,z*), the computational coordinates are calculated using 
 !     Newton's method for root finding to solve
 !
 !     x* = x(s,p,q),
 !     y* = y(s,p,q),
 !     z* = z(s,p,q)
 !
 !     for s, p, and q.
 ! =============================================================================================== !
  IMPLICIT NONE
  CLASS( MappedGeometry_3D ), INTENT(in) :: myGeom
  TYPE( Lagrange_3D ), INTENT(in)        :: interp
  REAL(prec), INTENT(in)                 :: x, y, z
  REAL(prec), INTENT(out)                :: s, p, q
  LOGICAL, INTENT(out), OPTIONAL         :: success
  ! LOCAL
  REAL(prec) :: dr(1:nDims), ds(1:nDims), A(1:nDims,1:nDims), Ainv(1:nDims,1:nDims)
  REAL(prec) :: thisX, thisY, thisZ, thisS, thisP, thisQ, resi, r0
  INTEGER    :: i 

     thisS = ZERO ! Initial guess is at the origin
     thisP = ZERO
     thisQ = ZERO
     
     IF( PRESENT(success) )THEN
        success = .FALSE.
     ENDIF
     
     DO i = 1, newtonMax
     
        ! Calculate the physical coordinate associated with the computational coordinate guess
        CALL myGeom % CalculateLocation( interp, thisS, thisP, thisQ, thisX, thisY, thisZ )
     
        ! Calculate the residual
        dr(1) = x - thisX
        dr(2) = y - thisY  
        dr(3) = z - thisZ
        resi = SQRT( DOT_PRODUCT( dr, dr ) )
     
        IF( resi < newtonTolerance )THEN
           s = thisS
           p = thisP
           q = thisQ
           IF( PRESENT(success) )THEN
              success = .TRUE.
           ENDIF
           RETURN
        ENDIF
        
        CALL myGeom % CalculateMetrics( interp, thisS, thisP, thisQ, A ) ! Calculate the covariant metric tensor
     
        Ainv = Invert_3x3( A ) ! Invert the covariant metric tensor.
                           ! This matrix is ivertable as long as the Jacobian is non-zero.
        
        ds = MATMUL( Ainv, dr ) ! calculate the correction in the computational coordinate
        thisS = thisS + ds(1)
        thisP = thisP + ds(2)
        thisQ = thisQ + ds(3)
     
     ENDDO
     
     ! Calculate the residual
     dr(1) = x - thisX
     dr(2) = y - thisY  
     dr(3) = z - thisZ
     resi = SQRT( DOT_PRODUCT( dr, dr ) )
     IF( resi < newtonTolerance )THEN
        s = thisS
        p = thisP
        q = thisQ
        IF( PRESENT(success) )THEN
           success = .TRUE.
        ENDIF
        RETURN
     ELSE
        s = fillValue
        p = fillValue
        q = fillValue
        PRINT*,'Module MappedGeometryClass_3D.f90 : S/R CalculateComputationalCoordinates :'
        PRINT*,'Search for coordinates failed. Final residual norm :', resi
        RETURN
    ENDIF
      
     
 END SUBROUTINE CalculateComputationalCoordinates_MappedGeometry_3D
!
!
!
 SUBROUTINE ScaleGeometry_MappedGeometry_3D( myGeom, interp, xScale, yScale, zScale )
 !
 !
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS( MappedGeometry_3D ), INTENT(inout) :: myGeom
   TYPE( Lagrange_3D ), INTENT(in)           :: interp
   REAL(prec), INTENT(in)                    :: xScale, yScale, zScale
   !
   
         myGeom % x = xScale*( myGeom % x )
         myGeom % y = yScale*( myGeom % y )
         myGeom % z = zScale*( myGeom % z )
         myGeom % xBound = xScale*( myGeom % xBound )
         myGeom % yBound = yScale*( myGeom % yBound )
         myGeom % zBound = zScale*( myGeom % zBound )

         myGeom % dxds = xScale*( myGeom % dxds )
         myGeom % dxdp = xScale*( myGeom % dxdp )
         myGeom % dxdq = xScale*( myGeom % dxdq )
         myGeom % dyds = yScale*( myGeom % dyds )
         myGeom % dydp = yScale*( myGeom % dydp )
         myGeom % dydq = yScale*( myGeom % dydq )
         myGeom % dzds = yScale*( myGeom % dzds )
         myGeom % dzdp = yScale*( myGeom % dzdp )
         myGeom % dzdq = yScale*( myGeom % dzdq )
          
         myGeom % J = xScale*yScale*zScale*( myGeom % J )

         ! Update the boundary metrics -- normals and normal lengths
         CALL myGeom % GenerateBoundaryMetrics( interp  )
         
 END SUBROUTINE ScaleGeometry_MappedGeometry_3D
!
!
!==================================================================================================!
!------------------------------------ File I/O Routines -------------------------------------------!
!==================================================================================================!
!
!
 SUBROUTINE WriteTecplot_MappedGeometry_3D( myGeom, filename )
 ! S/R WriteTecplot
 !  Description :
 ! 
 ! =============================================================================================== !
  IMPLICIT NONE
  CLASS( MappedGeometry_3D ), INTENT(in) :: myGeom
  CHARACTER(*), INTENT(in), OPTIONAL     :: filename  
  ! Local
  INTEGER :: iX, iY, iZ, nX, nY, nZ, fUnit
  REAL(prec) :: x, y, z, dxds, dxdp, dxdq, dyds, dydp, dydq, dzds, dzdp, dzdq, J
  
    CALL myGeom % GetNumberOfNodes( nX, nY, nZ )
    
    IF( PRESENT(filename) )THEN
       OPEN( UNIT   = NEWUNIT(fUnit), &
             FILE   = TRIM(filename)//'.tec', &
             FORM   = 'formatted', & 
             STATUS = 'REPLACE' )
    ELSE
       OPEN( UNIT   = NEWUNIT(fUnit), &
             FILE   = 'LocalGeometry.tec', &
             FORM   = 'formatted', & 
             STATUS = 'REPLACE' )
    ENDIF

    WRITE(fUnit,*) 'VARIABLES = "X", "Y", "Z", "Jacobian", "dxds", "dxdp", "dxdq", "dyds", "dydp" &
                               &"dydq", "dzds", "dzdp", "dzdq" '
    
    WRITE(fUnit,*)  'ZONE T="el0", I=',nX+1,', J=', nY+1,', K=',nZ+1,',F=POINT'
    DO iZ = 0, nZ
       DO iY = 0, nY
          DO iX = 0, nX
          
             CALL myGeom % GetPositionAtNode( x, y, z, iX, iY, iZ )
             CALL myGeom % GetJacobianAtNode( J, iX, iY, iZ )
             CALL myGeom % GetCovariantMetricsAtNode( dxds, dxdp, dxdq, &
                                                      dyds, dydp, dydq, &
                                                      dzds, dzdp, dzdq, &
                                                      iX, iY, iZ )

             WRITE(fUnit,*)  x, y, z, J, dxds, dxdp, dxdq, dyds, dydp, dydq, dzds, dzdp, dzdq
             
          ENDDO
       ENDDO
    ENDDO

    CLOSE(UNIT=fUnit)

    RETURN

 END SUBROUTINE WriteTecplot_MappedGeometry_3D
 
 
END MODULE MappedGeometryClass_3D
