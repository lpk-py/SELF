! MappedGeometryClass_2D.f90
! 
! Copyright 2015 Joseph Schoonover <schoonover.numerics@gmail.com>
! 
! MappedGeometryClass_2D.f90 is part of the Spectral Element Libraries in Fortran (SELF).
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
 
 
MODULE MappedGeometryClass_2D
! ========================================= Logs ================================================= !
!2016-05-11  Joseph Schoonover  schoonover.numerics@gmail.com 
!
! //////////////////////////////////////////////////////////////////////////////////////////////// ! 
! src/common/
USE ModelPrecision
USE ConstantsDictionary
USE ModelFlags
USE CommonRoutines
! src/interp/
USE Lagrange_2D_Class
! src/geom/
USE CurveClass_2D
USE VectorClass
USE GeometrySupportRoutines



IMPLICIT NONE



   TYPE MappedGeometry_2D
      INTEGER                             :: nS, nP, nMax
      TYPE(Vector), ALLOCATABLE, PRIVATE  :: nHat(:,:) 
      REAL(prec), ALLOCATABLE, PRIVATE    :: xBound(:,:)
      REAL(prec), ALLOCATABLE, PRIVATE    :: yBound(:,:) 
      REAL(prec), ALLOCATABLE, PRIVATE    :: x(:,:), y(:,:)
      REAL(prec), ALLOCATABLE, PRIVATE    :: J(:,:)    
      REAL(prec), ALLOCATABLE, PRIVATE    :: dxds(:,:), dxdp(:,:)
      REAL(prec), ALLOCATABLE, PRIVATE    :: dyds(:,:), dydp(:,:)

      CONTAINS

      PROCEDURE :: Build => Build_MappedGeometry_2D
      PROCEDURE :: Trash => Trash_MappedGeometry_2D
      
      PROCEDURE :: SetNumberOfNodes => SetNumberOfNodes_MappedGeometry_2D
      PROCEDURE :: GetNumberOfNodes => GetNumberOfNodes_MappedGeometry_2D
      PROCEDURE :: SetPositions => SetPositions_MappedGeometry_2D
      PROCEDURE :: GetPositions => GetPositions_MappedGeometry_2D
      PROCEDURE :: SetPositionAtNode => SetPositionAtNode_MappedGeometry_2D
      PROCEDURE :: GetPositionAtNode => GetPositionAtNode_MappedGeometry_2D
      PROCEDURE :: SetJacobian => SetJacobian_MappedGeometry_2D
      PROCEDURE :: GetJacobian => GetJacobian_MappedGeometry_2D
      PROCEDURE :: SetJacobianAtNode => SetJacobianAtNode_MappedGeometry_2D
      PROCEDURE :: GetJacobianAtNode => GetJacobianAtNode_MappedGeometry_2D
      PROCEDURE :: SetCovariantMetrics => SetCovariantMetrics_MappedGeometry_2D
      PROCEDURE :: GetCovariantMetrics => GetCovariantMetrics_MappedGeometry_2D
      PROCEDURE :: SetCovariantMetricsAtNode => SetCovariantMetricsAtNode_MappedGeometry_2D
      PROCEDURE :: GetCovariantMetricsAtNode => GetCovariantMetricsAtNode_MappedGeometry_2D
      PROCEDURE :: SetBoundaryLocation => SetBoundaryLocation_MappedGeometry_2D
      PROCEDURE :: GetBoundaryLocation => GetBoundaryLocation_MappedGeometry_2D
      PROCEDURE :: SetBoundaryLocationAtNode => SetBoundaryLocationAtNode_MappedGeometry_2D
      PROCEDURE :: GetBoundaryLocationAtNode => GetBoundaryLocationAtNode_MappedGeometry_2D
      PROCEDURE :: SetBoundaryNormalAtNode => SetBoundaryNormalAtNode_MappedGeometry_2D
      PROCEDURE :: GetBoundaryNormalAtNode => GetBoundaryNormalAtNode_MappedGeometry_2D

      PROCEDURE :: GenerateMesh => GenerateMesh_MappedGeometry_2D
      PROCEDURE :: GenerateMetrics => GenerateMetrics_MappedGeometry_2D
      PROCEDURE :: GenerateBoundaryMetrics => GenerateBoundaryMetrics_MappedGeometry_2D
      PROCEDURE :: ScaleGeometry => ScaleGeometry_MappedGeometry_2D
      PROCEDURE :: CalculateLocation => CalculateLocation_MappedGeometry_2D
      PROCEDURE :: CalculateMetrics => CalculateMetrics_MappedGeometry_2D
      PROCEDURE :: CalculateComputationalCoordinates => CalculateComputationalCoordinates_MappedGeometry_2D
      
      PROCEDURE :: WriteTecplot => WriteTecplot_MappedGeometry_2D
   END TYPE MappedGeometry_2D


 INTEGER, PRIVATE :: nDims = 2
 CONTAINS
!
!
!==================================================================================================!
!------------------------------- Manual Constructors/Destructors ----------------------------------!
!==================================================================================================!
!
!
 SUBROUTINE Build_MappedGeometry_2D( myGeom, interp, myCurves )
 ! S/R Build
 !
 ! =============================================================================================== !
 ! DECLARATIONS
  IMPLICIT NONE
  CLASS(MappedGeometry_2D), INTENT(out) :: myGeom
  TYPE(Lagrange_2D), INTENT(in)         :: interp
  TYPE(Curve_2D), INTENT(in)            :: myCurves(1:nQuadEdges)
  !LOCAL
  INTEGER :: nS, nP, i
   

      CALL interp % GetNumberOfNodes( nS, nP )
      CALL myGeom % SetNumberOfNodes( nS, nP )
       
      ! Allocate space
      ALLOCATE( myGeom % dxds(0:nS,0:nP), myGeom % dxdp(0:nS,0:nP) )
      ALLOCATE( myGeom % dyds(0:nS,0:nP), myGeom % dydp(0:nS,0:nP) )
      ALLOCATE( myGeom % J(0:nS,0:nP) )
      ALLOCATE( myGeom % x(0:nS,0:nP), myGeom % y(0:nS,0:nP) )
      ALLOCATE( myGeom % xBound(0:max(nS,nP),1:nQuadEdges) )
      ALLOCATE( myGeom % yBound(0:max(nS,nP),1:nQuadEdges) )
      ALLOCATE( myGeom % nHat(0:max(nS,nP),1:nQuadEdges) )
     
      DO i = 0,max(nS,nP)

         CALL myGeom % nHat(i,south) % Build( nDims ) ! Allocate space for the south nHats
         CALL myGeom % nHat(i,east)  % Build( nDims ) ! Allocate space for the east nHats
         CALL myGeom % nHat(i,north) % Build( nDims ) ! Allocate space for the north nHats
         CALL myGeom % nHat(i,west)  % Build( nDims ) ! Allocate space for the west nHats

      ENDDO
      
      ! Generate the mesh locations, and the mapping metrics

      CALL myGeom % GenerateMesh( interp, myCurves )
      CALL myGeom % GenerateMetrics( interp, myCurves )
 
 END SUBROUTINE Build_MappedGeometry_2D
!
!
!
 SUBROUTINE Trash_MappedGeometry_2D( myGeom )
 ! S/R Trash
 !
 ! =============================================================================================== !
 ! DECLARATIONS
  IMPLICIT NONE
  CLASS(MappedGeometry_2D), INTENT(inout)  :: myGeom
  ! LOCAL
  INTEGER :: i

      DEALLOCATE( myGeom % dxds, myGeom % dxdp )
      DEALLOCATE( myGeom % dyds, myGeom % dydp)
      DEALLOCATE( myGeom % J, myGeom % x, myGeom % y )
      DEALLOCATE( myGeom % xBound, myGeom % yBound )
     
      DO i = 0, myGeom % nMax

         CALL myGeom % nHat(i,south) % Trash( )
         CALL myGeom % nHat(i,east)  % Trash( )
         CALL myGeom % nHat(i,north) % Trash( )
         CALL myGeom % nHat(i,west)  % Trash( )

      ENDDO

      DEALLOCATE( myGeom % nHat )

 
 END SUBROUTINE Trash_MappedGeometry_2D
!
!
!==================================================================================================!
!--------------------------------------- Accessors ------------------------------------------------!
!==================================================================================================!
!
!
 SUBROUTINE SetNumberOfNodes_MappedGeometry_2D( myGeom, nS, nP )
 ! S/R SetNumberOfNodes
 !
 ! =============================================================================================== !
 ! DECLARATIONS
  IMPLICIT NONE
  CLASS(MappedGeometry_2D), INTENT(inout) :: myGeom
  INTEGER, INTENT(in)                     :: nS, nP
  
     myGeom % nS   = nS
     myGeom % nP   = nP
     myGeom % nMax = max(nS,nP)
     
 END SUBROUTINE SetNumberOfNodes_MappedGeometry_2D
!
!
!
 SUBROUTINE GetNumberOfNodes_MappedGeometry_2D( myGeom, nS, nP )
 ! S/R GetNumberOfNodes
 !
 ! =============================================================================================== !
 ! DECLARATIONS
  IMPLICIT NONE
  CLASS(MappedGeometry_2D), INTENT(in) :: myGeom
  INTEGER, INTENT(out)                 :: nS, nP
  
     nS = myGeom % nS
     nP = myGeom % nP
     
 END SUBROUTINE GetNumberOfNodes_MappedGeometry_2D
!
!
!
 SUBROUTINE SetPositions_MappedGeometry_2D( myGeom, x, y )
 ! S/R SetPositions
 !
 ! =============================================================================================== !
 ! DECLARATIONS
  IMPLICIT NONE
  CLASS(MappedGeometry_2D), INTENT(inout) :: myGeom
  REAL(prec), INTENT(in)                  :: x(0:myGeom % nS, 0:myGeom % nP)
  REAL(prec), INTENT(in)                  :: y(0:myGeom % nS, 0:myGeom % nP)
  
     myGeom % x = x
     myGeom % y = y
     
 END SUBROUTINE SetPositions_MappedGeometry_2D
!
!
!
 SUBROUTINE GetPositions_MappedGeometry_2D( myGeom, x, y )
 ! S/R GetPositions
 !
 ! =============================================================================================== !
 ! DECLARATIONS
  IMPLICIT NONE
  CLASS(MappedGeometry_2D), INTENT(in) :: myGeom
  REAL(prec), INTENT(out)              :: x(0:myGeom % nS, 0:myGeom % nP)
  REAL(prec), INTENT(out)              :: y(0:myGeom % nS, 0:myGeom % nP)
  
     x = myGeom % x
     y = myGeom % y
     
 END SUBROUTINE GetPositions_MappedGeometry_2D
!
!
!
 SUBROUTINE SetPositionAtNode_MappedGeometry_2D( myGeom, x, y, i, j )
 ! S/R SetPositionAtNode
 !
 ! =============================================================================================== !
 ! DECLARATIONS
  IMPLICIT NONE
  CLASS(MappedGeometry_2D), INTENT(inout) :: myGeom
  REAL(prec), INTENT(in)                  :: x, y
  INTEGER, INTENT(in)                     :: i, j
  
     myGeom % x(i,j) = x
     myGeom % y(i,j) = y
     
 END SUBROUTINE SetPositionAtNode_MappedGeometry_2D
!
!
!
 SUBROUTINE GetPositionAtNode_MappedGeometry_2D( myGeom, x, y, i, j )
 ! S/R GetPositionAtNode
 !
 ! =============================================================================================== !
 ! DECLARATIONS
  IMPLICIT NONE
  CLASS(MappedGeometry_2D), INTENT(in) :: myGeom
  REAL(prec), INTENT(out)              :: x, y
  INTEGER, INTENT(in)                  :: i, j
  
     x = myGeom % x(i,j)
     y = myGeom % y(i,j)
     
 END SUBROUTINE GetPositionAtNode_MappedGeometry_2D
!
!
!
 SUBROUTINE SetJacobian_MappedGeometry_2D( myGeom, J )
 ! S/R SetJacobian
 !
 ! =============================================================================================== !
 ! DECLARATIONS
  IMPLICIT NONE
  CLASS(MappedGeometry_2D), INTENT(inout) :: myGeom
  REAL(prec), INTENT(in)                  :: J(0:myGeom % nS, 0:myGeom % nP)
  
     myGeom % J = J
     
 END SUBROUTINE SetJacobian_MappedGeometry_2D
!
!
!
 SUBROUTINE GetJacobian_MappedGeometry_2D( myGeom, J )
 ! S/R GetJacobian
 !
 ! =============================================================================================== !
 ! DECLARATIONS
  IMPLICIT NONE
  CLASS(MappedGeometry_2D), INTENT(in) :: myGeom
  REAL(prec), INTENT(out)              :: J(0:myGeom % nS, 0:myGeom % nP)
  
     J = myGeom % J
     
 END SUBROUTINE GetJacobian_MappedGeometry_2D
!
!
!
 SUBROUTINE SetJacobianAtNode_MappedGeometry_2D( myGeom, J, iS, iP )
 ! S/R SetJacobianAtNode
 !
 ! =============================================================================================== !
 ! DECLARATIONS
  IMPLICIT NONE
  CLASS(MappedGeometry_2D), INTENT(inout) :: myGeom
  REAL(prec), INTENT(in)                  :: J
  INTEGER, INTENT(in)                     :: iS, iP
  
     myGeom % J(iS,iP) = J
     
 END SUBROUTINE SetJacobianAtNode_MappedGeometry_2D
!
!
!
 SUBROUTINE GetJacobianAtNode_MappedGeometry_2D( myGeom, J, iS, iP )
 ! S/R GetJacobianAtNode
 !
 ! =============================================================================================== !
 ! DECLARATIONS
  IMPLICIT NONE
  CLASS(MappedGeometry_2D), INTENT(in) :: myGeom
  REAL(prec), INTENT(out)              :: J
  INTEGER, INTENT(in)                  :: iS, iP
  
     J = myGeom % J(iS,iP)
     
 END SUBROUTINE GetJacobianAtNode_MappedGeometry_2D
!
!
!
 SUBROUTINE SetCovariantMetrics_MappedGeometry_2D( myGeom, dxds, dxdp, dyds, dydp )
 ! S/R SetCovariantMetrics
 !
 ! =============================================================================================== !
 ! DECLARATIONS
  IMPLICIT NONE
  CLASS(MappedGeometry_2D), INTENT(inout) :: myGeom
  REAL(prec), INTENT(in)                  :: dxds(0:myGeom % nS, 0:myGeom % nP)
  REAL(prec), INTENT(in)                  :: dxdp(0:myGeom % nS, 0:myGeom % nP)
  REAL(prec), INTENT(in)                  :: dyds(0:myGeom % nS, 0:myGeom % nP)
  REAL(prec), INTENT(in)                  :: dydp(0:myGeom % nS, 0:myGeom % nP)
  
     myGeom % dxds = dxds
     myGeom % dxdp = dxdp
     myGeom % dyds = dyds
     myGeom % dydp = dydp
     
 END SUBROUTINE SetCovariantMetrics_MappedGeometry_2D
!
!
!
 SUBROUTINE GetCovariantMetrics_MappedGeometry_2D( myGeom, dxds, dxdp, dyds, dydp )
 ! S/R GetCovariantMetrics
 !
 ! =============================================================================================== !
 ! DECLARATIONS
  IMPLICIT NONE
  CLASS(MappedGeometry_2D), INTENT(in) :: myGeom
  REAL(prec), INTENT(out)              :: dxds(0:myGeom % nS, 0:myGeom % nP)
  REAL(prec), INTENT(out)              :: dxdp(0:myGeom % nS, 0:myGeom % nP)
  REAL(prec), INTENT(out)              :: dyds(0:myGeom % nS, 0:myGeom % nP)
  REAL(prec), INTENT(out)              :: dydp(0:myGeom % nS, 0:myGeom % nP)
  
     dxds = myGeom % dxds
     dxdp = myGeom % dxdp
     dyds = myGeom % dyds
     dydp = myGeom % dydp
     
 END SUBROUTINE GetCovariantMetrics_MappedGeometry_2D
!
!
!
 SUBROUTINE SetCovariantMetricsAtNode_MappedGeometry_2D( myGeom, dxds, dxdp, dyds, dydp, iS, iP )
 ! S/R SetCovariantMetricsAtNode
 !
 ! =============================================================================================== !
 ! DECLARATIONS
  IMPLICIT NONE
  CLASS(MappedGeometry_2D), INTENT(inout) :: myGeom
  REAL(prec), INTENT(in)                  :: dxds, dxdp, dyds, dydp
  INTEGER, INTENT(in)                     :: iS, iP
  
     myGeom % dxds(iS,iP) = dxds
     myGeom % dxdp(iS,iP) = dxdp
     myGeom % dyds(iS,iP) = dyds
     myGeom % dydp(iS,iP) = dydp
     
 END SUBROUTINE SetCovariantMetricsAtNode_MappedGeometry_2D
!
!
!
 SUBROUTINE GetCovariantMetricsAtNode_MappedGeometry_2D( myGeom, dxds, dxdp, dyds, dydp, iS, iP )
 ! S/R GetCovariantMetricsAtNode
 !
 ! =============================================================================================== !
 ! DECLARATIONS
  IMPLICIT NONE
  CLASS(MappedGeometry_2D), INTENT(in) :: myGeom
  REAL(prec), INTENT(out)              :: dxds, dxdp, dyds, dydp
  INTEGER, INTENT(in)                  :: iS, iP
  
     dxds = myGeom % dxds(iS,iP)
     dxdp = myGeom % dxdp(iS,iP)
     dyds = myGeom % dyds(iS,iP)
     dydp = myGeom % dydp(iS,iP)
     
 END SUBROUTINE GetCovariantMetricsAtNode_MappedGeometry_2D
!
!
!
 SUBROUTINE SetBoundaryLocation_MappedGeometry_2D( myGeom, x, y, iBound )
 ! S/R SetBoundaryLocation
 !
 ! =============================================================================================== !
 ! DECLARATIONS
  IMPLICIT NONE
  CLASS(MappedGeometry_2D), INTENT(inout) :: myGeom
  INTEGER, INTENT(in)                     :: iBound
  REAL(prec), INTENT(in)                  :: x(0:myGeom % nMax)
  REAL(prec), INTENT(in)                  :: y(0:myGeom % nMax)
  
     myGeom % xBound(:,iBound) = x
     myGeom % yBound(:,iBound) = y
     
 END SUBROUTINE SetBoundaryLocation_MappedGeometry_2D
!
!
!
 SUBROUTINE GetBoundaryLocation_MappedGeometry_2D( myGeom, x, y, iBound )
 ! S/R GetBoundaryLocation
 !
 ! =============================================================================================== !
 ! DECLARATIONS
  IMPLICIT NONE
  CLASS(MappedGeometry_2D), INTENT(in) :: myGeom
  INTEGER, INTENT(in)                  :: iBound
  REAL(prec), INTENT(out)              :: x(0:myGeom % nMax)
  REAL(prec), INTENT(out)              :: y(0:myGeom % nMax)
  
     x = myGeom % xBound(:,iBound)
     y = myGeom % yBound(:,iBound)
     
 END SUBROUTINE GetBoundaryLocation_MappedGeometry_2D
!
!
!
 SUBROUTINE SetBoundaryLocationAtNode_MappedGeometry_2D( myGeom, x, y, i, iBound )
 ! S/R SetBoundaryLocationAtNode
 !
 ! =============================================================================================== !
 ! DECLARATIONS
  IMPLICIT NONE
  CLASS(MappedGeometry_2D), INTENT(inout) :: myGeom
  REAL(prec), INTENT(in)                  :: x, y
  INTEGER, INTENT(in)                     :: i, iBound
  
     myGeom % xBound(i,iBound) = x
     myGeom % yBound(i,iBound) = y
     
 END SUBROUTINE SetBoundaryLocationAtNode_MappedGeometry_2D
!
!
!
 SUBROUTINE GetBoundaryLocationAtNode_MappedGeometry_2D( myGeom, x, y, i, iBound )
 ! S/R GetBoundaryLocationAtNode
 !
 ! =============================================================================================== !
 ! DECLARATIONS
  IMPLICIT NONE
  CLASS(MappedGeometry_2D), INTENT(in) :: myGeom
  REAL(prec), INTENT(out)              :: x, y
  INTEGER, INTENT(in)                  :: i, iBound
  
     x = myGeom % xBound(i,iBound)
     y = myGeom % yBound(i,iBound)
     
 END SUBROUTINE GetBoundaryLocationAtNode_MappedGeometry_2D
!
!
!
 SUBROUTINE SetBoundaryNormalAtNode_MappedGeometry_2D( myGeom, dir, i, iBound )
 ! S/R SetBoundaryNormalAtNode
 !
 ! =============================================================================================== !
 ! DECLARATIONS
  IMPLICIT NONE
  CLASS(MappedGeometry_2D), INTENT(inout) :: myGeom
  REAL(prec), INTENT(in)                  :: dir(1:2)
  INTEGER, INTENT(in)                     :: i, iBound
  
     CALL myGeom % nHat(i,iBound) % SetDirection( dir(1:2) )

     CALL myGeom % nHat(i,iBound) % SetLength( SQRT(DOT_PRODUCT(dir,dir)) )
     
 END SUBROUTINE SetBoundaryNormalAtNode_MappedGeometry_2D
!
!
!
 SUBROUTINE GetBoundaryNormalAtNode_MappedGeometry_2D( myGeom, dir, length, i, iBound )
 ! S/R GetBoundaryNormalAtNode
 !
 ! =============================================================================================== !
 ! DECLARATIONS
  IMPLICIT NONE
  CLASS(MappedGeometry_2D), INTENT(in) :: myGeom
  REAL(prec), INTENT(out)              :: dir(1:2), length
  INTEGER, INTENT(in)                  :: i, iBound
  
     CALL myGeom % nHat(i,iBound) % GetDirection( dir(1:2) )

     CALL myGeom % nHat(i,iBound) % GetLength( length )
     
 END SUBROUTINE GetBoundaryNormalAtNode_MappedGeometry_2D
!
!
!==================================================================================================!
!--------------------------------- Type Specific Routines -----------------------------------------!
!==================================================================================================!
!
!
 SUBROUTINE GenerateMesh_MappedGeometry_2D( myGeom, interp, theCurves )
 ! S/R GenerateMesh
 !
 ! =============================================================================================== !
 ! DECLARATIONS
  IMPLICIT NONE
  CLASS( MappedGeometry_2D ), INTENT(inout) :: myGeom
  TYPE( Lagrange_2D ), INTENT(in)           :: interp
  TYPE( Curve_2D ), INTENT(in)              :: theCurves(1:4)
  ! Local
  INTEGER    :: iS, iP, nS, nP
  REAL(prec) :: s, p, x(1:2)
   
      CALL interp % GetNumberOfNodes( nS, nP )
      
      DO iS = 0, nS
 
         DO iP = 0,nP
      
            ! Get the interpolation point where we need to generate the mesh
            CALL interp % GetNode( iS, iP, s, p ) 
           
            ! Mesh location is computed via transfinite interpolation
            x = TransfiniteInterpolation( theCurves, s, p )
            
            CALL myGeom % SetPositionAtNode( x(1), x(2), iS, iP )
            
         ENDDO ! iP

      ENDDO ! iS
      
      ! Do the boundary locations
      DO iS = 0, nS
      
         CALL interp % GetNode( iS, 0, s, p )
         p = -ONE  ! south boundary
         
         x = TransfiniteInterpolation( theCurves, s, p )
         
         ! Set the south boundary
         CALL myGeom % SetBoundaryLocationAtNode( x(1), x(2), iS, 1 )
         
         p = ONE  ! north boundary
         
         x = TransfiniteInterpolation( theCurves, s, p )
         
         ! Set the north boundary
         CALL myGeom % SetBoundaryLocationAtNode( x(1), x(2), iS, 3 )
      
      ENDDO
      
      ! Do the boundary locations
      DO iP = 0, nP
      
         CALL interp % GetNode( 0, iP, s, p )
         s = -ONE  ! west boundary
         
         x = TransfiniteInterpolation( theCurves, s, p )
         
         ! Set the west boundary
         CALL myGeom % SetBoundaryLocationAtNode( x(1), x(2), iP, 4 )
         
         s = ONE  ! east boundary
         
         x = TransfiniteInterpolation( theCurves, s, p )
         
         ! Set the east boundary
         CALL myGeom % SetBoundaryLocationAtNode( x(1), x(2), iP, 2 )
      
      ENDDO

 END SUBROUTINE GenerateMesh_MappedGeometry_2D
!
!
!
 SUBROUTINE GenerateMetrics_MappedGeometry_2D( myGeom, interp, theCurves )
 ! S/R GenerateMetrics
 !
 ! =============================================================================================== !
 ! DECLARATIONS
  IMPLICIT NONE
  CLASS( MappedGeometry_2D ), INTENT(inout) :: myGeom
  TYPE( Lagrange_2D ), INTENT(in)           :: interp
  TYPE( Curve_2D ), INTENT(in)              :: theCurves(1:4)
  ! Local
  INTEGER    :: iS, iP, nS, nP
  REAL(prec) :: s, p, covT(1:2,1:2), J, signJ
   
      CALL interp % GetNumberOfNodes( nS, nP )
      
      DO iS = 0, nS
 
         DO iP = 0,nP
      
            ! Get the interpolation point where we need to generate the mesh
            CALL interp % GetNode( iS, iP, s, p ) 
            
            ! Mesh location is computed via transfinite interpolation
            covT = TransfiniteInterpolationMetrics( theCurves, s, p )
            
            CALL myGeom % SetCovariantMetricsAtNode( covT(1,1), & ! dxds
                                                     covT(1,2), & ! dxdp
                                                     covT(2,1), & ! dyds
                                                     covT(2,2), & ! dydp 
                                                     iS, iP )
            
            J = Determinant( covT, 2 )
            
            CALL myGeom % SetJacobianAtNode( J, iS, iP )

         ENDDO ! iP

      ENDDO ! iS
      
      ! Do the boundary locations
      DO iS = 0, nS
      
         CALL interp % GetNode( iS, 0, s, p )
         p = -ONE  ! south boundary
         
         covT = TransfiniteInterpolationMetrics( theCurves, s, p )
         
         J = Determinant( covT, 2 )

         signJ = abs(J)/J                    

         ! Setting southern boundary normal                 -dyds    ,    dxds
         CALL myGeom % SetBoundaryNormalAtNode( -signJ*(/ -covT(2,1), covT(1,1) /), iS, south )
         CALL myGeom % nHat(iS, south) % Normalize( )
          
         p = ONE  ! north boundary
         
         covT = TransfiniteInterpolationMetrics( theCurves, s, p )
         
         J = Determinant( covT, 2 )

         signJ = abs(J)/J                    

         ! Setting northern boundary normal                -dyds    ,    dxds
         CALL myGeom % SetBoundaryNormalAtNode( signJ*(/ -covT(2,1), covT(1,1) /), iS, north )
         CALL myGeom % nHat(iS, north) % Normalize( )
         
      ENDDO
      
      ! Do the boundary locations
      DO iP = 0, nP
      
         CALL interp % GetNode( 0, iP, s, p )
         s = -ONE  ! west boundary
         
         covT = TransfiniteInterpolationMetrics( theCurves, s, p )
         
         J = Determinant( covT, 2 )

         signJ = abs(J)/J                    

         ! Setting western boundary normal                  dydp   ,   -dxdp
         CALL myGeom % SetBoundaryNormalAtNode( -signJ*(/ covT(2,2), -covT(1,2) /), iP, west )
         CALL myGeom % nHat(iP, west) % Normalize( )
         
         s = ONE  ! east boundary
         
         covT = TransfiniteInterpolationMetrics( theCurves, s, p )
         
         J = Determinant( covT, 2 )

         signJ = abs(J)/J                    

         ! Setting eastern boundary normal                  dydp   ,   -dxdp
         CALL myGeom % SetBoundaryNormalAtNode( signJ*(/ covT(2,2), -covT(1,2) /), iP, east )
         CALL myGeom % nHat(iP, east) % Normalize( )
      
      ENDDO

 END SUBROUTINE GenerateMetrics_MappedGeometry_2D
!
!
! 
  SUBROUTINE GenerateBoundaryMetrics_MappedGeometry_2D( myGeom, interp  )
 ! S/R GenerateBoundaryMetrics
 !
 ! =============================================================================================== !
 ! DECLARATIONS
  IMPLICIT NONE
  CLASS( MappedGeometry_2D ), INTENT(inout) :: myGeom
  TYPE( Lagrange_2D ), INTENT(in)           :: interp
  ! Local
  INTEGER    :: iS, iP, nS, nP
  REAL(prec) :: s, p, covT(1:2,1:2), J, signJ
   
      CALL interp % GetNumberOfNodes( nS, nP )
      
      ! Do the boundary locations
      DO iS = 0, nS
      
         CALL interp % GetNode( iS, 0, s, p )
         p = -ONE  ! south boundary
         
         CALL myGeom % CalculateMetrics( interp, s, p, covT )
         
         J = Determinant( covT, 2 )

         signJ = abs(J)/J                    

         ! Setting southern boundary normal                 -dyds    ,    dxds
         CALL myGeom % SetBoundaryNormalAtNode( -signJ*(/ -covT(2,1), covT(1,1) /), iS, south )
         CALL myGeom % nHat(iS, south) % Normalize( )
          
         p = ONE  ! north boundary
         
         CALL myGeom % CalculateMetrics( interp, s, p, covT )
         
         J = Determinant( covT, 2 )

         signJ = abs(J)/J                    

         ! Setting northern boundary normal                -dyds    ,    dxds
         CALL myGeom % SetBoundaryNormalAtNode( signJ*(/ -covT(2,1), covT(1,1) /), iS, north )
         CALL myGeom % nHat(iS, north) % Normalize( )
         
      ENDDO
      
      ! Do the boundary locations
      DO iP = 0, nP
      
         CALL interp % GetNode( 0, iP, s, p )
         s = -ONE  ! west boundary
         
         CALL myGeom % CalculateMetrics( interp, s, p, covT )
         
         J = Determinant( covT, 2 )

         signJ = abs(J)/J                    

         ! Setting western boundary normal                  dydp   ,   -dxdp
         CALL myGeom % SetBoundaryNormalAtNode( -signJ*(/ covT(2,2), -covT(1,2) /), iP, west )
         CALL myGeom % nHat(iP, west) % Normalize( )
         
         s = ONE  ! east boundary
         
         CALL myGeom % CalculateMetrics( interp, s, p, covT )
         
         J = Determinant( covT, 2 )

         signJ = abs(J)/J                    

         ! Setting eastern boundary normal                  dydp   ,   -dxdp
         CALL myGeom % SetBoundaryNormalAtNode( signJ*(/ covT(2,2), -covT(1,2) /), iP, east )
         CALL myGeom % nHat(iP, east) % Normalize( )
      
      ENDDO

 END SUBROUTINE GenerateBoundaryMetrics_MappedGeometry_2D
!
!
!
 SUBROUTINE CalculateLocation_MappedGeometry_2D( myGeom, interp, s, p, x, y )
 ! S/R CalculateLocation
 !  Description :
 ! 
 !     This subroutine calculates the physical coordinates (x,y) given the computational coordinates
 !     (s,p).
 ! =============================================================================================== !
  IMPLICIT NONE
  CLASS( MappedGeometry_2D ), INTENT(in) :: myGeom
  TYPE( Lagrange_2D ), INTENT(in)        :: interp
  REAL(prec), INTENT(in)                 :: s, p
  REAL(prec), INTENT(out)                :: x, y
  
  
     x = interp % EvaluateInterpolant( s, p, myGeom % x )
     y = interp % EvaluateInterpolant( s, p, myGeom % y )
  
 END SUBROUTINE CalculateLocation_MappedGeometry_2D
!
!
!
 SUBROUTINE CalculateMetrics_MappedGeometry_2D( myGeom, interp, s, p, covT )
 ! S/R CalculateMetrics
 !  Description :
 ! 
 !     This subroutine calculates the covariant metric terms given the computational coordinates
 !     (s,p).
 !     covT(1,1) -> dxds
 !     covT(1,2) -> dxdp
 !     covT(2,1) -> dyds
 !     covT(2,2) -> dydp
 ! =============================================================================================== !
  IMPLICIT NONE
  CLASS( MappedGeometry_2D ), INTENT(in) :: myGeom
  TYPE( Lagrange_2D ), INTENT(in)        :: interp
  REAL(prec), INTENT(in)                 :: s, p
  REAL(prec), INTENT(out)                :: covT(1:2,1:2)

     covT(1,1:2) = interp % EvaluateDerivative( s, p, myGeom % x )
     covT(2,1:2) = interp % EvaluateDerivative( s, p, myGeom % y )
       
 END SUBROUTINE CalculateMetrics_MappedGeometry_2D
!
!
!
 SUBROUTINE CalculateComputationalCoordinates_MappedGeometry_2D( myGeom, interp, x, y, s, p, success )
 ! S/R CalculateComputationalCoordinates
 !  Description :
 ! 
 !     Given the physical coordinates (x*,y*), the computational coordinates are calculated using 
 !     Newton's method for root finding to solve
 !
 !     x* = x(s,p),
 !     y* = y(s,p)
 !
 !     for s and p.
 ! =============================================================================================== !
  IMPLICIT NONE
  CLASS( MappedGeometry_2D ), INTENT(in) :: myGeom
  TYPE( Lagrange_2D ), INTENT(in)        :: interp
  REAL(prec), INTENT(in)                 :: x, y
  REAL(prec), INTENT(out)                :: s, p
  LOGICAL, INTENT(out), OPTIONAL         :: success
  ! LOCAL
  REAL(prec) :: dr(1:2), ds(1:2), A(1:2,1:2), Ainv(1:2,1:2)
  REAL(prec) :: thisX, thisY, thisS, thisP, resi
  INTEGER    :: i 

     thisS = ZERO ! Initial guess is at the origin
     thisP = ZERO
     
     IF( PRESENT(success) )THEN
        success = .FALSE.
     ENDIF
     
     DO i = 1, newtonMax
     
        ! Calculate the physical coordinate associated with the computational coordinate guess
        CALL myGeom % CalculateLocation( interp, thisS, thisP, thisX, thisY )
     
        ! Calculate the residual
        dr(1) = x - thisX
        dr(2) = y - thisY  
        resi = SQRT( DOT_PRODUCT( dr, dr ) )
     
        IF( resi < newtonTolerance )THEN
           s = thisS
           p = thisP
           IF( PRESENT(success) )THEN
              success = .TRUE.
           ENDIF
           RETURN
        ENDIF
        
        CALL myGeom % CalculateMetrics( interp, thisS, thisP, A ) ! Calculate the covariant metric tensor
     
        Ainv = Invert( A ) ! Invert the covariant metric tensor.
                           ! This matrix is ivertable as long as the Jacobian is non-zero.
        
        ds = MATMUL( Ainv, dr ) ! calculate the correction in the computational coordinate
        thisS = thisS + ds(1)
        thisP = thisP + ds(2)
     
     ENDDO
     
     ! Calculate the residual
     dr(1) = x - thisX
     dr(2) = y - thisY  
     resi = SQRT( DOT_PRODUCT( dr, dr ) )
     IF( resi < newtonTolerance )THEN
        s = thisS
        p = thisP
        IF( PRESENT(success) )THEN
           success = .TRUE.
        ENDIF
        RETURN
     ELSE
        s = fillValue
        p = fillValue
        PRINT*,'Module MappedGeometryClass_2D.f90 : S/R CalculateComputationalCoordinates :'
        PRINT*,'Search for coordinates failed. Final residual norm :', resi
        RETURN
    ENDIF
      
     
 END SUBROUTINE CalculateComputationalCoordinates_MappedGeometry_2D
!
!
!
 SUBROUTINE ScaleGeometry_MappedGeometry_2D( myGeom, interp, xScale, yScale )
 !
 !
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS( MappedGeometry_2D ), INTENT(inout) :: myGeom
   TYPE( Lagrange_2D ), INTENT(in)           :: interp
   REAL(prec), INTENT(in)                    :: xScale, yScale
   !
   
         myGeom % x = xScale*( myGeom % x )
         myGeom % y = yScale*( myGeom % y )
         myGeom % xBound = xScale*( myGeom % xBound )
         myGeom % yBound = yScale*( myGeom % yBound )

         myGeom % dxds = xScale*( myGeom % dxds )
         myGeom % dxdp = xScale*( myGeom % dxdp )
         myGeom % dyds = yScale*( myGeom % dyds )
         myGeom % dydp = yScale*( myGeom % dydp )
          
         myGeom % J = xScale*yScale*( myGeom % J )

         ! Update the boundary metrics -- normals and normal lengths
         CALL myGeom % GenerateBoundaryMetrics( interp  )
         
 END SUBROUTINE ScaleGeometry_MappedGeometry_2D
!
!
!==================================================================================================!
!------------------------------------ File I/O Routines -------------------------------------------!
!==================================================================================================!
!
!
 SUBROUTINE WriteTecplot_MappedGeometry_2D( myGeom, filename )
 ! S/R WriteTecplot
 !  Description :
 ! 
 ! =============================================================================================== !
  IMPLICIT NONE
  CLASS( MappedGeometry_2D ), INTENT(in) :: myGeom
  CHARACTER(*), INTENT(in), OPTIONAL     :: filename  
  ! Local
  INTEGER :: iX, iY, nX, nY, fUnit
  REAL(prec) :: x, y, dxds, dxdp, dyds, dydp, J
  
    CALL myGeom % GetNumberOfNodes( nX, nY )
    
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

    WRITE(fUnit,*) 'VARIABLES = "X", "Y", "Jacobian", "dxds", "dxdp", "dyds", "dydp" '
    
    WRITE(fUnit,*)  'ZONE T="el0", I=',nX+1,', J=', nY+1,',F=POINT'
    DO iX = 0, nX
       DO iY = 0, nY
       
          CALL myGeom % GetPositionAtNode( x, y, iX, iY )
          CALL myGeom % GetJacobianAtNode( J, iX, iY )
          CALL myGeom % GetCovariantMetricsAtNode( dxds, dxdp, dyds, dydp, iX, iY )

          WRITE(fUnit,*)  x, y, J, dxds, dxdp, dyds, dydp

       ENDDO
    ENDDO

    CLOSE(UNIT=fUnit)

    RETURN

 END SUBROUTINE WriteTecplot_MappedGeometry_2D
 
 
END MODULE MappedGeometryClass_2D
