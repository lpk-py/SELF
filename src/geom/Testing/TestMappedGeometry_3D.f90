PROGRAM TestMappedGeometry_3D
! 
! TestMappedGeometry_3D.f90 (25/01/2016, new with v 2.1)
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
!
! Module History
!
! (v 2.1 - 25 January, 2016)
!
!
!
! ================================================================================================ !


 ! src/common/
 USE ModelPrecision
 USE ConstantsDictionary
 ! src/interp/
 USE Lagrange_3D_Class
 ! src/nodal/
 USE NodalStorage_3D_Class
 ! src/geom/
 USE SurfaceClass_3D
 USE GeometryParamsClass
 USE MappedGeometryClass_3D
 
 IMPLICIT NONE
 
 TYPE( NodalStorage_3D )   :: nodal
 TYPE( MappedGeometry_3D ) :: geometry, plottingGeometry
 TYPE( Surface_3D )        :: boundarySurfs(1:nHexFaces)
 TYPE( GeometryParams )    :: params
 REAL(prec)                :: s, p, q
 INTEGER                   :: i
 LOGICAL                   :: successful
 
 
 
    ! Read in the parameters
    CALL params % Build( )
    
    ! Build an interpolant
    CALL nodal % Build( N = params % polyDeg, &
                        M = params % polyDeg, &
                        P = params % polyDeg, &
                        quadrature = GAUSS_LOBATTO,  &
                        approxForm = CG )
                        
    ! Build the boundary curves
    CALL ConstructSurfaces( boundarySurfs, nodal % interp )
    
    ! Build the Geometry
    CALL geometry % Build( nodal % interp, boundarySurfs )
    
    ! Write to a tecplot file
    CALL geometry % WriteTecplot( )
    
    ! Test out the computational coordinate solver
    CALL geometry % CalculateComputationalCoordinates( nodal % interp, HALF, HALF, -0.5_prec, &
                                                       s, p, q, successful )
    
    IF( successful )THEN
       PRINT*,'Computational coordinate search successful '
       PRINT*, 'Coordinates (s,p,q) : (', s,',',p,',',q,')'
    ENDIF
    
    CALL geometry % Trash( )
    DO i = 1, nHexFaces
       CALL boundarySurfs(i) % Trash()
    ENDDO
 CONTAINS
 
 SUBROUTINE ConstructSurfaces( boundSurfaces, interp )
    TYPE( Surface_3D ), INTENT(out) :: boundSurfaces(1:nHexFaces)
    TYPE( Lagrange_3D ), INTENT(in) :: interp
    !
    REAL(prec) :: r(0:interp % nS, 0:interp % nP, 1:3)
    REAL(prec) :: s(0:interp % nS), p(0:interp % nP), q(0:interp % nQ)
    REAL(prec) :: xb, yb, zb, xt, yt, zt
    INTEGER    :: nS, nP, nQ, iS, iP, iQ
    
       CALL interp % GetNodes( s, p, q )
       CALL interp % GetNumberOfNodes( nS, nP, nQ )
       
       r = TopSurface( s, p, nS, nP )
       CALL boundSurfaces(top) % Build( r(:,:,1), r(:,:,2), r(:,:,3), s, p )

       r = BottomSurface( s, p, nS, nP )
       CALL boundSurfaces(bottom) % Build( r(:,:,1), r(:,:,2), r(:,:,3), s, p )
       
       ! The remaining sides (south,north,west,east) can be generated however we like. Here, we 
       ! choose to build an extrusion mesh where these remaining four sides are built by linearly 
       ! interpolating between the top and bottom surfaces. In this implementation, it is assumed
       ! that Gauss-Lobatto quadrature is used. The main purpose of Gauss-Lobatto quadrature is that
       ! the interpolation nodes include the boundaries of the computational domain (whereas Gauss
       ! quadrature does not).
       
       ! South Face
       DO iQ = 0, nQ
          DO iS = 0, nS
             r(iS,iQ,1) = HALF*(s(iS)+ONE) ! x
             r(iS,iQ,2) = ZERO             ! y
             
             CALL boundSurfaces(bottom) % GetPositionAtNode( xb, yb, zb, iS, 0 ) ! south-bottom
             CALL boundSurfaces(top)    % GetPositionAtNode( xt, yt, zt, iS, 0 ) ! south-top
             r(iS,iQ,3) = HALF*(zt-zb)*(q(iQ) + ONE) + zb
          ENDDO
       ENDDO
       CALL boundSurfaces(south) % Build( r(:,:,1), r(:,:,2), r(:,:,3), s, q )
       
       ! North Face
       DO iQ = 0, nQ
          DO iS = 0, nS
             r(iS,iQ,1) = HALF*(s(iS)+ONE) ! x
             r(iS,iQ,2) = ONE             ! y
             
             CALL boundSurfaces(bottom) % GetPositionAtNode( xb, yb, zb, iS, nP ) ! north-bottom
             CALL boundSurfaces(top)    % GetPositionAtNode( xt, yt, zt, iS, nP ) ! north-top
             r(iS,iQ,3) = HALF*(zt-zb)*(q(iQ) + ONE) + zb
          ENDDO
       ENDDO
       CALL boundSurfaces(north) % Build( r(:,:,1), r(:,:,2), r(:,:,3), s, q )
       
       ! West Face
       DO iQ = 0, nQ
          DO iP = 0, nP
             r(iP,iQ,1) = ZERO             ! x
             r(iP,iQ,2) = HALF*(p(iP)+ONE) ! y
             
             CALL boundSurfaces(bottom) % GetPositionAtNode( xb, yb, zb, 0, iP ) ! west-bottom
             CALL boundSurfaces(top)    % GetPositionAtNode( xt, yt, zt, 0, iP ) ! west-top
             r(iP,iQ,3) = HALF*(zt-zb)*(q(iQ) + ONE) + zb
          ENDDO
       ENDDO
       CALL boundSurfaces(west) % Build( r(:,:,1), r(:,:,2), r(:,:,3), p, q )
       
       ! East Face
       DO iQ = 0, nQ
          DO iP = 0, nP
             r(iP,iQ,1) = ONE              ! x
             r(iP,iQ,2) = HALF*(p(iP)+ONE) ! y
             
             CALL boundSurfaces(bottom) % GetPositionAtNode( xb, yb, zb, nS, iP ) ! north-bottom
             CALL boundSurfaces(top)    % GetPositionAtNode( xt, yt, zt, nS, iP ) ! north-top
             r(iP,iQ,3) = HALF*(zt-zb)*(q(iQ) + ONE) + zb
          ENDDO
       ENDDO
       CALL boundSurfaces(east) % Build( r(:,:,1), r(:,:,2), r(:,:,3), p, q )
       
 END SUBROUTINE ConstructSurfaces

 FUNCTION TopSurface( s, p, N, M ) RESULT( r )
    INTEGER    :: N, M
    REAL(prec) :: s(0:N), p(0:M)
    REAL(prec) :: r(0:N,0:M,1:3)
    !LOCAL
    INTEGER :: i, j
       
       DO j = 0, M
          DO i = 0, N
             r(i,j,1) = HALF*(s(i) + ONE)
             r(i,j,2) = HALF*(p(j) + ONE)
             r(i,j,3) = ZERO
          ENDDO
       ENDDO
       
 END FUNCTION TopSurface
! 
!
!
 FUNCTION BottomSurface( s, p, N, M ) RESULT( r )
    INTEGER    :: N, M
    REAL(prec) :: s(0:N), p(0:M)
    REAL(prec) :: r(0:N,0:M,1:3)
    !LOCAL
    INTEGER :: i, j
       
       DO j = 0, M
          DO i = 0, N
             r(i,j,1) = HALF*(s(i) + ONE)
             r(i,j,2) = HALF*(p(j) + ONE)
             r(i,j,3) = -( 0.1_prec + exp( -HALF*( (r(i,j,1)-HALF)**2 + (r(i,j,2)-HALF)**2 )/(0.3_prec**2)  ) )
          ENDDO
       ENDDO
       
 END FUNCTION BottomSurface
! 
!
!
 
END PROGRAM TestMappedGeometry_3D
