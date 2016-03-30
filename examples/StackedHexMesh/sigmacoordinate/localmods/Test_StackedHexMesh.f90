! Test_StackedHexMesh.f90 ( new with v2.1 - 29 March 2016)
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

PROGRAM Test_StackedHexMesh

USE ModelPrecision
USE ModelFlags
USE ConstantsDictionary

USE NodalStorage_2D_Class
USE NodalStorage_3D_Class

USE QuadMeshClass
USE StackedHexMeshClass

USE MeshParamsClass

IMPLICIT NONE

   TYPE(NodalStorage_2D) :: cg2dStorage
   TYPE(NodalStorage_3D) :: dg3dStorage

   TYPE(QuadMesh)        :: mesh2d
   TYPE(StackedHexMesh)  :: mesh3d

   TYPE(MeshParams)      :: params
   
   REAL(prec), ALLOCATABLE :: h(:,:,:)

      CALL params % Build( )


      CALL cg2dStorage % Build( N = params % geomPolyDeg, &
                                quadrature = GAUSS_LOBATTO, &
                                approxForm = CG  )

      CALL dg3dStorage % Build( N = params % polyDeg, &
                                quadrature = GAUSS, &
                                approxForm = DG  )

      IF( TRIM( params % SpecMeshFile ) == nada )THEN
         PRINT*,' Loading default mesh.'
         CALL  mesh2d % LoadDefaultMesh( cg2dStorage % interp, &
                                         params % nXelem, &
                                         params % nYelem )
      ELSE
      ! Builds the lateral mesh
         PRINT*, 'Reading mesh from '//trim(params % SpecMeshFile)//'.'
         CALL mesh2d % ReadSpecMeshFile( cg2dStorage % interp, params % SpecMeshFile )
      ENDIF


      ALLOCATE( h(0:params % geomPolyDeg, 0:params % geomPolyDeg, 1:mesh2d % nElems) )
      h = ZERO

      CALL SetBathymetry( mesh2d, h, params % geomPolyDeg )

      CALL mesh2d % ScaleTheMesh( cg2dStorage % interp, &
                                  params % xScale, &
                                  params % yScale )     

      CALL mesh3d % TerrainFollowingMesh( dg3dStorage % interp, mesh2d, &
                                          cg2dStorage % interp, &
                                          params % nZElem, h  )

      CALL mesh2d % WriteTecplot( 'mesh2d' )
      CALL mesh3d % WriteTecplot( 'mesh3d' )


      ! ///// Cleanup ///// !
      CALL mesh2d % Trash( )
      CALL mesh3d % Trash( )
      CALL cg2dStorage % Trash( )
      CALL dg3dStorage % Trash( )
      DEALLOCATE( h )

CONTAINS

 SUBROUTINE SetBathymetry( mesh, h, N )
 !
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   TYPE(QuadMesh), INTENT(in) :: mesh
   INTEGER, INTENT(in)        :: N
   REAL(prec), INTENT(inout)  :: h(0:N,0:N,1:mesh % nElems)
   ! Local
   INTEGER :: i, j, iEl, nEl
   REAL(prec) :: x, y

      nEl = mesh % nElems
      DO iEl = 1, nEl
         DO j = 0, N
            DO i = 0, N
               CALL mesh % GetPositionAtNode( iEl, x, y, i, j )
               h(i,j,iEl) = exp( -( (x-HALF)**2 + (y-HALF)**2 ) )
            ENDDO
         ENDDO
      ENDDO

 END SUBROUTINE SetBathymetry
      

END PROGRAM Test_StackedHexMesh
