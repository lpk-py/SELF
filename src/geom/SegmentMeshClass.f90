! SegmentMeshClass.f90
! 
! Copyright 2015 Joseph Schoonover <schoonover.numerics@gmail.com>
! 
! SegmentMeshClass.f90 is part of the Spectral Element Libraries in Fortran (SELF).
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
 
MODULE SegmentMeshClass
! ========================================= Logs ================================================= !
!2016-05-11  Joseph Schoonover  schoonover.numerics@gmail.com 
!
! //////////////////////////////////////////////////////////////////////////////////////////////// ! 
! src/common/
USE ModelPrecision
USE ConstantsDictionary
USE ModelFlags
! src/interp/
USE Lagrange_1D_Class
! src/geom/
USE MappedGeometryClass_1D


IMPLICIT NONE


! All properties are left as public in version 2.1 in order to allow the programmer the ability to
! have direct access to the attributes of the mesh class. This way, it can later be determined if
! additional accessor routines would be necessary.
 
    TYPE SegmentMesh 
       INTEGER                                :: nElems
       TYPE( MappedGeometry_1D ), ALLOCATABLE :: elements(:)
       INTEGER                                :: sideMap(1:2)
       CONTAINS

       PROCEDURE :: Build => Build_SegmentMesh
       PROCEDURE :: Trash => Trash_SegmentMesh
       
       PROCEDURE :: SetNumberOfElements => SetNumberOfElements_SegmentMesh
       PROCEDURE :: GetNumberOfElements => GetNumberOfElements_SegmentMesh
      
       ! MappedGeometry_2D Wrapper Routines
       PROCEDURE :: SetNumberOfInternalNodes => SetNumberOfInternalNodes_SegmentMesh
       PROCEDURE :: GetNumberOfInternalNodes => GetNumberOfInternalNodes_SegmentMesh
       PROCEDURE :: SetPositions => SetPositions_SegmentMesh
       PROCEDURE :: GetPositions => GetPositions_SegmentMesh
       PROCEDURE :: SetPositionAtNode => SetPositionAtNode_SegmentMesh
       PROCEDURE :: GetPositionAtNode => GetPositionAtNode_SegmentMesh
       PROCEDURE :: SetJacobian => SetJacobian_SegmentMesh
       PROCEDURE :: GetJacobian => GetJacobian_SegmentMesh
       PROCEDURE :: SetJacobianAtNode => SetJacobianAtNode_SegmentMesh
       PROCEDURE :: GetJacobianAtNode => GetJacobianAtNode_SegmentMesh
       PROCEDURE :: SetBoundaryLocation => SetBoundaryLocation_SegmentMesh
       PROCEDURE :: GetBoundaryLocation => GetBoundaryLocation_SegmentMesh

       PROCEDURE :: ScaleTheMesh => ScaleTheMesh_SegmentMesh
       PROCEDURE :: LoadDefaultMesh => LoadDefaultMesh_SegmentMesh
       
       PROCEDURE :: WriteTecplot => WriteTecplot_SegmentMesh
       
    END TYPE SegmentMesh

 INTEGER, PRIVATE, PARAMETER    :: nXElemDefault = 5
 INTEGER, PRIVATE, PARAMETER    :: BoundaryFlagDefault = NO_NORMAL_FLOW
 REAL(prec), PRIVATE, PARAMETER :: dXDefault = ONE/(nXElemDefault)


 CONTAINS
!
!
!==================================================================================================!
!------------------------------- Manual Constructors/Destructors ----------------------------------!
!==================================================================================================!
!
!
 SUBROUTINE Build_SegmentMesh( mySegmentMesh, nElems, nS )
 ! S/R Build
 ! 
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS(SegmentMesh), INTENT(out) :: mySegmentMesh
   INTEGER, INTENT(in)             :: nElems, nS

      mySegmentMesh % sideMap(1:2) = (/ 0, nS /)
      mySegmentMesh % nElems = nElems

      ALLOCATE( mySegmentMesh % elements(1:nElems) )

 END SUBROUTINE Build_SegmentMesh
!
!
!
 SUBROUTINE Trash_SegmentMesh( mySegmentMesh )
 ! S/R Trash
 ! 
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS(SegmentMesh), INTENT(inout) :: mySegmentMesh
  ! LOCAL
   INTEGER :: iEl

      DO iEl = 1, mySegmentMesh % nElems
         CALL mySegmentMesh % elements(iEl) % TRASH( )
      ENDDO

      DEALLOCATE( mySegmentMesh % elements )
      
 END SUBROUTINE Trash_SegmentMesh
!
!
!==================================================================================================!
!--------------------------------------- Accessors ------------------------------------------------!
!==================================================================================================!
!
! 
 SUBROUTINE SetNumberOfElements_SegmentMesh( mySegmentMesh, nElems )
 ! S/R SetNumberOfElements
 ! 
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS(SegmentMesh), INTENT(inout) :: mySegmentMesh
   INTEGER, INTENT(in)               :: nElems
   
      mySegmentMesh % nElems = nElems
      
 END SUBROUTINE SetNumberOfElements_SegmentMesh
!
!
!
 FUNCTION GetNumberOfElements_SegmentMesh( mySegmentMesh ) RESULT( nElems )
 ! FUNCTION GetNumberOfElements
 ! 
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS(SegmentMesh) :: mySegmentMesh
   INTEGER            :: nElems
   
      nElems = mySegmentMesh % nElems
      
 END FUNCTION GetNumberOfElements_SegmentMesh
!
! ---------------------------- MappedGeometryClass_1D Wrappers ----------------------------------- !
!
SUBROUTINE SetNumberOfInternalNodes_SegmentMesh( mySegmentMesh, iEl, nS )
 ! S/R SetNumberOfInternalNodes
 !
 ! =============================================================================================== !
 ! DECLARATIONS
  IMPLICIT NONE
  CLASS(SegmentMesh), INTENT(inout) :: mySegmentMesh
  INTEGER, INTENT(in)               :: iEl
  INTEGER, INTENT(in)               :: nS
  
     CALL mySegmentMesh % elements(iEl) % SetNumberOfNodes( nS )
     
 END SUBROUTINE SetNumberOfInternalNodes_SegmentMesh
!
!
!
 FUNCTION GetNumberOfInternalNodes_SegmentMesh( mySegmentMesh, iEl ) RESULT( nS )
 ! FUNCTION GetNumberOfInternalNodes
 !
 ! =============================================================================================== !
 ! DECLARATIONS
  IMPLICIT NONE
  CLASS(SegmentMesh) :: mySegmentMesh
  INTEGER            :: iEl
  INTEGER            :: nS
  
     nS = mySegmentMesh % elements(iEl) % GetNumberOfNodes( )
     
 END FUNCTION GetNumberOfInternalNodes_SegmentMesh
!
!
!
 SUBROUTINE SetPositions_SegmentMesh( mySegmentMesh, iEl, x )
 ! S/R SetPositions
 !
 ! =============================================================================================== !
 ! DECLARATIONS
  IMPLICIT NONE
  CLASS(SegmentMesh), INTENT(inout) :: mySegmentMesh
  INTEGER, INTENT(in)               :: iEl
  REAL(prec), INTENT(in)            :: x(0:mySegmentMesh % elements(iEl) % nS)
  
     CALL mySegmentMesh % elements(iEl) % SetPositions( x )
     
 END SUBROUTINE SetPositions_SegmentMesh
!
!
!
 FUNCTION GetPositions_SegmentMesh( mySegmentMesh, iEl ) RESULT( x )
 ! FUNCTION GetPositions
 !
 ! =============================================================================================== !
 ! DECLARATIONS
  IMPLICIT NONE
  CLASS(SegmentMesh) :: mySegmentMesh
  INTEGER            :: iEl
  REAL(prec)         :: x(0:mySegmentMesh % elements(iEl) % nS)
  
     x = mySegmentMesh % elements(iEl) % GetPositions( )
     
 END FUNCTION GetPositions_SegmentMesh
!
!
!
 SUBROUTINE SetPositionAtNode_SegmentMesh( mySegmentMesh, iEl, i, x )
 ! S/R SetPositionAtNode
 !
 ! =============================================================================================== !
 ! DECLARATIONS
  IMPLICIT NONE
  CLASS(SegmentMesh), INTENT(inout) :: mySegmentMesh
  INTEGER, INTENT(in)               :: iEl
  REAL(prec), INTENT(in)            :: x
  INTEGER, INTENT(in)               :: i
  
     CALL mySegmentMesh % elements(iEl) % SetPositionAtNode( i, x )
     
 END SUBROUTINE SetPositionAtNode_SegmentMesh
!
!
!
 FUNCTION GetPositionAtNode_SegmentMesh( mySegmentMesh, iEl, i ) RESULT( x )
 ! FUNCTION GetPositionAtNode
 !
 ! =============================================================================================== !
 ! DECLARATIONS
  IMPLICIT NONE
  CLASS(SegmentMesh) :: mySegmentMesh
  INTEGER            :: iEl
  REAL(prec)         :: x
  INTEGER            :: i
  
     x = mySegmentMesh % elements(iEl) % GetPositionAtNode( i )
     
 END FUNCTION GetPositionAtNode_SegmentMesh
!
!
!
 SUBROUTINE SetJacobian_SegmentMesh( mySegmentMesh, iEl, J )
 ! S/R SetJacobian
 !
 ! =============================================================================================== !
 ! DECLARATIONS
  IMPLICIT NONE
  CLASS(SegmentMesh), INTENT(inout) :: mySegmentMesh
  INTEGER, INTENT(in)               :: iEl
  REAL(prec), INTENT(in)            :: J(0:mySegmentMesh % elements(iEl) % nS)
  
     CALL mySegmentMesh % elements(iEl) % SetJacobian( J )
     
 END SUBROUTINE SetJacobian_SegmentMesh
!
!
!
 FUNCTION GetJacobian_SegmentMesh( mySegmentMesh, iEl ) RESULT( J )
 ! FUNCTION GetJacobian
 !
 ! =============================================================================================== !
 ! DECLARATIONS
  IMPLICIT NONE
  CLASS(SegmentMesh) :: mySegmentMesh
  INTEGER            :: iEl
  REAL(prec)         :: J(0:mySegmentMesh % elements(iEl) % nS)
  
     J = mySegmentMesh % elements(iEl) % GetJacobian( )
     
 END FUNCTION GetJacobian_SegmentMesh
!
!
!
 SUBROUTINE SetJacobianAtNode_SegmentMesh( mySegmentMesh, iEl, iS, J )
 ! S/R SetJacobianAtNode
 !
 ! =============================================================================================== !
 ! DECLARATIONS
  IMPLICIT NONE
  CLASS(SegmentMesh), INTENT(inout) :: mySegmentMesh
  INTEGER, INTENT(in)               :: iEl
  REAL(prec), INTENT(in)            :: J
  INTEGER, INTENT(in)               :: iS
  
     CALL mySegmentMesh % elements(iEl) % SetJacobianAtNode( iS, J )
     
 END SUBROUTINE SetJacobianAtNode_SegmentMesh
!
!
!
 FUNCTION GetJacobianAtNode_SegmentMesh( mySegmentMesh, iEl, iS ) RESULT( J )
 ! FUNCTION GetJacobianAtNode 
 !
 ! =============================================================================================== !
 ! DECLARATIONS
  IMPLICIT NONE
  CLASS(SegmentMesh) :: mySegmentMesh
  INTEGER            :: iEl
  REAL(prec)         :: J
  INTEGER            :: iS
  
     J = mySegmentMesh % elements(iEl) % GetJacobianAtNode( iS )
     
 END FUNCTION GetJacobianAtNode_SegmentMesh
!
!
!
 SUBROUTINE SetBoundaryLocation_SegmentMesh( mySegmentMesh, iEl, iBound, x )
 ! S/R SetBoundaryLocation
 !
 ! =============================================================================================== !
 ! DECLARATIONS
  IMPLICIT NONE
  CLASS(SegmentMesh), INTENT(inout) :: mySegmentMesh
  INTEGER, INTENT(in)               :: iEl
  INTEGER, INTENT(in)               :: iBound
  REAL(prec), INTENT(in)            :: x
  
     CALL mySegmentMesh % elements(iEl) % SetBoundaryLocation( iBound, x )
     
 END SUBROUTINE SetBoundaryLocation_SegmentMesh
!
!
!
 FUNCTION GetBoundaryLocation_SegmentMesh( mySegmentMesh, iEl, iBound ) RESULT( x )
 ! FUNCTION GetBoundaryLocation
 !
 ! =============================================================================================== !
 ! DECLARATIONS
  IMPLICIT NONE
  CLASS(SegmentMesh) :: mySegmentMesh
  INTEGER            :: iEl
  INTEGER            :: iBound
  REAL(prec)         :: x
  
     x = mySegmentMesh % elements(iEl) % GetBoundaryLocation( iBound )
     
 END FUNCTION GetBoundaryLocation_SegmentMesh
!
!
!==================================================================================================!
!--------------------------------- Type Specific Routines -----------------------------------------!
!==================================================================================================!
!
!
 SUBROUTINE ScaleTheMesh_SegmentMesh( mySegmentMesh, xScale )
 ! S/R ScaleTheMesh
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS( SegmentMesh ), INTENT(inout) :: mySegmentMesh
   REAL(prec), INTENT(in)              :: xScale
   ! LOCAL
   INTEGER :: nElems, iEl
   
      nElems = mySegmentMesh % GetNumberOfElements( )

      DO iEl = 1, nElems
         CALL mySegmentMesh % elements(iEl) % ScaleGeometry( xScale )
      ENDDO

 END SUBROUTINE ScaleTheMesh_SegmentMesh
!
!
!
 SUBROUTINE LoadDefaultMesh_SegmentMesh( mySegmentMesh, interp, nXelem )
 ! S/R LoadDefaultMesh
 !  
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS( SegmentMesh ), INTENT(inout)  :: mySegmentMesh
   TYPE( Lagrange_1D ), INTENT(in)      :: interp
   INTEGER, INTENT(in)                  :: nXelem
   ! Local
   REAL(prec) :: dxElem
   REAL(prec) :: x1, x2
   INTEGER :: iEl, nS

      dxElem = ONE/real(nXElem,prec)
      CALL interp % GetNumberOfNodes( nS )

      CALL mySegmentMesh % Build( nXElem, nS ) 
      
      DO iEl = 1, nXElem
         x1 = real(iEl-1,prec)*dxElem
         x2 = x1 + dxElem
         CALL mySegmentMesh % elements(iEl) % Build( interp, x1, x2 )
      ENDDO
   
 END SUBROUTINE LoadDefaultMesh_SegmentMesh
!
!
!==================================================================================================!
!--------------------------------- Mesh File I/O Routines -----------------------------------------!
!==================================================================================================!
!
!
 SUBROUTINE WriteTecplot_SegmentMesh( mySegmentMesh, filename )
 ! S/R WriteTecplot
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS(SegmentMesh), INTENT(inout)     :: mySegmentMesh
   CHARACTER(*), INTENT(in), OPTIONAL    :: filename  
   ! Local
   INTEGER :: iS, nS, iEl, fUnit
   REAL(prec) :: x, J

      nS = mySegmentMesh % elements(1) % GetNumberOfNodes( )

      IF( PRESENT(filename) )THEN
         OPEN( UNIT=NEWUNIT(fUnit), &
               FILE= TRIM(filename)//'.curve', &
               FORM='formatted')
      ELSE
         OPEN( UNIT=NEWUNIT(fUnit), &
               FILE= 'mesh.curve', &
               FORM='formatted')
      ENDIF
      
      WRITE(fUnit,*) '#J'
   
      DO iEl = 1, mySegmentMesh % nElems
         DO iS = 0,nS
            x = mySegmentMesh % GetPositionAtNode( iEl, iS )
            J = mySegmentMesh % GetJacobianAtNode( iEl, iS )
            WRITE(fUnit,*) x, J  
         ENDDO
      ENDDO
    
      CLOSE(UNIT=fUnit)

 END SUBROUTINE WriteTecplot_SegmentMesh

END MODULE SegmentMeshClass
