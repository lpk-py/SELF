! MappedGeometryClass_1D.f90 ( new with v2.1 - 25 March 2016)
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
! o  (ver 2.1) March 2016
!
! ========================================= Logs ================================================= !
! yyyy-mm-dd  Joe  <joe@clay>
!
!   
! 
! //////////////////////////////////////////////////////////////////////////////////////////////// !


MODULE MappedGeometryClass_1D

! src/common/
USE ModelPrecision
USE ConstantsDictionary
USE ModelFlags
USE CommonRoutines
! src/interp/
USE Lagrange_1D_Class



IMPLICIT NONE



   TYPE MappedGeometry_1D
      INTEGER                             :: nS
      REAL(prec), ALLOCATABLE, PRIVATE    :: xBound(:)
      REAL(prec), ALLOCATABLE, PRIVATE    :: x(:)
      REAL(prec), ALLOCATABLE, PRIVATE    :: J(:)    

      CONTAINS

      PROCEDURE :: Build => Build_MappedGeometry_1D
      PROCEDURE :: Trash => Trash_MappedGeometry_1D
      
      PROCEDURE :: SetNumberOfNodes => SetNumberOfNodes_MappedGeometry_1D
      PROCEDURE :: GetNumberOfNodes => GetNumberOfNodes_MappedGeometry_1D
      PROCEDURE :: SetPositions => SetPositions_MappedGeometry_1D
      PROCEDURE :: GetPositions => GetPositions_MappedGeometry_1D
      PROCEDURE :: SetPositionAtNode => SetPositionAtNode_MappedGeometry_1D
      PROCEDURE :: GetPositionAtNode => GetPositionAtNode_MappedGeometry_1D
      PROCEDURE :: SetJacobian => SetJacobian_MappedGeometry_1D
      PROCEDURE :: GetJacobian => GetJacobian_MappedGeometry_1D
      PROCEDURE :: SetJacobianAtNode => SetJacobianAtNode_MappedGeometry_1D
      PROCEDURE :: GetJacobianAtNode => GetJacobianAtNode_MappedGeometry_1D
      PROCEDURE :: SetBoundaryLocation => SetBoundaryLocation_MappedGeometry_1D
      PROCEDURE :: GetBoundaryLocation => GetBoundaryLocation_MappedGeometry_1D

      PROCEDURE :: GenerateMesh => GenerateMesh_MappedGeometry_1D
      PROCEDURE :: GenerateMetrics => GenerateMetrics_MappedGeometry_1D
      PROCEDURE :: ScaleGeometry => ScaleGeometry_MappedGeometry_1D
      PROCEDURE :: CalculateLocation => CalculateLocation_MappedGeometry_1D
      PROCEDURE :: CalculateMetrics => CalculateMetrics_MappedGeometry_1D
      PROCEDURE :: CalculateComputationalCoordinates => CalculateComputationalCoordinates_MappedGeometry_1D
      
      PROCEDURE :: WriteTecplot => WriteTecplot_MappedGeometry_1D
      
   END TYPE MappedGeometry_1D


 INTEGER, PRIVATE :: nDims = 1
 CONTAINS
!
!
!==================================================================================================!
!------------------------------- Manual Constructors/Destructors ----------------------------------!
!==================================================================================================!
!
!
 SUBROUTINE Build_MappedGeometry_1D( myGeom, interp, x1, x2 )
 ! S/R Build
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS(MappedGeometry_1D), INTENT(out) :: myGeom
   TYPE(Lagrange_1D), INTENT(in)         :: interp
   REAL(prec), INTENT(in)                :: x1, x2
   !LOCAL
   INTEGER :: nS
   

      CALL interp % GetNumberOfNodes( nS )
      CALL myGeom % SetNumberOfNodes( nS )
       
      ! Allocate space
      ALLOCATE( myGeom % J(0:nS) )
      ALLOCATE( myGeom % x(0:nS) )
      ALLOCATE( myGeom % xBound(1:2) )
     
      ! Generate the mesh locations, and the mapping metrics

      CALL myGeom % GenerateMesh( interp, x1, x2 )
      CALL myGeom % GenerateMetrics( interp )
 
 END SUBROUTINE Build_MappedGeometry_1D
!
!
!
 SUBROUTINE Trash_MappedGeometry_1D( myGeom )
 ! S/R Trash
 !
 ! =============================================================================================== !
 ! DECLARATIONS
  IMPLICIT NONE
  CLASS(MappedGeometry_1D), INTENT(inout)  :: myGeom

      DEALLOCATE( myGeom % J, myGeom % x )
      DEALLOCATE( myGeom % xBound )

 END SUBROUTINE Trash_MappedGeometry_1D
!
!
!==================================================================================================!
!--------------------------------------- Accessors ------------------------------------------------!
!==================================================================================================!
!
!
 SUBROUTINE SetNumberOfNodes_MappedGeometry_1D( myGeom, nS )
 ! S/R SetNumberOfNodes
 !
 ! =============================================================================================== !
 ! DECLARATIONS
  IMPLICIT NONE
  CLASS(MappedGeometry_1D), INTENT(inout) :: myGeom
  INTEGER, INTENT(in)                     :: nS
  
     myGeom % nS   = nS
     
 END SUBROUTINE SetNumberOfNodes_MappedGeometry_1D
!
!
!
 FUNCTION GetNumberOfNodes_MappedGeometry_1D( myGeom ) RESULT( nS )
 ! FUNCTION GetNumberOfNodes
 !
 ! =============================================================================================== !
 ! DECLARATIONS
  IMPLICIT NONE
  CLASS(MappedGeometry_1D) :: myGeom
  INTEGER                  :: nS
  
     nS = myGeom % nS
     
 END FUNCTION GetNumberOfNodes_MappedGeometry_1D
!
!
!
 SUBROUTINE SetPositions_MappedGeometry_1D( myGeom, x )
 ! S/R SetPositions
 !
 ! =============================================================================================== !
 ! DECLARATIONS
  IMPLICIT NONE
  CLASS(MappedGeometry_1D), INTENT(inout) :: myGeom
  REAL(prec), INTENT(in)                  :: x(0:myGeom % nS)
  
     myGeom % x = x
     
 END SUBROUTINE SetPositions_MappedGeometry_1D
!
!
!
 FUNCTION GetPositions_MappedGeometry_1D( myGeom ) RESULT( x )
 ! FUNCTION GetPositions
 !
 ! =============================================================================================== !
 ! DECLARATIONS
  IMPLICIT NONE
  CLASS(MappedGeometry_1D) :: myGeom
  REAL(prec)               :: x(0:myGeom % nS)
  
     x = myGeom % x
     
 END FUNCTION GetPositions_MappedGeometry_1D
!
!
!
 SUBROUTINE SetPositionAtNode_MappedGeometry_1D( myGeom, i, x )
 ! S/R SetPositionAtNode
 !
 ! =============================================================================================== !
 ! DECLARATIONS
  IMPLICIT NONE
  CLASS(MappedGeometry_1D), INTENT(inout) :: myGeom
  REAL(prec), INTENT(in)                  :: x
  INTEGER, INTENT(in)                     :: i
  
     myGeom % x(i) = x
     
 END SUBROUTINE SetPositionAtNode_MappedGeometry_1D
!
!
!
 FUNCTION GetPositionAtNode_MappedGeometry_1D( myGeom, i ) RESULT( x )
 ! FUNCTION GetPositionAtNode
 !
 ! =============================================================================================== !
 ! DECLARATIONS
  IMPLICIT NONE
  CLASS(MappedGeometry_1D) :: myGeom
  INTEGER                  :: i
  REAL(prec)               :: x
  
     x = myGeom % x(i)
     
 END FUNCTION GetPositionAtNode_MappedGeometry_1D
!
!
!
 SUBROUTINE SetJacobian_MappedGeometry_1D( myGeom, J )
 ! S/R SetJacobian
 !
 ! =============================================================================================== !
 ! DECLARATIONS
  IMPLICIT NONE
  CLASS(MappedGeometry_1D), INTENT(inout) :: myGeom
  REAL(prec), INTENT(in)                  :: J(0:myGeom % nS)
  
     myGeom % J = J
     
 END SUBROUTINE SetJacobian_MappedGeometry_1D
!
!
!
 FUNCTION GetJacobian_MappedGeometry_1D( myGeom ) RESULT( J )
 ! FUNCTION GetJacobian
 !
 ! =============================================================================================== !
 ! DECLARATIONS
  IMPLICIT NONE
  CLASS(MappedGeometry_1D) :: myGeom
  REAL(prec)               :: J(0:myGeom % nS)
  
     J = myGeom % J
     
 END FUNCTION GetJacobian_MappedGeometry_1D
!
!
!
 SUBROUTINE SetJacobianAtNode_MappedGeometry_1D( myGeom, iS, J )
 ! S/R SetJacobianAtNode
 !
 ! =============================================================================================== !
 ! DECLARATIONS
  IMPLICIT NONE
  CLASS(MappedGeometry_1D), INTENT(inout) :: myGeom
  REAL(prec), INTENT(in)                  :: J
  INTEGER, INTENT(in)                     :: iS
  
     myGeom % J(iS) = J
     
 END SUBROUTINE SetJacobianAtNode_MappedGeometry_1D
!
!
!
 FUNCTION GetJacobianAtNode_MappedGeometry_1D( myGeom, iS ) RESULT( J )
 ! FUNCTION GetJacobianAtNode
 !
 ! =============================================================================================== !
 ! DECLARATIONS
  IMPLICIT NONE
  CLASS(MappedGeometry_1D) :: myGeom
  INTEGER                  :: iS
  REAL(prec)               :: J
  
     J = myGeom % J(iS)
     
 END FUNCTION GetJacobianAtNode_MappedGeometry_1D
!
!
!
 SUBROUTINE SetBoundaryLocation_MappedGeometry_1D( myGeom, iBound, x )
 ! S/R SetBoundaryLocation
 !
 ! =============================================================================================== !
 ! DECLARATIONS
  IMPLICIT NONE
  CLASS(MappedGeometry_1D), INTENT(inout) :: myGeom
  INTEGER, INTENT(in)                     :: iBound
  REAL(prec), INTENT(in)                  :: x
  
     myGeom % xBound(iBound) = x
     
 END SUBROUTINE SetBoundaryLocation_MappedGeometry_1D
!
!
!
 FUNCTION GetBoundaryLocation_MappedGeometry_1D( myGeom, iBound ) RESULT( x )
 ! FUNCTION GetBoundaryLocation
 !
 ! =============================================================================================== !
 ! DECLARATIONS
  IMPLICIT NONE
  CLASS(MappedGeometry_1D) :: myGeom
  INTEGER                  :: iBound
  REAL(prec)               :: x
  
     x = myGeom % xBound(iBound)
     
 END FUNCTION GetBoundaryLocation_MappedGeometry_1D
!
!
!==================================================================================================!
!--------------------------------- Type Specific Routines -----------------------------------------!
!==================================================================================================!
!
!
 SUBROUTINE GenerateMesh_MappedGeometry_1D( myGeom, interp, x1, x2 )
 ! S/R GenerateMesh
 !
 ! =============================================================================================== !
 ! DECLARATIONS
  IMPLICIT NONE
  CLASS( MappedGeometry_1D ), INTENT(inout) :: myGeom
  TYPE( Lagrange_1D ), INTENT(in)           :: interp
  REAL(prec), INTENT(in)                    :: x1, x2
  ! Local
  INTEGER    :: iS, nS
  REAL(prec) :: s, x
   
      CALL interp % GetNumberOfNodes( nS )
      
      DO iS = 0, nS
         CALL interp % GetNode( iS, s ) 
         x = HALF*(x2-x1)*(s + ONE) + x1
         CALL myGeom % SetPositionAtNode( iS, x )
      ENDDO 
      
      CALL myGeom % SetBoundaryLocation( left,  x1 )
      CALL myGeom % SetBoundaryLocation( right, x2 )

 END SUBROUTINE GenerateMesh_MappedGeometry_1D
!
!
!
 SUBROUTINE GenerateMetrics_MappedGeometry_1D( myGeom, interp )
 ! S/R GenerateMetrics
 !
 ! =============================================================================================== !
 ! DECLARATIONS
  IMPLICIT NONE
  CLASS( MappedGeometry_1D ), INTENT(inout) :: myGeom
  TYPE( Lagrange_1D ), INTENT(in)           :: interp
  ! Local
  INTEGER    :: iS, nS
  REAL(prec) :: dxds, s
   
      CALL interp % GetNumberOfNodes( nS )
      
      DO iS = 0, nS
         CALL interp % GetNode( iS, s )
         dxds = interp % EvaluateDerivative( s, myGeom % x )
         CALL myGeom % SetJacobianAtNode( iS, dxds )
      ENDDO ! iS

 END SUBROUTINE GenerateMetrics_MappedGeometry_1D
!
!
! 
 FUNCTION CalculateLocation_MappedGeometry_1D( myGeom, interp, s ) RESULT( x )
 ! FUNCTION CalculateLocation
 !  Description :
 ! 
 !     This subroutine calculates the physical coordinates (x) given the computational coordinates
 !     (s).
 ! =============================================================================================== !
  IMPLICIT NONE
  CLASS( MappedGeometry_1D ) :: myGeom
  TYPE( Lagrange_1D )        :: interp
  REAL(prec)                 :: s
  REAL(prec)                 :: x
  
     x = interp % EvaluateInterpolant( s, myGeom % x )
  
 END FUNCTION CalculateLocation_MappedGeometry_1D
!
!
!
 FUNCTION CalculateMetrics_MappedGeometry_1D( myGeom, interp, s ) RESULT( J )
 ! FUNCTION CalculateMetrics
 !  Description :
 ! 
 ! =============================================================================================== !
  IMPLICIT NONE
  CLASS( MappedGeometry_1D ) :: myGeom
  TYPE( Lagrange_1D )        :: interp
  REAL(prec)                 :: s
  REAL(prec)                 :: J

     J = interp % EvaluateDerivative( s, myGeom % x )

 END FUNCTION CalculateMetrics_MappedGeometry_1D
!
!
!
 SUBROUTINE CalculateComputationalCoordinates_MappedGeometry_1D( myGeom, interp, x, s, success )
 ! S/R CalculateComputationalCoordinates
 !  Description :
 ! 
 !     Given the physical coordinates (x*), the computational coordinates are calculated using 
 !     Newton's method for root finding to solve
 !
 !     x* = x(s)
 !
 !     for s.
 ! =============================================================================================== !
  IMPLICIT NONE
  CLASS( MappedGeometry_1D ), INTENT(in) :: myGeom
  TYPE( Lagrange_1D ), INTENT(in)        :: interp
  REAL(prec), INTENT(in)                 :: x
  REAL(prec), INTENT(out)                :: s
  LOGICAL, INTENT(out), OPTIONAL         :: success
  ! LOCAL
  REAL(prec) :: dr, ds, J
  REAL(prec) :: thisX, thisS, resi
  INTEGER    :: i 

     thisS = ZERO ! Initial guess is at the origin
     
     IF( PRESENT(success) )THEN
        success = .FALSE.
     ENDIF
     
     DO i = 1, newtonMax
     
        ! Calculate the physical coordinate associated with the computational coordinate guess
        thisX = myGeom % CalculateLocation( interp, thisS )
     
        ! Calculate the residual
        dr = x - thisX 
        resi = abs(dr)
     
        IF( resi < newtonTolerance )THEN
           s = thisS
           IF( PRESENT(success) )THEN
              success = .TRUE.
           ENDIF
           RETURN
        ENDIF
        
        J = myGeom % CalculateMetrics( interp, thisS ) ! Calculate the Jacobian
        
        ds = dr/J ! calculate the correction in the computational coordinate
        thisS = thisS + ds
     
     ENDDO
     
     ! Calculate the residual
     dr = x - thisX 
     resi = ABS(dr)
     IF( resi < newtonTolerance )THEN
        s = thisS
        IF( PRESENT(success) )THEN
           success = .TRUE.
        ENDIF
        RETURN
     ELSE
        s = fillValue
        PRINT*,'Module MappedGeometryClass_1D.f90 : S/R CalculateComputationalCoordinates :'
        PRINT*,'Search for coordinates failed. Final residual norm :', resi
        RETURN
    ENDIF
      
     
 END SUBROUTINE CalculateComputationalCoordinates_MappedGeometry_1D
!
!
!
 SUBROUTINE ScaleGeometry_MappedGeometry_1D( myGeom, xScale )
 !
 !
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS( MappedGeometry_1D ), INTENT(inout) :: myGeom
   REAL(prec), INTENT(in)                    :: xScale
   
         myGeom % x = xScale*( myGeom % x )
         myGeom % xBound = xScale*( myGeom % xBound )
         myGeom % J = xScale*( myGeom % J )
         
 END SUBROUTINE ScaleGeometry_MappedGeometry_1D
!
!
!==================================================================================================!
!------------------------------------ File I/O Routines -------------------------------------------!
!==================================================================================================!
!
!
 SUBROUTINE WriteTecplot_MappedGeometry_1D( myGeom, filename )
 ! S/R WriteTecplot
 !  Description :
 ! 
 ! =============================================================================================== !
   IMPLICIT NONE
   CLASS( MappedGeometry_1D ), INTENT(in) :: myGeom
   CHARACTER(*), INTENT(in), OPTIONAL     :: filename  
   ! Local
   INTEGER :: iS, nS, fUnit
  
      nS = myGeom % GetNumberOfNodes( )
    
      IF( PRESENT(filename) )THEN
         OPEN( UNIT   = NEWUNIT(fUnit), &
               FILE   = TRIM(filename)//'.curve', &
               FORM   = 'formatted', & 
               STATUS = 'REPLACE' )
      ELSE
         OPEN( UNIT   = NEWUNIT(fUnit), &
               FILE   = 'LocalGeometry.curve', &
               FORM   = 'formatted', & 
               STATUS = 'REPLACE' )
      ENDIF

      WRITE(fUnit,*) '#J '
    
      DO iS = 0, nS
         WRITE(fUnit,*)  myGeom % x(iS), myGeom % J(iS)
      ENDDO

    CLOSE(UNIT=fUnit)

    RETURN

 END SUBROUTINE WriteTecplot_MappedGeometry_1D
 
 
END MODULE MappedGeometryClass_1D
