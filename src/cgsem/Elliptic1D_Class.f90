! Elliptic1D_Class.f90 ( new with v2.1 - 25 March 2016)
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
! o  (ver 1.0) April 2015
! o  (ver 2.1) March 2016
!
! ========================================= Logs ================================================= !
! yyyy-mm-dd  Joe  <joe@clay>
!    This module defines the data structure used for implementing the Continuous Galerkin
!    Spectral Element Method. 
!
!    The module  is set up so that we solve d/dx( F ) = Q
!
!    For the heat equation : F = k du/dx, Q = du/dt
!
!    With this "generic flux" (F), and source (Q), it should be a trivial task to solve a different
!    system with CGSEM through the modification of F and Q. Here, the module is set up to solve
!    the heat-equation using the trapezoidal rule for time integration.
!
!   
! 
! //////////////////////////////////////////////////////////////////////////////////////////////// !
MODULE Elliptic1D_Class

! src/common/
USE ModelPrecision
USE ConstantsDictionary
USE ModelFlags
! src/interp/
USE Lagrange_1D_Class
! src/geom/
USE MappedGeometryClass_1D
USE SegmentMeshClass
!

IMPLICIT NONE

   TYPE Elliptic1D
      INTEGER                  :: nS, nElem
      TYPE( NodalStorage_1D )  :: cgStorage
      TYPE( SegmentMesh )      :: mesh
      REAL(prec), ALLOCATABLE  :: sol(:,:)
      REAL(prec), ALLOCATABLE  :: source(:,:)
      REAL(prec)               :: boundaryFlux(1:2)
      TYPE( Elliptic1DParams ) :: params
      
      CONTAINS
      
      
   END TYPE Elliptic1D

CONTAINS
!
!
!==================================================================================================!
!------------------------------- Manual Constructors/Destructors ----------------------------------!
!==================================================================================================!
!
!
 SUBROUTINE Build_Elliptic1D( thisElliptic )
 ! S/R Build
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS( Elliptic1D ), INTENT(inout) :: thisElliptic
   ! Local
   INTEGER :: nS, nEl
   
      CALL thisElliptic % params % Build( )
      nS  = thisElliptic % params % polyDeg
      nEl = thisElliptic % params % nElems
        
      thisElliptic % nElems = nEl
      thisElliptic % nS     = nS
      
      CALL thisElliptic % cgStorage % Build( nS, GAUSS_LOBATTO, CG )
      CALL thisElliptic % mesh % LoadDefaultMesh( thisElliptic % cgStorage % interp, &
                                                  nElems )
      
      ALLOCATE( thisElliptic % sol(0:nS,1:nEl), thisElliptic % source(0:nS,1:nEl) )
      thisElliptic % sol    = ZERO
      thisElliptic % source = ZERO
      thisElliptic % boundaryFlux = ZERO
      
 END SUBROUTINE Build_Elliptic1D
!
!
! 
 SUBROUTINE Trash_Elliptic1D( thisElliptic )
 ! S/R Trash
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS( Elliptic1D ), INTENT(inout) :: thisElliptic
   
      CALL thisElliptic % cgStorage % Trash( )
      CALL thisElliptic % mesh % Trash( )
      DEALLOCATE( thisElliptic % sol, thisElliptic % source )

 END SUBROUTINE Trash_Elliptic1D
!
!
!==================================================================================================!
!--------------------------------------- Accessors ------------------------------------------------!
!==================================================================================================!
!
!

!
!
!==================================================================================================!
!--------------------------------- Type Specific Routines -----------------------------------------!
!==================================================================================================!
!
!
 FUNCTION StiffnessWeight( x ) RESULT( K )
 ! FUNCTION StiffnessWeight
 !
 ! This function returns the "K" used in the computation of 
 ! the operator  " d/dx( K(x) d/dx ) ."
 ! =============================================================================================== !
 ! DECLARATION
   IMPLICIT NONE
   real(prec) :: x,K

      K = ONE

 END FUNCTION StiffnessWeight
!
!
!
 SUBROUTINE SetBoundaryConditions_Elliptic1D( thisElliptic, sL, sR )
 ! S/R SetBoundaryConditions
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS( Elliptic1D ), INTENT(inout) :: thisElliptic
   REAL(prec), INTENT(in)             :: sL, sR
   ! LOCAL
   INTEGER :: nS, nEl
   
      nS  = thisElliptic % nS
      nEl = thisElliptic % nEl
      
      ! Left-most boundary
      thisElliptic % sol(0,1) = sL

      ! Right-most boundary
      thisElliptic % sol(nS,nEl) = sR

 END SUBROUTINE SetBoundaryConditions_Elliptic1D
!
!
!
 SUBROUTINE Mask_Elliptic1D( thisElliptic, thisArray )
 ! S/R Mask
 !
 ! This subroutine masks out 'thisArray' along the shared element boundaries. 
 ! Arbitrarily, we choose to mask the element to the right (in 1-D) of
 ! the shared node.
 !
 ! For local operations, an "UNMASK" routine is needed - this is included 
 ! below.
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS( Elliptic1D ), INTENT(in) :: thisElliptic
   REAL(prec), INTENT(inout)       :: thisArray(0:thisElliptic % nS, 1:thisElliptic % nEl)
   ! LOCAL
   INTEGER :: iEl
 
      DO iEl = 1, (thisElliptic % nEl - 1)
         thisArray(0,iEl+1) = ZERO
      ENDDO


 END SUBROUTINE Mask_Elliptic1D
!
!
!
 SUBROUTINE Unmask_Elliptic1D( thisElliptic, thisArray )
 ! S/R Unmask
 !
 ! This subroutine unmasks 'thisArray' along the nodes shared by two elements.
 ! A copy of the solution retained in the "MASK" operation is used to unmask.
 !
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS( Elliptic1D ), INTENT(in) :: thisElliptic
   REAL(prec), INTENT(inout)       :: thisArray(0:thisElliptic % nS, 1:thisElliptic % nEl)
   ! LOCAL
   INTEGER :: iEl, nS, nEl
 
      nS  = thisElliptic % nS
      nEl = thisElliptic % nEl
      
      DO iEl = 1, (nEl - 1)
         thisArray(0, iEl+1) = thisArray(nS,iEl)
      ENDDO


 END SUBROUTINE Unmask_Elliptic1D
!
!
!
 SUBROUTINE GlobalSum_Elliptic1D( thisElliptic, thisArray )
 ! S/R GlobalSum
 !
 ! In the continuous galerkin spectral element method, the equation at the nodes shared 
 ! by two elements involve the sum of the equations implied by each element. This routine
 ! computes the sum of 'thisArray' at the shared nodes.
 !
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS( Elliptic1D ), INTENT(in) :: thisElliptic
   REAL(prec), INTENT(inout)       :: thisArray(0:thisElliptic % nS, 1:thisElliptic % nEl)
   ! LOCAL
   INTEGER    :: iEl, nEl, nS
   REAL(prec) :: temp
 
      nS  = thisElliptic % nS
      nEl = thisElliptic % nEl
      
      DO iEl = 1, (nEl - 1)
         temp = thisArray(0, iEl+1) + thisArray(nS,iEl)
         thisArray(0, iEl+1) = temp
         thisArray(nS, iEl)  = temp
      ENDDO


 END SUBROUTINE GlobalSum_Elliptic1D
 
END MODULE Elliptic1D_Class
