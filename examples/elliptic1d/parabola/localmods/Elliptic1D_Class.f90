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
! src/nodal/
USE NodalStorage_1D_Class
! src/geom/
USE MappedGeometryClass_1D
USE SegmentMeshClass
!src/iterativesolve/
USE IterativeSolvers_Class
!src/cgsem/
USE Elliptic1DParams_Class


IMPLICIT NONE

   TYPE, EXTENDS(ConjugateGradient) :: Elliptic1D
      INTEGER                  :: nS, nElems
      TYPE( NodalStorage_1D )  :: cgStorage
      TYPE( SegmentMesh )      :: mesh
      REAL(prec), ALLOCATABLE  :: cgsol(:,:)
      REAL(prec), ALLOCATABLE  :: source(:,:)
      REAL(prec), ALLOCATABLE  :: d2Mat(:,:,:) ! matrix for d/dx( k(x) d/dx ) operator
      REAL(prec)               :: boundaryFlux(1:2)
      REAL(prec)               :: dirichletMask(1:2)
      TYPE( Elliptic1DParams ) :: params
      
      CONTAINS
      
      PROCEDURE :: Build => Build_Elliptic1D
      PROCEDURE :: Trash => Trash_Elliptic1D

      PROCEDURE :: FillSource          => FillSource_Elliptic1D
      PROCEDURE :: DirichletConditions => DirichletConditions_Elliptic1D
     
      PROCEDURE :: Mask           => Mask_Elliptic1D
      PROCEDURE :: UnMask         => UnMask_Elliptic1D
      PROCEDURE :: UnMaskSolution => UnMaskSolution_Elliptic1D
      PROCEDURE :: GlobalSum      => GlobalSum_Elliptic1D
      PROCEDURE :: ApplyOperator  => ApplyOperator_Elliptic1D

      PROCEDURE :: MatrixAction => MatrixAction_Elliptic1D 
      PROCEDURE :: Residual     => Residual_Elliptic1D
      PROCEDURE :: CopyTo       => CopyTo_Elliptic1D
      PROCEDURE :: CopyFrom     => CopyFrom_Elliptic1D

      PROCEDURE :: WriteTecplot => WriteTecplot_Elliptic1D

   END TYPE Elliptic1D

CONTAINS
!
!
!==================================================================================================!
!------------------------------- Manual Constructors/Destructors ----------------------------------!
!==================================================================================================!
!
!
 SUBROUTINE Build_Elliptic1D( this )
 ! S/R Build
 !
 ! 
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS( Elliptic1D ), INTENT(inout) :: this
   ! Local
   INTEGER :: row, col, m, iEl, nS, nEl, maxIter, nDOF
   REAL(prec) :: x, sK, J, lmc, lmr
   REAL(prec), ALLOCATABLE :: w(:), dMat(:,:)
   
      CALL this % params % Build( )
      nS  = this % params % polyDeg
      nEl = this % params % nElems
      maxIter = this % params % maximumIterates
      nDOF = (nS+1)*nEl      

      this % nElems = nEl
      this % nS     = nS
      this % maxIters  = maxIter
      this % nDOF      = nDOF
      this % tolerance = this % params % tolerance

      CALL this % cgStorage % Build( nS, GAUSS_LOBATTO, CG )
      CALL this % mesh % LoadDefaultMesh( this % cgStorage % interp, nEl )

      ALLOCATE( this % sol(1:nDOF), this % resi(0:maxIter) )
      ALLOCATE( this % cgsol(0:nS,1:nEl), this % source(0:nS,1:nEl) )
      this % sol           = ZERO
      this % resi          = ZERO 
      this % cgsol         = ZERO
      this % source        = ZERO
      this % boundaryFlux  = ZERO
      this % dirichletMask = ZERO

      CALL this % FillSource( )
      ! In 1-D, for a fixed mesh and stiffness coefficient, the second-derivative matrix
      ! for each element can be constructed and stored once. This way, to compute the elliptic
      ! matrix action, only one matrix-vector multiply is needed for each element. The alternative
      ! is to apply two matrix-vector operations with the 1st derivative matrices.
      !
      ! Here, we opt to build and store the second-derivative operator for each element.

      ALLOCATE( w(0:nS), dMat(0:nS,0:nS) )
      ALLOCATE( this % d2Mat(0:nS,0:nS,1:nEl) )
      CALL this % cgStorage % GetQuadratureWeights( w )
      CALL this % cgStorage % GetDerivativeMatrix( dMat )

      DO iEl = 1, nEl  
         DO col = 0, nS
            DO row = 0, nS
               this % d2Mat(row,col,iEl) = ZERO ! Initialize the sum to zero
               DO m = 0, nS 

                  x = this % mesh % GetPositionAtNode( iEl, m ) 
                  sK = StiffnessWeight( x ) ! calculate the stiffness weight, "K" in d/dx( K(x) d/dx )
                  J = this % mesh % GetJacobianAtNode( iEl, m )

                  lmr = dMat(m,row) 
                  lmc = dMat(m,col) 

                  this % d2Mat(row,col,iEl) = this % d2Mat(row,col,iEl) - (sK*w(m)/J)*lmr*lmc

               ENDDO 
            ENDDO 
         ENDDO 
      ENDDO 

      DEALLOCATE( w, dMat )

 END SUBROUTINE Build_Elliptic1D
!
!
! 
 SUBROUTINE Trash_Elliptic1D( this )
 ! S/R Trash
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS( Elliptic1D ), INTENT(inout) :: this
   
      CALL this % cgStorage % Trash( )
      CALL this % mesh % Trash( )
      DEALLOCATE( this % sol, this % source, this % d2Mat )

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
 SUBROUTINE FillSource_Elliptic1D( thisElliptic )
 !
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS( Elliptic1D ), INTENT(inout) :: thisElliptic
   ! LOCAL
   INTEGER    :: iS, iEl, nS, nEl
   REAL(prec) :: x

      nS  = thisElliptic % nS
      nEl = thisElliptic % nElems

      DO iEl = 1, nEl
         DO iS = 0, nS
            x = thisElliptic % mesh % GetPositionAtNode(iEl,iS)
            thisElliptic % source(iS,iEl) = ONE
         ENDDO
      ENDDO
 
 END SUBROUTINE FillSource_Elliptic1D
!
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
   REAL(prec) :: x, K

      K = ONE

 END FUNCTION StiffnessWeight
!
!
!
 SUBROUTINE DirichletConditions_Elliptic1D( thisElliptic )
 ! S/R DirichletConditions
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS( Elliptic1D ), INTENT(inout) :: thisElliptic
   ! LOCAL
   INTEGER :: nS, nEl
   
      nS  = thisElliptic % nS
      nEl = thisElliptic % nElems
      
      IF( thisElliptic % dirichletMask(1) == ZERO )THEN
         ! Left-most boundary
         thisElliptic % cgsol(0,1) = ZERO
      ENDIF

      IF( thisElliptic % dirichletMask(2) == ZERO )THEN
         ! Right-most boundary
         thisElliptic % cgsol(nS,nEl) = ZERO
      ENDIF

 END SUBROUTINE DirichletConditions_Elliptic1D
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
   REAL(prec), INTENT(inout)       :: thisArray(0:thisElliptic % nS, 1:thisElliptic % nElems)
   ! LOCAL
   INTEGER :: iEl, nEl, nS
 
      nEl = thisElliptic % nElems
      nS  = thisElliptic % nS

      DO iEl = 1, (nEl - 1)
         thisArray(0,iEl+1) = ZERO
      ENDDO

      ! If dirichlet boundary conditions are applied at either boundary, the masking routine 
      ! (assuming it is taking in a residual for an iterative solver) needs to force the 
      ! incoming array to zero at those boundary points. Here, this is achieved by using the
      ! "dirichletMask" attribute.
      thisArray(0,1)    = thisArray(0,1)*thisElliptic % dirichletMask(1)
      thisArray(nS,nEl) = thisArray(nS,nEl)*thisElliptic % dirichletMask(2) 

 END SUBROUTINE Mask_Elliptic1D
!
!
!
 SUBROUTINE UnMask_Elliptic1D( thisElliptic, thisArray )
 ! S/R UnMask
 !
 ! This subroutine unmasks 'thisArray' along the nodes shared by two elements.
 ! A copy of the solution retained in the "MASK" operation is used to unmask.
 !
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS( Elliptic1D ), INTENT(in) :: thisElliptic
   REAL(prec), INTENT(inout)       :: thisArray(0:thisElliptic % nS, 1:thisElliptic % nElems)
   ! LOCAL
   INTEGER :: iEl, nS, nEl
 
      nS  = thisElliptic % nS
      nEl = thisElliptic % nElems
      
      DO iEl = 1, (nEl - 1)
         thisArray(0, iEl+1) = thisArray(nS,iEl)
      ENDDO


 END SUBROUTINE UnMask_Elliptic1D
!
!
!
 SUBROUTINE UnMaskSolution_Elliptic1D( thisElliptic )
 ! S/R UnMaskSolution
 !
 ! This subroutine unmasks the solution attribute and applies the Dirichlet boundary conditions
 ! if appropriate.
 !
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS( Elliptic1D ), INTENT(inout) :: thisElliptic
      
      CALL thisElliptic % UnMask( thisElliptic % cgsol )
      CALL thisElliptic % DirichletConditions( )


 END SUBROUTINE UnMaskSolution_Elliptic1D
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
   REAL(prec), INTENT(inout)       :: thisArray(0:thisElliptic % nS, 1:thisElliptic % nElems)
   ! LOCAL
   INTEGER    :: iEl, nEl, nS
   REAL(prec) :: temp
 
      nS  = thisElliptic % nS
      nEl = thisElliptic % nElems
      
      DO iEl = 1, (nEl - 1)
         temp = thisArray(0, iEl+1) + thisArray(nS,iEl)
         thisArray(0, iEl+1) = temp
         thisArray(nS, iEl)  = temp
      ENDDO

 END SUBROUTINE GlobalSum_Elliptic1D
!
!
!
 SUBROUTINE ApplyOperator_Elliptic1D( thisElliptic, u, Lu )
 ! S/R ApplyOperator
 !
 ! This subroutine loops over all of the elements an applies the internal matrix operations
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS( Elliptic1D ), INTENT(in) :: thisElliptic
   REAL(prec), INTENT(in)          :: u(0:thisElliptic % nS, 1:thisElliptic % nElems)
   REAL(prec), INTENT(out)         :: Lu(0:thisElliptic % nS, 1:thisElliptic % nElems) 
   ! LOCAL
   INTEGER :: iEl, nS, nEl
   REAL(prec) :: locSol(0:thisElliptic % nS), locDmat(0:thisElliptic % nS, 0:thisElliptic % nS)

      nS = thisElliptic % nS
      nEl = thisElliptic % nElems

      DO iEl = 1, nEl ! loop over the elements
         locSol  = u(0:nS,iEl)
         locDmat = thisElliptic % d2Mat(0:nS,0:nS,iEl)
         Lu(0:nS,iEl) = MATMUL( locDmat, locSol )
      ENDDO ! iEl, loop over the elements
        
      
 END SUBROUTINE ApplyOperator_Elliptic1D
!
!
! 
 FUNCTION MatrixAction_Elliptic1D( this, s )
 !
 ! 
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS( Elliptic1D ) :: this
   REAL(prec)          :: s(1:this % nDOF)
   REAL(prec)          :: MatrixAction_Elliptic1D(1:this % nDOF)
   ! Local
   REAL(prec) :: Au(0:this % nS, 1:this % nElems)
   integer    :: iEl, iS, nEl, nS, k

      nEl = this % nElems
      nS  = this % nS

      CALL this % CopyFrom( s )

      CALL this % UnMaskSolution( )

      CALL this % ApplyOperator( this % cgsol, Au )

      CALL this % GlobalSum( Au )

      CALL this % Mask( Au )
 
      DO iEl = 1, nEl
         DO iS = 0, nS
            k = iS + 1 + (nS+1)*(iEl-1)
            MatrixAction_Elliptic1D(k) = Au(iS,iEl)
         ENDDO
      ENDDO
 

 END FUNCTION MatrixAction_Elliptic1D
!
!
!
 FUNCTION Residual_Elliptic1D( this, s )
 ! Residual
 ! 
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS( Elliptic1D ) :: this
   REAL(prec)          :: s(1:this % nDOF)
   REAL(prec)          :: Residual_Elliptic1D(1:this % nDOF)
   ! Local
   REAL(prec) :: Au(1:this % nDOF)
   REAL(prec) :: r(0:this % nS, 1:this % nElems)
   REAL(prec) :: w(0:this % nS), J
   integer    :: iEl, iS, nEl, nS, k

      nEl = this % nElems
      nS  = this % nS
      CALL this % cgStorage % GetQuadratureWeights( w )

      Au = this % MatrixAction( s )

      DO iEl = 1, nEl
         DO iS = 0, nS
            k = iS + 1 + (nS+1)*(iEl-1)
            J = this % mesh % GetJacobianAtNode( iEl, iS )
            r(iS,iEl) = this % source(iS,iEl)*J*w(iS) - Au(k)
         ENDDO
      ENDDO

      ! If the dirichletMask is set to zero at either boundary, then the dirichlet boundary
      ! conditions are enforced. If it is set to one at either boundary, then a boundary
      ! flux should be added to the matrix action for the associated degree of freedom. To
      ! implement such a Neumann boundary condition, the boundary flux is added to the 
      ! source for the first and last nodes if the dirichletMask is set to one.
      r(0,1)    = r(0,1) + this % boundaryFlux(1)*this % dirichletMask(1)
      r(nS,nEl) = r(nS,nEl) - this % boundaryFlux(2)*this % dirichletMask(2)

      CALL this % Mask( r )

      DO iEl = 1, nEl
         DO iS = 0, nS
            k = iS + 1 + (nS+1)*(iEl-1)
            Residual_Elliptic1D(k) = r(iS,iEl)
         ENDDO
      ENDDO

 END FUNCTION Residual_Elliptic1D
!
!
!
 SUBROUTINE CopyFrom_Elliptic1D( this, sol )
 ! S/R CopyFrom
 !
 !  Copies FROM sol to the Elliptic1D class
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS( Elliptic1D ), INTENT(inout) :: this
   REAL(prec), INTENT(in)             :: sol(1:this % nDOF)
   ! Local
   INTEGER :: iS, iEl, k, nS, nEl
   
      nEl = this % nElems
      nS  = this % nS

      DO iEl = 1, nEl
         DO iS = 0, nS
            k = iS + 1 + (nS+1)*(iEl-1)
            this % cgsol(iS,iEl) = sol(k)
         ENDDO
      ENDDO
  
 END SUBROUTINE CopyFrom_Elliptic1D 
!
!
!
 SUBROUTINE CopyTo_Elliptic1D( this, sol )
 ! S/R CopyTo
 !
 !   Copies TO sol from the Elliptic1D class
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS( Elliptic1D ), INTENT(in) :: this
   REAL(prec), INTENT(out)         :: sol(1:this % nDOF)
   ! Local
   INTEGER :: iS, iEl, k, nS, nEl
   
      nEl = this % nElems
      nS  = this % nS

      DO iEl = 1, nEl
         DO iS = 0, nS
            k = iS + 1 + (nS+1)*(iEl-1)
            sol(k) = this % cgsol(iS,iEl) 
         ENDDO
      ENDDO
  
 END SUBROUTINE CopyTo_Elliptic1D 
!
!
!==================================================================================================!
!-------------------------------------- File I/O Routines -----------------------------------------!
!==================================================================================================!
!
!
 SUBROUTINE WriteTecplot_Elliptic1D( thisElliptic )
 ! S/R WriteTecplot
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS( Elliptic1D ), INTENT(inout) :: thisElliptic
   ! Local
   INTEGER    :: iS, iEl, nS, nEl, fUnit
   REAL(prec) :: x

      nEl = thisElliptic % nElems
      nS  = thisElliptic % nS

      CALL thisElliptic % UnMaskSolution( )

      OPEN( UNIT = NewUnit(fUnit), &
            FILE = 'elliptic.curve', &
            FORM = 'FORMATTED' )

      WRITE(fUnit,*)'#solution'
      DO iEl = 1, nEl
         DO iS = 0, nS
            x = thisElliptic % mesh % GetPositionAtNode(iEl,iS)
            WRITE(fUnit,*) x, thisElliptic % cgsol(iS,iEl)
         ENDDO
      ENDDO

      WRITE(fUnit,*)'#source'
      DO iEl = 1, nEl
         DO iS = 0, nS
            x = thisElliptic % mesh % GetPositionAtNode(iEl,iS)
            WRITE(fUnit,*) x, thisElliptic % source(iS,iEl)
         ENDDO
      ENDDO

      CLOSE(fUnit)
   
 END SUBROUTINE WriteTecplot_Elliptic1D

END MODULE Elliptic1D_Class
