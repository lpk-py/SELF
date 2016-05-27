! BarotropicShelfWaves_Class.f90 ( new with v2.1 - 25 March 2016)
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
!    With this "generic flux" (F), and rhsVar (Q), it should be a trivial task to solve a different
!    system with CGSEM through the modification of F and Q. Here, the module is set up to solve
!    the heat-equation using the trapezoidal rule for time integration.
!
!   
! 
! //////////////////////////////////////////////////////////////////////////////////////////////// !

MODULE BarotropicShelfWaves_Class

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

!src/cgsem/
USE BarotropicShelfWavesParams_Class


IMPLICIT NONE

   TYPE BarotropicShelfWaves
      INTEGER                  :: nS, nElems
      INTEGER                  :: nEigenPairs, maxIters, nMaxToSolve
      REAL(prec)               :: tolerance
      TYPE( NodalStorage_1D )  :: cgStorage
      TYPE( SegmentMesh )      :: mesh
      REAL(prec), ALLOCATABLE  :: eigenvalues(:)
      REAL(prec), ALLOCATABLE  :: growthRates(:)
      REAL(prec), ALLOCATABLE  :: eigenfunctions(:,:,:)
      REAL(prec), ALLOCATABLE  :: cgSol(:,:)
      REAL(prec), ALLOCATABLE  :: rhsVar(:,:)
      REAL(prec), ALLOCATABLE  :: d2Mat(:,:,:) ! matrix for d/dx( k(x) d/dx ) operator within each element
      REAL(prec), ALLOCATABLE  :: wFunc(:,:)   ! Sturm-Liouville weight function
      REAL(prec), ALLOCATABLE  :: H(:,:)       ! The shelf profile
      REAL(prec), ALLOCATABLE  :: G(:,:)       ! The forcing function due to steepening
      REAL(prec)               :: boundaryFlux(1:2)
      REAL(prec)               :: dirichletMask(1:2)
      TYPE( BarotropicShelfWavesParams ) :: params
      
      CONTAINS
      
      PROCEDURE :: Build => Build_BarotropicShelfWaves
      PROCEDURE :: Trash => Trash_BarotropicShelfWaves

      PROCEDURE :: DirichletConditions => DirichletConditions_BarotropicShelfWaves
     
      PROCEDURE :: Mask           => Mask_BarotropicShelfWaves
      PROCEDURE :: UnMask         => UnMask_BarotropicShelfWaves
      PROCEDURE :: UnMaskSolution => UnMaskSolution_BarotropicShelfWaves
      PROCEDURE :: GlobalSum      => GlobalSum_BarotropicShelfWaves
      PROCEDURE :: ApplyOperator  => ApplyOperator_BarotropicShelfWaves

      PROCEDURE :: StiffnessWeight
      PROCEDURE :: ShelfProfile          => ShelfProfile_BarotropicShelfWaves
      PROCEDURE :: ShelfMod              => ShelfMod_BarotropicShelfWaves
      PROCEDURE :: FillForcingTerm       => FillForcingTerm_BarotropicShelfWaves
      PROCEDURE :: StiffnessMatrixAction => StiffnessMatrixAction_BarotropicShelfWaves 
      PROCEDURE :: Residual              => Residual_BarotropicShelfWaves
      PROCEDURE :: VectorProduct         => VectorProduct_BarotropicShelfWaves
      PROCEDURE :: HBarInnerProduct      => HBarInnerProduct_BarotropicShelfWaves
      PROCEDURE :: HBarNormalize         => HBarNormalize_BarotropicShelfWaves
      PROCEDURE :: ComputeGrowthRates    => ComputeGrowthRates_BarotropicShelfWaves
      PROCEDURE :: CGInvert              => CGInvert_BarotropicShelfWaves
      PROCEDURE :: MatrixAction          => MatrixAction_BarotropicShelfWaves
      PROCEDURE :: ArnoldiProcess        => ArnoldiProcess_BarotropicShelfWaves
      PROCEDURE :: IRAM                  => IRAM_BarotropicShelfWaves
      PROCEDURE :: RecomputeEigenValues  => RecomputeEigenValues_BarotropicShelfWaves

      PROCEDURE :: WriteTecplot => WriteTecplot_BarotropicShelfWaves

   END TYPE BarotropicShelfWaves

CONTAINS
!
!
!==================================================================================================!
!------------------------------- Manual Constructors/Destructors ----------------------------------!
!==================================================================================================!
!
!
 SUBROUTINE Build_BarotropicShelfWaves( this )
 ! S/R Build
 !
 ! 
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS( BarotropicShelfWaves ), INTENT(inout) :: this
   ! Local
   INTEGER :: row, col, m, iEl, nS, nEl, maxIter, nDOF, nEpairs
   REAL(prec) :: x, sK, J, lmc, lmr, f0
   REAL(prec), ALLOCATABLE :: w(:), dMat(:,:), Hinv(:)
   
      CALL this % params % Build( )

      f0      = this % params % f0
      nS      = this % params % polyDeg
      nEl     = this % params % nElems
      nEpairs = this % params % nEpairs
      maxIter = this % params % maximumIterates
      nDOF    = (nS+1)*nEl      

      this % nEigenpairs = nEpairs
      this % nElems      = nEl
      this % nS          = nS
      this % maxIters    = maxIter
      this % tolerance   = this % params % tolerance
      this % nMaxToSolve = this % params % nMaxToSolve

      CALL this % cgStorage % Build( nS, GAUSS_LOBATTO, CG )

      CALL this % mesh % LoadDefaultMesh( this % cgStorage % interp, nEl )
      CALL this % mesh % ScaleTheMesh( this % params % xScale )

      ALLOCATE( this % eigenvalues(1:nEpairs), this % eigenFunctions(0:nS,1:nEl,1:nEpairs) )
      ALLOCATE( this % growthRates(1:nEpairs) )
      ALLOCATE( this % cgsol(0:nS,1:nEl), this % rhsVar(0:nS,1:nEl) )
      this % eigenvalues      = ZERO
      this % eigenfunctions   = ZERO
      this % growthRates      = ZERO
      this % cgsol            = ZERO
      this % rhsVar           = ZERO
      this % boundaryFlux     = ZERO
      this % dirichletMask    = ZERO
      this % dirichletMask(2) = ONE

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
                  sK = this % StiffnessWeight( x ) ! calculate the stiffness weight, "K" in d/dx( K(x) d/dx )
                  J = this % mesh % GetJacobianAtNode( iEl, m )

                  lmr = dMat(m,row) 
                  lmc = dMat(m,col) 

                  this % d2Mat(row,col,iEl) = this % d2Mat(row,col,iEl) - (sK*w(m)/J)*lmr*lmc

               ENDDO 
            ENDDO
         ENDDO 
      ENDDO 

      ALLOCATE( this % wFunc(0:nS,1:nEl), this % H(0:nS,1:nEl), this % G(0:nS,1:nEl) )
      ALLOCATE( Hinv(0:nS) )

      DO iEl = 1, nEl
         DO m = 0, nS
            x = this % mesh % GetPositionAtNode( iEl, m ) 
            this % H(m,iEl) = this % ShelfProfile( x )
         ENDDO
      ENDDO
     
      DO iEl = 1, nEl

         Hinv = ONE/this % H(0:nS,iEl)
         this % wFunc(:,iEl) = MATMUL( dMat, Hinv )
  
         DO m = 0, nS
            J = this % mesh % GetJacobianAtNode( iEl, m )
            this % wFunc(m,iEl) = f0*this % wFunc(m,iEl)*w(m)
           ! this % wFunc(m,iEl) = w(m)*J
         ENDDO

      ENDDO
         
      DEALLOCATE( w, dMat, Hinv )

      CALL this % FillForcingTerm( )

 END SUBROUTINE Build_BarotropicShelfWaves
!
!
! 
 SUBROUTINE Trash_BarotropicShelfWaves( this )
 ! S/R Trash
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS( BarotropicShelfWaves ), INTENT(inout) :: this
   
      CALL this % cgStorage % Trash( )
      CALL this % mesh % Trash( )
      DEALLOCATE( this % cgsol, this % rhsVar, this % d2Mat )
      DEALLOCATE( this % eigenvalues, this % eigenfunctions )
      DEALLOCATE( this % growthrates, this % G )

 END SUBROUTINE Trash_BarotropicShelfWaves
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
 FUNCTION ShelfProfile_BarotropicShelfWaves( this, x ) RESULT( h )
 ! ShelfProfile
 !
 ! =============================================================================================== !
   IMPLICIT NONE
   CLASS( BarotropicShelfWaves ) :: this
   REAL(prec)                    :: x
   REAL(prec)                    :: h
   ! Local
   REAL(prec) :: L1, L2, hcoast, hshelf, hdeep, xSlope

      L1     = this % params % Linner
      L2     = this % params % Lslope
      xSlope = this % params % xSlope
      hcoast = this % params % hcoast
      hshelf = this % params % hShelf
      hDeep  = this % params % hDeep

      h = hcoast + (hshelf - hcoast)*tanh( x/L1 ) + &
          HALF*(hdeep - hshelf)*( tanh((x-xSlope)/L2) + ONE )


 END FUNCTION ShelfProfile_BarotropicShelfWaves
!
!
!
 FUNCTION ShelfMod_BarotropicShelfWaves( this, x ) RESULT( h )
 ! ShelfMod
 !
 ! =============================================================================================== !
   IMPLICIT NONE
   CLASS( BarotropicShelfWaves ) :: this
   REAL(prec)                    :: x
   REAL(prec)                    :: h
   ! Local
   REAL(prec) :: L1, L2, hcoast, hshelf, hdeep, xSlope
   REAL(prec) :: dx, dL
   REAL(prec) :: s, ds, sech2slope, sech2inner

      L1     = this % params % Linner
      L2     = this % params % Lslope
      xSlope = this % params % xSlope
      hcoast = this % params % hcoast
      hshelf = this % params % hShelf
      hDeep  = this % params % hDeep
      dX     = this % params % dxSlope
      dL     = this % params % dLinner

      ! Note that hyperbolic secant is not an intrinsic function. We can, however, approximate
      ! it using the definition of the derivative - introduces O( ds**2 )error
      s  = ( x - xSlope )/L2
      ds = (10.0_prec**(-7))!*s
      sech2slope = HALF*(tanh(s+ds) - tanh(s-ds))/ds

      s  = ( x )/L1
      ds = (10.0_prec**(-7))!*s
      sech2inner = HALF*(tanh(s+ds) - tanh(s-ds))/ds

      h = HALF*(hDeep - hShelf)*dX/L2*sech2slope - &
          (hShelf - hcoast)*dL*x/(L1*L1)*sech2inner


 END FUNCTION ShelfMod_BarotropicShelfWaves
!
!
!
 SUBROUTINE FillForcingTerm_BarotropicShelfWaves( this )
 ! S/R FillForcingTerm
 !
 ! =============================================================================================== !
 ! DECLARATIONS 
   IMPLICIT NONE
   CLASS( BarotropicShelfWaves ), INTENT(inout) :: this
   ! Local
   REAL(prec) :: x, J, vbar, f0
   REAL(prec) :: gOnH(0:this % nS, 1:this % nElems)
   REAL(prec) :: DgOnHDx(0:this % nS, 1:this % nElems)
   REAL(prec) :: dMat(0:this % nS,0:this % nS)
   INTEGER    :: nS, nEl, iS, iEl


      nEl  = this % nElems
      nS   = this % nS
      vbar = this % params % vbar
      f0   = this % params % f0

      CALL this % cgStorage % GetDerivativeMatrix( dMat )

      DO iEl = 1, nEl
         DO iS = 0, nS

            x = this % mesh % GetPositionAtNode( iEl, iS )
            gOnH(iS,iEl) = this % ShelfMod( x )/this % H(iS,iEl)

         ENDDO
        
         DgOnHDx(0:nS, iEl) = MATMUL( dMat, gOnH(0:nS,iEl) )
 
         DO iS = 0, nS
            J = this % mesh % GetJacobianAtNode( iEl, iS )
            DgOnHDx(iS,iEl) = DgOnHDx(iS,iEl)/J
         ENDDO

      ENDDO

      this % G = vbar*( vbar*DgOnHDx + f0*gOnH )


 END SUBROUTINE FillForcingTerm_BarotropicShelfWaves
!
!
!
 FUNCTION StiffnessWeight( this, x ) RESULT( K )
 ! FUNCTION StiffnessWeight
 !
 ! This function returns the "K" used in the computation of 
 ! the operator  " d/dx( K(x) d/dx ) ."
 ! =============================================================================================== !
 ! DECLARATION
   IMPLICIT NONE
   CLASS( BarotropicShelfWaves ) :: this
   REAL(prec)                    :: x, K

      K = ONE/this % ShelfProfile( x )

 END FUNCTION StiffnessWeight
!
!
!
 SUBROUTINE DirichletConditions_BarotropicShelfWaves( thisElliptic )
 ! S/R DirichletConditions
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS( BarotropicShelfWaves ), INTENT(inout) :: thisElliptic
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

 END SUBROUTINE DirichletConditions_BarotropicShelfWaves
!
!
!
 SUBROUTINE Mask_BarotropicShelfWaves( thisElliptic, thisArray )
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
   CLASS( BarotropicShelfWaves ), INTENT(in) :: thisElliptic
   REAL(prec), INTENT(inout)       :: thisArray(0:thisElliptic % nS, 1:thisElliptic % nElems)
   ! LOCAL
   INTEGER :: iEl, nEl, nS
 
      nEl = thisElliptic % nElems
      nS  = thisElliptic % nS

      DO iEl = 1, nEl-1
         thisArray(0,iEl+1) = ZERO
      ENDDO

      ! If dirichlet boundary conditions are applied at either boundary, the masking routine 
      ! (assuming it is taking in a residual for an iterative solver) needs to force the 
      ! incoming array to zero at those boundary points. Here, this is achieved by using the
      ! "dirichletMask" attribute.
      thisArray(0,1)    = thisArray(0,1)*thisElliptic % dirichletMask(1)
      thisArray(nS,nEl) = thisArray(nS,nEl)*thisElliptic % dirichletMask(2) 

 END SUBROUTINE Mask_BarotropicShelfWaves
!
!
!
 SUBROUTINE UnMask_BarotropicShelfWaves( thisElliptic, thisArray )
 ! S/R UnMask
 !
 ! This subroutine unmasks 'thisArray' along the nodes shared by two elements.
 ! A copy of the solution retained in the "MASK" operation is used to unmask.
 !
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS( BarotropicShelfWaves ), INTENT(in) :: thisElliptic
   REAL(prec), INTENT(inout)                 :: thisArray(0:thisElliptic % nS, 1:thisElliptic % nElems)
   ! LOCAL
   INTEGER :: iEl, nS, nEl
 
      nS  = thisElliptic % nS
      nEl = thisElliptic % nElems
      
      DO iEl = 1, nEl-1
         thisArray(0,iEl+1) = thisArray(nS,iEl)
      ENDDO


 END SUBROUTINE UnMask_BarotropicShelfWaves
!
!
!
 SUBROUTINE UnMaskSolution_BarotropicShelfWaves( thisElliptic )
 ! S/R UnMaskSolution
 !
 ! This subroutine unmasks the solution attribute and applies the Dirichlet boundary conditions
 ! if appropriate.
 !
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS( BarotropicShelfWaves ), INTENT(inout) :: thisElliptic
      
      CALL thisElliptic % UnMask( thisElliptic % cgsol )
      CALL thisElliptic % DirichletConditions( )


 END SUBROUTINE UnMaskSolution_BarotropicShelfWaves
!
!
!
 SUBROUTINE GlobalSum_BarotropicShelfWaves( thisElliptic, thisArray )
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
   CLASS( BarotropicShelfWaves ), INTENT(in) :: thisElliptic
   REAL(prec), INTENT(inout)                 :: thisArray(0:thisElliptic % nS, 1:thisElliptic % nElems)
   ! LOCAL
   INTEGER    :: iEl, nEl, nS
   REAL(prec) :: temp
 
      nS  = thisElliptic % nS
      nEl = thisElliptic % nElems
      
      DO iEl = 2, nEl
         temp = thisArray(0, iEl) + thisArray(nS,iEl-1)
         thisArray(0, iEl)   = temp
         thisArray(nS,iEl-1) = temp
      ENDDO

 END SUBROUTINE GlobalSum_BarotropicShelfWaves
!
!
!
 SUBROUTINE ApplyOperator_BarotropicShelfWaves( thisElliptic, u, Lu )
 ! S/R ApplyOperator
 !
 ! This subroutine loops over all of the elements an applies the internal matrix operations
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS( BarotropicShelfWaves ), INTENT(in) :: thisElliptic
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
        
      

 END SUBROUTINE ApplyOperator_BarotropicShelfWaves
!
!
! 
 FUNCTION StiffnessMatrixAction_BarotropicShelfWaves( this, s ) RESULT( As )
 !
 ! 
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS( BarotropicShelfWaves ) :: this
   REAL(prec)                    :: s(0:this % nS, 1:this % nElems)
   REAL(prec)                    :: As(0:this % nS, 1:this % nElems)

      CALL this % UnMask( s )

      CALL this % ApplyOperator( s, As )

      CALL this % GlobalSum( As )

      CALL this % Mask( As )
 
 END FUNCTION StiffnessMatrixAction_BarotropicShelfWaves
!
!
!
 FUNCTION Residual_BarotropicShelfWaves( this, s ) RESULT( r )
 ! Residual
 ! 
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS( BarotropicShelfWaves ) :: this
   REAL(prec)                    :: s(0:this % nS, 1:this % nElems)
   REAL(prec)                    :: r(0:this % nS, 1:this % nElems)
   ! Local
   REAL(prec) :: As(0:this % nS, 1:this % nElems)
   integer    :: nEl, nS

      nEl = this % nElems
      nS  = this % nS

      ! This call calculates the 2nd Derivative Operator to s
      As = this % StiffnessMatrixAction( s )

      r = this % wFunc*this % rhsVar ! This is the inner-product weight applied to solution array

      CALL this % GlobalSum( r )
      CALL this % Mask( r )
      r = r - As

      ! If the dirichletMask is set to zero at either boundary, then the dirichlet boundary
      ! conditions are enforced. If it is set to one at either boundary, then a boundary
      ! flux should be added to the matrix action for the associated degree of freedom. To
      ! implement such a Neumann boundary condition, the boundary flux is added to the 
      ! rhsVar for the first and last nodes if the dirichletMask is set to one.
     ! r(0,1)    = r(0,1) + this % boundaryFlux(1)*this % dirichletMask(1)
     ! r(nS,nEl) = r(nS,nEl) - this % boundaryFlux(2)*this % dirichletMask(2)

     ! CALL this % Mask( r )

 END FUNCTION Residual_BarotropicShelfWaves
!
!
!
 FUNCTION VectorProduct_BarotropicShelfWaves( this, u, v ) RESULT( uDotV )
 ! VectorProduct
 ! 
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS( BarotropicShelfWaves ) :: this
   REAL(prec)                    :: u(0:this % nS, 1:this % nElems)
   REAL(prec)                    :: v(0:this % nS, 1:this % nElems)
   REAL(prec)                    :: uDotV
   ! Local
   integer    :: iEl

      CALL this % Mask( u )
      CALL this % Mask( v )
      uDotV = ZERO

      DO iEl = 1, this % nElems
         uDotV = uDotV + DOT_PRODUCT( u(:,iEl), v(:,iEl) )
      ENDDO


 END FUNCTION VectorProduct_BarotropicShelfWaves
!
!
!
 FUNCTION HBarInnerProduct_BarotropicShelfWaves( this, u, v ) RESULT( uDotV )
 ! HBarInnerProduct
 ! 
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS( BarotropicShelfWaves ) :: this
   REAL(prec)                    :: u(0:this % nS, 1:this % nElems)
   REAL(prec)                    :: v(0:this % nS, 1:this % nElems)
   REAL(prec)                    :: uDotV
   ! Local
   INTEGER    :: iEl
   REAL(prec) :: uWeighted(0:this % nS, 1:this % nElems)

      CALL this % Mask( u )
      CALL this % Mask( v )
      uDotV = ZERO
      uWeighted = u*this % wFunc

      DO iEl = 1, this % nElems
         uDotV = uDotV + DOT_PRODUCT( uWeighted(:,iEl), v(:,iEl) )
      ENDDO


 END FUNCTION HBarInnerProduct_BarotropicShelfWaves
!
!
!
 SUBROUTINE HBarNormalize_BarotropicShelfWaves( this )
 ! S/R HBarNormalize
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS( BarotropicShelfWaves ), INTENT(inout) :: this
   ! LOCAL
   REAL(prec) :: s(0:this % nS, 1:this % nElems)
   REAL(prec) :: Ws(0:this % nS, 1:this % nElems)
   REAL(prec) :: sMag
   INTEGER    :: i, m

      m = this % nEigenPairs
      DO i = 1, m

         s    = this % eigenfunctions(:,:,i)
         Ws   = s*this % wFunc
  
         sMag = SQRT( ABS(this % VectorProduct( s, Ws )) )
         
         this % eigenfunctions(:,:,i) = s/sMag

         CALL this % UnMask( this % eigenfunctions(:,:,i) )

      ENDDO 

 END SUBROUTINE HBarNormalize_BarotropicShelfWaves
!
!
!
 SUBROUTINE ComputeGrowthRates_BarotropicShelfWaves( this )
 ! S/R ComputeGrowthRates
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS( BarotropicShelfWaves ), INTENT(inout) :: this
   ! LOCAL
   REAL(prec) :: s(0:this % nS, 1:this % nElems)
   REAL(prec) :: G(0:this % nS, 1:this % nElems)
   REAL(prec) :: c
   INTEGER    :: i, m

      m = this % nEigenPairs
      G = this % G

      DO i = 1, m

         s = this % eigenfunctions(:,:,i)
         c = this % eigenvalues(i)
         
         this % growthRates(i) = this % VectorProduct( s, G )*c

         CALL this % UnMask( this % eigenfunctions(:,:,i) )

      ENDDO 

 END SUBROUTINE ComputeGrowthRates_BarotropicShelfWaves
!
!
!
 SUBROUTINE CGInvert_BarotropicShelfWaves( this, ioerr )
 !  S/R CGInvert
 !
 !  This subroutine solves the system Ax = b using the un-preconditioned conjugate gradient method.
 !  The matrix action and residual routines are supplied by a non-abstracted type-extension of
 !  ConjugateGradient. These routines should return an array indexed from 1 to nDOF. Thus,
 !  in addition to a MatrixAction and Residual, the user should map their data-structure to a 1-D 
 !  array.  
 !
 !  On output ioerr is set to an error checking flag. 
 !  If ioerr ==  0, the method converged within the maximum number of iterations.
 !     ioerr == -1, the method did not converge within the maximum number of iterations.
 !     ioerr == -2, something that is not caught by the current construct happened.
 !  
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS( BarotropicShelfWaves ), INTENT(inout) :: this
   INTEGER, INTENT(out)                         :: ioerr
   ! LOCAL
   INTEGER    :: iter
   INTEGER    :: nIt
   REAL(prec) :: TOL
   REAL(prec) :: r(0:this % nS, 1:this % nElems)
   REAL(prec) :: v(0:this % nS, 1:this % nElems)
   REAL(prec) :: z(0:this % nS, 1:this % nElems)
   REAL(prec) :: rNorm
   REAL(prec) :: a, b, num, den, r0
   
      ioerr = -2
      nIt = this % maxIters
      TOL = this % tolerance
      
      r = this % Residual( this % cgsol )
      r0 = sqrt( this % VectorProduct( r, r ) ) 
 
      DO iter = 1,nIt ! Loop over the CG iterates

         num = this % VectorProduct( r, r )

         IF( iter == 1) THEN
            v = r
         ELSE
            b = num/den
            v = r + b*v
         ENDIF

         z = this % StiffnessMatrixAction( v )
         a = num/this % VectorProduct( v, z )
         this % cgsol = this % cgsol + a*v
         r = r - a*z
         
         den = num

         rNorm = sqrt( this % VectorProduct(r,r) )
        ! this % resi(iter) = rNorm
         IF( sqrt(rNorm)/r0 < TOL ) then
           EXIT
           ioerr = 0
         ENDIF

      ENDDO ! iter, loop over the CG iterates 

      IF( SQRT(rNorm) > TOL )THEN
         PRINT*, 'MODULE IterativeSolvers : ConjugateGradient failed to converge '
         PRINT*, 'Last L-2 residual : ', sqrt(rNorm)
         ioerr=-1
      ENDIF   

 END SUBROUTINE CGInvert_BarotropicShelfWaves
!
!
!
 FUNCTION MatrixAction_BarotropicShelfWaves( this, s ) RESULT( Ms )
 !  S/R MatrixAction
 !
 !    This subroutine calculates the matrix action for the discretized Sturm-Liouville shelf
 !    problem.   
 !                 c As = W s  ===>  c s = A^{-1}Ws
 !
 !    Given a variable "s", the matrix action is obtained by calculating Ws and applying A^{-1},
 !    where A is the "StiffnessMatrix" Operator.
 !  
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS( BarotropicShelfWaves ) :: this
   REAL(prec)                    :: s(0:this % nS, 1:this % nElems)
   REAL(prec)                    :: Ms(0:this % nS, 1:this % nElems)
   ! Local
   INTEGER :: ioerr, iEl, iS
   REAL(prec) :: w(0:this % nS), J
                    
      this % rhsVar = s    ! Set the solution variable that the inner-product weight will be applied to
      this % cgSol  = s ! Initial guess for the solution v = inverse(A)*W*s

      CALL this % CGInvert( ioerr )

      Ms = this % cgSol

      CALL this % UnMask( Ms )

 END FUNCTION MatrixAction_BarotropicShelfWaves
!
!
!
 SUBROUTINE ArnoldiProcess_BarotropicShelfWaves( this, U, H, mStart, earlyTermination )
 !  S/R ArnoldiProcess
 !
 !     This subroutine computes the iterates of the Arnoldi process starting from mStart and
 !     terminating at mEnd. Reorthogonalization is done to ensure that the basis vectors in the
 !     "columns" of U.
 !
 !     If the process terminates early, the integer "earlyTermination" is set to equal "j", the
 !     step in the process where h_{j+1,j} = 0. Otherwise, earlyTermination = -1.
 !
 !     ** The "basis function" U(:,:,mStart) is assumed to be set and normallized before being
 !        passed into this subroutine.
 !
 ! =============================================================================================== !
 !
   IMPLICIT NONE
   CLASS( BarotropicShelfWaves ), INTENT(inout) :: this
   INTEGER, INTENT(in)                          :: mStart
   REAL(prec), INTENT(inout)                    :: U(0:this % nS, 1:this % nElems, 1:this % nMaxToSolve+1)
   REAL(prec), INTENT(inout)                    :: H(1:this % nMaxToSolve+1,1:this % nMaxToSolve)

   INTEGER, INTENT(out)                         :: earlyTermination
   ! Local
   INTEGER :: i, j
   REAL(prec) :: w(0:this % nS,1:this % nElems), v(0:this % nS, 1:this % nElems)
   REAL(prec) :: orthoErr

      earlyTermination = -1

      DO j = mStart, this % nMaxToSolve
         
         v = U(:,:,j)
         w = this % MatrixAction( v )
         U(:,:,j+1) = w

         ! The new basis vector is obtained by multiplying the previous basis vector by the matrix
         ! and orthogonalizing wrt to all of the previous basis vectors using a Gram-Schmidt process.
         DO i = 1, j
            v      = U(:,:,i)
            H(i,j) = this % VectorProduct( v, w )
            w      = w - H(i,j)*v
         ENDDO

         ! Reorthogonalize -- to do this, we first compute the vector product of "w"  with each of 
         ! the previous basis functions. This additional component (that arises from roundoff error)
         ! is removed from the new basis vector.
         DO i = 1, j
            v        = U(:,:,i)
            orthoErr = this % VectorProduct( v, w )
            w        = w - orthoErr*v
            H(i,j)   = H(i,j) + orthoErr
         ENDDO

         H(j+1,j) = sqrt( this % VectorProduct( w, w ) )
       
         IF( H(j+1,j) <= this % tolerance  )THEN
            PRINT*,'S/R ArnoldiProcess : Early Termination at iterate', j
            earlyTermination = j
            EXIT
         ENDIF

         U(:,:,j+1) = w/H(j+1,j)

      ENDDO

!      DO j = 1, this % nMaxToSolve + 1
!         PRINT*, H(j,:)
!      ENDDO
!      STOP
 END SUBROUTINE ArnoldiProcess_BarotropicShelfWaves
!
!
!
 SUBROUTINE IRAM_BarotropicShelfWaves( this )
 ! S/R IRAM
 !
 !  This subroutine calculates the dominant eigenpairs of the "Shelf-Wave Sturm-Liouville" problem
 !  using the Implicitly Restarted Arnoldi Method. The number of desired eigenpairs is set by the
 !  parameter "nEigenPairs" in "runtime.params".
 !
 ! =============================================================================================== !
 ! 
   IMPLICIT NONE
   CLASS( BarotropicShelfWaves ), INTENT(inout) :: this
   ! Local
   REAL(prec) :: U(0:this % nS, 1:this % nElems, 1:this % nMaxToSolve+1)
   REAL(prec) :: w(0:this % nS)
   REAL(prec) :: Unew(0:this % nS, 1:this % nElems, 1:this % nMaxToSolve)
   REAL(prec) :: H(1:this % nMaxToSolve+1,1:this % nMaxToSolve)
   REAL(prec) :: Q(1:this % nMaxToSolve,1:this % nMaxToSolve)
   REAL(prec) :: R(1:this % nMaxToSolve,1:this % nMaxToSolve)
   REAL(prec) :: Qt(1:this % nMaxToSolve,1:this % nMaxToSolve)
   REAL(prec) :: Hsquare(1:this % nMaxToSolve,1:this % nMaxToSolve)
   REAL(prec) :: hEvalsR(1:this % nMaxToSolve), hEvalsI(1:this % nMaxToSolve)
   REAL(prec) :: hEvecs(1:this % nMaxToSolve,1:this % nMaxToSolve)
   REAL(prec) :: hShifted(1:this % nMaxToSolve,1:this % nMaxToSolve)
   REAL(prec) :: temp1(1:this % nMaxToSolve), temp2(1:this % nMaxToSolve)
   REAL(prec) :: work(1:this % nMaxToSolve,1:this % nMaxToSolve)
   REAL(prec) :: lapackWork(1:this % nMaxToSolve*5)
   REAL(prec) :: r0, resi, tolerance, x, L
   INTEGER    :: terminationIndex, lwork, ioerr
   INTEGER    :: m, k, j, i, iter, row, mstart, maxIters, iEl, iS, nEl, nS, col

      nEl       = this % nElems
      nS        = this % nS
      m         = this % nEigenPairs ! The number of desired eigenpairs
      k         = this % nMaxToSolve ! The dimension of the subspace used to calculate the desired eigenpairs
      maxIters  = this % maxIters
      tolerance = this % tolerance
      lwork     = 5*k

      CALL this % cgStorage % GetQuadratureWeights( w )
      mStart = 1
      H = ZERO
      U = ZERO
      DO iEl = 1, nEl
         DO iS = 0, nS
            x = this % mesh % GetPositionAtNode( iEl, iS )
            L = this % params % xScale
            U(iS,iEl,1) = tanh( TWO*x/L )
         ENDDO
      ENDDO
      Unew = U(:,:,1:k)
      U(:,:,1) = U(:,:,1)/SQRT( this % VectorProduct(Unew(:,:,1),Unew(:,:,1)) )

      DO iter = 1, maxIters
         
         ! For iterates beyond the first on, the first m-columns are retained after deflation
         ! and the remaining columns are recomputed to refine the estimate of the small eigenvalues
         ! This aids in further deflation of spectrum associated with small eigenvalues.

         CALL this % ArnoldiProcess( U, H, mStart, terminationIndex )
         mStart = m
         ! For now we assume that terminationIndex = -1
         r0 = H(k+1,k)
         resi = H(m+1,m)

         Hsquare = H(1:k,1:k)

         ! Calculate the eigenvalues and eigenvectors of "Hsquare"
         CALL DGEEV( 'N', 'V', k, Hsquare, k, &
                     hEvalsR, hEvalsI, &
                     work, k, &
                     hEvecs, k, &
                     lapackWork, lwork, ioerr )

         CALL SortEigenPairs( hEvalsR, hEvalsI, hEvecs, k )

         ! If the residual is smaller than the specified tolerance, then we can compute the 
         ! eigenfunctions, and exit the subroutine.
         IF( resi <= tolerance )THEN
            DO iEl = 1, nEl
               DO iS = 0, nS
                  DO col = 1, m
                     this % eigenfunctions(iS,iEl,col) = DOT_PRODUCT( U(iS,iEl,1:k), hEvecs(:,col) )
                  ENDDO
               ENDDO
            ENDDO
            DO col = 1, m
               this % eigenfunctions(:,:,col) =  this % eigenfunctions(:,:,col)/&
                                           sqrt( this % VectorProduct( this % eigenfunctions(:,:,col), &
                                                                        this % eigenfunctions(:,:,col) ) )
               CALL this % UnMask( this % eigenfunctions(:,:,col) )
            ENDDO
            this % eigenValues = hEvalsR(1:m) 
            CALL this % RecomputeEigenValues( )
            RETURN
         ENDIF

         ! Otherwise, we attempt to deflate the portion of the spectrum associated with the smaller
         ! eigenvalues using a QR algorithm with shifts
        
         DO i = m+1, k
            R = ZERO
            Q = ZERO
            Hshifted = H(1:k,1:k)
            DO row = 1, k
               Hshifted(row,row) = Hshifted(row,row) - hEvalsR(i)
            ENDDO

            ! Now, an orthonormal basis for the shifted Hessenburg matrix is constructed using
            ! a QR decomposition.
            Q(:,1) = Hshifted(:,1)
            temp1  = Q(:,1)
            R(1,1) = SQRT( DOT_PRODUCT(temp1,temp1) )
            Q(:,1) = Q(:,1)/R(1,1)
!            print*, hEvalsR(i)
!            STOP
            DO col = 2, k
            
               temp1 = Hshifted(:,col)
               Q(:,col) = temp1
               DO j = 1, col-1
            
                  temp2    = Q(:,j)
                  R(j,col) = DOT_PRODUCT(temp1,temp2)
                  Q(:,col) = Q(:,col) - R(j,col)*temp2
 
               ENDDO

               temp1      = Q(:,col)
               R(col,col) = SQRT( DOT_PRODUCT(temp1,temp1) )
               Q(:,col)   = Q(:,col)/R(col,col)

            ENDDO
            
            ! And reassemble Hsquare with the deflated basis
            Hsquare = MATMUL( R, Q )
            DO row = 1, k
               Hsquare(row,row) = Hsquare(row,row) + hEvalsR(i)
            ENDDO

         ENDDO

         ! Project the subspace onto the "small-eigenvalue-deflated" basis
         Unew = ZERO
         DO iEl = 1, nEl
            DO iS = 0, nS
               DO col = 1, m
                  Unew(iS,iEl,col) = DOT_PRODUCT( U(iS,iEl,1:k), Q(:,col) )
               ENDDO
            ENDDO
         ENDDO
   
         U(:,:,1:m) = Unew(:,:,1:m)
         H = ZERO
         H(1:m,1:m) = Hsquare(1:m,1:m)
!         DO i = 1,m
!             PRINT*, H(i,1:m)
!         ENDDO

!         PRINT*, iter
      ENDDO

 END SUBROUTINE IRAM_BarotropicShelfWaves
!
!
!
 SUBROUTINE RecomputeEigenValues_BarotropicShelfWaves( this )
 ! S/R RecomputeEigenValues
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS( BarotropicShelfWaves ), INTENT(inout) :: this
   ! Local
   REAL(prec) :: lambda, l
   REAL(prec) :: s(0:this % nS, 1:this % nElems)
   REAL(prec) :: As(0:this % nS, 1:this % nElems)
   REAL(prec) :: Ws(0:this % nS, 1:this % nElems)
   REAL(prec) :: sMag
   INTEGER    :: i, m

!      PRINT*, '============================================'
!      PRINT*, '========== RecomputeEigenValues ============'
!      PRINT*, '============================================'

!      m = this % nEigenPairs
!      DO i = 1, m

!         s    = this % eigenfunctions(:,:,i)
!         sMag = this % VectorProduct( s, s )
!         Ms   = this % MatrixAction( s )
  
!         lambda = this % VectorProduct( s, Ms )/sMag
!         l      = this % eigenvalues(i)
!         PRINT*, ' eigenvalue relative change = ', abs( lambda - l )/abs( l )

!         this % eigenvalues(i) = l

!      ENDDO 

!      PRINT*, '============================================'

      PRINT*, '============================================'
      PRINT*, '========== RecomputeEigenValues ============'
      PRINT*, '============================================'

      CALL this % HBarNormalize( )

      m = this % nEigenPairs
      DO i = 1, m

         s    = this % eigenfunctions(:,:,i)
         As   = this % StiffnessMatrixAction( s )
         Ws   = s*this % wFunc
  
         lambda = this % VectorProduct( s, Ws )/this % VectorProduct( s, As )
         l      = this % eigenvalues(i)
         PRINT*, ' eigenvalue relative change = ', abs( lambda - l )/abs( l )

         this % eigenvalues(i) = l

      ENDDO 

      PRINT*, '============================================'


 END SUBROUTINE RecomputeEigenValues_BarotropicShelfWaves
!
!
!
 SUBROUTINE SortEigenPairs( eReal, eImag, eVecs, N )
 ! S/R SortEigenPairs
 !
 !   Sorts eigenpairs in order of increasing absolute value of real part of eigenvalues. 
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   INTEGER, INTENT(in)  :: N
   REAL(prec), INTENT(inout) :: eReal(1:N), eImag(1:N)
   REAL(prec), INTENT(inout) :: eVecs(1:N,1:N)
   ! LOCAL
   INTEGER :: i, j
   REAL(prec) :: tempEreal, tempEimag
   REAL(prec) :: tempEvec(1:N)


    DO i = 2,  N
       j = i
       DO WHILE( j > 1 )
          IF(  ABS(eReal(j-1)) < ABS(eReal(j)) )THEN
             !Swap outArray(j) outArray(j-1)
             tempEreal  = eReal(j)
             eReal(j)   = eReal(j-1)
             eReal(j-1) = tempEreal 

             tempEimag  = eImag(j)
             eImag(j)   = eImag(j-1)
             eImag(j-1) = tempEimag

             tempEvec       = eVecs(1:N,j)
             eVecs(1:N,j)   = eVecs(1:N,j-1)
             eVecs(1:N,j-1) = tempEvec  
        
             j = j-1
          ELSE
             EXIT
          ENDIF
       ENDDO

    ENDDO

 END SUBROUTINE SortEigenPairs
!
!
!==================================================================================================!
!-------------------------------------- File I/O Routines -----------------------------------------!
!==================================================================================================!
!
!
 SUBROUTINE WriteTecplot_BarotropicShelfWaves( thisElliptic )
 ! S/R WriteTecplot
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS( BarotropicShelfWaves ), INTENT(inout) :: thisElliptic
   ! Local
   INTEGER      :: iS, iEl, nS, nEl, fUnit, i
   REAL(prec)   :: x
   CHARACTER(2) :: eigenChar

      nEl = thisElliptic % nElems
      nS  = thisElliptic % nS

      CALL thisElliptic % UnMaskSolution( )

      OPEN( UNIT = NewUnit(fUnit), &
            FILE = 'eigenfunctions.curve', &
            FORM = 'FORMATTED' )

      DO i = 1, thisElliptic % nEigenPairs

         WRITE(eigenChar,'(I2.2)')i

         WRITE(fUnit,*)'#e-funct'//eigenChar
         DO iEl = 1, nEl
            DO iS = 0, nS
              x = thisElliptic % mesh % GetPositionAtNode(iEl,iS)
              WRITE(fUnit,*) x, thisElliptic % eigenfunctions(iS,iEl,i)
            ENDDO
         ENDDO

      ENDDO 

      CLOSE(fUnit)

      OPEN( UNIT = NewUnit(fUnit), &
            FILE = 'shelf.curve', &
            FORM = 'FORMATTED' )

         WRITE(fUnit,*)'#bathymetry'
         DO iEl = 1, nEl
            DO iS = 0, nS
              x = thisElliptic % mesh % GetPositionAtNode(iEl,iS)
              WRITE(fUnit,*) x, thisElliptic % ShelfProfile(x)
            ENDDO
         ENDDO


      CLOSE(fUnit)

      OPEN( UNIT = NewUnit(fUnit), &
            FILE = 'forcing.curve', &
            FORM = 'FORMATTED' )

         WRITE(fUnit,*)'#vts'
         DO iEl = 1, nEl
            DO iS = 0, nS
              x = thisElliptic % mesh % GetPositionAtNode(iEl,iS)
              WRITE(fUnit,*) x, thisElliptic % G(iS,iEl)
            ENDDO
         ENDDO


      CLOSE(fUnit)

      OPEN( UNIT = NewUnit(fUnit), &
            FILE = 'eigenvalues.curve', &
            FORM = 'FORMATTED' )

      WRITE(fUnit,*)'#eigenvalues'

      DO i = 1, thisElliptic % nEigenPairs
         WRITE(fUnit,*) i, thisElliptic % eigenvalues(i)
      ENDDO

      CLOSE(fUnit)

            OPEN( UNIT = NewUnit(fUnit), &
            FILE = 'growthrates.curve', &
            FORM = 'FORMATTED' )

      WRITE(fUnit,*)'#eigenvalues'

      DO i = 1, thisElliptic % nEigenPairs
         WRITE(fUnit,*) i, thisElliptic % growthrates(i)
      ENDDO

      CLOSE(fUnit)
   
 END SUBROUTINE WriteTecplot_BarotropicShelfWaves

END MODULE BarotropicShelfWaves_Class
