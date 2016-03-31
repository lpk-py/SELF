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
!   ** This version of Elliptic1D_Class is modified to store multiple solution, source, d2Mat, 
!      boundaryFlux, and dirichletMask arrays so that it can be used to vertically integrate
!      hydrostatic balance and the incompressibility condition.
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
! src/hydrostaticprimitive/
USE StackedHexMeshClass
USE HydrostaticParams_Class



IMPLICIT NONE

   TYPE EllipticStorage
      INTEGER                  :: nLayers, nS
      TYPE( SegmentMesh )      :: mesh
      REAL(prec), ALLOCATABLE  :: cgsol(:,:)
      REAL(prec), ALLOCATABLE  :: source(:,:)
      REAL(prec), ALLOCATABLE  :: d2Mat(:,:,:) ! matrix for d/dx( k(x) d/dx ) operator
      REAL(prec)               :: boundaryFlux(1:2)
      REAL(prec)               :: dirichletMask(1:2)

      CONTAINS

      PROCEDURE :: Build => Build_EllipticStorage
      PROCEDURE :: Trash => Trash_EllipticStorage

      PROCEDURE :: FillSource          => FillSource_EllipticStorage
      PROCEDURE :: DirichletConditions => DirichletConditions_EllipticStorage

      PROCEDURE :: Mask           => Mask_EllipticStorage
      PROCEDURE :: UnMask         => UnMask_EllipticStorage
      PROCEDURE :: UnMaskSolution => UnMaskSolution_EllipticStorage
      PROCEDURE :: GlobalSum      => GlobalSum_EllipticStorage
      PROCEDURE :: ApplyOperator  => ApplyOperator_EllipticStorage

      PROCEDURE :: MatrixAction  => MatrixAction_EllipticStorage
      PROCEDURE :: Residual      => Residual_EllipticStorage
      PROCEDURE :: VectorProduct => VectorProduct_EllipticStorage

      PROCEDURE :: ConjugateGradient => ConjugateGradient_EllipticStorage

   END TYPE EllipticStorage

   TYPE :: Elliptic1D
      INTEGER                              :: maxIters  ! Maximum number of iterates
      REAL(prec)                           :: tolerance ! tolerance for convergence
      REAL(prec), ALLOCATABLE              :: resi(:)   ! An array for storing the residual at each iterate
      INTEGER                              :: nHElems, nS
      TYPE( NodalStorage_1D )              :: cgStorage
      TYPE( EllipticStorage ), ALLOCATABLE :: ellStorage(:,:,:)
      
      CONTAINS
      
      PROCEDURE :: Build => Build_Elliptic1D
      PROCEDURE :: Trash => Trash_Elliptic1D

!      PROCEDURE :: WriteTecplot => WriteTecplot_Elliptic1D

      PROCEDURE :: SolveAll => SolveAll_Elliptic1D
   END TYPE Elliptic1D

CONTAINS
!
!
!==================================================================================================!
!------------------------------- Manual Constructors/Destructors ----------------------------------!
!==================================================================================================!
!
!
 SUBROUTINE Build_EllipticStorage( this, nS, nLayers, cgStorage, z )
 !
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS( EllipticStorage ), INTENT(inout) :: this
   INTEGER, INTENT(in)                     :: nS, nLayers
   TYPE( NodalStorage_1D )                 :: cgStorage
   REAL(prec), INTENT(in)                  :: z(1:2, 1:nLayers)
   ! Local
   INTEGER :: row, col, m, iEl
   REAL(prec), ALLOCATABLE :: w(:), dMat(:,:)
   REAL(prec) :: x, J, lmc, lmr

      this % nLayers = nLayers
      this % nS      = nS
      ALLOCATE( this % cgsol(0:nS,1:nLayers), this % source(0:nS,1:nLayers) )
      this % cgsol         = ZERO
      this % source        = ZERO
      this % boundaryFlux  = ZERO
      this % dirichletMask = ZERO
      CALL this % mesh % Build( nLayers, nS ) 
      
      DO iEl = 1, nLayers
         CALL this % mesh % elements(iEl) % Build( cgStorage % interp, z(1,iEl), z(2,iEl) )
      ENDDO


      CALL this % FillSource( )
      ! In 1-D, for a fixed mesh and stiffness coefficient, the second-derivative matrix
      ! for each element can be constructed and stored once. This way, to compute the elliptic
      ! matrix action, only one matrix-vector multiply is needed for each element. The alternative
      ! is to apply two matrix-vector operations with the 1st derivative matrices.
      !
      ! Here, we opt to build and store the second-derivative operator for each element.

      ALLOCATE( w(0:nS), dMat(0:nS,0:nS) )
      ALLOCATE( this % d2Mat(0:nS,0:nS,1:nLayers) )
      CALL cgStorage % GetQuadratureWeights( w )
      CALL cgStorage % GetDerivativeMatrix( dMat )

      DO iEl = 1, nLayers
         DO col = 0, nS
            DO row = 0, nS
               this % d2Mat(row,col,iEl) = ZERO ! Initialize the sum to zero
               DO m = 0, nS 

                  x = this % mesh % GetPositionAtNode( iEl, m ) 
                  J = this % mesh % GetJacobianAtNode( iEl, m )

                  lmr = dMat(m,row) 
                  lmc = dMat(m,col) 

                  this % d2Mat(row,col,iEl) = this % d2Mat(row,col,iEl) - (w(m)/J)*lmr*lmc

               ENDDO 
            ENDDO 
         ENDDO 
      ENDDO 

      DEALLOCATE( w, dMat )

 END SUBROUTINE Build_EllipticStorage
!
!
!
 SUBROUTINE Build_Elliptic1D( this, nS, stackedHMesh, params )
 ! S/R Build
 !
 ! 
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS( Elliptic1D ), INTENT(inout)    :: this
   INTEGER, INTENT(in)                   :: nS
   TYPE( StackedHexMesh ), INTENT(in)    :: stackedHMesh
   TYPE( HydrostaticParams ), INTENT(in) :: params
   ! Local
   INTEGER    :: iEl, nHElems, nLayers, maxIter, elID, i, j, iLayer
   REAL(prec) ::  x, y, z(1:2,1:stackedHMesh % nLayers)
   
   
      nHElems = stackedHMesh % nHElems
      nLayers = stackedHMesh % nLayers

      maxIter = params % maximumIterates

      this % nHElems   = nHElems
      this % nS        = nS
      this % maxIters  = maxIter
      this % tolerance = params % tolerance

      CALL this % cgStorage % Build( nS, GAUSS_LOBATTO, CG )
      

      ALLOCATE( this % resi(0:maxIter) )
      this % resi = ZERO 

      ALLOCATE( this % ellStorage(0:nS,0:nS,1:nHElems) ) 

      DO iEl = 1, nHElems
         DO j = 0, nS
            DO i = 0, nS

               DO iLayer = 1, nLayers
                  ! Get the element ID for this layer
                  elID = stackedHMesh % stackMap(iEl, iLayer)
                  CALL stackedHMesh % GetBoundaryLocationAtNode( elID, x, y, z(1,iLayer), i, j, bottom )
                  CALL stackedHMesh % GetBoundaryLocationAtNode( elID, x, y, z(2,iLayer), i, j, top )
               ENDDO

               CALL this % ellStorage(i,j,iEl) % Build( nS, nLayers, this % cgStorage, z )
             
            ENDDO
         ENDDO
      ENDDO


 END SUBROUTINE Build_Elliptic1D
!
!
!
  SUBROUTINE Trash_EllipticStorage( this )
 ! S/R Trash
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS( EllipticStorage ), INTENT(inout) :: this
   
      CALL this % mesh % Trash( )
      DEALLOCATE( this % cgsol, this % source, this % d2Mat )

 END SUBROUTINE Trash_EllipticStorage
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
   ! Local
   INTEGER :: iEl, i, j, nHElems, nS
     
      nHElems = this % nHElems
      nS      = this % nS
 
      CALL this % cgStorage % Trash( )
      DEALLOCATE( this % resi )

      DO iEl = 1, nHElems
         DO j = 0, nS
            DO i = 0, nS
               CALL this % ellStorage(i,j,iEl) % Trash( )
            ENDDO
         ENDDO
      ENDDO

      DEALLOCATE( this % ellStorage )
      
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
 SUBROUTINE FillSource_EllipticStorage( thisElliptic )
 !
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS( EllipticStorage ), INTENT(inout) :: thisElliptic
   ! LOCAL
   INTEGER    :: iS, iEl, nS, nEl
   REAL(prec) :: x

      nS  = thisElliptic % nS
      nEl = thisElliptic % nLayers

      DO iEl = 1, nEl
         DO iS = 0, nS
            x = thisElliptic % mesh % GetPositionAtNode(iEl,iS)
            thisElliptic % source(iS,iEl) = ZERO
         ENDDO
      ENDDO
 
 END SUBROUTINE FillSource_EllipticStorage
!
!
!
 SUBROUTINE DirichletConditions_EllipticStorage( thisElliptic )
 ! S/R DirichletConditions
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS( EllipticStorage ), INTENT(inout) :: thisElliptic
   ! LOCAL
   INTEGER :: nS, nEl
   
      nS  = thisElliptic % nS
      nEl = thisElliptic % nLayers
      
      IF( thisElliptic % dirichletMask(1) == ZERO )THEN
         ! Left-most boundary
         thisElliptic % cgsol(0,1) = ZERO
      ENDIF

      IF( thisElliptic % dirichletMask(2) == ZERO )THEN
         ! Right-most boundary
         thisElliptic % cgsol(nS,nEl) = ZERO
      ENDIF

 END SUBROUTINE DirichletConditions_EllipticStorage
!
!
!
 SUBROUTINE Mask_EllipticStorage( thisElliptic, thisArray )
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
   CLASS( EllipticStorage ), INTENT(in) :: thisElliptic
   REAL(prec), INTENT(inout)            :: thisArray(0:thisElliptic % nS, 1:thisElliptic % nLayers)
   ! LOCAL
   INTEGER :: iEl, nEl, nS
 
      nEl = thisElliptic % nLayers
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

 END SUBROUTINE Mask_EllipticStorage
!
!
!
 SUBROUTINE UnMask_EllipticStorage( thisElliptic, thisArray )
 ! S/R UnMask
 !
 ! This subroutine unmasks 'thisArray' along the nodes shared by two elements.
 ! A copy of the solution retained in the "MASK" operation is used to unmask.
 !
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS( EllipticStorage ), INTENT(in) :: thisElliptic
   REAL(prec), INTENT(inout)            :: thisArray(0:thisElliptic % nS, 1:thisElliptic % nLayers)
   ! LOCAL
   INTEGER :: iEl, nS, nEl
 
      nS  = thisElliptic % nS
      nEl = thisElliptic % nLayers
      
      DO iEl = 1, (nEl - 1)
         thisArray(0, iEl+1) = thisArray(nS,iEl)
      ENDDO


 END SUBROUTINE UnMask_EllipticStorage
!
!
!
 SUBROUTINE UnMaskSolution_EllipticStorage( thisElliptic )
 ! S/R UnMaskSolution
 !
 ! This subroutine unmasks the solution attribute and applies the Dirichlet boundary conditions
 ! if appropriate.
 !
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS( EllipticStorage ), INTENT(inout) :: thisElliptic
      
      CALL thisElliptic % UnMask( thisElliptic % cgsol )
      CALL thisElliptic % DirichletConditions( )

 END SUBROUTINE UnMaskSolution_EllipticStorage
!
!
!
 SUBROUTINE GlobalSum_EllipticStorage( thisElliptic, thisArray )
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
   CLASS( EllipticStorage ), INTENT(in) :: thisElliptic
   REAL(prec), INTENT(inout)            :: thisArray(0:thisElliptic % nS, 1:thisElliptic % nLayers)
   ! LOCAL
   INTEGER    :: iEl, nEl, nS
   REAL(prec) :: temp
 
      nS  = thisElliptic % nS
      nEl = thisElliptic % nLayers
      
      DO iEl = 1, (nEl - 1)
         temp = thisArray(0, iEl+1) + thisArray(nS,iEl)
         thisArray(0, iEl+1) = temp
         thisArray(nS, iEl)  = temp
      ENDDO

 END SUBROUTINE GlobalSum_EllipticStorage
!
!
!
 SUBROUTINE ApplyOperator_EllipticStorage( thisElliptic, u, Lu )
 ! S/R ApplyOperator
 !
 ! This subroutine loops over all of the elements an applies the internal matrix operations
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS( EllipticStorage ), INTENT(in) :: thisElliptic
   REAL(prec), INTENT(in)               :: u(0:thisElliptic % nS, 1:thisElliptic % nLayers)
   REAL(prec), INTENT(out)              :: Lu(0:thisElliptic % nS, 1:thisElliptic % nLayers) 
   ! LOCAL
   INTEGER :: iEl, nS, nEl
   REAL(prec) :: locSol(0:thisElliptic % nS), locDmat(0:thisElliptic % nS, 0:thisElliptic % nS)

      nS = thisElliptic % nS
      nEl = thisElliptic % nLayers

      DO iEl = 1, nEl ! loop over the elements
         locSol  = u(0:nS,iEl)
         locDmat = thisElliptic % d2Mat(0:nS,0:nS,iEl)
         Lu(0:nS,iEl) = MATMUL( locDmat, locSol )
      ENDDO ! iEl, loop over the elements
        
      
 END SUBROUTINE ApplyOperator_EllipticStorage
!
!
! 
 FUNCTION MatrixAction_EllipticStorage( this, s ) RESULT( As )
 !
 ! 
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS( EllipticStorage ) :: this
   REAL(prec)               :: s(0:this % nS, 1:this % nLayers)
   REAL(prec)               :: As(0:this % nS, 1:this % nLayers)

      CALL this % UnMask( s )

      CALL this % ApplyOperator( this % cgsol, As )

      CALL this % GlobalSum( As )

      CALL this % Mask( As )
 

 END FUNCTION MatrixAction_EllipticStorage
!
!
!
 FUNCTION Residual_EllipticStorage( this, cgStorage, s ) RESULT( r )
 ! Residual
 ! 
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS( EllipticStorage ) :: this
   TYPE( NodalStorage_1D )  :: cgStorage
   REAL(prec)               :: s(0:this % nS, 1:this % nLayers)
   REAL(prec)               :: r(0:this % nS, 1:this % nLayers)
   ! Local
   REAL(prec) :: As(0:this % nS, 1:this % nLayers)
   REAL(prec) :: w(0:this % nS), J
   INTEGER    :: iEl, iS, nEl, nS

      nEl = this % nLayers
      nS  = this % nS
      CALL cgStorage % GetQuadratureWeights( w )

      As = this % MatrixAction( s )

      DO iEl = 1, nEl
         DO iS = 0, nS
            J = this % mesh % GetJacobianAtNode( iEl, iS )
            r(iS,iEl) = this % source(iS,iEl)*J*w(iS) - As(iS,iEl)
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


 END FUNCTION Residual_EllipticStorage
!
!
!
 FUNCTION VectorProduct_EllipticStorage( this, u, v ) RESULT( uDotv )
 !
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE 
   CLASS( EllipticStorage ) :: this
   REAL(prec)               :: u(0:this % nS, 1:this % nLayers)
   REAL(prec)               :: v(0:this % nS, 1:this % nLayers)
   REAL(prec)               :: uDotv
   ! LOCAL
   INTEGER :: iEl, nEl, nS

      nEl = this % nLayers
      nS  = this % nS 

      uDotv = ZERO
      CALL this % Mask( u )
      CALL this % Mask( v )

      DO iEl = 1, nEl 
         uDotv = uDotv + DOT_PRODUCT( u(0:nS,iEl), v(0:nS,iEl) )
      ENDDO

 END FUNCTION VectorProduct_EllipticStorage
!
!
!
 SUBROUTINE ConjugateGradient_EllipticStorage( this, cgStorage, maxIter, tol, resi, ioerr )
 !  S/R Solve
 !
 !  This subroutine solves the system Ax = b using the un-preconditioned conjugate gradient method.
 !  
 !
 !  On output ioerr is set to an error checking flag. 
 !  If ioerr ==  0, the method converged within the maximum number of iterations.
 !     ioerr == -1, the method did not converge within the maximum number of iterations.
 !     ioerr == -2, something that is not caught by the current construct happened.
 !  
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS( EllipticStorage ), INTENT(inout) :: this
   TYPE( NodalStorage_1D ), INTENT(in)     :: cgStorage
   INTEGER, INTENT(in)                     :: maxIter
   REAL(prec), INTENT(in)                  :: tol
   REAL(prec), INTENT(out)                 :: resi(0:maxIter)
   INTEGER, INTENT(out)                    :: ioerr
   ! LOCAL
   INTEGER    :: iter
   REAL(prec) :: r(0:this % nS, 1:this % nLayers)
   REAL(prec) :: v(0:this % nS, 1:this % nLayers)
   REAL(prec) :: z(0:this % nS, 1:this % nLayers)
   REAL(prec) :: rNorm
   REAL(prec) :: a, b, num, den, r0
   
      ioerr = -2
      
      r = this % Residual( cgStorage, this % cgsol )
      r0 = sqrt( this % VectorProduct( r, r ) ) 
      
      resi = ZERO
      resi(0) = r0
 
      IF( r0 < tol ) then
         ioerr = 0
         RETURN
      ENDIF

      DO iter = 1, maxIter ! Loop over the CG iterates

         num = this % VectorProduct( r, r )

         IF( iter == 1) THEN
            v = r
         ELSE
            b = num/den
            v = r + b*v
         ENDIF

         z = this % MatrixAction( v )
         a = num/this % VectorProduct( v, z )
         this % cgsol = this % cgsol + a*v
         r = r - a*z
         
         den = num

         rNorm = sqrt( this % VectorProduct(r,r) )
         resi(iter) = rNorm
         IF( sqrt(rNorm)/r0 < tol ) then
           EXIT
           ioerr = 0
         ENDIF

      ENDDO ! iter, loop over the CG iterates 


      IF( SQRT(rNorm) > TOL )THEN
         PRINT*, 'MODULE IterativeSolvers : ConjugateGradient failed to converge '
         PRINT*, 'Last L-2 residual : ', sqrt(rNorm)
         ioerr=-1
      ENDIF   

 END SUBROUTINE ConjugateGradient_EllipticStorage
!
!
! 
 SUBROUTINE SolveAll_Elliptic1D( this, iEl )
 ! S/R SolveAll
 !
 !    This subroutine solves the vertical elliptic equations for all of the quadrature points in 
 !    the horizontal element stack "iEl"
 ! 
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE 
   CLASS( Elliptic1D ), INTENT(inout) :: this
   INTEGER, INTENT(in)                :: iEl
   ! Local
   INTEGER ::  i, j, nS, ioerr

      nS  = this % nS

      DO j = 0, nS
         DO i = 0, nS
            CALL this % ellStorage(i,j,iEl) % ConjugateGradient( this % cgStorage, &
                                                                 this % maxIters, &
                                                                 this % tolerance, &
                                                                 this % resi, & 
                                                                 ioerr )
         ENDDO
      ENDDO

 END SUBROUTINE SolveAll_Elliptic1D
!
!
!==================================================================================================!
!-------------------------------------- File I/O Routines -----------------------------------------!
!==================================================================================================!
!
!

END MODULE Elliptic1D_Class
 
