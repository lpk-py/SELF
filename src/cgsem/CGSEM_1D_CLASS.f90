MODULE CGSEM_1D_CLASS
! CGSEM_1D_CLASS.f90
! 
!    This module defines the data structure used for implementing the Continuous Galerkin
!    Spectral Element Method. All data is stored globally, instead of element by element.
!    This facilitates the use of LAPACK and BLAS routines for matrix operations.
!
!    This version of the CGSEM_1D_CLASS does not implement over-integration. Instead the 
!    interpolation nodes are equivalent to the quadrature nodes which can result in aliasing
!    errors; though these errors are spectrally accurate, in some problems this may require
!    high order methods to reduce them sufficiently.
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
! Module History 
! 
! o  April 2015
!    
!
! 
! List of Module Procedures
!
!
!
!
!  
! =======================================================================================

USE COMMONDATA 
USE MATRIXROUTINES

USE LEGENDRE
USE LAGRANGE_1D_CLASS

USE NODAL_STORAGE_1D_CLASS

USE MESH_1D_CLASS
USE ELEMENT_1D_CLASS
USE MAPPEDGEOM_1D_CLASS

USE CGSEM_1D_FEM_PRECONDITIONER


IMPLICIT NONE

   TYPE SHARED_NODE_PTR
      integer :: eLeft, eRight ! element ID of the elements to the left and right
      integer :: nodeLeft, nodeRight ! nodeLeft is the local node ID of element(eLeft) that is shared with the current element
                                     ! nodeRight is the local node ID of element(eRight) that is shared with the current element
   END TYPE SHARED_NODE_PTR
   
   TYPE CGSEM_1D
      integer                            :: N, K ! polynomial degree, and number of elements
      TYPE(NODAL_STORAGE_1D)             :: nodalStorage ! Contains solution interpolant and first derivative matrix
      TYPE(MESH_1D)                      :: mesh
      TYPE(SHARED_NODE_PTR), allocatable :: p(:)
      real(prec), allocatable            :: d2Mat(:,:,:) ! matrix for d/dx( k(x) d/dx ) operator
      real(prec), allocatable            :: sol(:,:)   ! Solution storage (0:N,1:K)
      real(prec), allocatable            :: RHS(:,:)   ! Tendency for implicit solvers (0:N,1:K)

      CONTAINS

      ! Manual Constructors/Destructors
      PROCEDURE :: BUILD => BUILD_CGSEM_1D
      PROCEDURE :: TRASH => TRASH_CGSEM_1D

      ! Problem specific
      PROCEDURE :: SET_BOUNDARY_CONDITIONS 
      PROCEDURE :: TRAPEZOIDAL_INTEGRATION => TRAPEZOIDAL_INTEGRATION_CGSEM1D

      ! Global procedure
      PROCEDURE :: MASK
      PROCEDURE :: UNMASK
      PROCEDURE :: GLOBAL_SUM
      PROCEDURE :: VECTOR_PRODUCT

      ! Local procedure
      PROCEDURE :: LAPLACIAN
      PROCEDURE :: MATRIX_ACTION
      PROCEDURE :: RESIDUAL

      ! Iterative solvers
      PROCEDURE :: CONJUGATE_GRADIENT => CONJUGATE_GRADIENT_CGSEM1D
 
      ! File I/O
      PROCEDURE :: WRITE_TECPLOT => WRITE_TECPLOT_CGSEM1D


   END TYPE CGSEM_1D


   CONTAINS


 SUBROUTINE BUILD_CGSEM_1D( myCGSEM, N, K, x )
 ! S/R BUILD_CGSEM_1D
 !
 ! Build the CGSEM_1D data structure with polynomial degree N and with K elements
 ! x(0:K) defines the corner nodes of the elements
 !
 ! ============================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS( CGSEM_1D ), intent(out) :: myCGSEM
   integer, intent(in)            :: N, K
   real(prec), intent(in)         :: x(0:K)
   ! LOCAL
   integer    :: iEl, row, col, m
   real(prec) :: lmr, lmc, w, J, sK, thisX

      myCGSEM % N = N
      myCGSEM % K = K

      CALL myCGSEM % nodalStorage % BUILD( N, CG )

      CALL myCGSEM % mesh % BUILD( myCGSEM % nodalStorage, x, K )


      ALLOCATE( myCGSEM % p(1:K), myCGSEM % d2Mat(0:N,0:N,1:K) )

      ALLOCATE( myCGSEM % sol(0:N,1:K), myCGSEM % RHS(0:N,1:K) )

      myCGSEM % sol = ZERO ! Initialize

      myCGSEM % RHS = ZERO ! Initialize

     ! Build the node pointers

     do iEl = 1, K

        myCGSEM % p(iEl) % eLeft = iEl
        myCGSEM % p(iEl) % nodeLeft = N

        myCGSEM % p(iEl) % eRight = iEl + 1
        myCGSEM % p(iEl) % nodeRight = 0


     enddo ! iEl, loop over the elements

     ! Set up the operator d/dx( K(x) d/dx). Relies on specification of k(x) in FUNCTION STIFFNESS_WEIGHT
      
     do iEl = 1,K ! Loop over the elements  
 
        do col = 0, N ! Loop over the columns of the matrix

           do row = 0, N ! Loop over the rows of the matrix

              myCGSEM % d2Mat(row,col,iEl) = ZERO ! Initialize the sum to zero

              do m = 0, N ! Loop over the quadrature points (here, the same as the interpolation points)

              
                 thisX = myCGSEM % mesh % elements(iEl) % geometry % x(m) ! Get the physical location

                 sK = STIFFNESS_WEIGHT( thisX ) ! calculate the stiffness weight, "K" in d/dx( K(x) d/dx )

                 J = myCGSEM % mesh % elements(iEl) % geometry % J(m) ! Get the Jacobian of the coordinate transformation

                 w = myCGSEM % nodalStorage % qWeight(m) ! Get the quadrature weight

                 lmr = myCGSEM % nodalStorage % dMat(m,row) ! Derivative matrix

                 lmc = myCGSEM % nodalStorage % dMat(m,col) ! Derivative matrix

                 myCGSEM % d2Mat(row,col,iEl) = myCGSEM % d2Mat(row,col,iEl) - (sK*w/J)*lmr*lmc


              enddo ! m, loop over the quadrature points

           enddo ! row
 
        enddo ! col

     enddo ! iEl

 END SUBROUTINE BUILD_CGSEM_1D
!
!
!
 SUBROUTINE TRASH_CGSEM_1D( myCGSEM )
 ! S/R BUILD_CGSEM_1D
 !
 ! Build the CGSEM_1D data structure with polynomial degree N and with K elements
 ! x(0:K) defines the corner nodes of the elements
 !
 ! ============================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS( CGSEM_1D ), intent(inout) :: myCGSEM
   ! LOCAL
   integer :: iEl


      CALL myCGSEM % nodalStorage % TRASH( )

      CALL myCGSEM % mesh % TRASH( )


      DEALLOCATE( myCGSEM % p, myCGSEM % d2Mat )

      DEALLOCATE( myCGSEM % sol, myCGSEM % RHS )




   
 END SUBROUTINE TRASH_CGSEM_1D
!
!
!
!================================================================================!
!=======================       PROBLEM SPECIFIC       ===========================!
!================================================================================!
! This includes specification of the time integrator, stiffness function, and
! additional sources.
!
!
 SUBROUTINE TRAPEZOIDAL_INTEGRATION_CGSEM1D( myCGSEM, t, dt, nIterMax, cgTol, myFEM )
 ! S/R TRAPEZOIDAL_INTEGRATION_CGSEM1D
 !
 ! Uses the trapezoidal rule to integrate the heat equation in time.
 ! This subroutine depends on the preconditioned conjugate gradient 
 ! method subroutine provided below (S/R CONJUGATE_GRADIENT_SEM1D)
 !
 ! =============================================================================!
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS( CGSEM_1D ), intent(inout)        :: myCGSEM
   real(prec), intent(in)                  :: t, dt    ! The time at the beginning of this time interval and the time step size
   integer, intent(in)                     :: nIterMax ! The maximum number of iterations for conjugate gradient method
   real(prec), intent(in)                  :: cgTol    ! The residual tolerance for the conjugate gradient method
   TYPE( FEM1D_PRECONDITIONER), intent(in) :: myFEM    ! The finite element preconditioner data structure used for the conjugate gradient method

   
      CALL myCGSEM % MATRIX_ACTION( ONE, dt, myCGSEM % sol, myCGSEM % RHS ) ! Compute the right hand side of the discrete heat equation

      CALL myCGSEM % SET_BOUNDARY_CONDITIONS( t + dt) ! Fill in the boundary conditions at the next time level for the implicit side
 
      CALL myCGSEM % CONJUGATE_GRADIENT( dt, nIterMax, cgTol, myFEM ) ! Invert the system with the Conjugate Gradient method

      CALL myCGSEM % UNMASK( myCGSEM % sol )
 END SUBROUTINE TRAPEZOIDAL_INTEGRATION_CGSEM1D
!
!
!
 FUNCTION STIFFNESS_WEIGHT( x ) RESULT( K )
 ! 
 ! This function returns the "K" used in the computation of 
 ! the operator  " d/dx( K(x) d/dx ) ."
 !
 !
 ! ==============================================================================!
 ! DECLARATION
   IMPLICIT NONE
   real(prec) :: x,K


      K = ONE

 END FUNCTION STIFFNESS_WEIGHT
!
!
!
 SUBROUTINE SET_BOUNDARY_CONDITIONS( myCGSEM, t )
 !
 ! The subroutine does as the name suggests - it sets the boundary conditions.
 ! Note that this is only appropriate for Dirichlet boundary conditions.
 ! The time "t" is passed in the anticipation that one may have time-varying
 ! boundary conditions at some point.
 !
 ! ==============================================================================!
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS( CGSEM_1D ), intent(inout) :: myCGSEM
   real(prec), intent(in)           :: t

      ! Left-most boundary
      myCGSEM % sol(0, 1) = ZERO

      ! Right-most boundary
      myCGSEM % sol(myCGSEM % N, myCGSEM % K) = ZERO

 END SUBROUTINE SET_BOUNDARY_CONDITIONS
!
!
!
!================================================================================!
!=======================      GLOBAL OPERATIONS       ===========================!
!================================================================================!
!
!
!
 SUBROUTINE MASK( myCGSEM, thisArray )
 ! S/R MASK 
 !
 ! This subroutine masks out 'thisArray' along the shared element boundaries. 
 ! Arbitrarily, we choose to mask the element to the right (in 1-D) of
 ! the shared node.
 !
 ! For local operations, an "UNMASK" routine is needed - this is included 
 ! below.
 !
 ! ========================================================================================= !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS( CGSEM_1D ), intent(in) :: myCGSEM
   real(prec), intent(inout)     :: thisArray(0:myCGSEM % N, 1:myCGSEM % K)
   ! LOCAL
   integer :: iEl, eRight, nRight
 
      do iEl = 1, (myCGSEM % K - 1) ! Loop over all of the elements except the last one
  
         eRight = myCGSEM % p(iEl) % eRight ! Get  the element ID to the right of this edge

         nRight = myCGSEM % p(iEl) % nodeRight ! Get the local node ID of the shared global node in that element

         
         ! Mask out that duplicate entry in the array
         thisArray(nRight, eRight) = ZERO
    
      enddo ! iEl, loop over the elements


 END SUBROUTINE MASK
!
!
!
 SUBROUTINE UNMASK( myCGSEM, thisArray )
 ! S/R UNMASK 
 !
 ! This subroutine unmasks 'thisArray' along the nodes shared by two elements.
 ! A copy of the solution retained in the "MASK" operation is used to unmask.
 !
 !
 ! ========================================================================================= !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS( CGSEM_1D ), intent(in) :: myCGSEM
   real(prec), intent(inout)     :: thisArray(0:myCGSEM % N, 1:myCGSEM % K)
   ! LOCAL
   integer :: iEl, eRight, nRight, eLeft, nLeft
 
      do iEl = 1, (myCGSEM % K - 1) ! Loop over all of the elements except the last one
  
         eLeft = myCGSEM % p(iEl) % eLeft ! Get  the element ID to the left of this edge

         nLeft = myCGSEM % p(iEl) % nodeLeft ! Get the local node ID of the shared global node in that element

         eRight = myCGSEM % p(iEl) % eRight ! Get  the element ID to the right of this edge

         nRight = myCGSEM % p(iEl) % nodeRight ! Get the local node ID of the shared global node in that element

         ! Copy the solution into the masked out solution
         thisArray(nRight, eRight) = thisArray(nLeft,eLeft)
    
      enddo ! iEl, loop over the elements


 END SUBROUTINE UNMASK
!
!
!
 SUBROUTINE GLOBAL_SUM( myCGSEM, thisArray )
 ! S/R GLOBAL_SUM 
 !
 ! In the continuous galerkin spectral element method, the equation at the nodes shared 
 ! by two elements involve the sum of the equations implied by each element. This routine
 ! computes the sum of 'thisArray' at the shared nodes.
 !
 !
 ! ========================================================================================= !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS( CGSEM_1D ), intent(in) :: myCGSEM
   real(prec), intent(inout)     :: thisArray(0:myCGSEM % N, 1:myCGSEM % K)
   ! LOCAL
   integer    :: iEl, eRight, nRight, eLeft, nLeft
   real(prec) :: temp
 
      do iEl = 1, (myCGSEM % K - 1) ! Loop over all of the elements except the last one
  
         eLeft = myCGSEM % p(iEl) % eLeft ! Get  the element ID to the left of this edge

         nLeft = myCGSEM % p(iEl) % nodeLeft ! Get the local node ID of the shared global node in that element

         eRight = myCGSEM % p(iEl) % eRight ! Get  the element ID to the right of this edge

         nRight = myCGSEM % p(iEl) % nodeRight ! Get the local node ID of the shared global node in that element

         ! Add the contributions from both elements
         temp = thisArray(nRight, eRight) + thisArray(nLeft,eLeft)

         ! Store the resulting sum in each copy
         thisArray(nRight, eRight) = temp

         thisArray(nLeft, eLeft) = temp
         
    
      enddo ! iEl, loop over the elements


 END SUBROUTINE GLOBAL_SUM
!
!
!
!================================================================================!
!=======================      MATRIX OPERATIONS       ===========================!
!================================================================================!
!
!
!
 SUBROUTINE LAPLACIAN( myCGSEM, u, Lu )
 ! S/R LAPLACIAN
 !
 ! This subroutine uses computes the matrix action of the weighted laplacian :
 !  Lu = d/dx( K(x) du/dx ) where K(x) is set in the FUNCTION STIFFNESS_WEIGHT
 !
 ! ========================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS( CGSEM_1D ), intent(in) :: myCGSEM
   real(prec), intent(in)        :: u(0:myCGSEM % N, 1:myCGSEM % K)
   real(prec), intent(out)       :: Lu(0:myCGSEM % N, 1:myCGSEM % K) 
   ! LOCAL
   integer :: iEl, N, K

      N = myCGSEM % N
      K = myCGSEM % K

      do iEl = 1, K ! loop over the elements

         Lu(0:N,iEl) = MATVECMUL( myCGSEM % d2Mat(0:N,0:N,iEl), u(0:N,iEl), N, N )

      enddo ! iEl, loop over the elements
        
 END SUBROUTINE LAPLACIAN
!
!
! 
 SUBROUTINE MATRIX_ACTION( myCGSEM, s, dt, u, Au )
 !
 ! Computes the matrix action of the implicit matrix for the
 ! heat-equation using the implicit trapezoidal rule.
 !
 !
 !  [ J*w - dt*(G,)/2 ]U^{n+1} = [J*w + dt*(G,)/2]U^{n}
 !
 !  For the left-hand-side, s < 0 
 !  For the right-hand-side, s > 0
 ! ==========================================================!
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS( CGSEM_1D ), intent(inout) :: myCGSEM
   real(prec), intent(in)           :: dt, s
   real(prec), intent(inout)        :: u(0:myCGSEM % N, 1:myCGSEM % K)
   real(prec), intent(out)          :: Au(0:myCGSEM % N, 1:myCGSEM % K)
   ! Local
   integer :: iEl, iNode, K, N
   real(prec) :: w, J
   
      K = myCGSEM % K 
      N = myCGSEM % N

     ! unmask the solution array
     CALL myCGSEM % UNMASK( u ) 

     ! Apply the operator d/dx( K(x) d/dx )U
     CALL myCGSEM % LAPLACIAN( u, Au )

     ! Compute the matrix action for the left-hand-side of the trapezoidal rule
     do iEl = 1, K ! Loop over the elements
        do iNode = 0, N ! Loop over the nodes

           w = myCGSEM % nodalStorage % qWeight(iNode) ! Get the quadrature weight

           J = myCGSEM % mesh % elements(iEl) % geometry % J(iNode) ! Get the Jacobian of the coordinate transformation

           Au(iNode, iEl) = J*w*u(iNode, iEl) + HALF*s*dt*Au(iNode, iEl)

        enddo
     enddo ! iEl, loop over the elements

     ! Add the contributions to the edges together
     CALL myCGSEM % GLOBAL_SUM( Au )

     CALL myCGSEM  % MASK( Au )
     ! Mask the solution
     CALL myCGSEM % MASK( u )

     ! For the left-hand-side only, the homogeneous Dirichlet Boundary conditions should have no influence
     if( s < ZERO )then

        ! Mask out the copied contributions from the shared edges
     !   CALL myCGSEM  % MASK( Au )
          
        ! Zero out the matrix action where the Dirichlet boundary conditions are applied
        Au(0, 1) = ZERO
        Au(N, K) = ZERO

     endif

 END SUBROUTINE MATRIX_ACTION
!
!
!
 SUBROUTINE RESIDUAL( myCGSEM, dt, r  )
 !
 ! To compute the solution to the system
 !
 !  [ J*w - dt*(G,) ]U^{n+1} = [J*w + dt*(G,)]U^{n}
 !                 AU = B
 !
 ! for U^{n+1}, the pre-conditioned conjugate gradient algorithm 
 ! is employed (below). For the PCG algorithm, it is necessary to 
 ! calculate the residual " r = B - AU  ".
 ! 
 !  
 ! ==========================================================!
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS( CGSEM_1D ), intent(inout) :: myCGSEM
   real(prec), intent(in)           :: dt
   real(prec), intent(inout)        :: r(0:myCGSEM % N, 1:myCGSEM % K)
   ! Local
   integer :: iEl, iNode, K, N
   
      K = myCGSEM % K 
      N = myCGSEM % N

     ! Calculate the left-hand-side of the system
     CALL myCGSEM % MATRIX_ACTION( -ONE, dt, myCGSEM % sol, r )

     ! Compute the residual
     do iEl = 1, K ! Loop over the elements
        do iNode = 0, N ! Loop over the nodes

           r(iNode, iEl) = myCGSEM % RHS(iNode, iEl) - r(iNode, iEl) 

        enddo
     enddo ! iEl, loop over the elements

     ! Mask the residuals
     CALL myCGSEM % MASK( r )

     ! Set the residual to zero where the Dirichlet boundary conditions are applied
     r(0, 1) = ZERO
     r(N, K) = ZERO

 END SUBROUTINE RESIDUAL
!
!
!
  SUBROUTINE VECTOR_PRODUCT( myCGSEM, u, v, uDotv )
 !
 ! Computes the dot product of two "vectors" whose data are
 ! stored in the "local node, element-ID" format
 !
 ! ==========================================================!
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS( CGSEM_1D ), intent(in) :: myCGSEM
   real(prec), intent(inout)     :: u(0:myCGSEM % N, 1:myCGSEM % K)
   real(prec), intent(inout)     :: v(0:myCGSEM % N, 1:myCGSEM % K)
   real(prec), intent(out)       :: uDotv
   ! Local
   integer :: iEl, iNode, K, N
   
      K = myCGSEM % K 
      N = myCGSEM % N

     ! mask the arrays
     CALL myCGSEM % MASK( u )
     CALL myCGSEM % MASK( v )
 
     uDotv = ZERO

     ! Compute the dot product
     do iEl = 1, K ! Loop over the elements

           uDotV = uDotV + DOT_PRODUCT( u(0:N,iEl), v(0:N,iEl) )

     enddo ! iEl, loop over the elements

     ! unmask the arrays
    ! CALL myCGSEM % UNMASK( u )

    ! CALL myCGSEM % UNMASK( v )

 END SUBROUTINE VECTOR_PRODUCT
!
!
!
!================================================================================!
!=======================       ITERATIVE SOLVER       ===========================!
!================================================================================!
!
!
!
 SUBROUTINE CONJUGATE_GRADIENT_CGSEM1D( myCGSEM, dt, nIt, TOL, myFEM )
 !
 ! Inverts the system implied by the routine matrix action and residual 
 ! using the conjugate gradient algorithm.
 !
 ! =============================================================================!
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS( CGSEM_1D ), intent(inout)          :: myCGSEM
   TYPE( FEM1D_PRECONDITIONER ), intent(in)  :: myFEM   ! Pre-conditioner
   real(prec), intent(in)                    :: dt
   integer, intent(in)                       :: nIt
   real(prec), intent(in)                    :: TOL
   ! LOCAL
   integer    :: nEl, nDeg, iEl
   integer    :: iter
   real(prec) :: r(0:myCGSEM % N, 1:myCGSEM % K)
   real(prec) :: v(0:myCGSEM % N, 1:myCGSEM % K)
   real(prec) :: z(0:myCGSEM % N, 1:myCGSEM % K)
   real(prec) :: Ad(0:myCGSEM % N, 1:myCGSEM % K)
   real(prec) :: rNorm
   real(prec) :: a, b, num, den, r0
   

      nEl = myCGSEM % K
      nDeg = myCGSEM % N
      
      ! Calculate the residual
      CALL myCGSEM % RESIDUAL( dt, r )
      
      CALL myCGSEM % VECTOR_PRODUCT( r, r, r0 )

      ! Apply the preconditioner
      CALL myFEM % SOLVE( r, z, 1.5_prec )  ! Trying SORfac of 1.5 (should be between 1.0 and 2  for convergence )
      
      v = z ! Copy 

      CALL myCGSEM % VECTOR_PRODUCT( r, z, num ) ! numerator   r (DOT) (H^(-1)r )

      do iter = 1,nIt ! Loop over the PCG iterates
   
         ! Compute Ad matrix-vector product
         CALL myCGSEM % MATRIX_ACTION( -ONE, dt, v, z ) ! z = A*v
         
         ! Compute the search-direction magnitude
         CALL myCGSEM % VECTOR_PRODUCT( v, z, den) ! denominator

         a = num/den

         ! Update the solution guess
         myCGSEM % sol(0:nDeg, 1:nEl) = myCGSEM % sol(0:nDeg, 1:nEl) + a*v(0:nDeg, 1:nEl)

          
         ! Update the residual
         r(0:nDeg, 1:nEl) = r(0:nDeg, 1:nEl) - a*z(0:nDeg, 1:nEl)
         
         CALL myCGSEM % VECTOR_PRODUCT( r, r, rNorm )
    
       ! print*, iter, MAXVAL(rmaxes)
         if( rNorm/r0 < TOL ) then
            EXIT
         endif   

         ! Apply the preconidionter to the residual
         CALL myFEM % SOLVE( r, z, 1.5_prec ) ! solves H d2 = r, to get d2 = H^(-1) r

         den = num ! r(DOT)[ H^(-1)r ] = r(DOT)d

         ! Calculate the change in the search direction
         CALL myCGSEM % VECTOR_PRODUCT( r, z, num )
      
         ! Update the search direction
         b = num/den

         v(0:nDeg,1:nEl) = z(0:nDeg,1:nEl) + b*v(0:nDeg,1:nEl)

         
      enddo ! iter, loop over the PCG iterates 


      if( rNorm > TOL ) then
         print*, 'MODULE CGSEM_1D_CLASS : CONJUGATE_GRADIENT failed to converge '
         print*, 'Last L-2 residual : ', sqrt(rNorm)/sqrt(r0)

         !do iEl = 1, nEl
         !   print*, r(:,iEl)
         !enddo
      endif   

 END SUBROUTINE CONJUGATE_GRADIENT_CGSEM1D
!
!
!
!================================================================================!
!=======================           FILE I/O           ===========================!
!================================================================================!
!
!
!
 SUBROUTINE WRITE_TECPLOT_CGSEM1D( myCGsem, Tx, nPlot, plotInterp, filename )
 !
 !
 ! 
 ! ==============================================================================!
 ! DECLARATIONS
   CLASS( CGSEM_1D ), intent(in)    :: myCGSEM
   integer, intent(in)              :: nPlot
   real(prec), intent(in)           :: Tx(0:nPlot,0:myCGSEM % N)
   TYPE( LAG_INTERP1D ), intent(in) :: plotInterp
   character(*), intent(in)         :: filename
   ! LOCAL
   real(prec)   :: x(0:nPlot), p(0:nPlot)
   integer      :: iEl, iX
   character(5) :: zoneID

    open( unit=2, file= trim(filename)//'.curve', form='formatted',status='replace')

    write(2,*) '#heat'
    

    do iEl = 1, myCGSEM % K

       ! Interpolate the solutions onto a uniform mesh
       CALL COARSE_TO_FINE_1D(myCGsem % nodalStorage % interp,&
                                      myCGsem % mesh % elements(iEl) % geometry % x,&
                                      plotInterp, x, Tx )


       CALL COARSE_TO_FINE_1D(myCGsem % nodalStorage % interp,&
                                      myCGsem % sol(:,iEl),&
                                      plotInterp, p, Tx )
 
     

     !   write(zoneID,'(I5.5)') iEl


 !       write(2,*) 'ZONE T="el'//trim(zoneID)//'", I=',nPlot+1,',F=POINT'

        do iX = 0, nPlot-1

           write(2,'(F15.7,2x,F15.7)') x(iX), p(iX)

        enddo
        

    enddo

    write(2,'(F15.7,2x,F15.7)') x(nPlot), p(nPlot)

    close(unit=2)


 END SUBROUTINE WRITE_TECPLOT_CGSEM1D
!
!
!
END MODULE CGSEM_1D_CLASS
