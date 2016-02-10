MODULE CGSM2D_CLASS
!    
! May 2015
!
!    This module defines a data structure and accompanPing routines which can be used
!    to solve elliptic equations using the single-domain nodal galerkin spectral method.
!    This module is meant to be a precursor to the multidomain continuous galerkin 
!    spectral element method. 
!
!
!
!
! ========================================================================================== !

USE COMMONDATA
USE COMMONROUTINES
USE MATRIXROUTINES
USE RUN_PARAMS_CLASS

USE NODAL_STORAGE_2D_CLASS

USE GEOMETRY_BASICS
USE VECTOR_CLASS
USE MAPPEDGEOM_2D_CLASS

USE DIAGONAL_2D_PRECONDITIONER


IMPLICIT NONE

      TYPE CGSM2D
         TYPE( NODAL_STORAGE_2D )    :: cgStorage
         TYPE( MAPPEDGEOM_2D )       :: geometry
         integer                     :: boundaryFlags(1:4) 
         real(prec), allocatable     :: solution(:,:)
         real(prec), allocatable     :: nu(:,:)      ! The " diffusion coefficient "
         real(prec), allocatable     :: source(:,:)     
         TYPE( RUN_PARAMS )          :: params
         TYPE( DIAG_PRECONDITIONER ) :: diag 

         CONTAINS

         PROCEDURE :: BUILD => BUILD_CGSM2D
         PROCEDURE :: TRASH => TRASH_CGSM2D
       
         PROCEDURE :: SET_DIAG_COEFFS

         PROCEDURE :: APPLY_BOUNDARY_MASK
         PROCEDURE :: SET_BOUNDARY_CONDITION
         PROCEDURE :: CALC_GRADIENT
         PROCEDURE :: MAPPED_LAPLACIAN
         PROCEDURE :: MATRIX_ACTION
         PROCEDURE :: RESIDUAL

         PROCEDURE :: VECTOR_PRODUCT
         PROCEDURE :: CONJUGATE_GRADIENT_CGSM2D

         PROCEDURE :: WRITE_TECPLOT => WRITE_TECPLOT_CGSM2D


      END TYPE CGSM2D


      

CONTAINS

 SUBROUTINE BUILD_CGSM2D( myCGSM, curves )
 ! S/R BUILD_CGSM2D
 !
 ! 
 ! ================================================ !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS( CGSM2D ), intent(inout) :: myCGSM
   TYPE( CURVE_2D ), intent(in)   :: curves(1:4)
   ! LOCAL
   integer :: nS, nP, iS, iP

      nS = myCGSM % params % polyDeg
      nP = nS

      ! Build the nodal storage structure
      CALL myCGSM % cgStorage % BUILD( nS, nP, CG )

      ! Build the geometry
      CALL myCGSM % geometry % BUILD( curves, myCGSM % cgStorage )


      ! Default all of the boundary flags to Dirichlet
      myCGSM % boundaryFlags = DIRICHLET
      myCGSM % boundaryFlags(1) = HOMOGENEOUS_NEUMANN
      myCGSM % boundaryFlags(3) = HOMOGENEOUS_NEUMANN
      myCGSM % boundaryFlags(4) = HOMOGENEOUS_NEUMANN
      


      ! Allocate space for the solution, source term, and diffusivity coefficient
      ALLOCATE( myCGSM % solution(0:nS,0:nP), myCGSM  % nu(0:nS,0:nP), myCGSM % source(0:nS,0:nP) )

      ! Initialize values for the solution, source term, and diffusivity coefficient
      myCGSM % solution = ZERO

      CALL myCGSM % SET_BOUNDARY_CONDITION( )   

      myCGSM % nu = ONE
   
      myCGSM % source = ZERO

      CALL myCGSM % diag % BUILD( nP )


     CALL myCGSM % SET_DIAG_COEFFS( nP )

 END SUBROUTINE BUILD_CGSM2D
!
!
!
 SUBROUTINE TRASH_CGSM2D( myCGSM )
 ! S/R TRASH_CGSM2D
 !
 ! 
 ! ================================================ !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS( CGSM2D ), intent(inout) :: myCGSM



      ! Trash the nodal storage structure
      CALL myCGSM % cgStorage % TRASH( )

      ! Trash the geometry
      CALL myCGSM % geometry % TRASH( )

      ! Trash the diagonal preconditioner
      CALL myCGSM % diag % TRASH( )

      ! Deallocate space for the solution, source term, and diffusivity coefficient
      DEALLOCATE( myCGSM % solution, myCGSM  % nu, myCGSM % source )


 END SUBROUTINE TRASH_CGSM2D
!
!
!
!================================================================================!
!=======================   PROBLEM SPECIFIC ROUTINES    =========================!
!================================================================================!
!
!
!
 SUBROUTINE SET_DIAG_COEFFS( myCGSM, polyDeg )
 ! S/R SET_DIAG_COEFFS
 !
 !   Sets the coefficients for the diagonal preconditioner.
 !   An impulse response function is used to obtain a matrix column of the laplacian
 !   operator. 
 !
 ! ============================================================================= !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS( CGSM2D ), intent(inout) :: myCGSM
   integer, intent(in)            :: polyDeg
   ! LOCAL
   integer :: iS, iP
   real(prec) :: IRF(0:polyDeg, 0:polyDeg), aIRF(0:polyDeg,0:polyDeg)
 

      do iP = 0, polyDeg
         do iS = 0, polyDeg

            ! Set the impulse response function
            IRF = ZERO
            IRF(iS,iP) = ONE
            
            ! Calculate the Laplacian of the impulse response function
            CALL myCGSM % MATRIX_ACTION(IRF, aIRF )
            

            ! The diagonal contribution is the matrix action at iS,iP
            myCGSM % diag % coeffs(iS,iP) = aIRF(iS,iP)

         enddo
      enddo
 



 END SUBROUTINE SET_DIAG_COEFFS
!
!
!================================================================================!
!=======================       SUPPORT ROUTINES       ===========================!
!================================================================================!
!
!
!
 SUBROUTINE APPLY_BOUNDARY_MASK( myCGSM, u )
 ! S/R APPLY_BOUNDARY_MASK
 !
 !
 ! ============================================================================= !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS( CGSM2D ), intent(in) :: myCGSM
   real(prec), intent(inout)   :: u(0:myCGSM % cgStorage % nS, &
                                    0:myCGSM % cgStorage % nP )
   !LOCAL
   integer :: nS, nP

      nS = myCGSM % cgStorage % nS
      nP = myCGSM % cgStorage % nP
   
      if( myCGSM % boundaryFlags(1) == DIRICHLET )then ! Set the south wall to zero

         u(:,0) = ZERO

      endif

      if( myCGSM % boundaryFlags(2) == DIRICHLET )then ! Set the east wall to zero

         u(nS,:) = ZERO

      endif

      if( myCGSM % boundaryFlags(3) == DIRICHLET )then ! Set the north wall to zero

         u(:,nP) = ZERO

      endif

      if( myCGSM % boundaryFlags(4) == DIRICHLET )then ! Set the west wall to zero

         u(0,:) = ZERO

      endif


 END SUBROUTINE APPLY_BOUNDARY_MASK
!
!
!
 SUBROUTINE SET_BOUNDARY_CONDITION( myCGSM )
 ! S/R SET_BOUNDARY_CONDITION
 !
 !
 ! ============================================================================= !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS( CGSM2D ), intent(inout) :: myCGSM
   !LOCAL
   integer :: nS, nP, iS, iP
   real(prec) :: x, y, r

      nS = myCGSM % cgStorage % nS
      nP = myCGSM % cgStorage % nP
   
      if( myCGSM % boundaryFlags(1) == DIRICHLET )then ! Set the south wall to zero

         myCGSM % solution(:,0) = ZERO

      endif

      if( myCGSM % boundaryFlags(2) == DIRICHLET )then ! Set the east wall to zero

         do iP = 0, nP
       
            x = myCGSM % geometry % x(nS,iP)

            myCGSM % solution(nS,iP) = 0.5_prec*x  ! V0*x, far field potential

         enddo

      endif

      if( myCGSM % boundaryFlags(3) == DIRICHLET )then ! Set the north wall to zero

         
         myCGSM % solution(:,nP) =  ZERO
      endif

      if( myCGSM % boundaryFlags(4) == DIRICHLET )then ! Set the west wall to zero

         myCGSM % solution(0,:) =  ZERO

      endif
 END SUBROUTINE SET_BOUNDARY_CONDITION
!
!
!
 SUBROUTINE MAPPED_LAPLACIAN( myCGSM, u, Lu )
 ! S/R MAPPED_LAPLACIAN
 !
 !
 ! ============================================================================= !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS( CGSM2D ), intent(in) :: myCGSM
   real(prec), intent(in)      :: u(0:myCGSM % cgStorage % nS, &
                                    0:myCGSM % cgStorage % nP )
   real(prec), intent(out)     :: Lu(0:myCGSM % cgStorage % nS, &
                                     0:myCGSM % cgStorage % nP )
   ! LOCAL
   integer    :: iS, iP, nS, nP
   real(prec) :: F1(0:myCGSM % cgStorage % nS, &
                    0:myCGSM % cgStorage % nP )
   real(prec) :: F2(0:myCGSM % cgStorage % nS, &
                    0:myCGSM % cgStorage % nP )
   real(prec) :: dF1ds(0:myCGSM % cgStorage % nS, &
                       0:myCGSM % cgStorage % nP )
   real(prec) :: dF2dp(0:myCGSM % cgStorage % nS, &
                       0:myCGSM % cgStorage % nP )
   real(prec) :: A, B, C, D
   real(prec) :: duds(0:myCGSM % cgStorage % nS, &
                       0:myCGSM % cgStorage % nP )
   real(prec) ::  dudp(0:myCGSM % cgStorage % nS, &
                       0:myCGSM % cgStorage % nP )

   nS = myCGSM % cgStorage % nS
   nP = myCGSM % cgStorage % nP

   ! First calculate the contravariant fluxes
      

      do iP = 0, nP ! Loop over the second computational direction

         ! Calculate the derivative in the first computational direction
         duds(0:nS,iP) = MATVECMUL( myCGSM % cgStorage % dMatX(0:nS,0:nS), &
                                    u(0:nS,iP), nS, nS )

      enddo ! iP, loop over the second computational direction



      do iS = 0, nS ! Loop over the first computational direction

         ! Calculate the derivative in the second computational direction
         dudp(iS,0:nP) = MATVECMUL( myCGSM % cgStorage % dMatY(0:nP,0:nP), &
                                    u(iS,0:nP), nP, nP )

      enddo ! iS, loop over the first computational direction

      ! Now apply the metric terms
      do iP = 0, nP ! Loop over the second computational direction
         do iS = 0, nS ! Loop over the first computational direction

   
            A = myCGSM % geometry % dydp(iS,iP)**2 + myCGSM % geometry % dxdp(iS,iP)**2 

            B = myCGSM % geometry % dydp(iS,iP)*myCGSM % geometry % dyds(iS,iP) + &
                myCGSM % geometry % dxdp(iS,iP)*myCGSM % geometry % dxds(iS,iP) 

            C = myCGSM % geometry % dyds(iS,iP)**2 + myCGSM % geometry % dxds(iS,iP)**2 

            D = myCGSM % nu(iS,iP)/myCGSM % geometry % J(iS,iP)

            ! Contravariant flux calculation
            F1(iS,iP) = D*( A*duds(iS,iP) - B*dudp(iS,iP) )*myCGSM % cgStorage % qWeightX(iS)

            F2(iS,iP) = D*( C*dudp(iS,iP) - B*duds(iS,iP) )*myCGSM % cgStorage % qWeightY(iP)
        

         enddo ! iS
      enddo ! iP


   ! Now calculate the divergence of the flux

   do iP = 0, nP ! loop over the second computational direction

      dF1ds(0:nS,iP) = MATVECMUL_TRANSPOSE( myCGSM % cgStorage % dMatX(0:nS,0:nS), F1(0:nS,iP), nS )  

   enddo ! iP

   do iS = 0, nS ! loop over the second computational direction

      dF2dp(iS,0:nP) = MATVECMUL_TRANSPOSE( myCGSM % cgStorage % dMatY(0:nP,0:nP), F2(iS,0:nP), nP )  

   enddo ! iP

   ! Now calculate the Laplacian

   do iP = 0, nP
      do iS = 0, nS

         Lu(iS,iP) = -( dF1ds(iS,iP)*myCGSM % cgStorage % qWeightY(iP) +&
                        dF2dp(iS,iP)*myCGSM % cgStorage % qWeightX(iS) )

      enddo ! iS
   enddo ! iP
       

 END SUBROUTINE MAPPED_LAPLACIAN 
!
!
!
 SUBROUTINE CALC_GRADIENT( myCGSM, u, F1, F2, nS, nP )
 ! S/R CALC_GRADIENT
 !
 !
 ! ============================================================================= !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS( CGSM2D ), intent(in) :: myCGSM
   real(prec), intent(in)      :: u(0:myCGSM % cgStorage % nS, &
                                    0:myCGSM % cgStorage % nP )
   integer, intent(in)         :: nS, nP
   real(prec), intent(out)     :: F1(0:nS,0:nP), F2(0:nS,0:nP)
   ! LOCAL
   integer    :: iS, iP
   real(prec) :: D
   real(prec) :: duds(0:nS,0:nP), dudp(0:nS,0:nP) 
      

      do iP = 0, nP ! Loop over the second computational direction

         ! Calculate the derivative in the first computational direction
         duds(0:nS,iP) = MATVECMUL( myCGSM % cgStorage % dMatX(0:nS,0:nS), &
                                    u(0:nS,iP), nS, nS )

      enddo ! iP, loop over the second computational direction



      do iS = 0, nS ! Loop over the first computational direction

         ! Calculate the derivative in the second computational direction
         dudp(iS,0:nP) = MATVECMUL( myCGSM % cgStorage % dMatY(0:nP,0:nP), &
                                    u(iS,0:nP), nP, nP )

      enddo ! iS, loop over the first computational direction

      ! Now apply the metric terms
      do iP = 0, nP ! Loop over the second computational direction
         do iS = 0, nS ! Loop over the first computational direction

            D = ONE/myCGSM % geometry % J(iS,iP)

            ! Contravariant flux calculation
            F1(iS,iP) = D*( myCGSM % geometry % dydp(iS,iP)*duds(iS,iP) - myCGSM % geometry % dyds(iS,iP)*dudp(iS,iP) )

            F2(iS,iP) = D*( -myCGSM % geometry % dxdp(iS,iP)*duds(iS,iP) + myCGSM % geometry % dxds(iS,iP)*dudp(iS,iP) ) 
        

         enddo ! iS
      enddo ! iP

 END SUBROUTINE CALC_GRADIENT
!
!
!
 SUBROUTINE MATRIX_ACTION( myCGSM, u, Au )
 !
 ! 
 ! ==========================================================!
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS( CGSM2D), intent(inout) :: myCGSM
   real(prec), intent(inout)       :: u(0:myCGSM % cgStorage % nS, &
                                        0:myCGSM % cgStorage % nP )
   real(prec), intent(out)         :: Au(0:myCGSM % cgStorage % nS, &
                                         0:myCGSM % cgStorage % nP )
   ! Local
   
     
     CALL myCGSM % MAPPED_LAPLACIAN( u, Au )


     CALL myCGSM % APPLY_BOUNDARY_MASK( Au )

 END SUBROUTINE MATRIX_ACTION
!
!
!
 SUBROUTINE RESIDUAL( myCGSM, r )
 !
 ! 
 ! ==========================================================!
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS( CGSM2D ), intent(inout) :: myCGSM
   real(prec), intent(out)         :: r(0:myCGSM % cgStorage % nS, &
                                        0:myCGSM % cgStorage % nP )
   ! LOCAL
   integer    :: iS, iP, nS, nP
   real(prec) :: w1, w2, J
   real(prec) :: Au(0:myCGSM % cgStorage % nS, &
                    0:myCGSM % cgStorage % nP )
   real(prec) :: u(0:myCGSM % cgStorage % nS, &
                    0:myCGSM % cgStorage % nP )
   
      nS = myCGSM % cgStorage % nS 
      nP = myCGSM % cgStorage % nP

      u = myCGSM % solution
     
      CALL myCGSM % MATRIX_ACTION( u, Au )

     do iP = 0, nP
        do iS = 0, nS

           w1 = myCGSM % cgStorage % qWeightX(iS)
           w2 = myCGSM % cgStorage % qWeightY(iP)

           J = myCGSM % geometry % J(iS,iP)

           r(iS,iP) = -Au(iS,iP) + myCGSM % source(iS,iP)*J*w1*w2
           

        enddo
     enddo

     CALL myCGSM % APPLY_BOUNDARY_MASK( r )


 END SUBROUTINE RESIDUAL
!
!
!
  SUBROUTINE VECTOR_PRODUCT( myCGSM, u, v, uDotv )
 !
 ! Computes the dot product of two "vectors" whose data are
 ! stored in the "(iS,iP)" format
 !
 ! ==========================================================!
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS( CGSM2D ), intent(in)  :: myCGSM
   real(prec), intent(inout)     :: u(0:myCGSM % cgStorage % nS, &
                                      0:myCGSM % cgStorage % nP )
   real(prec), intent(inout)     :: v(0:myCGSM % cgStorage % nS, &
                                      0:myCGSM % cgStorage % nP )
   real(prec), intent(out)       :: uDotv
   ! Local
   integer :: iS, iP, nS, nP
   
      nS = myCGSM % cgStorage % nS 
      nP = myCGSM % cgStorage % nP

     ! mask the arrays
     CALL myCGSM % APPLY_BOUNDARY_MASK( u )
     CALL myCGSM % APPLY_BOUNDARY_MASK( v )
 
     uDotv = ZERO

     ! Compute the dot product
     do iP = 0, nP

           uDotV = uDotV + DOT_PRODUCT( u(0:nS,iP), v(0:nS,iP) )

     enddo ! iEl, loop over the elements


 END SUBROUTINE VECTOR_PRODUCT
!
!
!
 SUBROUTINE CONJUGATE_GRADIENT_CGSM2D( myCGSM, resi )
 !
 ! Inverts the system implied by the routine matrix action and residual 
 ! using the conjugate gradient algorithm.
 !
 ! =============================================================================!
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS( CGSM2D ), intent(inout)    :: myCGSM
   real(prec), intent(out), optional :: resi(0:myCGSM % params % cgIterMax)
   ! LOCAL
   integer    :: nS, nP, iEl
   integer    :: iter
   integer    :: nIt
   real(prec) :: TOL
   real(prec) :: r(0:myCGSM % cgStorage % nS, &
                   0:myCGSM % cgStorage % nP )
   real(prec) :: v(0:myCGSM % cgStorage % nS, &
                   0:myCGSM % cgStorage % nP )
   real(prec) :: z(0:myCGSM % cgStorage % nS, &
                   0:myCGSM % cgStorage % nP )
   real(prec) :: Ad(0:myCGSM % cgStorage % nS, &
                    0:myCGSM % cgStorage % nP )
   real(prec) :: rNorm
   real(prec) :: a, b, num, den, r0
   
   

      nS = myCGSM % cgStorage % nS 
      nP = myCGSM % cgStorage % nP
      
      nIt = myCGSM % params % cgIterMax
      TOL = myCGSM % params % cgTol

      ! Calculate the residual
      CALL myCGSM % RESIDUAL( r )
      
      CALL myCGSM % VECTOR_PRODUCT( r, r, r0 )
      
      if( PRESENT( resi ) )then
         resi = ZERO
         resi(0) = r0
      endif

      ! Apply the preconditioner
      CALL myCGSM % diag % SOLVE( r, z )  ! 
     
      CALL myCGSM % APPLY_BOUNDARY_MASK( z )
     
      v = z ! Copy 

      CALL myCGSM % VECTOR_PRODUCT( r, z, num ) ! numerator   r (DOT) (H^(-1)r )

      do iter = 1,nIt ! Loop over the PCG iterates
 
         ! Compute Ad matrix-vector productz
         CALL myCGSM % MATRIX_ACTION( v, z ) ! z = A*v

         ! Compute the search-direction magnitude
         CALL myCGSM % VECTOR_PRODUCT( v, z, den) ! denominator

         a = num/den

         ! Update the solution guess
         myCGSM % solution(0:nS,0:nP) = myCGSM % solution(0:nS,0:nP) + a*v(0:nS, 0:nP)

          
         ! Update the residual
         r(0:nS,0:nP) = r(0:nS,0:nP) - a*z(0:nS,0:nP)
         
         CALL myCGSM % VECTOR_PRODUCT( r, r, rNorm )

         if( PRESENT( resi ) )then
            resi(iter) = sqrt(rNorm)
         endif

        ! print*, rNorm
         if( sqrt(rNorm) < TOL ) then
            EXIT
         endif   

         ! Apply the preconidionter to the residual
         CALL myCGSM % diag % SOLVE( r, z )! solves H z = r, to get z = H^(-1) r

         den = num ! r(DOT)[ H^(-1)r ] = r(DOT)d

         ! Calculate the change in the search direction
         CALL myCGSM % VECTOR_PRODUCT( r, z, num )
      
         ! Update the search direction
         b = num/den

         v(0:nS,0:nP) = z(0:nS,0:nP) + b*v(0:nS,0:nP)

         
      enddo ! iter, loop over the PCG iterates 


      if( sqrt(rNorm) > TOL ) then
         print*, 'MODULE CGSM2D_CLASS : CONJUGATE_GRADIENT failed to converge '
         print*, 'Last L-2 residual : ', sqrt(rNorm)

      endif   

 END SUBROUTINE CONJUGATE_GRADIENT_CGSM2D
!
!
!================================================================================!
!=======================       FILE I/O      ===========================!
!================================================================================!
!
!
!
 SUBROUTINE WRITE_TECPLOT_CGSM2D( myCGSM, nPlot, nOld, plotInterp, Tmat, filename )
 ! WRITE_TECPLOT_CGSM2D
 !  Description :
 ! 
 !
 !    Subroutine dependencies :
 !    (NONE)
 !    
 !
 !  Input :
 !    type(LAG_INTERP2D) :: myPoly
 !
 !    real(prec) :: fAtNodes(:,:) - the interpolation nodes function values.
 !
 !    character(*) :: filename - name of output file
 !
 !
 ! 
 !------------------------------------------------------------------------
  IMPLICIT NONE
  CLASS( CGSM2D ), intent(in)           :: myCGSM
  integer, intent(in)                   :: nOld, nPlot
  real(prec), intent(in)                :: Tmat(0:nPlot, 0:nOld)
  TYPE( LAG_INTERP2D ), intent(in)      :: plotInterp
  character(*), intent(in)              :: filename
  !LOCAL
  real(prec)  :: x(0:nPlot,0:nPlot), y(0:nPlot,0:nPlot)!, depth(0:nPlot,0:nPlot)
  real(prec)  :: u(0:nPlot,0:nPlot)
  real(prec)  :: J(0:nPlot,0:nPlot)
  integer :: iX, iY, iZ, iEl
  character(len=5) :: zoneID


    open( unit=2, file= trim(filename)//'.tec', form='formatted',status='replace')

    write(2,*) 'VARIABLES = "X", "Y", "Potential", "Jac"'
    


       ! Interpolate the solutions onto a uniform mesh
       CALL COARSE_TO_FINE_2D(myCGSM % cgStorage % interp,&
                                      myCGSM % geometry % x,&
                                      plotInterp, x, Tmat, Tmat)


       CALL COARSE_TO_FINE_2D(myCGSM % cgStorage % interp,&
                                      myCGSM %  geometry % y,&
                                      plotInterp, y, Tmat, Tmat)

       CALL COARSE_TO_FINE_2D(myCGSM % cgStorage % interp,&
                                      myCGSM % solution,&
                                      plotInterp, u, Tmat, Tmat)
 
        CALL COARSE_TO_FINE_2D(myCGSM % cgStorage % interp,&
                                      myCGSM % geometry % J,&
                                      plotInterp, J, Tmat, Tmat)

        write(zoneID,'(I5.5)') 0


        write(2,*) 'ZONE T="el'//trim(zoneID)//'", I=',nPlot+1,', J=', nPlot+1,',F=POINT'

        do iY = 0, nPlot
           do iX = 0, nPlot

              write (2,*)  x( iX, iY ), y( iX, iY ), u(iX,iY), J(iX,iY) 

           enddo
        enddo


    close(unit=2)

 END SUBROUTINE WRITE_TECPLOT_CGSM2D

END MODULE CGSM2D_CLASS
