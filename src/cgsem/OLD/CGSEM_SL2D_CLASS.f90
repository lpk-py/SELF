MODULE CGSEM_SL2D_CLASS
!    
! July 2015
!
!    This module defines a data structure and accompanPing routines which can be used
!    to solve elliptic equations using the multi-domain nodal continuous galerkin spectral element method (CGSEM).
!
!
!
!
! ========================================================================================== !

! The model settings
USE MODEL_PRECISION
USE CONSTANTS_DICTIONARY
USE MODEL_FLAGS
!
! Basic Routines
USE COMMONROUTINES
USE MATRIXROUTINES
USE RUN_PARAMS_CLASS
!
! Storage container modules
USE NODAL_STORAGE_2D_CLASS
!
! Geometry routines
USE QUADMESH_CLASS
USE EDGE_CLASS
USE CORNERNODE_2D_CLASS
USE GEOMETRY_BASICS
USE VECTOR_CLASS
USE MAPPEDGEOM_2D_CLASS
!
!
USE CGSEM2D_CLASS
USE DIAGONAL_2DSEM_PRECONDITIONER


IMPLICIT NONE

      TYPE, EXTENDS( CGSEM2D ) :: CGSEM_SL2D
      !   ****** INHERITED FROM CGSEM2D ****** !
      !   TYPE( NODAL_STORAGE_2D )        :: cgStorage
      !   TYPE( QUADMESH )                :: mesh           ! Boundary condition flags are stored in the secondary element ID's of the edges.
      !   real(prec), allocatable         :: solution(:,:,:)
      !   real(prec), allocatable         :: source(:,:,:)     
      !   TYPE( RUN_PARAMS )              :: params
      !   TYPE( DIAG_SEM_PRECONDITIONER ) :: diag 
          
         integer                            :: nFree      ! number of degrees of freedom
         integer, allocatable               :: arrayMap(:,:) ! (1:nFree,1:3), maps from the array to the solution mesh.
         real(prec), allocatable            :: eVals(:,:)   ! The eigenvalues
         real(prec), allocatable            :: eModes(:,:,:,:) ! eigenvectors dimension (0:nS,0:nP,1:nEl,1:nFree)

         CONTAINS

         PROCEDURE :: BUILD => BUILD_CGSEM_SL2D
         PROCEDURE :: TRASH => TRASH_CGSEM_SL2D
       
       !  PROCEDURE :: SET_DIAG_COEFFS ( inherited )

         PROCEDURE :: COUNT_DEGREES_OF_FREEDOM
         PROCEDURE :: GET_ARRAY_MAP
         PROCEDURE :: BUILD_MATRIX_FROM_ACTION
         PROCEDURE :: BUILD_WEIGHT_MATRIX

         PROCEDURE :: CALC_EIGENMODES

         PROCEDURE :: ARRAY_TO_SOLUTION
         PROCEDURE :: WRITE_TECPLOT_MODES => WRITE_TECPLOT_CGSEM_SL2D
         PROCEDURE :: WRITE_EVALS_TO_CURVE
         
       !  PROCEDURE :: MASK => MASK_CGSEM2D
       !  PROCEDURE :: UNMASK => UNMASK_CGSEM2D
       !  PROCEDURE :: GLOBAL_SUM => GLOBALSUM_CGSEM2D

       !  PROCEDURE :: MATRIX_ACTION
       !  PROCEDURE :: RESIDUAL

       !  PROCEDURE :: VECTOR_PRODUCT
       !  PROCEDURE :: CONJUGATE_GRADIENT_CGSEM2D

       !  PROCEDURE :: WRITE_TECPLOT => WRITE_TECPLOT_CGSEM2D
        

      END TYPE CGSEM_SL2D


      

CONTAINS

 SUBROUTINE BUILD_CGSEM_SL2D( myCGSEM )
 ! S/R BUILD_CGSEM_SL2D
 !
 ! 
 ! ================================================ !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS( CGSEM_SL2D ), intent(inout) :: myCGSEM
   ! LOCAL
   integer :: nS, nP, iS, iP, nF
   character(40) :: meshFile
   character(10) :: pickupIterChar
   real(prec), allocatable :: dummyMask(:,:,:)

      CALL myCGSEM % params % BUILD( )

      nS = myCGSEM % params % polyDeg
      nP = nS

      ! Build the nodal storage structure
      CALL myCGSEM % cgStorage % BUILD( nS, nP, GAUSS_LOBATTO, CG )

      ! Build the geometry
      meshFile = 'mesh'

     ! Builds the lateral mesh
      if( TRIM(myCGSEM % params % ISMmeshFile) == ' ' )then
 
         CALL myCGSEM % mesh % READ_PICKUP( trim(meshfile), myCGSEM % cgStorage  ) 

      else

         CALL myCGSEM % mesh % READ_SPECMESH_FILE( myCGSEM % params, myCGSEM % cgStorage )

      endif

      ! Allocate space for the solution, source term, and diffusivity coefficient
      ALLOCATE( myCGSEM % solution(0:nS, 0:nP, 1:myCGSEM % mesh % nElems), &
                myCGSEM % source(0:nS, 0:nP, 1:myCGSEM % mesh % nElems), &
                dummyMask(0:nS, 0:nP, 1:myCGSEM % mesh % nElems) )

        

      ! Initialize values for the solution and source term
    
      myCGSEM % solution = ZERO

      myCGSEM % source = ZERO
     

      CALL myCGSEM % diag % BUILD( nP, myCGSEM % mesh % nElems )

      CALL myCGSEM % SET_DIAG_COEFFS( nP )


      ! Setting up the Sturm-Liouville matrices
      CALL myCGSEM % COUNT_DEGREES_OF_FREEDOM( dummyMask )

      nF = myCGSEM % nFree

      ! Get the array map
      ALLOCATE( myCGSEM % arrayMap(1:nF,1:3) )
      
      CALL myCGSEM % GET_ARRAY_MAP( dummyMask )
 

      ALLOCATE( myCGSEM % eVals(1:nF,1:3), myCGSEM % eModes(0:nS, 0:nP, 1:myCGSEM % mesh % nElems,1:nF) )

      

 END SUBROUTINE BUILD_CGSEM_SL2D
!
!
!
 SUBROUTINE TRASH_CGSEM_SL2D( myCGSEM )
 ! S/R TRASH_CGSEM2D
 !
 ! 
 ! ================================================ !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS( CGSEM_SL2D ), intent(inout) :: myCGSEM



      ! Trash the nodal storage structure
      CALL myCGSEM % cgStorage % TRASH( )

      ! Trash the geometry
      CALL myCGSEM % mesh % TRASH( )

      ! Trash the diagonal preconditioner
      CALL myCGSEM % diag % TRASH( )

      ! Deallocate space for the solution, source term, and diffusivity coefficient
      DEALLOCATE( myCGSEM % solution, myCGSEM % source )
      DEALLOCATE( myCGSEM % eVals, myCGSEM % eModes )
      DEALLOCATE( myCGSEM % arrayMap )
 


 END SUBROUTINE TRASH_CGSEM_SL2D
!
!
!
!================================================================================!
!=======================   PROBLEM SPECIFIC ROUTINES    =========================!
!================================================================================!
!
!
!

!
!
!================================================================================!
!=======================       SUPPORT ROUTINES       ===========================!
!================================================================================!
!
!
!
 FUNCTION ROBIN_FUNCTION_CGSEM_SL2D( u, x, y, nHat, params ) RESULT( rFunc )
 ! FUNCTION ROBIN_FUNCTION_CGSEM_SL2D
 !
 !   For a problem with Robin type boundary conditions, the flux is given
 !   as a function of the solution. That function is specified here.
 !
 !
 ! ============================================================================= !
 ! DECLARATIONS
   IMPLICIT NONE
   real(prec)         :: u, x, y, nHat(1:2)
   TYPE( RUN_PARAMS ) :: params
   real(prec) :: rFunc
 
      rFunc = ZERO
 

 END FUNCTION ROBIN_FUNCTION_CGSEM_SL2D
!
!
!
  FUNCTION FLUX_X( u, dudx, dudy, a, b, c ) RESULT( fx )
 ! FUNCTION FLUX_X
 !
 !   In the problem Div( F ) = s, the x-component of the flux is calculated here.
 !
 !   The x-component of F is given by
 !   
 !       Fx = (a*dudx + b*dudy)/c
 ! ============================================================================= !
 ! DECLARATIONS
   IMPLICIT NONE
   real(prec) :: u, dudx, dudy, fx, a ,b ,c


      fx = (a*dudx + b*dudy)/c


 END FUNCTION FLUX_X
!
!
!
  FUNCTION FLUX_Y( u, dudx, dudy, a, b, c ) RESULT( fy )
 ! FUNCTION FLUX_Y
 !
 !   In the problem Div( F ) = s, the y-component of the flux is calculated here.
 !
 !   The y-component of F is given by
 !   
 !       Fy = (a*dudx + b*dudy)/c
 ! ============================================================================= !
 ! DECLARATIONS
   IMPLICIT NONE
   real(prec)         :: u, dudx, dudy, fy, a ,b ,c

    

      fy = (a*dudx + b*dudy)/c


 END FUNCTION FLUX_Y
!
!
!
SUBROUTINE FLUX_DIVERGENCE( geometry, cgStorage, u, Lu, nS, nP, params )
 ! S/R FLUX_DIVERGENCE
 !
 !
 ! ============================================================================= !
 ! DECLARATIONS
   IMPLICIT NONE
   integer, intent(in)                  :: nS, nP
   TYPE( MAPPEDGEOM_2D ), intent(in)    :: geometry
   TYPE( NODAL_STORAGE_2D ), intent(in) :: cgStorage
   real(prec), intent(in)               :: u(0:nS, 0:nP)
   real(prec), intent(out)              :: Lu(0:nS, 0:nP)
   TYPE( RUN_PARAMS ), intent(in)       :: params
   ! LOCAL
   integer    :: iS, iP
   real(prec) :: F1(0:nS,0:nP)
   real(prec) :: F2(0:nS,0:nP)
   real(prec) :: dF1ds(0:nS, 0:nP)
   real(prec) :: dF2dp(0:nS, 0:nP)
   real(prec) :: dudx(0:nS, 0:nP)
   real(prec) :: dudy(0:nS, 0:nP)
   real(prec) :: J, dxds, dxdp, dyds, dydp, Fx, Fy

      

      CALL CALC_GRADIENT( geometry, cgStorage, u, dudx, dudy, nS, nP ) ! calculate the gradient in physical coordinates

      ! Now apply the metric terms
      do iP = 0, nP ! Loop over the second computational direction
         do iS = 0, nS ! Loop over the first computational direction

            J = geometry % J(iS,iP)

            dxds = geometry % dxds(iS,iP)
            dxdp = geometry % dxdp(iS,iP)
            dyds = geometry % dyds(iS,iP)
            dydp = geometry % dydp(iS,iP)
            

            ! Calculate the conservative flux
            Fx = FLUX_X( u(iS,iP), dudx(iS,iP), dudy(iS,iP), ONE, ZERO, ONE )
            Fy = FLUX_Y( u(iS,iP), dudx(iS,iP), dudy(iS,iP), ZERO, ONE, ZERO )

            ! Contravariant flux calculation
            F1(iS,iP) = ( Fx*dydp - Fy*dxdp )*cgStorage % qWeightX(iS)

            F2(iS,iP) = ( Fy*dxds - Fx*dyds )*cgStorage % qWeightY(iP)
        

         enddo ! iS
      enddo ! iP


   ! Now calculate the divergence of the flux

   do iP = 0, nP ! loop over the second computational direction

      dF1ds(0:nS,iP) = MATVECMUL_TRANSPOSE( cgStorage % dMatX(0:nS,0:nS), F1(0:nS,iP), nS )  

   enddo ! iP

   do iS = 0, nS ! loop over the second computational direction

      dF2dp(iS,0:nP) = MATVECMUL_TRANSPOSE( cgStorage % dMatY(0:nP,0:nP), F2(iS,0:nP), nP )  

   enddo ! iP

   ! Now calculate the Laplacian

   do iP = 0, nP
      do iS = 0, nS

         Lu(iS,iP) = -( dF1ds(iS,iP)*cgStorage % qWeightY(iP) +&
                        dF2dp(iS,iP)*cgStorage % qWeightX(iS) )

      enddo ! iS
   enddo ! iP
       

 END SUBROUTINE FLUX_DIVERGENCE 
!
!
!
 SUBROUTINE MATRIX_ACTION_CGSEMSL2D( myCGSEM, u, Au )
 !
 ! 
 ! ==========================================================!
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS( CGSEM_SL2D ), intent(inout) :: myCGSEM
   real(prec), intent(inout)       :: u(0:myCGSEM % cgStorage % nS, &
                                        0:myCGSEM % cgStorage % nP, &
                                        1:myCGSEM % mesh % nElems )
   real(prec), intent(out)         :: Au(0:myCGSEM % cgStorage % nS, &
                                         0:myCGSEM % cgStorage % nP, &
                                         1:myCGSEM % mesh % nElems  )
   ! Local
   integer :: iEl, nS, nP, iS, iP, iEdge, nEdges
   integer    :: e1, s1, nB
   real(prec) :: temp(0:myCGSEM % cgStorage % nS, &
                       0:myCGSEM % cgStorage % nP )
   real(prec) :: Atemp(0:myCGSEM % cgStorage % nS, &
                       0:myCGSEM % cgStorage % nP )
   real(prec) :: nHat(1:2), length, rhs, w1, w2, x, y


     nS = myCGSEM % cgStorage % nS
     nP = myCGSEM % cgStorage % nP
     nEdges = myCGSEM % mesh % nEdges
     
     CALL myCGSEM % UNMASK( u ) ! unmask the solution

!$OMP PARALLEL PRIVATE( temp, Atemp )
!$OMP DO

     do iEl = 1, myCGSEM  % mesh % nElems

        temp(0:nS,0:nP) = u(0:nS,0:nP,iEl)


        CALL FLUX_DIVERGENCE( myCGSEM % mesh % elements(iEl) % geometry, &
                              myCGSEM % cgStorage, &
                              temp, &
                              Atemp, nS, nP, myCGSEM % params )

        Au(0:nS,0:nP,iEl) = Atemp(0:nS,0:nP)

     enddo
!$OMP END DO
!$OMP FLUSH(Au)
!$OMP END PARALLEL


     ! Add the contributions from the shared edges and corners
     CALL myCGSEM % GLOBAL_SUM( Au )

     ! Mask the edges in preparation for residual calculation

     CALL myCGSEM % MASK( Au )

     ! Re-mask the solution

     CALL myCGSEM % MASK( u )


 END SUBROUTINE MATRIX_ACTION_CGSEMSL2D
!
!
!
 SUBROUTINE COUNT_DEGREES_OF_FREEDOM( myCGSEM, dummy )
 ! SUBROUTINE COUNT_DEGREES_OF_FREEDOM
 ! 
 ! This subroutine counts the degrees of freedom in the elliptic solver
 ! by applying a mask to dummy solution array of all one's. The sum of all
 ! the entries in this array is equal to the number of degrees of freedom
 ! since the masked out points are masked to zero.
 !
 ! Additionaly a mapping array of dimension (1:nFree,1:3)
 ! is stored which holds the element ID and local quadrature node IDs.
 !
 !
 ! This information is used to construct and store the matrix implied by the CGSEM
 ! algorithm, and to map from the "1-D" array to the solution storage array.
 !
 ! 
 ! ==================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS( CGSEM_SL2D ), intent(inout) :: myCGSEM
   real(prec), intent(out)            :: dummy(0:myCGSEM % cgStorage % nS, &
                                              0:myCGSEM % cgStorage % nP, &
                                              1:myCGSEM % mesh % nElems )
   ! LOCAL
   integer :: iS, iP, iEl, nS, nP, nEl, nFree

       nS = myCGSEM % cgStorage % nS
       nP = myCGSEM % cgStorage % nP
       nEl = myCGSEM % mesh % nElems 
 

       ! Set all of the dummy values to ONE

       dummy = ONE

       ! Mask out the dummy values to obtain the unique solution locations

       CALL myCGSEM % MASK( dummy )

       ! Loop over all of the points and count the number of degrees of freedom
       nFree = 0
       
       do iEl = 1, nEl
          do iP = 0, nP
             do iS = 0, nS

                if( dummy(iS,iP,iEl) == ONE )then
                   ! Note that "dummy" is converted to an integer here
                   nFree = nFree + dummy(iS,iP,iEl)
                endif


             enddo
          enddo
       enddo
   
       myCGSEM % nFree = nFree

 END SUBROUTINE COUNT_DEGREES_OF_FREEDOM
!
!
!
 SUBROUTINE GET_ARRAY_MAP( myCGSEM, dummy )
 ! SUBROUTINE GET_ARRAY_MAP
 ! 
 ! 
 !
 ! 
 ! ==================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS( CGSEM_SL2D ), intent(inout) :: myCGSEM
   real(prec), intent(in)             :: dummy(0:myCGSEM % cgStorage % nS, &
                                              0:myCGSEM % cgStorage % nP, &
                                              1:myCGSEM % mesh % nElems )
   ! LOCAL
   integer :: iS, iP, iEl, nS, nP, nEl, nFree

       nS = myCGSEM % cgStorage % nS
       nP = myCGSEM % cgStorage % nP
       nEl = myCGSEM % mesh % nElems 

       ! Mask out the dummy values to obtain the unique solution locations

       nFree = 0
       do iEl = 1, nEl
          do iP = 0, nP
             do iS = 0, nS

                if( dummy(iS,iP,iEl) == ONE )then

                   ! Note that "dummy" is converted to an integer here
                   nFree = nFree + dummy(iS,iP,iEl)

                   myCGSEM % arrayMap(nFree,1) = iEl ! Store the element ID
                   myCGSEM % arrayMap(nFree,2) = iS  ! Store the first computational direction node ID
                   myCGSEM % arrayMap(nFree,3) = iP  ! Store the second computational direction node ID

                endif

             enddo
          enddo
       enddo

 END SUBROUTINE GET_ARRAY_MAP
!
!
!
SUBROUTINE BUILD_MATRIX_FROM_ACTION( myCGSEM, dummy, A )
 ! SUBROUTINE BUILD_MATRIX_FROM_ACTION
 ! 
 ! This subroutine uses calls to the MATRIX_ACTION routine to build and store
 ! the matrix implied the CGSEM algorithm. An impulse function (1 at a node and zero
 ! at all others) is used to compute the columns of A
 !
 ! 
 ! ==================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS( CGSEM_SL2D ), intent(inout) :: myCGSEM
   real(prec), intent(out)            :: A(1:myCGSEM % nFree, 1:myCGSEM % nFree)
   real(prec), intent(in)             :: dummy(0:myCGSEM % cgStorage % nS, &
                                              0:myCGSEM % cgStorage % nP, &
                                              1:myCGSEM % mesh % nElems )
   ! LOCAL
   real(prec) :: IRF(0:myCGSEM % cgStorage % nS, &
                     0:myCGSEM % cgStorage % nP, &
                     1:myCGSEM % mesh % nElems )
   real(prec) :: AIRF(0:myCGSEM % cgStorage % nS, &
                      0:myCGSEM % cgStorage % nP, &
                      1:myCGSEM % mesh % nElems )
   integer :: iS, iP, iEl, nS, nP, nEl, jS, jP, jEl
   integer :: col, row

       nS = myCGSEM % cgStorage % nS
       nP = myCGSEM % cgStorage % nP
       nEl = myCGSEM % mesh % nElems 
 

       col = 0

       ! Loop over the "columns" of A
       do jEl = 1, nEl
          do jP = 0, nP
             do jS = 0, nS 

                if( dummy(jS,jP,jEl) == ONE )then ! then this is a unique point and not duplicated from a shared edge or node.

                   col = col + 1 ! increment the column

                   ! Build the impulse function
                   IRF = ZERO
                   IRF(jS,jP,jEl) = ONE

                   CALL myCGSEM % MATRIX_ACTION( IRF, AIRF ) ! Do not include robin boundary conditions
                   
                   row = 0
                   ! Loop over the "rows" of A
                   do iEl = 1, nEl
                      do iP = 0, nP
                         do iS = 0, nS
            
                            if( dummy(iS,iP,iEl) == ONE ) then
                            
                               row = row + 1 ! increment the row

                               A(row,col) = AIRF(iS,iP,iEl)

                            endif
 
                         enddo
                      enddo
                   enddo

                endif 
 
             enddo
          enddo
       enddo                              
                   
                                   



 END SUBROUTINE BUILD_MATRIX_FROM_ACTION
!
!
!
SUBROUTINE BUILD_WEIGHT_MATRIX( myCGSEM, dummy, W )
 ! SUBROUTINE BUILD_WEIGHT_MATRIX
 ! 
 ! 
 !
 ! 
 ! ==================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS( CGSEM_SL2D ), intent(inout) :: myCGSEM
   real(prec), intent(out)            :: W(1:myCGSEM % nFree, 1:myCGSEM % nFree)
   real(prec), intent(in)             :: dummy(0:myCGSEM % cgStorage % nS, &
                                              0:myCGSEM % cgStorage % nP, &
                                              1:myCGSEM % mesh % nElems )
   ! LOCAL
   real(prec) :: WsolForm(0:myCGSEM % cgStorage % nS, &
                          0:myCGSEM % cgStorage % nP, &
                          1:myCGSEM % mesh % nElems )
   real(prec) :: s
   integer :: iS, iP, iEl, nS, nP, nEl, jS, jP, jEl, nEdges, iEdge
   integer :: e1, s1, nB
   integer :: col, row
   real(prec) :: J, w1, w2, x, y, length, nHat(1:2), rhs

       nS = myCGSEM % cgStorage % nS
       nP = myCGSEM % cgStorage % nP
       nEl = myCGSEM % mesh % nElems 
       nEdges = myCGSEM % mesh % nEdges
 

       col = 0
       W = ZERO


      ! Loop over the "columns" of A
       do jEl = 1, nEl
          do jP = 0, nP
             do jS = 0, nS 

                w1 = myCGSEM % cgStorage % qWeightX(jS)
                w2 = myCGSEM % cgStorage % qWeightY(jP)
                J = myCGSEM % mesh % elements(jEl) % geometry % J(jS,jP)

                s = myCGSEM % source(jS,jP,jEl)

                WsolForm(jS,jP,jEl) = s*J*w1*w2

 
             enddo
          enddo
       enddo

     ! Add in ROBIN boundary conditions if an eigenvalue multiplies the boundary condition

     ! Loop over the edges and search for Robin boundary conditions
     do iEdge = 1, nEdges

        if( myCGSEM % mesh % edges(iEdge) % elementIDs(2) == ROBIN )then 

           ! Gather the local addressing information
           e1 = myCGSEM % mesh % edges(iEdge) % elementIDs(1)
           s1 = myCGSEM % mesh % edges(iEdge) % elementSides(1)

           nB = myCGSEM % mesh % sideMap(s1) ! get the local quadrature node ID for this boundary

           if( s1 == 1 .OR. s1 == 3 )then ! south or north boundary (respectively)

               do iS = 0, nS

                  ! Get nHat
                  CALL myCGSEM % mesh % elements(e1) % geometry % GET_NHAT_ATPOINT( s1, iS, nhat, length ) 
 
                  w1 = myCGSEM % cgStorage % qWeightX(iS)

                  x = myCGSEM % mesh % elements(e1) % geometry % x(iS, nB)
                  y = myCGSEM % mesh % elements(e1) % geometry % y(iS, nB)

                  rhs = ROBIN_FUNCTION_CGSEM_SL2D( myCGSEM % solution(iS, nB, e1), &
                                                   x, y, nhat, myCGSEM % params )

                   WsolForm(iS, nB, e1) =  WsolForm(iS, nB, e1) + rhs*w1*length

               enddo


           elseif( s1 == 2 .OR. s1 == 4 )then ! east or west boundary (respectively)

              do iP = 0, nP

                  ! Get nHat
                  CALL myCGSEM % mesh % elements(e1) % geometry % GET_NHAT_ATPOINT( s1, iP, nhat, length ) 
 
                  w2 = myCGSEM % cgStorage % qWeightY(iP)

                  x = myCGSEM % mesh % elements(e1) % geometry % x(nB, iP)
                  y = myCGSEM % mesh % elements(e1) % geometry % y(nB, iP)

                  rhs = ROBIN_FUNCTION_CGSEM_SL2D( myCGSEM % solution(nB, iP, e1), &
                                                   x, y, nhat, myCGSEM % params )

                   WsolForm(nB, iP, e1) =  WsolForm(nB, iP, e1) + rhs*w2*length

               enddo

           endif
 
        endif

     enddo ! iEdge, loop over the edges


       CALL myCGSEM % GLOBAL_SUM( WsolForm )

       CALL myCGSEM % MASK( WsolForm )

       ! Loop over the "columns" of A
       do jEl = 1, nEl
          do jP = 0, nP
             do jS = 0, nS 

                if( dummy(jS,jP,jEl) == ONE )then ! then this is a unique point and not duplicated from a shared edge or node.

                   col = col + 1 ! increment the column

                   W(col,col) = WsolForm(jS,jP,jEl)

                endif 
 
             enddo
          enddo
       enddo                              
                   
                                   



 END SUBROUTINE BUILD_WEIGHT_MATRIX
!
!
!
SUBROUTINE CALC_EIGENMODES( myCGSEM ) 
 ! S/R CALC_EIGENMODES
 !
 !  Calculates the Eigenmodes and eigenvalues for the Sturm-Liouville system 
 !
 !
 ! ==============================================================================!
 ! DECLARATION
   IMPLICIT NONE
    CLASS( CGSEM_SL2D ), intent(inout) :: myCGSEM
    ! LOCAL
    real(prec) :: A(1:myCGSEM % nFree,1:myCGSEM % nFree), W(1:myCGSEM % nFree,1:myCGSEM % nFree)
    real(prec) :: ReEv(1:myCGSEM % nFree)
    real(prec) :: ImEv(1:myCGSEM % nFree)
    real(prec) :: denEV(1:myCGSEM % nFree)
    real(prec) :: lVec(1:myCGSEM % nFree, 1:myCGSEM % nFree)
    real(prec) :: rVec(1:myCGSEM % nFree, 1:myCGSEM % nFree)
    real(prec) :: eVals(1:myCGSEM % nFree, 1:3)
    real(prec) :: work(1:8*(myCGSEM % nFree+1))
    real(prec) :: norm
    integer    :: info, iNode, iMode, iFree,nFree
    real(prec) :: dummy(0:myCGSEM % cgStorage % nS, &
                        0:myCGSEM % cgStorage % nP, &
                        1:myCGSEM % mesh % nElems )
  

       ! Set all of the dummy values to ONE

       dummy = ONE

       ! Mask out the dummy values to obtain the unique solution locations

       CALL myCGSEM % MASK( dummy )

      nFree = myCGSEM % nFree

      ! Build and store the mass+stiffness and weight matrices
      CALL myCGSEM % BUILD_MATRIX_FROM_ACTION( dummy, A ) ! Build the  "mass + stiffness" matrix
      CALL myCGSEM %  BUILD_WEIGHT_MATRIX( dummy, W ) ! Build the "weight" matrix

      if( prec==sp ) then !single precision

         CALL SGGEV('N','V', nFree, A(1:nFree,1:nFree), nFree, &      ! Sturm-Liouvill matrix
                                    W(1:nFree,1:nFree), nFree, &      ! Weight matrix
                                    ReEv(1:nFree), ImEv(1:nFree), & ! Real and complex component of eigenvalues (scaled) 
                                    denEv(1:nFree),&                ! Scale for eigenvalues
                                    lVec(1:nFree,1:nFree), nFree, &   ! left eigenvectors
                                    rvec(1:nFree,1:nFree), nFree, &
                                    work(1:8*(nFree+1)), 8*(nFree+1), &
                                    info )
                                     
        
      elseif( prec==dp )then ! double precision


         CALL DGGEV('N','V', nFree, A(1:nFree,1:nFree), nFree, &      ! Sturm-Liouvill matrix
                                    W(1:nFree,1:nFree), nFree, &      ! Weight matrix
                                    ReEv(1:nFree), ImEv(1:nFree), & ! Real and complex component of eigenvalues (scaled) 
                                    denEv(1:nFree),&                ! Scale for eigenvalues
                                    lVec(1:nFree,1:nFree), nFree, &   ! left eigenvectors
                                    rvec(1:nFree,1:nFree), nFree, &
                                    work(1:8*(nFree+1)), 8*(nFree+1), &
                                    info )
                                     
 
      else

         print*, 'MODULE CGSEM_SL2D_CLASS :: S/R CALC_EIGENMODES '
         print*, 'Real precision ', prec
         print*, 'Not supported under LAPACK. STOPPING!'
         
         STOP

      endif

      


      ! Sort and store the eigenvectors and eigen-modes
      if( info == 0 )then ! successful

         do iFree = 1, nFree


            eVals(iFree,1) = ReEv(iFree)/denEv(iFree)
            eVals(iFree,2) = ImEv(iFree)/denEv(iFree)
            eVals(iFree,3) = denEv(iFree)


         enddo

   
         ! And sort the eigenvalues and eigenfunctions in order of decreasing eigenvalue
         CALL SORT_EVALS( rVec, eVals, nFree )

         ! Map the eigenvectors to the solution mesh.

         CALL myCGSEM % ARRAY_TO_SOLUTION( rVec )

         ! Fill in the eigenvalues
         myCGSEM % eVals = eVals

      elseif( info > 0 )then

         print*, 'MODULE CGSEM_SL2D_CLASS :: S/R CALC_EIGENMODES '
         print*, 'S/R CALC_EIGENMODES '
         print*, 'LAPACK (S/D)GGEV QZ iteration failed. No eigenvectors have been calculated. STOPPING!'
         STOP
         
            
      endif

 END SUBROUTINE CALC_EIGENMODES
!
!
!
SUBROUTINE SORT_EVALS( eVecs, eVals, nFree )
 ! 
 ! Sorts the eigenvalues from largest to smallest, and swaps the eigenvectors as well.
 !
 ! A bubble sort algorithm is used
 !
 ! ================================================================================= !
 ! DECLARATIONS
   integer, intent(in) :: nFree 
   real(prec), intent(inout) :: eVecs(1:nFree,1:nFree), eVals(1:nFree,1:3)
   !LOCAL
   integer :: iS, jS, kS
   integer :: N, K
   real(prec) :: thisVal(1:3), eVec(1:nFree)


      do iS = nFree, 1, -1 ! 
 
         do jS = 2, iS ! 

            if( abs(eVals(jS-1,1)) > abs(eVals(jS,1)) )then ! swap the eigenvalues and eigenvectors !

               thisVal = eVals(jS-1,1:3) ! Grab temporary storage
               
               eVals(jS-1,1:3) = eVals(jS,1:3)

               eVals(jS,1:3) = thisVal

               ! Swap the eigenvectors

               eVec(1:nFree) = eVecs(1:nFree,jS-1)
                 
               eVecs(1:nFree,jS-1) = eVecs(1:nFree,jS)

               eVecs(1:nFree,jS) = eVec(1:nFree)

            endif


         enddo

      enddo ! iS, Loop over the sort-list entries, find the maximum
      
   


 END SUBROUTINE SORT_EVALS 
!
!
!
 SUBROUTINE ARRAY_TO_SOLUTION( myCGSEM, array )
 ! SUBROUTINE ARRAY_TO_SOLUTION
 ! 
 ! Takes the "array(1:nFree)" and uses the array map to reshape into the 
 ! solution form, "solForm(0:nS,0:nP,1:nEl)"
 ! 
 ! ==================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS( CGSEM_SL2D ), intent(inout) :: myCGSEM
   real(prec), intent(in)             :: array(1: myCGSEM % nFree, 1:myCGSEM % nFree)
   ! LOCAL
   integer :: iMode, iNode, iS, iP, iEl, nFree

      nFree = myCGSEM % nFree
      
      do iMode = 1, nFree ! Loop over the modes
        
         ! Mask the mode out to enforce Dirichlet boundary conditions  
         CALL myCGSEM % MASK( myCGSEM % eModes(:,:,:,iMode) )

         do iNode = 1, nFree ! Loop over the nodes
      
      
            iEl = myCGSEM % arrayMap(iNode,1) ! Obtain the element ID
            iS = myCGSEM % arrayMap(iNode,2)  ! Obtain the first computational direction quadrature node ID
            iP = myCGSEM % arrayMap(iNode,3)  ! Obtain the second computational direction quadrature node ID 
           
            ! Copy the solution over to the solution form
            myCGSEM % eModes(iS,iP,iEl,iMode) = array(iNode, iMode)

         enddo
 
          ! Unmask the solution
          CALL myCGSEM % UNMASK( myCGSEM % eModes(:,:,:,iMode) )

      enddo
 
     


 END SUBROUTINE ARRAY_TO_SOLUTION
!
!
!================================================================================!
!=======================       FILE I/O      ===========================!
!================================================================================!
!
!
!
 SUBROUTINE WRITE_TECPLOT_CGSEM_SL2D( myCGSEM, nPlot, nOld, plotInterp, Tmat, nModes )
 ! WRITE_TECPLOT_CGSEM_SL2D
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
  CLASS( CGSEM_SL2D ), intent(in)       :: myCGSEM
  integer, intent(in)                   :: nOld, nPlot, nModes
  real(prec), intent(in)                :: Tmat(0:nPlot, 0:nOld)
  TYPE( LAG_INTERP2D ), intent(in)      :: plotInterp
  !LOCAL
  real(prec)  :: x(0:nPlot,0:nPlot), y(0:nPlot,0:nPlot)!, depth(0:nPlot,0:nPlot)
  real(prec)  :: u(0:nPlot,0:nPlot)
  real(prec)  :: J(0:nPlot,0:nPlot)
  integer :: iX, iY, iZ, iEl, iMode
  character(len=5) :: zoneID, modeChar


    do iMode = 1, nModes

       write(modeChar,'(I5.5)') iMode

       open( unit=2, file= 'EigenMode'//modeChar//'.tec', form='formatted',status='replace')

       write(2,*) 'VARIABLES = "X", "Y", "EigenMode"'
    
       do iEL = 1, myCGSEM % mesh % nElems

          ! Interpolate the solutions onto a uniform mesh
          CALL COARSE_TO_FINE_2D(myCGSEM % cgStorage % interp,&
                                      myCGSEM % mesh % elements(iEl) % geometry % x,&
                                      plotInterp, x, Tmat, Tmat)


          CALL COARSE_TO_FINE_2D(myCGSEM % cgStorage % interp,&
                                      myCGSEM % mesh % elements(iEl) %  geometry % y,&
                                      plotInterp, y, Tmat, Tmat)

          CALL COARSE_TO_FINE_2D(myCGSEM % cgStorage % interp,&
                                      myCGSEM % eModes(:,:,iEl,iMode),&
                                      plotInterp, u, Tmat, Tmat)
 
    

           write(zoneID,'(I5.5)') iEl


           write(2,*) 'ZONE T="el'//trim(zoneID)//'", I=',nPlot+1,', J=', nPlot+1,',F=POINT'

           do iY = 0, nPlot
              do iX = 0, nPlot

                 write (2,*)  x( iX, iY ), y( iX, iY ), u(iX,iY)

              enddo
           enddo

       enddo

       close(unit=2)

    enddo


 END SUBROUTINE WRITE_TECPLOT_CGSEM_SL2D
!
!
!
 SUBROUTINE WRITE_EVALS_TO_CURVE( myCGSEM )
 ! S/R WRITE_EVALS_TO_CURVE
 !
 !  This subroutine writes the eigenvalues to three separate curve files:
 !   real, complex, and denominator.
 !
 ! ======================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS( CGSEM_SL2D ), intent(in) :: myCGSEM
   ! LOCAL
   integer :: iMode, nFree, nUnit

      nFree = myCGSEM % nFree

      open( unit=NEWUNIT(nUnit), file='evals-real.curve' )

      write( nUnit, * ) '#realpart-eigenvalue' 

      do iMode = 1, nFree

         write( nUnit, '(I5, E19.7)' ) iMode, myCGSEM % eVals(iMode,1)

      enddo

      close(nUnit)


      open( unit=NEWUNIT(nUnit), file='evals-complex.curve' )

      write( nUnit, * ) '#imagpart-eigenvalue' 

      do iMode = 1, nFree

         write( nUnit, '(I5, E19.7)' ) iMode, myCGSEM % eVals(iMode,2)

      enddo

      close(nUnit)


      open( unit=NEWUNIT(nUnit), file='evals-den.curve' )

      write( nUnit, * ) '#denominator-eigenvalue' 

      do iMode = 1, nFree

         write( nUnit, '(I5, E19.7)' ) iMode, myCGSEM % eVals(iMode,3)

      enddo

      close(nUnit)


 END SUBROUTINE WRITE_EVALS_TO_CURVE

END MODULE CGSEM_SL2D_CLASS
