MODULE DGSEM_BAROTROPIC_CLASS
! DGSEM_BAROTROPIC_CLASS.f90
! 
! Module History 
! 
! o  February 2014
!    
!
! 
! List of Module Procedures
!
! ------ CONSTRUCTION ROUTINES -----
!
!
!
! ------ EVALUATION ROUTINES -------
!
!
! ------ DERIVATIVE ROUTINES --------
!
!
! ----- MESH REFINEMENT ROUTINES ------
!
!
!
!
!  
! =======================================================================================

! "The basics"
USE COMMONDATA 
USE LEGENDRE
USE MATRIXROUTINES
USE RUN_PARAMS_CLASS

!
!
! Interpolation modules
USE LAGRANGE_1D_CLASS
USE LAGRANGE_2D_CLASS
!
!
! Storage Container modules
USE GHOST_ELEMENT_CLASS
USE NODAL_STORAGE_2D_CLASS
USE DG2D_SOLUTION_STORAGE_CLASS
USE DG2D_RELAXATION_FIELD_CLASS

! Geometry Modules
USE MAPPEDGEOM_2D_CLASS
USE BARO_MAPPEDGEOM_CLASS

! Mesh Primitives
USE GEOMETRY_BASICS
USE EDGE_CLASS
USE BARO_QUADELEMENT_CLASS  !****
USE BARO_QUADMESH_CLASS     !****
 
! Additional modules for shallow water equations
USE BAROTROPIC_PROPERTIES_CLASS !****


! Parallel Libraries
USE OMP_LIB


IMPLICIT NONE




    ! Two-dimensional polynomial is the tensor
    ! product of two one-dimensional polynomials

    TYPE DGSEM_BAROTROPIC
      integer                                    :: nEq, nGhosts
      type( DG2D_SOLUTION ), allocatable         :: sol(:)
      type( BAROTROPIC_PROPERTIES), allocatable  :: diags(:) 
      type( NODAL_STORAGE_2D )                   :: dGStorage
      type( BARO_QUADMESH )                      :: mesh
      type( GHOST_ELEMENT), allocatable          :: ghostEl(:)
      type( RUN_PARAMS )                         :: params
      type( DG2D_RELAXATION_FIELD ), allocatable :: relax(:)

      CONTAINS

      ! Manual constructors/Destructors
      PROCEDURE :: BUILD => BUILD_DGSEM_BAROTROPIC
      PROCEDURE :: TRASH => TRASH_DGSEM_BAROTROPIC
      ! Accessors

      
       ! Type Specific Routines
       PROCEDURE :: GLOBAL_TIME_DER_BAROTROPIC
       PROCEDURE :: MAPPED_DGSEM_BAROTROPIC_RK3
       PROCEDURE :: UPDATE_GHOSTS
       PROCEDURE :: EDGE_FLUXES

       PROCEDURE ::  DIAGNOSE_INTEGRATED_ENERGIES

       PROCEDURE :: WRITE_TECPLOT => WRITE_BAROTROPIC_TECPLOT
       PROCEDURE :: WRITE_TECPLOT_ASIS => WRITE_BAROTROPIC_TECPLOT_ASIS
       PROCEDURE :: WRITE_PICKUP => WRITE_BAROTROPIC_PICKUP 

    END TYPE DGSEM_BAROTROPIC



 CONTAINS
!
!
!=========================================================================!
!------------------- NODAL DG CONSTRUCTION ROUTINES ----------------------!
!=========================================================================!
!
!
 SUBROUTINE BUILD_DGSEM_BAROTROPIC( myDGSEM, nEq, initialize, theMPICOMM, myRank )
 ! S/R BUILD_DGSEM_BAROTROPIC
 !  Desription:
 !    
 !
 !    Subroutine dependencies :
 !    
 !
 !  Input :
 ! 
 !    
 !
 !  Output :
 !   
 !
 !------------------------------------------------------------------------
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS(DGSEM_BAROTROPIC), intent(inout) :: myDGSEM
   logical, intent(in)                    :: initialize
   integer, intent(in)         :: nEq, theMPICOMM, myRank
   !LOCAL
   integer :: iEl, nGhosts, globElID, pID2, recStart, recEnd, nX, nY
   integer :: rStartRel, rEndRel
   character(10) :: pickupIterChar
   character(5)  :: pIDChar
   character(40) :: meshFile


  

   myDGSEM % nEq = nEq

   CALL myDGSEM % params % BUILD( )
   nX = myDGSEM % params % polyDeg
   nY = nX

   write(pickupIterChar,'(I10.10)') myDGSEM % params % iterInit
   write(pIDChar,'(I5.5)') myRank
 
   CALL myDGSEM % dGStorage % BUILD( nX, nY, DG )

   meshFile = 'mesh'//'.'//pIDChar

   ! Builds the lateral mesh

 !  CALL myDGSEM % mesh % BUILD( trim(meshFile), myDGSEM % dGStorage )
   CALL myDGSEM % mesh % READ_PICKUP( trim(meshfile), myDGSEM % dGStorage  ) 
  
   ALLOCATE( myDGSEM % sol(1:myDGSEM % mesh % nElems) )

   ALLOCATE( myDGSEM % diags(1:myDGSEM % mesh % nElems) )

   ALLOCATE( myDGSEM % relax(1:myDGSEM % mesh % nElems) )


   ! If the iter0  (initial iterate) is negative this signifies that no initial condition is given
   ! Default to a "zero for all" initial condition
   if( INITIALIZE )then

      do iEl = 1, myDGSEM % mesh % nElems

         CALL myDGSEM % sol(iEl) % BUILD( nX, nY, nEq ) 

         CALL myDGSEM % relax(iEl) % BUILD( nX, nY, nEq ) 

         CALL myDGSEM % diags(iEl) % BUILD( nX, nY )
      
      enddo

   else

      recStart = 0
      rStartRel = 0
      print*, 'State.'//pIDChar//'.'//pickupIterChar
      
      do iEl = 1, myDGSEM % mesh % nElems

         CALL myDGSEM % sol(iEl) % BUILD( nX, nY, nEq ) 

         CALL myDGSEM % sol(iEl) % READ_PICKUP( nX, nY, nEq, recStart, recEnd, 'State.'//pIDChar//'.'//pickupIterChar ) 
         recStart = recEnd

         CALL myDGSEM % relax(iEl) % BUILD( nX, nY, nEq ) 

         CALL myDGSEM % relax(iEl) % READ_PICKUP( nX, nY, nEq, rStartRel, rEndRel, 'State.'//pIDChar//'.'//pickupIterChar ) 
         rStartRel = rEndRel

         CALL myDGSEM % diags(iEl) % BUILD( nX, nY )
      

      enddo

   endif
  

   ! Count the number of boundary edges and allocate space for the ghost elements

   nGhosts = 0

   do iEl = 1, myDGSEM % mesh % nEdges

      ! Get proc ID for the process which holds the secondary element
      globElID = myDGSEM % mesh % edges(iEl) % elementIDs(2)

      if( globElID <= NO_NORMAL_FLOW )then
     
         pID2 = myRank

      else

         pID2 = myDGSEM % mesh % glob2LocElemMap( globElID, 2 )

      endif

      

      if( globElID <= NO_NORMAL_FLOW .OR. pID2 /= myRank )then

         nGhosts = nGhosts + 1

      endif

   enddo

 
   myDGSEM % nGhosts = nGhosts
   
   ALLOCATE( myDGSEM % ghostEl(1:nGhosts) )
   
   ! Build ghost element
   do iEl = 1, nGhosts
 
      CALL myDGSEM % ghostEl(iEl) % BUILD( max(nX,nY), nEq )
 
   enddo
 

   nGhosts = 0
   do iEl = 1, myDGSEM % mesh % nEdges

      ! Get proc ID for the process which holds the secondary element
      globElID = myDGSEM % mesh % edges(iEl) % elementIDs(2)

      if( globElID <= NO_NORMAL_FLOW )then
     
         pID2 = myRank

      else

         pID2 = myDGSEM % mesh % glob2LocElemMap( globElID, 2 )

      endif

      if( globElID <= NO_NORMAL_FLOW .OR. pID2 /= myRank )then

         nGhosts = nGhosts + 1

        ! Associate a ghost element address with an internal element ID
         myDGSEM % ghostEl(nGhosts) % intElemID = myDGSEM % mesh % edges(iEl) % elementIDs(1)


         myDGSEM % ghostEl(nGhosts) % intEdgeID = myDGSEM % mesh % edges(iEl) % elementSides(1)

       ! and an external element ID
         myDGSEM % ghostEl(nGhosts) % extElemID = globElID ! Ghost element now uses the element ID it associates with
        
         myDGSEM % ghostEl(nGhosts) % extEdgeID = myDGSEM % mesh % edges(iEl) % elementSides(2)

       ! Now change the external element ID to an index for the ghost array

         myDGSEM % mesh % edges(iEl) % elementIDs(2) = -nGhosts ! Use a negative value as a flag to to use a ghost value

       
      endif

   enddo


   
 END SUBROUTINE BUILD_DGSEM_BAROTROPIC
!
!
!
 SUBROUTINE TRASH_DGSEM_BAROTROPIC( myDGSEM )
 ! S/R TRASH_DGSEM_BAROTROPIC
 !  Desription:
 !  
 !
 !    Subroutine dependencies :
 !   
 !
 !  Input :
 ! 
 !  
 !
 !  Output :
 !    
 !
 !------------------------------------------------------------------------
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS(DGSEM_BAROTROPIC), intent(inout) :: myDGSEM
   ! LOCAL
   integer :: iEl
 

   do iEl = 1, myDGSEM % mesh % nElems

      CALL myDGSEM % sol(iEl) % TRASH( )

      CALL myDGSEM % diags(iEl) % TRASH( )

      CALL myDGSEM % relax(iEl) % TRASH( )
  
   enddo

   CALL myDGSEM % dGStorage % TRASH( )

   CALL myDGSEM % mesh % TRASH( )

   DEALLOCATE( myDGSEM % sol )

   DEALLOCATE( myDGSEM % diags ) 

   DEALLOCATE( myDGSEM % relax )


   do iEl = 1, myDGSEM % nGhosts

      CALL myDGSEM % ghostEl(iEl) % TRASH( )

   enddo

   DEALLOCATE( myDGSEM % ghostEl )


 END SUBROUTINE TRASH_DGSEM_BAROTROPIC
!
!
!=========================================================================!
!-------------------------- ACCESSOR ROUTINES ----------------------------!
!=========================================================================!
!
!

SUBROUTINE MAPPED_DGSEM_BAROTROPIC_RK3( myDGSEM, tn, theMPICOMM, myRank )
 ! S/R MAPPED_DGSEM_BAROTROPIC_RK3( 3rd order Runge-Kutta)
 ! 
 ! Description :
 !  
 !
 !    Subroutine Dependencies :
 !    
 !
 ! Input :
 !     real(prec) :: f(0:myInterp%nNodes)
 ! 
 ! Output :
 !     real(prec) :: f(0:myInterp%nNodes)
 !
 !------------------------------------------------------------------------
 ! DECLARATIONS
   CLASS(DGSEM_BAROTROPIC), intent(inout) :: myDGSEM
   real(prec), intent(in)                 :: tn
   integer, intent(in)                    :: theMPICOMM, myRank
   ! LOCAL
   real(prec) :: t, dSdt, dt
   real(prec) :: G2D(0:myDGSEM % dGStorage % nS,&
                     0:myDGSEM % dGStorage % nP,&
                     1:myDGSEM % nEq,1:myDGSEM % mesh % nElems) 
   integer    :: m, iP, jP, nX, nY, nZ, nEq, iEqn, iEl


     nEq = myDGSEM % nEq 

     dt = myDGSEM % params % dt

     CALL myDGSEM% dgStorage % GET_N_NODES( nX, nY )



      G2D = 0.0_prec

    
      do m = 1,3 ! Loop over RK3 steps

         t = tn + rk3_b(m)*dt

         ! Calculate the tendency
         CALL myDGSEM % GLOBAL_TIME_DER_BAROTROPIC( t, theMPICOMM, myRank )
        

         do iEl = 1, myDGSEM % mesh % nElems ! Loop over all of the elements
 

            do iP = 0,nX ! Loop over the quadrature points (x)
 
               do jP = 0,nY ! Loop over the quadrature points (y)

                  do iEqn = 1,nEq ! Loop over the equations

                     CALL myDGSEM % sol(iEl) % GET_TENDENCY_ATNODE( iEqn, iP, jP, dSdt )

                     G2D(iP,jP,iEqn, iEl) = rk3_a(m)*G2D(iP,jP,iEqn, iEl) + dSdt

                     myDGSEM % sol(iEl) % solInterior(iP,jP,iEqn) = &
                               myDGSEM % sol(iEl) % solInterior(iP,jP,iEqn) + &
                               rk3_g(m)*dt*G2D(iP,jP,iEqn,iEl)

                   enddo ! iEqn, loop over the equations

                enddo ! jP, loop over the quadrature points (y) 

             enddo ! iP, loop over the quadrature points (x)

         enddo ! iEl, loop over all of the elements
 
      enddo ! m, loop over the RK3 steps
   
   
       
 END SUBROUTINE MAPPED_DGSEM_BAROTROPIC_RK3
!
!
!
 SUBROUTINE GLOBAL_TIME_DER_BAROTROPIC( myDGSEM, tn, theMPICOMM, myRank ) 
 ! S/R GLOBAL_TIME_DER_BAROTROPIC
 ! 
 ! Calculates the tendency terms for the shallow water system
 !
 ! Description :
 !  
 !
 !    Subroutine Dependencies :
 !    
 !
 ! Input :
 !
 ! 
 ! Output :
 !
 !
 !------------------------------------------------------------------------
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS(DGSEM_BAROTROPIC), intent(inout) :: myDGSEM
   real(prec), intent(in)        :: tn
   integer, intent(in)           :: theMPICOMM, myRank
   ! Local
   integer :: iEl, iEdge, nEq, iGhost, iError
   real(prec) :: cpuStart, cpuEnd

      nEq = myDGSEM % nEq

! Open the parallel block
!$OMP PARALLEL


!$OMP DO
!
! For each element, obtain the solution at the element edges via interpolation
!
    !   print*, 'Rank', myRank,': CALC_SWSOL_AT_BOUNDS : Begin'
      do iEl = 1, myDGSEM % mesh % nElems

         ! Use the interior solution to update solution at cell faces
         CALL CALC_SWSOL_AT_BOUNDS( myDGSEM % sol(iEl),&
                                    myDGSEM % dGStorage ) 


      enddo 

    !   print*, 'Rank', myRank,': CALC_SWSOL_AT_BOUNDS : End'
!$OMP END DO 

      
      CALL MPI_BARRIER( theMPICOMM, iError )

     ! The only stage where communication between processes is required.
     ! Update the Ghost Points
    !   print*, 'Rank', myRank,': UPDATE_GHOSTS : Begin'

       CALL myDGSEM % UPDATE_GHOSTS( tn, myDGSEM % dGStorage % nS, &
                                       myDGSEM % mesh % nElems, &
                                       myDGSEM % nGhosts, &
                                       nEq, theMPICOMM, myRank )   
      
     !  print*, 'Rank', myRank,': UPDATE_GHOSTS : End'
!
!$OMP DO
!
! For each edge, compute the element fluxes
!
   !   print*, 'Rank', myRank,': EDGE_FLUXES : Begin'

      do iEdge = 1, myDGSEM % mesh % nEdges


         CALL myDGSEM % EDGE_FLUXES( tn, iEdge, myRank, myDGSEM % params )

    
      enddo ! iEdge, cycle over the edges

    !  print*, 'Rank', myRank,': EDGE_FLUXES : End'

     ! Note** : the flux calculations also update the boundary values to equal half of the neighboring values


     !  print*, 'Rank', myRank,': DIAGNOSTICS : Begin'
      do iEl = 1,myDGSEM % mesh % nElems

         
       ! Calculate the vorticity
         
         CALL myDGSEM % diags(iEl) % COMPUTE_RELATIVE_VORTICITY( myDGSEM % dGStorage, &
                                                                 myDGSEM % sol(iEl), &
                                                                 myDGSEM % mesh % elements(iEl) % geometry, & 
                                                                 myDGSEM % params % f0, &
                                                                 myDGSEM % params % betaX, &
                                                                 myDGSEM % params % betaY, &
                                                                 myDGSEM % params % MODEL_FORMULATION, & 
                                                                 myDGSEM % dGStorage % nS, &
                                                                 myDGSEM % dGStorage % nP) 

      enddo

      ! print*, 'Rank', myRank,': DIAGNOSTICS : End'
      
!$OMP END DO
!
!$OMP DO
!
! For each element, calculate the time derivative
!
     !  print*, 'Rank', myRank,': MAPPED_DGSW_TIMEDER : Begin'
      do iEl = 1, myDGSEM % mesh % nElems

         ! Calculate the tendency
         CALL MAPPED_DGSW_TIMEDER( tn, nEq, myDGSEM % dGStorage % nS, &
                                   myDGSEM % dGStorage % nP, &
                                   myDGSEM % dGStorage, & 
                                   myDGSEM % mesh % elements(iEl) % geometry, &  
                                   myDGSEM % sol(iEl), &
                                   myDGSEM % diags(iEl) , &
                                   myDGSEM % params, &
                                   myDGSEM % relax(iEl) )

      enddo

      ! print*, 'Rank', myRank,': MAPPED_DGSW_TIMEDER : End'
!$OMP END DO
!
! Close the parallel block
!
!$OMP END PARALLEL
 
 END SUBROUTINE GLOBAL_TIME_DER_BAROTROPIC
!
!
! 
 SUBROUTINE EDGE_FLUXES( myDGSEM, tn, iEdge, myRank, locParams )
 ! S/R EDGE_FLUXES
 ! 
 ! Description :
 !  
 !
 !    Subroutine Dependencies :
 !    
 !
 ! Input :
 !
 ! 
 ! Output :
 !
 !
 !------------------------------------------------------------------------
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS( DGSEM_BAROTROPIC ), intent(inout) :: myDGSEM
   integer, intent(in)                      :: iEdge
   integer, intent(in)                      :: myRank
   real(prec), intent(in)                   :: tn
   type( RUN_PARAMS ), intent(in)           :: locParams
   
   ! Local
   integer :: k, iNode, iEqn
   integer :: e1, s1, e2, s2, iError
   real(prec) :: flux(1:myDGSEM % nEq)
   real(prec) :: inState(1:myDGSEM % nEq), extState(1:myDGSEM % nEq)
   real(prec) :: wavespeed, h 
   real(preC) :: nHat(1:2), nHatLength 
    

      e1 = myDGSEM % mesh % edges(iEdge) % elementIDs(1) ! This is a global element ID
      ! Obtain the local element ID
              
      e1 = myDGSEM % mesh % glob2LocElemMap(e1,3)
 
      s1 = myDGSEM % mesh % edges(iEdge) % elementSides(1)
 
      e2 = myDGSEM % mesh % edges(iEdge) % elementIDs(2)  ! This is a global element ID, or a ghost array index

      s2 = abs( myDGSEM % mesh % edges(iEdge) % elementSides(2) )

      if( e2 > 0 )then ! this is an interior edge

         k = myDGSEM % mesh % edges(iEdge) % start - myDGSEM % mesh % edges(iEdge) % inc
         e2 = myDGSEM % mesh % glob2LocElemMap(e2,3)

         do iNode = 0, myDGSEM % dGStorage % nS ! Loop over the nodes

           ! Obtain the local element ID for the secondary element

            

            CALL myDGSEM % sol(e1) % GET_SOLBOUND_ATPOINT( s1, iNode, inState ) ! Get the interior state

            CALL myDGSEM % sol(e2) % GET_SOLBOUND_ATPOINT( s2, k, extState ) ! Get the exterior state

            h = 0.5_prec*(myDGSEM % mesh % elements(e1) % geometry % dBound(iNode, s1) +&
                myDGSEM % mesh % elements(e2) % geometry % dBound(k, s2)) ! get the depth on the primary element boundary (s1)
           ! waveSpeed = sqrt( g*h ) ! Get the wave speed at the boundary
            
            CALL myDGSEM % mesh % elements(e1) % geometry % GET_NHAT_ATPOINT( s1, iNode, nHat, nHatLength ) ! Get nHat

            ! Calculate the RIEMANN flux
            flux = RIEMANN_LAXFRIEDRICH( inState, extState, myDGSEM % nEq, h, nHat, locParams)*nHatLength

            ! Store the flux for the elements which share this edge
            CALL myDGSEM % sol(e1) % SET_FLUXBOUND_ATPOINT ( s1, iNode, flux )

            CALL myDGSEM % sol(e2) % SET_FLUXBOUND_ATPOINT ( s2, k, -flux )
 
           ! Reset the the solutions at this boundary point as the average of the two states
           ! This is meant to communicate the neighbor information to the vorticity calculation
           ! CALL myDGSEM % sol(e1) % SET_SOLBOUND_ATPOINT( s1, iNode, (inState+extState)*0.5_prec )
        
           ! CALL myDGSEM % sol(e2) % SET_SOLBOUND_ATPOINT( s2, k, (inState+extState)*0.5_prec )

            k = k + myDGSEM % mesh % edges(iEdge) % inc


         enddo ! iNode, loop over the nodes


      else ! this edge is a boundary edge, we need to obtain ghost values


          e2 = abs(e2) ! the array address for the ghost values
  
          k = myDGSEM % mesh % edges(iEdge) % start - myDGSEM % mesh % edges(iEdge) % inc

         do iNode = 0, myDGSEM % dGStorage % nS ! loop over the nodes on this edge
       

            CALL myDGSEM % sol(e1) % GET_SOLBOUND_ATPOINT( s1, iNode, inState ) ! Get the interior state

            CALL myDGSEM % mesh % elements(e1) % geometry % GET_NHAT_ATPOINT( s1, iNode, nHat, nHatLength ) ! Get nHat

            CALL myDGSEM % ghostEl(e2) % GET_EXTSTATE( k, myDGSEM % nEq, extState ) ! Get the external state for the boundary

             h = myDGSEM % mesh % elements(e1) % geometry % dBound(iNode, s1) ! get the depth on the primary element boundary (s1)
          !  waveSpeed = sqrt( g*h ) ! Get the wave speed at the boundary

            ! Calculate the RIEMANN flux
            flux = RIEMANN_LAXFRIEDRICH( inState, extState, myDGSEM % nEq, h, nHat, locParams)*nHatLength

            CALL myDGSEM % sol(e1) % SET_FLUXBOUND_ATPOINT ( s1, iNode, flux )

            ! Reset the the solutions at this boundary point as the average of the two states
           ! This is meant to communicate the neighbor information to the vorticity calculation
          !  CALL myDGSEM % sol(e1) % SET_SOLBOUND_ATPOINT( s1, iNode, (inState+extState)*0.5_prec )

   
            k = k + myDGSEM % mesh % edges(iEdge) % inc

         enddo ! iNode, loop over the nodes on this edge


      endif ! Choosing the edge type
 
 END SUBROUTINE EDGE_FLUXES
!
!
!
 SUBROUTINE MAPPED_DGSW_TIMEDER( tn, nEq, nX, nY, dGStorage, geometry, dGSol, dG2DProps, locParams, locRelax ) 
 ! S/R MAPPED_DGSW_TIMEDER
 ! 
 ! Description :
 !  
 !
 !    Subroutine Dependencies :
 !    
 !
 ! Input :
 !
 ! 
 ! Output :
 !
 !
 !------------------------------------------------------------------------
 ! DECLARATIONS
   IMPLICIT NONE
   TYPE( NODAL_STORAGE_2D ), intent(in)    :: dGStorage
   TYPE( BARO_MAPPEDGEOM_2D ), intent(in)    :: geometry
   TYPE( BAROTROPIC_PROPERTIES), intent(in)  :: dG2DProps
   TYPE( DG2D_SOLUTION ), intent(inout)      :: dGSol
   TYPE( RUN_PARAMS ), intent(in)            :: locParams
   TYPE( DG2D_RELAXATION_FIELD ), intent(in) :: locRelax
   real(prec), intent(in)                    :: tn
   integer, intent(in)                       :: nX, nY, nEq
   ! LOCAL
   real(prec) :: dxds, dxdp, dyds, dydp
   real(prec) :: h
   real(prec) :: xF(1:nEq), yF(1:nEq), sol(1:nEq)
   real(prec) :: pContFlux(0:nY,1:nEq), sContFlux(0:nX,1:nEq)
   real(prec) :: pContFluxDer(0:nY,1:nEq), sContFluxDer(0:nX,1:nEq)
   real(prec) :: fL(1:nEq), fR(1:nEq)
   integer    :: iX, iY, iEq




     do iY = 0,nY ! Loop over the y-points

        do iX = 0,nX ! Loop over the x-points

           ! Get the metric terms to calculate the contravariant flux
           CALL geometry % GET_METRICS_ATPOINT( iX, iY, dxds, dxdp, dyds, dydp )


           ! Get the solution at x,y
           CALL dgSol % GET_SOL_ATNODE( iX, iY, sol )

           ! Get the depth - needed for wave speed calculation 
           h = geometry % depth(iX,iY)


           ! Calculate the x and y fluxes
           xF = XSWFLUX(  tn, nEq, locParams % g, h, sol, locParams % MODEL_FORMULATION )

           yF = YSWFLUX(  tn, nEq, locParams % g, h, sol, locParams % MODEL_FORMULATION )

           !And now the contravariant flux
           sContFlux(iX,1:nEq) = dydp*xF(1:nEq) - dxdp*yF(1:nEq)


        enddo ! iX, Loop over the x-points

        ! Get the numerical fluxes at the boundaries (should be calculated before
        ! this routine is called )
        
        ! Get the "west" flux at iY
        CALL dgSol % GET_FLUXBOUND_ATPOINT( 4, iY, fL )

        ! Get the "east" flux at iY
        CALL dgSol % GET_FLUXBOUND_ATPOINT( 2, iY, fR )

        ! At this y-level, calculate the DG-advective derivative
         sContFluxDer(0:nX,1:nEq) = DG_SYSTEMDER(  nEq, nX, dGStorage % dMatX, dGStorage % qWeightX, &
                                                   fL, fR, sContFlux,  dGStorage % lagWest, &
                                                   dGStorage % lagEast  )

         
         do iX = 0, nX ! Loop over the x-points
            do iEq = 1, nEq  ! Loop over the number of equations

               dGSol % solTendency(iX, iY, iEq) = -sContFluxDer(iX,iEq)

            enddo ! Loop over the number of equations
         enddo ! Loop over the x-points

      
     enddo ! iY  


     do iX = 0,nX ! Loop over the x-points

        do iY = 0,nY ! Loop over the y-points

           ! Get the metric terms to calculate the contravariant flux
           CALL geometry % GET_METRICS_ATPOINT( iX, iY, dxds, dxdp, dyds, dydp )


           CALL dgSol % GET_SOL_ATNODE( iX, iY, sol )

           ! Get the depth - needed for wave speed calculation 
           h = geometry % depth(iX,iY)


           ! Calculate the x and y fluxes
           xF = XSWFLUX(  tn, nEq, locParams % g, h, sol, locParams % MODEL_FORMULATION ) 

           yF = YSWFLUX(  tn, nEq, locParams % g, h, sol, locParams % MODEL_FORMULATION )


           !And now the contravariant flux

           pContFlux(iY,1:nEq) = -dyds*xF(1:nEq) + dxds*yF(1:nEq)

        enddo ! iY, Loop over the y-points

        ! Get the numerical fluxes at the boundaries (should be calculated before
        ! this routine is called )
        
        ! Get the "south" flux at iX
        CALL dgSol % GET_FLUXBOUND_ATPOINT( 1, iX, fL )

        ! Get the "north" flux at iX
        CALL dgSol % GET_FLUXBOUND_ATPOINT( 3, iX, fR )


        ! At this x-level, calculate the y-DG-advective derivative
         pContFluxDer(0:nY,1:nEq) = DG_SYSTEMDER(  nEq, nY, dGStorage % dMatY, dGStorage % qWeightY, &
                                                   fL, fR, pContFlux, dGStorage % lagSouth, &
                                                   dGStorage % lagNorth  )

         
         do iY = 0, nY ! Loop over the y-points
            do iEq = 1, nEq  ! Loop over the number of equations

                dGSol % solTendency(iX, iY, iEq) = ( dGSol % solTendency(iX, iY, iEq)- &
                                                     pContFluxDer(iY,iEq) )/geometry % J(iX,iY)

            enddo ! Loop over the number of equations
         enddo ! Loop over the x-points

      
     enddo ! iY

     ! COMPUTE THE SOURCE TERMS
     do iX = 0, nX  ! Loop over the x-points
        do iY = 0, nY ! Loop over the y-points


           dGSol % solTendency(iX, iY, 1:nEq) = dGSol % solTendency(iX, iY, 1:nEq)+&
                                                SWSOURCE( tn, nEq, &
                                                          dG2DProps % vorticity(iX,iY),&
                                                          geometry % x(iX,iY), geometry % y(iX,iY),&
                                                          locParams % g, geometry % depth(iX,iY),&
                                                          geometry % dDepthdX(iX,iY), &
                                                          geometry % dDepthdY(iX,iY), &
                                                          dGSol % solInterior(iX, iY, 1:nEq),&
                                                          locParams % MODEL_FORMULATION, &
                                                          locParams % lDrag,&
                                                          locRelax % solField(iX,iY,1:nEq),&
                                                          locRelax % rDrag(iX,iY) )

        enddo ! iY, Loop over the y-points
     enddo ! iX, Loop over the x-points



 END SUBROUTINE MAPPED_DGSW_TIMEDER
!
!
!
FUNCTION DG_SYSTEMDER(  nEq, nP, dMat, qWei, lFlux, rFlux, intFlux, lagLeft, lagRight  ) RESULT( tendency )
 ! FUNCTION DG_SYSTEMDER ( Discontinous Galerking time DERivative)
 ! 
 ! Description :
 !  
 !
 !    Subroutine Dependencies :
 !    
 !
 ! Input :
 !
 ! 
 ! Output :
 !
 !
 !------------------------------------------------------------------------
 ! DECLARATIONS
   IMPLICIT NONE
   integer    :: nP, nEq
   real(prec) :: dMat(0:nP,0:nP)
   real(prec) :: qWei(0:nP)
   real(prec) :: lFlux(1:nEq)
   real(prec) :: rFlux(1:nEq)
   real(prec) :: intFlux(0:nP,1:nEq)
   real(prec) :: lagLeft(0:nP), lagRight(0:nP)
   real(prec) :: tendency(0:nP,1:nEq)
   ! Local 
   integer :: iEqn, kP

      tendency = 0.0_prec
      
      ! Apply the derivative operator
      do iEqn = 1, nEq

         tendency(0:nP,iEqn) = MATVECMUL( dMat, intflux(0:nP,iEqn),nP, nP )


      enddo


      do kP = 0, nP ! Loop over the quadrature points

         do iEqn = 1, nEq  ! Loop over the equations

            tendency(kP,iEqn) = tendency(kP,iEqn) + ( rFlux(iEqn)*lagRight(kP) + lFlux(iEqn)*lagLeft(kP) )/qWei(kP)

         enddo ! iEqn, cycle over the equations

      enddo ! kP, loop over the quadrature points

      
                                                  

 END FUNCTION DG_SYSTEMDER
!
!
!
 SUBROUTINE CALC_SWSOL_AT_BOUNDS( dGSol, dGStorage ) 
 ! SUBOUTINE CALC_SWSOL_AT_BOUNDS
 !  Desription:
 !   On input myMappedNDG % DGSol % solInterior contains updated interior solution 
 !   values. On output myMappedNDG % DGSol % solBound is updated by performing vector
 !   dot product between myMappedNDG % DGSol % solInterior and myMappedNDG % DGStor % lag*
 !   where lag* is the lagrange interpolating polynomial evaluated at the appropriate
 !   boundary.
 !
 !  Procedure dependencies:
 !   FUNCTION INTERP_TO_POINT_1D, MODULE LAGRANGE_1D.f90
 !
 !
 !  Input :
 ! 
 !    
 !
 !  Output :
 !   
 !
 !------------------------------------------------------------------------
 ! DECLARATIONS
   IMPLICIT NONE
   TYPE( DG2D_SOLUTION ), intent(inout)   :: dGSol
   TYPE( NODAL_STORAGE_2D ), intent(in) :: dGStorage
   ! LOCAL
   integer :: nEq, nX, nY
   integer :: iX, iY, iEq
   


       CALL dgSol % GET_NEQN( nEq )

       CALL dgStorage % GET_N_NODES( nX, nY )

       do iY = 0, nY ! Loop over y - points
       
          do iEq = 1, nEq ! Loop over the number of equations (-1 to exclude bathymetry)

             ! Setting the "west boundary"
              dgSol % solBound(iY,iEq,4) = &
                           DOT_PRODUCT( dgStorage % lagWest, &
                                      dgSol % solInterior(:,iY,iEq) )

            ! Setting the "east boundary"
              dgSol % solBound(iY,iEq,2) = &
                           DOT_PRODUCT( dgStorage % lagEast, &
                                      dgSol % solInterior(:,iY,iEq) )

          enddo ! iEq, loop over the number of equations

       enddo ! iY, loop over y-points
         

       do iX = 0, nX ! Loop over x - points
       
          do iEq = 1, nEq ! Loop over the number of equations (-1 to exclude bathymetry)

             ! Setting the "south boundary"
             dgSol % solBound(iX,iEq,1)  = &
                           DOT_PRODUCT( dgStorage % lagSouth, &
                                      dgSol % solInterior(iX,:,iEq) )

            ! Setting the "north boundary"
             dgSol % solBound(iX,iEq,3)  = &
                           DOT_PRODUCT( dgStorage % lagNorth, &
                                      dgSol % solInterior(iX,:,iEq) )

          enddo ! iEq, loop over the number of equations

       enddo ! iX, loop over x-points


 END SUBROUTINE CALC_SWSOL_AT_BOUNDS
!
!
!
 SUBROUTINE UPDATE_GHOSTS( myDGSEM, t, nDeg, nEl, nGhosts, nEq, theMPICOMM, myRank ) 
 ! SUBOUTINE UPDATE_GHOSTS
 !  Desription:
 !   
 !
 !  Procedure dependencies:
 !   FUNCTION INTERP_TO_POINT_1D, MODULE LAGRANGE_1D.f90
 !
 !
 !  Input :
 ! 
 !    
 !
 !  Output :
 !   
 !
 !------------------------------------------------------------------------
 ! DECLARATIONS
   IMPLICIT NONE
   INCLUDE 'mpif.h'
   CLASS( DGSEM_BAROTROPIC ), intent(inout) :: myDGSEM
   real(prec), intent(in)                   :: t 
   integer, intent(in)                      :: nGhosts, nDeg, nEl, nEq
   integer, intent(in)                      :: theMPICOMM, myRank
   ! LOCAL
   integer :: iGhost, iNode, e1, e2, s1, s2
   integer :: source, dest, locElID, iEq
   real(prec) :: nHat(1:2), nHatLength, intState(1:nEq),x ,y , h
   real(prec) :: boundData(0:nDeg,1:nEq), incomingData(0:nDeg,1:nEq) 
   integer :: iError, tag, MPI_PREC
   integer :: stat(MPI_STATUS_SIZE) 
   integer :: formulation
   
   if( prec == dp )then
     MPI_PREC = MPI_DOUBLE
   else
     MPI_PREC = MPI_FLOAT
   endif


      ! Loop over the ghost points
      
         do iGhost = 1, nGhosts
         
            e1 = myDGSEM % ghostEl(iGhost) % intElemID ! e1 is a global ID

            e1 = myDGSEM % mesh % glob2LocElemMap(e1,3) ! e1 is now a local ID

            s1 = myDGSEM % ghostEl(iGhost) % intEdgeID

            e2 = myDGSEM % ghostEl(iGhost) % extElemID ! external is a global ID or a boundary condition flag
              

            if( e2 <= NO_NORMAL_FLOW ) then ! then this truly is a physical boundary,

               do iNode = 0, nDeg

                  CALL myDGSEM % mesh % elements(e1) % geometry % GET_NHAT_ATPOINT( s1, iNode, nHat, nHatLength ) ! Get nHat

                  intState(1:nEq) = myDGSEM % sol(e1) % solBound(iNode,1:nEq,s1)

                  ! Get the boundary point locations
                  CALL myDGSEM % mesh % elements(e1) % geometry % GET_BOUNDARY_CURVE_ATPOINT( x, y, s1, iNode )

                  !
                  h = myDGSEM % mesh % elements(e1) % geometry % dbound( iNode, s1 )

                  myDGSEM % ghostEl(ighost) % extState(iNode,1:nEq) = GET_EXTERNAL_STATE( nEq, nHat, x, y, h, t,&
                                                                                          e2, intState, myDGSEM % params)

               enddo


            ! for parallel situations the possibility exists for e2 > 0. In this case
            ! a lookup table is used which relates the global element ID to a process ID
            ! and local element ID. from this information, the boundary data is exchanged
            else

               ! e2 contains the global element ID of the neighboring element
               ! we retrieve the process ID and local element ID from the glob2loc list
               
               dest  = myDGSEM % mesh % glob2LocElemMap(e2,2) ! process ID which has the element where we wish to 
                                                          ! send our boundary data             
               
               ! Here, myRank process sends it's boundary information to the dest process
           
               ! grab the boundary data
               boundData = myDGSEM % sol(e1) % solBound(0:nDeg,1:nEq, s1) 
                
               ! the tag for the message will be the global element ID
               tag = myDGSEM % ghostEl(iGhost) % intElemID 

               CALL MPI_SEND(boundData, (nDeg+1)*nEq, MPI_PREC, dest, tag, theMPICOMM, iError)
!if(myRank==1) then
!print*, tag, dest
!endif
            
            endif

         enddo

         CALL MPI_BARRIER( theMPICOMM, iError )
!if(myRank==3) then
!print*, 'receiving'
!endif
     


            do iGhost = 1, nGhosts
         
               e1 = myDGSEM % ghostEl(iGhost) % intElemID

               s1 = myDGSEM % ghostEl(iGhost) % intEdgeID

               e2 = myDGSEM % ghostEl(iGhost) % extElemID


               if( e2 > NO_NORMAL_FLOW ) then ! then we exchange data

                  ! e2 contains the global element ID of the neighboring element
                  ! we retrieve the process ID and local element ID from the glob2loc list
              
                  source  = myDGSEM % mesh % glob2LocElemMap(e2,2) ! process ID which we exchange info with                
            
                
                  ! the tag for the message will be the secondary element ID
                  
                  ! This tag should match the tag in the send statement
                  tag = myDGSEM % ghostEl(iGhost) % extElemID
!if(myRank==0) then
!print*, tag, source
!endif
                               
                  CALL MPI_RECV( incomingData,& ! where to store data
                                 (nDeg+1)*nEq, & ! how much data
                                 MPI_PREC, & ! type of data
                                 source, tag,& ! source Process, message tag
                                 theMPICOMM, stat, iError )  !  communicator
           
                  
                  myDGSEM % ghostEl(iGhost) % extState(0:nDeg,1:nEq) = incomingData

               endif
            
            enddo

           ! CALL MPI_BARRIER( theMPICOMM, iError )
!call MPI_FINALIZE( theMPICOMM )
! stop
 END SUBROUTINE UPDATE_GHOSTS

!
!
!
FUNCTION RIEMANN_LAXFRIEDRICH( inState, outState, nEqn,h, nHat, locParams ) RESULT( numFlux )
 ! FUNCTION RIEMANN_LAXFRIEDRICH ( )
 ! 
 ! Description :
 !  A Riemann solver for computing the numerical flux associated with the
 !  linear shallow water equations
 !
 !    Subroutine Dependencies :
 !    (None)
 !
 ! Input :
 !
 ! 
 ! Output :
 !
 !
 !------------------------------------------------------------------------
 ! DECLARATIONS
   IMPLICIT NONE
   integer    :: nEqn
   real(prec) :: inState(1:nEqn)
   real(prec) :: outState(1:nEqn)
   real(prec) :: h
   real(prec) :: nHat(1:2) ! normal direction
   real(prec) :: numFlux(1:nEqn)
   TYPE( RUN_PARAMS ) :: locParams
   ! LOCAL
   real(prec) :: uL, uR, vL, vR, pL, pR, uOut, uIn
   real(prec) :: gamm1, gamm2, locG
   real(prec) :: jump(1:nEqn), aS(1:nEqn)
   real(prec) :: uNorm, fac


      jump = outState - inState
    
      locG = locParams % g
      uL = inState(1)
      vL = inState(2)
      pL = inState(3)

      uR = outState(1)
      vR = outState(2)
      pR = outState(3)


      uOut =  uR*nHat(1) + vR*nHat(2)

      uIn = uL*nHat(1) + vL*nHat(2)


      if( locParams % MODEL_FORMULATION == LINEAR ) then

         fac = sqrt(locG*h)

         numFlux(1) = 0.5_prec*( locG*(pR + pL) + fac*(uIn-uOut) )*nHat(1)
         numFlux(2) = 0.5_prec*( locG*(pR + pL) + fac*(uIn-uOut) )*nHat(2)
         numFlux(3) = 0.5_prec*( h*(uIn + uOut) + fac*(pL - pR) )

      elseif( locParams % MODEL_FORMULATION == SKEW_SYMMETRIC )then
       
      
         fac = max( abs(0.25_prec*( uOut + sqrt( uOut*uOut + 16.0_prec*locg*(pR + h) ) )), &
                    abs(0.25_prec*( uOut - sqrt( uOut*uOut + 16.0_prec*locg*(pR + h) ) )), &
                    abs(0.25_prec*( uIn + sqrt( uIn*uIn + 16.0_prec*locg*(pL + h) ) )), &
                    abs(0.25_prec*( uIn - sqrt( uIn*uIn + 16.0_prec*locg*(pL + h) ) )) )
      

         aS(1) = ( 0.5_prec*(uR*uR + vR*vR) + locg*pR + &
                   0.5_prec*(uL*uL + vL*vL) + locg*pL )*nHat(1)

         aS(2) = ( 0.5_prec*(uR*uR + vR*vR) + locg*pR + &
                   0.5_prec*(uL*uL + vL*vL) + locg*pL )*nHat(2)

         aS(3) = ( uIn*(pL + h) + uOut*(pR + h) )

         numFlux = 0.5_prec*( aS - fac*jump ) 

      else ! CONSERVATIVE

         fac = max( abs(0.5_prec*( uOut + sqrt( uOut*uOut + 2.0_prec*locg*(pR*pR*pR) ) )/pR ), &
                    abs(0.5_prec*( uOut - sqrt( uOut*uOut + 2.0_prec*locg*(pR*pR*pR) ) )/pR ), &
                    abs(0.5_prec*( uIn + sqrt( uIn*uIn + 2.0_prec*locg*(pL*pL*pL) ) )/pL ), &
                    abs(0.5_prec*( uIn - sqrt( uIn*uIn + 2.0_prec*locg*(pL*pL*pL) ) )/pL ) )

         aS(1) = (uOut/pR)*uR + 0.5_prec*locg*pR*pR*nHat(1) +&
                 (uIn/pL)*uL + 0.5_prec*locg*pL*pL*nHat(1)

         aS(2) = (uOut/pR)*vR + 0.5_prec*locg*pR*pR*nHat(2) +&
                 (uIn/pL)*vL + 0.5_prec*locg*pL*pL*nHat(2)

         aS(3) = uOut + uIn

         numFlux = 0.5_prec*( aS - fac*jump ) 

      endif     
      


 END FUNCTION RIEMANN_LAXFRIEDRICH
!
!
! 
 FUNCTION XSWFLUX(  tn, nEqn, g, h, solAtX, formulation ) RESULT( fx )
 ! FUNCTION XSWFLUX
 ! 
 ! Description :
 !  
 !
 !    Subroutine Dependencies :
 !    
 !
 ! Input :
 !
 ! 
 ! Output :
 !
 !
 !------------------------------------------------------------------------
 ! DECLARATIONS
   IMPLICIT NONE
   integer    :: nEqn, formulation
   real(prec) :: h, g
   real(prec) :: solAtX(1:nEqn)
   real(prec) :: tn
   real(prec) :: fx(1:nEqn)
   real(prec) :: KE

 

      if( formulation == LINEAR ) then

         fx(1) = g*solAtX(3)

         fx(2) = ZERO

         fx(3) = h*solAtX(1)
         
      elseif( formulation == SKEW_SYMMETRIC ) then

         KE = (solAtX(1)*solAtX(1) + solAtX(2)*solAtX(2))*0.5_prec
     
         fx(1) = g*solAtX(3) + KE  !u-momentum flux =  barotropic pressure + kinetic energy
    
         fx(2) = ZERO  !v-momentum flux =0

         fx(3) = (h + solAtX(3))*solAtX(1) ! pressure flux (c^2 + P)*u

      else

         fx(1) = solAtX(1)*solAtX(1)/(solAtX(3)) + g*0.5_prec*solAtX(3)*solAtX(3)

         fx(2) = solAtX(1)*solAtX(2)/(solAtX(3))

         fx(3) = solAtX(1)

      endif
                                                  

 END FUNCTION XSWFLUX
!
!
!
 FUNCTION YSWFLUX(  tn, nEqn, g, h, solAtX, formulation ) RESULT( fy )
 ! FUNCTION YSWFLUX
 ! 
 ! Description :
 !  Computes the y-flux for the nonlinear shallow water equation in 1-D
 !  the Discontinuous Galerkin approximation
 !
 !    Subroutine Dependencies :
 !    
 !
 ! Input :
 !
 ! 
 ! Output :
 !
 !
 !------------------------------------------------------------------------
 ! DECLARATIONS
   IMPLICIT NONE
   integer    :: nEqn, formulation
   real(prec) :: h, g 
   real(prec) :: solAtX(1:nEqn)
   real(prec) :: tn
   real(prec) :: fy(1:nEqn)
   real(prec) :: KE


      if( formulation == LINEAR ) then

         fy(1) = ZERO

         fy(2) = g*solAtX(3)

         fy(3) = h*solAtX(2)
         
      elseif( formulation == SKEW_SYMMETRIC ) then

         ! Compute the kinetic energy
         KE = (solAtX(1)*solAtX(1) + solAtX(2)*solAtX(2))*0.5_prec
     
         fy(1) = 0.0_prec !u-momentum flux = 0
    
         fy(2) = g*solAtX(3) + KE !v-momentum flux = barotropic pressure + kinetic energy

         fy(3) = (h + solAtX(3))*solAtX(2) ! pressure flux (h + eta)*v

      else

         fy(1) = solAtX(2)*solAtX(1)/(solAtX(3)) 

         fy(2) = solAtX(2)*solAtX(2)/(solAtX(3)) + g*0.5_prec*solAtX(3)*solAtX(3)

         fy(3) = solAtX(2)

      endif
     

                                                  

 END FUNCTION YSWFLUX
!
!
!
 FUNCTION SWSOURCE(  tn, nEqn, vorticity, x, y, g, h, hx, hy, solAtX, formulation, lDrag, rField, rDrag ) RESULT( q )
 ! FUNCTION SWSOURCE ( Shallow Water Source terms)
 ! 
 ! Description :
 !  Computes source terms in the Shallow water equations, including :
 !      Vorticity times hat{z} cross U
 !      Net Evap minus precipitation
 ! 
 !   ** IF MODEL_FORMULATION == CONSERVATIVE, 
 !      
 ! 
 !    Subroutine Dependencies :
 !    
 !
 ! Input :
 !
 ! 
 ! Output :
 !
 !
 !------------------------------------------------------------------------
 ! DECLARATIONS
   IMPLICIT NONE
   integer    :: nEqn, formulation
   real(prec) :: vorticity, g, lDrag
   real(prec) :: solAtX(1:nEqn)
   real(prec) :: tn, x, y, h, hx, hy
   real(prec) :: q(1:nEqn)
   real(prec) :: rField(1:nEqn), rDrag ! relaxation field and associated "drag"
   ! LOCAL
   real(prec) :: omega, x0, y0, r, amp 
  

     if( formulation == SKEW_SYMMETRIC .OR. formulation == LINEAR ) then


           q(1) = ZERO! vorticity*solAtX(2) - lDrag*solAtX(1)/h - &
                  !(solAtX(1) - rField(1))*rDrag !u-momentum source = f*v, or (f+zeta)*v 
    
           q(2) = ZERO!-vorticity*solAtX(1) - lDrag*solAtX(2)/h -&
                   !(solAtX(2) - rField(2))*rDrag !v-momentum source = -f*u, or -(f+zeta)*u
 
           q(3) = ZERO  ! pressure source    

     else

           q(1) = vorticity*solAtX(2) + g*solAtX(3)*hx - &
                  lDrag*solAtX(1) -& !
                  (solAtX(1) - rField(1))*rDrag !u-momentum source = f*v, or (f+zeta)*v 
    
           q(2) = -vorticity*solAtX(1) + g*solAtX(3)*hy -&
                   lDrag*solAtX(2) -&
                   (solAtX(2) - rField(2))*rDrag !v-momentum source = -f*u, or -(f+zeta)*u
 
           q(3) = ZERO  ! pressure source   

     endif
    

 END FUNCTION SWSOURCE
!
!
!
 FUNCTION GET_EXTERNAL_STATE( nEqn, nHat, x, y, h, t, bcFlag, intState, locParams) RESULT( extState )
 ! S/R GET_EXTERNAL_STATE
 ! 
 ! Description :
 !  Computes the External state using ''whichBoundFlag'' and bcFlag
 !
 !    Subroutine Dependencies :
 !    (NONE)
 !
 ! Input :
 !
 ! 
 ! Output :
 !
 !
 !------------------------------------------------------------------------
 ! DECLARATIONS
   IMPLICIT NONE
   integer    :: nEqn
   real(prec) :: intState(1:nEqn)
   real(prec) :: x,y,t,h
   real(prec) :: nHat(1:2)
   integer    :: bcFlag, formulation
   real(prec) :: extState(1:nEqn)
   TYPE( RUN_PARAMS ) :: locParams
! LOCAL
   real(prec) :: planRot, phase

   real(prec) :: k, d, xr, xloc, yloc,c 
  

       if( bcFlag == NO_NORMAL_FLOW ) then

          extState(1) = (nHat(2)*nHat(2) -nHat(1)*nHat(1))*intState(1) - 2.0_prec*nHat(1)*nHat(2)*intState(2) ! u velocity
          extState(2) = (nHat(1)*nHat(1) -nHat(2)*nHat(2))*intState(2) - 2.0_prec*nHat(1)*nHat(2)*intState(1) ! v velocity
          extState(3) = intState(3) ! barotropic pressure

       elseif ( bcFlag == PRESCRIBED ) then
 


             k = sqrt(2.0_prec)/2.0_prec
             c = sqrt(locparams % g* h)
             

             extState(3) = locparams % A*sin( 2.0_prec*pi*locparams % omega *t)
           
             extState(1) = -locparams % A*sin( 2.0_prec*pi*locparams % omega *t)/c

             extState(2) = ZERO

             

       
       elseif( bcFlag == RADIATION ) then


             extState(1) = ZERO ! 0.5_prec*(intState(1) + (nHat(2)*nHat(2) -nHat(1)*nHat(1))*intState(1)) -&
                             ! 2.0_prec*nHat(1)*nHat(2)*intState(2) ! u velocity
             extState(2) = ZERO ! 0.5_prec*(intState(2) + (nHat(1)*nHat(1) -nHat(2)*nHat(2))*intState(2)) -&
                             ! 2.0_prec*nHat(1)*nHat(2)*intState(1) ! v velocity
             extState(3) = ZERO !-intState(3)


       else
         ! The other conditional will be implemented in the Spectral Element Method...
         ! The "neighbor" attribute will be added adn will be passed to this routine
         ! Using the neighbor attribute, one can extract the external state from the 
         ! neighboring element solution at the boundary point.

         print*, 'FUNCTION : GET_EXTERNAL_STATE : Invalid bcflag =', bcflag
         print*, 'Stopping.'
         STOP

       endif

        
 END FUNCTION GET_EXTERNAL_STATE
!
!
!=========================================================================!
!-------------------------- DIAGNOSTIC ROUTINES --------------------------!
!=========================================================================!
!
! 
 SUBROUTINE DIAGNOSE_INTEGRATED_ENERGIES( myDGsem, totEnergy )
! S/R DIAGNOSE_INTEGRATED_ENERGIES
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
  CLASS( DGSEM_BAROTROPIC ), intent(in) :: myDGsem
  real(prec), intent(out)               :: totEnergy
  !LOCAL
  integer :: iX, iY, iEl
  real(prec) :: KE, PE, H, J
  character(len=5) :: zoneID


   ! open( unit=2, file= trim(filename)//'.curve', &
   !       form='formatted', status=trim(filestat) )

   
   totEnergy = 0.0_prec

    do iEl = 1, myDGsem % mesh % nElems


      do iY = 0, myDGSEM % dGStorage % nS

         do iX = 0, myDGSEM % dGStorage % nS

           if( myDGSEM % params % MODEL_FORMULATION == SKEW_SYMMETRIC ) then


              KE = 0.5_prec*(myDGsem % sol(iEl) % solInterior(iX,iY,1)*&
                             myDGsem % sol(iEl) % solInterior(iX,iY,1) + &
                             myDGsem % sol(iEl) % solInterior(iX,iY,2)*&
                             myDGsem % sol(iEl) % solInterior(iX,iY,2) )

              H = myDGsem % mesh % elements(iEl) % geometry % depth(iX,iY) + &
                  myDGsem % sol(iEl) % solInterior(iX,iY,3) 

              J = myDGSEM % mesh % elements(iEl) % geometry % J(iX,iY)
 
              PE = 0.5_prec*myDGsem % sol(iEl) % solInterior(iX,iY,3)*&
                   myDGsem % sol(iEl) % solInterior(iX,iY,3)*&
                   myDGSEM % params % g

             !(g*eta^2 + (h+eta)*KE) = total energy
             totEnergy  = totEnergy + ( PE + H*KE )*&
                          myDGSEM % dGStorage % qWeightX(iX)*&
                          myDGSEM % dGStorage % qWeightY(iY) 

          elseif( myDGSEM % params % MODEL_FORMULATION == CONSERVATIVE )then ! assumed conservative

             KE = 0.5_prec*(myDGsem % sol(iEl) % solInterior(iX,iY,1)*&
                             myDGsem % sol(iEl) % solInterior(iX,iY,1) + &
                             myDGsem % sol(iEl) % solInterior(iX,iY,2)*&
                             myDGsem % sol(iEl) % solInterior(iX,iY,2) )
              
             H = myDGsem % sol(iEl) % solInterior(iX,iY,3) -&
                 myDGsem % mesh % elements(iEl) % geometry % depth(iX,iY)

             PE = 0.5_prec*H*H*myDGSEM % params % g  
          
             H =  myDGsem % sol(iEl) % solInterior(iX,iY,3)
 
             J = myDGSEM % mesh % elements(iEl) % geometry % J(iX,iY)

             ! (g*eta^2 + (h+eta)*KE)
             totEnergy = totEnergy + (PE + KE/H)*J*&
                         myDGSEM % dGStorage % qWeightX(iX)*&
                         myDGSEM % dGStorage % qWeightY(iY)

          else

             totEnergy = ZERO

          endif
          

         enddo

      enddo

    enddo



 END SUBROUTINE DIAGNOSE_INTEGRATED_ENERGIES
!
!
!
!
!
!=========================================================================!
!-------------------------- FILE I/O ROUTINES ----------------------------!
!=========================================================================!
!
!
 SUBROUTINE WRITE_BAROTROPIC_TECPLOT( myDGsem, Tx, Ty, nOld, nPlot, plotInterp, filename )
! S/R WRITE_BAROTROPIC_TECPLOT
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
  CLASS( DGSEM_BAROTROPIC ), intent(in) :: myDGsem
  integer, intent(in)                   :: nOld, nPlot
  real(prec), intent(in)                :: Tx(0:nPlot, 0:nOld), Ty(0:nPlot,0:nOld) 
  TYPE( LAG_INTERP2D ), intent(in)      :: plotInterp
  character(*), intent(in)              :: filename
  !LOCAL
  real(prec)  :: x(0:nPlot,0:nPlot), y(0:nPlot,0:nPlot)!, depth(0:nPlot,0:nPlot)
  real(prec)  :: u(0:nPlot,0:nPlot)
  real(prec)  :: v(0:nPlot,0:nPlot)
  real(prec)  :: eta(0:nPlot,0:nPlot)
  real(prec)  :: dep(0:nPlot,0:nPlot), locG
  integer :: iX, iY, iZ, iEl
  character(len=5) :: zoneID


    open( unit=2, file= trim(filename)//'.tec', form='formatted',status='replace')

    write(2,*) 'VARIABLES = "X", "Y", "U", "V", "Eta" '
    
    locG = myDGSEM % params % g

    do iEl = 1, myDGsem % mesh % nElems

       ! Interpolate the solutions onto a uniform mesh
       CALL COARSE_TO_FINE_2D(myDGsem % dGStorage % interp,&
                                      myDGsem % mesh % elements(iEl) % geometry % x,&
                                      plotInterp, x, Tx, Ty)


       CALL COARSE_TO_FINE_2D(myDGsem % dGStorage % interp,&
                                      myDGsem % mesh % elements(iEl) % geometry% y,&
                                      plotInterp, y, Tx, Ty)

       CALL COARSE_TO_FINE_2D(myDGsem % dGStorage % interp,&
                                      myDGsem % sol(iEl) % solInterior(:,:,1),&
                                      plotInterp, u, Tx, Ty)
 
       CALL COARSE_TO_FINE_2D(myDGsem % dGStorage % interp,&
                                      myDGsem % sol(iEl) % solInterior(:,:,2),&
                                      plotInterp, v, Tx, Ty)

       CALL COARSE_TO_FINE_2D(myDGsem % dGStorage % interp,&
                                      myDGsem % sol(iEl) % solInterior(:,:,3),&
                                      plotInterp, eta, Tx, Ty)

       if( myDGSEM % params % MODEL_FORMULATION == CONSERVATIVE ) then

          CALL COARSE_TO_FINE_2D(myDGsem % dGStorage % interp,&
                                 myDGsem % mesh % elements(iEl) % geometry % depth,&
                                 plotInterp, dep, Tx, Ty)

       endif


        write(zoneID,'(I5.5)') myDGsem % mesh % elements(iEl) % globElID


        write(2,*) 'ZONE T="el'//trim(zoneID)//'", I=',nPlot+1,', J=', nPlot+1,',F=POINT'!, ZONETYPE=FEQUADRILATERAL'

        if( myDGSEM % params % MODEL_FORMULATION == SKEW_SYMMETRIC .OR. &
            myDGSEM % params % MODEL_FORMULATION == LINEAR )then

           do iY = 0, nPlot
              do iX = 0, nPlot

                 write (2,*)  x( iX, iY ), y( iX, iY ), u(iX,iY), v(iX,iY), eta(iX,iY)

              enddo
           enddo

        elseif( myDGSEM % params % MODEL_FORMULATION == CONSERVATIVE )then
 
           do iY = 0, nPlot
              do iX = 0, nPlot

                 write (2,*)  x( iX, iY ), y( iX, iY ), &
                              u(iX,iY)/(eta(iX,iY)),&
                              v(iX,iY)/(eta(iX,iY)),&
                              eta(iX,iY) - dep(iX,iY)

              enddo
           enddo

        endif


    enddo


    close(unit=2)

    

 END SUBROUTINE WRITE_BAROTROPIC_TECPLOT
!
!
!
 SUBROUTINE WRITE_BAROTROPIC_TECPLOT_ASIS( myDGsem, filename )
! S/R WRITE_BAROTROPIC_TECPLOT_ASIS
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
  CLASS( DGSEM_BAROTROPIC ), intent(in) :: myDGsem
  character(*), intent(in)              :: filename
  !LOCAL
  real(prec)  :: locG, x, y, u, v, p, h
  integer :: iX, iY, iZ, iEl, nX, nY
  character(len=5) :: zoneID


    open( unit=2, file= trim(filename)//'.tec', form='formatted', status='replace')

    write(2,*) 'VARIABLES = "X", "Y", "U", "V", "Eta" '
     
    locG = myDGSEM % params % g

    nX = myDGSEM % dGStorage % nS
    nY = myDGSEM % dGStorage % nP

    if( myDGSEM % params % MODEL_FORMULATION == SKEW_SYMMETRIC .OR. &
        myDGSEM % params % MODEL_FORMULATION == LINEAR )then


       do iEl = 1, myDGsem % mesh % nElems


          write(zoneID,'(I5.5)') myDGsem % mesh % elements(iEl) % globElID


          write(2,*) 'ZONE T="el'//trim(zoneID)//'", I=',nX+1,', J=', nY+1,',F=POINT'!, ZONETYPE=FEQUADRILATERAL'

        
          101 FORMAT(E13.7,1x,E13.7,1x,E13.7,1x,E13.7,1x,E13.7,1x)
          do iY = 0, nY
             do iX = 0, nX

                x = myDGSEM % mesh % elements(iEl) % geometry % x(iX,iY)
                y = myDGSEM % mesh % elements(iEl) % geometry % y(iX,iY)
                h = myDGSEM % mesh % elements(iEl) % geometry % depth(iX,iY)

                u = myDGSEM % sol(iEl) % solInterior(iX,iY,1)
                v = myDGSEM % sol(iEl) % solInterior(iX,iY,2)
                p = myDGSEM % sol(iEl) % solInterior(iX,iY,3)
                 

                write (2,*)  x, y, u, v, p

             enddo
          enddo

       enddo

     elseif( myDGSEM % params % MODEL_FORMULATION == CONSERVATIVE )then
 
        do iEl = 1, myDGsem % mesh % nElems
           write(zoneID,'(I5.5)') myDGsem % mesh % elements(iEl) % globElID


          write(2,*) 'ZONE T="el'//trim(zoneID)//'", I=',nX+1,', J=', nY+1,',F=POINT'!, ZONETYPE=FEQUADRILATERAL'


           do iY = 0, nY
              do iX = 0, nX

                 x = myDGSEM % mesh % elements(iEl) % geometry % x(iX,iY)
                 y = myDGSEM % mesh % elements(iEl) % geometry % y(iX,iY)
                 h = myDGSEM % mesh % elements(iEl) % geometry % depth(iX,iY)

                 u = myDGSEM % sol(iEl) % solInterior(iX,iY,1)
                 v = myDGSEM % sol(iEl) % solInterior(iX,iY,2)
                 p = myDGSEM % sol(iEl) % solInterior(iX,iY,3)
                 
               !  u = myDGSEM % mesh % elements(iEl) % geometry % dDepthdx(iX,iY)
               !  v = myDGSEM % mesh % elements(iEl) % geometry % dDepthdy(iX,iY)
                 write (2,*)  x, y, u/(p), v/(p), p-h
               ! write(2,*) x, y, u, v, h
              enddo
           enddo
       enddo

    endif


    close(unit=2)

    

 END SUBROUTINE WRITE_BAROTROPIC_TECPLOT_ASIS
!
!
!
 SUBROUTINE WRITE_BAROTROPIC_RELAXFIELDS_ASIS( myDGsem, filename )
! S/R WRITE_BAROTROPIC_RELAXFIELDS_ASIS
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
  CLASS( DGSEM_BAROTROPIC ), intent(in) :: myDGsem
  character(*), intent(in)              :: filename
  !LOCAL
  real(prec)  :: locG, x, y, u, v, p, h
  integer :: iX, iY, iZ, iEl, nX, nY
  character(len=5) :: zoneID


    open( unit=2, file= trim(filename)//'.tec', form='formatted', status='replace')

    write(2,*) 'VARIABLES = "X", "Y", "U", "V", "Eta", "InvTimeScale" '
     
    locG = myDGSEM % params % g

    nX = myDGSEM % dGStorage % nS
    nY = myDGSEM % dGStorage % nP

    if( myDGSEM % params % MODEL_FORMULATION == SKEW_SYMMETRIC .OR. &
        myDGSEM % params % MODEL_FORMULATION == LINEAR )then
       do iEl = 1, myDGsem % mesh % nElems


          write(zoneID,'(I5.5)') myDGsem % mesh % elements(iEl) % globElID


          write(2,*) 'ZONE T="el'//trim(zoneID)//'", I=',nX+1,', J=', nY+1,',F=POINT'!, ZONETYPE=FEQUADRILATERAL'

        
          101 FORMAT(E13.7,1x,E13.7,1x,E13.7,1x,E13.7,1x,E13.7,1x)
          do iY = 0, nY
             do iX = 0, nX

                x = myDGSEM % mesh % elements(iEl) % geometry % x(iX,iY)
                y = myDGSEM % mesh % elements(iEl) % geometry % y(iX,iY)
                h = myDGSEM % mesh % elements(iEl) % geometry % depth(iX,iY)

                u = myDGSEM % relax(iEl) % solField(iX,iY,1)
                v = myDGSEM % relax(iEl) % solField(iX,iY,2)
                p = myDGSEM % relax(iEl) % solField(iX,iY,3)
                 

                write (2,*)  x, y, u, v, p, myDGSEM % relax(iEl) % rDrag(iX,iY)

             enddo
          enddo

       enddo

     elseif( myDGSEM % params % MODEL_FORMULATION == CONSERVATIVE )then
 
        do iEl = 1, myDGsem % mesh % nElems
           write(zoneID,'(I5.5)') myDGsem % mesh % elements(iEl) % globElID


          write(2,*) 'ZONE T="el'//trim(zoneID)//'", I=',nX+1,', J=', nY+1,',F=POINT'!, ZONETYPE=FEQUADRILATERAL'


           do iY = 0, nY
              do iX = 0, nX

                 x = myDGSEM % mesh % elements(iEl) % geometry % x(iX,iY)
                 y = myDGSEM % mesh % elements(iEl) % geometry % y(iX,iY)
                 h = myDGSEM % mesh % elements(iEl) % geometry % depth(iX,iY)

                 u = myDGSEM % relax(iEl) % solField(iX,iY,1)
                 v = myDGSEM % relax(iEl) % solField(iX,iY,2)
                 p = myDGSEM % relax(iEl) % solField(iX,iY,3)
              
                 write (2,*)  x, y, u, v, p, myDGSEM % relax(iEl) % rDrag(iX,iY)
               
              enddo
           enddo
       enddo

    endif


    close(unit=2)

    

 END SUBROUTINE WRITE_BAROTROPIC_RELAXFIELDS_ASIS
!
!
!
 SUBROUTINE WRITE_BAROTROPIC_PICKUP( myDGSEM, filename, iter )
 ! S/R WRITE_BAROTROPIC_PICKUP
 !  Description :
 ! 
 !
 !    Subroutine dependencies :
 !    (NONE)
 !    
 !
 !  Input :
 !    
 !
 !
 ! 
 !------------------------------------------------------------------------
 IMPLICIT NONE
  CLASS( DGSEM_BAROTROPIC ), intent(in) :: myDGSEM
  character(*), intent(in)              :: filename
  integer, intent(in)                   :: iter
  ! LOCAL
  character(10) :: iterChar
  integer :: iEl, recStart, recEnd
  integer :: rStartRel, rEndRel


     write(iterChar,'(I10.10)') iter

     recStart = 0

     rStartRel = 0

     do iEl = 1, myDGSEM % mesh % nElems
        
        CALL myDGSEM % sol(iEl) % WRITE_PICKUP( myDGSEM % dgStorage % nS, &
                                                myDGSEM % dgStorage % nP, &
                                                myDGSEM % nEq, &
                                                recStart, recEnd, filename//'.'//iterChar ) 

        recStart = recEnd

        CALL myDGSEM % relax(iEl) % WRITE_PICKUP( myDGSEM % dgStorage % nS, &
                                                  myDGSEM % dgStorage % nP, &
                                                  myDGSEM % nEq, &
                                                  rStartRel, rEndRel, filename//'.'//iterChar ) 

        rStartRel = rEndRel

     enddo

 END SUBROUTINE WRITE_BAROTROPIC_PICKUP 
 
 END MODULE DGSEM_BAROTROPIC_CLASS



