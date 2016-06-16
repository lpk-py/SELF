PROGRAM ShallowWater_InitialConditions



! src/common/
USE ModelPrecision
! src/geom/
USE QuadMeshClass
! src/highend/shallowwater
USE SWParamsClass
USE ConservativeShallowWaterClass

 IMPLICIT NONE

 TYPE( ShallowWater ) :: mysw

 INTEGER :: forceNoPickup
 
    CALL mysw % Build( forceNoPickup )
 

    CALL Bathymetry( mysw )   
    CALL InitialCondition( mysw )
    CALL mysw % GlobalTimeDerivative( ZERO, 0 )

    
    CALL mysw % WritePickup( 0 ) 
    CALL mysw % WriteTecplot( 'ShallowWater.init' )
    CALL mysw % mesh % WriteTecplot( )
    
    CALL mysw % Trash( )

 CONTAINS
 
 SUBROUTINE Bathymetry( myDGSEM )
 ! S/R Bathymetry
 !
 !    This subroutine sets the bathymetry as a flat bottom with a unit depth
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   TYPE( ShallowWater ), INTENT(inout) :: myDGSEM
   ! Local
   INTEGER    :: iS, iP, iEl, nS, nP, nEl
   REAL(prec) :: x, y, h0, h1, L, dL, yc, Ly, Ls, xs
   REAL(prec) :: h(0:myDGSEM % nS,0:myDGSEM % nP,1:3)
   
      h0 = myDGSEM % params % hmin
      h1 = myDGSEM % params % hMax
      L  = myDGSEM % params % Lshelf
      dL = myDGSEM % params % dL
      yc = myDGSEM % params % steepeningCenter
      Ly = myDGSEM % params % steepeningZoneLength   

      CALL myDGSEM % dgStorage % GetNumberOfNodes( nS, nP )
      nEl = myDGSEM % mesh % nElems
      h = ZERO
      
      DO iEl = 1, nEl
         DO iP = 0, nP
            DO iS = 0, nS
            
              CALL myDGSEM % mesh % GetPositionAtNode( iEl, x, y, iS, iP )

               Ls = L - HALF*dL*( tanh((y-yc)/Ly) + ONE )
               xs = Ls - L

               h(iS,iP,1) = h0 + (h1-h0)*tanh( (x-xs)/Ls )
               
            ENDDO
         ENDDO
      
         
         CALL myDGSEM % SetBathymetry( iEl, h )
         CALL myDGSEM % CalculateBathymetryAtBoundaries( iEl )
         CALL myDGSEM % CalculateBathymetryGradient( iEl )
      ENDDO  
      
 END SUBROUTINE Bathymetry
!
!
!
 SUBROUTINE InitialCondition( myDGSEM )
 ! S/R InitialCondition
 !
 !    This subroutine (in v2.1) constructs a geostrophic dipole.
 !
 ! ==============================================================================================  !
 ! DECLARATIONS
   TYPE( ShallowWater ), INTENT(inout) :: myDGSEM
   ! Local
   INTEGER    :: iS, iP, iEl, nS, nP, nEl, bID, iEdge, e1, e2, s1
   REAL(prec) :: x, y, x0, y0, x1, y1, dn, L, Ls, ys, c0, c1
   REAL(prec) :: dpdx, dpdy, vmax, g
   REAL(prec) :: f(0:myDGSEM % nS,0:myDGSEM % nP)
   REAL(prec) :: sol(0:myDGSEM % nS,0:myDGSEM % nP,1:3)
   REAL(prec) :: rsol(0:myDGSEM % nS,0:myDGSEM % nP,1:3)
   REAL(prec) :: tScale(0:myDGSEM % nS, 0:myDGSEM % nP)
   REAL(prec) :: state(0:myDGSEM % nS, 1:3)
   REAL(prec) :: f0, h(1:3)
   
      vmax = myDGSEM % params % vMax
      x0   = myDGSEM % params % x0
      y0   = myDGSEM % params % y0
      x1   = myDGSEM % params % x1
      y1   = myDGSEM % params % y1
      L    = myDGSEM % params % L
      f0   = myDGSEM % params % f0
      g    = myDGSEM % params % g
      ys   = TWO*myDGSEM % params % yScale
      Ls   = myDGSEM % params % Lsponge

      dn = f0*vmax*L/g

      CALL myDGSEM % dgStorage % GetNumberOfNodes( nS, nP )
      nEl = myDGSEM % mesh % nElems
      sol = ZERO

      rsol = ZERO
      DO iEl = 1, nEl
         tScale = myDGSEM % params % rFacMax
         DO iP = 0, nP
            DO iS = 0, nS
            
               CALL myDGSEM % mesh % GetPositionAtNode( iEl, x, y, iS, iP )
               CALL myDGSEM % GetBathymetryAtNode( iEl, iS, iP, h )
               f(iS,iP) = f0   
               c0 = (x-x0)**2 + (y-y0)**2    
               c1 = (x-x1)**2 + (y-y1)**2   
               ! Setting the dipole free surface height
               sol(iS,iP,3) = dn*( exp( -HALF*c0/(L*L) ) - &
                                   exp( -HALF*c1/(L*L) ) )

               tScale(iS,iP) = ZERO!tScale(iS,iP)*( tanh( (y - ys)/Ls ) + ONE )

            ENDDO
         ENDDO
         CALL myDGSEM % SetRelaxationField( iEl, rsol )
         CALL myDGSEM % SetRelaxationTimeScale( iEl, tScale )
         CALL myDGSEM % SetSolution( iEl, sol )
         CALL myDGSEM % SetPlanetaryVorticity( iEl, f )
      ENDDO
      
      CALL myDGSEM % GlobalTimeDerivative( ZERO, 0 ) ! Calculate the pressure gradient
      
      DO iEl = 1, nEl
         DO iP = 0, nP
            DO iS = 0, nS
            
               CALL myDGSEM % GetTendencyAtNodeWithVarID( iEl, iS, iP, 1, dpdx )
               dpdx = -dpdx
               CALL myDGSEM % GetTendencyAtNodeWithVarID( iEl, iS, iP, 2, dpdy )
               dpdy = -dpdy
                       
               ! Setting the velocity via geostrophy
               sol(iS,iP,1) = -dpdy/f0
               sol(iS,iP,2) = dpdx/f0
               
            ENDDO
         ENDDO
         
         CALL myDGSEM % SetVelocity( iEl, sol(:,:,1), sol(:,:,2) )
         CALL myDGSEM % CalculateSolutionAtBoundaries( iEl )
      ENDDO
      
      ! Set the inflow boundary conditions. To do this, we cycle over the boundary edges, and check
      ! for the edges with the INFLOW boundary flag in the secondary element ID,'
      ! For these edges, we copy the boundary solution into the "externalState" attribute
               
!      DO bID = 1, myDGSEM % nBoundaryEdges

!         iEdge = myDGSEM % boundaryEdgeIDs(bID)
!         CALL myDGSEM % mesh % GetEdgeSecondaryElementID( iEdge, e2 )

!         IF( e2 == PRESCRIBED )THEN
!            CALL myDGSEM % mesh % GetEdgePrimaryElementID( iEdge, e1 )
!            CALL myDGSEM % mesh % GetEdgePrimaryElementSide( iEdge, s1 )
!            CALL myDGSEM % GetBoundarySolutionAtBoundary( e1, s1, state )

!            myDGSEM % prescribedState(:,:,bID) = state
!         ENDIF

!      ENDDO

 END SUBROUTINE InitialCondition
!
! 
 
END PROGRAM ShallowWater_InitialConditions
