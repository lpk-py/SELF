PROGRAM ConservativeShallowWater_InitialConditions



! src/common/
USE ModelPrecision
! src/geom/
USE QuadMeshClass
! src/highend/shallowwater
USE CSWParamsClass
USE ConservativeShallowWaterClass

 IMPLICIT NONE

 TYPE( ShallowWater ) :: mysw

 INTEGER :: iEl, forceNoPickup
 
    forceNoPickup = 1
    CALL mysw % Build( forceNoPickup )
 

    CALL Bathymetry( mysw )   
    CALL InitialCondition( mysw )

    
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
   REAL(prec) :: x, y, L, yc, dLx, dLy, xc, w, h0, h1
   REAL(prec) :: h(0:myDGSEM % nS,0:myDGSEM % nP,1:3)
   
      CALL myDGSEM % dgStorage % GetNumberOfNodes( nS, nP )
      nEl = myDGSEM % mesh % nElems
      
      h = ZERO
      
      L = myDGSEM % params % shelfwidth
      dLx = myDGSEM % params % dLx
      dLy = myDGSEM % params % dLy
      yc = myDGSEM % params % yC
      h0 = myDGSEM % params % h0
      h1 = myDGSEM % params % h1
      
      DO iEl = 1, nEl
         DO iP = 0, nP
            DO iS = 0, nS
            
               CALL myDGSEM % mesh % GetPositionAtNode( iEl, x, y, iS, iP )
               
               w  = L - HALF*dLx*( tanh( (y-yc)/dLy ) + ONE )
               xc = w - L
                
               h(iS,iP,1) = h0 + (h1-h0)*( tanh( (x-xc)/w ) )
               
            ENDDO
         ENDDO
         
         CALL myDGSEM % SetBathymetry( iEl, h )
      ENDDO
      
 END SUBROUTINE Bathymetry
!
!
!
 SUBROUTINE InitialCondition( myDGSEM )
 ! S/R InitialCondition
 !
 !    This subroutine (in v2.1) constructs a geostrophic dipole by setting the free surface height
 !    and calculating the velocity from geostrophy.
 !
 ! ==============================================================================================  !
 ! DECLARATIONS
   TYPE( ShallowWater ), INTENT(inout) :: myDGSEM
   ! Local
   INTEGER    :: iS, iP, iEl, nS, nP, nEl
   REAL(prec) :: f0, x, y
   REAL(prec) :: f(0:myDGSEM % nS,0:myDGSEM % nP)
   REAL(prec) :: sol(0:myDGSEM % nS,0:myDGSEM % nP,1:3)
   REAL(prec) :: dpdx, dpdy, h(1:3)
   REAL(prec) :: dEta, x0, L
   
      dEta = myDGSEM % params % dEta
      x0   = myDGSEM % params % jetCenter
      L    = myDGSEM % params % jetWidth
   
      CALL myDGSEM % dgStorage % GetNumberOfNodes( nS, nP )
      nEl = myDGSEM % mesh % nElems
      f0 = myDGSEM % params % f0
      
      sol = ZERO
      
      DO iEl = 1, nEl
         DO iP = 0, nP
            DO iS = 0, nS
            
               CALL myDGSEM % mesh % GetPositionAtNode( iEl, x, y, iS, iP )
               CALL myDGSEM % GetBathymetryAtNode( iEl, iS, iP, h )
               f(iS,iP) = f0         
               
               ! Setting the dipole free surface height
               sol(iS,iP,3) = h(1) +  HALF*dEta*( tanh( (x-x0)/L ) )
               
            ENDDO
         ENDDO
         
         CALL myDGSEM % SetSolution( iEl, sol )
         CALL myDGSEM % SetPlanetaryVorticity( iEl, f )
      ENDDO
      
      ! Calculating the tendency gives the pressure gradient
      CALL myDGSEM % GlobalTimeDerivative( ZERO )
      
      ! Use the tendency to calculate the velocity by invoking geostrophic balance  
      DO iEl = 1, nEl
         DO iP = 0, nP
            DO iS = 0, nS
               ! Get the (negative) pressure gradient
               CALL myDGSEM % GetTendencyAtNodeWithVarID( iEl, iS, iP, 1, dpdx )
               CALL myDGSEM % GetTendencyAtNodeWithVarID( iEl, iS, iP, 2, dpdy )
               ! Flip the sign on the negative pressure gradient
               dpdx = -dpdx
               dpdy = -dpdy
               
               ! Set the "zonal" velocity
               sol(iS,iP,1) = -dpdy/f0
               ! Set the "meridional" velocity
               sol(iS,iP,2) = dpdx/f0
            ENDDO
         ENDDO
         CALL myDGSEM % SetVelocity( iEl, sol(:,:,1), sol(:,:,2) )
      ENDDO
               
               
 END SUBROUTINE InitialCondition
!
!
! 
 
END PROGRAM ConservativeShallowWater_InitialConditions
