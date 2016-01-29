PROGRAM ShallowWater_InitialConditions



! src/common/
USE ModelPrecision
! src/geom/
USE QuadMeshClass
! src/highend/shallowwater
USE SWParamsClass
USE ShallowWaterClass

 IMPLICIT NONE

 TYPE( ShallowWater ) :: mysw

 INTEGER :: iEl, forceNoPickup
 
    CALL mysw % Build( forceNoPickup )
 

    CALL Bathymetry( mysw )   
    CALL InitialCondition( mysw )
    CALL mysw % GlobalTimeDerivative( ZERO )

    
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
   REAL(prec) :: x, y
   REAL(prec) :: h(0:myDGSEM % nS,0:myDGSEM % nP,1:3)
   REAL(prec) :: x0, L1, L2, H1, H2, Ly, y0, xc, dx, L
   
      x0 = myDGSEM % params % shelfCenter     
      L1 = myDGSEM % params % UpstreamWidth   ! L1
      L2 = myDGSEM % params % DownstreamWidth ! L2
      H1 = myDGSEM % params % shelfDepth      ! H1
      H2 = myDGSEM % params % offshoreDepth   ! H2
      Ly = myDGSEM % params % transitionScale ! Ly
      y0 = myDGSEM % params % yTransition
      
      dx = L2 - L1
      CALL myDGSEM % dgStorage % GetNumberOfNodes( nS, nP )
      nEl = myDGSEM % mesh % nElems
      h = ZERO
      
      DO iEl = 1, nEl
         DO iP = 0, nP
            DO iS = 0, nS
            
              CALL myDGSEM % mesh % GetPositionAtNode( iEl, x, y, iS, iP )
               xc = x0 - HALF*dx*( ONE + tanh( (y - y0)/Ly ) )
               L  = L1 + HALF*(L2-L1)*( ONE + tanh( (y-y0)/Ly ) )
               h(iS,iP,1) = H1 + HALF*(H2-H1)*( ONE + tanh( (x-xc)/L ) )
               
            ENDDO
         ENDDO
         
         CALL myDGSEM % SetBathymetry( iEl, h )
         CALL myDGSEM % CalculateBathymetryAtBoundaries( iEl )
         
      ENDDO
      
 END SUBROUTINE Bathymetry
!
!
!
 SUBROUTINE InitialCondition( myDGSEM )
 ! S/R InitialCondition
 !
 !    This subroutine (in v2.1) constructs a geostrophic jet.
 !
 ! ==============================================================================================  !
 ! DECLARATIONS
   TYPE( ShallowWater ), INTENT(inout) :: myDGSEM
   ! Local
   INTEGER    :: iS, iP, iEl, nS, nP, nEl
   REAL(prec) :: x, y
   REAL(prec) :: f(0:myDGSEM % nS,0:myDGSEM % nP)
   REAL(prec) :: sol(0:myDGSEM % nS,0:myDGSEM % nP,1:3)
   REAL(prec) :: dpdx, dpdy
   REAL(prec) :: f0, g, dn, vm, Lj, x0
   
      
      
      f0 = myDGSEM % params % f0
      g  = myDGSEM % params % g
      vm = myDGSEM % params % maxV
      Lj = myDGSEM % params % jetWidth
      x0 = myDGSEM % params % jetCenter
      dn = TWO*f0*vm*Lj/g
      
      CALL myDGSEM % dgStorage % GetNumberOfNodes( nS, nP )
      nEl = myDGSEM % mesh % nElems
      sol = ZERO
      
      DO iEl = 1, nEl
         DO iP = 0, nP
            DO iS = 0, nS
            
               CALL myDGSEM % mesh % GetPositionAtNode( iEl, x, y, iS, iP )
               f(iS,iP) = f0         
               ! Setting the dipole free surface height
               sol(iS,iP,3) = HALF*dn*( ONE + tanh( (x-x0)/Lj ) )
               
            ENDDO
         ENDDO
         
         CALL myDGSEM % SetSolution( iEl, sol )
         CALL myDGSEM % SetPlanetaryVorticity( iEl, f )
      ENDDO
      
      CALL myDGSEM % GlobalTimeDerivative( ZERO ) ! Calculate the pressure gradient
      
      DO iEl = 1, nEl
         DO iP = 0, nP
            DO iS = 0, nS
            
               CALL myDGSEM % GetTendencyAtNodeWithVarID( iEl, iS, iP, 1, dpdx )
               dpdx = -dpdx
               CALL myDGSEM % GetTendencyAtNodeWithVarID( iEl, iS, iP, 2, dpdy )
               dpdy = -dpdy
                       
               ! Setting the velocity via geostrophy
               sol(iS,iP,1) = -dpdy/f0 ! u
               sol(iS,iP,2) = dpdx/f0
               
            ENDDO
         ENDDO
         
         CALL myDGSEM % SetVelocity( iEl, sol(:,:,1), sol(:,:,2) )
         
         CALL myDGSEM % CalculateSolutionAtBoundaries( iEl ) 
     
         CALL myDGSEM % CalculateRelativeVorticity( iEl )
      ENDDO
      

               
 END SUBROUTINE InitialCondition
!
!
! 
 
END PROGRAM ShallowWater_InitialConditions
