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
   REAL(prec) :: x, y, xs, ys, h0
   REAL(prec) :: h(0:myDGSEM % nS,0:myDGSEM % nP,1:3)
   
      xs = myDGSEM % params % xScale
      ys = myDGSEM % params % yScale
      h0 = myDGSEM % params % h0
      CALL myDGSEM % dgStorage % GetNumberOfNodes( nS, nP )
      nEl = myDGSEM % mesh % nElems
      h = ZERO
      
      DO iEl = 1, nEl
         DO iP = 0, nP
            DO iS = 0, nS
            
              CALL myDGSEM % mesh % GetPositionAtNode( iEl, x, y, iS, iP )

               h(iS,iP,1) = h0
               
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
   REAL(prec) :: x, y, x0, y0, y1, dn, r0, r1, L
   REAL(prec) :: dpdx, dpdy
   REAL(prec) :: f(0:myDGSEM % nS,0:myDGSEM % nP)
   REAL(prec) :: sol(0:myDGSEM % nS,0:myDGSEM % nP,1:3)
   REAL(prec) :: tScale(0:myDGSEM % nS, 0:myDGSEM % nP)
   REAL(prec) :: f0, h(1:3)
   
      x0 = myDGSEM % params % x0
      y0 = myDGSEM % params % y0
      y1 = myDGSEM % params % y1
      dn = myDGSEM % params % eta0
      f0 = myDGSEM % params % f0
      CALL myDGSEM % dgStorage % GetNumberOfNodes( nS, nP )
      nEl = myDGSEM % mesh % nElems
      sol = ZERO
      L = myDGSEM %  params % L
      tScale = myDGSEM % params % tScale
      DO iEl = 1, nEl
         DO iP = 0, nP
            DO iS = 0, nS
            
               CALL myDGSEM % mesh % GetPositionAtNode( iEl, x, y, iS, iP )
               CALL myDGSEM % GetBathymetryAtNode( iEl, iS, iP, h )
               f(iS,iP) = f0         
               ! Setting the dipole free surface height
               r0 = HALF*( (x-x0)**2 + (y-y0)**2 )/(L**2)
               r1 = HALF*( (x-x0)**2 + (y-y1)**2 )/(L**2)
               sol(iS,iP,3) = h(1) + dn*( exp(-r0) )! - exp(-r1) ) 
               
            ENDDO
         ENDDO
         CALL myDGSEM % SetRelaxationField( iEl, sol )
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
               sol(iS,iP,1) = -dpdy/f0 ! u
               sol(iS,iP,2) = dpdx/f0
               
            ENDDO
         ENDDO
         
         CALL myDGSEM % SetVelocity( iEl, sol(:,:,1), sol(:,:,2) )
         
      ENDDO
      

               
 END SUBROUTINE InitialCondition
!
! 
 
END PROGRAM ShallowWater_InitialConditions
