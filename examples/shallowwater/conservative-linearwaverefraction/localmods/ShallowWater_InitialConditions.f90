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
   REAL(prec) :: x, y, xs, ys
   REAL(prec) :: h(0:myDGSEM % nS,0:myDGSEM % nP,1:3)
   
      xs = myDGSEM % params % xScale
      ys = myDGSEM % params % yScale

      CALL myDGSEM % dgStorage % GetNumberOfNodes( nS, nP )
      nEl = myDGSEM % mesh % nElems
      h = ZERO
      
      DO iEl = 1, nEl
         DO iP = 0, nP
            DO iS = 0, nS
            
              CALL myDGSEM % mesh % GetPositionAtNode( iEl, x, y, iS, iP )

               h(iS,iP,1) = 1000.0_prec*( ONE - 0.8_prec*exp(-( (x-HALF*xs)**2 + (y-HALF*ys)**2 )/((0.1_prec*xs)**2) ) )
               
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
   REAL(prec) :: x, y, xs, ys
   REAL(prec) :: f(0:myDGSEM % nS,0:myDGSEM % nP)
   REAL(prec) :: sol(0:myDGSEM % nS,0:myDGSEM % nP,1:3)
   REAL(prec) :: f0, h(1:3)
   
      xs = myDGSEM % params % xScale
      ys = myDGSEM % params % yScale
      f0 = myDGSEM % params % f0
      CALL myDGSEM % dgStorage % GetNumberOfNodes( nS, nP )
      nEl = myDGSEM % mesh % nElems
      sol = ZERO
      
      DO iEl = 1, nEl
         DO iP = 0, nP
            DO iS = 0, nS
            
               CALL myDGSEM % mesh % GetPositionAtNode( iEl, x, y, iS, iP )
               CALL myDGSEM % GetBathymetryAtNode( iEl, iS, iP, h )
               f(iS,iP) = f0         
               ! Setting the dipole free surface height
           
               sol(iS,iP,3) = h(1) + 0.4_prec*exp(-( (x-0.75_prec*xs)**2 + (y-HALF*ys)**2 )/((0.05_prec*xs)**2) ) 
               
            ENDDO
         ENDDO
         
         CALL myDGSEM % SetSolution( iEl, sol )
         CALL myDGSEM % SetPlanetaryVorticity( iEl, f )
      ENDDO
      
!      CALL myDGSEM % GlobalTimeDerivative( ZERO ) ! Calculate the pressure gradient
      
!      DO iEl = 1, nEl
!         DO iP = 0, nP
!            DO iS = 0, nS
            
!               CALL myDGSEM % GetTendencyAtNodeWithVarID( iEl, iS, iP, 1, dpdx )
!               dpdx = -dpdx
!               CALL myDGSEM % GetTendencyAtNodeWithVarID( iEl, iS, iP, 2, dpdy )
!               dpdy = -dpdy
                       
!               ! Setting the velocity via geostrophy
!               sol(iS,iP,1) = -dpdy/f0 ! u
!               sol(iS,iP,2) = dpdx/f0
               
!            ENDDO
!         ENDDO
         
!         CALL myDGSEM % SetVelocity( iEl, sol(:,:,1), sol(:,:,2) )
         
!         CALL myDGSEM % CalculateSolutionAtBoundaries( iEl ) 
     
!         CALL myDGSEM % CalculateRelativeVorticity( iEl )
!      ENDDO
      

               
 END SUBROUTINE InitialCondition
!
!
! 
 
END PROGRAM ShallowWater_InitialConditions
