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
   REAL(prec) :: h(0:myDGSEM % nS,0:myDGSEM % nP,1:3), h0
   
      CALL myDGSEM % dgStorage % GetNumberOfNodes( nS, nP )
      nEl = myDGSEM % mesh % nElems
      h0 = myDGSEM % params % h0
      h = ZERO
      
      DO iEl = 1, nEl
         DO iP = 0, nP
            DO iS = 0, nS
   
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
 !    This subroutine (in v2.1) constructs a geostrophic dipole.
 !
 ! ==============================================================================================  !
 ! DECLARATIONS
   TYPE( ShallowWater ), INTENT(inout) :: myDGSEM
   ! Local
   INTEGER    :: iS, iP, iEl, nS, nP, nEl
   REAL(prec) :: f0, x, y
   REAL(prec) :: f(0:myDGSEM % nS,0:myDGSEM % nP)
   REAL(prec) :: sol(0:myDGSEM % nS,0:myDGSEM % nP,1:3)
   REAL(prec) :: dpdx, dpdy, g
   REAL(prec) :: eta0, x0, y0, x1, y1, L, delta, R0, R1, KE, u0, u1, B
   
      eta0 = myDGSEM % params % eta0
      x0   = myDGSEM % params % x0
      y0   = myDGSEM % params % y0
      x1   = myDGSEM % params % x1
      y1   = myDGSEM % params % y1
      delta = myDGSEM % params % delta
      u0 = myDGSEM % params % u0
      u1 = myDGSEM % params % u1
      R0 = myDGSEM % params % R0
      R1 = myDGSEM % params % R1
      
      
      
      CALL myDGSEM % dgStorage % GetNumberOfNodes( nS, nP )
      nEl = myDGSEM % mesh % nElems
      f0 = myDGSEM % params % f0
      g  = myDGSEM % params % g
      
      sol = ZERO
      
      DO iEl = 1, nEl
         DO iP = 0, nP
            DO iS = 0, nS
            
               CALL myDGSEM % mesh % GetPositionAtNode( iEl, x, y, iS, iP )
               f(iS,iP) = f0         
               ! Setting the dipole free surface height
               sol(iS,iP,3) = eta0*( exp( -HALF*((x-x0)**2 + (y-y0)**2)/R0**2 ) - &
                                     exp( -HALF*((x-x1)**2 + (y-y1)**2)/R0**2 ) )
               
            ENDDO
         ENDDO
         
         CALL myDGSEM % SetSolution( iEl, sol )
         CALL myDGSEM % SetPlanetaryVorticity( iEl, f )
      ENDDO
      
      CALL myDGSEM % GlobalTimeDerivative( ZERO )
      
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
         
         CALL myDGSEM % SetSolution( iEl, sol )
      ENDDO
      

               
 END SUBROUTINE InitialCondition
!
!
! 
 
END PROGRAM ShallowWater_InitialConditions
