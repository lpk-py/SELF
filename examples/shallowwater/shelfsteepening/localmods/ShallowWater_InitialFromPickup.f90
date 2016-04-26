PROGRAM ShallowWater_InitialFromPickup



! src/common/
USE ModelPrecision
! src/geom/
USE QuadMeshClass
! src/highend/shallowwater
USE SWParamsClass
USE ConservativeShallowWaterClass

 IMPLICIT NONE

 TYPE( ShallowWater ) :: mysw

    CALL mysw % Build( )
 

    CALL Bathymetry( mysw )   

    
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
 
END PROGRAM ShallowWater_InitialFromPickup
