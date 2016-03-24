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

 INTEGER :: iEl
 
    CALL mysw % Build( 0 )
 

    CALL Bathymetry( mysw )  
    CALL SetPlanetaryRotation( mysw )
    CALL SetMassSources( mysw ) 
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
      ENDDO
      
 END SUBROUTINE Bathymetry
!
!
!
 SUBROUTINE SetPlanetaryRotation( myDGSEM )
 ! S/R InitialCondition
 !
 !    
 !
 ! ==============================================================================================  !
 ! DECLARATIONS
   TYPE( ShallowWater ), INTENT(inout) :: myDGSEM
   ! Local
   INTEGER    :: iS, iP, iEl, nS, nP, nEl
   REAL(prec) :: x, y, x0, y0, Lx, Ly, f0, betaX, betaY
   REAL(prec) :: f(0:myDGSEM % nS,0:myDGSEM % nP)


      CALL myDGSEM % dgStorage % GetNumberOfNodes( nS, nP )
      nEl = myDGSEM % mesh % nElems
      f0  = myDGSEM % params % f0
      betaX = myDGSEM % params % betaX
      betaY = myDGSEM % params % betaY
      Lx = myDGSEM % params % xScale
      Ly = myDGSEM % params % yScale 

      ! In this example, the parameter "f0" is the planetary rotation at the center of the mesh
      ! Mesh center
      x0 = 0.5_prec*Lx
      y0 = 0.5_prec*Ly
      f = ZERO

      DO iEl = 1, nEl
         DO iP = 0, nP
            DO iS = 0, nS

               CALL myDGSEM % mesh % GetPositionAtNode( iEl, x, y, iS, iP )
               f(iS,iP) = f0 + betaY*(y-y0) + betaX*(x-x0)

            ENDDO
         ENDDO
         CALL myDGSEM % SetPlanetaryVorticity( iEl, f )
      ENDDO

 END SUBROUTINE SetPlanetaryRotation
!
!
!
 SUBROUTINE SetMassSources( myDGSEM )
 ! S/R InitialCondition
 !
 !    
 !
 ! ==============================================================================================  !
 ! DECLARATIONS
   TYPE( ShallowWater ), INTENT(inout) :: myDGSEM
   ! Local
   INTEGER    :: iS, iP, iEl, nS, nP, nEl
   REAL(prec) :: x, y, x0, y0, x1, y1, Lx, Ly, Ls, wEk
   REAL(prec) :: s(0:myDGSEM % nS,0:myDGSEM % nP,1:3)
   REAL(prec) :: st(0:myDGSEM % nS,0:myDGSEM % nP)


      CALL myDGSEM % dgStorage % GetNumberOfNodes( nS, nP )
      nEl = myDGSEM % mesh % nElems

      wEk = myDGSEM % params % wEk
      Lx = myDGSEM % params % xScale
      Ly = myDGSEM % params % yScale
      Ls = myDGSEM % params % sourceLFac*sqrt((Lx**2 + Ly**2)*0.5_prec)

      ! In this example, there is one source and one sink
      ! The source is set in the southeast corner of the mesh
      ! and the sink is set at the northeast corner of the mesh.
      ! For low rossby number, theory predicts a northward flowing western
      ! boundary current, provided there is an opposing torque to balance
      ! planetary vorticity advection in the wbc. For this model, bottom
      ! friction provides that sink.
      
      ! Source location
      x0 = 0.75_prec*Lx
      y0 = 0.25_prec*Ly
      ! Sink Location
      x1 = x0
      y1 = 0.75*Ly

      ! Set the "ramp-up" timescale as the inertial period
      st = 1.0_prec/(myDGSEM % params % f0)
      s = ZERO
      
      DO iEl = 1, nEl
         DO iP = 0, nP
            DO iS = 0, nS

               ! The sources and sinks are imposed as Gaussian bumps with half-width of Ls
               ! The units on the source or sink are length/time so that it represents an Ekman 
               ! pumping velocity. The geostrophic velocity scale in the interior of the basin is
               ! related to the ekman pumping scale via
               !
               !   U ~ g*wEk/(f0^2*Ls)
               !   For g ~ 10 m/s^2, f0 ~ 10^(-4) 1/s, Ls ~ 100 km
               !   U ~ wEK*10^(4)
               !
               !   Here, the Ekman pumping velocity, wEk = 10^(-6)
               !
               CALL myDGSEM % mesh % GetPositionAtNode( iEl, x, y, iS, iP )
               s(iS,iP,3) = wEk*( exp( -HALF*( (x-x0)**2 + (y-y0)**2 )/Ls**2  ) - &
                                  exp( -HALF*( (x-x1)**2 + (y-y1)**2 )/Ls**2  ) )

            ENDDO
         ENDDO
         CALL myDGSEM % SetRelaxationField( iEl, s )
         CALL myDGSEM % SetRelaxationTimeScale( iEl, st )
      ENDDO

 END SUBROUTINE SetMassSources 
!
!
!
 SUBROUTINE InitialCondition( myDGSEM )
 ! S/R InitialCondition
 !
 !    
 !
 ! ==============================================================================================  !
 ! DECLARATIONS
   TYPE( ShallowWater ), INTENT(inout) :: myDGSEM
   ! Local
   INTEGER    :: iS, iP, iEl, nS, nP, nEl
   REAL(prec) :: x, y
   REAL(prec) :: sol(0:myDGSEM % nS,0:myDGSEM % nP,1:3)


      CALL myDGSEM % dgStorage % GetNumberOfNodes( nS, nP )
      nEl = myDGSEM % mesh % nElems
      
      sol = ZERO
      
      DO iEl = 1, nEl
         CALL myDGSEM % SetSolution( iEl, sol )
      ENDDO
      
               
 END SUBROUTINE InitialCondition
!
!
! 
 
END PROGRAM ShallowWater_InitialConditions
