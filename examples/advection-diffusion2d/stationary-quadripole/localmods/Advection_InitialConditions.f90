PROGRAM Advection_InitialConditions



! src/common/
USE ModelPrecision
! src/geom/
USE QuadMeshClass
! src/highend/advection2d
USE AdvectionParamsClass
USE AdvectionClass

 IMPLICIT NONE

 TYPE( Advection ) :: myAdv

 INTEGER :: iEl
 
    CALL myAdv % Build( 0 )
 

    CALL RelaxationFields( myAdv )   
    CALL InitialCondition( myAdv )
    CALL InitVelocity( myAdv )
    
    DO iEl = 1, myAdv % mesh % nElems
       CALL myAdv % CalculateConcentrationGradient( iEl )
    ENDDO
    
    CALL myAdv % WritePickup( 0 ) 
    CALL myAdv % WriteTecplot( 'Advection.init' )
    CALL myAdv % mesh % WriteTecplot( )
    
    CALL myAdv % Trash( )

 CONTAINS
 
 SUBROUTINE RelaxationFields( myDGSEM )
 ! S/R RelaxationFields
 !
 !    This subroutine (in v2.1) constructs the intial conditions for the tracer fields
 !
 ! ==============================================================================================  !
 ! DECLARATIONS
   TYPE( Advection ), INTENT(inout) :: myDGSEM
   ! Local
   INTEGER    :: iS, iP, iEl, nS, nP, nEl
   REAL(prec) :: x, y, ys, x1, x2, x3, x4, y1, y2, y3, y4, delta
   REAL(prec) :: m1, m2, m3, m4, tScale
   REAL(prec) :: sol(0:myDGSEM % nS,0:myDGSEM % nP,1:myDGSEM % nEq)
   REAL(prec) :: rfac(0:myDGSEM % nS,0:myDGSEM % nP,1:myDGSEM % nEq)

      CALL myDGSEM % dgStorage % GetNumberOfNodes( nS, nP )
      nEl = myDGSEM % mesh % nElems
      
      sol = ZERO
      
      delta = 0.01_prec
      x1 = 0.25_prec
      x2 = 0.75_prec
      x3 = x2
      x4 = x1
      
      y1 = 0.25_prec
      y2 = y1
      y3 = 0.75_prec
      y4 = y3
      tScale = 0.01_prec

      ! Use the tendency to calculate the velocity by invoking geostrophic balance  
      DO iEl = 1, nEl
         DO iP = 0, nP
            DO iS = 0, nS

               CALL myDGSEM % mesh % GetPositionAtNode( iEl, x, y, iS, iP )
               
               m1 = exp( -HALF*( (x-x1-10.0_prec*delta)**2 + (y-y1)**2 )/(delta**2) )
               m2 = exp( -HALF*( (x-x2-10.0_prec*delta)**2 + (y-y2)**2 )/(delta**2) )
               m3 = exp( -HALF*( (x-x3-10.0_prec*delta)**2 + (y-y3)**2 )/(delta**2) )
               m4 = exp( -HALF*( (x-x4-10.0_prec*delta)**2 + (y-y4)**2 )/(delta**2) )

               sol(iS,iP,1) = ONE*tScale
               sol(iS,iP,2) = ONE*tScale
               sol(iS,iP,3) = ONE*tScale
               sol(iS,iP,4) = ONE*tScale

               rFac(iS,iP,1) = m1/tScale
               rFac(iS,iP,2) = m2/tScale
               rFac(iS,iP,3) = m3/tScale
               rFac(iS,iP,4) = m4/tScale
               
            ENDDO
         ENDDO
         CALL myDGSEM % SetRelaxationField( iEl, sol )
         CALL myDGSEM % SetRelaxationFactor( iEl, rFac )
      ENDDO
               
               
 END SUBROUTINE RelaxationFields
!
!
!
 SUBROUTINE InitialCondition( myDGSEM )
 ! S/R InitialCondition
 !
 !    This subroutine (in v2.1) constructs the intial conditions for the tracer fields
 !
 ! ==============================================================================================  !
 ! DECLARATIONS
   TYPE( Advection ), INTENT(inout) :: myDGSEM
   ! Local
   INTEGER    :: iS, iP, iEl, nS, nP, nEl
   REAL(prec) :: x, y
   REAL(prec) :: sol(0:myDGSEM % nS,0:myDGSEM % nP,1:myDGSEM % nEq)

      CALL myDGSEM % dgStorage % GetNumberOfNodes( nS, nP )
      nEl = myDGSEM % mesh % nElems
      
      sol = ZERO
      
      ! Use the tendency to calculate the velocity by invoking geostrophic balance  
      DO iEl = 1, nEl
        ! DO iP = 0, nP
        !    DO iS = 0, nS

        !       CALL myDGSEM % mesh % GetPositionAtNode( iEl, x, y, iS, iP )
        !       sol(iS,iP,:) = exp( -(x-0.2_prec)**2/(2.0_prec*0.001_prec) )
 
!            ENDDO
!         ENDDO

               CALL myDGSEM % SetSolution( iEl, sol )
      ENDDO
               
               
 END SUBROUTINE InitialCondition
!
!
! 
  SUBROUTINE InitVelocity( myDGSEM )
 ! S/R InitVelocity
 !
 !    This subroutine (in v2.1) constructs the intial conditions for the tracer fields
 !
 ! ==============================================================================================  !
 ! DECLARATIONS
   TYPE( Advection ), INTENT(inout) :: myDGSEM
   ! Local
   INTEGER    :: iS, iP, iEl, nS, nP, nEl
   REAL(prec) :: x, y, ys, x1, x2, x3, x4, y1, y2, y3, y4
   REAL(prec) :: theta1, theta2, theta3, theta4
   REAL(prec) :: m1, m2, m3, m4, delta
   REAL(prec) :: u(0:myDGSEM % nS,0:myDGSEM % nP)
   REAL(prec) :: v(0:myDGSEM % nS,0:myDGSEM % nP)

      CALL myDGSEM % dgStorage % GetNumberOfNodes( nS, nP )
      nEl = myDGSEM % mesh % nElems
      ys = myDGSEM % params % yScale

      delta = 0.1_prec
      x1 = 0.25_prec
      x2 = 0.75_prec
      x3 = x2
      x4 = x1
      
      y1 = 0.25_prec
      y2 = y1
      y3 = 0.75_prec
      y4 = y3
      
      ! Use the tendency to calculate the velocity by invoking geostrophic balance  
      DO iEl = 1, nEl
         DO iP = 0, nP
            DO iS = 0, nS

               CALL myDGSEM % mesh % GetPositionAtNode( iEl, x, y, iS, iP )
               theta1 = atan2( y-y1, x-x1 )
               theta2 = atan2( y-y2, x-x2 )
               theta3 = atan2( y-y3, x-x3 )
               theta4 = atan2( y-y4, x-x4 )
               
               m1 = exp( -HALF*( (x-x1)**2 + (y-y1)**2 )/(delta**2) )
               m2 = exp( -HALF*( (x-x2)**2 + (y-y2)**2 )/(delta**2) )
               m3 = exp( -HALF*( (x-x3)**2 + (y-y3)**2 )/(delta**2) )
               m4 = exp( -HALF*( (x-x4)**2 + (y-y4)**2 )/(delta**2) )

               u(iS,iP) = ( m1*sin(theta1) - m2*sin(theta2) + m3*sin(theta3) - m4*sin(theta4) )
               v(iS,iP) = ( -m1*cos(theta1) + m2*cos(theta2) - m3*cos(theta3) + m4*cos(theta4) )

               ! Add image velocity fields to satisfy no-normal flow
               m1 = exp( -HALF*( (x+x1)**2 + (y-y1)**2 )/(delta**2) )       ! Left boundary, lower
               m2 = exp( -HALF*( (x-(ONE+x2))**2 + (y-y2)**2 )/(delta**2) ) ! Right boundary, lower
               m3 = exp( -HALF*( (x-(ONE+x3))**2 + (y-y3)**2 )/(delta**2) ) ! Right boundary, upper
               m4 = exp( -HALF*( (x+x4)**2 + (y-y4)**2 )/(delta**2) )       ! Left boundary, upper
               
               u(iS,iP) = u(iS,iP) - ( m1*sin(theta1) - m2*sin(theta2) + m3*sin(theta3) - m4*sin(theta4) )
               v(iS,iP) = v(iS,iP) - ( -m1*cos(theta1) + m2*cos(theta2) - m3*cos(theta3) + m4*cos(theta4) )

               ! Add image velocity fields to satisfy no-normal flow
               m1 = exp( -HALF*( (x-x1)**2 + (y+y1)**2 )/(delta**2) ) ! Lower left
               m2 = exp( -HALF*( (x-x2)**2 + (y+y2)**2 )/(delta**2) ) ! Lower right
               m3 = exp( -HALF*( (x-x3)**2 + (y-(ONE+y3))**2 )/(delta**2) ) ! Upper right
               m4 = exp( -HALF*( (x-x4)**2 + (y-(ONE+y4))**2 )/(delta**2) ) ! Lower left
               
               u(iS,iP) = u(iS,iP) - ( m1*sin(theta1) - m2*sin(theta2) + m3*sin(theta3) - m4*sin(theta4) )
               v(iS,iP) = v(iS,iP) - ( -m1*cos(theta1) + m2*cos(theta2) - m3*cos(theta3) + m4*cos(theta4) )
 
            ENDDO
         ENDDO
               CALL myDGSEM % SetVelocity( iEl, u, v )
      ENDDO
               
               
 END SUBROUTINE InitVelocity
END PROGRAM Advection_InitialConditions
