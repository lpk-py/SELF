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
   REAL(prec) :: x, y
   REAL(prec) :: sol(0:myDGSEM % nS,0:myDGSEM % nP,1:myDGSEM % nEq)

      CALL myDGSEM % dgStorage % GetNumberOfNodes( nS, nP )
      nEl = myDGSEM % mesh % nElems
      
      sol = ZERO
      
      ! Use the tendency to calculate the velocity by invoking geostrophic balance  
      DO iEl = 1, nEl
         CALL myDGSEM % SetRelaxationField( iEl, sol )
         CALL myDGSEM % SetRelaxationFactor( iEl, sol )
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
        !       sol(iS,iP,:) = ZERO!exp( -(x-0.2_prec)**2/(2.0_prec*0.001_prec) )
 
        !    ENDDO
        ! ENDDO

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
   REAL(prec) :: x, y, ys
   REAL(prec) :: u(0:myDGSEM % nS,0:myDGSEM % nP)
   REAL(prec) :: v(0:myDGSEM % nS,0:myDGSEM % nP)

      CALL myDGSEM % dgStorage % GetNumberOfNodes( nS, nP )
      nEl = myDGSEM % mesh % nElems

      ys = myDGSEM % params % yScale
      
      ! Use the tendency to calculate the velocity by invoking geostrophic balance  
      DO iEl = 1, nEl
         DO iP = 0, nP
            DO iS = 0, nS

               CALL myDGSEM % mesh % GetPositionAtNode( iEl, x, y, iS, iP )
               v(iS,iP) = ZERO 
               u(iS,iP) = tanh( (y-ys*HALF)/(0.1_prec*ys)  )
               
            ENDDO
         ENDDO
               CALL myDGSEM % SetVelocity( iEl, u, v )
      ENDDO
               
               
 END SUBROUTINE InitVelocity
END PROGRAM Advection_InitialConditions
