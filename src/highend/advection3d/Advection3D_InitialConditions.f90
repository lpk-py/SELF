PROGRAM Advection3D_InitialConditions



! src/common/
USE ModelPrecision
! src/geom/
USE HexMeshClass
! src/highend/advection2d
USE AdvectionParamsClass
USE Advection3DClass

 IMPLICIT NONE

 TYPE( Advection ) :: myAdv

 INTEGER :: iEl
 
    CALL myAdv % Build( 0 )
 

    CALL RelaxationFields( myAdv )   
    CALL InitialCondition( myAdv )
    CALL InitVelocity( myAdv )
    
    
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
   INTEGER    :: iS, iP, iEl, nS, nP, nQ, nEl
   REAL(prec) :: x, y
   REAL(prec) :: sol(0:myDGSEM % nS,0:myDGSEM % nP,0:myDGSEM % nQ, 1:myDGSEM % nEq)

      CALL myDGSEM % dgStorage % GetNumberOfNodes( nS, nP, nQ )
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
   INTEGER    :: iS, iP, iQ, iEl, nS, nP, nQ, nEl
   REAL(prec) :: x, y, z
   REAL(prec) :: sol(0:myDGSEM % nS,0:myDGSEM % nP, 0:myDGSEM % nQ, 1:myDGSEM % nEq)

      CALL myDGSEM % dgStorage % GetNumberOfNodes( nS, nP, nQ )
      nEl = myDGSEM % mesh % nElems
      
      sol = ZERO
      
      ! Use the tendency to calculate the velocity by invoking geostrophic balance  
      DO iEl = 1, nEl
         DO iQ = 0, nQ
            DO iP = 0, nP
               DO iS = 0, nS

                  CALL myDGSEM % mesh % GetPositionAtNode( iEl, x, y, z, iS, iP, iQ )
                  sol(iS,iP,iQ,:) = exp( -(x**2 + y**2 + z**2)/(0.2_prec**2) )
               ENDDO
            ENDDO
         ENDDO

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
   INTEGER    :: iS, iP, iEl, nS, nP, nQ, nEl
   REAL(prec) :: x, y, ys
   REAL(prec) :: u(0:myDGSEM % nS,0:myDGSEM % nP, 0:myDGSEM % nQ)
   REAL(prec) :: v(0:myDGSEM % nS,0:myDGSEM % nP, 0:myDGSEM % nQ)
   REAL(prec) :: w(0:myDGSEM % nS,0:myDGSEM % nP, 0:myDGSEM % nQ)

      CALL myDGSEM % dgStorage % GetNumberOfNodes( nS, nP, nQ )
      nEl = myDGSEM % mesh % nElems

      ys = myDGSEM % params % yScale
      
      ! Use the tendency to calculate the velocity by invoking geostrophic balance  
      DO iEl = 1, nEl

         u = ONE
         v = ONE
         w = ONE
         CALL myDGSEM % SetVelocity( iEl, u, v, w )

      ENDDO
                              
 END SUBROUTINE InitVelocity

END PROGRAM Advection3D_InitialConditions
