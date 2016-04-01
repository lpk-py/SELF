PROGRAM HydrostaticPrimitives_InitialConditions



! src/common/
USE ModelPrecision
! src/geom/
USE StackedHexMeshClass
! src/highend/advection2d
USE HydrostaticParams_Class
USE HydrostaticPrimitivesClass

 IMPLICIT NONE

 TYPE( HydrostaticPrimitive ) :: myHyd
 
    CALL myHyd % Build( 0 )
   
    CALL InitialCondition( myHyd )
    
    CALL myHyd % GlobalTimeDerivative( ZERO, 0 ) 
    
    CALL myHyd % WritePickup( 0 ) 
    CALL myHyd % WriteTecplot( 'Hydrostatic.init' )
    CALL myHyd % mesh % WriteTecplot( )
    CALL myHyd % mesh % WriteMeshFile( GAUSS, LEGENDRE_BASIS )
    
    CALL myHyd % Trash( )

 CONTAINS
 

 SUBROUTINE InitialCondition( myDGSEM )
 ! S/R InitialCondition
 !
 !    This subroutine (in v2.1) constructs the intial conditions for the tracer fields
 !
 ! ==============================================================================================  !
 ! DECLARATIONS
   TYPE( HydrostaticPrimitive ), INTENT(inout) :: myDGSEM
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
                  sol(iS,iP,iQ,1) = -x ! convergent zonal flow
                  sol(iS,iP,iQ,3) = z
               ENDDO
            ENDDO
         ENDDO

         CALL myDGSEM % SetSolution( iEl, sol )
      ENDDO
               
               
 END SUBROUTINE InitialCondition
!
!
!
END PROGRAM HydrostaticPrimitives_InitialConditions
