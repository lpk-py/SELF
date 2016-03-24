PROGRAM Burgers_InitialConditions



! src/common/
USE ModelPrecision
! src/highend/burgers1D
USE BurgersParamsClass
USE BurgersClass

 IMPLICIT NONE

 TYPE( Burgers1D ) :: myBur

 INTEGER :: iEl
 
    CALL myBur % Build( 0 )
   
    CALL InitialCondition( myBur )

    CALL myBur % WritePickup( 0 ) 
    CALL myBur % WriteTecplot( 'Burgers.init' )
    
    CALL myBur % Trash( )

 CONTAINS
 
 SUBROUTINE InitialCondition( myDGSEM )
 ! S/R InitialCondition
 !
 !    This subroutine (in v2.1) constructs the intial conditions for the tracer fields
 !
 ! ==============================================================================================  !
 ! DECLARATIONS
   TYPE( Burgers1D ), INTENT(inout) :: myDGSEM
   ! Local
   INTEGER    :: iS, iEl, nS,  nEl
   REAL(prec) :: x
   REAL(prec) :: sol(0:myDGSEM % nS,1)

      CALL myDGSEM % dgStorage % GetNumberOfNodes( nS )
      nEl = myDGSEM % nElems
      
      sol = ZERO
      
      ! Use the tendency to calculate the velocity by invoking geostrophic balance  
      DO iEl = 1, nEl
         DO iS = 0, nS
            x = myDGSEM % x(iS,iEl)
            sol(iS,1) = -(TWO/10.0_prec)*x + ONE
         ENDDO
         CALL myDGSEM % SetSolution( iEl, sol )
      ENDDO
               
               
 END SUBROUTINE InitialCondition
!
!
! 
END PROGRAM Burgers_InitialConditions
