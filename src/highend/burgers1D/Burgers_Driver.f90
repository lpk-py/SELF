PROGRAM Burgers_Driver



! src/common/
USE ModelPrecision
! src/highend/burgers1D
USE BurgersParamsClass
USE BurgersClass

 IMPLICIT NONE

 TYPE( Burgers1D ) :: myBur

 INTEGER :: iter0, nT, dFreq
 INTEGER :: iT
 REAL(prec) :: tn, deltaT
 
 CHARACTER(len=10) :: iterChar
 
    CALL myBur % Build( )
 
    ! ////////////////////////////////// Time Stepping /////////////////////////////////////////// !
       
      iter0 = myBur % params % iterInit
      nT = myBur % params % nTimeSteps
      dFreq = myBur % params % dumpFreq
      deltaT = myBur % params % dt

      WRITE( iterChar, '(I10.10)') iter0   
      CALL myBur % GlobalTimeDerivative( ZERO )   
      CALL myBur % WriteTecplot( 'Burgers.'//iterChar )
      
      DO iT = iter0, iter0+nT-1 ! Loop over time-steps

          tn = real( iT, prec )*deltaT ! time at the beginning of this timestep          
          CALL myBur % ForwardStepRK3( tn ) ! Forward Step

          IF( mod( iT+1, dFreq) == 0 )THEN

             WRITE(iterChar,'(I10.10)') iT+1
             CALL myBur % WriteTecplot( 'Burgers.'//iterChar )

          ENDIF

       ENDDO
    
    CALL myBur % WritePickup( iter0+nT ) 
    
    CALL myBur % Trash( )


END PROGRAM Burgers_Driver
