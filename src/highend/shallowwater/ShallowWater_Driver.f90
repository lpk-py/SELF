PROGRAM ShallowWater_Driver



! src/common/
USE ModelPrecision
! src/geom/
USE QuadMeshClass
! src/highend/shallowwater
USE SWParamsClass
USE ShallowWaterClass

 IMPLICIT NONE

 TYPE( ShallowWater ) :: mysw

 INTEGER :: iter0, nT, dFreq
 INTEGER :: iT
 REAL(prec) :: tn, deltaT
 
 CHARACTER(len=10) :: iterChar
 
    CALL mysw % Build( )
 
    ! ////////////////////////////////// Time Stepping /////////////////////////////////////////// !
       
      iter0 = mysw % params % iterInit
      nT = mysw % params % nTimeSteps
      dFreq = mysw % params % dumpFreq
      deltaT = mysw % params % dt

      WRITE( iterChar, '(I10.10)') iter0   
      CALL mysw % GlobalTimeDerivative( ZERO )   
      CALL mysw % WriteTecplot( 'ShallowWater.'//iterChar )
      
      DO iT = iter0, iter0+nT ! Loop over time-steps

          tn = real( iT, prec )*deltaT ! time at the beginning of this timestep          
          CALL mysw % ForwardStepRK3( tn ) ! Forward Step

          IF( mod( iT+1, dFreq) == 0 )THEN

             WRITE(iterChar,'(I10.10)') iT+1
             CALL mysw % WriteTecplot( 'ShallowWater.'//iterChar )

          ENDIF

       ENDDO
    
    CALL mysw % WritePickup( iT ) 
    CALL mysw % mesh % WriteTecplot( )
    
    CALL mysw % Trash( )


END PROGRAM ShallowWater_Driver
