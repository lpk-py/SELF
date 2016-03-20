PROGRAM Advection3D_Driver



! src/common/
USE ModelPrecision
USE Timing
! src/geom/
USE HexMeshClass
! src/highend/shallowwater
USE AdvectionParamsClass
USE Advection3DClass

 IMPLICIT NONE

 TYPE( Advection ) :: myAdv

 INTEGER :: iter0, nT, dFreq
 INTEGER :: iT
 REAL(prec) :: tn, deltaT
 
 CHARACTER(len=10) :: iterChar
 
    CALL myAdv % Build( )
 
    ! ////////////////////////////////// Time Stepping /////////////////////////////////////////// !
       
      iter0 = myAdv % params % iterInit
      nT = myAdv % params % nTimeSteps
      dFreq = myAdv % params % dumpFreq
      deltaT = myAdv % params % dt

      WRITE( iterChar, '(I10.10)') iter0   
      CALL myAdv % GlobalTimeDerivative( ZERO )   

      CALL myAdv % clocks % StartThisTimer( 0 )
         CALL myAdv % WriteTecplot( 'Advection.'//iterChar )
      CALL myAdv % clocks % StopThisTimer( 0 )
      CALL myAdv % clocks % AccumulateTimings( ) 


      DO iT = iter0, iter0+nT-1 ! Loop over time-steps

          tn = real( iT, prec )*deltaT ! time at the beginning of this timestep          
          CALL myAdv % ForwardStepRK3( tn ) ! Forward Step

          IF( mod( iT+1, dFreq) == 0 )THEN

             WRITE(iterChar,'(I10.10)') iT+1
             CALL myAdv % clocks % StartThisTimer( 0 )
                CALL myAdv % WriteTecplot( 'Advection.'//iterChar )
             CALL myAdv % clocks % StopThisTimer( 0 )
             CALL myAdv % clocks % AccumulateTimings( ) 

          ENDIF

       ENDDO
    
    CALL myAdv % WritePickup( iT ) 
    CALL myAdv % mesh % WriteTecplot( )
    
    CALL myAdv % Trash( )


END PROGRAM Advection3D_Driver
