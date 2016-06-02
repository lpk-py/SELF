PROGRAM BSW_Integrator


 USE ModelPrecision
 USE ModelFlags
 USE ConstantsDictionary

 USE AlongShelfClass
 USE BarotropicShelfWaves_Class

 IMPLICIT NONE

 TYPE( AlongShelf1D ) :: mybsw

 INTEGER :: iter0, nT, dFreq
 INTEGER :: iT
 REAL(prec) :: tn, deltaT
 CHARACTER(10) :: iterChar

    CALL mybsw % Build( )

    CALL mybsw % shelfwaves % IRAM( )

    CALL mybsw % shelfwaves % ComputeGrowthRates( )

    CALL mybsw % shelfwaves % WriteTecplot( )

    ! ////////////////////////////////// Time Stepping /////////////////////////////////////////// !
       
      iter0 = mybsw % shelfwaves % params % iterInit
      nT = mybsw % shelfwaves % params % nTimeSteps
      dFreq = mybsw % shelfwaves % params % dumpFreq
      deltaT = mybsw % shelfwaves % params % dt

      WRITE( iterChar, '(I10.10)') iter0   
      CALL mybsw % GlobalTimeDerivative( ZERO )   
      CALL mybsw % WriteTecplot( 'StreamFunction.'//iterChar )
      
      DO iT = iter0, iter0+nT-1 ! Loop over time-steps

          tn = real( iT, prec )*deltaT ! time at the beginning of this timestep          
          CALL mybsw % ForwardStepRK3( tn ) ! Forward Step

          IF( mod( iT+1, dFreq) == 0 )THEN

             WRITE(iterChar,'(I10.10)') iT+1
             CALL mybsw % WriteTecplot( 'StreamFunction.'//iterChar )

          ENDIF

       ENDDO
    
    CALL mybsw % WritePickup( iter0+nT )

    CALL mybsw % Trash( )

END PROGRAM BSW_Integrator
