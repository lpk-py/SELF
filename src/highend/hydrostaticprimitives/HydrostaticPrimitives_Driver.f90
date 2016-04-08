PROGRAM HydrostaticPrimitives_Driver



! src/common/
USE ModelPrecision
USE Timing
! src/geom/
USE HexMeshClass
! src/highend/shallowwater
USE HydrostaticParams_Class
USE HydrostaticPrimitivesClass

 IMPLICIT NONE

 TYPE( HydrostaticPrimitive ) :: myHyd

 INTEGER :: iter0, nT, dFreq
 INTEGER :: iT
 REAL(prec) :: tn, deltaT

    
      CALL myHyd % Build( InitializeOnly = .FALSE. )

       iter0 = myHyd % params % iterInit
       nT = myHyd % params % nTimeSteps
       dFreq = myHyd % params % dumpFreq
       deltaT = myHyd % params % dt
   ! ////////////////////////////////// Time Stepping //////////////////////////////////////////// !

      CALL myHyd % clocks % StartThisTimer( fIOTimer )
         CALL myHyd % WriteTecplot( iter0 )
      CALL myHyd % clocks % StopThisTimer( fIOTimer )


      DO iT = iter0, iter0+nT-1 ! Loop over time-steps

          tn = real( iT, prec )*deltaT ! time at the beginning of this timestep          
          CALL myHyd % ForwardStepRK3( tn ) ! Forward Step

          IF( mod( iT+1, dFreq) == 0 )THEN

             CALL myHyd % clocks % StartThisTimer( fIOTimer )
                CALL myHyd % WriteTecplot( iT+1 )
             CALL myHyd % clocks % StopThisTimer( fIOTimer )

          ENDIF

       ENDDO

    ! //////////////////////////////////////////////////////////////////////////////////////////// !

    CALL myHyd % WritePickup( iter0+nT ) 
    CALL myHyd % mesh % WriteTecplot( )
    
    CALL myHyd % Trash( )


END PROGRAM HydrostaticPrimitives_Driver
