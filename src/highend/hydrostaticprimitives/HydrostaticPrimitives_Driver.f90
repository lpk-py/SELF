! HydrostaticPrimitives_Driver.f90
! 
! Copyright 2015 Joseph Schoonover <schoonover.numerics@gmail.com>
! 
! HydrostaticPrimitives_Driver.f90 is part of the Spectral Element Libraries in Fortran (SELF).
! 
! Permission is hereby granted, free of charge, to any person obtaining a copy of this software 
! and associated documentation files (the "Software"), to deal in the Software without restriction, 
! including without limitation the rights to use, copy, modify, merge, publish, distribute, 
! sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is 
! furnished to do so, subject to the following conditions: 
! 
! The above copyright notice and this permission notice shall be included in all copies or  
! substantial portions of the Software. 
! 
! THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING 
! BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND 
! NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, 
! DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, 
! OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE. 
!
! //////////////////////////////////////////////////////////////////////////////////////////////// !
 
 
PROGRAM HydrostaticPrimitives_Driver
! ========================================= Logs ================================================= !
!2016-05-11  Joseph Schoonover  schoonover.numerics@gmail.com 
!
! //////////////////////////////////////////////////////////////////////////////////////////////// ! 


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
