! Advection3D_Driver.f90
! 
! Copyright 2015 Joseph Schoonover <schoonover.numerics@gmail.com>
! 
! Advection3D_Driver.f90 is part of the Spectral Element Libraries in Fortran (SELF).
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
 
 
PROGRAM Advection3D_Driver
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
      CALL myAdv % GlobalTimeDerivative( ZERO, 1 )   

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
