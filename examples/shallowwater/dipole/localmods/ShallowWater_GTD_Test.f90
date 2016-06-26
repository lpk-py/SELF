! ShallowWater_GTD_Test.f90
! 
! Copyright 2015 Joe <joe@clay>
! 
! ShallowWater_GTD_Test.f90 is part of the Spectral Element Libraries in Fortran (SELF).
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
 
 
PROGRAM ShallowWater_GTD_Test


! src/common/
USE ModelPrecision
!USE Timing
! src/highend/shallowwater
USE ConservativeShallowWaterClass

 IMPLICIT NONE

 TYPE( ShallowWater ) :: mysw
 TYPE( MultiTimers )  :: timers
 INTEGER :: i, n

 
   CALL mysw % Build( )
   n = mysw % params % nTimesteps

   !CALL timers % Build( )

   !CALL timers % AddTimer( 'GlobalTimeDerivative', 1 )

   !!!$OMP PARALLEL
   !CALL timers % StartThisTimer( 1 ) 
   DO i = 1, n
      CALL mysw % GlobalTimeDerivative( ZERO, 0 )
   ENDDO
   !CALL timers % StopThisTimer( 1 )
   !!!$OMP END PARALLEL 
 
   !CALL timers % Write_MultiTimers( )

   !CALL timers % Trash( )
   CALL mysw % Trash( )


END PROGRAM ShallowWater_GTD_Test
