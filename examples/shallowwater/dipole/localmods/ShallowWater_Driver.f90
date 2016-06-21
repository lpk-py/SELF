! ShallowWater_Driver.f90
! 
! Copyright 2015 Joe <joe@clay>
! 
! ShallowWater_Driver.f90 is part of the Spectral Element Libraries in Fortran (SELF).
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
 
 
PROGRAM ShallowWater_Driver



! src/common/
USE ModelPrecision
! src/geom/
USE QuadMeshClass
! src/highend/shallowwater
USE SWParamsClass
USE ConservativeShallowWaterClass

 IMPLICIT NONE

 TYPE( ShallowWater ) :: mysw

 INTEGER :: iter0, nT, dFreq
 INTEGER :: iT, fUnit, nDumps, k
 REAL(prec) :: tn, deltaT
 REAL(prec), ALLOCATABLE :: KE(:), PE(:)
 
 CHARACTER(len=10) :: iterChar
 
    CALL mysw % Build( )
 
    ! ////////////////////////////////// Time Stepping /////////////////////////////////////////// !
       
      iter0 = mysw % params % iterInit
      nT = mysw % params % nTimeSteps
      dFreq = mysw % params % dumpFreq
      deltaT = mysw % params % dt

      nDumps = (nT)/dFreq
      ALLOCATE( KE(0:nDumps), PE(0:nDumps) )
      KE = ZERO
      PE = ZERO
      
      WRITE( iterChar, '(I10.10)') iter0   
      CALL mysw % GlobalTimeDerivative( ZERO, 0 ) 
      CALL mysw % WriteTecplot( 'ShallowWater.'//iterChar )
      k = 0
      CALL mysw % IntegrateEnergies( KE(k), PE(k) ) 
      
      DO iT = iter0, iter0+nT-1 ! Loop over time-steps

          tn = real( iT, prec )*deltaT ! time at the beginning of this timestep          
          CALL mysw % ForwardStepRK3( tn ) ! Forward Step

          IF( mod( iT+1, dFreq) == 0 )THEN
             k = k+1
             WRITE(iterChar,'(I10.10)') iT+1
             CALL mysw % WriteTecplot( 'ShallowWater.'//iterChar )

 !            CALL mysw % IntegrateEnergies( KE(k), PE(k) )

          ENDIF

       ENDDO
       
       ! ///////////// Energy File I/O ////////////// !

!       WRITE( iterChar, '(I10.10)') iter0 
!       OPEN( UNIT = NewUnit(fUnit), &
!            FILE = 'Energies.'//iterChar//'.curve', &
!            FORM = 'FORMATTED' )
            
!       WRITE( fUnit, * )'#KE'
!       k = 0
!       WRITE( fUnit, * ) real( iter0, prec )*deltaT, KE(0)
!       DO iT = iter0, iter0+nT-1 ! Loop over time-steps
!          IF( mod( iT+1, dFreq) == 0 )THEN
!             k = k+1
!             tn = real( iT+1, prec )*deltaT
!             WRITE( fUnit, * ) tn, KE(k)
!          ENDIF
!       ENDDO
!       WRITE( fUnit, * )'#PE'
!       k = 0
!       WRITE( fUnit, * ) real( iter0, prec )*deltaT, PE(0)
!       DO iT = iter0, iter0+nT-1 ! Loop over time-steps
!          IF( mod( iT+1, dFreq) == 0 )THEN
!             k = k+1
!             tn = real( iT+1, prec )*deltaT
!             WRITE( fUnit, * ) tn, PE(k)
!          ENDIF
!       ENDDO
!       CLOSE( fUnit )
       
    CALL mysw % WritePickup(iter0+nT) 
    CALL mysw % mesh % WriteTecplot( )
    
    CALL mysw % Trash( )


END PROGRAM ShallowWater_Driver
