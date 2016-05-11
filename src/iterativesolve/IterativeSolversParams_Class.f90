! IterativeSolversParams_Class.f90
! 
! Copyright 2015 Joseph Schoonover <schoonover.numerics@gmail.com>
! 
! IterativeSolversParams_Class.f90 is part of the Spectral Element Libraries in Fortran (SELF).
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
 
 
MODULE IterativeSolversParams_Class
! ========================================= Logs ================================================= !
!2016-05-11  Joseph Schoonover  schoonover.numerics@gmail.com 
!
! //////////////////////////////////////////////////////////////////////////////////////////////// ! 

 USE ModelPrecision
 USE CommonRoutines
 USE ConstantsDictionary
 USE ModelFlags
 

 IMPLICIT NONE


    TYPE IterativeSolversParams
       ! MODEL_FORM
       INTEGER       :: MaximumIterates
       REAL(prec)    :: tolerance
       INTEGER       :: nDOF


       CONTAINS

       PROCEDURE :: Build => BuildParams

    END TYPE IterativeSolversParams 
 

 CONTAINS


 SUBROUTINE BuildParams( thisParam )
 ! S/R BuildParams
 !
 !
 ! ==============================================!
 ! DECLARATIONS
   CLASS( IterativeSolversParams ), intent(out) :: thisParam
   ! LOCAL
   INTEGER :: nUnit
   ! MODEL_FORM
   INTEGER    :: MaximumIterates
   REAL(prec) :: tolerance
   INTEGER    :: nDOF


      ! Set the namelists
      NAMELIST / SolverCriteria / MaximumIterates, tolerance, nDOF
      
 
      ! Set the default parameters
      ! SolverCriteria
      MaximumIterates = 1000   ! Max number of conjugate gradient iterations
      tolerance = 10.0_prec**(-8) ! conjugate gradient residual tolerance
      nDOF = 1

      ! Reading in the namelist FILE

      OPEN( UNIT = NEWUNIT(nUnit), FILE = 'itsolve.params')
         READ( UNIT = nUnit, NML = SolverCriteria )
      CLOSE( UNIT = nUnit ) 

      ! Sanity check - PRINT the results of the namelist READ to screen
      WRITE( UNIT = *, NML = SolverCriteria )

      ! Fill in the data structure
      thisParam % MaximumIterates = MaximumIterates
      thisParam % tolerance       = tolerance
      thisParam % nDOF            = nDOF
      
      


 END SUBROUTINE BuildParams

END MODULE IterativeSolversParams_Class
