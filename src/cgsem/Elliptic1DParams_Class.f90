!  Elliptic1DParams_Class.f90
!  
!  Copyright 2015 Joseph Schoonover <schoonover.numerics@gmail.com>
!  
!  Elliptic1DParams_Class.f90 is part of the Spectral Element Libraries in Fortran (SELF).
!
!  Permission is hereby granted, free of charge, to any person obtaining a copy of this software
!  and associated documentation files (the "Software"), to deal in the Software without restriction,
!  including without limitation the rights to use, copy, modify, merge, publish, distribute, 
!  sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is 
!  furnished to do so, subject to the following conditions:
!
!  The above copyright notice and this permission notice shall be included in all copies or 
!  substantial portions of the Software.
!
!  THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING
!  BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
!  NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM,
!  DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, 
!  OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
!  
! //////////////////////////////////////////////////////////////////////////////////////////////// !
 
MODULE Elliptic1DParams_Class
! ========================================= Logs ================================================= !
! yyyy-mm-dd  Joe  <joe@clay>
! 
! //////////////////////////////////////////////////////////////////////////////////////////////// !


 USE ModelPrecision
 USE CommonRoutines
 USE ConstantsDictionary
 USE ModelFlags
 

 IMPLICIT NONE


    TYPE Elliptic1DParams
       ! SolverCriteria
       INTEGER       :: MaximumIterates
       REAL(prec)    :: tolerance

       ! SpaceManagement
       INTEGER       :: polyDeg
       INTEGER       :: nElems
       INTEGER       :: nPlot
       REAL(prec)    :: dxPlot
       REAL(prec)    :: xScale
      


       CONTAINS

       PROCEDURE :: Build => BuildParams_Elliptic1D

    END TYPE Elliptic1DParams 
 

 CONTAINS


 SUBROUTINE BuildParams_Elliptic1D( thisParam )
 ! S/R BuildParams
 !
 !
 ! ==============================================!
 ! DECLARATIONS
   CLASS( Elliptic1DParams ), intent(out) :: thisParam
   ! LOCAL
   INTEGER :: nUnit
   
       ! SolverCriteria
       INTEGER       :: MaximumIterates
       REAL(prec)    :: tolerance

       ! SpaceManagement
       INTEGER       :: polyDeg
       INTEGER       :: nElems
       INTEGER       :: nPlot
!       REAL(prec)    :: dxPlot
       REAL(prec)    :: xScale

      ! Set the namelists
      NAMELIST / SolverCriteria / MaximumIterates, tolerance

      NAMELIST / SpaceManagement / polyDeg, nElems, nPlot, xScale
 
      ! Set the default parameters
      ! SolverCriteria
      MaximumIterates = 1000   ! Max number of conjugate gradient iterations
      tolerance = 10.0_prec**(-8) ! conjugate gradient residual tolerance

      ! SpaceManagement 
      polyDeg = 5
      nElems  = 5
      nPlot = 10
      xScale = ONE

      ! Reading in the namelist FILE

      OPEN( UNIT = NEWUNIT(nUnit), FILE = 'runtime.params')
         READ( UNIT = nUnit, NML = SolverCriteria )
         READ( UNIT = nUnit, NML = SpaceManagement )
      CLOSE( UNIT = nUnit ) 

      ! Sanity check - PRINT the results of the namelist READ to screen
      WRITE( UNIT = *, NML = SolverCriteria )
      WRITE( UNIT = *, NML = SpaceManagement )

      ! Fill in the data structure
      ! SolverCriteria
      thisParam % MaximumIterates = MaximumIterates
      thisParam % tolerance       = tolerance
      
      ! SpaceManagement
      thisParam % polyDeg = polyDeg
      thisParam % nElems  = nElems
      thisParam % nPlot   = nPlot
      thisParam % dxPlot  = 2.0_prec/REAL(nPlot,prec)
      thisParam % xScale  = xScale


 END SUBROUTINE BuildParams_Elliptic1D

END MODULE Elliptic1DParams_Class
