! StommelParams_Class.f90
! 
! Copyright 2015 Joseph Schoonover <schoonover.numerics@gmail.com>
! 
! StommelParams_Class.f90 is part of the Spectral Element Libraries in Fortran (SELF).
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
 
 
MODULE StommelParams_Class
! ========================================= Logs ================================================= !
!2016-05-13  Joseph Schoonover  schoonover.numerics@gmail.com 
!
! //////////////////////////////////////////////////////////////////////////////////////////////// ! 

 USE ModelPrecision
 USE CommonRoutines
 USE ConstantsDictionary
 USE ModelFlags
 

 IMPLICIT NONE


    TYPE StommelParams
       ! SolverCriteria
       INTEGER       :: MaximumIterates
       INTEGER       :: mInnerIters
       INTEGER       :: pcIterates
       REAL(prec)    :: tolerance
       REAL(prec)    :: pcTolerance
       INTEGER       :: pcDegree

       ! SpaceManagement
       CHARACTER(50) :: SpecMeshFile
       INTEGER       :: polyDeg
       INTEGER       :: nXElem
       INTEGER       :: nYElem
       INTEGER       :: nPlot
       REAL(prec)    :: dxPlot
       REAL(prec)    :: xScale
       REAL(prec)    :: yScale

       ! StommelParameters
       REAL(prec)    :: cDrag
       REAL(prec)    :: f0
       REAL(prec)    :: betaX
       REAL(prec)    :: betaY 
      


       CONTAINS

       PROCEDURE :: Build => BuildParams_Stommel

    END TYPE StommelParams 
 

 CONTAINS


 SUBROUTINE BuildParams_Stommel( thisParam )
 ! S/R BuildParams
 !
 !
 ! ==============================================!
 ! DECLARATIONS
   CLASS( StommelParams ), intent(out) :: thisParam
   ! LOCAL
   INTEGER :: nUnit
   
       ! SolverCriteria
       INTEGER       :: MaximumIterates
       INTEGER       :: mInnerIters
       INTEGER       :: pcIterates
       REAL(prec)    :: tolerance
       REAL(prec)    :: pcTolerance
       INTEGER       :: pcDegree

       ! SpaceManagement
       CHARACTER(50) :: SpecMeshFile
       INTEGER       :: polyDeg
       INTEGER       :: nXElem
       INTEGER       :: nYElem
       INTEGER       :: nPlot
!       REAL(prec)    :: dxPlot
       REAL(prec)    :: xScale
       REAL(prec)    :: yScale

       ! StommelParameters
       REAL(prec)    :: cDrag
       REAL(prec)    :: f0
       REAL(prec)    :: betaX
       REAL(prec)    :: betaY

      ! Set the namelists
      NAMELIST / SolverCriteria / MaximumIterates, mInnerIters, pcIterates, tolerance, pcTolerance, pcDegree

      NAMELIST / SpaceManagement / SpecMeshFile, polyDeg, nXElem, nYElem, nPlot, xScale, yScale
 
      NAMELIST / StommelParameters / cDrag, f0, betaX, betaY

      ! Set the default parameters
      ! SolverCriteria
      MaximumIterates = 500   ! Max number of conjugate gradient iterations
      mInnerIters     = 50
      pcIterates      = 50     !
      tolerance       = 10.0_prec**(-8) ! conjugate gradient residual tolerance
      pcTolerance     = 10.0_prec**(-4) 
      pcDegree        = 1
      ! SpaceManagement
      SpecMeshFile = nada 
      polyDeg      = 5
      nXElem       = 5
      nYElem       = 5
      nPlot        = 10
      xScale       = 1.0_prec*10.0_prec**(6)
      yScale       = 3.0_prec*10.0_prec**(6)
      ! StommelParameters
      cDrag = 0.01_prec
      f0    = 1.0_prec*10.0_prec**(-4)
      betaX = ZERO
      betaY = 1.0_prec*10.0_prec**(-11)

      ! Reading in the namelist FILE

      OPEN( UNIT = NEWUNIT(nUnit), FILE = 'runtime.params')
         READ( UNIT = nUnit, NML = SolverCriteria )
         READ( UNIT = nUnit, NML = SpaceManagement )
         READ( Unit = nUnit, NML = StommelParameters )
      CLOSE( UNIT = nUnit ) 

      ! Sanity check - PRINT the results of the namelist READ to screen
      WRITE( UNIT = *, NML = SolverCriteria )
      WRITE( UNIT = *, NML = SpaceManagement )
      WRITE( UNIT = *, NML = StommelParameters )

      ! Fill in the data structure
      ! SolverCriteria
      thisParam % MaximumIterates = MaximumIterates
      thisParam % mInnerIters     = mInnerIters
      thisParam % pcIterates      = pcIterates
      thisParam % tolerance       = tolerance
      thisParam % pcTolerance     = pcTolerance
      thisParam % pcDegree        = 1
      
      ! SpaceManagement
      thisParam % SpecMeshFile = SpecMeshFile
      thisParam % polyDeg = polyDeg
      thisParam % nXElem  = nXElem
      thisParam % nYElem  = nYElem
      thisParam % nPlot   = nPlot
      thisParam % dxPlot  = 2.0_prec/REAL(nPlot,prec)
      thisParam % xScale  = xScale
      thisParam % yScale  = yScale

      ! StommelParameters
      thisParam % cDrag = cDrag
      thisParam % f0    = f0
      thisParam % betaX = betaX
      thisParam % betaY = betaY


 END SUBROUTINE BuildParams_Stommel

END MODULE StommelParams_Class
