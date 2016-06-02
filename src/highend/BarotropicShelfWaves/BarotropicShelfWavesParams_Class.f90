MODULE BarotropicShelfWavesParams_Class
!
! BarotropicShelfWavesParamsClass.f90
!
! Joe Schoonover
!
! schoonover.numerics@gmail.com
!
!
! This MODULE provides a data structure with a build routine for reading in a namelist FILE
! that is USEd for setting run-time PARAMETERs for the SEM software.
! 
! =========================================================================================== !
!



 USE ModelPrecision
 USE CommonRoutines
 USE ConstantsDictionary
 USE ModelFlags
 

 IMPLICIT NONE

   INTEGER, PARAMETER :: nMaxModes = 5

    TYPE BarotropicShelfWavesParams
       ! SolverCriteria
       INTEGER       :: MaximumIterates
       INTEGER       :: nEpairs
       INTEGER       :: nMaxToSolve
       REAL(prec)    :: tolerance
       ! TimeManagement
       REAL(prec)    :: dt
       INTEGER       :: iterInit
       INTEGER       :: nTimeSteps
       INTEGER       :: dumpFreq
       ! SpaceManagement
       INTEGER       :: polyDeg
       INTEGER       :: nElems
       INTEGER       :: nPlot
       REAL(prec)    :: dxPlot
       REAL(prec)    :: xScale
       REAL(prec)    :: dy
       INTEGER       :: nYElems
       ! Physical
       REAL(prec)    :: f0
       REAL(prec)    :: Linner
       REAL(prec)    :: Lslope
       REAL(prec)    :: xSlope
       REAL(prec)    :: hcoast
       REAL(prec)    :: hShelf
       REAL(prec)    :: hDeep
       ! ShelfPerturbation
       REAL(prec)    :: vbar
       REAL(prec)    :: dLinner
       REAL(prec)    :: dXslope
       REAL(prec)    :: steepeningZoneLength
       REAL(prec)    :: steepeningCenter
       ! FileIO
       INTEGER :: nModes2Write
       INTEGER :: modes2Write(1:nMaxModes)

       CONTAINS

       PROCEDURE :: Build => BuildParams_BarotropicShelfWaves

    END TYPE BarotropicShelfWavesParams 
 

 CONTAINS


 SUBROUTINE BuildParams_BarotropicShelfWaves( thisParam )
 ! S/R BuildParams
 !
 !
 ! ==============================================!
 ! DECLARATIONS
   CLASS( BarotropicShelfWavesParams ), intent(out) :: thisParam
   ! LOCAL
   INTEGER :: nUnit
   
       ! SolverCriteria
       INTEGER       :: MaximumIterates
       INTEGER       :: nEpairs
       INTEGER       :: nMaxToSolve
       REAL(prec)    :: tolerance
       ! TimeManagement
       REAL(prec)    :: dt
       INTEGER       :: iterInit
       INTEGER       :: nTimeSteps
       INTEGER       :: dumpFreq
       ! SpaceManagement
       INTEGER       :: polyDeg
       INTEGER       :: nElems
       INTEGER       :: nPlot
       REAL(prec)    :: dxPlot
       REAL(prec)    :: xScale
       REAL(prec)    :: dy
       INTEGER       :: nYElems
       ! Physical
       REAL(prec)    :: f0
       REAL(prec)    :: Linner
       REAL(prec)    :: Lslope
       REAL(prec)    :: xSlope
       REAL(prec)    :: hcoast
       REAL(prec)    :: hShelf
       REAL(prec)    :: hDeep
       ! ShelfPerturbation
       REAL(prec)    :: vbar
       REAL(prec)    :: dLinner
       REAL(prec)    :: dXslope
       REAL(prec)    :: steepeningZoneLength
       REAL(prec)    :: steepeningCenter
       ! FileIO
       INTEGER :: nModes2Write
       INTEGER :: modes2Write(1:nMaxModes)
       
       !LOCAL
      INTEGER :: i

      ! Set the namelists
      NAMELIST / SolverCriteria / MaximumIterates, nEpairs, nMaxToSolve, tolerance
      NAMELIST / TimeManagement / dt, iterInit, nTimeSteps, dumpFreq
      NAMELIST / SpaceManagement / polyDeg, nElems, nPlot, xScale, dy, nYElems
      NAMELIST / Physical / f0, Linner, Lslope, xSlope, hcoast, &
                            hShelf, hDeep
      NAMELIST / ShelfPerturbation / vbar, dLinner, dXslope, steepeningZoneLength, steepeningCenter
      NAMELIST / FileIO / nModes2Write, modes2Write
 
      ! Set the default parameters
      ! SolverCriteria
      MaximumIterates = 1000   ! Max number of conjugate gradient iterations
      nEpairs         = 10
      nMaxToSolve     = 15
      tolerance       = 10.0_prec**(-8) ! conjugate gradient residual tolerance

      ! TimeManagement
      dt = ZERO
      iterInit = 0
      nTimeSteps = 0
      dumpFreq = 1
      ! SpaceManagement 
      polyDeg = 10
      nElems  = 10
      nPlot   = 10
      xScale  = 2.0D6
      dy      = 5.0_prec*10.0_prec**(3)
      nYElems = 200

      ! Physical
      f0     = 1.0D-4
      Linner = 1.0D5
      Lslope = 5.0D4
      xSlope = Linner
      hcoast = 50.0_prec
      hShelf = 2000.0_prec
      hDeep  = 4000.0_prec

      ! ShelfPerturbation
      vbar    = ONE
      dLinner = 5.0D4
      dXslope = 5.0D4
      steepeningZoneLength = 200.0D3
      steepeningCenter     = 500.0D3

      ! FileIO
      nModes2Write   = nMaxModes
      DO i = 1, nMaxModes 
         modes2Write(i) = i
      ENDDO

      ! Reading in the namelist FILE

      OPEN( UNIT = NEWUNIT(nUnit), FILE = 'runtime.params')
         READ( UNIT = nUnit, NML = SolverCriteria )
         READ( UNIT = nUnit, NML = TimeManagement )
         READ( UNIT = nUnit, NML = SpaceManagement )
         READ( UNIT = nUnit, NML = Physical )
         READ( UNIT = nUnit, NML = ShelfPerturbation )
         READ( UNIT = nUnit, NML = FileIO )
      CLOSE( UNIT = nUnit ) 

      ! Sanity check - PRINT the results of the namelist READ to screen
      WRITE( UNIT = *, NML = SolverCriteria )
      WRITE( UNIT = *, NML = TimeManagement )
      WRITE( UNIT = *, NML = SpaceManagement )
      WRITE( UNIT = *, NML = Physical )
      WRITE( UNIT = *, NML = ShelfPerturbation )
      WRITE( UNIT = *, NML = FileIO )

      ! Fill in the data structure
      ! SolverCriteria
      thisParam % MaximumIterates = MaximumIterates
      thisParam % nEpairs         = nEpairs
      thisParam % nMaxToSolve     = nMaxToSolve
      thisParam % tolerance       = tolerance
      
      ! TimeManagement
      thisParam % dt = dt
      thisParam % iterInit = iterInit
      thisParam % nTimeSteps = nTimeSteps
      thisParam % dumpFreq = dumpFreq

      ! SpaceManagement
      thisParam % polyDeg = polyDeg
      thisParam % nElems  = nElems
      thisParam % nPlot   = nPlot
      thisParam % dxPlot  = 2.0_prec/REAL(nPlot,prec)
      thisParam % xScale  = xScale
      thisParam % dy      = dy
      thisParam % nYElems = nYElems

      ! Physical
      thisParam % f0 = f0
      thisParam % Linner = Linner
      thisParam % Lslope = Lslope
      thisParam % xSlope = xSlope
      thisParam % hcoast = hcoast
      thisParam % hShelf = hShelf
      thisParam % hDeep  = hDeep

      ! ShelfPerturbation
      thisParam % vbar    = vbar
      thisParam % dLinner = dLinner
      thisParam % dXSlope = dXSlope
      thisParam % steepeningZoneLength = steepeningZoneLength
      thisParam % steepeningCenter = steepeningCenter

      ! FileIO
      thisParam % nModes2Write = nModes2Write
      thisParam % modes2Write  = modes2Write

 END SUBROUTINE BuildParams_BarotropicShelfWaves

END MODULE BarotropicShelfWavesParams_Class
