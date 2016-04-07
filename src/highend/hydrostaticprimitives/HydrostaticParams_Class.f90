MODULE HydrostaticParams_Class
!
! HydrostaticParams_Class.f90
!
! 
! =========================================================================================== !
!



 USE ModelPrecision
 USE CommonRoutines
 USE ConstantsDictionary
 USE ModelFlags
 

 IMPLICIT NONE


    TYPE HydrostaticParams
       ! IterativeSolver
       INTEGER       :: MaximumIterates
       INTEGER       :: pcIterates
       REAL(prec)    :: tolerance
       REAL(prec)    :: pcTolerance
      ! TimeManagement
       REAL(prec)    :: dt
       INTEGER       :: iterInit
       INTEGER       :: nTimeSteps
       INTEGER       :: dumpFreq
       ! SpaceManagement
       CHARACTER(50) :: SpecMeshFile
       INTEGER       :: polyDeg
       INTEGER       :: geomPolyDeg
       INTEGER       :: nCutoff
       INTEGER       :: nXElem
       INTEGER       :: nYElem
       INTEGER       :: nZElem 
       INTEGER       :: nPlot
       REAL(prec)    :: dxPlot
       REAL(prec)    :: xScale
       REAL(prec)    :: yScale
       REAL(prec)    :: zScale
       ! PhysicalParameters
       REAL(prec)    :: g
       REAL(prec)    :: f0
      
       CONTAINS

       PROCEDURE :: Build => Build_HydrostaticParams

    END TYPE HydrostaticParams 
 

 CONTAINS


 SUBROUTINE Build_HydrostaticParams( thisParam )
 ! S/R Build
 !
 !
 ! ==============================================!
 ! DECLARATIONS
   CLASS( HydrostaticParams ), intent(out) :: thisParam
   ! LOCAL
   INTEGER :: nUnit
       ! IterativeSolver
       INTEGER       :: MaximumIterates
       INTEGER       :: pcIterates
       REAL(prec)    :: tolerance
       REAL(prec)    :: pcTolerance
       ! TimeManagement
       REAL(prec)    :: dt
       INTEGER       :: iterInit
       INTEGER       :: nTimeSteps
       INTEGER       :: dumpFreq
       ! SpaceManagement
       CHARACTER(50) :: SpecMeshFile
       INTEGER       :: geomPolyDeg
       INTEGER       :: polyDeg
       INTEGER       :: nCutoff
       INTEGER       :: nXElem
       INTEGER       :: nYElem
       INTEGER       :: nZElem
       INTEGER       :: nPlot
       REAL(prec)    :: xScale
       REAL(prec)    :: yScale
       REAL(prec)    :: zScale
       ! PhysicalParameters
       REAL(prec)    :: g
       REAL(prec)    :: f0
       
      NAMELIST / IterativeSolver / MaximumIterates, pcIterates, tolerance, pcTolerance
      NAMELIST / TimeManagement / dt, iterInit, nTimeSteps, dumpFreq
      NAMELIST / SpaceManagement / SpecMeshFile, geomPolyDeg, polyDeg, nCutoff, nXElem, nYElem, nZElem, &
                                   nPlot, xScale, yScale, zScale
      NAMELIST / PhysicalParameters / g, f0
      
      ! IterativeSolver
      MaximumIterates = 500   ! Max number of conjugate gradient iterations
      pcIterates      = 50     !
      tolerance       = 10.0_prec**(-8) ! conjugate gradient residual tolerance
      pcTolerance     = 10.0_prec**(-5)
      ! TimeManagement
      dt = ZERO
      iterInit = 0
      nTimeSteps = 0
      dumpFreq = 1
      ! SpaceManagement
      SpecMeshFile = nada
      geomPolyDeg = 5
      polyDeg = 5
      nCutoff = 3
      nXElem = 2
      nYElem = 2
      nZElem = 2
      nPlot = 10
      xScale = ONE
      yScale = ONE 
      zScale = ONE
      ! PhysicalParameters
      g = 9.806_prec
      f0 = ZERO
      
      ! Reading in the namelist FILE

      OPEN( UNIT = NEWUNIT(nUnit), FILE = 'hydrostatic.params')
         READ( UNIT = nUnit, NML = IterativeSolver )
         READ( UNIT = nUnit, NML = TimeManagement )
         READ( UNIT = nUnit, NML = SpaceManagement )
         READ( UNIT = nUnit, NML = PhysicalParameters )
      CLOSE( UNIT = nUnit ) 

      ! Sanity check - PRINT the results of the namelist READ to screen
      WRITE( UNIT = *, NML = IterativeSolver )
      WRITE( UNIT = *, NML = TimeManagement )
      WRITE( UNIT = *, NML = SpaceManagement )
      WRITE( UNIT = *, NML = PhysicalParameters )
      
      ! IterativeSolver
      thisParam % MaximumIterates = MaximumIterates
      thisParam % pcIterates      = pcIterates
      thisParam % tolerance       = tolerance
      thisParam % pcTolerance     = pcTolerance
      ! TimeManagement
      thisParam % dt         = dt
      thisParam % iterInit   = iterInit
      thisParam % nTimeSteps = nTimeSteps
      thisParam % dumpFreq   = dumpFreq
      ! SpaceManagement
      thisParam % SpecMeshFile = SpecMeshFile
      thisParam % geomPolyDeg  = geomPolyDeg
      thisParam % polyDeg      = polyDeg
      thisParam % nCutoff      = nCutoff
      thisParam % nXElem       = nXElem
      thisParam % nYElem       = nYElem 
      thisParam % nZElem       = nZElem
      thisParam % nPlot        = nPlot
      thisParam % dxPlot       = 2.0_prec/REAL(nPlot,prec)
      thisParam % xScale       = xScale
      thisParam % yScale       = yScale
      thisParam % zScale       = zScale
      ! Physical Parameters
      thisParam % g  = g
      thisParam % f0 = f0
      
      
 END SUBROUTINE Build_HydrostaticParams

END MODULE HydrostaticParams_Class
