MODULE AdvectionParamsClass
!
! AdvectionParamsClass.f90
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


    TYPE AdvectionParams
      ! TracerSetup
       INTEGER       :: nTracers
       REAL(prec)    :: kappa
      ! TimeManagement
       REAL(prec)    :: dt
       INTEGER       :: iterInit
       INTEGER       :: nTimeSteps
       INTEGER       :: dumpFreq
       ! SpaceManagement
       CHARACTER(50) :: SpecMeshFile
       INTEGER       :: polyDeg
       INTEGER       :: nXElem
       INTEGER       :: nYElem
       INTEGER       :: geomPolyDeg
       INTEGER       :: nPlot
       REAL(prec)    :: dxPlot
       REAL(prec)    :: xScale
       REAL(prec)    :: yScale
       ! VelocityField
    !   REAL(prec)    :: 
      
       CONTAINS

       PROCEDURE :: Build => Build_AdvectionParams

    END TYPE AdvectionParams 
 

 CONTAINS


 SUBROUTINE Build_AdvectionParams( thisParam )
 ! S/R Build
 !
 !
 ! ==============================================!
 ! DECLARATIONS
   CLASS( AdvectionParams ), intent(out) :: thisParam
   ! LOCAL
   INTEGER :: nUnit
       ! TracerSetup
       INTEGER       :: nTracers
       REAL(prec)    :: kappa
       ! TimeManagement
       REAL(prec)    :: dt
       INTEGER       :: iterInit
       INTEGER       :: nTimeSteps
       INTEGER       :: dumpFreq
       ! SpaceManagement
       CHARACTER(50) :: SpecMeshFile
       INTEGER       :: polyDeg
       INTEGER       :: nXElem
       INTEGER       :: nYElem
       INTEGER       :: nPlot
       REAL(prec)    :: xScale
       REAL(prec)    :: yScale

       
      NAMELIST / TracerSetup / nTracers, kappa
      NAMELIST / TimeManagement / dt, iterInit, nTimeSteps, dumpFreq
      NAMELIST / SpaceManagement / SpecMeshFile, polyDeg, nXElem, nYElem, nPlot, xScale, yScale
      
      ! TracerSetup
      nTracers = 2
      kappa    = ZERO
      ! TimeManagement
      dt = ZERO
      iterInit = 0
      nTimeSteps = 0
      dumpFreq = 1
      ! SpaceManagement
      SpecMeshFile = nada
      polyDeg = 5
      nXElem = 5
      nYElem = 5
      nPlot = 10
      xScale = ONE
      yScale = ONE 
      
      ! Reading in the namelist FILE

      OPEN( UNIT = NEWUNIT(nUnit), FILE = 'advection.params')
         READ( UNIT = nUnit, NML = TracerSetup )
         READ( UNIT = nUnit, NML = TimeManagement )
         READ( UNIT = nUnit, NML = SpaceManagement )
      CLOSE( UNIT = nUnit ) 

      ! Sanity check - PRINT the results of the namelist READ to screen
      WRITE( UNIT = *, NML = TracerSetup )
      WRITE( UNIT = *, NML = TimeManagement )
      WRITE( UNIT = *, NML = SpaceManagement )
      
      ! MODEL_FORM
      thisParam % nTracers = nTracers
      thisParam % kappa    = kappa
      ! TIME_MANAGEMENT
      thisParam % dt = dt
      thisParam % iterInit = iterInit
      thisParam % nTimeSteps = nTimeSteps
      thisParam % dumpFreq = dumpFreq
      ! SPACE_MANAGEMENT (Default to isoparametric elements - if overintegration is USEd, default is to double the number of points )
      ! Plotting is defaulted to 10 evenly spaced points
      thisParam % SpecMeshFile = SpecMeshFile
      thisParam % polyDeg = polyDeg
      thisParam % nXElem = nXElem
      thisParam % nYElem = nYElem 
      thisParam % nPlot = nPlot
      thisParam % dxPlot = 2.0_prec/REAL(nPlot,prec)
      thisParam % xScale = xScale
      thisParam % yScale = yScale

      
      
 END SUBROUTINE Build_AdvectionParams

END MODULE AdvectionParamsClass
