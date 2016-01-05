MODULE RunParamsClass
!
! RunParamsClass.f90
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


    TYPE RunParams
       ! MODEL_FORM
       INTEGER       :: ModelFormulation
       INTEGER       :: GeoMetricsForm
       INTEGER       :: quadrature
       INTEGER       :: cgIterMax
       REAL(prec)    :: cgTol
       CHARACTER(40) :: ISMmeshFile
       INTEGER       :: nModes
       ! TIME_MANAGEMENT
       REAL(prec)    :: dt
       INTEGER       :: iterInit
       INTEGER       :: nTimeSteps
       INTEGER       :: dumpFreq
       ! SPACE_MANAGEMENT
       INTEGER       :: polyDeg
       INTEGER       :: geomPolyDeg
       INTEGER       :: nPlot
       REAL(prec)    :: dxPlot
       REAL(prec)    :: xScale
       REAL(prec)    :: yScale
       ! PHYSICAL_PARAMETERS
       REAL(prec)    :: g 
       REAL(prec)    :: f0
       REAL(prec)    :: betaX 
       REAL(prec)    :: betaY
       REAL(prec)    :: lDrag
       REAL(prec)    :: R
       REAL(prec)    :: N2
       REAL(prec)    :: rho0
      


       CONTAINS

       PROCEDURE :: Build => BuildParams

    END TYPE RunParams 
 

 CONTAINS


 SUBROUTINE BuildParams( thisParam )
 ! S/R BuildParams
 !
 !
 ! ==============================================!
 ! DECLARATIONS
   CLASS( RunParams ), intent(out) :: thisParam
   ! LOCAL
   INTEGER :: nUnit
   ! MODEL_FORM
       CHARACTER(20) :: ModelFormulation
       INTEGER       :: GeoMetricsForm
       CHARACTER(20) :: quadrature
       INTEGER       :: cgIterMax
       REAL(prec)    :: cgTol
       INTEGER       :: nModes
       ! TIME_MANAGEMENT
       REAL(prec)    :: dt
       INTEGER       :: iterInit
       INTEGER       :: nTimeSteps
       INTEGER       :: dumpFreq
       ! SPACE_MANAGEMENT
       INTEGER       :: polyDeg
       INTEGER       :: geomPolyDeg
       INTEGER       :: nPlot
       REAL(prec)    :: dxPlot
       REAL(prec)    :: xScale
       REAL(prec)    :: yScale
       ! PHYSICAL_PARAMETERS
       REAL(prec)    :: g 
       REAL(prec)    :: f0
       REAL(prec)    :: betaX 
       REAL(prec)    :: betaY
       REAL(prec)    :: lDrag
       REAL(prec)    :: R
       REAL(prec)    :: N2
       REAL(prec)    :: rho0


      ! Set the namelists
      NAMELIST / MODEL / ModelFormulation, GeoMetricsForm, quadrature, cgIterMax, cgTol, ISMmeshFile, nModes
      
      NAMELIST / TIME_MANAGEMENT / dt, iterInit, nTimeSteps, dumpFreq

      NAMELIST / SPACE_MANAGEMENT / polyDeg, polyDegOI, geomPolyDeg, nPlot, xScale, yScale
 
      NAMELIST / PHYSICAL_PARAMETERS / g, f0, betaX, betaY, lDrag, R, N2, rho0
 
      ! Set the default PARAMETERs 
      ! MODEL
      ModelFormulation = 'Linear'
      GeoMetricsForm = CLASSIC_CONTRAVARIANT
      quadrature = 'Gauss Lobatto'
      cgIterMax = 1000   ! Max number of conjugate gradient iterations
      cgTol = 10.0_prec**(-8) ! conjugate gradient residual tolerance

      ! TIME_MANAGEMENT
      dt = ZERO
      iterInit = 0
      nTimeSteps = 0
      dumpFreq = 1
      
      ! SPACE_MANAGEMENT (Default to isoparametric elements - if overintegration is USEd, default is to double the number of points )
      ! Plotting is defaulted to 10 evenly spaced points
      polyDeg = 5
      geomPolyDeg = 5
      nPlot = 10
      xScale = ONE
      yScale = ONE 

      ! PHYSICAL_PARAMETERS - midlatitude PARAMETERs on the planet earth with no beta and no bottom drag
      g = 9.81_prec        ! m/sec^2
      f0 = 10.0_prec**(-4) ! 1/sec
      betaX = ZERO         ! 1/(m*sec)
      betaY = ZERO         ! 1/(m*sec)
      lDrag = ZERO         ! 1/sec
      R = ZERO
      N2 = 10.0_prec**(-4) ! 1/sec^2
      rho0 = 999.98_prec   ! kg/m^3
 


     
      

      ! Reading in the namelist FILE

      OPEN( UNIT = NEWUNIT(nUnit), FILE = 'runtime.params')
         READ( UNIT = nUnit, NML = MODEL )
         READ( UNIT = nUnit, NML = TIME_MANAGEMENT )
         READ( UNIT = nUnit, NML = SPACE_MANAGEMENT )
         READ( UNIT = nUnit, NML = PHYSICAL_PARAMETERS )
      CLOSE( UNIT = nUnit ) 

      ! Sanity check - PRINT the results of the namelist READ to screen
      WRITE( UNIT = *, NML = MODEL )
      WRITE( UNIT = *, NML = TIME_MANAGEMENT )
      WRITE( UNIT = *, NML = SPACE_MANAGEMENT )
      WRITE( UNIT = *, NML = PHYSICAL_PARAMETERS )


      ! Fill in the data structure
      thisParam % ModelFormulation = GetFlagForChar( ModelFormulation )
      thisParam % GeoMetricsForm = GeoMetricsForm
      thisParam % quadrature = GetFlagForChar( quadrature )
      thisParam % cgIterMax = cgIterMax
      thisParam % cgTol = cgTol
      
      ! TIME_MANAGEMENT
      thisParam % dt = dt
      thisParam % iterInit = iterInit
      thisParam % nTimeSteps = nTimeSteps
      thisParam % dumpFreq = dumpFreq
      
      ! SPACE_MANAGEMENT (Default to isoparametric elements - if overintegration is USEd, default is to double the number of points )
      ! Plotting is defaulted to 10 evenly spaced points
      thisParam % polyDeg = polyDeg
      thisParam % nPlot = nPlot
      thisParam % dxPlot = 2.0_prec/REAL(nPlot,prec)
      thisParam % xScale = xScale
      thisParam % yScale = yScale

      ! PHYSICAL_PARAMETERS - midlatitude PARAMETERs on the planet earth with no beta and no bottom drag
      thisParam % g = g          ! m/sec^2
      thisParam % f0 = f0        ! 1/sec
      thisParam % betaX = betaX  ! 1/(m*sec)
      thisParam % betaY = betaY  ! 1/(m*sec)
      thisParam % lDrag = lDrag  ! 1/sec
      thisParam % R = R
      thisParam % N2 = N2        ! 1/sec^2
      thisParam % rho0 = rho0

 END SUBROUTINE BuildParams

END MODULE RunParamsClass
