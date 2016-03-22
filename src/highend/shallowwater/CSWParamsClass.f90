MODULE CSWParamsClass
!
! CSWParamsClass.f90
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


    TYPE SWParams
      ! MODEL_FORM
       INTEGER       :: ModelFormulation
      ! TIME_MANAGEMENT
       REAL(prec)    :: dt
       INTEGER       :: iterInit
       INTEGER       :: nTimeSteps
       INTEGER       :: dumpFreq
       ! SPACE_MANAGEMENT
       CHARACTER(50) :: SpecMeshFile
       INTEGER       :: polyDeg
       INTEGER       :: nXElem
       INTEGER       :: nYElem
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
       REAL(prec)    :: linearDrag
       ! JET_PARAMETERS
       REAL(prec)    :: dEta
       REAL(prec)    :: jetCenter
       REAL(prec)    :: jetWidth
       ! Bathymetry
       REAL(prec)    :: shelfWidth
       REAL(prec)    :: dLx
       REAL(prec)    :: dLy
       REAL(prec)    :: yc
       REAL(prec)    :: h0
       REAL(prec)    :: h1
      
       CONTAINS

       PROCEDURE :: Build => Build_SWParams

    END TYPE SWParams 
 

 CONTAINS


 SUBROUTINE Build_SWParams( thisParam )
 ! S/R Build
 !
 !
 ! ==============================================!
 ! DECLARATIONS
   CLASS( SWParams ), intent(out) :: thisParam
   ! LOCAL
   INTEGER :: nUnit
       ! MODEL_FORM
       CHARACTER(50) :: ModelFormulation
       ! TIME_MANAGEMENT
       REAL(prec)    :: dt
       INTEGER       :: iterInit
       INTEGER       :: nTimeSteps
       INTEGER       :: dumpFreq
       ! SPACE_MANAGEMENT
       CHARACTER(50) :: SpecMeshFile
       INTEGER       :: polyDeg
       INTEGER       :: nXElem
       INTEGER       :: nYElem
       INTEGER       :: nPlot
       REAL(prec)    :: dxPlot
       REAL(prec)    :: xScale
       REAL(prec)    :: yScale
       ! PHYSICAL_PARAMETERS
       REAL(prec)    :: g 
       REAL(prec)    :: f0
       REAL(prec)    :: betaX 
       REAL(prec)    :: betaY
       REAL(prec)    :: linearDrag
       ! JET_PARAMETERS
       REAL(prec)    :: dEta
       REAL(prec)    :: jetCenter
       REAL(prec)    :: jetWidth
       ! Bathymetry
       REAL(prec)    :: shelfWidth
       REAL(prec)    :: dLx
       REAL(prec)    :: dLy
       REAL(prec)    :: yc
       REAL(prec)    :: h0
       REAL(prec)    :: h1
       
      NAMELIST / MODEL_FORM / ModelFormulation
      NAMELIST / TIME_MANAGEMENT / dt, iterInit, nTimeSteps, dumpFreq
      NAMELIST / SPACE_MANAGEMENT / SpecMeshFile, polyDeg, nXElem, nYElem, nPlot, xScale, yScale
      NAMELIST / PHYSICAL_PARAMETERS / g, f0, betaX, betaY, linearDrag
      NAMELIST / JET_PARAMETERS / dEta, jetCenter, jetWidth
      NAMELIST / BATHYMETRY / shelfWidth, dLx, dLy, yc, h0, h1
      
      ! MODEL_FORM
      ModelFormulation = 'LINEAR'
      ! TIME_MANAGEMENT
      dt = ZERO
      iterInit = 0
      nTimeSteps = 0
      dumpFreq = 1
      ! SPACE_MANAGEMENT
      SpecMeshFile = nada
      polyDeg = 5
      nXElem = 5
      nYElem = 5
      nPlot = 10
      xScale = ONE
      yScale = ONE 
      ! PHYSICAL_PARAMETERS
      g = ONE        ! m/sec^2
      f0 = 10.0_prec ! 1/sec
      betaX = ZERO         ! 1/(m*sec)
      betaY = ZERO         ! 1/(m*sec)
      linearDrag = ZERO         ! 1/sec
       ! JET_PARAMETERS
      dEta = HALF
      jetCenter = 150000.0_prec 
      jetWidth  = 100000.0_prec
       ! Bathymetry
      shelfWidth = 300000.0_prec
      dLx = 20000.0_prec
      dLy = 1000000.0_prec
      yc  = 1000000.0_prec
      h0 = 500.0_prec
      h1 = 3000.0_prec
      ! Reading in the namelist FILE

      OPEN( UNIT = NEWUNIT(nUnit), FILE = 'cshallowwater.params')
         READ( UNIT = nUnit, NML = MODEL_FORM )
         READ( UNIT = nUnit, NML = TIME_MANAGEMENT )
         READ( UNIT = nUnit, NML = SPACE_MANAGEMENT )
         READ( UNIT = nUnit, NML = PHYSICAL_PARAMETERS )
         READ( UNIT = nUnit, NML = JET_PARAMETERS )
         READ( UNIT = nUnit, NML = BATHYMETRY )
      CLOSE( UNIT = nUnit ) 

      ! Sanity check - PRINT the results of the namelist READ to screen
      WRITE( UNIT = *, NML = MODEL_FORM )
      WRITE( UNIT = *, NML = TIME_MANAGEMENT )
      WRITE( UNIT = *, NML = SPACE_MANAGEMENT )
      WRITE( UNIT = *, NML = PHYSICAL_PARAMETERS )
      WRITE( UNIT = *, NML = JET_PARAMETERS )
      WRITE( UNIT = *, NML = BATHYMETRY )
      
      ! MODEL_FORM
      thisParam % ModelFormulation = GetFlagForChar( ModelFormulation )
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
      ! PHYSICAL_PARAMETERS - midlatitude PARAMETERs on the planet earth with no beta and no bottom drag
      thisParam % g = g          ! m/sec^2
      thisParam % f0 = f0        ! 1/sec
      thisParam % betaX = betaX  ! 1/(m*sec)
      thisParam % betaY = betaY  ! 1/(m*sec)
      thisParam % linearDrag = linearDrag  ! 1/sec
       ! JET_PARAMETERS
      thisParam % dEta = dEta
      thisParam % jetCenter = jetCenter 
      thisParam % jetWidth = jetWidth 
       ! Bathymetry
      thisParam % shelfWidth = shelfWidth
      thisParam % dLx = dLx
      thisParam % dLy = dLy
      thisParam % yc  = yc
      thisParam % h0 = h0
      thisParam % h1 = h1
      
      
 END SUBROUTINE Build_SWParams

END MODULE CSWParamsClass
