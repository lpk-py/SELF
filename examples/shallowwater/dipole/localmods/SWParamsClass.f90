MODULE SWParamsClass
!
! SWParamsClass.f90
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
       INTEGER       :: nCutoff
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
       ! DipoleParameters
       REAL(prec)    :: vMax
       REAL(prec)    :: x0
       REAL(prec)    :: x1
       REAL(prec)    :: y0
       REAL(prec)    :: y1
       REAL(prec)    :: L
       ! SpongeParameters
       REAL(prec)    :: Lsponge
       REAL(prec)    :: rFacMax
       ! Shelf
       REAL(prec)    :: hMin
       REAL(prec)    :: hMax
       REAL(prec)    :: Lshelf
       REAL(prec)    :: dL
       REAL(prec)    :: steepeningCenter
       REAL(prec)    :: steepeningZoneLength
            
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
       INTEGER       :: nCutoff
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
       ! JetParameters
       REAL(prec)    :: vMax
       REAL(prec)    :: x0
       REAL(prec)    :: x1
       REAL(prec)    :: y0
       REAL(prec)    :: y1
       REAL(prec)    :: L
       ! SpongeParameters
       REAL(prec)    :: Lsponge
       REAL(prec)    :: rFacMax
       ! Shelf
       REAL(prec)    :: hMin
       REAL(prec)    :: hMax
       REAL(prec)    :: Lshelf
       REAL(prec)    :: dL
       REAL(prec)    :: steepeningCenter
       REAL(prec)    :: steepeningZoneLength
       
      NAMELIST / ModelForm / ModelFormulation, nCutoff
      NAMELIST / TimeManagement / dt, iterInit, nTimeSteps, dumpFreq
      NAMELIST / SpaceManagement / SpecMeshFile, polyDeg, nXElem, nYElem, nPlot, xScale, yScale
      NAMELIST / PhysicalParameters / g, f0, betaX, betaY, linearDrag
      NAMELIST / DipoleParameters / vMax, x0, x1, y0, y1, L
      NAMELIST / SpongeParameters / Lsponge, rFacMax
      NAMELIST / ShelfParameters / hMin, hMax, Lshelf, dL, steepeningCenter, steepeningZoneLength
      
      ! MODEL_FORM
      ModelFormulation = 'LINEAR'
      nCutoff = 3
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
      xScale = 10.0_prec**(6)
      yScale = 10.0_prec**(6) 
      ! PHYSICAL_PARAMETERS
      g = 9.81        ! m/sec^2
      f0 = 10.0_prec**(-4) ! 1/sec
      betaX = ZERO         ! 1/(m*sec)
      betaY = ZERO         ! 1/(m*sec)
      linearDrag = ZERO         ! 1/sec
      ! DipoleParameters
      vMax = 1.0_prec
      x0   = HALF*xScale
      y0   = (HALF+0.05_prec)*yScale
      x1   = HALF*xScale
      y1   = (HALF-0.05_prec)*yScale
      L    = 0.05_prec*xScale
      ! SpongeParameter
      Lsponge = L
      rFacMax = f0
      ! ShelfParameters
      hMin   = 200.0_prec
      hMax   = 3000.0_prec
      Lshelf = L
      dL     = ZERO
      steepeningCenter = 1000.0_prec*10.0_prec**(3)
      steepeningZoneLength = 500.0_prec*10.0_prec**(3)

      OPEN( UNIT = NEWUNIT(nUnit), FILE = 'shallowwater.params')
         READ( UNIT = nUnit, NML = ModelForm )
         READ( UNIT = nUnit, NML = TimeManagement )
         READ( UNIT = nUnit, NML = SpaceManagement )
         READ( UNIT = nUnit, NML = PhysicalParameters )
         READ( UNIT = nUnit, NML = DipoleParameters )
         READ( UNIT = nUnit, NML = SpongeParameters )
         READ( UNIT = nUnit, NML = ShelfParameters )
      CLOSE( UNIT = nUnit ) 

      ! Sanity check - PRINT the results of the namelist READ to screen
         WRITE( UNIT = *, NML = ModelForm )
         WRITE( UNIT = *, NML = TimeManagement )
         WRITE( UNIT = *, NML = SpaceManagement )
         WRITE( UNIT = *, NML = PhysicalParameters )
         WRITE( UNIT = *, NML = DipoleParameters )
         WRITE( UNIT = *, NML = SpongeParameters )
         WRITE( UNIT = *, NML = ShelfParameters )
      
      ! MODEL_FORM
      thisParam % ModelFormulation = GetFlagForChar( ModelFormulation )
      thisParam % nCutoff          = nCutoff
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
      ! DipoleParameters
      thisParam % vMax = vMax
      thisParam % x0   = x0
      thisParam % x1   = x1
      thisParam % y0   = y0
      thisParam % y1   = y1
      thisParam % L    = L
      ! SpongeParameter
      thisParam % Lsponge = Lsponge
      thisParam % rFacMax = rFacMax
      ! ShelfParameters
      thisParam % hMin   = hMin
      thisParam % hMax   = hMax
      thisParam % Lshelf = Lshelf
      thisParam % dL     = dL
      thisParam % steepeningCenter = steepeningCenter
      thisParam % steepeningZoneLength = steepeningZoneLength
      
 END SUBROUTINE Build_SWParams

END MODULE SWParamsClass
