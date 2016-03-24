MODULE BurgersParamsClass
!
! BurgersParamsClass.f90
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


    TYPE BurgersParams
      ! TimeManagement
       REAL(prec)    :: dt
       INTEGER       :: iterInit
       INTEGER       :: nTimeSteps
       INTEGER       :: dumpFreq
       ! SpaceManagement
       INTEGER       :: polyDeg
       INTEGER       :: nCutoff
       INTEGER       :: nElems
       REAL(prec)    :: dx
       INTEGER       :: nPlot
       REAL(prec)    :: dxPlot
       ! VelocityField
    !   REAL(prec)    :: 
      
       CONTAINS

       PROCEDURE :: Build => Build_BurgersParams

    END TYPE BurgersParams 
 

 CONTAINS


 SUBROUTINE Build_BurgersParams( thisParam )
 ! S/R Build
 !
 !
 ! ==============================================!
 ! DECLARATIONS
   CLASS( BurgersParams ), intent(out) :: thisParam
   ! LOCAL
   INTEGER :: nUnit
       ! TimeManagement
       REAL(prec)    :: dt
       INTEGER       :: iterInit
       INTEGER       :: nTimeSteps
       INTEGER       :: dumpFreq
       ! SpaceManagement
       INTEGER       :: polyDeg
       INTEGER       :: nCutoff
       INTEGER       :: nElems
       REAL(prec)    :: dx
       INTEGER       :: nPlot

      NAMELIST / TimeManagement / dt, iterInit, nTimeSteps, dumpFreq
      NAMELIST / SpaceManagement / polyDeg, nCutoff, nElems, dx, nPlot

      ! TimeManagement
      dt = ZERO
      iterInit = 0
      nTimeSteps = 0
      dumpFreq = 1
      ! SpaceManagement
      polyDeg = 15
      nCutoff = 8
      nElems = 5
      dx = 2.0_prec
      nPlot = 20
      
      ! Reading in the namelist FILE

      OPEN( UNIT = NEWUNIT(nUnit), FILE = 'burgers.params')
         READ( UNIT = nUnit, NML = TimeManagement )
         READ( UNIT = nUnit, NML = SpaceManagement )
      CLOSE( UNIT = nUnit ) 

      ! Sanity check - PRINT the results of the namelist READ to screen
      WRITE( UNIT = *, NML = TimeManagement )
      WRITE( UNIT = *, NML = SpaceManagement )
      

      ! TIME_MANAGEMENT
      thisParam % dt = dt
      thisParam % iterInit = iterInit
      thisParam % nTimeSteps = nTimeSteps
      thisParam % dumpFreq = dumpFreq
      ! SPACE_MANAGEMENT (Default to isoparametric elements - if overintegration is USEd, default is to double the number of points )
      ! Plotting is defaulted to 10 evenly spaced points
      thisParam % polyDeg = polyDeg
      thisParam % nCutoff = nCutoff
      thisParam % nElems  = nElems
      thisParam % dx      = dx
      thisParam % nPlot = nPlot
      thisParam % dxPlot = 2.0_prec/REAL(nPlot,prec)

      
      
 END SUBROUTINE Build_BurgersParams

END MODULE BurgersParamsClass
