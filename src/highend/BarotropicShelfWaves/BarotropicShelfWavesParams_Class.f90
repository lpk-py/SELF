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


    TYPE BarotropicShelfWavesParams
       ! SolverCriteria
       INTEGER       :: MaximumIterates
       INTEGER       :: nEpairs
       INTEGER       :: nMaxToSolve
       REAL(prec)    :: tolerance
       ! SpaceManagement
       INTEGER       :: polyDeg
       INTEGER       :: nElems
       INTEGER       :: nPlot
       REAL(prec)    :: dxPlot
       REAL(prec)    :: xScale
       ! Physical
       REAL(prec)    :: f0
       REAL(prec)    :: innerSlopeLength
       REAL(prec)    :: outerSlopeLength
       REAL(prec)    :: shelfwidth
       REAL(prec)    :: hmin
       REAL(prec)    :: shelfDepth
       REAL(prec)    :: offshoreDepth


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
       ! SpaceManagement
       INTEGER       :: polyDeg
       INTEGER       :: nElems
       INTEGER       :: nPlot
       REAL(prec)    :: dxPlot
       REAL(prec)    :: xScale
       ! Physical
       REAL(prec)    :: f0
       REAL(prec)    :: innerSlopeLength
       REAL(prec)    :: outerSlopeLength
       REAL(prec)    :: shelfwidth
       REAL(prec)    :: hmin
       REAL(prec)    :: shelfDepth
       REAL(prec)    :: offshoreDepth


      ! Set the namelists
      NAMELIST / SolverCriteria / MaximumIterates, nEpairs, nMaxToSolve, tolerance
      NAMELIST / SpaceManagement / polyDeg, nElems, nPlot, xScale
      NAMELIST / Physical / f0, innerSlopeLength, outerSlopeLength, shelfWidth, hmin, &
                            shelfDepth, offshoreDepth
 
      ! Set the default parameters
      ! SolverCriteria
      MaximumIterates = 1000   ! Max number of conjugate gradient iterations
      nEpairs         = 10
      nMaxToSolve     = 15
      tolerance       = 10.0_prec**(-8) ! conjugate gradient residual tolerance

      ! SpaceManagement 
      polyDeg = 10
      nElems  = 10
      nPlot   = 10
      xScale  = 2.0E6

      ! Physical
      f0 = 1.0E-4
      innerSlopeLength = 1.0E5
      outerSlopeLength = 1.0E5
      shelfwidth       = 10.0_prec*xScale
      hmin             = 50.0_prec
      shelfDepth       = 2050.0_prec
      offshoreDepth    = shelfDepth
      ! Reading in the namelist FILE

      OPEN( UNIT = NEWUNIT(nUnit), FILE = 'runtime.params')
         READ( UNIT = nUnit, NML = SolverCriteria )
         READ( UNIT = nUnit, NML = SpaceManagement )
         READ( UNIT = nUnit, NML = Physical )
      CLOSE( UNIT = nUnit ) 

      ! Sanity check - PRINT the results of the namelist READ to screen
      WRITE( UNIT = *, NML = SolverCriteria )
      WRITE( UNIT = *, NML = SpaceManagement )
      WRITE( UNIT = *, NML = Physical )

      ! Fill in the data structure
      ! SolverCriteria
      thisParam % MaximumIterates = MaximumIterates
      thisParam % nEpairs         = nEpairs
      thisParam % nMaxToSolve     = nMaxToSolve
      thisParam % tolerance       = tolerance
      
      ! SpaceManagement
      thisParam % polyDeg = polyDeg
      thisParam % nElems  = nElems
      thisParam % nPlot   = nPlot
      thisParam % dxPlot  = 2.0_prec/REAL(nPlot,prec)
      thisParam % xScale  = xScale

      ! Physical
      thisParam % f0 = f0
      thisParam % innerSlopeLength = innerSlopeLength
      thisParam % outerSlopeLength = outerSlopeLength
      thisParam % shelfwidth       = shelfwidth
      thisParam % hmin             = hmin
      thisParam % shelfDepth       = shelfDepth
      thisParam % offshoreDepth    = offshoreDepth

 END SUBROUTINE BuildParams_BarotropicShelfWaves

END MODULE BarotropicShelfWavesParams_Class
