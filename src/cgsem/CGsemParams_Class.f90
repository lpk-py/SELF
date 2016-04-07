MODULE CGsemParams_Class
!
! CGsemParamsClass.f90
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


    TYPE CGsemParams
       ! SolverCriteria
       INTEGER       :: MaximumIterates
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
      


       CONTAINS

       PROCEDURE :: Build => BuildParams_CGsem

    END TYPE CGsemParams 
 

 CONTAINS


 SUBROUTINE BuildParams_CGsem( thisParam )
 ! S/R BuildParams
 !
 !
 ! ==============================================!
 ! DECLARATIONS
   CLASS( CGsemParams ), intent(out) :: thisParam
   ! LOCAL
   INTEGER :: nUnit
   
       ! SolverCriteria
       INTEGER       :: MaximumIterates
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

      ! Set the namelists
      NAMELIST / SolverCriteria / MaximumIterates, pcIterates, tolerance, pcTolerance, pcDegree

      NAMELIST / SpaceManagement / SpecMeshFile, polyDeg, nXElem, nYElem, nPlot, xScale
 
      ! Set the default parameters
      ! SolverCriteria
      MaximumIterates = 500   ! Max number of conjugate gradient iterations
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
      xScale       = ONE

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


 END SUBROUTINE BuildParams_CGsem

END MODULE CGsemParams_Class
