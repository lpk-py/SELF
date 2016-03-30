MODULE Elliptic1DParams_Class
!
! Elliptic1DParamsClass.f90
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


    TYPE Elliptic1DParams
       ! SolverCriteria
       INTEGER       :: MaximumIterates
       REAL(prec)    :: tolerance

       ! SpaceManagement
       INTEGER       :: polyDeg
       INTEGER       :: nElems
       INTEGER       :: nPlot
       REAL(prec)    :: dxPlot
       REAL(prec)    :: xScale
      


       CONTAINS

       PROCEDURE :: Build => BuildParams_Elliptic1D

    END TYPE Elliptic1DParams 
 

 CONTAINS


 SUBROUTINE BuildParams_Elliptic1D( thisParam )
 ! S/R BuildParams
 !
 !
 ! ==============================================!
 ! DECLARATIONS
   CLASS( Elliptic1DParams ), intent(out) :: thisParam
   ! LOCAL
   INTEGER :: nUnit
   
       ! SolverCriteria
       INTEGER       :: MaximumIterates
       REAL(prec)    :: tolerance

       ! SpaceManagement
       INTEGER       :: polyDeg
       INTEGER       :: nElems
       INTEGER       :: nPlot
!       REAL(prec)    :: dxPlot
       REAL(prec)    :: xScale

      ! Set the namelists
      NAMELIST / SolverCriteria / MaximumIterates, tolerance

      NAMELIST / SpaceManagement / polyDeg, nElems, nPlot, xScale
 
      ! Set the default parameters
      ! SolverCriteria
      MaximumIterates = 1000   ! Max number of conjugate gradient iterations
      tolerance = 10.0_prec**(-8) ! conjugate gradient residual tolerance

      ! SpaceManagement 
      polyDeg = 5
      nElems  = 5
      nPlot = 10
      xScale = ONE

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
      thisParam % tolerance       = tolerance
      
      ! SpaceManagement
      thisParam % polyDeg = polyDeg
      thisParam % nElems  = nElems
      thisParam % nPlot   = nPlot
      thisParam % dxPlot  = 2.0_prec/REAL(nPlot,prec)
      thisParam % xScale  = xScale


 END SUBROUTINE BuildParams_Elliptic1D

END MODULE Elliptic1DParams_Class
