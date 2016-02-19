
MODULE IterativeSolversParams_Class


 USE ModelPrecision
 USE CommonRoutines
 USE ConstantsDictionary
 USE ModelFlags
 

 IMPLICIT NONE


    TYPE IterativeSolversParams
       ! SolverCriteria
       INTEGER       :: MaximumIterates
       REAL(prec)    :: tolerance
       INTEGER       :: nDOF
       ! FV1D
       REAL(prec)    :: xL, xR, dx, kappa
       INTEGER       :: nTimingCycles

       CONTAINS

       PROCEDURE :: Build => BuildParams

    END TYPE IterativeSolversParams 
 

 CONTAINS


 SUBROUTINE BuildParams( thisParam )
 ! S/R BuildParams
 !
 !
 ! ==============================================!
 ! DECLARATIONS
   CLASS( IterativeSolversParams ), intent(out) :: thisParam
   ! LOCAL
   INTEGER :: nUnit
   ! SolverCriteria
   INTEGER    :: MaximumIterates
   REAL(prec) :: tolerance
   ! FV1D
   REAL(prec) :: xL, xR, dx, kappa
   INTEGER    :: nTimingCycles

      ! Set the namelists
      NAMELIST / SolverCriteria / MaximumIterates, tolerance
      
      NAMELIST/ FV1D / xL, xR, dx, kappa, nTimingCycles
 
      ! Set the default parameters
      ! SolverCriteria
      MaximumIterates = 1000   ! Max number of conjugate gradient iterations
      tolerance = 10.0_prec**(-8) ! conjugate gradient residual tolerance

      ! FV1D
      xL    = ZERO
      xR    = ONE
      dx    = 0.02_prec
      kappa = ONE
      nTimingCycles = 1
      
      ! Reading in the namelist FILE

      OPEN( UNIT = NEWUNIT(nUnit), FILE = 'itsolve.params')
         READ( UNIT = nUnit, NML = SolverCriteria )
         READ( UNIT = nUnit, NML = FV1D )
      CLOSE( UNIT = nUnit ) 

      ! Sanity check - PRINT the results of the namelist READ to screen
      WRITE( UNIT = *, NML = SolverCriteria )

      ! Fill in the data structure
      thisParam % MaximumIterates = MaximumIterates
      thisParam % tolerance       = tolerance
      
      thisParam % xL = xL
      thisParam % xR = xR
      thisParam % dx = dx
      thisParam % kappa = kappa
      thisParam % nTimingCycles = nTimingCycles
  
      thisParam % nDOF = INT( (thisParam % xR - thisParam % xL)/thisParam % dx ) 
      


 END SUBROUTINE BuildParams

END MODULE IterativeSolversParams_Class
