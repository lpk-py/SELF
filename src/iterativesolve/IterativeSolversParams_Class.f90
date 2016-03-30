
MODULE IterativeSolversParams_Class


 USE ModelPrecision
 USE CommonRoutines
 USE ConstantsDictionary
 USE ModelFlags
 

 IMPLICIT NONE


    TYPE IterativeSolversParams
       ! MODEL_FORM
       INTEGER       :: MaximumIterates
       REAL(prec)    :: tolerance
       INTEGER       :: nDOF


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
   ! MODEL_FORM
   INTEGER    :: MaximumIterates
   REAL(prec) :: tolerance
   INTEGER    :: nDOF


      ! Set the namelists
      NAMELIST / SolverCriteria / MaximumIterates, tolerance, nDOF
      
 
      ! Set the default parameters
      ! SolverCriteria
      MaximumIterates = 1000   ! Max number of conjugate gradient iterations
      tolerance = 10.0_prec**(-8) ! conjugate gradient residual tolerance
      nDOF = 1

      ! Reading in the namelist FILE

      OPEN( UNIT = NEWUNIT(nUnit), FILE = 'itsolve.params')
         READ( UNIT = nUnit, NML = SolverCriteria )
      CLOSE( UNIT = nUnit ) 

      ! Sanity check - PRINT the results of the namelist READ to screen
      WRITE( UNIT = *, NML = SolverCriteria )

      ! Fill in the data structure
      thisParam % MaximumIterates = MaximumIterates
      thisParam % tolerance       = tolerance
      thisParam % nDOF            = nDOF
      
      


 END SUBROUTINE BuildParams

END MODULE IterativeSolversParams_Class
