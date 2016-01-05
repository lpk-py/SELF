MODULE GeometryParamsClass
!
! GeometryParamsClass.f90
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


    TYPE GeometryParams
       ! SPACE_MANAGEMENT
       CHARACTER(50) :: SpecMeshFile
       INTEGER       :: polyDeg
       INTEGER       :: geomPolyDeg
       INTEGER       :: nPlot
       REAL(prec)    :: dxPlot
       REAL(prec)    :: xScale
       REAL(prec)    :: yScale
      
       CONTAINS

       PROCEDURE :: Build => Build_GeometryParams

    END TYPE GeometryParams 
 

 CONTAINS


 SUBROUTINE Build_GeometryParams( thisParam )
 ! S/R Build
 !
 !
 ! ==============================================!
 ! DECLARATIONS
   CLASS( GeometryParams ), intent(out) :: thisParam
   ! LOCAL
   INTEGER :: nUnit
       ! SPACE_MANAGEMENT
       CHARACTER(50) :: SpecMeshFile
       INTEGER       :: polyDeg
       INTEGER       :: nPlot
       REAL(prec)    :: dxPlot
       REAL(prec)    :: xScale
       REAL(prec)    :: yScale


      NAMELIST / SPACE_MANAGEMENT / SpecMeshFile, polyDeg, nPlot, xScale, yScale
 
      SpecMeshFile = nada
      polyDeg = 5
      nPlot = 10
      xScale = ONE
      yScale = ONE 

      ! Reading in the namelist FILE

      OPEN( UNIT = NEWUNIT(nUnit), FILE = 'geometry.params')
         READ( UNIT = nUnit, NML = SPACE_MANAGEMENT )
      CLOSE( UNIT = nUnit ) 

      ! Sanity check - PRINT the results of the namelist READ to screen
      WRITE( UNIT = *, NML = SPACE_MANAGEMENT )
      
      ! SPACE_MANAGEMENT (Default to isoparametric elements - if overintegration is USEd, default is to double the number of points )
      ! Plotting is defaulted to 10 evenly spaced points
      thisParam % SpecMeshFile = SpecMeshFile
      thisParam % polyDeg = polyDeg
      thisParam % nPlot = nPlot
      thisParam % dxPlot = 2.0_prec/REAL(nPlot,prec)
      thisParam % xScale = xScale
      thisParam % yScale = yScale


 END SUBROUTINE Build_GeometryParams

END MODULE GeometryParamsClass
