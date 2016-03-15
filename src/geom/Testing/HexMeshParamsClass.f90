MODULE HexMeshParamsClass
!
! HexMeshParamsClass.f90
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


    TYPE HexMeshParams
       ! SPACE_MANAGEMENT
       INTEGER       :: polyDeg
       INTEGER       :: geomPolyDeg
       INTEGER       :: nPlot
       INTEGER       :: nXElem
       INTEGER       :: nYElem 
       INTEGER       :: nZElem 
       REAL(prec)    :: dxPlot
       REAL(prec)    :: xScale
       REAL(prec)    :: yScale
       REAL(prec)    :: zScale
      
       CONTAINS

       PROCEDURE :: Build => Build_GeometryParams

    END TYPE HexMeshParams 
 

 CONTAINS


 SUBROUTINE Build_GeometryParams( thisParam )
 ! S/R Build
 !
 !
 ! ==============================================!
 ! DECLARATIONS
   CLASS( HexMeshParams ), intent(out) :: thisParam
   ! LOCAL
   INTEGER :: nUnit
       ! SPACE_MANAGEMENT
       INTEGER       :: polyDeg
       INTEGER       :: geomPolyDeg
       INTEGER       :: nPlot
       INTEGER       :: nXElem
       INTEGER       :: nYElem 
       INTEGER       :: nZElem 
       REAL(prec)    :: dxPlot
       REAL(prec)    :: xScale
       REAL(prec)    :: yScale
       REAL(prec)    :: zScale


      NAMELIST / SPACE_MANAGEMENT / polyDeg, geomPolyDeg, nPlot, nXElem, nYElem, nZElem, xScale, yScale, zScale
 
      geomPolyDeg = 1
      polyDeg = 5
      nPlot = 10
      nXElem = 2
      nYElem = 2
      nZElem = 2
      xScale = ONE
      yScale = ONE
      zScale = ONE

      ! Reading in the namelist FILE

      OPEN( UNIT = NEWUNIT(nUnit), FILE = 'hexmesh.params')
         READ( UNIT = nUnit, NML = SPACE_MANAGEMENT )
      CLOSE( UNIT = nUnit ) 

      ! Sanity check - PRINT the results of the namelist READ to screen
      WRITE( UNIT = *, NML = SPACE_MANAGEMENT )
      
      ! SPACE_MANAGEMENT (Default to isoparametric elements - if overintegration is USEd, default is to double the number of points )
      ! Plotting is defaulted to 10 evenly spaced points
      thisParam % geomPolyDeg = geomPolyDeg
      thisParam % polyDeg = polyDeg
      thisParam % nPlot = nPlot
      thisParam % nXElem = nXElem 
      thisParam % nYElem = nYElem 
      thisParam % nZElem = nZElem 
      thisParam % dxPlot = 2.0_prec/REAL(nPlot,prec)
      thisParam % xScale = xScale
      thisParam % yScale = yScale
      thisParam % zScale = zScale


 END SUBROUTINE Build_GeometryParams

END MODULE HexMeshParamsClass
