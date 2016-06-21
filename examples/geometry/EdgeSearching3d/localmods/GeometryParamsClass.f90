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
       ! SpaceManagement
       CHARACTER(100) :: SpecMeshFile
       INTEGER        :: nXElem
       INTEGER        :: nYElem
       INTEGER        :: nZElem
       INTEGER        :: polyDeg
       INTEGER        :: geomPolyDeg
       INTEGER        :: nPlot
       REAL(prec)     :: dxPlot
       REAL(prec)     :: xScale
       REAL(prec)     :: yScale
       REAL(prec)     :: zScale
      
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
       ! SpaceManagement
       CHARACTER(100) :: SpecMeshFile
       INTEGER        :: nXElem
       INTEGER        :: nYElem
       INTEGER        :: nZElem
       INTEGER       :: polyDeg
       INTEGER       :: geomPolyDeg
       INTEGER       :: nPlot
       REAL(prec)    :: dxPlot
       REAL(prec)    :: xScale
       REAL(prec)    :: yScale
       REAL(prec)    :: zScale


      NAMELIST / SpaceManagement / SpecMeshFile, nXElem, nYElem, nZElem, polyDeg,  &
                                   geomPolyDeg, nPlot, xScale, yScale, zScale
 
      SpecMeshFile = nada
      nXElem = 2
      nYElem = 2
      nZelem = 2
      polyDeg = 5
      geomPolyDeg = 5
      nPlot = 10
      xScale = ONE
      yScale = ONE 
      zScale = ONE
     
      

      ! Reading in the namelist FILE

      OPEN( UNIT = NEWUNIT(nUnit), FILE = 'geometry.params')
         READ( UNIT = nUnit, NML = SpaceManagement )
      CLOSE( UNIT = nUnit ) 

      ! Sanity check - PRINT the results of the namelist READ to screen
      WRITE( UNIT = *, NML = SpaceManagement )
      
      ! SpaceManagement (Default to isoparametric elements - if overintegration is USEd, default is to double the number of points )
      ! Plotting is defaulted to 10 evenly spaced points
      thisParam % SpecMeshFile = SpecMeshFile
      thisParam % nXElem = nXElem
      thisParam % nYElem = nYElem 
      thisParam % nZElem = nZElem
      thisParam % polyDeg = polyDeg
      thisParam % nPlot = nPlot
      thisParam % dxPlot = 2.0_prec/REAL(nPlot,prec)
      thisParam % xScale = xScale
      thisParam % yScale = yScale
      thisParam % zScale = zScale

 END SUBROUTINE Build_GeometryParams

END MODULE GeometryParamsClass
