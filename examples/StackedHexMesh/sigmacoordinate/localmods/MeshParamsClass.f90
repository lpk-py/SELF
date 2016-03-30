MODULE MeshParamsClass
!
! MeshParamsClass.f90
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


    TYPE MeshParams
      ! SpaceManagement
       CHARACTER(50) :: SpecMeshFile
       INTEGER       :: geomPolyDeg
       INTEGER       :: polyDeg
       INTEGER       :: nXElem
       INTEGER       :: nYElem
       INTEGER       :: nZElem
       INTEGER       :: nPlot
       REAL(prec)    :: dxPlot
       REAL(prec)    :: xScale
       REAL(prec)    :: yScale
       REAL(prec)    :: zScale
      
       CONTAINS

       PROCEDURE :: Build => Build_MeshParams

    END TYPE MeshParams 
 

 CONTAINS


 SUBROUTINE Build_MeshParams( thisParam )
 ! S/R Build
 !
 !
 ! ==============================================!
 ! DECLARATIONS
   CLASS( MeshParams ), intent(out) :: thisParam
   ! LOCAL
   INTEGER :: nUnit
       ! SpaceManagement
       CHARACTER(50) :: SpecMeshFile
       INTEGER       :: geomPolyDeg
       INTEGER       :: polyDeg
       INTEGER       :: nXElem
       INTEGER       :: nYElem
       INTEGER       :: nZElem
       INTEGER       :: nPlot
       REAL(prec)    :: dxPlot
       REAL(prec)    :: xScale
       REAL(prec)    :: yScale
       REAL(prec)    :: zScale

      NAMELIST / SpaceManagement / SpecMeshFile, geomPolyDeg, polyDeg, nXElem, nYElem, nZelem, &
                                   nPlot, xScale, yScale, zScale

      ! SpaceManagement
      SpecMeshFile = nada
      geomPolyDeg = 1
      polyDeg = 5
      nXElem = 5
      nYElem = 5
      nPlot = 10
      xScale = ONE
      yScale = ONE 
      zScale = ONE
      ! Reading in the namelist FILE

      OPEN( UNIT = NEWUNIT(nUnit), FILE = 'mesh.params')
         READ( UNIT = nUnit, NML = SpaceManagement )
      CLOSE( UNIT = nUnit ) 

      ! Sanity check - PRINT the results of the namelist READ to screen
      WRITE( UNIT = *, NML = SpaceManagement )

      ! SpaceManagement
      thisParam % SpecMeshFile = SpecMeshFile
      thisParam % geomPolyDeg = geomPolyDeg
      thisParam % polyDeg = polyDeg
      thisParam % nXElem = nXElem
      thisParam % nYElem = nYElem 
      thisParam % nZElem = nZElem
      thisParam % nPlot = nPlot
      thisParam % dxPlot = 2.0_prec/REAL(nPlot,prec)
      thisParam % xScale = xScale
      thisParam % yScale = yScale
      thisParam % zScale = zScale

      
      
 END SUBROUTINE Build_MeshParams

END MODULE MeshParamsClass
