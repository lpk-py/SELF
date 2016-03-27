PROGRAM TestSegmentMesh

USE ModelPrecision
USE ModelFlags

USE NodalStorage_1D_Class
USE SegmentMeshClass

IMPLICIT NONE
INTEGER, PARAMETER :: nS = 10
INTEGER, PARAMETER :: nEl = 10

TYPE(SegmentMesh)     :: thisMesh
TYPE(NodalStorage_1D) :: nodal

   CALL nodal % Build( nS, GAUSS_LOBATTO, CG )
   CALL thisMesh % LoadDefaultMesh( nodal % interp, nEl )
   
   CALL thisMesh % ScaleTheMesh( TWO )

   CALL thisMesh % WriteTecplot( )
   
   CALL thisMesh % Trash( )
   CALL nodal % Trash( )   
   
END PROGRAM TestSegmentMesh
