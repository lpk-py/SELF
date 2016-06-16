PROGRAM TestQuadMesh



 ! src/common/
 USE ModelPrecision
 USE ConstantsDictionary
 ! src/nodal/
 USE NodalStorage_2D_Class
 ! src/geom/
 USE QuadMeshClass
 USE GeometryParamsClass
 
 IMPLICIT NONE
 
 TYPE( NodalStorage_2D )   :: nodal
 TYPE( QuadMesh )          :: mesh
 TYPE( GeometryParams )    :: params
 REAL(prec)                :: s, p
 LOGICAL                   :: successful
 
 
 
    ! Read in the parameters
    CALL params % Build( )
    
    ! Build an interpolant
    CALL nodal % Build( N = params % polyDeg, &
                        M = params % polyDeg, &
                        quadrature = GAUSS_LOBATTO,  &
                        approxForm = CG )

    
    ! Build the Geometry
      IF( TRIM( params % SpecMeshFile ) == nada )THEN
         PRINT*,' Loading default mesh.'
         CALL mesh % LoadDefaultMesh( nodal % interp, &
                                       params % nXelem, &
                                       params % nYelem )
      ELSE
      ! Builds the lateral mesh
         PRINT*, 'Reading mesh from '//trim(params % SpecMeshFile)//'.'
         CALL mesh % ReadSpecMeshFile( nodal % interp, &
                                       params % SpecMeshFile )
      ENDIF

    ! Find the unique edges using hash-table. I know, this is redundant if a SpecMesh file is read.
    ! However, we can still do a fresh edge search in order to determine the performance of the 
    ! hash-table edge search.
    CALL mesh % ConstructEdges( )


    ! Scale the Geometry
    CALL mesh % ScaleTheMesh( nodal % interp, params % xScale, params % yScale )
    ! Write to a tecplot file
    CALL mesh % WriteTecplot( )
   
    CALL mesh % Trash( )
    

 
END PROGRAM TestQuadMesh
