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
 TYPE( QuadMesh )          :: linmesh
 TYPE( GeometryParams )    :: params
 REAL(prec)                :: s, p
 LOGICAL                   :: successful
 INTEGER, ALLOCATABLE      :: linToO(:,:,:,:), oToLin(:,:,:,:)
 
 
 
    ! Read in the parameters
    CALL params % Build( )
    
    ! Build an interpolant
    CALL nodal % Build( N = params % polyDeg, &
                        M = params % polyDeg, &
                        quadrature = GAUSS_LOBATTO,  &
                        approxForm = CG )

    
    ! Build the Geometry
    CALL mesh % LoadDefaultMesh( nodal % interp, 2, 2 )
    !CALL mesh % ReadSpecMeshFile( nodal % interp, params % SpecMeshFile )
    CALL linmesh % ConstructLinearMeshFromOther( mesh, linToO, oToLin )
    ! Scale the Geometry
!    CALL mesh % ScaleTheMesh( nodal % interp, params % xScale, params % yScale )
    ! Write to a tecplot file
    CALL mesh % WriteTecplot( )
    CALL linmesh % WriteTecplot( 'linmesh' )
   
    CALL mesh % Trash( )
    CALL linmesh % Trash( )
    
    DEALLOCATE( linToO, oToLin )
    ! Test the Default mesh
    !CALL mesh % LoadDefaultMesh( nodal % interp )
    !CALL mesh % WriteTecplot('DefaultMesh')
    !CALL mesh % Trash( )
    

 
END PROGRAM TestQuadMesh
