PROGRAM TestHexMesh


 USE ModelFlags 
 USE NodalStorage_3D_Class
 USE HexMeshParamsClass
 USE HexMeshClass

 TYPE( NodalStorage_3D ) :: nodal 
 TYPE( HexMeshParams ) :: params
 TYPE( HexMesh )       :: mesh
 


    CALL params % Build( )

    ! Build an interpolant
    CALL nodal % Build( N = params % polyDeg, &
                        M = params % polyDeg, &
                        P = params % polyDeg, &
                        quadrature = GAUSS_LOBATTO,  &
                        approxForm = CG )

    CALL mesh % LoadDefaultMesh( nodal % interp, &
                                 params % nXelem, &
                                 params % nYelem, &
                                 params % nZelem  )

    CALL mesh % WriteTecplot( )
    CALL mesh % WriteMeshFile( GAUSS_LOBATTO, LEGENDRE_BASIS )

    CALL mesh % Trash( )
    CALL nodal % Trash( )

END PROGRAM TestHexMesh
