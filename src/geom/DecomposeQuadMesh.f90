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
 REAL(prec), ALLOCATABLE   :: x(:,:), y(:,:), J(:,:), xc(:), yc(:), ws(:), wp(:)
 REAL(prec)                :: vol, swapFac
 INTEGER, ALLOCATABLE      :: mIDs(:)
 INTEGER                   :: nElems, nProc, nBlocks, iProc, iBlock, iEl, iter
 INTEGER                   :: iS, iP, nS, iSide, e2, nSwaps, mswap, nLocalSwaps, k 
 CHARACTER(5)              :: iterChar
 
      ! Read in the parameters
      CALL params % Build( )
      nS = params % polyDeg

      ALLOCATE( x(0:nS,0:nS), y(0:nS,0:nS), ws(0:nS), wp(0:nS), J(0:nS,0:nS) )

      ! Build an interpolant
      CALL nodal % Build( N = params % polyDeg, &
                          M = params % polyDeg, &
                          quadrature = GAUSS_LOBATTO,  &
                          approxForm = CG )

      CALL nodal % GetQuadratureWeights( ws, wp )

      ! Build the Geometry
      IF( params % SpecMeshFile == NADA )THEN
         CALL mesh % LoadDefaultMesh( nodal % interp, params % nXElem, params % nYElem )
      ELSE
         CALL mesh % ReadSpecMeshFile( nodal % interp, params % SpecMeshFile )
      ENDIF

      nElems = mesh % nElems
      nProc  = params % nProcesses
      ALLOCATE( mIDs(1:nElems) )
      nBlocks = nElems/nProc

      ALLOCATE( xc(1:nElems), yc(1:nElems) )
 
      ! Initialize the material ID's so that there is an equal number of elements per block
      iEl = 0
      DO iBlock = 1, nBlocks
         DO iProc = 1, nProc
            iEl = iEl + 1
            mIDs(iEl) = iProc
         ENDDO
      ENDDO
         
      iter = 0
      WRITE( iterChar, '(I5.5)' )iter
      CALL mesh % WriteMaterialTecplot( mIDs, 'procdecomp.'//iterChar )

      DO iEl = 1, nElems

         CALL mesh % GetPositions( iEl, x, y )
         CALL mesh % GetJacobian( iEl, J )
         xc(iEl) = ZERO
         yc(iEl) = ZERO
         vol     = ZERO

         DO iP = 0, nS
            DO iS = 0, nS

               vol = vol + J(iS,iP)*ws(iS)*wp(iP)
               xc(iEl)  = xc(iEl) + x(iS,iP)*J(iS,iP)*ws(iS)*wp(iP)
               yc(iEl)  = yc(iEl) + y(iS,iP)*J(iS,iP)*ws(iS)*wp(iP)

            ENDDO
         ENDDO

         xc(iEl) = xc(iEl)/vol
         yc(iEl) = yc(iEl)/vol

      ENDDO

      nSwaps = 0

      DO k = 1, params % nMaxSweeps
         nLocalSwaps = 0
         DO iEl = 1, nElems

            DO iSide = 1, 4
               CALL mesh % GetElementNeighbor( iEl, iSide, e2 )
               IF( e2 > 0 )THEN
                  swapFac = ( ( xc(e2) - xc(iEl) ) + ( yc(e2) - yc(iEl) ) )*( mIDs(e2) - mIDs(iEl) )
                  IF( swapFac < ZERO )THEN
                     nSwaps = nSwaps + 1
                     nLocalSwaps = nLocalSwaps + 1
                     mswap = mIDs(iEl)
                     mIDs(iEl) = mIDs(e2)
                     mIDs(e2)  = mswap

                     WRITE( iterChar, '(I5.5)' )nSwaps
                     CALL mesh % WriteMaterialTecplot( mIDs, 'procdecomp.'//iterChar )
                  ENDIF
               ENDIF
            ENDDO

         ENDDO
         IF( nLocalSwaps == 0 )THEN
            EXIT
         ENDIF
      ENDDO

      CALL mesh % WriteTecplot( )
      CALL mesh % Trash( )
      DEALLOCATE( x, y, ws, wp, J, xc, yc )
 
END PROGRAM TestQuadMesh
