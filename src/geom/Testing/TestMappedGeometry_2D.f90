PROGRAM TestMappedGeometry_2D



 ! src/common/
 USE ModelPrecision
 USE ConstantsDictionary
 ! src/interp/
 USE Lagrange_2D_Class
 ! src/nodal/
 USE NodalStorage_2D_Class
 ! src/geom/
 USE CurveClass_2D
 USE GeometryParamsClass
 USE MappedGeometryClass_2D
 
 IMPLICIT NONE
 
 TYPE( NodalStorage_2D )   :: nodal
 TYPE( MappedGeometry_2D ) :: geometry, plottingGeometry
 TYPE( Curve_2D )          :: boundaryCurves(1:4)
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
                        
    ! Build the boundary curves
    CALL ConstructCurves( boundaryCurves, nodal % interp )
    
    ! Build the Geometry
    CALL geometry % Build( nodal % interp, boundaryCurves )
    
    ! Write to a tecplot file
    CALL geometry % WriteTecplot( )
    
    ! Test out the computational coordinate solver
    CALL geometry % CalculateComputationalCoordinates( nodal % interp, ZERO, 1.25_prec, s, p, successful )
    
    IF( successful )THEN
       PRINT*,'Computational coordinate search successful '
       PRINT*, 'Coordinates (s,p) : (', s,',',p,')'
    ENDIF
    
    CALL geometry % Trash( )

 CONTAINS
 
 SUBROUTINE ConstructCurves( boundCurves, interp )
    TYPE( Curve_2D ), INTENT(out)   :: boundCurves(1:4)
    TYPE( Lagrange_2D ), INTENT(in) :: interp
    !
    REAL(prec) :: r(0:interp % nS, 1:2)
    REAL(prec) :: s(0:interp % nS), p(0:interp % nP)
    INTEGER    :: nS, nP
    
       CALL interp % GetNodes( s, p )
       CALL interp % GetNumberOfNodes( nS, nP )
       
       r = SouthCurve( s, nS )
       CALL boundCurves(1) % Build( r(:,1), r(:,2), s )
       
       r = EastCurve( p, nP )
       CALL boundCurves(2) % Build( r(:,1), r(:,2), p )
       
       r = NorthCurve( s, nS )
       CALL boundCurves(3) % Build( r(:,1), r(:,2), s )
       
       r = WestCurve( p, nP )
       CALL boundCurves(4) % Build( r(:,1), r(:,2), p )
    
 END SUBROUTINE ConstructCurves
 
 FUNCTION SouthCurve( s, N ) RESULT( r )
    INTEGER    :: N
    REAL(prec) :: s(0:N)
    REAL(prec) :: r(0:N,1:2)
    !LOCAL
    INTEGER :: i
    
       N = UBOUND(s,1)
       
       r(:,1) = HALF*(s+ONE) + HALF
       r(:,2) = ZERO
    
 END FUNCTION SouthCurve
!
 FUNCTION NorthCurve( s, N ) RESULT( r )
    INTEGER    :: N
    REAL(prec) :: s(0:N)
    REAL(prec) :: r(0:N,1:2)
    !LOCAL
    INTEGER :: i
    
       N = UBOUND(s,1)
       
       r(:,1) = -( HALF*(s+ONE) + HALF )
       r(:,2) = ZERO
       
 END FUNCTION NorthCurve
!
 FUNCTION EastCurve( s, N ) RESULT( r )
    INTEGER    :: N
    REAL(prec) :: s(0:N)
    REAL(prec) :: r(0:N,1:2)
    !LOCAL
    INTEGER :: i
    
       N = UBOUND(s,1)
       
       r(:,1) = (ONE+HALF)*cos( pi*(s+ONE)*HALF )
       r(:,2) = (ONE+HALF)*sin( pi*(s+ONE)*HALF )
    
 END FUNCTION EastCurve
!
 FUNCTION WestCurve( s, N ) RESULT( r )
    INTEGER    :: N
    REAL(prec) :: s(0:N)
    REAL(prec) :: r(0:N,1:2)
    !LOCAL
    INTEGER :: i
    
       N = UBOUND(s,1)
       
       r(:,1) = (HALF)*cos( pi*(s+ONE)*HALF )
       r(:,2) = (HALF)*sin( pi*(s+ONE)*HALF )
    
 END FUNCTION WestCurve
END PROGRAM TestMappedGeometry_2D
