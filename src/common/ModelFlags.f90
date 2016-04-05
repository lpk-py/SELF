MODULE ModelFlags

  !*************************************************************!
  ! ------------------- Edge/Face Numbering --------------------!
  ! ************************************************************!
  ! This block of flags pertains to the integer associates with !
  ! assigning boundary edge and face ordering of an element     !
  !                                                             !
  ! ------------------------------------------------------------!
    INTEGER, PARAMETER :: left   = 1
    INTEGER, PARAMETER :: right  = 1
    INTEGER, PARAMETER :: south  = 1
    INTEGER, PARAMETER :: east   = 2
    INTEGER, PARAMETER :: north  = 3
    INTEGER, PARAMETER :: west   = 4
    INTEGER, PARAMETER :: bottom = 5
    INTEGER, PARAMETER :: top    = 6
    
    INTEGER, PARAMETER :: nHexFaces  = 6
    INTEGER, PARAMETER :: nHexNodes  = 8
    INTEGER, PARAMETER :: nQuadEdges = 4
    INTEGER, PARAMETER :: nQuadNodes = 4

  !*************************************************************!
  ! ----------------- BOUNDARY CONDITION FLAGS -----------------!
  ! ************************************************************!
  ! This block of flags pertain to boundary condition flags     !
  ! that make it easy to reference which boundary condition     !
  ! to USE. These flags improve code readability and are        !
  ! easy to USE.                                                !
  !                                                             !
  ! For any boundary condition, it is imperative that the flag  !
  ! be negative.                                                !
  !                                                             !
  ! "Cornernodes" are an exception to this rule                 !
  !                                                             !
  ! The convention for any software within the                  !
  ! SCHOONER package is that the elements of a mesh and have    !
  ! a positive element ID. In the "edge" information for a mesh !
  ! the secondary element ID can be an actual element ID or set !
  ! to a boundary condition flag. A boundary condition check is !
  ! usually done by  passing through a conditional              !
  !                                                             !
  !           " if( secondaryElement < 0 )"                     !
  !                                                             !
  ! ------------------------------------------------------------!
  !==============================================!
  ! ------------- Node Specifiers -------------- !
  !==============================================!
   INTEGER, PARAMETER :: INTERIOR = 1
   INTEGER, PARAMETER :: BOUNDARY = 0
  !
  !==============================================!
  ! --------- Discontinuous Galerkin ----------- !
  !==============================================!
   INTEGER, PARAMETER :: NO_NORMAL_FLOW = -100
   INTEGER, PARAMETER :: RADIATION = -101
   INTEGER, PARAMETER :: PRESCRIBED = -102
   INTEGER, PARAMETER :: InflowOne = -103
   INTEGER, PARAMETER :: InflowTwo = -104
   INTEGER, PARAMETER :: SEA_FLOOR = -105
  !
  !==============================================!
  ! ---------- Continuous Galerkin ------------- !
  !==============================================!
   INTEGER, PARAMETER :: DIRICHLET = -200
   INTEGER, PARAMETER :: HOMOGENEOUS_NEUMANN = -201
   INTEGER, PARAMETER :: ROBIN = -202
   INTEGER, PARAMETER :: ROBIN_FORCED = -204
   INTEGER, PARAMETER :: INHOMOGENEOUS_NEUMANN = -203
   INTEGER, PARAMETER :: NEUMANN = -100 ! = NO_NORMAL_FLOW
   INTEGER, PARAMETER :: NEUMANN_WALL = -206
  !
  !==============================================!
  !
  !*************************************************************!
  ! ----------------- MODEL FORMULATION FLAGS ------------------!
  ! ************************************************************!
  ! This block of flags pertains to those which are USEd to     !
  ! change the formulation of a particular model, such as the   !
  ! "ShallowWater" class which offers up three different flavors!
  ! of the model. Additionally, multipurpose flags which are    !
  ! USEd for specifying the TYPE of quadrature are given here.  !
  !                                                             !
  ! ------------------------------------------------------------!
  !==============================================!
  ! --------------- Quadrature------------------ !
  !==============================================!
   INTEGER, PARAMETER :: GAUSS = 1
   INTEGER, PARAMETER :: GAUSS_LOBATTO = -1
   INTEGER, PARAMETER :: DG = 2000
   INTEGER, PARAMETER :: CG = 2001
   INTEGER, PARAMETER :: LEGENDRE_BASIS = 2100
   INTEGER, PARAMETER :: CHEBYSHEV_BASIS = 2101
  !
  !==============================================!
  ! ---------------- Geometry ------------------ !
  !==============================================!
   INTEGER, PARAMETER :: STANDARD_FLUX = 300       ! USEs "CLASSIC" calculation of contravariant basic vectors or Curl-Invariant (Kopriva, 2006)
   INTEGER, PARAMETER :: CROSS_PRODUCT_FLUX = 301  ! USEs Kopriva 2006 contravariant basis vectors and cross product identity for contravariant flux (MUST USE CURL_CONTRAVARIANT)
   INTEGER, PARAMETER :: TRANSFINITE_METRICS = 200 ! compute metric terms from transfinite interpolation
   INTEGER, PARAMETER :: INTERP_METRICS = 201      ! compute metric terms from interpolation of mesh locations
   INTEGER, PARAMETER :: CLASSIC_CONTRAVARIANT = 210
   INTEGER, PARAMETER :: CURL_CONTRAVARIANT = 211
  !
  !==============================================!
  ! --------------- Model forms ---------------- !
  !==============================================!
   INTEGER, PARAMETER :: LINEAR = 99               ! Solves the Linearized form of the shallow water equations
   INTEGER, PARAMETER :: NONLINEAR = 100
   INTEGER, PARAMETER :: SKEW_SYMMETRIC = 101      ! Solves the Skew-symmetric form of the shallow water equations
   INTEGER, PARAMETER :: CONSERVATIVE = 102        ! Solves the Conservative form of the shallow water equations
  !==============================================!
  !==============================================!
  ! --------- Special plotting flags ----------- !
  !==============================================!
   INTEGER, PARAMETER :: surf2din3d = 500
   INTEGER, PARAMETER :: keepItFlat = 501


! Misc. INTEGER and CHARACTER flag definitions
  INTEGER, PARAMETER      :: NONE = 0
  CHARACTER(1), PARAMETER :: nada = ' ' 


END MODULE ModelFlags
