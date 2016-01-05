MODULE ModelFlags



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
  !
  !==============================================!
  ! ---------- Continuous Galerkin ------------- !
  !==============================================!
   INTEGER, PARAMETER :: DIRICHLET = -100
   INTEGER, PARAMETER :: HOMOGENEOUS_NEUMANN = -101
   INTEGER, PARAMETER :: ROBIN = -102
   INTEGER, PARAMETER :: ROBIN_FORCED = -104
   INTEGER, PARAMETER :: INHOMOGENEOUS_NEUMANN = -103
   INTEGER, PARAMETER :: NEUMANN = -105
   INTEGER, PARAMETER :: NEUMANN_WALL = -106
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
