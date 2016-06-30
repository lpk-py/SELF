! Lagrange_3D_Class.f90 (v 3.0)
! 
! Copyright 2015 Joseph Schoonover <schoonover.numerics@gmail.com>, The Florida State University
!
! Lagrange_3D_Class.f90 is part of the Spectral Element Libraries in Fortran (SELF).
! 
! Licensed under the Apache License, Version 2.0 (the "License"); 
! You may obtain a copy of the License at 
!
! http://www.apache.org/licenses/LICENSE-2.0 
!
! Unless required by applicable law or agreed to in writing, software 
! distributed under the License is distributed on an "AS IS" BASIS, 
! WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. 
! See the License for the specific language governing permissions and  
! limitations under the License.
!
! //////////////////////////////////////////////////////////////////////////////////////////////// !

 
 
MODULE Lagrange_3D_Class
! ========================================= Logs ================================================= !
!  2016/06/29 : Joe Schoonover (schoonover.numerics@gmail.com) : NEW VERSION (3.0)
!     > Switched over to Apache 2.0 License. 
!     > Cleaned up comments and code formatting to enhance readability.
!     > Removed Get/Set...AtNode routines to encourage larger memory accesses in practice.
!     > Renamed attributes of the data-structure, so that higher-dimension interpolants can follow
!       a simple naming pattern.
!     > Added interpolation and derivative matrices to the data-structure.
!     > Changed interpolation and differentiation routines to yield large matrix-vector and
!       matrix-matrix products respectively.
!     > This version (3.0) implements 3D interpolation as tensor products of polynomial interpolants
!       of the same degree. Reduced flexibility also reduces memory footprint and allows for more
!       straightforward optimizations that are expected to increase the arithmetic intensity of the
!       software.
!
! //////////////////////////////////////////////////////////////////////////////////////////////// ! 

!src/common
USE ModelPrecision
USE ConstantsDictionary
USE CommonRoutines
! src/interp
USE InterpolationSupportRoutines

IMPLICIT NONE



   ! The Lagrange_3D structure and routines are constructed based on the idea that the 3-D 
   ! interpolation can be conducted via a tensor product of three 1-D interpolants. It is assumed
   ! that the polynomial degree in each computational direction ("s", "p", and "q") is the same.
   ! For this reason, only a 1-D array of nodes, and a single set of barycentric weights needs
   ! to be stored. Derivatives and Interpolations are achieved via matrix-matrix and
   ! matrix-vector multiplication respectively. The SELF-TechnicalDocumentation provides the
   ! details on how this is achieved.
   TYPE, PUBLIC :: Lagrange_3D
      INTEGER                 :: nS      ! number of nodal points where we have obersvations
      INTEGER                 :: nSo     ! number of nodal points where we want observations
      INTEGER                 :: nIso    ! number of rows of the interpolation matrix
      INTEGER                 :: nIs     ! number of columns of the interpolation matrix
      INTEGER                 :: nD      ! number of columns needed to collapse two of the three
                                         ! dimensions of a 3-D array to one
      REAL(prec), ALLOCATABLE :: s(:)    ! locations where we have obervations
      REAL(prec), ALLOCATABLE :: bWs(:)  ! barycentric weights
      REAL(prec), ALLOCATABLE :: so(:)   ! Locations where we want observations
      REAL(prec), ALLOCATABLE :: Ts(:,:) ! Interpolation matrix to get us from what we have to what
                                         ! we want, in each computational direction.
      REAL(prec), ALLOCATABLE :: Ds(:,:) ! Derivative matrix to calculate the derivative of an 
                                         ! interpolant at the interpolation nodes. 
                                         ! This derivative matrix is used to calculate the 
                                         ! derivative in all computational direction. Such use
                                         ! requires the polynomial degree to be the same in each
                                         ! computational direction. 
      CONTAINS
      
      !-------------!
      ! Constructors/Destructors
      PROCEDURE :: Build => Build_Lagrange_3D
      PROCEDURE :: Trash => Trash_Lagrange_3D
  
      ! Accessors
      PROCEDURE :: SetNodes                  => SetNodes_Lagrange_3D
      PROCEDURE :: GetNodes                  => GetNodes_Lagrange_3D
      PROCEDURE :: SetAlternateNodes         => SetAlternateNodes_Lagrange_3D
      PROCEDURE :: GetAlternateNodes         => GetAlternateNodes_Lagrange_3D
      PROCEDURE :: SetWeights                => SetWeights_Lagrange_3D
      PROCEDURE :: GetWeights                => GetWeights_Lagrange_3D
      PROCEDURE :: SetNumberOfNodes          => SetNumberOfNodes_Lagrange_3D
      PROCEDURE :: GetNumberOfNodes          => GetNumberOfNodes_Lagrange_3D
      PROCEDURE :: SetNumberOfAlternateNodes => SetNumberOfAlternateNodes_Lagrange_3D
      PROCEDURE :: GetNumberOfAlternateNodes => GetNumberOfAlternateNodes_Lagrange_3D

      ! Type-Specific
      PROCEDURE :: CalculateBarycentricWeights  => CalculateBarycentricWeights_Lagrange_3D
      PROCEDURE :: CalculateInterpolationMatrix => CalculateInterpolationMatrix_Lagrange_3D
      PROCEDURE :: CalculateDerivativeMatrix    => CalculateDerivativeMatrix_Lagrange_3D
      PROCEDURE :: LagrangePolynomials          => LagrangePolynomials_Lagrange_3D
      PROCEDURE :: ApplyInterpolationMatrix     => ApplyInterpolationMatrix_Lagrange_3D
      PROCEDURE :: ApplyDerivativeMatrix        => ApplyDerivativeMatrix_Lagrange_3D
      
    END TYPE Lagrange_3D

 INTEGER, PRIVATE :: nDim = 2
 
 CONTAINS
!
!
!==================================================================================================!
!------------------------------- Manual Constructors/Destructors ----------------------------------!
!==================================================================================================!
!
!
 SUBROUTINE Build_Lagrange_3D( myPoly, nS, nSo, s, so )
 ! S/R Build_Lagrange_3D 
 !
 !   A manual constructor for the Lagrange_3D class.
 ! 
 !   Usage :
 !      CALL myPoly % Build( nS, nSo, s, so )
 !
 !      Allocates memory and fills in data for the attributes of the Lagrange_3D class. 
 !
 !   Input : 
 !      myPoly       The Lagrange_3D data structure
 !      nS           The number of observations where we have data
 !      nSo          The number of observations that we want
 !      s(0:nS)      The node locations where we have data
 !      so(0:nSo)    The node locations where we want data
 !
 !   Output :
 !      myPoly       The Lagrange_3D data structure with its attributes filled in
 !
 ! =============================================================================================== !
 ! DECLARATIONS 
   IMPLICIT NONE
   CLASS(Lagrange_3D), INTENT(inout) :: myPoly
   INTEGER, INTENT(in)               :: nS, nSo
   REAL(prec), INTENT(in)            :: s(0:nS), so(0:nSo)
   ! Local
   INTEGER :: nIso, nIs
   
      ! Set the number of observations (those we have and those we want)
      myPoly % nS   = nS
      myPoly % nSo  = nSo
      myPoly % nIso = nSo*( (nSo+2)*(nSo+1) + 1 )
      myPoly % nIs  = nS*( (nS+2)*(nS+1) + 1) ! The last column index of the interpolation matrix is 
                                              ! determined by "packing" three separate interpolation 
                                              ! matrices (that begin indices at 0) for each 
                                              ! computational direction into a single matrix
      myPoly % nD = nS*(nS+2)

      nIso = myPoly % nIso
      nIs  = myPoly % nIs
      
      ! Allocate storage
      ALLOCATE( myPoly % s(0:nS), myPoly % bWs(0:nS) )
      ALLOCATE( myPoly % so(0:nSo), myPoly % Ts(0:nISo,0:nS) )
      ALLOCATE( myPoly % Ds(0:nS,0:nS) )
      
      ! Fill in the nodal locations
      myPoly % s  = s
      myPoly % so = so

      ! and calculate the barycentric weights for quick interpolation.
      CALL myPoly % CalculateBarycentricWeights( )

      ! Using the two nodal locations, we can construct the interpolation matrix. The interpolation
      ! matrix enables quick interpolation.
      CALL myPoly % CalculateInterpolationMatrix( )

      CALL myPoly % CalculateDerivativeMatrix( )
 
 END SUBROUTINE Build_Lagrange_3D
!
!
!
SUBROUTINE Trash_Lagrange_3D(myPoly)
 ! S/R Trash_Lagrange_3D 
 !
 !   A manual destructor for the Lagrange_3D class.
 ! 
 !   Usage :
 !      CALL myPoly % Trash( )
 !
 !      Deallocates memory for the Lagrange_3D class. 
 !
 !   Input : 
 !      myPoly       The Lagrange_3D data structure
 !
 !   Output :
 !      myPoly       An empty Lagrange_3D data structure.
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS(Lagrange_3D), INTENT(inout) :: myPoly

      DEALLOCATE( myPoly % s, myPoly % bWs )
      DEALLOCATE( myPoly % so, myPoly % Ts )
      DEALLOCATE( myPoly % Ds )

 END SUBROUTINE Trash_Lagrange_3D
!
!
!==================================================================================================!
!--------------------------------------- Accessors ------------------------------------------------!
! These routines are meant to be for convenience, so users do not have to make reference directly to 
! the attributes of the data-structure directly. Of course, if the programmer knows the attributes,
! nothing is currently stopping them from accessing the data directly; no data-hiding is implemented
! here.
!==================================================================================================!
!
!
 SUBROUTINE SetNodes_Lagrange_3D( myPoly, sInput )
 ! S/R SetNodes_Lagrange_3D 
 !  
 !   Uses "sInput" to assign the attribute "s" of the Lagrange_3D data structure. 
 !
 !   Usage :
 !      CALL myPoly % SetNodes( sInput )
 !
 !   Input :
 !      sInput       A REAL array on nodal locations where we have observations
 !
 !   Output :
 !      myPoly       The Lagrange_3D structure with the "s" attribute updated
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS(Lagrange_3D), INTENT(inout) :: myPoly
   REAL(prec), INTENT(in)            :: sInput(0:myPoly % ns)

      myPoly % s = sInput

 END SUBROUTINE SetNodes_Lagrange_3D
!
!
!
 SUBROUTINE GetNodes_Lagrange_3D( myPoly, sOutput )
 ! S/R GetNodes_Lagrange_3D
 !  
 !   Uses the Lagrange_3D data structure to report the locations where we have data ("s") in the
 !   REAL array "sOutput". 
 !
 !   Usage :
 !      CALL myPoly % GetNodes( sOutput )
 !
 !   Input :
 !      myPoly       The Lagrange_3D structure (with "s" assigned)
 !
 !   Output :
 !      sOutput      A REAL array on nodal locations where we have observations
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS(Lagrange_3D), INTENT(in) :: myPoly
   REAL(prec), INTENT(out)        :: sOutput(0:myPoly % ns)

      sOutput = myPoly % s 

 END SUBROUTINE GetNodes_Lagrange_3D
!
!
!
SUBROUTINE SetAlternateNodes_Lagrange_3D( myPoly, sInput )
 ! S/R SetAlternateNodes_Lagrange_3D 
 !  
 !   Uses "sInput" to assign the attribute "so" of the Lagrange_3D data structure. Recall that "so"
 !   refers to the locations where we want to have new observations, ie, it is the set of locations
 !   that we will interpolate to.
 !
 !   Usage :
 !      CALL myPoly % SetNodes( sInput )
 !
 !   Input :
 !      sInput       A REAL array on nodal locations where we have observations
 !
 !   Output :
 !      myPoly       The Lagrange_3D structure with the "so" attribute updated
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS(Lagrange_3D), INTENT(inout) :: myPoly
   REAL(prec), INTENT(in)            :: sInput(0:myPoly % ns)

      myPoly % so = sInput

 END SUBROUTINE SetAlternateNodes_Lagrange_3D
!
!
!
 SUBROUTINE GetAlternateNodes_Lagrange_3D( myPoly, sOutput )
 ! S/R GetAlternateNodes_Lagrange_3D
 !  
 !   Uses the Lagrange_3D data structure to report the locations where we want data ("so") in the
 !   REAL array "sOutput". Recall that "so" refers to the locations where we want to have new 
 !   observations, ie, it is the set of locations that we will interpolate to.
 !
 !   Usage :
 !      CALL myPoly % GetNodes( sOutput )
 !
 !   Input :
 !      myPoly       The Lagrange_3D structure (with "so" assigned)
 !
 !   Output :
 !      sOutput      A REAL array on nodal locations where we have observations
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS(Lagrange_3D), INTENT(in) :: myPoly
   REAL(prec), INTENT(out)        :: sOutput(0:myPoly % ns)

      sOutput = myPoly % so 

 END SUBROUTINE GetAlternateNodes_Lagrange_3D
!
!
!
 SUBROUTINE SetWeights_Lagrange_3D( myPoly, wInput )
 ! S/R SetWeights_Lagrange_3D
 !  
 !   Uses "wInput" to assign the attribute "bWs" of the Lagrange_3D data structure. Recall that 
 !   "bWs" refers to barycentric interpolation weights.
 !
 !   Usage :
 !      CALL myPoly % SetWeights( wInput )
 !
 !   Input :
 !      wIn          A REAL array on nodal locations where we have observations
 !
 !   Output :
 !      myPoly       The Lagrange_3D structure with the "bWs" attribute updated
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS(Lagrange_3D), INTENT(inout) :: myPoly
   REAL(prec), INTENT(in)            :: wInput(0:myPoly % ns)

      myPoly % bWs = wInput

 END SUBROUTINE SetWeights_Lagrange_3D
!
!
!
 SUBROUTINE GetWeights_Lagrange_3D( myPoly, wOutput )
 ! S/R GetWeights_Lagrange_3D
 !  
 !   Uses the Lagrange_3D data structure to report the barycentric interpolation weights ("bWs") in 
 !   the REAL array "wOutput". Recall that "bWs" refers to barycentric interpolation weights.
 !
 !   Usage :
 !      CALL myPoly % GetWeights( wOutput )
 !
 !   Input :
 !      myPoly       The Lagrange_3D structure, with the "bWs" attribute assigned
 !
 !   Output :
 !      wOuput       A REAL array containing the barycentric weights
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS(Lagrange_3D), INTENT(in) :: myPoly
   REAL(prec), INTENT(out)        :: wOutput(0:myPoly % nS)

      wOutput = myPoly % bWs 

 END SUBROUTINE GetWeights_Lagrange_3D
!
!
!
 SUBROUTINE SetNumberOfNodes_Lagrange_3D( myPoly, N )
 ! S/R SetNumberOfNodes_Lagrange_3D
 !  
 !    Sets the "nS" attribute of the Lagrange_3D data structure to N. The "nS" attribute refers
 !    to the number of nodes where we have observations.
 !
 !    Usage :
 !       CALL myPoly % SetNumberOfNodes( N )
 !
 !    Input :
 !       myPoly      Lagrange_3D data structure
 !       N           Integer, indicating the number of nodes that we have
 !
 !    Output :
 !       myPoly      Lagrange_3D data structure with the number of nodes that we have (nS ) assigned
 !  
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS(Lagrange_3D), INTENT(inout) :: myPoly
   INTEGER, INTENT(in)               :: N
    
      myPoly % nS = N 

 END SUBROUTINE SetNumberOfNodes_Lagrange_3D
!
!
!
 SUBROUTINE GetNumberOfNodes_Lagrange_3D( myPoly, N )
 ! S/R GetNumberOfNodes_Lagrange_3D
 !  
 !    Gets the "nS" attribute from the Lagrange_3D data structure. The "nS" attribute refers
 !    to the number of nodes where we have observations.
 !
 !    Usage :
 !       CALL myPoly % GetNumberOfNodes( N )
 !
 !    Input :
 !       myPoly      Lagrange_3D data structure
 !
 !    Output :
 !       N           Integer, indicating the number of nodes that we have
 !  
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS(Lagrange_3D), INTENT(in) :: myPoly
   INTEGER, INTENT(out)           :: N

      N = myPoly % nS 

 END SUBROUTINE GetNumberOfNodes_Lagrange_3D
!
!
!
 SUBROUTINE SetNumberOfAlternateNodes_Lagrange_3D( myPoly, N )
 ! S/R SetNumberOfAlternateNodes_Lagrange_3D
 !  
 !    Sets the "nSo" attribute of the Lagrange_3D data structure to N. The "nSo" attribute refers
 !    to the number of AlternateNodess where we want observations.
 !
 !    Usage :
 !       CALL myPoly % SetNumberOfAlternateNodes( N )
 !
 !    Input :
 !       myPoly      Lagrange_3D data structure
 !       N           Integer, indicating the number of nodes that we want
 !
 !    Output :
 !       myPoly      Lagrange_3D data structure with the number of nodes that we want (nSo) assigned
 !  
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS(Lagrange_3D), INTENT(inout) :: myPoly
   INTEGER, INTENT(in)               :: N
    
      myPoly % nSo = N 

 END SUBROUTINE SetNumberOfAlternateNodes_Lagrange_3D
!
!
!
 SUBROUTINE GetNumberOfAlternateNodes_Lagrange_3D( myPoly, N )
 ! S/R GetNumberOfAlternateNodes_Lagrange_3D
 !  
 !    Gets the "nSo" attribute from the Lagrange_3D data structure. The "nSo" attribute refers
 !    to the number of nodes where we want observations.
 !
 !    Usage :
 !       CALL myPoly % GetNumberOfAlternateNodes( N )
 !
 !    Input :
 !       myPoly      Lagrange_3D data structure
 !
 !    Output :
 !       N           Integer, indicating the number of AlternateNodes that we have
 !  
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS(Lagrange_3D), INTENT(in) :: myPoly
   INTEGER, INTENT(out)           :: N

      N = myPoly % nso 

 END SUBROUTINE GetNumberOfAlternateNodes_Lagrange_3D
!
!
!==================================================================================================!
!--------------------------------- Type Specific Routines -----------------------------------------!
!==================================================================================================!
!
!
 SUBROUTINE CalculateBarycentricWeights_Lagrange_3D( myPoly )
 ! S/R CalculateBarycentricWeights_Lagrange_3D
 !  
 !    Calculates the barycentric weights from the interpolation nodes and stores them in 
 !    the "bWs" attribute. This routine should be called after the interpolation nodes (s) have 
 !    been assigned.
 !
 !    Usage :
 !       CALL myPoly % CalculateBarycentricWeights( )
 !
 !    Input :
 !       myPoly (% s)     The Lagrange_3D data structure is sent in with the interpolation nodes
 !                        already filled in.
 !
 !    Output :
 !       myPoly (% bWs)   The Lagrange_3D data structure is returned with the barycentric weights
 !                        filled in.
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS(Lagrange_3D), INTENT(inout) :: myPoly
   ! Local
   REAL(prec) :: s(0:myPoly % nS)
   INTEGER    :: N

      N = myPoly % nS
      s = myPoly % s

      myPoly % bWs = BarycentricWeights( s, N )
 
 END SUBROUTINE CalculateBarycentricWeights_Lagrange_3D
!
!
!
 SUBROUTINE CalculateInterpolationMatrix_Lagrange_3D( myPoly )  
 ! S/R CalculateInterpolationMatrix_Lagrange_3D 
 ! 
 !    The interpolation of a function from one set of points (myPoly % s) to another (myPoly % so)
 !    can be written as {f}_j = sum_{i}( f_i l_i(so_j) ). The sum (for each j) represents a matrix
 !    vector product where the factor "l_i(so_j)" is the interpolation matrix. This subroutine
 !    fills in the interpolation matrix.
 !
 !    Usage :
 !       CALL myPoly % CalculateInterpolationMatrix( )
 !
 !    Input :
 !       myPoly (% s) (% bWs) (% so)   Lagrange_3D data structure with the nodes, barycentric 
 !                                     weights, and alternate nodes filled in
 !
 !    Output
 !       myPoly (% Ts)                 Lagrange_3D data structure with the interpolation matrix
 !                                     filled in.
 ! 
 ! 
 ! =============================================================================================== !
   IMPLICIT NONE
   CLASS(Lagrange_3D), INTENT(inout) :: myPoly
   ! LOCAL
   REAL(prec) :: temp1, temp2
   REAL(prec) :: s(0:myPoly % nS), w(0:myPoly % nS)
   REAL(prec) :: so(0:myPoly % nSo), T(0:myPoly % nSo, 0:myPoly % nS)
   INTEGER    :: i, j, k, a, b, c, row, col,  N, Nnew
   LOGICAL    :: rowHasMatch 

      N    = myPoly % nS
      nNew = myPoly % nSo
      s    = myPoly % s
      w    = myPoly % bWs
      so   = myPoly % so

      T = InterpolationMatrix( s, w, so, N, nNew )

      DO c = 0, nNew
         DO b = 0, nNew
            DO a = 0, nNew

               row = a + (nNew+1)*( b + (nNew+1)*c )

               DO k = 0, N
                  DO j = 0, N
                     DO i = 0, N

                        col = i + (N+1)*( j + (N+1)*k )
                        myPoly % Ts(row,col) = T(a,i)*T(b,j)*T(c,k)

                     ENDDO
                  ENDDO
               ENDDO
            ENDDO
         ENDDO
      ENDDO

 END SUBROUTINE CalculateInterpolationMatrix_Lagrange_3D
!
!
!
 SUBROUTINE CalculateDerivativeMatrix_Lagrange_3D( myPoly )  
 ! S/R CalculateDerivativeMatrix_Lagrange_3D 
 !  
 !    Given nodal values of an interpolant, the derivative can be estimated at the interpolation 
 !    nodes using the summation
 !                             {f'}_j = sum_{i}( f_i l'_i(s_j) )
 !    The summation is a matrix-vector product, where the factor l'_i(s_j) is the derivative 
 !    matrix. This subroutine calculates the derivative matrix and stores it in the "Ds" attribute
 !    of myPoly for later use.
 ! 
 !     Usage :
 !        CALL myPoly % CalculateDerivativeMatrix( )
 !
 !     Input :
 !        myPoly (% s) (% bWs)    The Lagrange_3D data structure with the interpolation nodes ("s")
 !                                and barycentric weights ("bWs") filled in.
 !
 !     Output
 !        myPoly (% Ds)           The Lagrange_3D data structure with the derivative matrix filled
 !                                in.
 !
 ! =============================================================================================== !
   IMPLICIT NONE
   CLASS(Lagrange_3D), INTENT(inout) :: myPoly
   ! LOCAL
   REAL(prec) :: s(0:myPoly % nS), w(0:myPoly % nS)
   REAL(prec) :: dMat(0:myPoly % nS, 0:myPoly % nS)
   INTEGER    :: k, j, N

      N = myPoly % nS
      s = myPoly % s
      w = myPoly % bWs
      
      dMat = DerivativeMatrix( w, s, N )
  
      myPoly % Ds = dMat

 END SUBROUTINE CalculateDerivativeMatrix_Lagrange_3D
!
!
!
 FUNCTION LagrangePolynomials_Lagrange_3D( myPoly, sE ) RESULT( lAtS )  
 ! FUNCTION LagrangePolynomials_Lagrange_3D
 !  
 !    Given an evaluation location, this function returns each of the lagrange interpolating 
 !    polynomials evaluated at "sE". This is helpful if your program repeatedly requires 
 !    interpolation onto a given point.
 !
 !    This subroutine is from Alg. # on pg. # of D.A. Kopriva, 2011, "Implementing Spectral Element
 !    Methods for Scientists and Engineers"
 !
 !    Usage :
 !       lAtS = myPoly % LagrangePolynomials( sE )
 ! 
 !    Input :
 !       myPoly (% s) (% bWs)   Lagrange_3D Structure with the interpolation nodes and barycentric
 !                              weights filled in
 !   
 !       sE                     REAL(prec) value where we want to evaluate the Lagrange 
 !                              interpolating polynomials.
 !    Output :
 !       lAtS(0:myPoly % nS)    REAL(prec) array of the Lagrange interpolating polynomials evaluated
 !                              at sE
 !  
 ! =============================================================================================== !
   IMPLICIT NONE
   CLASS(Lagrange_3D) :: myPoly
   REAL(prec)         :: sE
   REAL(prec)         :: lAtS(0:myPoly % nS)
   ! LOCAL
   REAL(prec) :: temp1, temp2
   REAL(prec) :: s(0:myPoly % nS), w(0:myPoly % nS)
   INTEGER    :: j, N
   LOGICAL    :: xMatchesNode

      N = myPoly % nS
      s = myPoly % s
      w = myPoly % bWs

      xMatchesNode = .FALSE.

      DO j = 0, N
        
         lAtS(j) = ZERO

         IF( AlmostEqual(sE, s(j)) ) THEN
            lAtS(j) = ONE
            xMatchesNode = .TRUE.
         ENDIF 

     ENDDO

     IF( xMatchesNode )THEN 
        RETURN
     ENDIF

     temp1 = ZERO
     
     DO j = 0, N 
        temp2 = w(j)/(sE - s(j))
        lAtS(j) = temp2
        temp1 = temp1 + temp2
     ENDDO 
     
     DO j = 0, N 
        lAtS(j) = lAtS(j)/temp1
     ENDDO 
     

 END FUNCTION LagrangePolynomials_Lagrange_3D
!
!
!
 FUNCTION ApplyInterpolationMatrix_Lagrange_3D( myPoly, f ) RESULT( fNew )  
 ! FUNCTION ApplyInterpolationMatrix_Lagrange_3D
 !
 !   This function performs the matrix-vector multiply between the interpolation matrix 
 !   (myPoly % Ts) and the "vector" of nodal-values "f". The application of the interpolation 
 !   matrix maps "f" from the points "myPoly % s" to the points "myPoly % so".
 !
 !
 !   Usage :
 !      fNew = myPoly % ApplyInterpolationMatrix( f )
 !
 !   Input :
 !      myPoly (% Ts)               Lagrange_3D data structure that has the interpolation matrix 
 !                                  filled in.
 ! 
 !      f(:,:)                       A REAL(prec) 2-D array of nodal values of a function at the 
 !                                   points "myPoly % s".
 !
 !   Output :
 !      fNew(:,:)                   A REAL(prec) 2-D array of nodal values of a function at the 
 !                                  points "myPoly % so".
 !
 ! =============================================================================================== !
   IMPLICIT NONE
   CLASS(Lagrange_3D) :: myPoly
   REAL(prec)         :: f(0:myPoly % nS,0:myPoly % nS,0:myPoly % nS)
   REAL(prec)         :: fNew(0:myPoly % nSo,0:myPoly % nSo,0:myPoly % nSo)
   ! Local
   REAL(prec) :: fMapped(0:myPoly % nIs), fIntMapped(0:myPoly % nIso) 
   INTEGER    :: i, j, k, N, nNew, row
      
      N    = myPoly % nS
      nNew = myPoly % nSo
      
      DO k = 0, N
         DO j = 0, N
            DO i = 0, N
               row = i + (N+1)*( j + (N+1)*k )
               fMapped(row) = f(i,j,k)
            ENDDO
         ENDDO
      ENDDO

      ! Interpolation is achieved with a single matrix vector multiply
      fIntMapped = MATMUL( myPoly % Ts, f )

      DO k = 0, nNew
         DO j = 0, nNew
            DO i = 0, nNew
               row = i + (nNew+1)*( j + (nNew+1)*k )
               fNew(i,j,k) = fIntMapped(row)
            ENDDO
         ENDDO
      ENDDO

 END FUNCTION ApplyInterpolationMatrix_Lagrange_3D
!
!
!
 FUNCTION ApplyDerivativeMatrix_Lagrange_3D( myPoly, f ) RESULT( derF )  
 ! FUNCTION ApplyDerivativeMatrix_Lagrange_3D
 !
 !   This function performs the matrix-vector multiply between the Derivative matrix 
 !   (myPoly % Ds) and the "vector" of nodal-values "f". The application of the Derivative 
 !   matrix estimates the derivative of "f" at the nodal locations.
 !
 !   Usage :
 !      fNew = myPoly % ApplyDerivativeMatrix( f )
 !
 !   Input :
 !      myPoly (% Ds)               Lagrange_3D data structure that has the Derivative matrix 
 !                                  filled in.
 ! 
 !      f(0:myPoly % nS)            A REAL(prec) array of nodal values of a function at the points
 !                                  "myPoly % s".
 !
 !   Output :
 !      derF(:,:1:ndim)             A REAL(prec) array of nodal values of the derivative of the 
 !                                  interpolant at the tensor product of points "myPoly % s".
 !                                  derF(:,1) = dFds and derF(:,2) = dFdp. 
 !
 ! =============================================================================================== !
   IMPLICIT NONE
   CLASS(Lagrange_3D) :: myPoly
   REAL(prec)         :: f(0:myPoly % nS,0:myPoly % nS,0:myPoly % nS)
   REAL(prec)         :: derF(0:myPoly % nS,0:myPoly % nS,0:myPoly % nS, 1:ndim)
   ! LOCAL
   REAL(prec) :: floc1(0:myPoly % nS, 0:myPoly % nD)
   REAL(prec) :: floc2(0:myPoly % nS, 0:myPoly % nD)
   REAL(prec) :: floc3(0:myPoly % nS, 0:myPoly % nD)
   REAL(prec) :: dF1(0:myPoly % nS, 0:myPoly % nD)
   REAL(prec) :: dF2(0:myPoly % nS, 0:myPoly % nD)
   REAL(prec) :: dF3(0:myPoly % nS, 0:myPoly % nD)
   INTEGER    :: i, j, col

      N = myPoly % nS

      ! Condense the 3-D array to 2-D arrays to compute derivatives in each computational direction
      DO k = 0, N
         DO j = 0, N
            col = j + (N+1)*k
            floc1(0:N,col) = f(0:N,j,k)
            floc2(0:N,col) = f(j,0:N,k)
            floc3(0:N,col) = f(j,k,0:N)
         ENDDO
      ENDDO

      ! Matrix-Matrix multiply to compute derivatives
      dF1 = MATMUL( myPoly % Ds, floc1 )
      dF2 = MATMUL( myPoly % Ds, floc2 )
      dF3 = MATMUL( myPoly % Ds, floc3 )

      ! Remap derivatives to the 3-D arrays
      DO k = 0, N
         DO j = 0, N
            col = j + (N+1)*k
            derF(0:N,j,k,1) = dF1(0:N,col)
            derF(j,0:N,k,2) = dF2(0:N,col)
            derF(j,k,0:N,3) = dF3(0:N,col)
         ENDDO
      ENDDO

 END FUNCTION ApplyDerivativeMatrix_Lagrange_3D
!
!
!==================================================================================================!
!--------------------------------- File I/O Routines ----------------------------------------------!
!==================================================================================================!
!
!
END MODULE Lagrange_3D_Class
