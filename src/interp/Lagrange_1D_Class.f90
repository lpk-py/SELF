! Lagrange_1D_Class.f90
! 
! Copyright 2015 Joseph Schoonover <schoonover.numerics@gmail.com>, The Florida State University
!
! Lagrange_1D_Class.f90 is part of the Spectral Element Libraries in Fortran (SELF).
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

 
 
MODULE Lagrange_1D_Class
! ========================================= Logs ================================================= !
!  2016/06/29 : Joe Schoonover (schoonover.numerics@gmail.com)
!     > Switched over to Apache 2.0 License. 
!     > Cleaned up comments and code formatting to enhance readability.
!     > Removed Get/Set...AtNode routines to encourage larger memory accesses in practice.
!     > Renamed attributes of the data-structure, so that higher-dimension interpolants can follow
!       a simple naming pattern.
!     > Added interpolation and derivative matrices to the data-structure.
!
! //////////////////////////////////////////////////////////////////////////////////////////////// ! 

!src/common
USE ModelPrecision
USE ConstantsDictionary
USE CommonRoutines
! src/interp
USE InterpolationSupportRoutines

IMPLICIT NONE



   TYPE, PUBLIC :: Lagrange_1D
      INTEGER                 :: nS      ! number of nodal points where we have obersvations
      INTEGER                 :: nSo     ! number of nodal points where we want observations
      REAL(prec), ALLOCATABLE :: s(:)    ! locations where we have obervations
      REAL(prec), ALLOCATABLE :: bWs(:)  ! barycentric weights
      REAL(prec), ALLOCATABLE :: so(:)   ! Locations where we want observations
      REAL(prec), ALLOCATABLE :: Ts(:,:) ! Interpolation matrix to get us from what we have to what
                                         ! we want
      REAL(prec), ALLOCATABLE :: Ds(:,:) ! Derivative matrix to calculate the derivative of an 
                                         ! interpolant at the interpolation nodes ("s").

      CONTAINS
      
      !-------------!
      ! Constructors/Destructors
      PROCEDURE :: Build => Build_Lagrange_1D
      PROCEDURE :: Trash => Trash_Lagrange_1D
  
      ! Accessors
      PROCEDURE :: SetNodes                  => SetNodes_Lagrange_1D
      PROCEDURE :: GetNodes                  => GetNodes_Lagrange_1D
      PROCEDURE :: SetAlternateNodes         => SetAlternateNodes_Lagrange_1D
      PROCEDURE :: GetAlternateNodes         => GetAlternateNodes_Lagrange_1D
      PROCEDURE :: SetWeights                => SetWeights_Lagrange_1D
      PROCEDURE :: GetWeights                => GetWeights_Lagrange_1D
      PROCEDURE :: SetNumberOfNodes          => SetNumberOfNodes_Lagrange_1D
      PROCEDURE :: GetNumberOfNodes          => GetNumberOfNodes_Lagrange_1D
      PROCEDURE :: SetNumberOfAlternateNodes => SetNumberOfAlternateNodes_Lagrange_1D
      PROCEDURE :: GetNumberOfAlternateNodes => GetNumberOfAlternateNodes_Lagrange_1D

      ! Type-Specific
      PROCEDURE :: CalculateBarycentricWeights  => CalculateBarycentricWeights_Lagrange_1D
      PROCEDURE :: CalculateInterpolationMatrix => CalculateInterpolationMatrix_Lagrange_1D
      PROCEDURE :: CalculateDerivativeMatrix    => CalculateDerivativeMatrix_Lagrange_1D
      PROCEDURE :: Interpolate                  => Interpolate_Lagrange_1D
      PROCEDURE :: LagrangePolynomials          => LagrangePolynomials_Lagrange_1D
      PROCEDURE :: Differentiate                => Differentiate_Lagrange_1D
      PROCEDURE :: ApplyInterpolationMatrix     => ApplyInterpolationMatrix_Lagrange_1D
      PROCEDURE :: ApplyDerivativeMatrix        => ApplyDerivativeMatrix_Lagrange_1D
      
    END TYPE Lagrange_1D

 
 CONTAINS
!
!
!==================================================================================================!
!------------------------------- Manual Constructors/Destructors ----------------------------------!
!==================================================================================================!
!
!
 SUBROUTINE Build_Lagrange_1D( myPoly, nS, nSo, s, so )
 ! S/R Build
 !
 !   A manual constructor for the Lagrange_1D class.
 ! 
 !   Usage :
 !      CALL myPoly % Build( nS, nSo, s, so )
 !
 !      Allocates memory and fills in data for the attributes of the Lagrange_1D class. 
 !
 !   Input : 
 !      myPoly       The Lagrange_1D data structure
 !      nS           The number of observations where we have data
 !      nSo          The number of observations that we want
 !      s(0:nS)      The node locations where we have data
 !      so(0:nSo)    The node locations where we want data
 !
 !   Output :
 !      myPoly       The Lagrange_1D data structure with its attributes filled in
 !
 ! =============================================================================================== !
 ! DECLARATIONS 
   IMPLICIT NONE
   CLASS(Lagrange_1D), INTENT(inout) :: myPoly
   INTEGER, INTENT(in)               :: nS, nSo
   REAL(prec), INTENT(in)            :: s(0:nS), so(0:nSo)
   
      ! Set the number of observations (those we have and those we want)
      myPoly % nS  = nS
      myPoly % nSo = nSo
      
      ! Allocate storage
      ALLOCATE( myPoly % s(0:nS), myPoly % bWs(0:nS) )
      ALLOCATE( myPoly % so(0:nSo), myPoly % Ts(0:nSo,0:nS) )
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
 
 END SUBROUTINE Build_Lagrange_1D
!
!
!
SUBROUTINE Trash_Lagrange_1D(myPoly)
 ! S/R Trash
 !
 !   A manual destructor for the Lagrange_1D class.
 ! 
 !   Usage :
 !      CALL myPoly % Trash( )
 !
 !      Deallocates memory for the Lagrange_1D class. 
 !
 !   Input : 
 !      myPoly       The Lagrange_1D data structure
 !
 !   Output :
 !      myPoly       An empty Lagrange_1D data structure.
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS(Lagrange_1D), INTENT(inout) :: myPoly

      DEALLOCATE( myPoly % s, myPoly % bWs )
      DEALLOCATE( myPoly % so, myPoly % Ts, myPoly % Ds )

 END SUBROUTINE Trash_Lagrange_1D
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
 SUBROUTINE SetNodes_Lagrange_1D( myPoly, sInput )
 ! S/R SetNodes
 !  
 !   Uses "sInput" to assign the attribute "s" of the Lagrange_1D data structure. 
 !
 !   Usage :
 !      CALL myPoly % SetNodes( sInput )
 !
 !   Input :
 !      sInput       A REAL array on nodal locations where we have observations
 !
 !   Output :
 !      myPoly       The Lagrange_1D structure with the "s" attribute updated
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS(Lagrange_1D), INTENT(inout) :: myPoly
   REAL(prec), INTENT(in)            :: sInput(0:myPoly % ns)

      myPoly % s = sInput

 END SUBROUTINE SetNodes_Lagrange_1D
!
!
!
 SUBROUTINE GetNodes_Lagrange_1D( myPoly, sOutput )
 ! S/R GetNodes_Lagrange_1D
 !  
 !   Uses the Lagrange_1D data structure to report the locations where we have data ("s") in the
 !   REAL array "sOutput". 
 !
 !   Usage :
 !      CALL myPoly % GetNodes( sOutput )
 !
 !   Input :
 !      myPoly       The Lagrange_1D structure (with "s" assigned)
 !
 !   Output :
 !      sOutput      A REAL array on nodal locations where we have observations
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS(Lagrange_1D), INTENT(in) :: myPoly
   REAL(prec), INTENT(out)        :: sOutput(0:myPoly % ns)

      sOutput = myPoly % s 

 END SUBROUTINE GetNodes_Lagrange_1D
!
!
!
SUBROUTINE SetAlternateNodes_Lagrange_1D( myPoly, sInput )
 ! S/R SetAlternateNodes
 !  
 !   Uses "sInput" to assign the attribute "so" of the Lagrange_1D data structure. Recall that "so"
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
 !      myPoly       The Lagrange_1D structure with the "so" attribute updated
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS(Lagrange_1D), INTENT(inout) :: myPoly
   REAL(prec), INTENT(in)            :: sInput(0:myPoly % ns)

      myPoly % so = sInput

 END SUBROUTINE SetAlternateNodes_Lagrange_1D
!
!
!
 SUBROUTINE GetAlternateNodes_Lagrange_1D( myPoly, sOutput )
 ! S/R GetAlternateNodes
 !  
 !   Uses the Lagrange_1D data structure to report the locations where we want data ("so") in the
 !   REAL array "sOutput". Recall that "so" refers to the locations where we want to have new 
 !   observations, ie, it is the set of locations that we will interpolate to.
 !
 !   Usage :
 !      CALL myPoly % GetNodes( sOutput )
 !
 !   Input :
 !      myPoly       The Lagrange_1D structure (with "so" assigned)
 !
 !   Output :
 !      sOutput      A REAL array on nodal locations where we have observations
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS(Lagrange_1D), INTENT(in) :: myPoly
   REAL(prec), INTENT(out)        :: sOutput(0:myPoly % ns)

      sOutput = myPoly % so 

 END SUBROUTINE GetAlternateNodes_Lagrange_1D
!
!
!
 SUBROUTINE SetWeights_Lagrange_1D( myPoly, wInput )
 ! S/R SetWeights
 !  
 !   Uses "wInput" to assign the attribute "bWs" of the Lagrange_1D data structure. Recall that 
 !   "bWs" refers to barycentric interpolation weights.
 !
 !   Usage :
 !      CALL myPoly % SetWeights( wInput )
 !
 !   Input :
 !      wIn          A REAL array on nodal locations where we have observations
 !
 !   Output :
 !      myPoly       The Lagrange_1D structure with the "bWs" attribute updated
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS(Lagrange_1D), INTENT(inout) :: myPoly
   REAL(prec), INTENT(in)            :: wInput(0:myPoly % ns)

      myPoly % bWs = wInput

 END SUBROUTINE SetWeights_Lagrange_1D
!
!
!
 SUBROUTINE GetWeights_Lagrange_1D( myPoly, wOutput )
 ! S/R GetWeights
 !  
 !   Uses the Lagrange_1D data structure to report the barycentric interpolation weights ("bWs") in 
 !   the REAL array "wOutput". Recall that "bWs" refers to barycentric interpolation weights.
 !
 !   Usage :
 !      CALL myPoly % GetWeights( wOutput )
 !
 !   Input :
 !      myPoly       The Lagrange_1D structure, with the "bWs" attribute assigned
 !
 !   Output :
 !      wOuput       A REAL array containing the barycentric weights
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS(Lagrange_1D), INTENT(in) :: myPoly
   REAL(prec), INTENT(out)        :: wOutput(0:myPoly % nS)

      wOutput = myPoly % bWs 

 END SUBROUTINE GetWeights_Lagrange_1D
!
!
!
 SUBROUTINE SetNumberOfNodes_Lagrange_1D( myPoly, N )
 ! S/R SetNumberOfNodes
 !  
 !    Sets the "nS" attribute of the Lagrange_1D data structure to N. The "nS" attribute refers
 !    to the number of nodes where we have observations.
 !
 !    Usage :
 !       CALL myPoly % SetNumberOfNodes( N )
 !
 !    Input :
 !       myPoly      Lagrange_1D data structure
 !       N           Integer, indicating the number of nodes that we have
 !
 !    Output :
 !       myPoly      Lagrange_1D data structure with the number of nodes that we have (nS ) assigned
 !  
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS(Lagrange_1D), INTENT(inout) :: myPoly
   INTEGER, INTENT(in)               :: N
    
      myPoly % nS = N 

 END SUBROUTINE SetNumberOfNodes_Lagrange_1D
!
!
!
 SUBROUTINE GetNumberOfNodes_Lagrange_1D( myPoly, N )
 ! S/R GetNumberOfNodes_Lagrange_1D
 !  
 !    Gets the "nS" attribute from the Lagrange_1D data structure. The "nS" attribute refers
 !    to the number of nodes where we have observations.
 !
 !    Usage :
 !       CALL myPoly % GetNumberOfNodes( N )
 !
 !    Input :
 !       myPoly      Lagrange_1D data structure
 !
 !    Output :
 !       N           Integer, indicating the number of nodes that we have
 !  
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS(Lagrange_1D), INTENT(in) :: myPoly
   INTEGER, INTENT(out)           :: N

      N = myPoly % nS 

 END SUBROUTINE GetNumberOfNodes_Lagrange_1D
!
!
!
 SUBROUTINE SetNumberOfAlternateNodes_Lagrange_1D( myPoly, N )
 ! S/R SetNumberOfAlternateNodes
 !  
 !    Sets the "nSo" attribute of the Lagrange_1D data structure to N. The "nSo" attribute refers
 !    to the number of AlternateNodess where we want observations.
 !
 !    Usage :
 !       CALL myPoly % SetNumberOfAlternateNodes( N )
 !
 !    Input :
 !       myPoly      Lagrange_1D data structure
 !       N           Integer, indicating the number of nodes that we want
 !
 !    Output :
 !       myPoly      Lagrange_1D data structure with the number of nodes that we want (nSo) assigned
 !  
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS(Lagrange_1D), INTENT(inout) :: myPoly
   INTEGER, INTENT(in)               :: N
    
      myPoly % nSo = N 

 END SUBROUTINE SetNumberOfAlternateNodes_Lagrange_1D
!
!
!
 SUBROUTINE GetNumberOfAlternateNodes_Lagrange_1D( myPoly, N )
 ! S/R GetNumberOfAlternateNodes
 !  
 !    Gets the "nSo" attribute from the Lagrange_1D data structure. The "nSo" attribute refers
 !    to the number of nodes where we want observations.
 !
 !    Usage :
 !       CALL myPoly % GetNumberOfAlternateNodes( N )
 !
 !    Input :
 !       myPoly      Lagrange_1D data structure
 !
 !    Output :
 !       N           Integer, indicating the number of AlternateNodess that we have
 !  
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS(Lagrange_1D), INTENT(in) :: myPoly
   INTEGER, INTENT(out)           :: N

      N = myPoly % nSo 

 END SUBROUTINE GetNumberOfAlternateNodes_Lagrange_1D
!
!
!==================================================================================================!
!--------------------------------- Type Specific Routines -----------------------------------------!
!==================================================================================================!
!
!
 SUBROUTINE CalculateBarycentricWeights_Lagrange_1D( myPoly )
 ! S/R CalculateBarycentricWeights
 !  
 !    Calculates the barycentric weights from the interpolation nodes and stores them in 
 !    the "bWs" attribute. This routine should be called after the interpolation nodes (s) have 
 !    been assigned.
 !
 !    This subroutine is from Alg. # on pg. # of D.A. Kopriva, 2011, "Implementing Spectral Element
 !    Methods for Scientists and Engineers"
 !
 !    Usage :
 !       CALL myPoly % CalculateBarycentricWeights( )
 !
 !    Input :
 !       myPoly (% s)     The Lagrange_1D data structure is sent in with the interpolation nodes
 !                        already filled in.
 !
 !    Output :
 !       myPoly (% bWs)   The Lagrange_1D data structure is returned with the barycentric weights
 !                        filled in.
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS(Lagrange_1D), INTENT(inout) :: myPoly
   ! Local
   REAL(prec) :: s(0:myPoly % nS)
   INTEGER    :: N

      N = myPoly % nS
      s = myPoly % s

      myPoly % bWs = BarycentricWeights( s, N )
 
 END SUBROUTINE CalculateBarycentricWeights_Lagrange_1D
!
!
!
 SUBROUTINE CalculateInterpolationMatrix_Lagrange_1D( myPoly )  
 ! S/R CalculateInterpolationMatrix 
 ! 
 !    The interpolation of a function from one set of points (myPoly % s) to another (myPoly % so)
 !    can be written as {f}_j = sum_{i}( f_i l_i(so_j) ). The sum (for each j) represents a matrix
 !    vector product where the factor "l_i(so_j)" is the interpolation matrix. This subroutine
 !    fills in the interpolation matrix.
 !
 !    This subroutine is from Alg. # on pg. # of D.A. Kopriva, 2011, "Implementing Spectral Element
 !    Methods for Scientists and Engineers"

 !    Usage :
 !       CALL myPoly % CalculateInterpolationMatrix( )
 !
 !    Input :
 !       myPoly (% s) (% bWs) (% so)   Lagrange_1D data structure with the nodes, barycentric 
 !                                     weights, and alternate nodes filled in
 !
 !    Output
 !       myPoly (% Ts)                 Lagrange_1D data structure with the interpolation matrix
 !                                     filled in.
 ! 
 ! =============================================================================================== !
   IMPLICIT NONE
   CLASS(Lagrange_1D), INTENT(inout) :: myPoly
   ! LOCAL
   REAL(prec) :: temp1, temp2
   REAL(prec) :: s(0:myPoly % nS), w(0:myPoly % nS)
   REAL(prec) :: so(0:myPoly % nSo), T(0:myPoly % nSo, 0:myPoly % nS)
   INTEGER    :: k, j,  N, Nnew
   LOGICAL    :: rowHasMatch 

      N    = myPoly % nS
      nNew = myPoly % nSo
      s    = myPoly % s
      w    = myPoly % bWs
      so   = myPoly % so

      T = InterpolationMatrix( s, w, so, N, nNew )

      myPoly % Ts = T

 END SUBROUTINE CalculateInterpolationMatrix_Lagrange_1D
!
!
!
 SUBROUTINE CalculateDerivativeMatrix_Lagrange_1D( myPoly )  
 ! S/R CalculateDerivativeMatrix 
 !  
 !    Given nodal values of an interpolant, the derivative can be estimated at the interpolation 
 !    nodes using the summation
 !                             {f'}_j = sum_{i}( f_i l'_i(s_j) )
 !    The summation is a matrix-vector product, where the factor l'_i(s_j) is the derivative 
 !    matrix. This subroutine calculates the derivative matrix and stores it in the "Ds" attribute
 !    of myPoly for later use.
 !
 !    This subroutine is from Alg. # on pg. # of D.A. Kopriva, 2011, "Implementing Spectral Element
 !    Methods for Scientists and Engineers"
 ! 
 !     Usage :
 !        CALL myPoly % CalculateDerivativeMatrix( )
 !
 !     Input :
 !        myPoly (% s) (% bWs)    The Lagrange_1D data structure with the interpolation nodes ("s")
 !                                and barycentric weights ("bWs") filled in.
 !
 !     Output
 !        myPoly (% Ds)           The Lagrange_1D data structure with the derivative matrix filled
 !                                in.
 !
 ! =============================================================================================== !
   IMPLICIT NONE
   CLASS(Lagrange_1D), INTENT(inout) :: myPoly
   ! LOCAL
   REAL(prec) :: s(0:myPoly % nS), w(0:myPoly % nS)
   REAL(prec) :: dMat(0:myPoly % nS, 0:myPoly % nS)
   INTEGER    :: k, j, N

      N = myPoly % nS
      s = myPoly % s
      w = myPoly % bWs
      
      dMat = DerivativeMatrix( w, s, N )

      myPoly % Ds = dMat

 END SUBROUTINE CalculateDerivativeMatrix_Lagrange_1D
!
!
!
 FUNCTION Interpolate_Lagrange_1D( myPoly, f, sE ) RESULT( interpF )  
 ! S/R Interpolate
 !  
 !    Interpolates the discrete array of data in f to the point sE using the Lagrange interpolating
 !    polynomials represented in myPoly.
 !    
 !    This subroutine is from Alg. # on pg. # of D.A. Kopriva, 2011, "Implementing Spectral Element
 !    Methods for Scientists and Engineers"
 !
 !    Usage :
 !       interpF = myPoly % Interpolate( f, sE )
 !
 !    Input :
 !       myPoly             Lagrange_1D data structure with the "s" and "bWs" filled in. "s"
 !                          corresponds to the nodal locations where "f" is specified.
 !
 !       f(0:myPoly % nS)   A single dimension REAL(prec) array of the function observations at
 !                          the locations specified by myPoly % s
 !
 !       sE                 A REAL(prec) value indicating where we want to interpolate the set of
 !                          data in "f".
 !
 !    Output :
 !       interpF            A REAL(prec) value that has the observations interpolated to sE.
 !
 ! =============================================================================================== !
   IMPLICIT NONE
   CLASS(Lagrange_1D) :: myPoly
   REAL(prec)         :: sE
   REAL(prec)         :: f(0:myPoly%ns)
   REAL(prec)         :: interpF
   ! LOCAL
   REAL(prec) :: num, den, t
   REAL(prec) :: w(0:myPoly % nS), s(0:myPoly % nS)
   INTEGER    :: j, N
 
     N = myPoly % nS
     w = myPoly % bWs
     s = myPoly % s
     
     num = ZERO
     den = ZERO

     DO j = 0, N

        IF( AlmostEqual(sE, s(j)) ) THEN 
          interpF = f(j)
          RETURN
        ELSE 
           
           t = w(j)/(sE - s(j))
           num = num + t*f(j)
           den = den + t

        ENDIF
        
     ENDDO

     interpF = num/den

 END FUNCTION Interpolate_Lagrange_1D
!
!
!
 FUNCTION Differentiate_Lagrange_1D(myPoly, f, sE ) RESULT(dInFdx)  
 ! FUNCTION Differentiate
 !  
 !    Evaluates the derivative of the interpolating polynomial described by the Lagrange_1D 
 !    data-structure "myPoly" and the nodal values "f" at the location "sE" 
 !
 !    This subroutine is from Alg. # on pg. # of D.A. Kopriva, 2011, "Implementing Spectral Element
 !    Methods for Scientists and Engineers"
 !
 !    Usage :
 !       interpF = myPoly % Differentiate( f, sE )
 !
 !    Input :
 !       myPoly             Lagrange_1D data structure with the "s" and "bWs" filled in. "s"
 !                          corresponds to the nodal locations where "f" is specified.
 !
 !       f(0:myPoly % nS)   A single dimension REAL(prec) array of the function observations at
 !                          the locations specified by myPoly % s
 !
 !       sE                 A REAL(prec) value indicating where we want to interpolate the set of
 !                          data in "f".
 !
 !    Output :
 !       dInFdx             A REAL(prec) value - the derivative of the interpolant at sE .
 !
 ! =============================================================================================== !
  IMPLICIT NONE
  CLASS(Lagrange_1D), INTENT(in) :: myPoly
  REAL(prec), INTENT(in)         :: sE
  REAL(prec), INTENT(in)         :: f(0:myPoly%ns)
  REAL(prec)                     :: dInFdx
  ! LOCAL
  REAL(prec) :: num, den, t, p
  REAL(prec) :: s(0:myPoly % nS), w(0:myPoly % nS)
  INTEGER    :: j, i, N 
  LOGICAL    :: atNode 

     N = myPoly % nS
     s = myPoly % s
     w = myPoly % bWs

     num = ZERO
     atNode = .FALSE.

     DO j = 0, N 

        IF( AlmostEqual(sE, s(j)) ) THEN 
           atNode = .TRUE.
           p = f(j)           
           den = -w(j)
           i = j
        ENDIF  
        
     ENDDO 

     IF( atNode ) THEN 

        DO j = 0, N 

           IF( .NOT.(j == i) ) THEN 
             num = num + w(j)*(p - f(j))/(sE - s(j))
           ENDIF

        ENDDO 

     ELSE !

        den = ZERO
        p = myPoly % Interpolate( f, sE ) 
        
        DO j = 0, N ! loop over the interpolation s
           t = w(j)/(sE - s(j))
           num = num + t*(p - f(j))/(sE - s(j))
           den = den + t
        ENDDO ! j, loop over the interpolation s

     ENDIF ! conditional, IF we're on an interpolating node

     dInFdx = num/den

 END FUNCTION Differentiate_Lagrange_1D
!
!
!
 FUNCTION LagrangePolynomials_Lagrange_1D( myPoly, sE ) RESULT( lAtS )  
 ! FUNCTION LagrangePolynomials
 !  
 !    Given an evaluation location, this function returns each of the lagrange interpolating 
 !    polynomials evaluated at "sE". This is helpful if your program repeatedly requires 
 !    interpolation onto a given point.
 !
 !    Usage :
 !       lAtS = myPoly % LagrangePolynomials( sE )
 ! 
 !    Input :
 !       myPoly (% s) (% bWs)   Lagrange_1D Structure with the interpolation nodes and barycentric
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
   CLASS(Lagrange_1D) :: myPoly
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
     

 END FUNCTION LagrangePolynomials_Lagrange_1D
!
!
!
 FUNCTION ApplyInterpolationMatrix_Lagrange_1D( myPoly, f ) RESULT( fNew )  
 ! ApplyInterpolationMatrix
 !
 !   This function performs the matrix-vector multiply between the interpolation matrix 
 !   (myPoly % Ts) and the "vector" of nodal-values "f". The application of the interpolation 
 !   matrix maps "f" from the points "myPoly % s" to the points "myPoly % so".
 !
 !   **This routine is provided as a template for higher-dimension interpolations. 
 !     It may be more beneficial (for 1-D problems) to simply call MATMUL in your main program.
 !
 !   Usage :
 !      fNew = myPoly % ApplyInterpolationMatrix( f )
 !
 !   Input :
 !      myPoly (% Ts)               Lagrange_1D data structure that has the interpolation matrix 
 !                                  filled in.
 ! 
 !      f(0:myPoly % nS)            A REAL(prec) array of nodal values of a function at the points
 !                                  "myPoly % s".
 !
 !   Output :
 !      fNew(0:myPoly % nSo)        A REAL(prec) array of nodal values of a function at the points
 !                                  "myPoly % so".
 !
 ! =============================================================================================== !
  IMPLICIT NONE
  CLASS(Lagrange_1D) :: myPoly
  REAL(prec)         :: f(0:myPoly % nS)
  REAL(prec)         :: fNew(0:myPoly % nSo)

     fNew = MATMUL( myPoly % Ts, f )

 END FUNCTION ApplyInterpolationMatrix_Lagrange_1D
!
!
!
 FUNCTION ApplyDerivativeMatrix_Lagrange_1D( myPoly, f ) RESULT( derF )  
 ! ApplyDerivativeMatrix
 !
 !   This function performs the matrix-vector multiply between the Derivative matrix 
 !   (myPoly % Ds) and the "vector" of nodal-values "f". The application of the Derivative 
 !   matrix estimates the derivative of "f" at the nodal locations.
 !
 !   **This routine is provided as a template for higher-dimensional derivatives. 
 !     It may be more beneficial (for 1-D problems) to simply call MATMUL in your main program.
 !
 !   Usage :
 !      fNew = myPoly % ApplyDerivativeMatrix( f )
 !
 !   Input :
 !      myPoly (% Ds)               Lagrange_1D data structure that has the Derivative matrix 
 !                                  filled in.
 ! 
 !      f(0:myPoly % nS)            A REAL(prec) array of nodal values of a function at the points
 !                                  "myPoly % s".
 !
 !   Output :
 !      derF(0:myPoly % nS)         A REAL(prec) array of nodal values of the derivative of the 
 !                                  interpolant at the points "myPoly % s".
 !
 ! =============================================================================================== !
  IMPLICIT NONE
  CLASS(Lagrange_1D) :: myPoly
  REAL(prec)         :: f(0:myPoly % nS)
  REAL(prec)         :: derF(0:myPoly % nS)

     derF = MATMUL( myPoly % Ds, f )

 END FUNCTION ApplyDerivativeMatrix_Lagrange_1D
!
!
!==================================================================================================!
!--------------------------------- File I/O Routines ----------------------------------------------!
!==================================================================================================!
!
!
 
!
!
!
END MODULE Lagrange_1D_Class
