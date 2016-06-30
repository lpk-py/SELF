! InterpolationSupportRoutines.f90 (v 3.0)
! 
! Copyright 2015 Joseph Schoonover <schoonover.numerics@gmail.com>, The Florida State University
!
! InterpolationSupportRoutines.f90 is part of the Spectral Element Libraries in Fortran (SELF).
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
 
 
MODULE InterpolationSupportRoutines

! src/common/
 USE ModelPrecision
 USE ConstantsDictionary
 USE CommonRoutines

 IMPLICIT NONE

 CONTAINS

 FUNCTION BarycentricWeights( s, N ) RESULT( w )
 !  BarycentricWeights
 !  
 !    Calculates the barycentric weights from the interpolation nodes s(0:N)
 !    The barycentric weights are the denominators of the lagrange interpolating polynomials
 !    that pass throught the points "s".
 !
 !    This subroutine is from Alg. # on pg. # of D.A. Kopriva, 2011, "Implementing Spectral Element
 !    Methods for Scientists and Engineers"
 !
 !    Usage :
 !       s = BarycentricWeights( s, N )
 ! 
 !    Input :
 !       s(0:N)       A REAL(prec) array of nodal locations where we have observations.
 !
 !       N            An integer specifying the polynomial degree.
 !
 !    Output :
 !       w(0:N)       A REAL(prec) array of the barycentric weights
 !    
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   INTEGER    :: N
   REAL(prec) :: s(0:N)
   REAL(prec) :: w(0:N)
   ! Local (to this function)
   INTEGER :: i, j
   REAL(prec) :: wlocal(0:N)
   
       wlocal(0:N) = ONE ! initialize the weights to 1

      ! Computes the product w_k = w_k*(x_k - x_j), k /= j
      DO j = 1,N ! loop over the interpolation nodes
         DO i = 0, j-1 ! loop to perform multiplication for weights

            wlocal(i) = wlocal(i)*( s(i) - s(j) )
            wlocal(j) = wlocal(j)*( s(j) - s(i) )

         ENDDO 
      ENDDO 
 
     ! Invert
      DO j = 0, N
        wlocal(j) = ONE/wlocal(j)
      ENDDO 
      w = wlocal
  
     RETURN

 END FUNCTION BarycentricWeights
!
!
!
 FUNCTION InterpolationMatrix( s, w, so, N, nNew ) RESULT( T )
 ! InterpolationMatrix 
 ! 
 !    This subroutine is from Alg. # on pg. # of D.A. Kopriva, 2011, "Implementing Spectral Element
 !    Methods for Scientists and Engineers"
 !
 ! 
 ! =============================================================================================== !
   IMPLICIT NONE
   INTEGER    :: N, nNew
   REAL(prec) :: s(0:N)
   REAL(prec) :: w(0:N)
   REAL(prec) :: so(0:nNew)
   REAL(prec) :: T(0:nNew,0:N)
   ! LOCAL
   REAL(prec) :: temp1, temp2
   INTEGER    :: row, col
   LOGICAL    :: rowHasMatch 


      DO row = 0, nNew ! loop over the new interpolation nodes ("so")

         rowHasMatch = .FALSE.
       
         DO col = 0, N ! loop over the old interpolation nodes ("s")

            T(row,col) = ZERO
           
            IF( AlmostEqual( so(row), s(col) ) )THEN
               rowHasMatch = .TRUE.
               T(row,col) = ONE
            ENDIF

         ENDDO 

         IF( .NOT.(rowHasMatch) )THEN 

            temp1 = ZERO

            DO col = 0, N ! loop over the old interpolation nodes ("s")         
               temp2 = w(col)/( so(row) - s(col) )
               T(row,col) = temp2
               temp1 = temp1 + temp2
            ENDDO 

            DO col = 0, N 
               T(row,col) = T(row,col)/temp1
            ENDDO

         ENDIF 

      ENDDO


 END FUNCTION InterpolationMatrix
!
!
!
 FUNCTION DerivativeMatrix( baryWeights, nodes, N ) RESULT( dMat )  
 !  DerivativeMatrix 
 !  
 !    This subroutine is from Alg. # on pg. # of D.A. Kopriva, 2011, "Implementing Spectral Element
 !    Methods for Scientists and Engineers"
 ! =============================================================================================== !
  IMPLICIT NONE
  INTEGER    :: N
  REAL(prec) :: baryWeights(0:N)
  REAL(prec) :: nodes(0:N)
  REAL(prec) :: dMat(0:N,0:N)
  ! LOCAL
  REAL(prec) :: x(0:N), w(0:N)
  INTEGER    :: row, col

      x = nodes
      w = baryweights

      DO row = 0, N ! loop over the interpolation nodes
         dMat(row,row) = ZERO
         DO col = 0, N ! loop over the interpolation nodes (again)
           
            IF( .NOT. (col == row) )THEN
               dMat(row,col) = w(col)/( w(row)*( x(row) - x(col) ) )
               dMat(row,row) = dMat(row,row) - dMat(row,col)
            ENDIF
        
         ENDDO 
      ENDDO 


 END FUNCTION DerivativeMatrix

END MODULE InterpolationSupportRoutines
