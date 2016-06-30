! Quadrature.f90 (v 3.0)
! 
! Copyright 2015 Joseph Schoonover <schoonover.numerics@gmail.com>, The Florida State University
!
! Quadrature.f90 is part of the Spectral Element Libraries in Fortran (SELF).
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

 
MODULE Quadrature
! ========================================= Logs ================================================= !
!
! ================================================================================================ !
! src/common/
USE ModelPrecision
USE ConstantsDictionary


IMPLICIT NONE

 PUBLIC  :: ChebyshevQuadrature, LegendreQuadrature
 PRIVATE :: ChebyshevGauss, ChebyshevGaussLobatto, LegendreGauss, LegendreGaussLobatto, &
            LegendreGaussRadau, LegendreQandL, LegendrePolynomial

 CONTAINS
!
! ================================================================================================ !
! -------------------------------------- Chebyshev ----------------------------------------------- !
! ================================================================================================ !
!
 SUBROUTINE ChebyshevQuadrature( N, quadType, nodes, weights )
 ! S/R ChebyshevQuadrature
 !
 !   Usage :
 !      CALL ChebyshevQuadrature( N, quadType, nodes, weights )
 !
 !   Input :
 !
 !   Output :
 !   
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   INTEGER, INTENT(in)     :: N
   REAL(prec), INTENT(out) :: nodes(0:N)
   REAL(prec), INTENT(out) :: weights(0:N)
   INTEGER, INTENT(in)     :: QuadType
   
   
      ! Compute the quadrature nodes and weights
      IF( QuadType  == GAUSS_LOBATTO )then ! Gauss Lobatto quadrature
         CALL ChebyshevGaussLobatto( N, nodes, weights )
      ELSEIF( QuadType == GAUSS )then  ! Gauss Quadrature
         CALL ChebyshevGauss( N, nodes, weights )
      ELSE
         PRINT*,'Module Chebyshev.f90 : S/R GenerateChebyshevQuadrature: Invalid form. Stopping'
         STOP
      ENDIF

 END SUBROUTINE ChebyshevQuadrature
!
!
!
 SUBROUTINE ChebyshevGauss(N, nodes, weights)  
 ! S/R ChebyshevGauss
 !
 ! =============================================================================================== !
   IMPLICIT NONE
   INTEGER       :: N
   REAL(prec)    :: nodes(0:N)
   REAL(prec)    :: weights(0:N)
   ! LOCAL
   INTEGER    :: j ! Loop counter 
   REAL(prec) :: nReal, jReal, den


      nReal = REAL(N, prec)
      den = nReal + ONE

      DO j = 0, N
         jReal = REAL(j, prec)
         weights(j) = pi/den
         nodes(j) = -cos( 0.5_prec*(2.0_prec*jReal + ONE)*weights(j) )
      ENDDO

 END SUBROUTINE ChebyshevGauss
!
!
!
 SUBROUTINE ChebyshevGaussLobatto(N, nodes, weights)  
 ! S/R ChebyshevGaussLobatto
 !
 ! =============================================================================================== !
   IMPLICIT NONE
   INTEGER       :: N
   REAL(prec)    :: nodes(0:N)
   REAL(prec)    :: weights(0:N)
   ! LOCAL
   INTEGER    :: j ! Loop counter 
   REAL(prec) :: nReal, jReal

      nReal = REAL(N, prec)

      DO j = 0, N
         jReal = REAL(j, prec)
         weights(j) = pi/nReal
         nodes(j) = -cos( jReal*weights(j) )
      ENDDO

      weights(0)  = weights(0)*HALF
      weights(N) = weights(N)*HALF 

 END SUBROUTINE ChebyshevGaussLobatto
!
! ================================================================================================ !
! ---------------------------------------- Legendre ---------------------------------------------- !
! ================================================================================================ !
!
 SUBROUTINE LegendreQuadrature( N, nodes, weights, QuadType )
 ! S/R LegendreQuadrature
 !
 !   This subroutine serves as a wrapper for the quadrature generation routines.
 !   By specifying the quadrature type (QuadType) as Gauss or Gauss_Lobatto, this
 !   routine manages the appropriate call.
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   INTEGER, INTENT(in)     :: N
   REAL(prec), INTENT(out) :: nodes(0:N)
   REAL(prec), INTENT(out) :: weights(0:N)
   INTEGER, INTENT(in)     :: QuadType
   
   
      ! Compute the quadrature nodes and weights
      IF( QuadType  == GAUSS_LOBATTO )then ! Gauss Lobatto quadrature
         CALL LegendreGaussLobatto( N, nodes, weights )
      ELSEIF( QuadType == GAUSS )then  ! Gauss Quadrature
         CALL LegendreGauss( N, nodes, weights )
      ELSE
         PRINT*,'Module Legendre.f90 : S/R GenerateLegendreQuadrature: Invalid form. Stopping'
         STOP
      ENDIF

 END SUBROUTINE LegendreQuadrature
!
!
!
SUBROUTINE LegendreGauss(N, nodes, weights)  
 ! S/R LegendreGauss
 !  Descrition :
 !    Calculates the roots of the N+1-th Legendre Polynomial and the associated
 !    quadrature weights. Currently a Newton's method is employed with an
 !    initial guess of the roots as the Chebyshev roots. 
 !     * From Kopriva(2009), "Implementing Spectral Methods for Partial 
 !       Differential Equations", Ch 3, Section 2.1, p. 63; For large N, 
 !       "better starting values would be use the asymptotic representations 
 !        for the roots"
 !    
 !    Subroutine dependencies :
 !    S/R LegendrePolynomial, MODULE LEGEDNRE.f90
 !    
 !
 !  INut :
 !    INTEGER :: N - (N+1) is the degree of the legendre polynomial for 
 !                    which we calculate the root
 !
 !  Output :
 !    REAL(prec) :: nodes(0:N)- the roots of the N Legendre polynomials
 !  
 !    REAL(prec) :: weights(0:N) - the associated quadrature weights
 ! 
 !------------------------------------------------------------------------
  IMPLICIT NONE
  INTEGER :: N
  REAL(prec)    :: nodes(0:N)
  REAL(prec)    :: weights(0:N)
  ! LOCAL
  REAL(prec)    :: lN1, dlN1  ! Legendre polynomial and derivative
  REAL(prec)    :: delta
  INTEGER :: j, kIt ! Loop counter 
 

     IF( N == 0 ) then

        nodes(0) = ZERO
        weights(0) = TWO

     ELSEIF( N == 1 ) then

        nodes(0) = -sqrt(ONE/3.0_prec)
        weights(0) = ONE

        nodes(1) = -nodes(0)
        weights(1) = weights(0)

     ELSE ! use Newton's method

        DO j = 0, ( (N+1)/2 ) ! Loop over the roots

           nodes(j) = -cos( (TWO*REAL(j,prec) + ONE)*pi/(TWO*REAL(N,prec) + ONE) )


           DO kIt = 1, kItMax ! Loop over the Newton's iterations

              CALL LegendrePolynomial(N+1, nodes(j), lN1, dlN1)

              delta = -lN1/dlN1

              nodes(j) = nodes(j) + delta
 
              IF( abs(delta) <= TOL*nodes(j) ) EXIT

           ENDDO ! kIt, loop over the Newton's iterations


           CALL LegendrePolynomial(N+1, nodes(j), lN1, dlN1)

           weights(j) = TWO/( (ONE - nodes(j)*nodes(j))*dlN1*dlN1 )
 
           weights(N - j) = weights(j) ! uses symmetry to assign weights

           nodes(N - j) = -nodes(j)
           

        ENDDO ! j, loop over all of the roots

     ENDIF ! conditional on whether to use newton's method

     
     IF( mod(REAL(N,prec),TWO) == ZERO)then ! odd number of roots - get the weight for xRoot=0.0
         
        CALL LegendrePolynomial(N+1, ZERO, lN1, dlN1)
 
        nodes(N/2) = ZERO

        weights(N/2) = 2.0/(dlN1*dlN1)

     ENDIF ! IF we are looking for an odd number of roots.

 END SUBROUTINE LegendreGauss
!
!
!
 SUBROUTINE LegendreGaussLobatto(N, nodes, weights)  
 ! S/R LegendreGaussLobatto
 !  Descrition :
 !    Calculates the roots of the N+1-th Legendre Polynomial and the associated
 !    quadrature weights. Currently a Newton's method is employed with an
 !    initial guess of the roots as the Chebyshev roots. 
 !     * From Kopriva(2009), "Implementing Spectral Methods for Partial 
 !       Differential Equations", Ch 3, Section 2.1, p. 63; For large N, 
 !       "better starting values would be use the asymptotic representations 
 !        for the roots"
 !    
 !    Subroutine dependencies :
 !    S/R LegendrePolynomial, MODULE LEGEDNRE.f90
 !    
 !
 !  INut :
 !    INTEGER :: N - (N+1) is the degree of the legendre polynomial for 
 !                    which we calculate the root
 !
 !  Output :
 !    REAL(prec) :: nodes(0:N)- the roots of the N Legendre polynomials
 !  
 !    REAL(prec) :: weights(0:N) - the associated quadrature weights
 ! 
 !------------------------------------------------------------------------
  IMPLICIT NONE
  INTEGER :: N
  REAL(prec)    :: nodes(0:N)
  REAL(prec)    :: weights(0:N)
  ! LOCAL
  REAL(prec)    :: delta, q, qprime, lN
  INTEGER :: j, kIt ! Loop counter 
 

     IF( N == 1 ) then

        nodes(0) = -ONE
        weights(0) = ONE

        nodes(1) = ONE
        weights(1) = ONE

     ELSE ! use Newton's method

        nodes(0) = -ONE
        weights(0) = TWO/(REAL(N,prec)*(REAL(N,prec) + ONE) )

        nodes(N) = ONE
        weights(N) = weights(0)

        DO j = 1, ( (N+1)/2 -1 ) ! Loop over the roots

           nodes(j) = -cos( (REAL(j,prec) + 0.25_prec)*pi/REAL(N,prec) - &
                              3.0_prec/(8.0_prec*REAL(N,prec)*pi*(REAL(j,prec) + 0.25_prec) ) )


           DO kIt = 1, kItMax ! Loop over the Newton's iterations

              CALL LegendreQandL(N, nodes(j), q, qprime, lN)

              delta = -q/qprime

              nodes(j) = nodes(j) + delta
 
              IF( abs(delta) <= TOL*nodes(j) ) EXIT

           ENDDO ! kIt, loop over the Newton's iterations


           CALL LegendreQandL(N, nodes(j), q, qprime, lN)

           weights(j) = TWO/( REAL(N,prec)*(REAL(N,prec) + ONE)*lN*lN )
 
           weights(N - j) = weights(j) ! uses symmetry to assign weights

           nodes(N - j) = -nodes(j)
           

        ENDDO ! j, loop over all of the roots

     ENDIF ! conditional on whether to use newton's method

     
     IF( mod(REAL(N,prec),TWO) == ZERO)then ! odd number of roots - get the weight for xRoot=0.0
         
        CALL LegendreQandL(N, ZERO, q, qprime, lN)
 
        nodes(N/2) = ZERO

        weights(N/2) = TWO/( REAL(N,prec)*(REAL(N,prec) + ONE)*lN*lN )

     ENDIF ! IF we are looking for an odd number of roots.

 END SUBROUTINE LegendreGaussLobatto
!
!
!
 SUBROUTINE LegendreGaussRadau(N, nodes, weights)  
 ! S/R LegendreGaussRadau
 !  Descrition :
 !    Calculates the roots of the N+1-th Legendre Polynomial and the associated
 !    quadrature weights. Currently a Newton's method is employed with an
 !    initial guess of the roots as the Chebyshev roots. 
 !     * From Kopriva(2009), "Implementing Spectral Methods for Partial 
 !       Differential Equations", Ch 3, Section 2.1, p. 63; For large N, 
 !       "better starting values would be use the asymptotic representations 
 !        for the roots"
 !    
 !    Subroutine dependencies :
 !    S/R LegendrePolynomial, MODULE LEGEDNRE.f90
 !    
 !
 !  INut :
 !    INTEGER :: N - (N+1) is the degree of the legendre polynomial for 
 !                    which we calculate the root
 !
 !  Output :
 !    REAL(prec) :: nodes(0:N)- the roots of the N Legendre polynomials
 !  
 !    REAL(prec) :: weights(0:N) - the associated quadrature weights
 ! 
 !------------------------------------------------------------------------
  IMPLICIT NONE
  INTEGER :: N
  REAL(prec)    :: nodes(0:N)
  REAL(prec)    :: weights(0:N)
  ! LOCAL
  REAL(prec)    :: delta, dLdx, lN, q, qprime
  INTEGER :: j, kIt ! Loop counter 


     nodes(0) = ONE
     weights(0) = TWO/(REAL(N,prec)*(REAL(N,prec)) )


     DO j = 1, N ! Loop over the roots

        nodes(j) =  -cos( HALF*pi*(TWO*REAL(j+1,prec)-ONE)/(TWO*REAL(N,prec)-ONE) )

        
        DO kIt = 1, kItMax ! Loop over the Newton's iterations

           CALL LegendrePolynomial(N+1, nodes(j), lN, dLdx)

           q = lN
           qprime = dLdx*(ONE + nodes(j)) - lN

           CALL LegendrePolynomial(N, nodes(j), lN, dLdx)

           q = (q + lN)
           qprime = (qprime + (ONE + nodes(j))*dLdx - lN)/( ONE + nodes(j) )

           delta = -q/qprime 

           nodes(j) = nodes(j) + delta
 
           IF( abs(delta) <= TOL*nodes(j) ) EXIT

        ENDDO ! kIt, loop over the Newton's iterations

        CALL LegendrePolynomial(N, nodes(j), lN, dLdx)

        weights(j) = ONE/( (ONE-nodes(j))*dLdx*dLdx )

        nodes(j) = -nodes(j)
           

     ENDDO ! j, loop over all of the roots

 

 END SUBROUTINE LegendreGaussRadau
!
!
!
 SUBROUTINE LegendrePolynomial(N, x, lAtX, dLdxAtX)
 ! S/R LegendrePolynomial
 !  Descrition :
 !    Evaluates the N-th Legendre polynomial basis function
 !    at some x  ( in [-1,1] )
 !    
 !    Subroutine dependencies :
 !    (NONE)
 !
 !  INut :
 !    INTEGER :: N - the degree of the legendre polynomial
 !
 !    REAL(prec) :: x - the location where we wish to evaluate the basis function
 !
 !  Output :
 !    REAL(prec) :: lAtx- the value of the legendre polynomial at x
 !  
 !    REAL(prec) :: dLdxAtX - the value of the derivative of the legendre polynomial at x
 ! 
 !------------------------------------------------------------------------
 ! DECLARATIONS
   IMPLICIT NONE
   INTEGER       :: N
   REAL(prec)    :: x
   REAL(prec)    :: lAtX, dLdxAtX
   ! LOCAL
   REAL(prec) :: lNm1, lNm2, dlNm1, dlNm2
   INTEGER :: i

      IF( N == 0 )then
 
         lAtX = ONE    ! Legendre Polynomial
         dLdxAtX = ZERO ! Derivative
      
      ELSEIF( N == 1)then


         lAtX = x       ! Legendre Polynomial
         dLdxAtX = ONE  ! Derivative

      ELSE  ! Then we turn to the recursive relation for higher order Legendre Polynomials
  
         lnM2 = ONE     ! Initializing the recursive loop

         lnM1 = x

         dlnM2 = ZERO
   
         dlnM1 = ONE
        

         DO i = 2,N ! Recursive relation for the legendre polynomials
        
             lAtX = ((TWO*REAL(i,prec) - ONE)*x*lnM1 -&
                    (REAL(i,prec) - 1.0)*lnM2)/(REAL(i,prec))

             dldxAtX = dlnM2 + (TWO*REAL(i,prec)-ONE)*lnM1

             lnM2 = lnM1

             lnM1 = lAtX

             dlnM2 = dlnM1

             dlnM1 = dldxAtX
 
         ENDDO ! i, loop over the legendre polynomial degrees.

      ENDIF

 END SUBROUTINE LegendrePolynomial
!
!
!
SUBROUTINE LegendreQandL(N, x, q, qprime, lN)
 ! S/R LegendreQandL
 !  Descrition :
 !    
 !    Subroutine dependencies :
 !    (NONE)
 !
 !  INut :
 !
 !  Output :
 ! 
 !------------------------------------------------------------------------
 ! DECLARATIONS
   IMPLICIT NONE
   INTEGER       :: N
   REAL(prec)    :: x
   REAL(prec)    :: lN, q, qprime
   ! LOCAL
   REAL(prec) :: lNm1, lNm2, dlNm1, dlNm2, dlN, lN1, dlN1
   INTEGER :: i

      lNm2 = ONE
  
      lNm1 = x
  
      dlNm2 = ZERO

      dlNm1 = ONE

      DO i = 2, N

         lN = (TWO*i - ONE)/(REAL(i,prec))*x*lNm1 - (REAL(i,prec) - ONE)/(REAL(i,prec))*lNm2

         dlN = dlNm2 + (TWO*REAL(i,prec) - ONE)*lNm1

         lNm2 = lNm1

         lNm1 = lN

         dlNm2 = dlNm1

         dlNm1 = dlN

      ENDDO

      i = N + 1

      lN1 = (TWO*i - ONE)/(REAL(i,prec))*x*lN - (REAL(i,prec) - ONE)/(REAL(i,prec))*lNm2

      dlN1 = dlNm2 + (TWO*REAL(i,prec) - ONE)*lNm1

      q = lN1 - lNm2

      qprime = dlN1 - dlNm2
      

 END SUBROUTINE LegendreQandL

END MODULE Quadrature
