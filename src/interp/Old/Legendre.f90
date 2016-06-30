! Legendre.f90
! 
! Copyright 2015 Joseph Schoonover <schoonover.numerics@gmail.com>
! 
! Legendre.f90 is part of the Spectral Element Libraries in Fortran (SELF).
! 
! Permission is hereby granted, free of charge, to any person obtaining a copy of this software 
! and associated documentation files (the "Software"), to deal in the Software without restriction, 
! including without limitation the rights to use, copy, modify, merge, publish, distribute, 
! sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is 
! furnished to do so, subject to the following conditions: 
! 
! The above copyright notice and this permission notice shall be included in all copies or  
! substantial portions of the Software. 
! 
! THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING 
! BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND 
! NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, 
! DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, 
! OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE. 
!
! //////////////////////////////////////////////////////////////////////////////////////////////// !
 
 
MODULE Legendre
! ========================================= Logs ================================================= !
!2016-05-11  Joseph Schoonover  schoonover.numerics@gmail.com 
!
! //////////////////////////////////////////////////////////////////////////////////////////////// ! 

USE ModelPrecision
USE ModelFlags
USE ConstantsDictionary


IMPLICIT NONE


 CONTAINS
 !=========================================================================!
 !-------------------- POLYNOMIAL EVALUATION ROUTINE ----------------------!
 !=========================================================================!
 !
 SUBROUTINE LegendrePolynomial(nP, x, lAtX, dLdxAtX)
 ! S/R LegendrePolynomial
 !  Description :
 !    Evaluates the nP-th Legendre polynomial basis function
 !    at some x  ( in [-1,1] )
 !    
 !    Subroutine dependencies :
 !    (NONE)
 !
 !  Input :
 !    INTEGER :: nP - the degree of the legendre polynomial
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
   INTEGER       :: nP
   REAL(prec)    :: x
   REAL(prec)    :: lAtX, dLdxAtX
   ! LOCAL
   REAL(prec) :: lNm1, lNm2, dlNm1, dlNm2
   INTEGER :: iP

      IF( nP == 0 )then
 
         lAtX = ONE    ! Legendre Polynomial
         dLdxAtX = ZERO ! Derivative
      
      ELSEIF( nP == 1)then


         lAtX = x       ! Legendre Polynomial
         dLdxAtX = ONE  ! Derivative

      ELSE  ! Then we turn to the recursive relation for higher order Legendre Polynomials
  
         lnM2 = ONE     ! Initializing the recursive loop

         lnM1 = x

         dlnM2 = ZERO
   
         dlnM1 = ONE
        

         DO iP = 2,nP ! Recursive relation for the legendre polynomials
        
             lAtX = ((TWO*REAL(iP,prec) - ONE)*x*lnM1 -&
                    (REAL(iP,prec) - 1.0)*lnM2)/(REAL(iP,prec))

             dldxAtX = dlnM2 + (TWO*REAL(iP,prec)-ONE)*lnM1

             lnM2 = lnM1

             lnM1 = lAtX

             dlnM2 = dlnM1

             dlnM1 = dldxAtX
 
         ENDDO ! iP, loop over the legendre polynomial degrees.

      ENDIF

 END SUBROUTINE LegendrePolynomial
!
!
!
SUBROUTINE LegendreQandL(nP, x, q, qprime, lN)
 ! S/R LegendreQandL
 !  Description :
 !    
 !    Subroutine dependencies :
 !    (NONE)
 !
 !  Input :
 !
 !  Output :
 ! 
 !------------------------------------------------------------------------
 ! DECLARATIONS
   IMPLICIT NONE
   INTEGER       :: nP
   REAL(prec)    :: x
   REAL(prec)    :: lN, q, qprime
   ! LOCAL
   REAL(prec) :: lNm1, lNm2, dlNm1, dlNm2, dlN, lNp1, dlNp1
   INTEGER :: iP

      lNm2 = ONE
  
      lNm1 = x
  
      dlNm2 = ZERO

      dlNm1 = ONE

      DO iP = 2, nP

         lN = (TWO*iP - ONE)/(REAL(iP,prec))*x*lNm1 - (REAL(iP,prec) - ONE)/(REAL(iP,prec))*lNm2

         dlN = dlNm2 + (TWO*REAL(iP,prec) - ONE)*lNm1

         lNm2 = lNm1

         lNm1 = lN

         dlNm2 = dlNm1

         dlNm1 = dlN

      ENDDO

      iP = nP + 1

      lNp1 = (TWO*iP - ONE)/(REAL(iP,prec))*x*lN - (REAL(iP,prec) - ONE)/(REAL(iP,prec))*lNm2

      dlNp1 = dlNm2 + (TWO*REAL(iP,prec) - ONE)*lNm1

      q = lNp1 - lNm2

      qprime = dlNp1 - dlNm2
      

 END SUBROUTINE LegendreQandL
 !
 !=========================================================================!
 !-------------------- QUADRATURE WEIGHTS AND NODES -----------------------!
 !=========================================================================!
 !
 SUBROUTINE GenerateLegendreQuadrature( nP, xNodes, weights, QuadType )
 !
 !   This subroutine serves as a wrapper for the quadrature generation routines.
 !   By specifying the quadrature type (QuadType) as Gauss or Gauss_Lobatto, this
 !   routine manages the appropriate call.
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   INTEGER, INTENT(in)     :: nP
   REAL(prec), INTENT(out) :: xNodes(0:nP)
   REAL(prec), INTENT(out) :: weights(0:nP)
   INTEGER, INTENT(in)     :: QuadType
   
   
      ! Compute the quadrature nodes and weights
      IF( QuadType  == GAUSS_LOBATTO )then ! Gauss Lobatto quadrature

         CALL LegendreGaussLobatto( nP, xNodes, weights )
   
      ELSEIF( QuadType == GAUSS )then  ! Gauss Quadrature

         CALL LegendreGauss( nP, xNodes, weights )

      ELSE
      
         PRINT*,'Module Legendre.f90 : S/R GenerateLegendreQuadrature: Invalid form. Stopping'
         STOP

      ENDIF


 END SUBROUTINE GenerateLegendreQuadrature
!
!
!
 SUBROUTINE LegendreGauss(nP, xNodes, weights)  
 ! S/R LegendreGauss
 !  Description :
 !    Calculates the roots of the nP+1-th Legendre Polynomial and the associated
 !    quadrature weights. Currently a Newton's method is employed with an
 !    initial guess of the roots as the Chebyshev roots. 
 !     * From Kopriva(2009), "Implementing Spectral Methods for Partial 
 !       Differential Equations", Ch 3, Section 2.1, p. 63; For large nP, 
 !       "better starting values would be use the asymptotic representations 
 !        for the roots"
 !    
 !    Subroutine dependencies :
 !    S/R LegendrePolynomial, MODULE LEGEDNRE.f90
 !    
 !
 !  Input :
 !    INTEGER :: nP - (nP+1) is the degree of the legendre polynomial for 
 !                    which we calculate the root
 !
 !  Output :
 !    REAL(prec) :: xNodes(0:nP)- the roots of the nP Legendre polynomials
 !  
 !    REAL(prec) :: weights(0:nP) - the associated quadrature weights
 ! 
 !------------------------------------------------------------------------
  IMPLICIT NONE
  INTEGER :: nP
  REAL(prec)    :: xNodes(0:nP)
  REAL(prec)    :: weights(0:nP)
  ! LOCAL
  REAL(prec)    :: lNp1, dlNp1  ! Legendre polynomial and derivative
  REAL(prec)    :: delta
  INTEGER :: jX, kIt ! Loop counter 
 

     IF( nP == 0 ) then

        xNodes(0) = ZERO
        weights(0) = TWO

     ELSEIF( nP == 1 ) then

        xNodes(0) = -sqrt(ONE/3.0_prec)
        weights(0) = ONE

        xNodes(1) = -xNodes(0)
        weights(1) = weights(0)

     ELSE ! use Newton's method

        DO jX = 0, ( (nP+1)/2 ) ! Loop over the roots

           xNodes(jX) = -cos( (TWO*REAL(jX,prec) + ONE)*pi/(TWO*REAL(nP,prec) + ONE) )


           DO kIt = 1, kItMax ! Loop over the Newton's iterations

              CALL LegendrePolynomial(nP+1, xNodes(jX), lNp1, dlNp1)

              delta = -lNp1/dlNp1

              xNodes(jX) = xNodes(jX) + delta
 
              IF( abs(delta) <= TOL*xNodes(jX) ) EXIT

           ENDDO ! kIt, loop over the Newton's iterations


           CALL LegendrePolynomial(nP+1, xNodes(jX), lNp1, dlNp1)

           weights(jX) = TWO/( (ONE - xNodes(jX)*xNodes(jX))*dlNp1*dlNp1 )
 
           weights(nP - jX) = weights(jX) ! uses symmetry to assign weights

           xNodes(nP - jX) = -xNodes(jX)
           

        ENDDO ! jX, loop over all of the roots

     ENDIF ! conditional on whether to use newton's method

     
     IF( mod(REAL(nP,prec),TWO) == ZERO)then ! odd number of roots - get the weight for xRoot=0.0
         
        CALL LegendrePolynomial(nP+1, ZERO, lNp1, dlNp1)
 
        xNodes(nP/2) = ZERO

        weights(nP/2) = 2.0/(dlNp1*dlNp1)

     ENDIF ! IF we are looking for an odd number of roots.

 END SUBROUTINE LegendreGauss
!
!
!
 SUBROUTINE LegendreGaussLobatto(nP, xNodes, weights)  
 ! S/R LegendreGaussLobatto
 !  Description :
 !    Calculates the roots of the nP+1-th Legendre Polynomial and the associated
 !    quadrature weights. Currently a Newton's method is employed with an
 !    initial guess of the roots as the Chebyshev roots. 
 !     * From Kopriva(2009), "Implementing Spectral Methods for Partial 
 !       Differential Equations", Ch 3, Section 2.1, p. 63; For large nP, 
 !       "better starting values would be use the asymptotic representations 
 !        for the roots"
 !    
 !    Subroutine dependencies :
 !    S/R LegendrePolynomial, MODULE LEGEDNRE.f90
 !    
 !
 !  Input :
 !    INTEGER :: nP - (nP+1) is the degree of the legendre polynomial for 
 !                    which we calculate the root
 !
 !  Output :
 !    REAL(prec) :: xNodes(0:nP)- the roots of the nP Legendre polynomials
 !  
 !    REAL(prec) :: weights(0:nP) - the associated quadrature weights
 ! 
 !------------------------------------------------------------------------
  IMPLICIT NONE
  INTEGER :: nP
  REAL(prec)    :: xNodes(0:nP)
  REAL(prec)    :: weights(0:nP)
  ! LOCAL
  REAL(prec)    :: delta, q, qprime, lN
  INTEGER :: jX, kIt ! Loop counter 
 

     IF( nP == 1 ) then

        xNodes(0) = -ONE
        weights(0) = ONE

        xNodes(1) = ONE
        weights(1) = ONE

     ELSE ! use Newton's method

        xNodes(0) = -ONE
        weights(0) = TWO/(REAL(nP,prec)*(REAL(nP,prec) + ONE) )

        xNodes(nP) = ONE
        weights(nP) = weights(0)

        DO jX = 1, ( (nP+1)/2 -1 ) ! Loop over the roots

           xNodes(jX) = -cos( (REAL(jX,prec) + 0.25_prec)*pi/REAL(nP,prec) - &
                              3.0_prec/(8.0_prec*REAL(nP,prec)*pi*(REAL(jX,prec) + 0.25_prec) ) )


           DO kIt = 1, kItMax ! Loop over the Newton's iterations

              CALL LegendreQandL(nP, xNodes(jX), q, qprime, lN)

              delta = -q/qprime

              xNodes(jX) = xNodes(jX) + delta
 
              IF( abs(delta) <= TOL*xNodes(jX) ) EXIT

           ENDDO ! kIt, loop over the Newton's iterations


           CALL LegendreQandL(nP, xNodes(jX), q, qprime, lN)

           weights(jX) = TWO/( REAL(nP,prec)*(REAL(nP,prec) + ONE)*lN*lN )
 
           weights(nP - jX) = weights(jX) ! uses symmetry to assign weights

           xNodes(nP - jX) = -xNodes(jX)
           

        ENDDO ! jX, loop over all of the roots

     ENDIF ! conditional on whether to use newton's method

     
     IF( mod(REAL(nP,prec),TWO) == ZERO)then ! odd number of roots - get the weight for xRoot=0.0
         
        CALL LegendreQandL(nP, ZERO, q, qprime, lN)
 
        xNodes(nP/2) = ZERO

        weights(nP/2) = TWO/( REAL(nP,prec)*(REAL(nP,prec) + ONE)*lN*lN )

     ENDIF ! IF we are looking for an odd number of roots.

 END SUBROUTINE LegendreGaussLobatto
!
!
!
 SUBROUTINE LegendreGaussRadau(nP, xNodes, weights)  
 ! S/R LegendreGaussRadau
 !  Description :
 !    Calculates the roots of the nP+1-th Legendre Polynomial and the associated
 !    quadrature weights. Currently a Newton's method is employed with an
 !    initial guess of the roots as the Chebyshev roots. 
 !     * From Kopriva(2009), "Implementing Spectral Methods for Partial 
 !       Differential Equations", Ch 3, Section 2.1, p. 63; For large nP, 
 !       "better starting values would be use the asymptotic representations 
 !        for the roots"
 !    
 !    Subroutine dependencies :
 !    S/R LegendrePolynomial, MODULE LEGEDNRE.f90
 !    
 !
 !  Input :
 !    INTEGER :: nP - (nP+1) is the degree of the legendre polynomial for 
 !                    which we calculate the root
 !
 !  Output :
 !    REAL(prec) :: xNodes(0:nP)- the roots of the nP Legendre polynomials
 !  
 !    REAL(prec) :: weights(0:nP) - the associated quadrature weights
 ! 
 !------------------------------------------------------------------------
  IMPLICIT NONE
  INTEGER :: nP
  REAL(prec)    :: xNodes(0:nP)
  REAL(prec)    :: weights(0:nP)
  ! LOCAL
  REAL(prec)    :: delta, dLdx, lN, q, qprime
  INTEGER :: jX, kIt ! Loop counter 


     xNodes(0) = ONE
     weights(0) = TWO/(REAL(nP,prec)*(REAL(nP,prec)) )


     DO jX = 1, nP ! Loop over the roots

        xNodes(jX) =  -cos( HALF*pi*(TWO*REAL(jX+1,prec)-ONE)/(TWO*REAL(nP,prec)-ONE) )

        
        DO kIt = 1, kItMax ! Loop over the Newton's iterations

           CALL LegendrePolynomial(nP+1, xNodes(jX), lN, dLdx)

           q = lN
           qprime = dLdx*(ONE + xNodes(jX)) - lN

           CALL LegendrePolynomial(nP, xNodes(jX), lN, dLdx)

           q = (q + lN)
           qprime = (qprime + (ONE + xNodes(jX))*dLdx - lN)/( ONE + xNodes(jX) )

           delta = -q/qprime 

           xNodes(jX) = xNodes(jX) + delta
 
           IF( abs(delta) <= TOL*xNodes(jX) ) EXIT

        ENDDO ! kIt, loop over the Newton's iterations

        CALL LegendrePolynomial(nP, xNodes(jX), lN, dLdx)

        weights(jX) = ONE/( (ONE-xNodes(jX))*dLdx*dLdx )

        xNodes(jX) = -xNodes(jX)
           

     ENDDO ! jX, loop over all of the roots

 

 END SUBROUTINE LegendreGaussRadau


END MODULE Legendre
