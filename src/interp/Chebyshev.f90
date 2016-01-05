MODULE Chebyshev
! Chebyshev.f90 (v2.1 - 12 Nov. 2015)
!
! schoonover.numerics@gmail.com
!
!
! 
! o  ( v2.0 - 17 July 2015 )
!
! This module contains auxiliary routines for generating the Chebyshev Gauss and Chebyshev Gauss-
! Lobatto integration nodes and quadrature.
!
! ( v2.1 - 12 Nov. 2015 )
! Additional comments were added and code format was switched to CaMeL for readability.
!    
! 
!
! 
! 
!
! ================================================================================================ !
USE ModelPrecision
USE ConstantsDictionary


IMPLICIT NONE


 CONTAINS

 !
 !=========================================================================!
 !-------------------- QUADRATURE WEIGHTS AND NODES -----------------------!
 !=========================================================================!
 !
 SUBROUTINE ChebyshevGauss(nP, xNodes, weights)  
 ! S/R ChebyshevGauss
 !  Description :
 !    
 !    
 !    Subroutine dependencies :
 !    
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
  INTEGER       :: nP
  REAL(prec)    :: xNodes(0:nP)
  REAL(prec)    :: weights(0:nP)
  ! LOCAL
  INTEGER    :: jX ! Loop counter 
  REAL(prec) :: nReal, jReal, den


     nReal = REAL(nP, prec)

     den = nReal + ONE

     DO jX = 0, nP

        jReal = REAL(jX, prec)

        weights(jX) = pi/den

        xNodes(jX) = -cos( 0.5_prec*(2.0_prec*jReal + ONE)*weights(jX) )

        
     ENDDO

 END SUBROUTINE ChebyshevGauss
!
!
!
 SUBROUTINE ChebyshevGaussLobatto(nP, xNodes, weights)  
 ! S/R ChebyshevGaussLobatto
 !  Description :
 !    
 !    
 !    Subroutine dependencies :
 !    
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
  INTEGER       :: nP
  REAL(prec)    :: xNodes(0:nP)
  REAL(prec)    :: weights(0:nP)
  ! LOCAL
  INTEGER    :: jX ! Loop counter 
  REAL(prec) :: nReal, jReal


     nReal = REAL(nP, prec)

     DO jX = 0, nP

        jReal = REAL(jX, prec)

        weights(jX) = pi/nReal

        xNodes(jX) = -cos( jReal*weights(jX) )
        
     ENDDO

     weights(0)  = weights(0)*HALF
     weights(nP) = weights(nP)*HALF 

 END SUBROUTINE ChebyshevGaussLobatto

END MODULE Chebyshev
