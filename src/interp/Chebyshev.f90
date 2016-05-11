! Chebyshev.f90
! 
! Copyright 2015 Joseph Schoonover <schoonover.numerics@gmail.com>
! 
! Chebyshev.f90 is part of the Spectral Element Libraries in Fortran (SELF).
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
 
 
MODULE Chebyshev
! ========================================= Logs ================================================= !
!2016-05-11  Joseph Schoonover  schoonover.numerics@gmail.com 
!
! //////////////////////////////////////////////////////////////////////////////////////////////// ! 
!
! This module contains auxiliary routines for generating the Chebyshev Gauss and Chebyshev Gauss-
! Lobatto integration nodes and quadrature.
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
