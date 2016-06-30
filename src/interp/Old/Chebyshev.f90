! Chebyshev.f90
! 
! Copyright 2015 Joseph Schoonover <schoonover.numerics@gmail.com>, The Florida State University
!
! Chebyshev.f90 is part of the Spectral Element Libraries in Fortran (SELF).
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

 
 
MODULE Chebyshev
! ========================================= Logs ================================================= !
!
! ================================================================================================ !
! src/common/
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
