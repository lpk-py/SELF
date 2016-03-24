! ModalCutoffFilter3D_Class.f90 ( new with v2.1 - 21 March 2016)
! 
! ====================================== LICENSE ================================================= !
!
!  This program is free software; you can redistribute it and/or modify
!  it under the terms of the GNU General Public License as published by
!  the Free Software Foundation; either version 2 of the License, or
!  (at your option) any later version.
!  
!  This program is distributed in the hope that it will be useful,
!  but WITHOUT ANY WARRANTY; without even the implied warranty of
!  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!  GNU General Public License for more details.
!  
!  You should have received a copy of the GNU General Public License
!  along with this program; if not, write to the Free Software
!  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
!  MA 02110-1301, USA.
! 
! ==================================== Module History ============================================ ! 
! 
! o  (ver 2.1) March 2016
!   This module provides a data structure for constructing and storing element-local filter matrices
!   that will be used for polynomial de-aliasing or as part of SGS parameterizations. 3-D arrays,
!   representative of the nodal values of an interpolant, can be passed in to the supplied routine
!   "ApplyFilter" to return a filtered form of the interpolant.  
!    In this version, a simple "cutoff" filter is used, and only the high-polynomial 
!    cutoff degrees are required. Only the Legendre basis is considered here.
!
!    This module was inspired by the paper
!  
!     D. Flad, A. Beck, and C. Munz, (2016) "Simulation of underresolved turbulent flows by adaptive 
!        filtering using the high order discontinuous Galerkin spectral element method", JCP, 313, 
!        pp. 1-12
!
!
!    The data-structure provided here can be used with a high-end routine in order to incorporate
!    "switch-based" filtering of the prognostic solution. For under-resolved simulation, such 
!    filtering can provide fine-tuned control of dissipation needed to remain stable.
!
! ========================================= Logs ================================================= !
! yyyy-mm-dd  Joe  <joe@clay>
!
!   
! 
! //////////////////////////////////////////////////////////////////////////////////////////////// !



MODULE ModalCutoffFilter3D_Class

USE ModelPrecision
USE ModelFlags
USE ConstantsDictionary
USE CommonRoutines

USE Legendre
USE Lagrange_3D_Class

USE NodalStorage_3D_Class


   TYPE ModalCutoffFilter3D
      INTEGER                 :: nS, nP, nQ, nDOF
      INTEGER                 :: nSco, nPco, nQco ! modal cutoff numbers
      REAL(prec), ALLOCATABLE :: filterMat(:,:)
      
      CONTAINS
  
      PROCEDURE :: Build => Build_ModalCutoffFilter3D
      PROCEDURE :: Trash => Trash_ModalCutoffFilter3D
      ! 
      PROCEDURE :: SetNumberOfNodes => SetNumberOfNodes_ModalCutoffFilter3D
      PROCEDURE :: GetNumberOfNodes => GetNumberOfNodes_ModalCutoffFilter3D


   END TYPE ModalCutoffFilter3D


CONTAINS
!
!
!==================================================================================================!
!------------------------------- Manual Constructors/Destructors ----------------------------------!
!==================================================================================================!
!
!
 SUBROUTINE Build_ModalCutoffFilter3D( thisFilter, thisStorage, sCutoff, pCutoff, qCutoff )
 ! S/R Build
 !
 ! This constructor allocates storage for the SE filter. The filter matrix is constructed by 
 ! multiplying three matrices together. The first matrix maps the nodal values of the interpolant
 ! to modal form (using the Legendre basis). The second matrix applies the modal truncation
 ! operator. The third matrix computes the nodal values associated with the filtered coefficients;
 ! the third matrix is the inverse of the first.
 ! 
 ! The filter matrix preserves the energy of the lower order modes (below the cutoff mode), and
 ! effectively sets the higher order modal energy to zero. 
 ! 
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS(ModalCutoffFilter3D), INTENT(inout) :: thisFilter
   TYPE(NodalStorage_3D), INTENT(in)             :: thisStorage
   INTEGER, INTENT(in)                           :: sCutoff, pCutoff, qCutoff
   ! Local
   REAL(prec) :: s(0:thisStorage % nS), p(0:thisStorage % nP), q(0:thisStorage % nQ)
   REAL(prec) :: ws(0:thisStorage % nS), wp(0:thisStorage % nP), wq(0:thisStorage % nQ)
   REAL(prec) :: Lnorm, Li, Lj, Lk, Ls
   REAL(prec), ALLOCATABLE :: Pfilt(:,:), V(:,:), VInv(:,:)
   INTEGER    :: nS, nP, nQ, i, j, k, a, b, c, row, col 
   
      CALL thisStorage % GetNumberOfNodes( nS, nP, nQ )
      CALL thisFilter % SetNumberOfNodes( nS, nP, nQ )
      nDOF = (nS+1)*(nP+1)*(nQ+1)
      CALL thisFilter % SetNDOF( nDOF )
      CALL thisFilter % SetCutoffs( sCutoff, pCutoff, qCutoff )

      ALLOCATE( thisFilter % filterMat(1:nDOF,1:nDOF) )

      ALLOCATE( Pfilt(1:nDOF,1:nDOF), V(1:nDOF,1:nDOF), VInv(1:nDOF,1:nDOF) )

      thisFilter % filterMat = ZERO

      CALL thisStorage % interp % GetNodes( s, p, q )
      CALL thisStorage % GetQuadratureWeights( ws, wp, wq )

      Pfilt = ZERO
      V     = ZERO
      VInv  = ZERO

      ! The Legendre-polynomials in each direction are evaluated at each of the computational points
      DO c = 0, nQ
         DO b = 0, nP
            DO a = 0, nS ! Loop over the Legendre Polynomials
               
               row = a + 1 + (b + c*(nQ+1))*(nS+1)
               Lnorm = ZERO

               IF( a <= sCutoff .OR. b <= pCutoff .OR. c <= qCutoff )THEN
                  Pfilt(row,row) = ONE
               ENDIF

               DO k = 0, nQ

                  CALL LegendrePolynomial(c, q(k), Lk, Ls)

                  DO j = 0, nP

                     CALL LegendrePolynomial(b, p(j), Lj, Ls)

                     DO i = 0, nS

                        col = i + 1 + (j + k*(nQ+1))*(nS+1)

                        CALL LegendrePolynomial(a, s(i), Li, Ls)

                        Lnorm = Lnorm + ( (Li*Lj*Lk)**2 )*(ws(i)*wp(j)*wq(k))
                        V(row,col) = (Li*Lj*Lk)*(ws(i)*wp(j)*wq(k))
                        ! The inverse of the modal projection matrix is easy to build since the Legendre basis
                        ! is an orthogonal basis
                        VInv(col,row) = Li*Lj*Lk

                     ENDDO
                  ENDDO
               ENDDO
               V(row,1:nDOF) = V(row,1:nDOF)/Lnorm
            ENDDO
         ENDDO
      ENDDO

      Pfilt = MATMUL( Pfilt, V )
      thisFilter % filterMat = MATMUL( VInv, Pfilt )
    
      DEALLOCATE( V, Pfilt, Vinv )

 END SUBROUTINE Build_ModalCutoffFilter3D
!
!
!
 SUBROUTINE Trash_ModalCutoffFilter3D( thisFilter )
 ! S/R Trash
 ! 
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS(ModalCutoffFilter3D), INTENT(inout) :: thisFilter
   
      DEALLOCATE( thisFilter % filterMat )

 END SUBROUTINE Trash_ModalCutoffFilter3D
!
!
!==================================================================================================!
!--------------------------------------- Accessors ------------------------------------------------!
!==================================================================================================!
!
! 
 SUBROUTINE SetNumberOfNodes_ModalCutoffFilter3D( thisFilter, nS, nP, nQ  )
 ! S/R SetNumberOfNodes
 !
 !    
 !
 ! =============================================================================================== !
 ! DECLARATIONS
    IMPLICIT NONE
    CLASS(ModalCutoffFilter3D), INTENT(inout) :: thisFilter
    INTEGER, INTENT(in)                   :: nS, nP, nQ

      thisFilter % nS = nS
      thisFilter % nP = nP
      thisFilter % nQ = nQ
       
 END SUBROUTINE SetNumberOfNodes_ModalCutoffFilter3D
!
!
!
 SUBROUTINE GetNumberOfNodes_ModalCutoffFilter3D( thisFilter, nS, nP, nQ  )
 ! S/R GetNumberOfNodes
 !
 !
 ! =============================================================================================== !
 ! DECLARATIONS
    IMPLICIT NONE
    CLASS(ModalCutoffFilter3D), INTENT(in) :: thisFilter
    INTEGER, INTENT(out)               :: nS, nP, nQ

      nS = thisFilter % nS
      nP = thisFilter % nP
      nQ = thisFilter % nQ

 END SUBROUTINE GetNumberOfNodes_ModalCutoffFilter3D
!
!
!==================================================================================================!
!------------------------------------- Type Specific ----------------------------------------------!
!==================================================================================================!
!
! 
 FUNCTION ApplyFilter_ModalCutoffFilter3D( thisFilter, sol ) RESULT(fSol)
 ! FUNCTION ApplyFilter
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS(ModalCutoffFilter3D) :: thisFilter
   REAL(prec)                     :: sol(0:thisFilter % nS, 0:thisFilter % nP, 0:thisFilter % nQ)
   REAL(prec)                     :: fSol(0:thisFilter % nS, 0:thisFilter % nP, 0:thisFilter % nQ)
   ! LOCAL 
   REAL(prec) :: locSol(1:thisFilter % nDOF)
   INTEGER    :: nS, nP, nQ, nDOF 

      CALL thisFilter % GetNumberOfNodes( nS, nP, nQ )
      CALL thisFilter % GetNDOF( nDOF )

      locSol = Map3Dto1D( sol, nS, nP, nQ, nDOF )

      locSol = MATMUL( thisFilter % filterMat, locSol )
 
      fSol   = Map1Dto3D( locSol, nS, nP, nQ, nDOF )
      
 END FUNCTION ApplyFilter_ModalCutoffFilter3D

END MODULE ModalCutoffFilter3D_Class
