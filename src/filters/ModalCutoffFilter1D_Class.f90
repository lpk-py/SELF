! ModalCutoffFilter1D_Class.f90 ( new with v2.1 - 21 March 2016)
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
!   that will be used for polynomial de-aliasing or as part of SGS parameterizations. 1-D arrays,
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



MODULE ModalCutoffFilter1D_Class

USE ModelPrecision
USE ModelFlags
USE ConstantsDictionary
USE CommonRoutines

USE Legendre
USE Lagrange_1D_Class

USE NodalStorage_1D_Class


   TYPE ModalCutoffFilter1D
      INTEGER                 :: nS
      INTEGER                 :: nSco ! modal cutoff 
      REAL(prec), ALLOCATABLE :: filterMat(:,:)
      
      CONTAINS
  
      PROCEDURE :: Build => Build_ModalCutoffFilter1D
      PROCEDURE :: Trash => Trash_ModalCutoffFilter1D
      ! 
      PROCEDURE :: SetNumberOfNodes => SetNumberOfNodes_ModalCutoffFilter1D
      PROCEDURE :: GetNumberOfNodes => GetNumberOfNodes_ModalCutoffFilter1D
      PROCEDURE :: SetCutoffs => SetCutoffs_ModalCutoffFilter1D
      PROCEDURE :: GetCutoffs => GetCutoffs_ModalCutoffFilter1D
      !
      PROCEDURE :: ApplyFilter => ApplyFilter_ModalCutoffFilter1D
   END TYPE ModalCutoffFilter1D


CONTAINS
!
!
!==================================================================================================!
!------------------------------- Manual Constructors/Destructors ----------------------------------!
!==================================================================================================!
!
!
 SUBROUTINE Build_ModalCutoffFilter1D( thisFilter, thisStorage, sCutoff )
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
   CLASS(ModalCutoffFilter1D), INTENT(inout) :: thisFilter
   TYPE(NodalStorage_1D), INTENT(in)              :: thisStorage
   INTEGER, INTENT(in)                            :: sCutoff
   ! Local
   REAL(prec) :: s(0:thisStorage % nS)
   REAL(prec) :: ws(0:thisStorage % nS)
   REAL(prec) :: Lnorm, L, Ls
   REAL(prec), ALLOCATABLE :: Pfilt(:,:), V(:,:), VInv(:,:)
   INTEGER    :: nS, row, col 
   
      CALL thisStorage % GetNumberOfNodes( nS )
      CALL thisFilter % SetNumberOfNodes( nS )

      CALL thisFilter % SetCutoffs( sCutoff )

      ALLOCATE( thisFilter % filterMat(0:nS,0:nS) )

      ALLOCATE( Pfilt(0:nS,0:nS), V(0:nS,0:nS), VInv(0:nS,0:nS) )

      thisFilter % filterMat = ZERO

      CALL thisStorage % interp % GetNodes( s )
      CALL thisStorage % GetQuadratureWeights( ws )

      Pfilt = ZERO
      V     = ZERO
      VInv  = ZERO

      ! The Legendre-polynomials in each direction are evaluated at each of the computational points
      DO row = 0, nS
         IF( row <= sCutoff )THEN
            Pfilt(row,row) = ONE
         ENDIF
         Lnorm = ZERO
         DO col = 0, nS

            CALL LegendrePolynomial(row, s(col), L, Ls)
            Lnorm = Lnorm + L*L*ws(col)
            V(row,col) = L*ws(col)
            ! The inverse of the modal projection matrix is easy to build since the Legendre basis
            ! is an orthogonal basis
            VInv(col,row) = L

         ENDDO
         V(row,0:nS) = V(row,0:nS)/Lnorm
      ENDDO

      Pfilt = MATMUL( Pfilt, V )
      thisFilter % filterMat = MATMUL( VInv, Pfilt )
      
    
      DEALLOCATE( V, Pfilt, Vinv )

 END SUBROUTINE Build_ModalCutoffFilter1D
!
!
!
 SUBROUTINE Trash_ModalCutoffFilter1D( thisFilter )
 ! S/R Trash
 ! 
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS(ModalCutoffFilter1D), INTENT(inout) :: thisFilter
   
      DEALLOCATE( thisFilter % filterMat )

 END SUBROUTINE Trash_ModalCutoffFilter1D
!
!
!==================================================================================================!
!--------------------------------------- Accessors ------------------------------------------------!
!==================================================================================================!
!
! 
 SUBROUTINE SetNumberOfNodes_ModalCutoffFilter1D( thisFilter, nS )
 ! S/R SetNumberOfNodes
 !
 !    
 !
 ! =============================================================================================== !
 ! DECLARATIONS
    IMPLICIT NONE
    CLASS(ModalCutoffFilter1D), INTENT(inout) :: thisFilter
    INTEGER, INTENT(in)                            :: nS

      thisFilter % nS = nS
       
 END SUBROUTINE SetNumberOfNodes_ModalCutoffFilter1D
!
!
!
 SUBROUTINE GetNumberOfNodes_ModalCutoffFilter1D( thisFilter, nS  )
 ! S/R GetNumberOfNodes
 !
 !
 ! =============================================================================================== !
 ! DECLARATIONS
    IMPLICIT NONE
    CLASS(ModalCutoffFilter1D), INTENT(in) :: thisFilter
    INTEGER, INTENT(out)                        :: nS

      nS = thisFilter % nS

 END SUBROUTINE GetNumberOfNodes_ModalCutoffFilter1D
!
!
!
 SUBROUTINE SetCutoffs_ModalCutoffFilter1D( thisFilter, sCutoff )
 !
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS(ModalCutoffFilter1D), INTENT(inout) :: thisFilter
   INTEGER, INTENT(in)                            :: sCutoff
   
      thisFilter % nSco = sCutoff 
      
 END SUBROUTINE SetCutoffs_ModalCutoffFilter1D
!
!
!
 SUBROUTINE GetCutoffs_ModalCutoffFilter1D( thisFilter, sCutoff )
 !
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS(ModalCutoffFilter1D), INTENT(in) :: thisFilter
   INTEGER, INTENT(out)                        :: sCutoff
   
      sCutoff = thisFilter % nSco
      
 END SUBROUTINE GetCutoffs_ModalCutoffFilter1D
!
!
!==================================================================================================!
!------------------------------------- Type Specific ----------------------------------------------!
!==================================================================================================!
!
! 
 FUNCTION ApplyFilter_ModalCutoffFilter1D( thisFilter, sol ) RESULT(fSol)
 ! FUNCTION ApplyFilter
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS(ModalCutoffFilter1D) :: thisFilter
   REAL(prec)                 :: sol(0:thisFilter % nS)
   REAL(prec)                 :: fSol(0:thisFilter % nS)

      fSol = MATMUL( thisFilter % filterMat, sol )
      
 END FUNCTION ApplyFilter_ModalCutoffFilter1D

END MODULE ModalCutoffFilter1D_Class