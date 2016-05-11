! RollOffFilter1D_Class.f90
! 
! Copyright 2015 Joseph Schoonover <schoonover.numerics@gmail.com>
! 
! RollOffFilter1D_Class.f90 is part of the Spectral Element Libraries in Fortran (SELF).
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

MODULE RollOffFilter1D_Class
! ========================================= Logs ================================================= !
! yyyy-mm-dd  Joe  <joe@clay>
!   This module provides a data structure for constructing and storing element-local filter matrices
!   that will be used for polynomial de-aliasing or as part of SGS parameterizations. 1-D arrays,
!   representative of the nodal values of an interpolant, can be passed in to the supplied routine
!   "ApplyFilter" to return a filtered form of the interpolant.  
!    In this version, a "Roll-off" filter is used. A polynomial "halfwidth" is required.
!    Only the Legendre basis is considered here.
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
! //////////////////////////////////////////////////////////////////////////////////////////////// !


USE ModelPrecision
USE ModelFlags
USE ConstantsDictionary
USE CommonRoutines

USE Legendre
USE Lagrange_1D_Class

USE NodalStorage_1D_Class


   TYPE RollOffFilter1D
      INTEGER                 :: nS
      INTEGER                 :: nSco ! modal cutoff 
      REAL(prec), ALLOCATABLE :: filterMat(:,:)
      
      CONTAINS
  
      PROCEDURE :: Build => Build_RollOffFilter1D
      PROCEDURE :: Trash => Trash_RollOffFilter1D
      ! 
      PROCEDURE :: SetNumberOfNodes => SetNumberOfNodes_RollOffFilter1D
      PROCEDURE :: GetNumberOfNodes => GetNumberOfNodes_RollOffFilter1D
      PROCEDURE :: SetCutoffs => SetCutoffs_RollOffFilter1D
      PROCEDURE :: GetCutoffs => GetCutoffs_RollOffFilter1D
      !
      PROCEDURE :: ApplyFilter => ApplyFilter_RollOffFilter1D
   END TYPE RollOffFilter1D


CONTAINS
!
!
!==================================================================================================!
!------------------------------- Manual Constructors/Destructors ----------------------------------!
!==================================================================================================!
!
!
 SUBROUTINE Build_RollOffFilter1D( thisFilter, thisStorage, sCutoff )
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
   CLASS(RollOffFilter1D), INTENT(inout) :: thisFilter
   TYPE(NodalStorage_1D), INTENT(in)              :: thisStorage
   INTEGER, INTENT(in)                            :: sCutoff
   ! Local
   REAL(prec) :: s(0:thisStorage % nS)
   REAL(prec) :: ws(0:thisStorage % nS)
   REAL(prec) :: Lnorm, L, Ls, r, alpha, sc
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

      sc = real(sCutoff,prec)
      alpha = log(TWO)/(sc*sc)
      ! The Legendre-polynomials in each direction are evaluated at each of the computational points
      DO row = 0, nS
         r = real(row,prec)
         Pfilt(row,row) = exp( -alpha*r**2 )
         
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

 END SUBROUTINE Build_RollOffFilter1D
!
!
!
 SUBROUTINE Trash_RollOffFilter1D( thisFilter )
 ! S/R Trash
 ! 
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS(RollOffFilter1D), INTENT(inout) :: thisFilter
   
      DEALLOCATE( thisFilter % filterMat )

 END SUBROUTINE Trash_RollOffFilter1D
!
!
!==================================================================================================!
!--------------------------------------- Accessors ------------------------------------------------!
!==================================================================================================!
!
! 
 SUBROUTINE SetNumberOfNodes_RollOffFilter1D( thisFilter, nS )
 ! S/R SetNumberOfNodes
 !
 !    
 !
 ! =============================================================================================== !
 ! DECLARATIONS
    IMPLICIT NONE
    CLASS(RollOffFilter1D), INTENT(inout) :: thisFilter
    INTEGER, INTENT(in)                            :: nS

      thisFilter % nS = nS
       
 END SUBROUTINE SetNumberOfNodes_RollOffFilter1D
!
!
!
 SUBROUTINE GetNumberOfNodes_RollOffFilter1D( thisFilter, nS  )
 ! S/R GetNumberOfNodes
 !
 !
 ! =============================================================================================== !
 ! DECLARATIONS
    IMPLICIT NONE
    CLASS(RollOffFilter1D), INTENT(in) :: thisFilter
    INTEGER, INTENT(out)                        :: nS

      nS = thisFilter % nS

 END SUBROUTINE GetNumberOfNodes_RollOffFilter1D
!
!
!
 SUBROUTINE SetCutoffs_RollOffFilter1D( thisFilter, sCutoff )
 !
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS(RollOffFilter1D), INTENT(inout) :: thisFilter
   INTEGER, INTENT(in)                            :: sCutoff
   
      thisFilter % nSco = sCutoff 
      
 END SUBROUTINE SetCutoffs_RollOffFilter1D
!
!
!
 SUBROUTINE GetCutoffs_RollOffFilter1D( thisFilter, sCutoff )
 !
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS(RollOffFilter1D), INTENT(in) :: thisFilter
   INTEGER, INTENT(out)                        :: sCutoff
   
      sCutoff = thisFilter % nSco
      
 END SUBROUTINE GetCutoffs_RollOffFilter1D
!
!
!==================================================================================================!
!------------------------------------- Type Specific ----------------------------------------------!
!==================================================================================================!
!
! 
 FUNCTION ApplyFilter_RollOffFilter1D( thisFilter, sol ) RESULT(fSol)
 ! FUNCTION ApplyFilter
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS(RollOffFilter1D) :: thisFilter
   REAL(prec)                 :: sol(0:thisFilter % nS)
   REAL(prec)                 :: fSol(0:thisFilter % nS)

      fSol = MATMUL( thisFilter % filterMat, sol )
      
 END FUNCTION ApplyFilter_RollOffFilter1D

END MODULE RollOffFilter1D_Class
