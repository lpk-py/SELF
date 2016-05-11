! RollOffFilter2D_Class.f90
! 
! Copyright 2015 Joseph Schoonover <schoonover.numerics@gmail.com>
! 
! RollOffFilter2D_Class.f90 is part of the Spectral Element Libraries in Fortran (SELF).
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
 
MODULE RollOffFilter2D_Class
! ========================================= Logs ================================================= !
! yyyy-mm-dd  Joe  <joe@clay>
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
! 
! //////////////////////////////////////////////////////////////////////////////////////////////// !


USE ModelPrecision
USE ModelFlags
USE ConstantsDictionary
USE CommonRoutines

USE Legendre
USE Lagrange_2D_Class

USE NodalStorage_2D_Class


   TYPE RollOffFilter2D
      INTEGER                 :: nS, nP, nDOF
      INTEGER                 :: nSco, nPco ! modal cutoff numbers
      REAL(prec), ALLOCATABLE :: filterMat(:,:)
      
      CONTAINS
  
      PROCEDURE :: Build => Build_RollOffFilter2D
      PROCEDURE :: Trash => Trash_RollOffFilter2D
      ! 
      PROCEDURE :: SetNDOF          => SetNDOF_RollOffFilter2D
      PROCEDURE :: GetNDOF          => GetNDOF_RollOffFilter2D
      PROCEDURE :: SetNumberOfNodes => SetNumberOfNodes_RollOffFilter2D
      PROCEDURE :: GetNumberOfNodes => GetNumberOfNodes_RollOffFilter2D

      PROCEDURE :: ApplyFilter => ApplyFilter_RollOffFilter2D
   END TYPE RollOffFilter2D


CONTAINS
!
!
!==================================================================================================!
!------------------------------- Manual Constructors/Destructors ----------------------------------!
!==================================================================================================!
!
!
 SUBROUTINE Build_RollOffFilter2D( thisFilter, thisStorage, sCutoff, pCutoff )
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
   CLASS(RollOffFilter2D), INTENT(inout) :: thisFilter
   TYPE(NodalStorage_2D), INTENT(in)     :: thisStorage
   INTEGER, INTENT(in)                   :: sCutoff, pCutoff
   ! Local
   REAL(prec) :: s(0:thisStorage % nS), p(0:thisStorage % nP)
   REAL(prec) :: ws(0:thisStorage % nS), wp(0:thisStorage % nP)
   REAL(prec) :: Lnorm, Li, Lj, Ls, alpha, r, sc
   REAL(prec), ALLOCATABLE :: Pfilt(:,:), V(:,:), VInv(:,:)
   INTEGER    :: nS, nP, i, j, a, b, row, col, nDOF 
   
      CALL thisStorage % GetNumberOfNodes( nS, nP )
      CALL thisFilter % SetNumberOfNodes( nS, nP )
      nDOF = (nS+1)*(nP+1)
      CALL thisFilter % SetNDOF( nDOF )
      thisFilter % nSco = sCutoff
      thisFilter % nPco = pCutoff

      ALLOCATE( thisFilter % filterMat(1:nDOF,1:nDOF) )

      ALLOCATE( Pfilt(1:nDOF,1:nDOF), V(1:nDOF,1:nDOF), VInv(1:nDOF,1:nDOF) )

      thisFilter % filterMat = ZERO

      CALL thisStorage % interp % GetNodes( s, p )
      CALL thisStorage % GetQuadratureWeights( ws, wp )

      Pfilt = ZERO
      V     = ZERO
      VInv  = ZERO

      sc = real(scutoff,prec)
      alpha = log(TWO)/(sc*sc)
      ! The Legendre-polynomials in each direction are evaluated at each of the computational points
      DO b = 0, nP
         DO a = 0, nS ! Loop over the Legendre Polynomials
               
            row = a + 1 + b*(nS+1)
            Lnorm = ZERO

            r = real(a**2+b**2,prec)
            Pfilt(row,row) = exp( -alpha*r )

            DO j = 0, nP

               CALL LegendrePolynomial(b, p(j), Lj, Ls)

               DO i = 0, nS

                  col = i + 1 + j*(nS+1)

                  CALL LegendrePolynomial(a, s(i), Li, Ls)

                  Lnorm = Lnorm + ( (Li*Lj)**2 )*(ws(i)*wp(j))
                  V(row,col) = (Li*Lj)*(ws(i)*wp(j))
                  ! The inverse of the modal projection matrix is easy to build since the Legendre basis
                  ! is an orthogonal basis
                  VInv(col,row) = Li*Lj

               ENDDO
            ENDDO
            V(row,1:nDOF) = V(row,1:nDOF)/Lnorm
         ENDDO
      ENDDO

      Pfilt = MATMUL( Pfilt, V )
      thisFilter % filterMat = MATMUL( VInv, Pfilt )
    
      DEALLOCATE( V, Pfilt, Vinv )

 END SUBROUTINE Build_RollOffFilter2D
!
!
!
 SUBROUTINE Trash_RollOffFilter2D( thisFilter )
 ! S/R Trash
 ! 
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS(RollOffFilter2D), INTENT(inout) :: thisFilter
   
      DEALLOCATE( thisFilter % filterMat )

 END SUBROUTINE Trash_RollOffFilter2D
!
!
!==================================================================================================!
!--------------------------------------- Accessors ------------------------------------------------!
!==================================================================================================!
!
!
 SUBROUTINE SetNDOF_RollOffFilter2D( thisFilter, nDOF  )
 ! S/R SetNDOF
 !
 !    
 !
 ! =============================================================================================== !
 ! DECLARATIONS
    IMPLICIT NONE
    CLASS(RollOffFilter2D), INTENT(inout) :: thisFilter
    INTEGER, INTENT(in)                   :: nDOF

      thisFilter % nDOF = nDOF
       
 END SUBROUTINE SetNDOF_RollOffFilter2D
!
!
!
 SUBROUTINE GetNDOF_RollOffFilter2D( thisFilter, nDOF  )
 ! S/R GetNDOF
 !
 !
 ! =============================================================================================== !
 ! DECLARATIONS
    IMPLICIT NONE
    CLASS(RollOffFilter2D), INTENT(in) :: thisFilter
    INTEGER, INTENT(out)                   :: nDOF

      nDOF = thisFilter % nDOF

 END SUBROUTINE GetNDOF_RollOffFilter2D
!
!
! 
 SUBROUTINE SetNumberOfNodes_RollOffFilter2D( thisFilter, nS, nP  )
 ! S/R SetNumberOfNodes
 !
 !    
 !
 ! =============================================================================================== !
 ! DECLARATIONS
    IMPLICIT NONE
    CLASS(RollOffFilter2D), INTENT(inout) :: thisFilter
    INTEGER, INTENT(in)                   :: nS, nP

      thisFilter % nS = nS
      thisFilter % nP = nP
       
 END SUBROUTINE SetNumberOfNodes_RollOffFilter2D
!
!
!
 SUBROUTINE GetNumberOfNodes_RollOffFilter2D( thisFilter, nS, nP  )
 ! S/R GetNumberOfNodes
 !
 !
 ! =============================================================================================== !
 ! DECLARATIONS
    IMPLICIT NONE
    CLASS(RollOffFilter2D), INTENT(in) :: thisFilter
    INTEGER, INTENT(out)                   :: nS, nP

      nS = thisFilter % nS
      nP = thisFilter % nP

 END SUBROUTINE GetNumberOfNodes_RollOffFilter2D
!
!
!==================================================================================================!
!------------------------------------- Type Specific ----------------------------------------------!
!==================================================================================================!
!
! 
 FUNCTION ApplyFilter_RollOffFilter2D( thisFilter, sol ) RESULT(fSol)
 ! FUNCTION ApplyFilter
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS(RollOffFilter2D) :: thisFilter
   REAL(prec)                     :: sol(0:thisFilter % nS, 0:thisFilter % nP)
   REAL(prec)                     :: fSol(0:thisFilter % nS, 0:thisFilter % nP)
   ! LOCAL 
   REAL(prec) :: locSol(1:thisFilter % nDOF)
   INTEGER    :: nS, nP, nDOF 

      CALL thisFilter % GetNumberOfNodes( nS, nP )
      CALL thisFilter % GetNDOF( nDOF )

      locSol = Map2Dto1D( sol, nS, nP, nDOF )

      locSol = MATMUL( thisFilter % filterMat, locSol )
 
      fSol   = Map1Dto2D( locSol, nS, nP, nDOF )
      
 END FUNCTION ApplyFilter_RollOffFilter2D

END MODULE RollOffFilter2D_Class
