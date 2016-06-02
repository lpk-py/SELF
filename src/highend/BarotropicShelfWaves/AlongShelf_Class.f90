! AlongShelf_Class.f90
! 
! Copyright 2015 Joseph Schoonover <schoonover.numerics@gmail.com>
! 
! AlongShelf_Class.f90 is part of the Spectral Element Libraries in Fortran (SELF).
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
 

MODULE AlongShelfClass
! ========================================= Logs ================================================= !
!2016-05-31  Joseph Schoonover  schoonover.numerics@gmail.com 
!
! //////////////////////////////////////////////////////////////////////////////////////////////// ! 
! src/common/
USE ModelPrecision
USE ConstantsDictionary
USE ModelFlags
USE CommonRoutines
! src/interp/
USE Lagrange_1D_Class
! src/nodal/
USE NodalStorage_1D_Class
USE DGSEMSolutionStorageClass_1D   
! src/highend/BarotropicShelfWaves/
USE BarotropicShelfWaves_Class
USE BarotropicShelfWavesParams_Class

! Parallel Libraries
USE OMP_LIB


IMPLICIT NONE
!
! Some parameters that are specific to this module


! ================================================================================================ !
!     The AlongShelf Class is provides an example of the 1-D DGSEM with an adaptive filter applied
! ================================================================================================ !  



    TYPE AlongShelf1D
      INTEGER                               :: nPlot, nS, nElems, nEq
      REAL(prec)                            :: dx, vbar
      REAL(prec), ALLOCATABLE               :: x(:,:)

      TYPE( BarotropicShelfWaves )          :: shelfwaves
      TYPE( NodalStorage_1D )               :: dGStorage
      TYPE( DGSEMSolution_1D ), ALLOCATABLE :: sol(:)
      REAL(prec), ALLOCATABLE               :: shelfShaper(:,:)
      REAL(prec), ALLOCATABLE               :: plMatS(:,:)

      CONTAINS

      PROCEDURE :: Build => Build_AlongShelf1D
      PROCEDURE :: Trash => Trash_AlongShelf1D
      PROCEDURE :: BuildMesh => BuildMesh_AlongShelf1D
      
      ! DGSEMSolution_2DStorage_Class Wrapper Routines
      PROCEDURE :: GetSolution => GetSolution_AlongShelf1D
      PROCEDURE :: SetSolution => SetSolution_AlongShelf1D
      PROCEDURE :: GetSolutionAtNode => GetSolutionAtNode_AlongShelf1D
      PROCEDURE :: SetSolutionAtNode => SetSolutionAtNode_AlongShelf1D
      PROCEDURE :: GetTendency => GetTendency_AlongShelf1D
      PROCEDURE :: SetTendency => SetTendency_AlongShelf1D
      PROCEDURE :: GetTendencyAtNode => GetTendencyAtNode_AlongShelf1D
      PROCEDURE :: SetTendencyAtNode => SetTendencyAtNode_AlongShelf1D
      PROCEDURE :: GetBoundarySolutionAtBoundary => GetBoundarySolutionAtBoundary_AlongShelf1D
      PROCEDURE :: SetBoundarySolutionAtBoundary => SetBoundarySolutionAtBoundary_AlongShelf1D
      PROCEDURE :: GetBoundaryFluxAtBoundary => GetBoundaryFluxAtBoundary_AlongShelf1D
      PROCEDURE :: SetBoundaryFluxAtBoundary => SetBoundaryFluxAtBoundary_AlongShelf1D
      
       ! Type Specific Routines
      PROCEDURE :: ShelfSteepeningShape => ShelfSteepeningShape_AlongShelf1D
      PROCEDURE :: FillInShelfShaper => FillInShelfShaper_AlongShelf1D
      PROCEDURE :: GlobalTimeDerivative => GlobalTimeDerivative_AlongShelf1D
      PROCEDURE :: ForwardStepRK3 => ForwardStepRK3_AlongShelf1D
      PROCEDURE :: EdgeFlux => EdgeFlux_AlongShelf1D
      PROCEDURE :: MappedTimeDerivative => MappedTimeDerivative_AlongShelf1D 
      PROCEDURE :: CalculateSolutionAtBoundaries => CalculateSolutionAtBoundaries_AlongShelf1D
      
      PROCEDURE :: CoarseToFine => CoarseToFine_AlongShelf1D
      PROCEDURE :: WriteTecplot => WriteTecplot_AlongShelf1D
      PROCEDURE :: WritePickup => WritePickup_AlongShelf1D
      PROCEDURE :: ReadPickup => ReadPickup_AlongShelf1D

    END TYPE AlongShelf1D

 CONTAINS
!
!
!==================================================================================================!
!------------------------------- Manual Constructors/Destructors ----------------------------------!
!==================================================================================================!
!
!
 SUBROUTINE Build_AlongShelf1D( myDGSEM )!, forceNoPickup )
 ! S/R Build
 ! 
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS(AlongShelf1D), INTENT(inout) :: myDGSEM
!   INTEGER, INTENT(in), OPTIONAL   :: forceNoPickup
   !LOCAL
   INTEGER                 :: iEl, nS, nPlot, nElems
   REAL(prec)              :: dx, nR, nC
   REAL(prec), ALLOCATABLE :: sNew(:)

      CALL myDGSEM % shelfwaves % Build( )
      
      nS     = myDGSEM % shelfwaves % params % polyDeg
      nElems = myDGSEM % shelfwaves % params % nYElems
      dx     = myDGSEM % shelfwaves % params % dy
 
      myDGSEM % nElems  = nElems
      myDGSEM % dx      = dx
      myDGSEM % nS      = nS
      myDGSEM % nPlot   = myDGSEM % shelfwaves % params % nPlot
      myDGSEM % nEq     = myDGSEM % shelfwaves % nEigenpairs
      myDGSEM % vbar    = myDGSEM % shelfwaves % params % vbar

      nPlot = myDGSEM % shelfwaves % params % nPlot
      
      ALLOCATE( sNew(0:nPlot) )
      CALL myDGSEM % dGStorage % Build( nS, GAUSS, DG )

      CALL myDGSEM % BuildMesh(  )
      
      ! Space is allocated for the solution-storage structure and an blank structure is built
   
      ALLOCATE( myDGSEM % sol(1:nElems), myDGSEM % shelfShaper(0:nS,1:nElems) )
      myDGSEM % shelfShaper = ZERO

      ! Build and initialize the solution storage structure to zero
      DO iEl = 1, myDGSEM % nElems
         CALL myDGSEM % sol(iEl) % Build( nS, myDGSEM % shelfwaves % nEigenpairs )
      ENDDO
      
      ! In the event that this is a pickup run, we'll READ the pickup file here for the solution 
      ! and the relaxation fields. A negative initial iterate can be used to specify to start 
      ! from zeros.
      
      IF( myDGSEM % shelfwaves % params % iterInit /= 0 )then
         
         ! This call reads the solution from a pickup file if desired
         CALL myDGSEM % ReadPickup( myDGSEM % shelfwaves % params % iterInit )

      ENDIF
      
      sNew = UniformPoints( -ONE, ONE, nPlot )
      ALLOCATE( myDGSEM % plMatS(0:nPlot,0:nS) )
      ! Build the plotting matrix
      CALL myDGSEM % dgStorage % interp % CalculateInterpolationMatrix( nPlot, sNew, myDGSEM % plMatS )
                                                                        
      DEALLOCATE(sNew)

      CALL myDGSEM % FillInShelfShaper( )
      
 END SUBROUTINE Build_AlongShelf1D
!
!
!
 SUBROUTINE Trash_AlongShelf1D( myDGSEM )
 ! S/R Trash
 ! 
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS(AlongShelf1D), INTENT(inout) :: myDGSEM
   ! LOCAL
   INTEGER :: iEl

     DO iEl = 1, myDGSEM % nElems
        CALL myDGSEM % sol(iEl) % Trash( )
     ENDDO

     CALL myDGSEM % dGStorage % Trash( )
     CALL myDGSEM % shelfwaves % Trash( )

     DEALLOCATE( myDGSEM % sol, myDGSEM % shelfShaper ) 
     DEALLOCATE( myDGSEM % plMatS )

 END SUBROUTINE Trash_AlongShelf1D
!
!
!
 SUBROUTINE BuildMesh_AlongShelf1D( myDGSEM )
 ! S/R BuildMesh
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS( AlongShelf1D ), INTENT(inout) :: myDGSEM
   ! LOCAL
   INTEGER :: iEl, nEl, i, nS
   REAL(prec) :: dx, x1, x2
   REAL(prec) :: s(0:myDGSEM % nS)
   
      PRINT*,'Module AlongShelf1D_Class.f90 : S/R BuildMesh :'
      PRINT*,'  Loading default mesh.'
      
      dx  = myDGSEM % dx
      nS  = myDGSEM % nS
      nEl = myDGSEM % nElems   

      ALLOCATE( myDGSEM % x(0:nS, 1:nEl) )
      CALL myDGSEM % dgStorage % GetNodes( s )
      x1 = ZERO
      DO iEl = 1, nEl
 
         x2 = x1 + dx

         DO i = 0, nS
            myDGSEM % x(i,iEl) = x1 + HALF*dx*( s(i) + ONE )
         ENDDO

         x1 = x2

      ENDDO

      PRINT *,'  Done!'

 END SUBROUTINE BuildMesh_AlongShelf1D
!
!
!==================================================================================================!
!--------------------------------------- Accessors ------------------------------------------------!
!==================================================================================================!
!
!
!
! ---------------------------------------- Solution ---------------------------------------------- !
! 
 SUBROUTINE GetSolution_AlongShelf1D( myDGSEM, iEl, theSolution  )
 ! S/R GetSolution
 !  
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS(AlongShelf1D), INTENT(in) :: myDGSEM
   INTEGER, INTENT(in)          :: iEl
   REAL(prec), INTENT(out)      :: theSolution(0:myDGSEM % nS, 1:myDGSEM % nEq)

      CALL myDGSEM % sol(iEl) % GetSolution( theSolution )

 END SUBROUTINE GetSolution_AlongShelf1D
!
!
!
 SUBROUTINE SetSolution_AlongShelf1D( myDGSEM, iEl, theSolution  )
 ! S/R SetSolution
 !  
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS(AlongShelf1D), INTENT(inout) :: myDGSEM
   INTEGER, INTENT(in)             :: iEl
   REAL(prec), INTENT(in)          :: theSolution(0:myDGSEM % nS, 1:myDGSEM % nEq)

      CALL myDGSEM % sol(iEl) % SetSolution( theSolution )

 END SUBROUTINE SetSolution_AlongShelf1D
!
!
!
 SUBROUTINE GetSolutionAtNode_AlongShelf1D( myDGSEM, iEl, i, theSolution  )
 ! S/R GetSolution
 !  
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS(AlongShelf1D), INTENT(in) :: myDGSEM
   INTEGER, INTENT(in)          :: iEl
   INTEGER, INTENT(in)          :: i
   REAL(prec), INTENT(out)      :: theSolution(1:myDGSEM % nEq)

      CALL myDGSEM % sol(iEl) % GetSolutionAtNode( i, theSolution )

 END SUBROUTINE GetSolutionAtNode_AlongShelf1D
!
!
!
 SUBROUTINE SetSolutionAtNode_AlongShelf1D( myDGSEM, iEl, i, theSolution  )
 ! S/R SetSolution
 !  
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS(AlongShelf1D), INTENT(inout) :: myDGSEM
   INTEGER, INTENT(in)             :: iEl
   INTEGER, INTENT(in)             :: i
   REAL(prec), INTENT(in)          :: theSolution(1:myDGSEM % nEq)

      CALL myDGSEM % sol(iEl) % SetSolutionAtNode( i, theSolution )

 END SUBROUTINE SetSolutionAtNode_AlongShelf1D
!
! ----------------------------------------- Tendency --------------------------------------------- !
!
 SUBROUTINE GetTendency_AlongShelf1D( myDGSEM, iEl, theTend  )
 ! S/R GetTendency
 !  
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS(AlongShelf1D), INTENT(in) :: myDGSEM
   INTEGER, INTENT(in)          :: iEl
   REAL(prec), INTENT(out)      :: theTend(0:myDGSEM % nS, 1:myDGSEM % nEq)

      CALL myDGSEM % sol(iEl) % GetTendency( theTend )

 END SUBROUTINE GetTendency_AlongShelf1D
!
!
!
 SUBROUTINE SetTendency_AlongShelf1D( myDGSEM, iEl, theTend  )
 ! S/R SetTendency
 !  
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS(AlongShelf1D), INTENT(inout) :: myDGSEM
   INTEGER, INTENT(in)             :: iEl
   REAL(prec), INTENT(in)          :: theTend(0:myDGSEM % nS, 1:myDGSEM % nEq)

      CALL myDGSEM % sol(iEl) % SetTendency( theTend )   

 END SUBROUTINE SetTendency_AlongShelf1D
!
!
!
 SUBROUTINE GetTendencyAtNode_AlongShelf1D( myDGSEM, iEl, i, theTend  )
 ! S/R GetTendencyAtNode
 !  
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS(AlongShelf1D), INTENT(in) :: myDGSEM
   INTEGER, INTENT(in)          :: iEl
   INTEGER, INTENT(in)          :: i
   REAL(prec), INTENT(out)      :: theTend(1:myDGSEM % nEq)

      CALL myDGSEM % sol(iEl) % GetTendencyAtNode( i, theTend )

 END SUBROUTINE GetTendencyAtNode_AlongShelf1D
!
!
!
 SUBROUTINE SetTendencyAtNode_AlongShelf1D( myDGSEM, iEl, i,  theTend  )
 ! S/R SetTendencyAtNode
 !  
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS(AlongShelf1D), INTENT(inout) :: myDGSEM
   INTEGER, INTENT(in)             :: iEl
   INTEGER, INTENT(in)             :: i
   REAL(prec), INTENT(in)          :: theTend(1:myDGSEM % nEq)

      CALL myDGSEM % sol(iEl) % SetTendencyAtNode( i, theTend )

 END SUBROUTINE SetTendencyAtNode_AlongShelf1D
!
! ------------------------------- Boundary Solution ---------------------------------------------- !
!
 SUBROUTINE GetBoundarySolutionAtBoundary_AlongShelf1D( myDGSEM, iEl, boundary, theSolution  )
 ! S/R GetBoundarySolution
 !  
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS(AlongShelf1D), INTENT(in) :: myDGSEM
   INTEGER, INTENT(in)          :: iEl
   INTEGER, INTENT(in)          :: boundary
   REAL(prec), INTENT(out)      :: theSolution(1:myDGSEM % nEq)

      CALL myDGSEM % sol(iEl) % GetBoundarySolutionAtBoundary( boundary, theSolution )

 END SUBROUTINE GetBoundarySolutionAtBoundary_AlongShelf1D
!
!
!
 SUBROUTINE SetBoundarySolutionAtBoundary_AlongShelf1D( myDGSEM, iEl, boundary, theSolution  )
 ! S/R SetBoundarySolution
 !  
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS(AlongShelf1D), INTENT(inout) :: myDGSEM
   INTEGER, INTENT(in)             :: iEl
   INTEGER, INTENT(in)             :: boundary
   REAL(prec), INTENT(in)          :: theSolution(1:myDGSEM % nEq)

      CALL myDGSEM % sol(iEl) % SetBoundarySolutionAtBoundary( boundary, theSolution )

 END SUBROUTINE SetBoundarySolutionAtBoundary_AlongShelf1D
!
! --------------------------------- Boundary Flux ------------------------------------------------ !
!
 SUBROUTINE GetBoundaryFluxAtBoundary_AlongShelf1D( myDGSEM, iEl, boundary, theFlux  )
 ! S/R GetBoundaryFlux
 !  
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS(AlongShelf1D), INTENT(in) :: myDGSEM
   INTEGER, INTENT(in)          :: iEl
   INTEGER, INTENT(in)          :: boundary
   REAL(prec), INTENT(out)      :: theFlux(1:myDGSEM % nEq)

      CALL myDGSEM % sol(iEl) % GetBoundaryFluxAtBoundary( boundary, theFlux )

 END SUBROUTINE GetBoundaryFluxAtBoundary_AlongShelf1D
!
!
!
 SUBROUTINE SetBoundaryFluxAtBoundary_AlongShelf1D( myDGSEM, iEl, boundary, theFlux  )
 ! S/R SetBoundaryFlux
 !  
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS(AlongShelf1D), INTENT(inout) :: myDGSEM
   INTEGER, INTENT(in)             :: iEl
   INTEGER, INTENT(in)             :: boundary
   REAL(prec), INTENT(in)          :: theFlux(1:myDGSEM % nEq)

      CALL myDGSEM % sol(iEl) % SetBoundaryFluxAtBoundary( boundary, theFlux )

 END SUBROUTINE SetBoundaryFluxAtBoundary_AlongShelf1D
!
!
!==================================================================================================!
!--------------------------------- Type Specific Routines -----------------------------------------!
!==================================================================================================!
!
!
 FUNCTION ShelfSteepeningShape_AlongShelf1D( myDGSEM, y ) RESULT( S )
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS(AlongShelf1D) :: myDGSEM
   REAL(prec)          :: y, S
   ! Local
   REAL(prec) :: Ly, yc

      Ly = myDGSEM % shelfwaves % params % steepeningZoneLength
      yc = myDGSEM % shelfwaves % params % steepeningCenter

      S = HALF*( tanh( (y-yc)/Ly ) + ONE )

 END FUNCTION ShelfSteepeningShape_AlongShelf1D
!
!
!
 SUBROUTINE FillInShelfShaper_AlongShelf1D( myDGSEM )
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS(AlongShelf1D), INTENT(inout) :: myDGSEM
   ! Local
   INTEGER    :: iEl, nEl, iS, nS
   REAL(prec), ALLOCATABLE :: dMat(:,:)
   REAL(prec) :: slopesteep(0:myDGSEM % nS)
   REAL(prec) :: y, dy

      nEl = myDGSEM % nElems
      nS  = myDGSEM % nS
      dy  = myDGSEM % dx
      CALL myDGSEM % dgStorage % interp % CalculateDerivativeMatrix( dMat )
      slopesteep = ZERO

      DO iEl = 1, nEl
         DO iS = 0, nS

            y = myDGSEM % x(iS,iEl)
            slopesteep(iS) = myDGSEM % ShelfSteepeningShape( y )

         ENDDO
         
         myDGSEM % shelfShaper(0:nS,iEl) = MATMUL( dMat, slopeSteep )*TWO/dy

      ENDDO

      DEALLOCATE( dMat )

 END SUBROUTINE FillInShelfShaper_AlongShelf1D
!
!
!
 SUBROUTINE CalculateSolutionAtBoundaries_AlongShelf1D( myDGSEM, iEl )
 ! S/R CalculateSolutionAtBoundaries
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS(AlongShelf1D), INTENT(inout) :: myDGSEM
   INTEGER, INTENT(in)             :: iEl
   
      CALL myDGSEM % sol(iEl) % CalculateSolutionAtBoundaries( myDGSEM % dgStorage )
      
 END SUBROUTINE CalculateSolutionAtBoundaries_AlongShelf1D
!
!
!
 SUBROUTINE ForwardStepRK3_AlongShelf1D( myDGSEM, tn )
 ! S/R ForwardStepRK3( 3rd order Runge-Kutta)
 ! 
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS(AlongShelf1D), INTENT(inout) :: myDGSEM
   REAL(prec), INTENT(in)          :: tn
   ! LOCAL
   REAL(prec) :: t, dt
   REAL(prec) :: G2D(0:myDGSEM % nS, 1:myDGSEM % nEq, 1:myDGSEM %  nElems) 
   REAL(prec) :: dSdt(0:myDGSEM % nS, 1:myDGSEM % nEq )
   INTEGER    :: m, iEl


     dt = myDGSEM % shelfwaves % params % dt
     G2D = ZERO
    
     DO m = 1,3 ! Loop over RK3 steps

        t = tn + rk3_b(m)*dt
        ! Calculate the tendency
        CALL myDGSEM % GlobalTimeDerivative( t )
        
        DO iEl = 1, myDGSEM % nElems ! Loop over all of the elements

           CALL myDGSEM % GetTendency( iEl, dSdt )
           G2D(:,:,iEl) = rk3_a(m)*G2D(:,:, iEl) + dSdt

           myDGSEM % sol(iEl) % solution = myDGSEM % sol(iEl) % solution + rk3_g(m)*dt*G2D(:,:,iEl)

         ENDDO ! iEl, loop over all of the elements
         
         
      ENDDO ! m, loop over the RK3 steps
   
   
       
 END SUBROUTINE ForwardStepRK3_AlongShelf1D
!
!
!
 SUBROUTINE GlobalTimeDerivative_AlongShelf1D( myDGSEM, tn ) 
 ! S/R GlobalTimeDerivative_AlongShelf1D
 ! 
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS(AlongShelf1D), INTENT(inout) :: myDGSEM
   REAL(prec), INTENT(in)          :: tn
   ! Local
   INTEGER :: iEl, iEdge

      ! CALL myDGSEM % LoadVelocityField( tn )
!$OMP PARALLEL


      ! Calculate the solution at the boundaries
!$OMP DO
      DO iEl = 1, myDGSEM % nElems
         CALL myDGSEM % CalculateSolutionAtBoundaries( iEl ) 
      ENDDO 
!$OMP END DO 
!$OMP FLUSH( myDGSEM )


!$OMP DO
      DO iEdge = 1, myDGSEM % nElems + 1
         CALL myDGSEM % EdgeFlux( iEdge, tn )
      ENDDO 
!$OMP END DO
!$OMP FLUSH( myDGSEM ) 
 

!$OMP DO
      DO iEl = 1, myDGSEM % nElems
         CALL myDGSEM % MappedTimeDerivative( iEl, tn )
      ENDDO
!$OMP END DO
!$OMP FLUSH( myDGSEM )
!$OMP END PARALLEL

 
 END SUBROUTINE GlobalTimeDerivative_AlongShelf1D
!
!
!
 SUBROUTINE EdgeFlux_AlongShelf1D( myDGSEM, iEdge, tn )
 ! S/R EdgeFlux
 ! 
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS( AlongShelf1D ), INTENT(inout) :: myDGSEM
   INTEGER, INTENT(in)               :: iEdge  
   REAL(prec), INTENT(in)            :: tn
   ! Local
   REAL(prec) :: flux(1:myDGSEM % nEq)
   REAL(prec) :: inState(1:myDGSEM % nEq)
   REAL(prec) :: exState(1:myDGSEM % nEq)
   REAL(prec) :: netV(1:myDGSEM % nEq)
   REAL(prec) :: nHat

      netV = myDGSEM % vbar - myDGSEM % shelfwaves % eigenvalues

      IF( iEdge == 1 )THEN! this is the left-most boundary "edge"
      
         CALL myDGSEM % GetBoundarySolutionAtBoundary( iEdge, 1, inState )
         exState = GetExternalState( tn, iEdge, myDGSEM % nElems, inState, myDGSEM % nEq )
         nHat = -ONE

         flux = RiemannSolver( inState, exState, netV, nHat, myDGSEM % nEq )

         ! Store the flux for the elements which share this edge
         CALL myDGSEM % SetBoundaryFluxAtBoundary( iEdge, 1, flux )

      ELSEIF( iEdge == myDGSEM % nElems + 1 )THEN  ! This is the right most boundary edge
 
         CALL myDGSEM % GetBoundarySolutionAtBoundary( iEdge-1, 2, inState )
         exState = GetExternalState( tn, iEdge, myDGSEM % nElems, inState, myDGSEM % nEq )
         nHat = ONE

         flux = RiemannSolver( inState, exState, netV, nHat, myDGSEM % nEq )

         ! Store the flux for the elements which share this edge
         CALL myDGSEM % SetBoundaryFluxAtBoundary( iEdge-1, 2, flux )

      ELSE ! this edge is a boundary edge

         CALL myDGSEM % GetBoundarySolutionAtBoundary( iEdge-1, 2, inState )
         CALL myDGSEM % GetBoundarySolutionAtBoundary( iEdge, 1, exState )
         
         nHat = ONE

         flux = RiemannSolver( inState, exState, netV, nHat, myDGSEM % nEq )

         ! Store the flux for the elements which share this edge
         CALL myDGSEM % SetBoundaryFluxAtBoundary( iEdge-1, 2, flux )
         CALL myDGSEM % SetBoundaryFluxAtBoundary( iEdge, 1, -flux )

      ENDIF ! Choosing the edge type
 
 END SUBROUTINE EdgeFlux_AlongShelf1D
!
!
!
 SUBROUTINE MappedTimeDerivative_AlongShelf1D( myDGSEM, iEl, tn ) 
 ! S/R MappedTimeDerivative_AlongShelf1D
 ! 
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS(AlongShelf1D), INTENT(inout) :: myDGSEM
   INTEGER, INTENT(in)                :: iEl
   REAL(prec), INTENT(in)             :: tn
   ! LOCAL
   REAL(prec) :: contFlux(0:myDGSEM % nS,1:myDGSEM % nEq)
   REAL(prec) :: contFluxDer(0:myDGSEM % nS,1:myDGSEM % nEq)
   REAL(prec) :: tend(0:myDGSEM % nS,1:myDGSEM % nEq)
   REAL(prec) :: dMatX(0:myDGSEM % nS,0:myDGSEM % nS)
   REAL(prec) :: qWeightX(0:myDGSEM % nS)
   REAL(prec) :: least(0:myDGSEM % nS)
   REAL(prec) :: lwest(0:myDGSEM % nS)
   REAL(prec) :: fL(1:myDGSEM % nEq), fR(1:myDGSEM % nEq), sol(1:myDGSEM % nEq)
   REAL(prec) :: netv(1:myDGSEM % nEq)
   REAL(prec) :: dx
   INTEGER    :: iS, nS, iEq, nEq

      CALL myDGSEM % dgStorage % GetNumberOfNodes( nS )
      CALL myDGSEM % dgStorage % GetDerivativeMatrix( dMatX )
      CALL myDGSEM % dgStorage % GetQuadratureWeights( qWeightX )
      CALL myDGSEM % dgStorage % GetEasternInterpolants( least )
      CALL myDGSEM % dgStorage % GetWesternInterpolants( lwest )
      nEq = myDGSEM % nEq
      dx = myDGSEM % dx

      contFluxDer = ZERO
      tend = ZERO

      DO iS = 0,nS ! Loop over the x-points

         ! Get the solution at x,y
         CALL myDGSEM % GetSolutionAtNode( iEl, iS, sol )
         netv(1:nEq) = myDGSEM % vbar - myDGSEM % shelfwaves % eigenvalues(1:nEq)
         contFlux(iS,:) = XFlux(  tn, netv, sol, nEq )

      ENDDO ! iS, Loop over the x-points

      ! Get the numerical fluxes at the boundaries (should be calculated before
      ! this routine is called )
        
      ! Get the "west" flux at iP
      CALL myDGSEM % GetBoundaryFluxAtBoundary( iEl, 1, fL )

      ! Get the "east" flux at iP
      CALL myDGSEM % GetBoundaryFluxAtBoundary( iEl, 2,  fR )

      ! At this y-level, calculate the DG-advective derivative

      DO iEq = 1, nEq
         contFluxDer(0:nS,iEq) = DGSystemDerivative( nS, dMatX, qWeightX, fL(iEq), fR(iEq), &
                                                     contFlux(:,iEq), lwest, least  )
      ENDDO

      tend(0:nS,1:nEq) = -TWO*contFluxDer(0:nS,1:nEq)/dx

      DO iEq = 1, nEq
         tend(0:nS,iEq) = tend(0:nS,iEq) + &
                          myDGSEM % shelfwaves % growthRates(iEq)*myDGSEM % shelfshaper(0:nS,iEl) 
      ENDDO

     CALL myDGSEM % SetTendency( iEl, tend )

 END SUBROUTINE MappedTimeDerivative_AlongShelf1D
!
!
!  
FUNCTION DGSystemDerivative(  nP, dMat, qWei, lFlux, rFlux, intFlux, lagLeft, lagRight  ) RESULT( tendency )
 ! FUNCTION DGSystemDerivative ( Discontinous Galerkin time Derivative)
 ! 
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   INTEGER    :: nP
   REAL(prec) :: dMat(0:nP,0:nP)
   REAL(prec) :: qWei(0:nP)
   REAL(prec) :: lFlux
   REAL(prec) :: rFlux
   REAL(prec) :: intFlux(0:nP)
   REAL(prec) :: lagLeft(0:nP), lagRight(0:nP)
   REAL(prec) :: tendency(0:nP)
   ! Local 
   INTEGER :: kP

      ! Apply the derivative operator
      tendency = MATMUL( dMat, intflux )

      DO kP = 0, nP ! Loop over the quadrature points
            tendency(kP) = tendency(kP) + ( rFlux*lagRight(kP) + lFlux*lagLeft(kP) )/qWei(kP)
      ENDDO ! kP, loop over the quadrature points
                                         
 END FUNCTION DGSystemDerivative
!
!
!
 FUNCTION RiemannSolver( inState, outState, netv, nHat, nEq ) RESULT( numFlux )
 ! FUNCTION RiemannSolver 
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   INTEGER    :: nEq
   REAL(prec) :: inState(1:nEq)
   REAL(prec) :: outState(1:nEq)
   REAL(prec) :: netv(1:nEq)
   REAL(prec) :: nHat
   REAL(prec) :: numFlux(1:nEq)
   
   ! LOCAL
   REAL(prec) :: fIn(1:nEq), fOut(1:nEq)
   REAL(prec) :: jump(1:nEq)
   REAL(prec) :: fac(1:nEq)


      jump = outState - inState
      fIn  = XFlux( ZERO, netv, inState, nEq )*nHat
      fOut = XFlux( ZERO, netv, outState, nEq )*nHat
      fac  = abs(netv)

      numFlux = HALF*( fIn + fOut - fac*jump )

 END FUNCTION RiemannSolver
!
!
! 
 FUNCTION XFlux( tn, netv, solAtX, nEq ) RESULT( fx )
 ! FUNCTION XFlux
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE 
   INTEGER    :: nEq
   REAL(prec) :: tn
   REAL(prec) :: netv(1:nEq)
   REAL(prec) :: solAtX(1:nEq)
   REAL(prec) :: fx(1:nEq)

      fx = netv*solAtX

 END FUNCTION XFlux
!
!
!                 
 FUNCTION GetExternalState( t, edgeID, nElems, intState, nEq ) RESULT( extState )
 ! S/R GetExternalState
 ! 
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   INTEGER    :: nEq
   REAL(prec) :: intState(1:nEq)
   REAL(prec) :: t
   INTEGER    :: edgeID, nElems
   REAL(prec) :: extState(1:nEq)

       extState = ZERO 

 END FUNCTION GetExternalState
!
!
!==================================================================================================!
!-------------------------------------- FILE I/O ROUTINES -----------------------------------------!
!==================================================================================================!
!
!
 SUBROUTINE CoarseToFine_AlongShelf1D( myDGSEM, iElx, iEly, x, y, psi, psimodes, bsf )
 ! CoarseToFine
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS( AlongShelf1D ), INTENT(in) :: myDGSEM
   INTEGER, INTENT(in)               :: iElx, iEly
   REAL(prec), INTENT(out)           :: x(0:myDGSEM % nPlot)
   REAL(prec), INTENT(out)           :: y(0:myDGSEM % nPlot)
   REAL(prec), INTENT(out)           :: psi(0:myDGSEM % nPlot,0:myDGSEM % nPlot)
   REAL(prec), INTENT(out)           :: psimodes(0:myDGSEM % nPlot,0:myDGSEM % nPlot,&
                                                 1:myDGSEM % shelfwaves % params % nModes2Write)
   REAL(prec), INTENT(out)           :: bsf(0:myDGSEM % nPlot)
   ! Local
   REAL(prec) :: localX(0:myDGSEM % nS), localY(0:myDGSEM % nS)
   REAL(prec) :: phi(0:myDGSEM % nS,1:myDGSEM % nEq)
   REAL(prec) :: phiPL(0:myDGSEM % nPlot)
   REAL(prec) :: F(0:myDGSEM % nPlot)
   INTEGER    :: iEq, iS, iP, i


      CALL myDGSEM % GetSolution( iEly, phi )
      localX = myDGSEM % shelfWaves % mesh % GetPositions( iElx )
      localY = myDGSEM % x(:,iEly)
      
      CALL myDGSEM % shelfwaves % cgStorage % interp % CoarseToFine( localX, &
                                                                     myDGSEM % shelfwaves % plMatS, &
                                                                     myDGSEM % nPlot, &
                                                                     x )
      
      CALL myDGSEM % dgStorage % interp % CoarseToFine( localY, &
                                                        myDGSEM % plMatS, &
                                                        myDGSEM % nPlot, &
                                                        y )

      CALL myDGSEM % shelfwaves % cgStorage % interp % CoarseToFine( myDGSEM % shelfwaves % backgroundStreamF(:,iElx), &
                                                                     myDGSEM % plMatS, &
                                                                     myDGSEM % nPlot, &
                                                                     bsf )
    
      psi = ZERO
      DO iEq = 1, myDGSEM % nEq

         CALL myDGSEM % dgStorage % interp % CoarseToFine( phi(:,iEq), &
                                                           myDGSEM % plMatS, &
                                                           myDGSEM % nPlot, &
                                                           phiPL )

         CALL myDGSEM % shelfwaves % cgStorage % interp % &
                        CoarseToFine( myDGSEM % shelfwaves % eigenfunctions(:,iElx,iEq), &
                                      myDGSEM % shelfwaves % plMatS, &
                                      myDGSEM % nPlot, &
                                      F )
         
         DO iP = 0, myDGSEM % nPlot
            DO iS = 0, myDGSEM % nPlot
               psi(iS,iP) = psi(iS,iP) + phiPL(iP)*F(iS)
            ENDDO
         ENDDO

      ENDDO

      DO iEq = 1, myDGSEM % shelfwaves % params % nModes2Write

         i = myDGSEM % shelfwaves % params % modes2Write(iEq)
         CALL myDGSEM % dgStorage % interp % CoarseToFine( phi(:,i), &
                                                           myDGSEM % plMatS, &
                                                           myDGSEM % nPlot, &
                                                           phiPL )

         CALL myDGSEM % shelfwaves % cgStorage % interp % &
                        CoarseToFine( myDGSEM % shelfwaves % eigenfunctions(:,iElx,i), &
                                      myDGSEM % shelfwaves % plMatS, &
                                      myDGSEM % nPlot, &
                                      F )
         DO iP = 0, myDGSEM % nPlot
            DO iS = 0, myDGSEM % nPlot
               psimodes(iS,iP,iEq) = phiPL(iP)*F(iS)
            ENDDO
         ENDDO

      ENDDO

      
 END SUBROUTINE CoarseToFine_AlongShelf1D
!
!
!
 SUBROUTINE WriteTecplot_AlongShelf1D( myDGSEM, filename )
 ! S/R WriteTecplot
 !  
 ! =============================================================================================== !
 ! DECLARATIONS
  IMPLICIT NONE
  CLASS( AlongShelf1D ), INTENT(in)  :: myDGsem
  CHARACTER(*), INTENT(in), OPTIONAL :: filename
  !LOCAL
  REAL(prec)  :: x(0:myDGSEM % nPlot), y(0:myDGSEM % nPlot)
  REAL(prec)  :: psi(0:myDGSEM % nPlot,0:myDGSEM % nPlot)
  REAL(prec)  :: psims(0:myDGSEM % nPlot,0:myDGSEM % nPlot,&
                       1:myDGSEM % shelfwaves % params % nModes2Write)
  REAL(prec)  :: bsf(0:myDGSEM % nPlot)
  REAL(prec)  :: h, hp
  INTEGER     :: iS, iP, iEly, iElx, fUnit, nPlot, i
  CHARACTER(5) :: zoneID
  CHARACTER(2) :: modeID
  CHARACTER(200) :: modeTags

      nPlot = myDGSEM % nPlot
    
    
      IF( PRESENT(filename) )THEN
         OPEN( UNIT=NEWUNIT(fUnit), &
               FILE= TRIM(filename)//'.tec', &
               FORM='formatted', &
               STATUS='replace')
      ELSE
         OPEN( UNIT=NEWUNIT(fUnit), &
               FILE= 'StreamFunction.tec', &
               FORM='formatted', &
               STATUS='replace')  
      ENDIF
    
      modeTags = ''
      DO i = 1, myDGSEM % shelfwaves % params % nModes2Write
         WRITE(modeID,'(I2.2)') myDGSEM % shelfwaves % params % modes2Write(i)
         modeTags = TRIM(modeTags)//', "Psi'//modeID//'"'
      ENDDO
     !PRINT*, TRIM(modeTags)
      WRITE(fUnit,*) 'VARIABLES = "X", "Y", "h", "PsiBar", "Psi"'//TRIM(modeTags)
  
      DO iEly = 1, myDGSEM % nElems

         DO iElx = 1, myDGSEM % shelfwaves % nElems/2
  
            CALL myDGSEM % CoarseToFine( iElx, iEly, x, y, psi, psims, bsf )
            WRITE(zoneID,'(I5.5)') iElx + (iEly - 1)*( myDGSEM % shelfwaves % nElems )
            WRITE(fUnit,*) 'ZONE T="el'//trim(zoneID)//'", I=',nPlot+1,', J=', nPlot+1,',F=POINT'

            DO iP = 0, nPlot
               DO iS = 0, nPlot
                  h = myDGSEM % shelfwaves % ShelfProfile( x(iS) )
                  hp = myDGSEM % shelfwaves % ShelfMod( x(iS) )*&
                       myDGSEM % ShelfSteepeningShape( y(iP) )
                  WRITE (fUnit,*)  x(iS), y(iP), h + 0.5_prec*hp, bsf(iS), psi(iS,iP), psims(iS,iP,:)
               ENDDO
            ENDDO

         ENDDO
    ENDDO

    CLOSE(UNIT=fUnit)

 END SUBROUTINE WriteTecplot_AlongShelf1D
!
!
!
 SUBROUTINE WritePickup_AlongShelf1D( myDGSEM, iter )
 ! S/R WritePickup
 ! 
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS( AlongShelf1D ), INTENT(in) :: myDGSEM
   INTEGER, INTENT(in)               :: iter
  ! LOCAL
   REAL(prec)    :: sol(0:myDGSEM % nS,1:myDGSEM % nEq)
   CHARACTER(10) :: iterChar
   INTEGER       :: iEl
   INTEGER       :: thisRec, fUnit
   INTEGER       :: nS

     nS = myDGSEM % nS
     
      WRITE(iterChar,'(I10.10)') iter

      OPEN( UNIT=NEWUNIT(fUnit), &
            FILE='AlongShelf1D.'//iterChar//'.pickup', &
            FORM='unformatted',&
            ACCESS='direct',&
            STATUS='replace',&
            ACTION='WRITE',&
            CONVERT='big_endian',&
            RECL=prec*(nS+1)*(myDGSEM % nEq) )

      thisRec = 1 
      DO iEl = 1, myDGSEM % nElems
        
         CALL myDGSEM % GetSolution( iEl, sol )
         WRITE( fUnit, REC=thisRec )sol
         thisRec = thisRec+1

     ENDDO

     CLOSE(UNIT=fUnit)

 END SUBROUTINE WritePickup_AlongShelf1D
!
!
!
  SUBROUTINE ReadPickup_AlongShelf1D( myDGSEM, iter )
 ! S/R ReadPickup
 ! 
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS( AlongShelf1D ), INTENT(inout) :: myDGSEM
   INTEGER, INTENT(in)                  :: iter
  ! LOCAL
   REAL(prec)    :: sol(0:myDGSEM % nS,1:myDGSEM % nEq)
   CHARACTER(10) :: iterChar
   INTEGER       :: iEl
   INTEGER       :: thisRec, fUnit
   INTEGER       :: nS

      nS = myDGSEM % nS
     
      WRITE(iterChar,'(I10.10)') iter

      OPEN( UNIT=NEWUNIT(fUnit), &
            FILE='AlongShelf1D.'//iterChar//'.pickup', &
            FORM='unformatted',&
            ACCESS='direct',&
            STATUS='old',&
            ACTION='READ',&
            CONVERT='big_endian',&
            RECL=prec*(nS+1)*(myDGSEM % nEq) )
     
      thisRec = 1
      DO iEl = 1, myDGSEM % nElems
        
         READ( fUnit, REC=thisRec )sol 
         thisRec = thisRec+1
         CALL myDGSEM % SetSolution( iEl, sol )
        
      ENDDO

      CLOSE(UNIT=fUnit)

 END SUBROUTINE ReadPickup_AlongShelf1D
!
!
! 
 END MODULE AlongShelfClass



