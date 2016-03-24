! BurgersClass.f90 ( new with v2.1 - 21 March 2016)
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
!
! ========================================= Logs ================================================= !
! yyyy-mm-dd  Joe  <joe@clay>
!
!   
! 
! //////////////////////////////////////////////////////////////////////////////////////////////// !


MODULE BurgersClass
! BurgersClass.f90
!
! schoonover.numerics@gmail.com
!
! o (ver 2.1) March 2016
!
!
! 
!  
! ================================================================================================ !

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
!
USE RollOffFilter1D_Class
! src/highend/burgers1D/
USE BurgersParamsClass
! Nocturnal Aviation classes and extensions
!USE FTTimerClass
!USE TIMING


! Parallel Libraries
USE OMP_LIB


IMPLICIT NONE
!
! Some parameters that are specific to this module


! ================================================================================================ !
!     The Burgers Class is provides an example of the 1-D DGSEM with an adaptive filter applied
! ================================================================================================ !  



    TYPE Burgers1D
      INTEGER                               :: nPlot, nS, nElems
      REAL(prec)                            :: dx
      REAL(prec), ALLOCATABLE               :: x(:,:)

      TYPE( RollOffFilter1D )               :: modalFilter
      REAL(prec), ALLOCATABLE               :: E1(:), E2(:), lim(:,:)
      REAL(prec), ALLOCATABLE               :: dE1(:), dE2(:)
      REAL(prec)                            :: filtFac

      TYPE( NodalStorage_1D )               :: dGStorage
      TYPE( DGSEMSolution_1D ), ALLOCATABLE :: sol(:)
      TYPE( BurgersParams )                 :: params
      REAL(prec), ALLOCATABLE               :: plMatS(:,:)

      CONTAINS

      PROCEDURE :: Build => Build_Burgers1D
      PROCEDURE :: Trash => Trash_Burgers1D
      PROCEDURE :: BuildMesh => BuildMesh_Burgers1D
      
      ! DGSEMSolution_2DStorage_Class Wrapper Routines
      PROCEDURE :: GetSolution => GetSolution_Burgers1D
      PROCEDURE :: SetSolution => SetSolution_Burgers1D
      PROCEDURE :: GetSolutionAtNode => GetSolutionAtNode_Burgers1D
      PROCEDURE :: SetSolutionAtNode => SetSolutionAtNode_Burgers1D
      PROCEDURE :: GetTendency => GetTendency_Burgers1D
      PROCEDURE :: SetTendency => SetTendency_Burgers1D
      PROCEDURE :: GetTendencyAtNode => GetTendencyAtNode_Burgers1D
      PROCEDURE :: SetTendencyAtNode => SetTendencyAtNode_Burgers1D
      PROCEDURE :: GetBoundarySolutionAtBoundary => GetBoundarySolutionAtBoundary_Burgers1D
      PROCEDURE :: SetBoundarySolutionAtBoundary => SetBoundarySolutionAtBoundary_Burgers1D
      PROCEDURE :: GetBoundaryFluxAtBoundary => GetBoundaryFluxAtBoundary_Burgers1D
      PROCEDURE :: SetBoundaryFluxAtBoundary => SetBoundaryFluxAtBoundary_Burgers1D
      
       ! Type Specific Routines
      PROCEDURE :: GlobalTimeDerivative => GlobalTimeDerivative_Burgers1D
      PROCEDURE :: ForwardStepRK3 => ForwardStepRK3_Burgers1D
      PROCEDURE :: EdgeFlux => EdgeFlux_Burgers1D
      PROCEDURE :: MappedTimeDerivative => MappedTimeDerivative_Burgers1D 
      PROCEDURE :: CalculateSolutionAtBoundaries => CalculateSolutionAtBoundaries_Burgers1D
      PROCEDURE :: DoTheAdaptiveFiltering => DoTheAdaptiveFiltering_Burgers1D
      
      PROCEDURE :: CoarseToFine => CoarseToFine_Burgers1D
      PROCEDURE :: WriteTecplot => WriteTecplot_Burgers1D
      PROCEDURE :: WritePickup => WritePickup_Burgers1D
      PROCEDURE :: ReadPickup => ReadPickup_Burgers1D

    END TYPE Burgers1D


 INTEGER, PARAMETER, PRIVATE :: nBurgerEq = 1
 CONTAINS
!
!
!==================================================================================================!
!------------------------------- Manual Constructors/Destructors ----------------------------------!
!==================================================================================================!
!
!
 SUBROUTINE Build_Burgers1D( myDGSEM, forceNoPickup )
 ! S/R Build
 ! 
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS(Burgers1D), INTENT(inout) :: myDGSEM
   INTEGER, INTENT(in), OPTIONAL   :: forceNoPickup
   !LOCAL
   INTEGER :: iEl, nS, nPlot, nElems
   REAL(prec)              :: dx, nR, nC
   REAL(prec), ALLOCATABLE ::sNew(:)

      CALL myDGSEM % params % Build( )
     
      
      nS     = myDGSEM % params % polyDeg
      nElems = myDGSEM % params % nElems
      dx     = myDGSEM % params % dx
 
      myDGSEM % nElems  = nElems
      myDGSEM % dx      = dx
      myDGSEM % nS      = nS
      myDGSEM % nPlot   = myDGSEM % params % nPlot
      
      nPlot = myDGSEM % params % nPlot
      
      ALLOCATE( sNew(0:nPlot) )
      CALL myDGSEM % dGStorage % Build( nS, GAUSS, DG )
 
      CALL myDGSEM % modalFilter % Build( myDGSEM % dgStorage, &
                                          myDGSEM % params % nCutoff )
      nC = real( myDGSEM % params % nCutoff, prec )
      nR = real( myDGSEM % nS, prec )
      myDGSEM % filtFac = ( (nC + ONE)**(-2/3) - TWO**(-2/3) )/( (nR+ONE)**(-2/3) - (nC+ONE)**(-2/3) )

      CALL myDGSEM % BuildMesh(  )
      
      ! Space is allocated for the solution-storage structure and an blank structure is built
   
      ALLOCATE( myDGSEM % sol(1:nElems) )
      ALLOCATE( myDGSEM % E1(1:nElems), myDGSEM % E2(1:nElems), myDGSEM % lim(1:nElems,1:3) )
      ALLOCATE( myDGSEM % dE1(1:nElems), myDGSEM % dE2(1:nElems) )

      ! Build and initialize the solution storage structure to zero
      DO iEl = 1, myDGSEM % nElems
         CALL myDGSEM % sol(iEl) % Build( nS, nBurgerEq )
      ENDDO
      
      ! In the event that this is a pickup run, we'll READ the pickup file here for the solution 
      ! and the relaxation fields. A negative initial iterate can be used to specify to start 
      ! from zeros.
      
      IF( .NOT. PRESENT(forceNoPickup) )then
         
         ! This call reads the solution from a pickup file if desired
         CALL myDGSEM % ReadPickup( myDGSEM % params % iterInit )

      ENDIF
      
      
      nPlot = myDGSEM % params % nPlot
      myDGSEM % nPlot = nPlot
      
      sNew = UniformPoints( -ONE, ONE, nPlot )
      ALLOCATE( myDGSEM % plMatS(0:nPlot,0:nS) )
      ! Build the plotting matrix
      CALL myDGSEM % dgStorage % interp % CalculateInterpolationMatrix( nPlot, sNew, myDGSEM % plMatS )
                                                                        
      DEALLOCATE(sNew)
      
 END SUBROUTINE Build_Burgers1D
!
!
!
 SUBROUTINE Trash_Burgers1D( myDGSEM )
 ! S/R Trash
 ! 
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS(Burgers1D), INTENT(inout) :: myDGSEM
   ! LOCAL
   INTEGER :: iEl

     DO iEl = 1, myDGSEM % nElems
        CALL myDGSEM % sol(iEl) % Trash( )
     ENDDO

     CALL myDGSEM % dGStorage % Trash( )
     CALL myDGSEM % modalFilter % Trash( )

     DEALLOCATE( myDGSEM % E1, myDGSEM % E2, myDGSEM % dE1, myDGSEM % dE2 )
     DEALLOCATE( myDGSEM % sol ) 
     DEALLOCATE( myDGSEM % plMatS )

 END SUBROUTINE Trash_Burgers1D
!
!
!
 SUBROUTINE BuildMesh_Burgers1D( myDGSEM )
 ! S/R BuildMesh
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS( Burgers1D ), INTENT(inout) :: myDGSEM
   ! LOCAL
   INTEGER :: iEl, nEl, i, nS
   REAL(prec) :: dx, x1, x2
   REAL(prec) :: s(0:myDGSEM % nS)
   
      PRINT*,'Module Burgers1DClass.f90 : S/R BuildQuadMesh :'
      PRINT*,' Loading default mesh.'
      
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



 END SUBROUTINE BuildMesh_Burgers1D
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
 SUBROUTINE GetSolution_Burgers1D( myDGSEM, iEl, theSolution  )
 ! S/R GetSolution
 !  
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS(Burgers1D), INTENT(in) :: myDGSEM
   INTEGER, INTENT(in)          :: iEl
   REAL(prec), INTENT(out)      :: theSolution(0:myDGSEM % nS, 1:nBurgerEq)

      CALL myDGSEM % sol(iEl) % GetSolution( theSolution )

 END SUBROUTINE GetSolution_Burgers1D
!
!
!
 SUBROUTINE SetSolution_Burgers1D( myDGSEM, iEl, theSolution  )
 ! S/R SetSolution
 !  
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS(Burgers1D), INTENT(inout) :: myDGSEM
   INTEGER, INTENT(in)             :: iEl
   REAL(prec), INTENT(in)          :: theSolution(0:myDGSEM % nS, 1:nBurgerEq)

      CALL myDGSEM % sol(iEl) % SetSolution( theSolution )

 END SUBROUTINE SetSolution_Burgers1D
!
!
!
 SUBROUTINE GetSolutionAtNode_Burgers1D( myDGSEM, iEl, i, theSolution  )
 ! S/R GetSolution
 !  
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS(Burgers1D), INTENT(in) :: myDGSEM
   INTEGER, INTENT(in)          :: iEl
   INTEGER, INTENT(in)          :: i
   REAL(prec), INTENT(out)      :: theSolution(1:nBurgerEq)

      CALL myDGSEM % sol(iEl) % GetSolutionAtNode( i, theSolution )

 END SUBROUTINE GetSolutionAtNode_Burgers1D
!
!
!
 SUBROUTINE SetSolutionAtNode_Burgers1D( myDGSEM, iEl, i, theSolution  )
 ! S/R SetSolution
 !  
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS(Burgers1D), INTENT(inout) :: myDGSEM
   INTEGER, INTENT(in)             :: iEl
   INTEGER, INTENT(in)             :: i
   REAL(prec), INTENT(in)          :: theSolution(1:nBurgerEq)

      CALL myDGSEM % sol(iEl) % SetSolutionAtNode( i, theSolution )

 END SUBROUTINE SetSolutionAtNode_Burgers1D
!
! ----------------------------------------- Tendency --------------------------------------------- !
!
 SUBROUTINE GetTendency_Burgers1D( myDGSEM, iEl, theTend  )
 ! S/R GetTendency
 !  
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS(Burgers1D), INTENT(in) :: myDGSEM
   INTEGER, INTENT(in)          :: iEl
   REAL(prec), INTENT(out)      :: theTend(0:myDGSEM % nS, 1:nBurgerEq)

      CALL myDGSEM % sol(iEl) % GetTendency( theTend )

 END SUBROUTINE GetTendency_Burgers1D
!
!
!
 SUBROUTINE SetTendency_Burgers1D( myDGSEM, iEl, theTend  )
 ! S/R SetTendency
 !  
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS(Burgers1D), INTENT(inout) :: myDGSEM
   INTEGER, INTENT(in)             :: iEl
   REAL(prec), INTENT(in)          :: theTend(0:myDGSEM % nS, 1:nBurgerEq)

      CALL myDGSEM % sol(iEl) % SetTendency( theTend )   

 END SUBROUTINE SetTendency_Burgers1D
!
!
!
 SUBROUTINE GetTendencyAtNode_Burgers1D( myDGSEM, iEl, i, theTend  )
 ! S/R GetTendencyAtNode
 !  
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS(Burgers1D), INTENT(in) :: myDGSEM
   INTEGER, INTENT(in)          :: iEl
   INTEGER, INTENT(in)          :: i
   REAL(prec), INTENT(out)      :: theTend(1:nBurgerEq)

      CALL myDGSEM % sol(iEl) % GetTendencyAtNode( i, theTend )

 END SUBROUTINE GetTendencyAtNode_Burgers1D
!
!
!
 SUBROUTINE SetTendencyAtNode_Burgers1D( myDGSEM, iEl, i,  theTend  )
 ! S/R SetTendencyAtNode
 !  
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS(Burgers1D), INTENT(inout) :: myDGSEM
   INTEGER, INTENT(in)             :: iEl
   INTEGER, INTENT(in)             :: i
   REAL(prec), INTENT(in)          :: theTend(1:nBurgerEq)

      CALL myDGSEM % sol(iEl) % SetTendencyAtNode( i, theTend )

 END SUBROUTINE SetTendencyAtNode_Burgers1D
!
! ------------------------------- Boundary Solution ---------------------------------------------- !
!
 SUBROUTINE GetBoundarySolutionAtBoundary_Burgers1D( myDGSEM, iEl, boundary, theSolution  )
 ! S/R GetBoundarySolution
 !  
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS(Burgers1D), INTENT(in) :: myDGSEM
   INTEGER, INTENT(in)          :: iEl
   INTEGER, INTENT(in)          :: boundary
   REAL(prec), INTENT(out)      :: theSolution(1:nBurgerEq)

      CALL myDGSEM % sol(iEl) % GetBoundarySolutionAtBoundary( boundary, theSolution )

 END SUBROUTINE GetBoundarySolutionAtBoundary_Burgers1D
!
!
!
 SUBROUTINE SetBoundarySolutionAtBoundary_Burgers1D( myDGSEM, iEl, boundary, theSolution  )
 ! S/R SetBoundarySolution
 !  
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS(Burgers1D), INTENT(inout) :: myDGSEM
   INTEGER, INTENT(in)             :: iEl
   INTEGER, INTENT(in)             :: boundary
   REAL(prec), INTENT(in)          :: theSolution(1:nBurgerEq)

      CALL myDGSEM % sol(iEl) % SetBoundarySolutionAtBoundary( boundary, theSolution )

 END SUBROUTINE SetBoundarySolutionAtBoundary_Burgers1D
!
! --------------------------------- Boundary Flux ------------------------------------------------ !
!
 SUBROUTINE GetBoundaryFluxAtBoundary_Burgers1D( myDGSEM, iEl, boundary, theFlux  )
 ! S/R GetBoundaryFlux
 !  
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS(Burgers1D), INTENT(in) :: myDGSEM
   INTEGER, INTENT(in)          :: iEl
   INTEGER, INTENT(in)          :: boundary
   REAL(prec), INTENT(out)      :: theFlux(1:nBurgerEq)

      CALL myDGSEM % sol(iEl) % GetBoundaryFluxAtBoundary( boundary, theFlux )

 END SUBROUTINE GetBoundaryFluxAtBoundary_Burgers1D
!
!
!
 SUBROUTINE SetBoundaryFluxAtBoundary_Burgers1D( myDGSEM, iEl, boundary, theFlux  )
 ! S/R SetBoundaryFlux
 !  
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS(Burgers1D), INTENT(inout) :: myDGSEM
   INTEGER, INTENT(in)             :: iEl
   INTEGER, INTENT(in)             :: boundary
   REAL(prec), INTENT(in)          :: theFlux(1:nBurgerEq)

      CALL myDGSEM % sol(iEl) % SetBoundaryFluxAtBoundary( boundary, theFlux )

 END SUBROUTINE SetBoundaryFluxAtBoundary_Burgers1D
!
!
!==================================================================================================!
!--------------------------------- Type Specific Routines -----------------------------------------!
!==================================================================================================!
!
!
 SUBROUTINE CalculateSolutionAtBoundaries_Burgers1D( myDGSEM, iEl )
 ! S/R CalculateSolutionAtBoundaries
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS(Burgers1D), INTENT(inout) :: myDGSEM
   INTEGER, INTENT(in)             :: iEl
   
      CALL myDGSEM % sol(iEl) % CalculateSolutionAtBoundaries( myDGSEM % dgStorage )
      
 END SUBROUTINE CalculateSolutionAtBoundaries_Burgers1D
!
!
!
 SUBROUTINE ForwardStepRK3_Burgers1D( myDGSEM, tn )
 ! S/R ForwardStepRK3( 3rd order Runge-Kutta)
 ! 
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS(Burgers1D), INTENT(inout) :: myDGSEM
   REAL(prec), INTENT(in)          :: tn
   ! LOCAL
   REAL(prec) :: t, dt
   REAL(prec) :: G2D(0:myDGSEM % nS, 1:nBurgerEq, 1:myDGSEM %  nElems) 
   REAL(prec) :: dSdt(0:myDGSEM % nS, 1:nBurgerEq )
   INTEGER    :: m, iEl


     dt = myDGSEM % params % dt
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
         

         !$OMP PARALLEL
         !$OMP DO
         DO iEl = 1, myDGSEM % nElems
            CALL myDGSEM % DoTheAdaptiveFiltering( iEl, m ) 
         ENDDO 
         !$OMP END DO 
         !$OMP FLUSH( myDGSEM )
         !$OMP END PARALLEL
         
      ENDDO ! m, loop over the RK3 steps
   
   
       
 END SUBROUTINE ForwardStepRK3_Burgers1D
!
!
!
 SUBROUTINE GlobalTimeDerivative_Burgers1D( myDGSEM, tn ) 
 ! S/R GlobalTimeDerivative_Burgers1D
 ! 
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS(Burgers1D), INTENT(inout) :: myDGSEM
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

 
 END SUBROUTINE GlobalTimeDerivative_Burgers1D
!
!
! 
 SUBROUTINE DoTheAdaptiveFiltering_Burgers1D( myDGSEM, iEl, m )
 !
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS(Burgers1D), INTENT(inout) :: myDGSEM
   INTEGER, INTENT(inout)          :: iEl, m
   ! Local
   REAL(prec) :: sol(0:myDGSEM % nS, 1:nBurgerEq)
   REAL(prec) :: solBar(1:nBurgerEq)
   REAL(prec) :: solFiltered(0:myDGSEM % nS, 1:nBurgerEq)
   REAL(prec) :: qWeightX(0:myDGSEM % nS), dx, xi
   REAL(prec) :: E1prior, E2prior, E1, E2, dE2, dE1
   INTEGER    :: iS


      CALL myDGSEM % GetSolution( iEl, sol )
      CALL myDGSEM % dgStorage % GetQuadratureWeights( qWeightX )

      dx = myDGSEM % dx

      ! The first step in the process is to apply the filter to obtain the well resolved solution
      solFiltered(:,1) = myDGSEM % modalFilter % ApplyFilter( sol(:,1) )
      ! This filtered solution is subtracted from the full solution to obtain the marginally-resolved
      ! portion of the solution
      sol = sol - solFiltered ! sol now contains the under-resolved component

      ! Now that we have the two distinct components of the Legendre spectra, we want to calculate
      ! the energy in each component, and the change in the energy of each component from the 
      ! previous model state. 
      ! If the small scale (marginally resolved) exhibits a growth in energy, this should be balanced
      ! by a decay in the small scale energy. Aliasing errors may cause unphysical growth in the 
      ! energy associated with the marginally resolved. In this case, the solution is assigned to the
      ! filtered solution, effectively implying dissipation.
      E1prior = myDGSEM % E1(iEl)
      E2prior = myDGSEM % E2(iEl)

      E1 = ZERO
      E2 = ZERO

      solBar = ZERO
      DO iS = 0, myDGSEM % nS
         solBar = solBar + solFiltered(iS,:)*qWeightX(iS)
      ENDDO
      solBar = solBar/SUM(qWeightX)
    
      DO iS = 0, myDGSEM % nS
         E1 = E1 + (solFiltered(iS,1)-solBar(1))**2*qWeightX(iS)*HALF*dx
         E2 = E2 + (sol(iS,1))**2*qWeightX(iS)*HALF*dx
      ENDDO

      myDGSEM % E1(iEl) = E1 
      myDGSEM % E2(iEl) = E2
      
      dE1 = E1-E1prior
      dE2 = E2-E2prior

      myDGSEM % dE1(iEl) = dE1
      myDGSEM % dE2(iEl) = dE2

      xi = (E2/E1)

      IF( dE2 > ZERO .AND. abs(dE1)/dE2 > ONE )THEN ! The energy in the small scales is growing faster than the large scale is giving it up
                            
         myDGSEM % lim(iEl,m) = xi
   
      ELSEIF( dE2 < ZERO )THEN

         myDGSEM % lim(iEl,m) = xi

      ELSE
         
         IF( m > 1)THEN
            myDGSEM % lim(iEl,m) = myDGSEM % lim(iEl,m-1)
         ENDIF

      ENDIF


      IF( xi > 1.05_prec*myDGSEM % lim(iEl,m) )THEN
         CALL myDGSEM % SetSolution( iEl, solFiltered )
         PRINT*, 'Filtered!'
      ENDIF

 END SUBROUTINE DoTheAdaptiveFiltering_Burgers1D
!
!
!
 SUBROUTINE EdgeFlux_Burgers1D( myDGSEM, iEdge, tn )
 ! S/R EdgeFlux
 ! 
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS( Burgers1D ), INTENT(inout) :: myDGSEM
   INTEGER, INTENT(in)               :: iEdge  
   REAL(prec), INTENT(in)            :: tn
   ! Local
   REAL(prec) :: flux(1:nBurgerEq)
   REAL(prec) :: inState(1:nBurgerEq)
   REAL(prec) :: exState(1:nBurgerEq)
   REAL(prec) :: nHat


      IF( iEdge == 1 )THEN! this is the left-most boundary "edge"
      
         CALL myDGSEM % GetBoundarySolutionAtBoundary( iEdge, 1, inState )
         exState = GetExternalState( tn, iEdge, myDGSEM % nElems, inState )
         nHat = -ONE

         flux = RiemannSolver( inState, exState, nHat )

         ! Store the flux for the elements which share this edge
         CALL myDGSEM % SetBoundaryFluxAtBoundary( iEdge, 1, flux )

      ELSEIF( iEdge == myDGSEM % nElems + 1 )THEN  ! This is the right most boundary edge
 
         CALL myDGSEM % GetBoundarySolutionAtBoundary( iEdge-1, 2, inState )
         exState = GetExternalState( tn, iEdge, myDGSEM % nElems, inState )
         nHat = ONE

         flux = RiemannSolver( inState, exState, nHat )

         ! Store the flux for the elements which share this edge
         CALL myDGSEM % SetBoundaryFluxAtBoundary( iEdge-1, 2, flux )

      ELSE ! this edge is a boundary edge

         CALL myDGSEM % GetBoundarySolutionAtBoundary( iEdge-1, 2, inState )
         CALL myDGSEM % GetBoundarySolutionAtBoundary( iEdge, 1, exState )
         
         nHat = ONE

         flux = RiemannSolver( inState, exState, nHat )

         ! Store the flux for the elements which share this edge
         CALL myDGSEM % SetBoundaryFluxAtBoundary( iEdge-1, 2, flux )
         CALL myDGSEM % SetBoundaryFluxAtBoundary( iEdge, 1, -flux )

      ENDIF ! Choosing the edge type
 
 END SUBROUTINE EdgeFlux_Burgers1D
!
!
!
 SUBROUTINE MappedTimeDerivative_Burgers1D( myDGSEM, iEl, tn ) 
 ! S/R MappedTimeDerivative_Burgers1D
 ! 
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS(Burgers1D), INTENT(inout) :: myDGSEM
   INTEGER, INTENT(in)             :: iEl
   REAL(prec), INTENT(in)          :: tn
   ! LOCAL
   REAL(prec) :: contFlux(0:myDGSEM % nS,1:nBurgerEq)
   REAL(prec) :: contFluxDer(0:myDGSEM % nS,1:nBurgerEq)
   REAL(prec) :: tend(0:myDGSEM % nS,1:nBurgerEq)
   REAL(prec) :: dMatX(0:myDGSEM % nS,0:myDGSEM % nS)
   REAL(prec) :: qWeightX(0:myDGSEM % nS)
   REAL(prec) :: least(0:myDGSEM % nS)
   REAL(prec) :: lwest(0:myDGSEM % nS)
   REAL(prec) :: fL(1:nBurgerEq), fR(1:nBurgerEq), sol(1:nBurgerEq)
   REAL(prec) :: dx
   INTEGER    :: iS, nS

      CALL myDGSEM % dgStorage % GetNumberOfNodes( nS )
      CALL myDGSEM % dgStorage % GetDerivativeMatrix( dMatX )
      CALL myDGSEM % dgStorage % GetQuadratureWeights( qWeightX )
      CALL myDGSEM % dgStorage % GetEasternInterpolants( least )
      CALL myDGSEM % dgStorage % GetWesternInterpolants( lwest )

      dx = myDGSEM % dx

      DO iS = 0,nS ! Loop over the x-points

         ! Get the solution at x,y
         CALL myDGSEM % GetSolutionAtNode( iEl, iS, sol )
         contFlux(iS,:) = XFlux(  tn, sol )

      ENDDO ! iS, Loop over the x-points

      ! Get the numerical fluxes at the boundaries (should be calculated before
      ! this routine is called )
        
      ! Get the "west" flux at iP
      CALL myDGSEM % GetBoundaryFluxAtBoundary( iEl, 1, fL )

      ! Get the "east" flux at iP
      CALL myDGSEM % GetBoundaryFluxAtBoundary( iEl, 2,  fR )

      ! At this y-level, calculate the DG-advective derivative

      contFluxDer(0:nS,1) = DGSystemDerivative( nS, dMatX, qWeightX, fL(1), fR(1), &
                                                contFlux(:,1), lwest, least  )
     

      tend(0:nS,1) = -TWO*contFluxDer(0:nS,1)/dx

     CALL myDGSEM % SetTendency( iEl, tend )

 END SUBROUTINE MappedTimeDerivative_Burgers1D
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
 FUNCTION RiemannSolver( inState, outState, nHat ) RESULT( numFlux )
 ! FUNCTION RiemannSolver 
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   REAL(prec) :: inState(1:nBurgerEq)
   REAL(prec) :: outState(1:nBurgerEq)
   REAL(prec) :: nHat
   REAL(prec) :: numFlux(1:nBurgerEq)
   
   ! LOCAL
   REAL(prec) :: fIn(1:nBurgerEq), fOut(1:nBurgerEq)
   REAL(prec) :: jump(1:nBurgerEq)
   REAL(prec) :: fac(1:nBurgerEq)


      jump = outState - inState
      ! ** Starting off with one-way advection, to make sure the conversion to 1-D is ok. 22 March, 2016 ( 1:07 PM )
      fIn  = XFlux( ZERO, inState )*nHat
      fOut = XFlux( ZERO, outState )*nHat
      fac  = max( abs(inState), abs(outState) )

      numFlux = HALF*( fIn + fOut - fac*jump )

 END FUNCTION RiemannSolver
!
!
! 
 FUNCTION XFlux( tn, solAtX ) RESULT( fx )
 ! FUNCTION XFlux
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE 
   REAL(prec) :: tn
   REAL(prec) :: solAtX(1:nBurgerEq)
   REAL(prec) :: fx(1:nBurgerEq)

      fx = solAtX*solAtX*HALF

 END FUNCTION XFlux
!
!
!                 
 FUNCTION GetExternalState( t, edgeID, nElems, intState ) RESULT( extState )
 ! S/R GetExternalState
 ! 
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   REAL(prec) :: intState(1:nBurgerEq)
   REAL(prec) :: t
   INTEGER    :: edgeID, nElems
   REAL(prec) :: extState(1:nBurgerEq)

       extState = ZERO 

       IF( edgeID == 1 )THEN
          
          extState(1) = ONE 
             
       ELSEIF( edgeID == nElems + 1 )THEN

          extState(1) =  -ONE
           
       ENDIF

        
 END FUNCTION GetExternalState
!
!
!==================================================================================================!
!-------------------------------------- FILE I/O ROUTINES -----------------------------------------!
!==================================================================================================!
!
!
 SUBROUTINE CoarseToFine_Burgers1D( myDGSEM, iEl, x, u )
 ! CoarseToFine
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS( Burgers1D ), INTENT(in) :: myDGSEM
   INTEGER, INTENT(in)            :: iEl
   REAL(prec), INTENT(out)        :: x(0:myDGSEM % nPlot)
   REAL(prec), INTENT(out)        :: u(0:myDGSEM % nPlot)
   ! Local
   REAL(prec) :: localX(0:myDGSEM % nS)
   REAL(prec) :: sol(0:myDGSEM % nS,1:nBurgerEq)


      CALL myDGSEM % GetSolution( iEl, sol )
      localX = myDGSEM % x(:,iEl)
      
      CALL myDGSEM % dgStorage % interp % CoarseToFine( localX, &
                                                        myDGSEM % plMatS, &
                                                        myDGSEM % nPlot, &
                                                        x )
      
    
      CALL myDGSEM % dgStorage % interp % CoarseToFine( sol(:,1), &
                                                        myDGSEM % plMatS, &
                                                        myDGSEM % nPlot, &
                                                        u )
                                                        

      
 END SUBROUTINE CoarseToFine_Burgers1D
!
!
!
 SUBROUTINE WriteTecplot_Burgers1D( myDGSEM, filename )
 ! S/R WriteTecplot
 !  
 ! =============================================================================================== !
 ! DECLARATIONS
  IMPLICIT NONE
  CLASS( Burgers1D ), INTENT(in)  :: myDGsem
  CHARACTER(*), INTENT(in), OPTIONAL :: filename
  !LOCAL
  REAL(prec)  :: x(0:myDGSEM % nPlot)
  REAL(prec)  :: u(0:myDGSEM % nPlot)
  INTEGER     :: iS, iEl, fUnit, nPlot
  CHARACTER(len=5) :: zoneID

    nPlot = myDGSEM % nPlot
    
    IF( PRESENT(filename) )THEN
       OPEN( UNIT=NEWUNIT(fUnit), &
             FILE= TRIM(filename)//'.curve', &
             FORM='formatted', &
             STATUS='replace')
    ELSE
       OPEN( UNIT=NEWUNIT(fUnit), &
             FILE= 'Burgers1D.tec', &
             FORM='formatted', &
             STATUS='replace')  
    ENDIF
    
    WRITE(fUnit,*) '#solution'
 
    DO iEl = 1, myDGsem % nElems

      CALL myDGSEM % CoarseToFine( iEl, x, u )

      DO iS = 0, nPlot
         WRITE (fUnit,*)  x( iS ), u(iS)
      ENDDO
        
    ENDDO

    CLOSE(UNIT=fUnit)

 END SUBROUTINE WriteTecplot_Burgers1D
!
!
!
 SUBROUTINE WritePickup_Burgers1D( myDGSEM, iter )
 ! S/R WritePickup
 ! 
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS( Burgers1D ), INTENT(in) :: myDGSEM
   INTEGER, INTENT(in)               :: iter
  ! LOCAL
   REAL(prec)    :: sol(0:myDGSEM % nS,1:nBurgerEq)
   CHARACTER(10) :: iterChar
   INTEGER       :: iEl
   INTEGER       :: thisRec, fUnit
   INTEGER       :: nS

     nS = myDGSEM % nS
     
      WRITE(iterChar,'(I10.10)') iter

      OPEN( UNIT=NEWUNIT(fUnit), &
            FILE='Burgers1D.'//iterChar//'.pickup', &
            FORM='unformatted',&
            ACCESS='direct',&
            STATUS='replace',&
            ACTION='WRITE',&
            CONVERT='big_endian',&
            RECL=prec*(nS+1) )

      thisRec = 1 
      DO iEl = 1, myDGSEM % nElems
        
         CALL myDGSEM % GetSolution( iEl, sol )
         WRITE( fUnit, REC=thisRec )sol(:,1) 
         thisRec = thisRec+1

     ENDDO

     CLOSE(UNIT=fUnit)

 END SUBROUTINE WritePickup_Burgers1D
!
!
!
  SUBROUTINE ReadPickup_Burgers1D( myDGSEM, iter )
 ! S/R ReadPickup
 ! 
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS( Burgers1D ), INTENT(inout) :: myDGSEM
   INTEGER, INTENT(in)                  :: iter
  ! LOCAL
   REAL(prec)    :: sol(0:myDGSEM % nS,1:nBurgerEq)
   CHARACTER(10) :: iterChar
   INTEGER       :: iEl
   INTEGER       :: thisRec, fUnit
   INTEGER       :: nS

      nS = myDGSEM % nS
     
      WRITE(iterChar,'(I10.10)') iter

      OPEN( UNIT=NEWUNIT(fUnit), &
            FILE='Burgers1D.'//iterChar//'.pickup', &
            FORM='unformatted',&
            ACCESS='direct',&
            STATUS='old',&
            ACTION='READ',&
            CONVERT='big_endian',&
            RECL=prec*(nS+1) )
     
      thisRec = 1
      DO iEl = 1, myDGSEM % nElems
        
         READ( fUnit, REC=thisRec )sol(:,1) 
         thisRec = thisRec+1
         CALL myDGSEM % SetSolution( iEl, sol )
        
      ENDDO

      CLOSE(UNIT=fUnit)

 END SUBROUTINE ReadPickup_Burgers1D
!
!
! 
 END MODULE BurgersClass



