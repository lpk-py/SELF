! ShallowWaterDiagnostics.f90
! 
! Copyright 2015 Joseph Schoonover <schoonover.numerics@gmail.com>
! 
! ShallowWaterDiagnostics.f90 is part of the Spectral Element Libraries in Fortran (SELF).
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
 
 
MODULE ShallowWaterShallowWaterDiagnosticsClass
! ========================================= Logs ================================================= !
!2016-05-11  Joseph Schoonover  schoonover.numerics@gmail.com 
!
! //////////////////////////////////////////////////////////////////////////////////////////////// ! 

! src/common/
USE ModelPrecision
USE ConstantsDictionary
USE ModelFlags
USE CommonRoutines
! src/interp/
USE Legendre
USE Lagrange_2D_Class
! src/nodal/
USE NodalStorage_2D_Class
USE DGSEM_SolutionStorageClass_2D
! src/geometry/
USE EdgeClass
USE QuadElementClass  
USE QuadMeshClass    
! src/highend/shallowwater/
USE ShallowWaterClass
! Nocturnal Aviation classes and extensions
!USE FTTimerClass
!USE TIMING


! Parallel Libraries
USE OMP_LIB


IMPLICIT NONE
!
! Some parameters that are specific to this module





! ==============================================================================================!
!     The SHALLOWWATER_DIAGNOSTICS class contains a collection of
!     of diagnostics that are typically useful in understanding 
!     shallow-water physics
!     (So far) This includes 
!          
!           o       KE       :: Kinetic Energy
!           o       PE       :: Potential Energy
!           o       PV       :: Potential Vorticity 
! ==============================================================================================!

    TYPE ShallowWaterShallowWaterDiagnostics
      INTEGER                 :: nS, nP, nElems, nT
  
      REAL(prec), allocatable :: KE(:,:,:)
      REAL(prec), allocatable :: PE(:,:,:)
      REAL(prec), allocatable :: PV(:,:,:)

      ! Area integrated terms
      REAL(prec), allocatable :: intKE(:) 
      REAL(prec), allocatable :: intPE(:)
      REAL(prec), allocatable :: modelTime(:)  ! stores the times at which the solution is stored
      CONTAINS

      ! Manual constructors/Destructors
      PROCEDURE :: BUILD => BUILD_SHALLOWWATER_DIAGNOSTICS
      PROCEDURE :: TRASH => TRASH_SHALLOWWATER_DIAGNOSTICS
      ! Accessors

      
       ! Type Specific Routines
       PROCEDURE :: UPDATE_DIAGNOSTICS
       PROCEDURE :: COMPUTE_POTENTIAL_VORTICITY
       PROCEDURE :: DIAGNOSE_ENERGIES
       PROCEDURE :: DIAGNOSE_INTEGRATED_ENERGIES

       ! File I/O
       PROCEDURE :: WRITE_TECPLOT => WRITE_SHALLOWWATER_DIAGS_TECPLOT
       PROCEDURE :: WRITE_ENERGY_CURVE
 
    END TYPE SHALLOWWATER_DIAGNOSTICS
     

 CONTAINS
!
!
!==================================================================================================!
!------------------------------- Manual Constructors/Destructors ----------------------------------!
!==================================================================================================!
!
!
 SUBROUTINE Build_ShallowWaterDiagnostics( mySWdiags, params, nElems )
 ! S/R Build
 ! 
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS( ShallowWaterDiagnostics ), INTENT(out) :: mySWdiags
   TYPE( SWParams ), INTENT(in)      :: params
   INTEGER, INTENT(in)               :: nElems
   ! LOCAL
   INTEGER :: nS, nP, nT

      nS = params % polyDeg
      nP = nS

      mySWdiags % nS = nS
      mySWdiags % nP = nP
      mySWdiags % nElems = nElems

      nT = params % nTimeSteps/params % dumpFreq + 1! The number of time steps output is recorded
      mySWDiags % nT = nT
 
      ALLOCATE( mySWdiags % KE(0:nS,0:nP,1:nElems) )
      ALLOCATE( mySWdiags % PE(0:nS,0:nP,1:nElems) )
      ALLOCATE( mySWdiags % PV(0:nS,0:nP,1:nElems) )
      ALLOCATE( mySWdiags % intKE(0:nT) )
      ALLOCATE( mySWdiags % intPE(0:nT) )
      ALLOCATE( mySWdiags % modelTime(0:nT) )
   
      mySWdiags % KE = ZERO
      mySWdiags % PE = ZERO
      mySWdiags % PV = ZERO
      mySWdiags % intKE = ZERO
      mySWdiags % intPE = ZERO
      mySWdiags % modelTime = ZERO

 END SUBROUTINE Build_ShallowWaterDiagnostics
!
!
!
 SUBROUTINE Trash_ShallowWaterDiagnostics( mySWdiags )
 ! S/R Trash
 ! 
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS( ShallowWaterDiagnostics ), INTENT(inout) :: mySWdiags

      DEALLOCATE( mySWdiags % KE )
      DEALLOCATE( mySWdiags % PE )
      DEALLOCATE( mySWdiags % PV )
      DEALLOCATE( mySWdiags % intKE )
      DEALLOCATE( mySWdiags % intPE )
      DEALLOCATE( mySWdiags % modelTime )  

 END SUBROUTINE Trash_ShallowWaterDiagnostics
!
!
!=========================================================================!
!--------------------- TYPE-SPECIFIC ROUTINES ----------------------------!
!=========================================================================!
!
!
 SUBROUTINE UPDATE_DIAGNOSTICS( myDiag, myDGSEM, K, time )
 !
 !
 !
 ! ========================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS( SHALLOWWATER_DIAGNOSTICS ), intent(inout) :: myDiag
   TYPE( DGSEM_SHALLOWWATER ), intent(in)           :: myDGSEM
   INTEGER, intent(in)                              :: K
   REAL(prec), intent(in)                           :: time

      CALL myDiag % DIAGNOSE_ENERGIES( myDGSEM )


      CALL myDiag % DIAGNOSE_INTEGRATED_ENERGIES( K, time, &
                                                   myDGSEM % dgStorage, &
                                                   myDGSEM % mesh, &
                                                   myDGSEM % params )

      CALL myDiag % COMPUTE_POTENTIAL_VORTICITY( myDGSEM )

 END SUBROUTINE UPDATE_DIAGNOSTICS
!
!
!
 SUBROUTINE COMPUTE_POTENTIAL_VORTICITY( mydiag, myDGSEM ) 
 ! SUBROUTINE COMPUTE_POTENTIAL_VORTICITY
 ! 
 ! Description :
 !  
 ! 
 !   
 !    Subroutine Dependencies :
 !    
 !
 ! Input :
 !
 ! 
 ! Output :
 !
 !
 !------------------------------------------------------------------------
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS(SHALLOWWATER_DIAGNOSTICS), intent(inout) :: mydiag
   TYPE(DGSEM_SHALLOWWATER), intent(in)           :: myDGSEM
   ! LOCAL
   REAL(prec) :: f0, betaX, betaY
   INTEGER    :: formulation
   INTEGER    :: nS, nP, nEl
   INTEGER    :: iS, iP, iEl
   REAL(prec) :: u1(0:myDGSEM % dgStorage % nS,0:myDGSEM % dgStorage % nP)
   REAL(prec) :: u2(0:myDGSEM % dgStorage % nS,0:myDGSEM % dgStorage % nP)
   REAL(prec) :: u2North(0:myDGSEM % dgStorage % nS), u2South(0:myDGSEM % dgStorage % nS)
   REAL(prec) :: u1East(0:myDGSEM % dgStorage % nP), u1West(0:myDGSEM % dgStorage % nP)
   REAL(prec) :: pvort(0:myDGSEM % dgStorage % nS,0:myDGSEM % dgStorage % nP)
   REAL(prec) :: rvort(0:myDGSEM % dgStorage % nS,0:myDGSEM % dgStorage % nP)
   REAL(prec) :: u, v, H, normVec(1:2)
   REAL(prec) :: x, y, dxds, dxdp, dyds, dydp

      nS  = myDGSEM % dgStorage % nS
      nP  = myDGSEM % dgStorage % nP
      nEl = myDGSEM % mesh % nElems

      f0          = myDGSEM % params % f0
      betaX       = myDGSEM % params % betaX
      betaY       = myDGSEM % params % betaY
      formulation = myDGSEM % params % MODEL_FORMULATION

      

      DO iEl = 1, nEl

            ! Calculate the planetary vorticity
         DO iP = 0, nP
            DO iS = 0, nS                                   

               x = myDGSEM % mesh % elements(iEl) % geometry % x(iS,iP)
               y = myDGSEM % mesh % elements(iEl) % geometry % y(iS,iP)

               pvort(iS,iP) =  f0 + betaY*y +betaX*x ! Planetary contribution ONLY
                                          
            ENDDO
         ENDDO


         IF( formulation == LINEAR )THEN
      
            DO iP = 0, nP
               DO iS = 0, nS

                  H = myDGSEM % swAddons(iEl) % h(iS,iP)
                  mydiag % PV(iS,iP,iEl) = pvort(iS,iP)/H

               ENDDO
            ENDDO

         ELSE ! This is a nonlinear formulation         


            ! First calculate the relative vorticity

            ! Get the interior contravariant velocity
            DO iP = 0, nP ! Loop over the y-points

               DO iS = 0, nS ! Loop over the x-points

                  u = myDGSEM % sol(iEl) % solInterior(iS,iP,1)
                  v = myDGSEM % sol(iEl) % solInterior(iS,iP,2)
                  
                  dxds = myDGSEM % mesh % elements(iEl) % geometry % dxds(iS,iP)
                  dxdp = myDGSEM % mesh % elements(iEl) % geometry % dxdp(iS,iP)
                  dyds = myDGSEM % mesh % elements(iEl) % geometry % dyds(iS,iP)
                  dydp = myDGSEM % mesh % elements(iEl) % geometry % dydp(iS,iP)

            ! Compute the cross product of the velocities and the contravariant basis vector
                  u1(iS,iP) = v*dydp + u*dxdp

                  u2(iS,iP) = -(v*dyds + u*dxds)
  
               ENDDO ! iS loop over the x-points

            ENDDO ! iP, loop over the y-points


           ! Get the boundary values of the contravariant velocity
           ! Two ways to DO this : 
           ! (1) Use the already computed boundary solutions, and multiply by appropriate metrics to rotate velocities
           !  
           ! (2) Use u1 and u2 previously computed and then interpolate to the boundaries
           !
           ! Currently style (1) is used

     
 
           ! Get the east and west velocities
           DO iP = 0, nP

              u = myDGSEM % sol(iEl) % solBound(iP,1,2) 

              v = myDGSEM % sol(iEl) % solBound(iP,2,2)

              normVec(1) = DOT_PRODUCT(myDGSEM % mesh % elements(iEl) % geometry % dydp(:,iP),&
                                       myDGSEM % dGStorage % lagEast) !east

              normVec(2) = DOT_PRODUCT(myDGSEM % mesh % elements(iEl) % geometry % dxdp(:,iP),&
                                       myDGSEM % dGStorage % lagEast) !east

              u1East(iP) = v*normVec(1) + u*normVec(2)

              u = myDGSEM % sol(iEl) % solBound(iP,1,4) 

              v = myDGSEM % sol(iEl) % solBound(iP,2,4) 

              normVec(1) = DOT_PRODUCT( myDGSEM % mesh % elements(iEl) % geometry % dydp(:,iP),&
                                        myDGSEM % dGStorage % lagWest ) !west

              normVec(2) = DOT_PRODUCT( myDGSEM % mesh % elements(iEl) % geometry % dxdp(:,iP),&
                                        myDGSEM % dGStorage % lagWest ) !west

              u1West(iP) = v*normVec(1) + u*normVec(2)
 
              ! Compute du1/ds

              rvort(:,iP) = MATVECMUL( myDGSEM % dGStorage % dMatX, u1(:,iP), nS, nS )
  
           ENDDO ! Loop over all of the x -points


           ! Get the north and south velocities
           DO iS = 0, nS

              u = myDGSEM % sol(iEl) % solBound(iP,1,3) 
  
              v = myDGSEM % sol(iEl) % solBound(iP,2,3)

              normVec(1) = DOT_PRODUCT( myDGSEM % mesh % elements(iEl) % geometry % dyds(iS,:),&
                                        myDGSEM % dGStorage % lagNorth  ) !north

              normVec(2) = DOT_PRODUCT( myDGSEM % mesh % elements(iEl) % geometry % dxds(iS,:),&
                                        myDGSEM % dGStorage % lagNorth ) !north

              u2North(iS) = -(v*normVec(1) + u*normVec(2))

              u = myDGSEM % sol(iEl) % solBound(iP,1,1)

              v = myDGSEM % sol(iEl) % solBound(iP,2,1)

              normVec(1) = DOT_PRODUCT( myDGSEM % mesh % elements(iEl) % geometry % dyds(iS,:),&
                                        myDGSEM % dGStorage % lagSouth ) !west
  
              normVec(2) = DOT_PRODUCT( myDGSEM % mesh % elements(iEl) % geometry % dxds(iS,:),&
                                        myDGSEM % dGStorage % lagSouth ) !west

              u2South(iS) = -(v*normVec(1) + u*normVec(2))

              ! Compute du2/dp and add to vorticity

              rvort(iS,:) = rvort(iS,:) + MATVECMUL( myDGSEM % dGStorage % dMatY, u2(iS,:), nP, nP )
  
           ENDDO ! Loop over all of the x -points

           ! Calculate using the DG-derivative
            DO iP = 0, nP
               DO iS = 0, nS                                   

                  rvort(iS,iP) = ( rvort(iS,iP) + &                          ! Relative contribution
                                   ( myDGSEM % dgStorage % lagEast(iS)*u1East(iP) - &
                                     myDGSEM % dgStorage % lagWest(iS)*u1West(iP) )/myDGSEM % dgStorage % qWeightX(iS) + &
                                   ( myDGSEM % dgStorage % lagNorth(iP)*u2North(iS) - &
                                     myDGSEM % dgStorage % lagSouth(iP)*u2South(iS) )/myDGSEM % dgStorage % qWeightY(iP) )/&
                                 ( myDGSEM % mesh % elements(iEl) % geometry % J(iS,iP) )

               ENDDO
            ENDDO

            ! relative vorticity calculation done

            ! Now we compute the PV
            IF( formulation == CONSERVATIVE )then

               DO iP = 0, nP
                  DO iS = 0, nS
                     ! Here, the fluid thickness is given by the third solution variable
                     H = myDGSEM % Sol(iEl) % solInterior(iS,iP,iEl)

                     mydiag % PV(iS,iP,iEl) = (pvort(iS,iP) + rvort(iS,iP))/H

                  ENDDO
               ENDDO

            ELSEIF( formulation == SKEW_SYMMETRIC )then
 
               DO iP = 0, nP
                  DO iS = 0, nS
                     ! Here, the fluid thickness is given by the third solution variable 
                     ! plus the background depth
                     H = myDGSEM % Sol(iEl) % solInterior(iS,iP,iEl) + &
                         myDGSEM % swAddons(iEl) % h(iS,iP)

                     mydiag % PV(iS,iP,iEl) = (pvort(iS,iP) + rvort(iS,iP))/H

                  ENDDO
               ENDDO

            ENDIF

         ENDIF

      ENDDO ! Loop over the elements                 
                     

 END SUBROUTINE COMPUTE_POTENTIAL_VORTICITY
!
!
!
 SUBROUTINE DIAGNOSE_ENERGIES( mySWdiags, myDGSEM )
 ! S/R DIAGNOSE_ENERGIES
 ! 
 ! Calculates the kinetic and potential energy from the solution stored in myDGSEM
 !
 !
 ! 
 !------------------------------------------------------------------------
  IMPLICIT NONE
  CLASS(SHALLOWWATER_DIAGNOSTICS ), intent(inout) :: mySWdiags
  TYPE(DGSEM_SHALLOWWATER), intent(in)            :: myDGSEM
  !LOCAL
  INTEGER :: iS, iP, iEl
  REAL(prec) :: H,KE,PE


    

    DO iEl = 1, myDGSEM % mesh % nElems

      DO iP = 0, myDGSEM % dgStorage % nP

         DO iS = 0, myDGSEM % dgStorage % nS

             IF( myDGSEM % params % MODEL_FORMULATION == SKEW_SYMMETRIC .OR. &
                 myDGSEM % params % MODEL_FORMULATION == LINEAR ) then
  

                H = myDGsem % swAddons(iEl) % h(iS,iP) + &
                    myDGsem % sol(iEl) % solInterior(iS,iP,3)  ! Total thickness

                KE = HALF*H*(myDGsem % sol(iEl) % solInterior(iS,iP,1)*&
                             myDGsem % sol(iEl) % solInterior(iS,iP,1) + &
                             myDGsem % sol(iEl) % solInterior(iS,iP,2)*&
                             myDGsem % sol(iEl) % solInterior(iS,iP,2) )

                PE = HALF*myDGsem % sol(iEl) % solInterior(iS,iP,3)*&
                          myDGsem % sol(iEl) % solInterior(iS,iP,3)*&
                          myDGSEM % params % g

             ELSE

                H = myDGsem % sol(iEl) % solInterior(iS,iP,3) ! Total thickness
                KE = HALF*(myDGsem % sol(iEl) % solInterior(iS,iP,1)*&
                           myDGsem % sol(iEl) % solInterior(iS,iP,1) + &
                           myDGsem % sol(iEl) % solInterior(iS,iP,2)*&
                           myDGsem % sol(iEl) % solInterior(iS,iP,2) )/H
              
                H = H -myDGsem % swAddons(iEl) % h(iS,iP) ! Thickness anomaly

                PE = HALF*H*H*myDGSEM % params % g           

             ENDIF

             mySWdiags % KE(iS,iP,iEl) = KE
             mySWdiags % PE(iS,iP,iEl) = PE 

         ENDDO

      ENDDO

    ENDDO

   
 END SUBROUTINE DIAGNOSE_ENERGIES
!
!
!
 SUBROUTINE DIAGNOSE_INTEGRATED_ENERGIES( mySWdiags, iT, t, dgStor, mesh, params )
 ! S/R DIAGNOSE_INTEGRATED_ENERGIES
 ! 
 ! Calculates the area integrated kinetic and potential energies and assigns the 
 ! value to the "iT-th" address in the arrays 
 !
 !    intKE(iT)
 !    intPE(iT)
 !
 ! 
 !------------------------------------------------------------------------
  IMPLICIT NONE
  CLASS(SHALLOWWATER_DIAGNOSTICS ), intent(inout) :: mySWdiags
  INTEGER, intent(in)                             :: iT
  REAL(prec), intent(in)                          :: t
  TYPE( NODAL_STORAGE_2D ), intent(in)            :: dGStor
  TYPE( QUADMESH ), intent(in)                    :: mesh
  TYPE( RUN_PARAMS ), intent(in)                  :: params
  !LOCAL
  INTEGER :: iX, iY, iEl
  REAL(prec) :: intKE, intPE, J, w1, w2


    intKE = ZERO
    intPE = ZERO

    DO iEl = 1, mySWdiags % nElems

      DO iY = 0, mySWdiags % nP

         DO iX = 0, mySWdiags % nS


              J = mesh % elements(iEl) % geometry % J(iX,iY)
              w1 = dGStor % qWeightX(iX)
              w2 = dGStor % qWeightY(iY)
 
              intKE  = intKE + mySWdiags % KE(iX,iY,iEl)*J*w1*w2
              intPE  = intPE + mySWdiags % PE(iX,iY,iEl)*J*w1*w2
             

         ENDDO

      ENDDO

    ENDDO

    mySWdiags % intKE(iT) = intKE
    mySWdiags % intPE(iT) = intPE
   
 END SUBROUTINE DIAGNOSE_INTEGRATED_ENERGIES
!
!
!=========================================================================!
!-------------------------- FILE I/O ROUTINES ----------------------------!
!=========================================================================!
!
!
 SUBROUTINE WRITE_SHALLOWWATER_DIAGS_TECPLOT( mySWDiags, myDGsem, Tx, Ty, nOld, nPlot, plotInterp, filename )
! S/R WRITE_SHALLOWWATER_DIAGS_TECPLOT
 !  Description :
 ! 
 !
 !    Subroutine dependencies :
 !    (NONE)
 !    
 !
 !  Input :
 !    type(LAG_INTERP2D) :: myPoly
 !
 !    REAL(prec) :: fAtNodes(:,:) - the interpolation nodes function values.
 !
 !    CHARACTER(*) :: filename - name of output file
 !
 !
 ! 
 !------------------------------------------------------------------------
  IMPLICIT NONE
  CLASS( SHALLOWWATER_DIAGNOSTICS ), intent(in) :: mySWDiags
  TYPE( DGSEM_SHALLOWWATER ), intent(in)        :: myDGsem
  INTEGER, intent(in)                           :: nOld, nPlot
  REAL(prec), intent(in)                        :: Tx(0:nPlot, 0:nOld), Ty(0:nPlot,0:nOld) 
  TYPE( LAG_INTERP2D ), intent(in)              :: plotInterp
  CHARACTER(*), intent(in)                      :: filename
  !LOCAL
  REAL(prec)  :: x(0:nPlot,0:nPlot), y(0:nPlot,0:nPlot)!, depth(0:nPlot,0:nPlot)
  REAL(prec)  :: ke(0:nPlot,0:nPlot)
  REAL(prec)  :: pe(0:nPlot,0:nPlot)
  REAL(prec)  :: pv(0:nPlot,0:nPlot)
  INTEGER :: iX, iY, iZ, iEl
  CHARACTER(len=5) :: zoneID


    open( unit=2, file= trim(filename)//'.tec', form='formatted',status='replace')

    WRITE(2,*) 'VARIABLES = "X", "Y", "Kinetic Energy", "Potential Energy", "Potential Vorticity" '
    
    DO iEl = 1, myDGsem % mesh % nElems

       ! Interpolate the solutions onto a uniform mesh
       CALL COARSE_TO_FINE_2D(myDGsem % dGStorage % interp,&
                              myDGsem % mesh % elements(iEl) % geometry % x,&
                              plotInterp, x, Tx, Ty)


       CALL COARSE_TO_FINE_2D(myDGsem % dGStorage % interp,&
                              myDGsem % mesh % elements(iEl) % geometry % y,&
                              plotInterp, y, Tx, Ty)

       CALL COARSE_TO_FINE_2D(myDGsem % dGStorage % interp,&
                              mySWDiags % KE(:,:,iEl),&
                              plotInterp, ke, Tx, Ty)
 
       CALL COARSE_TO_FINE_2D(myDGsem % dGStorage % interp,&
                              mySWDiags % PE(:,:,iEl),&
                              plotInterp, pe, Tx, Ty)

       CALL COARSE_TO_FINE_2D(myDGsem % dGStorage % interp,&
                              mySWDiags % PV(:,:,iEl),&
                              plotInterp, pv, Tx, Ty)

        WRITE(zoneID,'(I5.5)') myDGsem % mesh % elements(iEl) % globElID


        WRITE(2,*) 'ZONE T="el'//trim(zoneID)//'", I=',nPlot+1,', J=', nPlot+1,',F=POINT'!, ZONETYPE=FEQUADRILATERAL'

        DO iY = 0, nPlot
           DO iX = 0, nPlot

              WRITE (2,*)  x( iX, iY ), y( iX, iY ), ke(iX,iY), pe(iX,iY), pv(iX,iY)

           ENDDO
        ENDDO

      

    ENDDO


    close(unit=2)

    

 END SUBROUTINE WRITE_SHALLOWWATER_DIAGS_TECPLOT
!
!
!
 SUBROUTINE WRITE_ENERGY_CURVE( mySWDiags )
 ! S/R WRITE_ENERGY_CURVE
 !
 !
 ! ================================================================= !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS( SHALLOWWATER_DIAGNOSTICS ), intent(in) :: mySWDiags
   ! LOCAL
   INTEGER :: iT, fUnit
   REAL(prec) :: initKEpPE

      open( unit = NEWUNIT(fUnit), file='Energy.curve' )

      WRITE(fUnit,*)'#Kinetic-Energy'
      DO iT = 0, mySWDiags % nT
         WRITE(fUnit,'(E17.7,1x,E17.7)') mySWDiags % modelTime(iT), mySWDiags % intKE(iT)
      ENDDO 

      WRITE(fUnit,*)'#Potential-Energy'
      DO iT = 0, mySWDiags % nT
         WRITE(fUnit,'(E17.7,1x,E17.7)') mySWDiags % modelTime(iT), mySWDiags % intPE(iT)
      ENDDO 

      initKEpPE = mySWDiags % intKE(0)+mySWDiags % intPE(0)
      WRITE(fUnit,*)'#Scaled-Total-Energy'
      DO iT = 0, mySWDiags % nT
         WRITE(fUnit,'(E17.7,1x,E17.7)')  mySWDiags % modelTime(iT),&
                                        ( mySWDiags % intKE(iT)+mySWDiags % intPE(iT) )/initKEpPE
      ENDDO 
    

      close(unit=fUnit)


 END SUBROUTINE WRITE_ENERGY_CURVE

 END MODULE ShallowWaterDiagnosticsClass



