PROGRAM TEST_CGSEM1D
!
! May 2015
!
!  This program is written in order to test the modules CGSEM_1D_CLASS.f90 and CGSEM_1D_FEM_PRECONDITIONER.f90 .
!  The heat equation is solved using the trapezoidal rule as implemented in the module CGSEM_1D_CLASS.f90. The 
!  implicit system is inverted using the Preconditioned-Conjugate Gradient Method with the Finite Element 
!  preconditioner.
!
!  The initial condition is chosen as
!
!   u = e^(-x^2) 
!  
!  on the interval [-8, 8]
!
!  One can control the number of elements by specifying the parameter "nEl" in the program here. The polynomial
!  degree can be set in the runtime.params file. Additionaly, the "RUN_PARAMS_CLASS" is modified to include
!  the max conjugate gradient iterations and cg tolerance - this way these parameters can be modified at runtime
!  instead of requiring recompilation.
!
! =============================================================================================================== !


USE COMMONDATA
USE RUN_PARAMS_CLASS

USE LAGRANGE_1D_CLASS

USE CGSEM_1D_CLASS
USE CGSEM_1D_FEM_PRECONDITIONER

IMPLICIT NONE

   integer, parameter  :: nEl = 3  ! The number of elements

   TYPE(RUN_PARAMS)           :: theParams
   TYPE(CGSEM_1D)             :: thisCGSEM
   TYPE(FEM1D_PRECONDITIONER) :: thisFEM
   TYPE(LAG_INTERP1D)         :: c2f

   real(prec)                 :: x(0:nEl), dx, tn, dt, thisX, l2Error
   real(prec), allocatable    :: plotNodes(:), c2fmat(:,:), error(:,:)

   integer :: iNode, iEl, iter0, iterF, iT, dFreq, nP
   character(10) :: iterChar
   character(5)  :: pDegChar


     ! Build the runtime parameters

     CALL theParams % BUILD( )

     write(pDegChar,'(I5.5)')theParams % polyDeg

     ! ================= GEOMETRY BUILD ================= !
     ! Build the element boundaries

     dx = 16.0_prec/real(nEl,prec)

     do iNode = 0,nEl

        x(iNode) = -8.0_prec + dx*real(iNode,prec)

     enddo


     ! ================================================== !


     ! Build the plot nodes (in computational space)
     nP = theParams % nPlot

     ALLOCATE( plotNodes( 0:nP ) )
     ALLOCATE( c2fmat( 0:nP, 0:theParams % polyDeg) )

     dx = 2.0_prec/real(nP,prec)

     do iNode = 0, nP

        plotNodes(iNode) = -ONE + dx*real(iNode,prec)

     enddo


     CALL c2f % BUILD_WITHNODES( plotNodes )


     ! Build the CGSEM data structure
     CALL thisCGSEM % BUILD( theParams % polyDeg, nEl, x(0:nel) )
    

     ! Build the preconditioner
     CALL thisFEM % BUILD( thisCGSEM % mesh, theParams % polyDeg )


     CALL PINTERP_MATRIX_1D( thisCGSEM % nodalStorage % interp , &
                              c2f, c2fmat )

     ! ================= INITIAL CONDITION ================= !
     ! Fill in the exact solution
     ALLOCATE( error(0:theParams % polyDeg, 1:nEl) )

     do iEl = 1,nEl

        do iNode = 0, theParams % polyDeg

           thisX = thisCGSEM % mesh % elements(iEl) % geometry % x(iNode)

           thisCGSEM % sol(iNode,iEl) = EXACT( thisX, ZERO )

           error(iNode, iEl) = EXACT( thisX, ZERO ) - thisCGSEM % sol(iNode, iEl) 
        enddo

     enddo

     ! ================================================== !

     ! Calculate the L-2 error
    
     CALL thisCGSEM % VECTOR_PRODUCT( error, error, l2Error )
 
     l2Error = sqrt(l2Error) 

     open( unit = 200, file = 'L2Error.'//pDegChar//'.curve', &
           form = 'FORMATTED', access= 'SEQUENTIAL', &
           status = 'REPLACE' )

     write(200,*) '#L2-error'
     write(200,'(F20.10,2x,F20.10)') ZERO, l2Error

     ! ================= TIME INTEGRATION ================= !

     iter0 = theParams % iterInit
     iterf = iter0 + theParams % nTimeSteps-1

     dt = theParams % dt
     dFreq = theParams % dumpFreq


     ! FILE I/O

     write(iterChar,'(I10.10)') 0

     CALL thisCGSEM % WRITE_TECPLOT(  c2fmat, nP, c2f, 'State.'//iterChar )
     CALL WRITE_EXACT( thisCGsem, c2fmat, nP, ZERO, c2f, 'Exact.'//iterChar )

      

     do iT = iter0, iterf

        tn = real(iT,prec)*dt ! the time at the beginning of this interval


        CALL thisCGSEM % TRAPEZOIDAL_INTEGRATION( tn, dt, &
                                                  theParams % cgIterMax, theParams % cgTol,&
                                                  thisFEM )

        ! FILE I/O
        if( mod( iT+1, dFreq) == 0 ) then

          write(iterChar,'(I10.10)') iT+1

          CALL thisCGSEM % WRITE_TECPLOT(  c2fmat, nP, c2f, 'State.'//iterChar )

          CALL WRITE_EXACT( thisCGsem, c2fmat, nP, tn+dt, c2f, 'Exact.'//iterChar )

          !Calculate the error
          do iEl = 1,nEl

             do iNode = 0, theParams % polyDeg

                thisX = thisCGSEM % mesh % elements(iEl) % geometry % x(iNode)

                error(iNode, iEl) = EXACT( thisX, tn+dt ) - thisCGSEM % sol(iNode, iEl) 
             enddo

          enddo

          CALL thisCGSEM % VECTOR_PRODUCT( error, error, l2Error )
 
          l2Error = sqrt(l2Error) 

          write(200,'(F20.10,2x,F20.10)') tn+dt, l2Error

        endif


     enddo !iT, loop over time


   ! DEALLOCATE AND CLEAR MEMORY !

   CALL thisCGSEM % TRASH( )
   CALL thisFEM % TRASH( )

   CALL c2f % TRASH( )

   DEALLOCATE( plotNodes, c2fmat, error )

   close(unit=200)
   

   CONTAINS


 FUNCTION EXACT (xi, t) RESULT( uAtxiT )
    real(prec) :: xi, t, uAtxiT

       uAtxiT = exp( -xi*xi/( 4.0_prec*t + ONE ) )/sqrt(4.0_prec*t + ONE)

 END FUNCTION EXACT

 SUBROUTINE WRITE_EXACT( myCGsem, Tx, nPlot, t, plotInterp, filename )
 !
 !
 ! 
 ! ==============================================================================!
 ! DECLARATIONS
   TYPE( CGSEM_1D ), intent(in)    :: myCGSEM
   integer, intent(in)              :: nPlot
   real(prec), intent(in)           :: Tx(0:nPlot,0:myCGSEM % N), t
   TYPE( LAG_INTERP1D ), intent(in) :: plotInterp
   character(*), intent(in)         :: filename
   ! LOCAL
   real(prec)   :: x(0:nPlot), p(0:nPlot)
   integer      :: iEl, iX
   character(5) :: zoneID

    open( unit=2, file= trim(filename)//'.curve', form='formatted',status='replace')

    write(2,*) '#heat'
    

    do iEl = 1, myCGSEM % K

       ! Interpolate the solutions onto a uniform mesh
       CALL COARSE_TO_FINE_1D(myCGsem % nodalStorage % interp,&
                                      myCGsem % mesh % elements(iEl) % geometry % x,&
                                      plotInterp, x, Tx )
 
     

     !   write(zoneID,'(I5.5)') iEl


 !       write(2,*) 'ZONE T="el'//trim(zoneID)//'", I=',nPlot+1,',F=POINT'

        do iX = 0, nPlot-1

           write(2,'(F15.7,2x,F15.7)') x(iX), EXACT( x(iX), t )

        enddo
        

    enddo

    write(2,'(F15.7,2x,F15.7)') x(nPlot), EXACT( x(nPlot), t )
    close(unit=2)


 END SUBROUTINE WRITE_EXACT
END PROGRAM TEST_CGSEM1D
