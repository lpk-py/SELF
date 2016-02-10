PROGRAM CGSEM2D_MAIN
!
! 
!
!
! =============================================================================== !
!
!
! The model settings
USE MODEL_PRECISION
USE CONSTANTS_DICTIONARY
USE MODEL_FLAGS
!  
 USE COMMONROUTINES
 USE RUN_PARAMS_CLASS
 USE MATRIXROUTINES

 USE LEGENDRE
 USE LAGRANGE_1D_CLASS
 USE LAGRANGE_2D_CLASS

 USE QUADMESH_CLASS
 USE QUADELEMENT_CLASS
 USE MAPPEDGEOM_2D_CLASS
 USE GEOMETRY_BASICS

 USE CGSEM2D_CLASS
 
IMPLICIT NONE


  TYPE( CGSEM2D )       :: thisCGSM
  TYPE( LAG_INTERP2D )  :: plotInterp
  

  real(prec), allocatable :: Tmat(:,:)
  real(prec), allocatable :: plotNodes(:)  
  real(prec), allocatable :: resi(:)
  real(prec)              :: thisX, thisY, error

  integer :: nPlot, pDeg, iS, iP, iEl

  logical :: fileExists

 
      ! Initial build
      CALL thisCGSM % BUILD( )

      nPlot = thisCGSM % params % nPlot

      pDeg = thisCGSM % params % polyDeg

      

      ! Allocate the quadrature to plot mesh matrix
      ALLOCATE( Tmat(0:nPlot,0:pDeg) )
 
      ALLOCATE( plotNodes(0:nPlot) )

      ALLOCATE( resi(0:thisCGSM % params % cgIterMax) )
      

      ! Build the uniform plotting nodes
      do iS = 0, nPlot

         plotNodes(iS) = -ONE + (thisCGSM % params % dxPlot)*real(iS, prec)

      enddo


      CALL plotInterp % BUILD_WITHNODES( plotNodes, plotNodes )


      CALL PINTERP_MATRIX_1D(thisCGSM % cgStorage% interp % xInterp, plotInterp % xInterp, Tmat )  

   
     ! Set the source 
     do iEl = 1, thisCGSM % mesh % nElems
        do iP = 0,pDeg
           do iS = 0,pDeg

              thisx = thisCGSM % mesh % elements(iEl) % geometry % x(iS,iP)
              thisy = thisCGSM % mesh % elements(iEl) % geometry % y(iS,iP)

              thisCGSM % source(iS,iP,iEl) = SOURCE( thisX, thisY )


           enddo
        enddo
     enddo

  
     CALL thisCGSM % CONJUGATE_GRADIENT_CGSEM2D( resi )

     CALL CALC_L2ERROR( thisCGSM, error )

     CALL thisCGSM % WRITE_TECPLOT( nPlot, pDeg, plotInterp, Tmat, 'State' )


     open(unit=200,file='resi.curve')
     write(200,*)'# residual'
     do iS = 0, thisCGSM % params % cgIterMax

        write(200,'(I6,2x,E18.8)') iS, resi(iS)

     enddo
     close(unit=200)


     INQUIRE(file='L2Error-4x4.curve', exist=fileExists)

     if( fileExists )then


        open(unit=200,file='L2Error-4x4.curve',STATUS='OLD', ACCESS='APPEND')
      
     else

        open(unit=200,file='L2Error-4x4.curve',STATUS='NEW')
        write(200,*)'# L2Error'

     endif

        write(200,'(I6,2x,E18.8)') pDeg, error

     close(unit=200)


   ! ------------------------------------------------------------------------ !
   ! Clean up memory

     DEALLOCATE( Tmat )
 
     DEALLOCATE( plotNodes  )

     DEALLOCATE( resi )

     CALL thisCGSM % TRASH( )

     CALL plotInterp % TRASH( )


 CONTAINS
!
!
!
 FUNCTION SOURCE( x, y ) RESULT( s )
 !
 !
 ! ===============================================
   IMPLICIT NONE
   real(prec) :: x, y
   real(prec) :: s

      s = 8.0_prec*pi*pi*sin(2.0_prec*pi*thisX)*sin(2.0_prec*pi*thisY)

 END FUNCTION SOURCE
!
!
!
 SUBROUTINE CALC_L2ERROR( myCGSEM, error )
 !
 ! 
 ! ================================================
 ! DECLARATIONS
   IMPLICIT NONE
   TYPE( CGSEM2D ), intent(inout) :: myCGSEM
   real(prec), intent(out)    :: error
   ! LOCAL
   integer :: iS, iP, nS, nP, iEl, nEl
   real(prec) :: x, y, J, w1, w2
   real(prec) :: errorfield(0:myCGSEM % cgStorage % nS, &
                            0:myCGSEM % cgStorage % nP, &
                            1:myCGSEM % mesh % nElems  )

      nS = myCGSEM % cgStorage % nS
      nP = myCGSEM % cgStorage % nP
      nEl = myCGSEM % mesh % nElems

      error = ZERO 

      do iEl = 1, nEl
         do iP = 0, nP
            do iS = 0, nS
 
               x = myCGSEM % mesh % elements(iEl) % geometry % x(iS,iP)
               y = myCGSEM % mesh % elements(iEl) % geometry % y(iS,iP)
               J = myCGSEM % mesh % elements(iEl) % geometry % J(iS,iP)
               w1 = myCGSEM % cgStorage % qWeightX(iS)
               w2 = myCGSEM % cgStorage % qWeightY(iP)              

               errorField(iS,iP,iEl) = ( myCGSEM % solution(iS,iP,iEl) + sin(2.0_prec*pi*x)*sin(2.0_prec*pi*y) )**2*J*w1*w2


            enddo
         enddo
      enddo

 
      CALL myCGSEM % MASK( errorField )

       do iEl = 1, nEl
         do iP = 0, nP
            do iS = 0, nS            

               error = error + errorField(iS,iP,iEl) 


            enddo
         enddo
      enddo
 
      error = sqrt(error)
      
 END SUBROUTINE CALC_L2ERROR





END PROGRAM CGSEM2D_MAIN
