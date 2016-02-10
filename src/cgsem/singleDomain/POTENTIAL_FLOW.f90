PROGRAM POTENTIAL_FLOW

! 



 USE COMMONDATA 
 USE COMMONROUTINES
 USE RUN_PARAMS_CLASS
 USE MATRIXROUTINES

 USE LEGENDRE
 USE LAGRANGE_1D_CLASS
 USE LAGRANGE_2D_CLASS

 USE MAPPEDGEOM_2D_CLASS
 USE GEOMETRY_BASICS

 USE CGSM2D_CLASS
 
IMPLICIT NONE

  real(prec), parameter :: rOut = 10.0_prec
  real(prec), parameter :: rIn = 0.5_prec

  TYPE( CGSM2D )        :: thisCGSM
  TYPE( CURVE_2D )      :: boundCurves(1:4)
  TYPE( LAG_INTERP2D )  :: plotInterp
  

  real(prec), allocatable :: Tmat(:,:), Lu(:,:), u(:,:), diags(:,:,:)
  real(prec), allocatable :: plotNodes(:)  
  real(prec), allocatable :: geomNodes(:), geomWeights(:)
  real(prec), allocatable :: x(:), y(:), resi(:)
  real(prec)              :: thisX, thisY, error, theta

  integer :: nPlot, pDeg, iS, iP

 
   ! Setup

      ! Read in the model parameters from runtime.params
      CALL thisCGSM % params % BUILD( )

      nPlot = thisCGSM % params % nPlot

      pDeg = thisCGSM % params % polyDeg

      ! Allocate the quadrature to plot mesh matrix
      ALLOCATE( Tmat(0:nPlot,0:pDeg), diags(0:pDeg,0:pDeg,1:3) )
 
      ALLOCATE( plotNodes(0:nPlot), geomNodes(0:pDeg), geomWeights(0:pDeg), u(0:pDeg,0:pDeg), Lu(0:pDeg,0:pDeg) )

      ALLOCATE( x(0:pDeg), y(0:pDeg), resi(0:thisCGSM % params % cgIterMax) )
      

      ! Build the uniform plotting nodes
      do iS = 0, nPlot

         plotNodes(iS) = -ONE + thisCGSM % params % dxPlot*real(iS, prec)

      enddo


      CALL plotInterp % BUILD_WITHNODES( plotNodes, plotNodes )


   ! Build the geometry

      CALL LEG_GAUSSLOBATTO_QUAD(pDeg, geomNodes, geomWeights) 

      
      do iS = 0, pDeg

         x(iS) = rIn + 0.5_prec*(rOut-rIn)*( geomNodes(iS) + ONE )

         y(iS) = ZERO

      enddo

      CALL boundCurves(1) % BUILD( x, y, geomNodes, pDeg)

      do iS = 0, pDeg

         theta = 0.5_prec*pi*( geomNodes(iS) + ONE )

         x(iS) = rOut*cos( theta )

         y(iS) = rOut*sin( theta )

      enddo
 
      CALL boundCurves(2) % BUILD( x, y, geomNodes, pDeg)

      do iS = 0, pDeg

         x(iS) = -rIn - 0.5_prec*(rOut-rIn)*( geomNodes(iS) + ONE )

         y(iS) = ZERO

      enddo

      CALL boundCurves(3) % BUILD( x, y, geomNodes, pDeg)

      do iS = 0, pDeg

         theta = 0.5_prec*pi*( geomNodes(iS) + ONE )

         x(iS) = rIn*cos( theta )

         y(iS) = rIn*sin( theta )

      enddo
   
      CALL boundCurves(4) % BUILD( x, y, geomNodes, pDeg)

      
   ! Build the data structure

      CALL thisCGSM % BUILD( boundCurves )

      CALL PINTERP_MATRIX_1D(thisCGSM % cgStorage% interp % xInterp, plotInterp % xInterp, Tmat )  

   ! Set the initial guess and the source terms
     

     ! Set the source 
  
     CALL thisCGSM % CONJUGATE_GRADIENT_CGSM2D( resi )

     !CALL CALC_L2ERROR( thisCGSM, error )
     ! Verify gradient and laplacian
      
     CALL thisCGSM % CALC_GRADIENT( thisCGSM % solution, diags(:,:,1), diags(:,:,2), pDeg, pDeg )


   ! Write to file
     do iP = 0, pDeg
        do iS = 0, pDeg

           diags(iS,iP,3) = ONE - 0.5_prec*( ONE - (diags(iS,iP,1)**2 + diags(iS,iP,2)**2) )

        enddo
     enddo
 

     

 
     CALL thisCGSM % WRITE_TECPLOT( nPlot, pDeg, plotInterp, Tmat, 'State' )

     CALL  WRITE_TECPLOT_DIAGS( thisCGSM, nPlot, pDeg, plotInterp, Tmat, diags, 3, 'Diags')



     open(unit=200,file='resi.curve')
     write(200,*)'# residual'
     do iS = 0, thisCGSM % params % cgIterMax

        write(200,'(I6,2x,E18.8)') iS, resi(iS)

     enddo
     close(unit=200)

    ! open(unit=200,file='L2Error.curve',ACCESS='APPEND')
    ! write(200,*)'# L2Error'

    !    write(200,'(I6,2x,E18.8)') pDeg, error

    ! close(unit=200)


   ! ------------------------------------------------------------------------ !
   ! Clean up memory

     DEALLOCATE( Tmat, diags )
 
     DEALLOCATE( plotNodes, geomNodes, geomWeights )

     DEALLOCATE( x, y, resi )

     CALL thisCGSM % TRASH( )

     CALL plotInterp % TRASH( )

     CALL boundCurves(1) % TRASH( )
     CALL boundCurves(2) % TRASH( )
     CALL boundCurves(3) % TRASH( )
     CALL boundCurves(4) % TRASH( )


 CONTAINS
!
!
!
 SUBROUTINE CALC_L2ERROR( myCGSM, error )
 !
 ! 
 ! ================================================
 ! DECLARATIONS
   IMPLICIT NONE
   TYPE( CGSM2D ), intent(in) :: myCGSM
   real(prec), intent(out)    :: error
   ! LOCAL
   integer :: iS, iP, nS, nP
   real(prec) :: x, y, J, w1, w2

      nS = myCGSM % cgStorage % nS
      nP = myCGSM % cgStorage % nP

      error = ZERO

      do iP = 0, nP
         do iS = 0, nS
 
            x = myCGSM % geometry % x(iS,iP)
            y = myCGSM % geometry % y(iS,iP)
            J = myCGSM % geometry % J(iS,iP)
            w1 = myCGSM % cgStorage % qWeightX(iS)
            w2 = myCGSM % cgStorage % qWeightY(iP)              

            error = error + ( myCGSM % solution(iS,iP) - cos(2.0_prec*pi*x)*sin(2.0_prec*pi*y) )**2*J*w1*w2

         enddo
      enddo

 END SUBROUTINE CALC_L2ERROR




 SUBROUTINE WRITE_TECPLOT_DIAGS( myCGSM, nPlot, nOld, plotInterp, Tmat, diags, nDiags, filename )
 ! WRITE_TECPLOT_CGSM2D
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
 !    real(prec) :: fAtNodes(:,:) - the interpolation nodes function values.
 !
 !    character(*) :: filename - name of output file
 !
 !
 ! 
 !------------------------------------------------------------------------
  IMPLICIT NONE
  TYPE( CGSM2D ), intent(in)            :: myCGSM
  integer, intent(in)                   :: nOld, nPlot, nDiags
  real(prec), intent(in)                :: diags(0:nOld,0:nOld,nDiags)
  real(prec), intent(in)                :: Tmat(0:nPlot, 0:nOld)
  TYPE( LAG_INTERP2D ), intent(in)      :: plotInterp
  character(*), intent(in)              :: filename
  !LOCAL
  real(prec)  :: x(0:nPlot,0:nPlot), y(0:nPlot,0:nPlot)!, depth(0:nPlot,0:nPlot)
  real(prec)  :: u(0:nPlot,0:nPlot), v(0:nPlot,0:nPlot), p(0:nPlot,0:nPlot)
  real(prec)  :: dxds(0:nPlot,0:nPlot), dxdp(0:nPlot,0:nPlot), dyds(0:nPlot,0:nPlot), dydp(0:nPlot,0:nPlot)
  integer :: iX, iY, iZ, iEl
  character(len=5) :: zoneID


    open( unit=2, file= trim(filename)//'.tec', form='formatted',status='replace')

    write(2,*) 'VARIABLES = "X", "Y", "U", "V", "Pressure", "dxds", "dxdp", "dyds","dydp"'
    


       ! Interpolate the solutions onto a uniform mesh
       CALL COARSE_TO_FINE_2D(myCGSM % cgStorage % interp,&
                                      myCGSM % geometry % x,&
                                      plotInterp, x, Tmat, Tmat)


       CALL COARSE_TO_FINE_2D(myCGSM % cgStorage % interp,&
                                      myCGSM %  geometry % y,&
                                      plotInterp, y, Tmat, Tmat)

       CALL COARSE_TO_FINE_2D(myCGSM % cgStorage % interp,&
                                      myCGSM %  geometry % dxds,&
                                      plotInterp, dxds, Tmat, Tmat)

       CALL COARSE_TO_FINE_2D(myCGSM % cgStorage % interp,&
                                      myCGSM %  geometry % dxdp,&
                                      plotInterp, dxdp, Tmat, Tmat)

       CALL COARSE_TO_FINE_2D(myCGSM % cgStorage % interp,&
                                      myCGSM %  geometry % dyds,&
                                      plotInterp, dyds, Tmat, Tmat)

       CALL COARSE_TO_FINE_2D(myCGSM % cgStorage % interp,&
                                      myCGSM %  geometry % dydp,&
                                      plotInterp, dydp, Tmat, Tmat)

       CALL COARSE_TO_FINE_2D(myCGSM % cgStorage % interp,&
                                      diags(:,:,1),&
                                      plotInterp, u, Tmat, Tmat)
 
        CALL COARSE_TO_FINE_2D(myCGSM % cgStorage % interp,&
                                      diags(:,:,2),&
                                      plotInterp, v, Tmat, Tmat)

        CALL COARSE_TO_FINE_2D(myCGSM % cgStorage % interp,&
                                      diags(:,:,3),&
                                      plotInterp, p, Tmat, Tmat)

        write(zoneID,'(I5.5)') 0


        write(2,*) 'ZONE T="el'//trim(zoneID)//'", I=',nPlot+1,', J=', nPlot+1,',F=POINT'

        do iY = 0, nPlot
           do iX = 0, nPlot

              write (2,*)  x( iX, iY ), y( iX, iY ), u(iX,iY), v(iX,iY), p(iX,iY), &
                           dxds(iX,iY), dxdp(iX,iY), dyds(iX,iY), dydp(iX,iY) 

           enddo
        enddo


    close(unit=2)

 END SUBROUTINE WRITE_TECPLOT_DIAGS

END PROGRAM POTENTIAL_FLOW
