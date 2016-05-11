! ConstantsDictionary.f90
! 
! Copyright 2015 Joseph Schoonover <schoonover.numerics@gmail.com>
! 
! ConstantsDictionary.f90 is part of the Spectral Element Libraries in Fortran (SELF).
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
 
 

 MODULE ConstantsDictionary

  USE ModelPrecision

  !*************************************************************!
  ! ------------------ MATHEMATICAL CONSTANTS ------------------!
  ! ************************************************************!
  !                                                             !
  ! ------------------------------------------------------------!
    REAL(prec), PARAMETER :: pi   = 4.0_prec*atan(1.0_prec)
    REAL(prec), PARAMETER :: ZERO = 0.0_prec
    REAL(prec), PARAMETER :: ONE  = 1.0_prec
    REAL(prec), PARAMETER :: TWO  = 2.0_prec
    REAL(prec), PARAMETER :: HALF = 0.5_prec
    REAL(prec), PARAMETER :: TOL  = epsilon(1.0_prec)


    INTEGER, PARAMETER :: TransFiniteMetrics = 200 ! compute metric terms from transfinite interpolation
    INTEGER, PARAMETER :: InterpMetrics = 201      ! compute metric terms from interpolation of mesh locations

    INTEGER, PARAMETER :: nUniBounds = 2
    INTEGER, PARAMETER :: nTriBounds = 6 
    INTEGER, PARAMETER :: kItMax = 50 ! Max iterations for Newton's method.
                                      ! USEd in : MODULE LEGENDRE.f90, S/R LEG_GAUSSQUAD
    REAL(prec), PARAMETER :: fillValue = -9999.99_prec
    INTEGER, PARAMETER    :: fillValueInt = -99999
  !*************************************************************!
  ! ----------------- ROOT FINDER CONSTANTS --------------------!
  ! ************************************************************!
  !                                                             !
  ! ------------------------------------------------------------!
    REAL(prec), PARAMETER :: newtonTolerance = 10.0**(-10)
    INTEGER, PARAMETER    :: newtonMax       = 100
  
  !*************************************************************!
  ! ----------------- TIME STEPPING CONSTANTS ------------------!
  ! ************************************************************!
  !                                                             !
  ! ------------------------------------------------------------!
  ! Runge-Kutta 3rd Order, low storage constants
    REAL(prec), PARAMETER :: rk3_a(1:3) = (/ 0.0_prec, -5.0_prec/9.0_prec, -153.0_prec/128.0_prec /)
    REAL(prec), PARAMETER :: rk3_b(1:3) = (/ 0.0_prec, 1.0_prec/3.0_prec, 3.0_prec/4.0_prec /)
    REAL(prec), PARAMETER :: rk3_g(1:3) = (/ 1.0_prec/3.0_prec, 15.0_prec/16.0_prec, 8.0_prec/15.0_prec /)

  !*************************************************************!
  ! ------------------- PHYSICAL CONSTANTS ---------------------!
  ! ************************************************************!
  !                                                             !
  ! ------------------------------------------------------------!
  ! Time conversion factors
    REAL(prec), PARAMETER   :: secondsToMinutes = 1.0_prec/60.0_prec                   ! conversion for seconds to minutes
    REAL(prec), PARAMETER   :: minutesToHours   = 1.0_prec/60.0_prec                   ! conversion for minutes to hours
    REAL(prec), PARAMETER   :: hoursToDays      = 1.0_prec/24.0_prec                   ! conversion for hours to days
    REAL(prec), PARAMETER   :: daysToMonths     = 12.0_prec/365.25_prec                ! conversion for days to months
    REAL(prec), PARAMETER   :: monthsToYears    = 1.0_prec/12.0_prec                   ! conversion for months to years
    REAL(prec), PARAMETER   :: daysToSeconds    = 86400.0_prec

  !*************************************************************!
  ! ------------------- Package defaults ---------------------!
  ! ************************************************************!
  !                                                             !
  ! ------------------------------------------------------------!
    INTEGER, PARAMETER :: defaultNameLength = 40

 END MODULE ConstantsDictionary
