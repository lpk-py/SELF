MODULE CommonRoutines
! 
! CommonRoutines.f90 (v2.1 - 12 Nov. 2015)
!
! schoonover.numerics@gmail.com
!
! This module provides routines for various operations that are broadly useful. It is recognized
! that there probably should be a more coherent structure to this end of the software, but this 
! has not be identified yet.
!
!
! ================================================================================================= !

USE ModelPrecision
USE ConstantsDictionary
USE ModelFlags

IMPLICIT NONE


   INTERFACE Invert
      MODULE PROCEDURE :: Invert_2x2
   END INTERFACE

CONTAINS


 FUNCTION AlmostEqual( a, b ) RESULT( AisB )
 !
 ! This FUNCTION takes in two REAL(prec) values (a and b) and
 ! returns a logical which tells if a is approximately b in 
 ! floating point aritmetic.
 !
 ! 
 ! Input : 
 !       REAL(prec) :: a
 !       REAL(prec) :: b
 ! 
 ! Output :
 !       logical :: AisB
 !
 !=============================================================
  IMPLICIT NONE
  REAL(prec) :: a, b
  logical :: AisB


    if( a == ZERO .OR. b == ZERO )then

       if( abs(a-b) <= epsilon(ONE) )then 
          AisB = .TRUE.
       else
          AisB = .FALSE.
       ENDif

    else

       if( (abs(a-b) <= epsilon(ONE)*abs(a)) .OR. (abs(a-b) <= epsilon(ONE)*abs(b)) )then
          AisB = .TRUE.
       else
          AisB = .FALSE.
       ENDif

    ENDif

 END FUNCTION AlmostEqual
!
!
!
 SUBROUTINE SortAndSum( myArray, low, high, arraysum)
 !
 ! 
 !
 !=============================================================
 INTEGER, INTENT(in)       :: low, high
 REAL(prec), INTENT(inout) :: myArray(low:high)
 REAL(prec), INTENT(out)   :: arraysum
 
 ! LOCAL
 INTEGER :: ind


    CALL SortArray( myArray, low, high )

    arraysum = ZERO

    DO ind = low, high

       arraysum = arraysum + myArray(ind)

    ENDDO
    

 END SUBROUTINE SortAndSum
!
!
!
 SUBROUTINE SortArray( myArray, low, high )
 !
 ! 
 !       
 !
 !=============================================================
 INTEGER, INTENT(in)       :: low, high
 REAL(prec), INTENT(inout) :: myArray(low:high)
 
 ! LOCAL
 REAL(prec) :: temp
 INTEGER :: locOfMin
 INTEGER :: ind


    DO ind = low, high-1

       locOfMin = MINLOC( abs(myArray(ind:high)),1 ) + low - 1 + ind

       CALL SwapIndices( myArray, low, high, ind, locOfMin )

    ENDDO
    

 END SUBROUTINE SortArray
!
!
!
 SUBROUTINE SwapIndices( myArray, low, high, ind1, ind2 )
 !=============================================================
 ! 
 !
 ! 
 ! Input : 
 !       
 ! 
 ! Output :
 !       
 !
 !=============================================================
 INTEGER, INTENT(in)       :: low, high
 REAL(prec), INTENT(inout) :: myArray(low:high)
 INTEGER, INTENT(in)       :: ind1, ind2
 ! LOCAL
 REAL(prec) :: temp


    temp = myArray(ind1)

    myArray(ind1) = myArray(ind2)

    myArray(ind2) = temp


 END SUBROUTINE SwapIndices
!
!
!
 SUBROUTINE ReverseArray( myArray, low, high )
 !=============================================================
 ! 
 !
 ! 
 ! Input : 
 !       
 ! 
 ! Output :
 !       
 !
 !=============================================================
 INTEGER, INTENT(in)       :: low, high
 REAL(prec), INTENT(inout) :: myArray(low:high)
 ! LOCAL
 REAL(prec) :: temp(low:high)
 INTEGER    :: i, j


    temp = myArray
    j = high

    DO i = low,high

       myArray(i) = temp(j)
       j = j-1

    ENDDO


 END SUBROUTINE ReverseArray
!
!
!
 INTEGER FUNCTION NewUnit(thisunit)
 ! This FUNCTION searches for a new FILE UNIT that is not currently in USE.
 !
 !=========================================================================!
 ! DECLARATIONS
  INTEGER, INTENT(out), optional :: thisunit
 ! LOCAL
  INTEGER, PARAMETER :: unitMin=100, unitMax=1000
  logical :: isopened
  INTEGER :: iUnit


     newunit=-1

     DO iUnit=unitMin, unitMax

        ! Check to see if this UNIT is opened
        INQUIRE(UNIT=iUnit,opened=isopened)

        if( .not. isopened )then

           newunit=iUnit
           EXIT

        ENDif

     ENDDO

     if (PRESENT(thisunit)) thisunit = newunit
 
END FUNCTION NewUnit
!
!
!
 FUNCTION GetFlagForChar( mychar ) RESULT( outFlag )
!
!
!
!================================================================ !
! DECLARATIONS
  CHARACTER(*) :: myChar
  INTEGER      :: outFlag

     if( UpperCase(trim(myChar)) == 'LINEAR' )then
 
         outFlag = LINEAR

     elseif( UpperCase(trim(myChar)) == 'NONLINEAR' )then
 
         outFlag = NONLINEAR

     elseif( UpperCase(trim(myChar)) == 'SKEW-SYMMETRIC' )then
 
         outFlag = SKEW_SYMMETRIC

     elseif( UpperCase(trim(myChar)) == 'CONSERVATIVE')then
 
         outFlag = CONSERVATIVE

     elseif( UpperCase(trim(myChar)) == 'GAUSS' )then
 
         outFlag = GAUSS

     elseif( UpperCase(trim(myChar)) == 'GAUSS LOBATTO'  )then
 
         outFlag = CONSERVATIVE
     ENDif

 END FUNCTION GetFlagForChar
!
!
!
 FUNCTION UpperCase ( myChar ) Result (myUpperChar)
 !
 !   Changes a string to all upper case letters
 ! 
 ! ============================================================== !
    IMPLICIT NONE
    CHARACTER(*)           :: myChar
    ! LOCAL
    CHARACTER(LEN(myChar))   :: myUpperChar
    INTEGER                  :: ic, i
    CHARACTER(28), PARAMETER :: cap = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ- '
    CHARACTER(28), PARAMETER :: low = 'abcdefghijklmnopqrstuvwxyz- '


       !   Capitalize each letter if it is lowecase
       myUpperChar = myChar

       DO i = 1, LEN_TRIM(myChar)

           ic = INDEX(low, myChar(i:i))
           if (ic > 0) myUpperChar(i:i) = cap(ic:ic)

       ENDDO

END FUNCTION UpperCase
!
!
!
 FUNCTION UniformPoints( a, b, N ) RESULT( xU )
 !
 !
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   REAL(prec) :: a, b
   INTEGER    :: N
   REAL(prec) :: xU(0:N)
   ! LOCAL
   REAL    :: dx
   INTEGER :: i
   
      dx = (b-a)/REAL(N,prec)
   
      DO i = 0,N
      
         xU(i) = a + dx*REAL(i,prec)
      
      ENDDO   
   
 END FUNCTION UniformPoints
!
!
!
  RECURSIVE FUNCTION Determinant( A, N ) RESULT( D )
 ! FUNCTION Determinant
 !
 !  Computes the determinant of the N X N matrix, A.
 !
 ! =================================================
 ! DECLARATIONS
   integer    :: N
   real(prec) :: A(1:N,1:N)
   real(prec) :: D
   ! LOCAL
   real(prec) :: M(1:N-1,1:N-1)
   integer    :: j


      if(  N == 2 )then
     !    print*, 'N=',N
         D = A(1,1)*A(2,2) - A(1,2)*A(2,1)
     !    print*, D
         RETURN

      else
        ! print*, 'N=',N
         D = ZERO

         do j = 1, N

            M = GetMinor( A, 1, j, N )
        !    print*, M(1,:)
        !    print*, M(2,:)

            D = D + A(1,j)*Determinant( M, N-1 )

         enddo

      endif
      
 
 END FUNCTION Determinant
!
!
!
 FUNCTION GetMinor( A, i, j, N ) RESULT( M )
 ! FUNCTION GetMinor
 !
 !  Returns the minor matrix of the N X N matrix A.
 !  The minor matrix of A is an (N-1)X(N-1) matrix.
 !
 ! =================================================
 ! DECLARATIONS
   integer    :: i, j, N
   real(prec) :: A(1:N,1:N)
   real(prec) :: M(1:N-1,1:N-1)
   ! LOCAL
   integer :: row, col
   integer :: thisRow, thisCol
      

      thisRow = 0
      

      do row = 1, N ! loop over the rows of A
  
         if( row /= i ) then

            thisRow = thisRow + 1      

            thisCol = 0
            do col = 1, N ! loop over the columns of A

               if( col /= j ) then

                  thisCol = thisCol + 1

                  M(thisRow,thisCol) = (-1.0_prec)**(i+j)*A(row, col)

               endif

            enddo ! col, loop over the columns of A

         endif 

      enddo ! row, loop over the rows of A


 END FUNCTION GetMinor
!
!
!
 FUNCTION Invert_2x2( A ) RESULT( Ainv )
 !
 ! =============================================================================================== !
   IMPLICIT NONE
   REAL(prec) :: A(1:2,1:2)
   REAL(prec) :: Ainv(1:2,1:2)
   ! LOCAL
   REAL(prec) :: detA
   
      detA = Determinant( A, 2 )
      
      Ainv(1,1) = A(2,2)/detA
      Ainv(2,2) = A(1,1)/detA
      Ainv(1,2) = -A(1,2)/detA
      Ainv(2,1) = -A(2,1)/detA

 END FUNCTION Invert_2x2
!
!
!
 FUNCTION Invert_3x3( A ) RESULT( Ainv )
 !
 ! =============================================================================================== !
   IMPLICIT NONE
   REAL(prec) :: A(1:3,1:3)
   REAL(prec) :: Ainv(1:3,1:3)
   ! LOCAL
   REAL(prec) :: detA
   REAL(prec) :: submat(1:2,1:2)
   REAL(prec) :: detSubmat
   
      detA = Determinant( A, 3 )
      
      ! Row 1 column 1 of inverse (use submatrix neglecting row 1 and column 1 of A)
      submat    =  A(2:3,2:3)
      detSubmat = Determinant( submat, 2 )
      Ainv(1,1) = detSubmat/detA
      
      ! Row 1 column 2 of inverse (use submatrix neglecting row 2 and column 1 of A)
      submat    =  A(1:3:2,2:3)
      detSubmat = Determinant( submat, 2 )
      Ainv(1,2) = -detSubmat/detA
      
      ! Row 1 column 3 of inverse (use submatrix neglecting row 3 and column 1 of A)
      submat    =  A(1:2,2:3)
      detSubmat = Determinant( submat, 2 )
      Ainv(1,3) = detSubmat/detA

      ! Row 2 column 1 of inverse (use submatrix neglecting row 1 and column 2 of A)
      submat    =  A(2:3,1:3:2)
      detSubmat = Determinant( submat, 2 )
      Ainv(2,1) = -detSubmat/detA
      
      ! Row 2 column 2 of inverse (use submatrix neglecting row 2 and column 2 of A)
      submat    =  A(1:3:2,1:3:2)
      detSubmat = Determinant( submat, 2 )
      Ainv(2,2) = -detSubmat/detA
      
      ! Row 2 column 3 of inverse (use submatrix neglecting row 3 and column 2 of A)
      submat    =  A(1:2,1:3:2)
      detSubmat = Determinant( submat, 2 )
      Ainv(2,3) = -detSubmat/detA
      
      ! Row 3 column 1 of inverse (use submatrix neglecting row 1 and column 3 of A)
      submat    =  A(2:3,1:2)
      detSubmat = Determinant( submat, 2 )
      Ainv(3,1) = detSubmat/detA
      
      ! Row 3 column 2 of inverse (use submatrix neglecting row 2 and column 3 of A)
      submat    =  A(1:3:2,1:2)
      detSubmat = Determinant( submat, 2 )
      Ainv(3,2) = -detSubmat/detA
      
      ! Row 3 column 3 of inverse (use submatrix neglecting row 3 and column 3 of A)
      submat    =  A(1:2,1:2)
      detSubmat = Determinant( submat, 2 )
      Ainv(3,3) = detSubmat/detA

 END FUNCTION Invert_3x3 
 
END MODULE COMMONROUTINES
