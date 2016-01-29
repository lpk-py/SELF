MODULE VectorClass
! VectorClass.f90 (v2.1 - 14 Dec. 2015) 
! 
! schoonover.numerics@gmail.com
! 
! o (ver 1.0) March 2014
! o (ver 2.1) Dec 2015
!
!
! 
!  
! ================================================================================================ !
USE ModelPrecision
USE ConstantsDictionary


IMPLICIT NONE

   TYPE Vector
       REAL(prec),ALLOCATABLE :: dir(:)
       REAL(prec)             :: length

       CONTAINS 

       PROCEDURE :: Build => Build_Vector
       PROCEDURE :: Trash => Trash_Vector

       ! ACCESSORS
       PROCEDURE :: SetDirection => SetDirection_Vector
       PROCEDURE :: GetDirection => GetDirection_Vector
       PROCEDURE :: SetLength => SetLength_Vector
       PROCEDURE :: GetLength => GetLength_Vector
       PROCEDURE :: GetDimension => GetDimension_Vector

       PROCEDURE :: Add => Add_Vector
       PROCEDURE :: Dot => Dot_Vector
       PROCEDURE :: Normalize => Normalize_Vector
       
   END TYPE Vector

   !INTERFACE OPERATOR(+)
   !   MODULE PROCEDURE Add_Vector
   !END INTERFACE

   !INTERFACE OPERATOR(*)
   !   MODULE PROCEDURE Dot_Vector
   !END INTERFACE

 CONTAINS
!
!
!==================================================================================================!
!------------------------------- Manual Constructors/Destructors ----------------------------------!
!==================================================================================================!
!
!
 SUBROUTINE Build_Vector( myVec, nDim, dir  )
 ! S/R Build
 !
 ! =============================================================================================== !
 ! DECLARATIONS
  IMPLICIT NONE
  CLASS(Vector), INTENT(out)       :: myVec
  INTEGER, INTENT(in)              :: nDim
  REAL(prec), INTENT(in), OPTIONAL :: dir(1:nDim) 

      ALLOCATE( myVec % dir(1:nDim) )
      myVec % dir = ZERO
      
      IF( PRESENT(dir) )THEN
         CALL myVec % SetDirection( dir )
      ENDIF
      
      CALL myVec % SetLength( DOT_PRODUCT( myvec % dir, myvec % dir ) )

 END SUBROUTINE Build_Vector
!
!
!
 SUBROUTINE Trash_Vector( myvec )
 ! S/R Trash
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS(Vector), INTENT(inout) :: myvec
   
      DEALLOCATE( myVec % dir )

 END SUBROUTINE Trash_Vector
!
!
!==================================================================================================!
!--------------------------------------- Accessors ------------------------------------------------!
!==================================================================================================!
!
!
 SUBROUTINE SetDirection_Vector( myvec, direc  )
 ! S/R SetDirection
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS(Vector), INTENT(inout) :: myvec
   REAL(prec), INTENT(in)       :: direc(1:)
   
      myvec % dir = direc

 END SUBROUTINE SetDirection_Vector
!
!
!
 SUBROUTINE GetDirection_Vector( myvec, direc )
 ! S/R GetDirection
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS(Vector), INTENT(in) :: myvec
   REAL(prec), INTENT(out)   :: direc(1:)
   
      direc = myvec % dir

 END SUBROUTINE GetDirection_Vector
!
!
!
 SUBROUTINE SetLength_Vector( myvec, leng )
 ! S/R SetLength
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS(Vector), INTENT(inout) :: myvec
   REAL(prec), INTENT(in)       :: leng

      myvec % length = leng

 END SUBROUTINE SetLength_Vector
!
!
!
 SUBROUTINE GetLength_Vector( myvec, leng )
 ! S/R GetLength
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS(Vector), INTENT(in) :: myvec
   REAL(prec), INTENT(out)   :: leng

      leng = myvec % length

 END SUBROUTINE GetLength_Vector
!
!
!
 SUBROUTINE GetDimension_Vector( myvec, dimen  )
 ! S/R GetDimension
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS(Vector), INTENT(in) :: myvec
   INTEGER, INTENT(out)      :: dimen

      dimen = SIZE( myvec % dir, 1 )

 END SUBROUTINE GetDimension_Vector
!
!
!==================================================================================================!
!--------------------------------- Type Specific Routines -----------------------------------------!
!==================================================================================================!
!
!
 FUNCTION Add_Vector( vec1, vec2 ) RESULT ( vec3 )
 ! FUNCTION Add_Vectors
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS(Vector), INTENT(in) :: vec1
   TYPE(Vector)               :: vec2
   TYPE(Vector)               :: vec3
   ! LOCAL
   INTEGER :: iDimen
   
      ! Currently no check is in place to ensure that all three vectors are
      ! the same dimension

      do iDimen = 1, SIZE(vec1 % dir, 1) ! Cycle over the vector dimensions
  
         vec3 % dir(iDimen) = vec1 % dir(iDimen) + vec2 % dir(iDimen) 

      enddo  ! iDim, cycle over the vector dimensions

      CALL vec3 % SetLength( sqrt( vec3 % Dot(vec3) ) )

 END FUNCTION Add_Vector
!
!
!
 FUNCTION Dot_Vector( vec1, vec2 ) RESULT ( vec1DotVec2 )
 ! FUNCTION Dot_Vector
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS(Vector), INTENT(in) :: vec1, vec2
   REAL(prec)               :: vec1DotVec2
   ! LOCAL
   INTEGER :: iDimen
   
      vec1DotVec2 = DOT_PRODUCT( vec1 % dir, vec2 % dir )

 END FUNCTION Dot_Vector
!
!
! 
 SUBROUTINE Normalize_Vector( myvec ) 
 ! S/R Normalize
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS(Vector), INTENT(inout)  :: myvec
   
      myvec % dir = myvec % dir/sqrt( Dot_Vector(myvec, myvec) )

 END SUBROUTINE Normalize_Vector

END MODULE VectorClass
