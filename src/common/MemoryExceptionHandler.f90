MODULE MemoryExceptionHandler
!
! MemoryExceptionHandler.f90
!
! Joe Schoonover
!
! schoonover.numerics@gmail.com
!
!
! This MODULE provides subroutines which process the status outputs of allocate and deallocate
! statements. During this process, a user-friendly message is printed to the screen and an
! action is decided upon.
!
! =========================================================================================== !
!
!
!

 USE ConstantsDictionary


 CONTAINS

 SUBROUTINE CheckForAllocationException( allocStat )
 ! S/R CheckForAllocationException
 !
 ! When an allocation statement is called, an optional "STAT" argument
 ! can be assigned to an integer. This integer, "allocStat" can be 
 ! sent to this subroutine for processing. Based on the value assigned
 ! to allocStat, a message is printed to the screen if there is an issue
 ! with memory allocation. Additionally, a call to the MainExceptionHandler
 ! is made to determine what the software should do.
 !
 ! ===================================================================== !
   IMPLICIT NONE
   INTEGER, INTENT(IN) :: allocStat
   
      
   
   
 END SUBROUTINE CheckForAllocationException


END MODULE MemoryExceptionHandler
