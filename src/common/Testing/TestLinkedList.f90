PROGRAM TestLinkedList
!
! TestLinkedList.f90
!
! schoonover.numerics@gmail.com
!
! This program tests the functionality of the LinkedListClass in a very simple example.
! The list is loaded up with integers from 1 to 50, and is printed to screen. 
!
! ************************************************************************************************ !
! On Nov. 12,2015, this program was tested with valgrind, from which the following output is given:
!
!==2154== 
!==2154== HEAP SUMMARY:
!==2154==     in use at exit: 0 bytes in 0 blocks
!==2154==   total heap usage: 71 allocs, 71 frees, 12,644 bytes allocated
!==2154== 
!==2154== All heap blocks were freed -- no leaks are possible
!==2154== 
!==2154== For counts of detected and suppressed errors, rerun with: -v
!==2154== ERROR SUMMARY: 0 errors from 0 contexts (suppressed: 0 from 0)
!
!
! This output indicates that there are no memory leaks.
! ************************************************************************************************ !
!
!
  USE LinkedListClass

  IMPLICIT NONE

  TYPE( LinkedList ) :: thisList
  
  integer :: i


     CALL thisList % Build( )


     do i = 1, 50

        CALL thisList % AddToList( i )
   
     enddo

     CALL thisList % PrintToScreen( )

     CALL thisList % MoveToHead( )
     DO i = 1,10
        CALL thisList % MoveToNext( )
     ENDDO
     
     CALL thisList % RemoveCurrent( )
     
     PRINT*, 'Removed 11th item....'
     
     CALL thisList % PrintToScreen( )
     
     CALL thisList % Trash( )

 
END PROGRAM TestLinkedList
