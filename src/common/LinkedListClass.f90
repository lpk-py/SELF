! LinkedListClass.f90
! 
! Copyright 2015 Joseph Schoonover <schoonover.numerics@gmail.com>
! 
! LinkedListClass.f90 is part of the Spectral Element Libraries in Fortran (SELF).
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
 
 

MODULE LinkedListClass
! ========================================= Logs ================================================= !
!2015-11-12  Joseph Schoonover  schoonover.numerics@gmail.com 
!  A linked-list data-structure with routines for construction, destruction, accession, and 
!  modification are provided. This provides a template for a linked list of other data-structures.
!  Here, the data is taken as an integer and we also provide an integer key which can be useful 
!  for accessing linked-list data for other data-structures.
!
!  The associated test program is found in ~/src/common/Testing/TestLinkedList.f90 . This program 
!  was tested with ValGrind ( on 12 Nov. 2015 ) to check for memory issues in this module. The 
!  output is given below :
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
! Which indicates that no memory leaks are possible and there are no other fatal errors regarding
! memory access.
! //////////////////////////////////////////////////////////////////////////////////////////////// ! 
 
 IMPLICIT NONE

   TYPE Record
      INTEGER                 :: listData
      INTEGER                 :: key
      TYPE( Record ), POINTER :: next

   END TYPE Record

   TYPE LinkedList      
      TYPE( Record ), POINTER :: head, tail, current

      CONTAINS

      PROCEDURE :: Build => Build_LinkedList
      PROCEDURE :: Trash => Trash_LinkedList
      
      
      PROCEDURE :: GetData => GetCurrentData_LinkedList
      PROCEDURE :: SetData => SetCurrentData_LinkedList
      PROCEDURE :: GetKey => GetCurrentKey_LinkedList
      PROCEDURE :: SetKey => SetCurrentKey_LinkedList
      
      PROCEDURE :: ListIsEmpty => ListIsEmpty_LinkedList
      PROCEDURE :: AddToList => AddToList_LinkedList
      PROCEDURE :: RemoveCurrent => RemoveCurrent_LinkedList

      PROCEDURE :: MoveToNext => MoveToNext_LinkedList
      PROCEDURE :: MoveToHead => MoveToHead_LinkedList
      PROCEDURE :: MoveToTail => MoveToTail_LinkedList
      PROCEDURE :: PrintToScreen => PrintToScreen_LinkedList

   END TYPE LinkedList


! Setting up some parameters pertaining to this module
 INTEGER, PARAMETER, PRIVATE :: keyInc   = 1 ! The default increment in the Record Key
 INTEGER, PARAMETER, PRIVATE :: keyStart = 1 ! The default starting Record key
 
 CONTAINS
 
! 
! ========================================================================== !
! ===================== Constructors and Destructors ======================= !
! ========================================================================== !
!
 SUBROUTINE Build_LinkedList( myList )
 ! S/R BUILD_LinkedList
 !
 !  This subroutine creates a linked-list by allocating space for the list
 !  head. The "next" attribute of the list head is nullified. 
 !  The tail is pointed to NULL and the current position in the list
 !  is pointed to the head. 
 ! 
 !   INPUT/OUTPUT : 
 !      CLASS( LinkedList) :: myList
 !
 ! ========================================================================== !
  IMPLICIT NONE
  CLASS( LinkedList ) :: myList

    
     myList % head => NULL( )

     ! Point the tail to null
     myList % tail => NULL()

     ! Set the current position to Null
     myList % current => NULL( )
  
 END SUBROUTINE Build_LinkedList
!
!
!  
 SUBROUTINE Trash_LinkedList( myList )
 ! S/R Trash_LinkedList
 ! 
 !  This subroutine traverses through the linked list and frees the memory
 !  that was allocated for the entries in the linked list.  
 !
 ! ========================================================================== !
  IMPLICIT NONE
  CLASS( LinkedList ) :: myList
  ! LOCAL
  TYPE( Record ), POINTER :: pNext

     ! Set the current position of the list to the head
     myList % current => myList % head
     
     ! Scroll through the list until the current position is nullified
     DO WHILE ( ASSOCIATED( myList % current ) )

        ! temporarily point to the next in the list
        pNext => myList % current % next 

        ! Deallocate memory pointed to by current position
        DEALLOCATE( myList % current ) 

        ! Update current position
        myList % current => pNext 

     ENDDO
  
 END SUBROUTINE Trash_LinkedList
! 
! ========================================================================== !
! =============================== Accessors ================================ !
! ========================================================================== !
!
 SUBROUTINE GetCurrentData_LinkedList( myList, outData )
 ! S/R GetCurrentData
 !
 !
 ! ========================================================================== !
 IMPLICIT NONE
  CLASS( LinkedList ) :: myList
  INTEGER             :: outData

     outData = myList % current % listData

 END SUBROUTINE GetCurrentData_LinkedList
!
!
!
 SUBROUTINE GetCurrentKey_LinkedList( myList, outKey )
 ! S/R GetCurrentKey
 !
 !
 ! ========================================================================== !
 IMPLICIT NONE
  CLASS( LinkedList ) :: myList
  INTEGER             :: outKey

     outKey = myList % current % key

 END SUBROUTINE GetCurrentKey_LinkedList
 !
 !
 !
 SUBROUTINE SetCurrentData_LinkedList( myList, inData )
 ! S/R SetCurrentData
 !
 !
 ! ========================================================================== !
 IMPLICIT NONE
  CLASS( LinkedList ) :: myList
  INTEGER             :: inData

     myList % current % listData = inData

 END SUBROUTINE SetCurrentData_LinkedList
!
!
!
 SUBROUTINE SetCurrentKey_LinkedList( myList, inKey )
 ! S/R SetCurrentKey
 !
 !
 ! ========================================================================== !
 IMPLICIT NONE
  CLASS( LinkedList ) :: myList
  INTEGER             :: inKey

     myList % current % key = inKey

 END SUBROUTINE SetCurrentKey_LinkedList
! 
! ========================================================================== !
! ======================= Data Structure Operations ======================== !
! ========================================================================== !
!
 FUNCTION ListIsEmpty_LinkedList( myList ) RESULT( TorF)
 ! Function ListIsEmpty
 !
 !  This function checks if "myList" is an empty list. A logical is returned
 !  
 ! ========================================================================== !
  IMPLICIT NONE
  CLASS( LinkedList ) :: myList
  LOGICAL             :: TorF

     TorF = .NOT.( ASSOCIATED( myList % head  ) )
     
 END FUNCTION ListIsEmpty_LinkedList
!
!
!
 SUBROUTINE AddToList_LinkedList( myList, inData, inKey )
 ! S/R AddToList
 !
 !  This subroutine adds an item to the linked list. If "inKey"
 !  is supplied, then this is filled in for the key data. Otherwise, the key is
 !  filled in as the previous records key plus one. The data is placed at the end
 !  of the linked list
 !  
 ! ========================================================================== !
  IMPLICIT NONE
  CLASS( LinkedList ) :: myList
  INTEGER             :: inData
  INTEGER, OPTIONAL   :: inKey
  ! LOCAL
  TYPE( Record ), POINTER :: previous
  INTEGER :: allocationStatus

     ! Check to see if this list is empty
     IF( myList % ListIsEmpty() )THEN
     
        ALLOCATE( myList % head, STAT = allocationStatus )
        IF( allocationStatus /=0 )THEN
           PRINT*, 'MODULE LinkedListClass.f90 : S/R AddToList : Memory not allocated for next entry in list.'
           ! An exception handler should be built to handle these problems
        ENDIF      
      
        ! Point the current position to the head
        myList % current => myList % head
        ! Set the data
        CALL myList % SetData( inData )
        
        IF( PRESENT(inKey) )THEN
           CALL myList % SetKey( inKey )
        ELSE
           CALL myList % SetKey( keyStart )
        ENDIF
        
        ! Point the next to null and the tail to current
        myList % current % next => NULL( )
        myList % tail => myList % current
        
     ELSE ! the list is not empty
    
        ! Then we allocate space for the next item in the list    
        ALLOCATE( myList % tail % next, STAT = allocationStatus )
        IF( allocationStatus /=0 )THEN
           PRINT*, 'MODULE LinkedListClass.f90 : S/R AddToList : Memory not allocated for next entry in list.'
           ! An exception handler should be built to handle these problems
        ENDIF      
        
        !  Temporarily hold onto the tail in case we need the key (if inKey is not present)  
        previous => myList % tail
        ! Reassign the tail
        myList % tail => myList % tail % next
        
        ! Set the current to the tail
        myList % current => myList % tail
  
        ! Fill in the data
        CALL myList % SetData( inData )
        
        ! Fill in the key information
        IF( PRESENT(inKey) )THEN
           CALL myList % SetKey( inKey )
        ELSE
           CALL myList % SetKey( previous % key + keyInc )
        ENDIF
        
        ! Point the next to null and the tail to current
        myList % current % next => NULL( )
        
     ENDIF

 END SUBROUTINE AddToList_LinkedList
!
!
!
 SUBROUTINE RemoveCurrent_LinkedList( myList )
 ! S/R RemoveCurrent
 !
 !  This subroutine removes the current item in the linked list and patches together the previous
 !  and next items.
 !  
 ! ===========================================================================================s==== !
  IMPLICIT NONE
  CLASS( LinkedList ) :: myList
  ! LOCAL
  TYPE( Record ), POINTER :: previous, pNext
  INTEGER                 :: currentKey, thisKey

     CALL myList % GetKey( currentKey )
     
     ! Check to see if this list is empty
     IF( myList % ListIsEmpty() )THEN
     
        PRINT*, 'Module LinkedListClass.f90 : S/R RemoveCurrent : List is empty. Nothing to remove'
        RETURN
        
     ELSE ! the list is not empty
    
        CALL myList % MoveToHead( ) ! Rewind the list
         
        ! Get the key for this list item
        CALL myList % GetKey( thisKey )
        
        ! Check if we are trying to remove the head of the list
        IF( thisKey == currentKey )THEN 
        
           ! temporarily point to the next in the list
           pNext => myList % current % next 
           
          ! Deallocate memory pointed to by current position
          DEALLOCATE( myList % current ) 

          ! Update current position
          myList % head => pNext ! Reset the head of the list
          
          RETURN
        ENDIF
        
        ! If the execution of the code has arrived here, then we are not removing the head of the 
        ! list. 
        ! Hang on to the head as the previous
        previous => myList % current
        CALL myList % MoveToNext( )
        
        DO WHILE( ASSOCIATED( myList % current ) )
        
           ! Get the key for this list item
           CALL myList % GetKey( thisKey )
           
           ! Check if we are trying to remove the head of the list
           IF( thisKey == currentKey )THEN 
           
              ! temporarily point to the next in the list
              pNext => myList % current % next 
            
              ! Patch the previous item to the next item
              previous % next => pNext
           
              IF( .NOT.ASSOCIATED(pNext)  )THEN
                myList % tail => previous
              ENDIF
              
             ! Deallocate memory pointed to by current position
             DEALLOCATE( myList % current ) 
          
             EXIT
           ELSE
           
              previous => myList % current
              CALL myList % moveToNext( )
              
           ENDIF
        
        ENDDO

     ENDIF

 END SUBROUTINE RemoveCurrent_LinkedList
!
!
!
 SUBROUTINE MoveToNext_LinkedList( myList )
 ! S/R MoveToNext
 !
 !  
 ! ========================================================================== !
  IMPLICIT NONE
  CLASS( LinkedList ) :: myList

     myList % current => myList % current % next

 END SUBROUTINE MoveToNext_LinkedList
!
!
!
 SUBROUTINE MoveToHead_LinkedList( myList )
 ! S/R MoveToHead
 !
 !  
 ! ========================================================================== !
  IMPLICIT NONE
  CLASS( LinkedList ) :: myList

     myList % current => myList % head

 END SUBROUTINE MoveToHead_LinkedList
!
!
!
 SUBROUTINE MoveToTail_LinkedList( myList )
 ! S/R MoveToTail
 !
 !  
 ! ========================================================================== !
  IMPLICIT NONE
  CLASS( LinkedList ) :: myList

     myList % current => myList % tail

 END SUBROUTINE MoveToTail_LinkedList
!
!
!
 
 SUBROUTINE GetCount_LinkedList( myList, numberOfListItems )
 ! S/R GetCount
 !
 !
 ! ===========================================================================================s==== !
  IMPLICIT NONE
  CLASS( LinkedList )  :: myList
  INTEGER, INTENT(out) :: numberOfListItems

     numberOfListItems = 0 ! Initialize the number of list items
     ! Check to see if this list is empty
     IF( myList % ListIsEmpty() )THEN
     
        PRINT*, 'Module LinkedListClass.f90 : S/R ListCount : List is empty.'
        RETURN
        
     ELSE ! the list is not empty
    
        CALL myList % MoveToHead( ) ! Rewind the list
        
        DO WHILE( ASSOCIATED( myList % current ) )
        
           numberOfListItems = numberOfListItems + 1
           CALL myList % moveToNext( )

        ENDDO

     ENDIF

 END SUBROUTINE GetCount_LinkedList
!
!
!
 SUBROUTINE PrintToScreen_LinkedList( myList )
 ! S/R PrintToScreen
 !
 !  
 ! ========================================================================== !
  IMPLICIT NONE
  CLASS( LinkedList ) :: myList

     myList % current => myList % head

        PRINT*, '          Data        Key'
     DO WHILE( ASSOCIATED( myList % current ) )
 
        PRINT*, myList % current % listData, myList % current % key

        CALL myList % MoveToNext()

     ENDDO

 END SUBROUTINE PrintToScreen_LinkedList
!
!
!
END MODULE LinkedListClass
