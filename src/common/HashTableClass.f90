MODULE HashTableClass
! A hash-table class is built using an array of pointers to linked lists
!
 
 USE ConstantsDictionary
 USE LinkedListClass

 IMPLICIT NONE


   TYPE HashTable      
      TYPE( LinkedList ), ALLOCATABLE :: list(:)

      CONTAINS

      PROCEDURE :: Build => Build_HashTable
      PROCEDURE :: Trash => Trash_HashTable
      PROCEDURE :: AddDataForKeys => AddDataForKeys_HashTable
      PROCEDURE :: ContainsKeys => ContainsKeys_HashTable
      PROCEDURE :: GetDataForKeys => GetDataForKeys_HashTable

   END TYPE HashTable

   

 CONTAINS
!
!
!==================================================================================================!
!------------------------------- Manual Constructors/Destructors ----------------------------------!
!==================================================================================================!
!
!
 SUBROUTINE Build_HashTable( myTable, N )
 ! S/R Build
 !
 !  
 ! =============================================================================================== !
  IMPLICIT NONE
  CLASS( HashTable ) :: myTable
  INTEGER            :: N
  ! LOCAL
  INTEGER :: i


     ALLOCATE( myTable % list(1:N) )

     DO i = 1,N
        CALL myTable % list(i) % Build( )
     ENDDO

 END SUBROUTINE Build_HashTable
!
!
!  
 SUBROUTINE Trash_HashTable( myTable )
 ! S/R Trash
 !
 !  
 ! =============================================================================================== !
  IMPLICIT NONE
  CLASS( HashTable ) :: myTable
  ! LOCAL
  INTEGER :: i, N


     N = SIZE( myTable % list )

     DO i = 1, N
        CALL myTable % list(i) % Trash( )
     ENDDO

     DEALLOCATE( myTable % list )
  
 END SUBROUTINE Trash_HashTable
!
!
!==================================================================================================!
!--------------------------------- Type Specific Routines  ----------------------------------------!
!==================================================================================================!
!
!
 SUBROUTINE AddDataForKeys_HashTable( myTable, inData, i, j )
 ! S/R AddDataForKeys
 !
 ! =============================================================================================== !
  IMPLICIT NONE
  CLASS( HashTable ) :: myTable
  INTEGER            :: inData, i, j
  
     CALL myTable % list(i) % AddToList( inData, j ) 

 END SUBROUTINE AddDataForKeys_HashTable
!
!
!
 FUNCTION ContainsKeys_HashTable( myTable, i, j ) RESULT( doesContain )
 ! FUNCTION ContainsKeys
 !
 ! =============================================================================================== !
  IMPLICIT NONE
  CLASS( HashTable ) :: myTable
  INTEGER            :: i, j
  LOGICAL            :: doesContain
  ! LOCAL
  INTEGER :: thisKey

     IF(  myTable % list(i) % ListIsEmpty() )THEN ! this list hasn't been started
        doesContain = .FALSE.
        RETURN
     ENDIF
 
     ! Rewind the list
     myTable % list(i) % current => myTable % list(i) % head

     DO WHILE( ASSOCIATED( myTable % list(i) % current ) )
    
        CALL myTable % list(i) % GetKey( thisKey )

        IF( thisKey == j )THEN ! This list already has this key
           doesContain =.TRUE.
           RETURN
        ENDIF

        ! other wise we move to the next element in the list
        CALL myTable % list(i) % MoveToNext( )

     ENDDO
 

 END FUNCTION ContainsKeys_HashTable
!
!
!
  SUBROUTINE GetDataForKeys_HashTable( myTable, outData, i, j )
 ! S/R GetDataForKeys
 !
 ! =============================================================================================== !
  CLASS( HashTable ) :: myTable
  INTEGER            :: outData, i, j
  ! LOCAL
  INTEGER :: thisData, thisKey
  
     IF(  myTable % list(i) % ListIsEmpty() )THEN ! this table entry is not pointing to a linked list
        PRINT*, 'MODULE HASHTABLE_CLASS : S/R GetDataForKeys :'
        PRINT*, 'List ', i,' is empty.'
        outData = fillValueInt
        RETURN
     ENDIF

     myTable % list(i) % current => myTable % list(i) % head

     DO WHILE( ASSOCIATED( myTable % list(i) % current ) ) ! this table entry does not contain the keys i,j

        ! Add to the list
        CALL myTable % list(i) % GetData( thisData )
        CALL myTable % list(i) % GetKey( thisKey )

        IF( thisKey == j )THEN
           outData = thisData
           RETURN
        ENDIF
        
        CALL myTable % list(i) % MoveToNext( )
      
     ENDDO
   
 END SUBROUTINE GetDataForKeys_HashTable

END MODULE HashTableClass
