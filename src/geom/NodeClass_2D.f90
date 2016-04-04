MODULE NodeClass_2D
! NodeClass_2D.f90
! 
! Module History 
! 
! o (ver 1.0)  April 2015
! o (ver 2.1)  Dec. 2015
! 
!
! This module defines the mesh primitive "node". The node contains a node type information, the 
! position, and a node-to-element connectivity list. Linked-List type structure and routines are
! provided for those applications which require dynamic memory allocation in which the number of 
! nodes is unknown a'priori (e.g. Mesh Generation). One can use the Node_2D data-structure on its own
! without the linked-list structure in cases where the number of nodes is known.
!
! =======================================================================================

USE ModelPrecision
USE ModelFlags
USE ConstantsDictionary
USE LinkedListClass


IMPLICIT NONE

   TYPE Node_2D
      INTEGER, PRIVATE           :: key
      INTEGER, PRIVATE           :: nodeType      ! An INTEGER flag for INTERIOR or BOUNDARY
      REAL(prec), PRIVATE        :: x, y
      TYPE( LinkedList ), PUBLIC :: nodeToElement ! A linked list of element IDs which share this node.
      TYPE( Node_2D ), POINTER   :: next

      CONTAINS

      PROCEDURE :: Build => Build_Node_2D
      PROCEDURE :: Trash => Trash_Node_2D

      PROCEDURE :: GetData => GetData_Node_2D
      PROCEDURE :: SetData => SetData_Node_2D
      PROCEDURE :: GetKey => GetKey_Node_2D
      PROCEDURE :: SetKey => SetKey_Node_2D
      PROCEDURE :: GetType => GetType_Node_2D
      PROCEDURE :: SetType => SetType_Node_2D
      PROCEDURE :: GetPosition => GetPosition_Node_2D
      PROCEDURE :: SetPosition => SetPosition_Node_2D
      
      PROCEDURE :: ScaleNodePosition => ScaleNodePosition_Node_2D
      
      
   END TYPE Node_2D

   TYPE NodeList_2D      
      TYPE( Node_2D ), POINTER :: head, tail, current

      CONTAINS

      PROCEDURE :: Build => Build_Node_2DList
      PROCEDURE :: Trash => Trash_Node_2DList
      
      PROCEDURE :: GetData => GetCurrentData_Node_2DList
      PROCEDURE :: SetData => SetCurrentData_Node_2DList
      PROCEDURE :: GetPosition => GetCurrentPosition_Node_2DList
      PROCEDURE :: SetPosition => SetCurrentPosition_Node_2DList
      PROCEDURE :: GetKey => GetCurrentKey_Node_2DList
      PROCEDURE :: SetKey => SetCurrentKey_Node_2DList
      
      PROCEDURE :: ListIsEmpty => ListIsEmpty_Node_2DList
      PROCEDURE :: AddToList => AddToList_Node_2DList
      PROCEDURE :: RemoveCurrent => RemoveCurrent_Node_2DList
      PROCEDURE :: MoveToHead => MoveToHead_Node_2DList
      PROCEDURE :: MoveToNext => MoveToNext_Node_2DList
      PROCEDURE :: MoveToTail => MoveToTail_Node_2DList
      PROCEDURE :: GetCount => GetCount_Node_2DList

   END TYPE NodeList_2D

! Setting up some parameters pertaining to this module
 INTEGER, PARAMETER, PRIVATE :: keyInc   = 1 ! The default increment in the Record Key
 INTEGER, PARAMETER, PRIVATE :: keyStart = 1 ! The default starting Record key

 CONTAINS
!
!
!==================================================================================================!
!------------------------------- Manual Constructors/Destructors ----------------------------------!
!==================================================================================================!
!
!
 SUBROUTINE Build_Node_2D( thisNode, x, y )
 ! S/R Build_Node_2D
 !  
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS( Node_2D ), INTENT(inout) :: thisNode
   REAL(prec), INTENT(in)       :: x, y

 
      CALL thisNode % SetPosition( x, y )
      
      CALL thisNode % SetType( INTERIOR )
      
      CALL thisNode % nodeToElement % Build( )
      
      thisNode % next => NULL(  )
      
 END SUBROUTINE Build_Node_2D
!
!
!
 SUBROUTINE Trash_Node_2D( thisNode )
 ! S/R Trash_Node_2D
 !  
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS( Node_2D ), INTENT(inout) :: thisNode

      CALL thisNode % nodeToElement % Trash( )

 END SUBROUTINE Trash_Node_2D
!
!
!
 SUBROUTINE Build_Node_2DList( myList )
 ! S/R BUILD_Node_2DList
 !
 !  This subroutine creates a linked-list by allocating space for the list
 !  head. The "next" attribute of the list head is nullified. 
 !  The tail is pointed to NULL and the current position in the list
 !  is pointed to the head. 
 ! 
 !   INPUT/OUTPUT : 
 !      CLASS( NodeList_2D) :: myList
 !
 ! =============================================================================================== !
 ! DECLARATIONS
  IMPLICIT NONE
  CLASS( NodeList_2D ) :: myList

    
     myList % head => NULL( )

     ! Point the tail to null
     myList % tail => NULL()

     ! Set the current position to Null
     myList % current => NULL( )
  
 END SUBROUTINE Build_Node_2DList
!
!
!  
 SUBROUTINE Trash_Node_2DList( myList )
 ! S/R Trash_Node_2DList
 ! 
 !  This subroutine traverses through the linked list and frees the memory
 !  that was allocated for the entries in the linked list.  
 !
 ! =============================================================================================== !
 ! DECLARATIONS
  IMPLICIT NONE
  CLASS( NodeList_2D ) :: myList
  ! LOCAL
  TYPE( Node_2D ), POINTER :: pNext

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
  
 END SUBROUTINE Trash_Node_2DList
!
!
!==================================================================================================!
!--------------------------------------- Accessors ------------------------------------------------!
!==================================================================================================!
!
!
 SUBROUTINE SetCurrentData_Node_2DList( myList, x, y, nodeType )
 ! S/R SetCurrentData
 !
 !
 ! ========================================================================== !
 IMPLICIT NONE
  CLASS( NodeList_2D )      :: myList
  REAL(prec), INTENT(in) :: x, y
  INTEGER, INTENT(in)    :: nodeType

     CALL myList % current % SetData( x,  y, nodeType )

 END SUBROUTINE SetCurrentData_Node_2DList
!
!
!
 SUBROUTINE GetCurrentData_Node_2DList( myList, x, y, nodeType )
 ! S/R GetCurrentData
 !
 !
 ! ========================================================================== !
 IMPLICIT NONE
  CLASS( NodeList_2D )       :: myList
  REAL(prec), INTENT(out) :: x, y
  INTEGER, INTENT(out)    :: nodeType

     CALL myList % current % GetData( x, y, nodeType )

 END SUBROUTINE GetCurrentData_Node_2DList
!
!
!
 SUBROUTINE SetCurrentPosition_Node_2DList( myList, x, y )
 ! S/R SetCurrentPosition
 !
 !
 ! ========================================================================== !
 IMPLICIT NONE
  CLASS( NodeList_2D )      :: myList
  REAL(prec), INTENT(in) :: x, y

     CALL myList % current % SetPosition( x, y )

 END SUBROUTINE SetCurrentPosition_Node_2DList
!
!
!
 SUBROUTINE GetCurrentPosition_Node_2DList( myList, x, y )
 ! S/R GetCurrentPosition
 !
 !
 ! ========================================================================== !
 IMPLICIT NONE
  CLASS( NodeList_2D )       :: myList
  REAL(prec), INTENT(out) :: x, y

     CALL myList % current % GetPosition( x, y )

 END SUBROUTINE GetCurrentPosition_Node_2DList
!
!
!
 SUBROUTINE SetCurrentKey_Node_2DList( myList, key )
 ! S/R SetCurrentKey
 !
 !
 ! ========================================================================== !
 IMPLICIT NONE
  CLASS( NodeList_2D )   :: myList
  INTEGER, INTENT(in) :: key

     CALL myList % current % SetKey( key )

 END SUBROUTINE SetCurrentKey_Node_2DList
!
!
!
 SUBROUTINE GetCurrentKey_Node_2DList( myList, key )
 ! S/R GetCurrentKey
 !
 !
 ! ========================================================================== !
  IMPLICIT NONE
  CLASS( NodeList_2D )    :: myList
  INTEGER, INTENT(out) :: key

     CALL myList % current % GetKey( key )

 END SUBROUTINE GetCurrentKey_Node_2DList
!
!
!
 SUBROUTINE SetData_Node_2D( thisNode, x, y, nodeType )
 ! S/R SetData
 !  
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS( Node_2D ), INTENT(inout) :: thisNode
   REAL(prec), INTENT(in)       :: x, y
   INTEGER, INTENT(in)          :: nodeType

      CALL thisNode % SetPosition( x, y )
      CALL thisNode % SetType( nodeType )

 END SUBROUTINE SetData_Node_2D
!
!
!
 SUBROUTINE GetData_Node_2D( thisNode, x, y, nodeType )
 ! S/R GetData
 !  
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS( Node_2D ), INTENT(in) :: thisNode
   REAL(prec), INTENT(out)   :: x, y
   INTEGER, INTENT(out)      :: nodeType

      CALL thisNode % GetPosition( x, y )
      CALL thisNode % GetType( nodeType )

 END SUBROUTINE GetData_Node_2D
!
! --------------------------------------- Node_2D Key ----------------------------------------------- !
!
 SUBROUTINE GetKey_Node_2D( thisNode, key )
 ! S/R GetKey
 !  
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS( Node_2D ), INTENT(in) :: thisNode
   INTEGER, INTENT(out)      :: key
   
      key = thisNode % key
   
 END SUBROUTINE GetKey_Node_2D
!
!
!
 SUBROUTINE SetKey_Node_2D( thisNode, key )
 ! S/R SetKey
 !  
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS( Node_2D ), INTENT(inout) :: thisNode
   INTEGER, INTENT(in)          :: key
   
      thisNode % key = key
   
 END SUBROUTINE SetKey_Node_2D
!
!
! --------------------------------------- Node_2D Type ---------------------------------------------- !
!
 SUBROUTINE GetType_Node_2D( thisNode, nodeType )
 ! S/R GetType_Node_2D
 !  
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS( Node_2D ), INTENT(in) :: thisNode
   INTEGER, INTENT(out)      :: nodeType
   
      nodeType = thisNode % nodeType
   
 END SUBROUTINE GetType_Node_2D
!
!
!
 SUBROUTINE SetType_Node_2D( thisNode, nodeType )
 ! S/R SetType_Node_2D
 !  
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS( Node_2D ), INTENT(inout) :: thisNode
   INTEGER, INTENT(in)          :: nodeType
   
      thisNode % nodeType = nodeType
   
 END SUBROUTINE SetType_Node_2D
!
! --------------------------------------- Position ----------------------------------------------- !
!
 SUBROUTINE GetPosition_Node_2D( thisNode, x, y )
 ! S/R GetPosition_Node_2D
 !  
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS( Node_2D ), INTENT(in) :: thisNode
   REAL(prec), INTENT(out)   :: x, y
   
      x = thisNode % x
      y = thisNode % y
   
 END SUBROUTINE GetPosition_Node_2D
!
!
!
 SUBROUTINE SetPosition_Node_2D( thisNode, x, y )
 ! S/R SetPosition_Node_2D
 !  
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS( Node_2D ), INTENT(inout) :: thisNode
   REAL(prec), INTENT(in)       :: x, y
   
      thisNode % x = x
      thisNode % y = y
   
 END SUBROUTINE SetPosition_Node_2D
!
!
!==================================================================================================!
!----------------------------------- Linked List Routines  ----------------------------------------!
!==================================================================================================!
!
!
 FUNCTION ListIsEmpty_Node_2DList( myList ) RESULT( TorF )
 ! Function ListIsEmpty
 !
 !  This function checks if "myList" is an empty list. A logical is returned
 !  
 ! =============================================================================================== !
 ! DECLARATIONS
  IMPLICIT NONE
  CLASS( NodeList_2D ) :: myList
  LOGICAL           :: TorF

     TorF = .NOT.( ASSOCIATED( myList % head  ) )
     
 END FUNCTION ListIsEmpty_Node_2DList
!
!
!
 SUBROUTINE AddToList_Node_2DList( myList, x, y, nodeType, inKey )
 ! S/R AddToList
 !
 !  This subroutine adds an item to the linked list. If "inKey"
 !  is supplied, then this is filled in for the key data. Otherwise, the key is
 !  filled in as the previous records key plus one. The data is placed at the end
 !  of the linked list
 !  
 ! =============================================================================================== !
 ! DECLARATIONS
  IMPLICIT NONE
  CLASS( NodeList_2D ) :: myList
  REAL(prec)        :: x, y
  INTEGER           :: nodeType
  INTEGER, OPTIONAL :: inKey
  ! LOCAL
  TYPE( Node_2D ), POINTER :: previous
  INTEGER                  :: allocationStatus

     ! Check to see if this list is empty
     IF( myList % ListIsEmpty() )THEN
     
        ALLOCATE( myList % head, STAT = allocationStatus )
        IF( allocationStatus /=0 )THEN
           PRINT*, 'MODULE NodeList_2DClass.f90 : S/R AddToList : Memory not allocated for next entry in list.'
           ! An exception handler should be built to handle these problems
        ENDIF      
      
        ! Point the current position to the head
        myList % current => myList % head
        ! Set the data
        CALL myList % SetData( x, y, nodeType )
        
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
           PRINT*, 'MODULE NodeList_2DClass.f90 : S/R AddToList : Memory not allocated for next entry in list.'
           ! An exception handler should be built to handle these problems
        ENDIF      
        
        !  Temporarily hold onto the tail in case we need the key (if inKey is not present)  
        previous => myList % tail
        ! Reassign the tail
        myList % tail => myList % tail % next
        
        ! Set the current to the tail
        myList % current => myList % tail
  
        ! Fill in the data
        CALL myList % SetData( x, y, nodeType )
        
        ! Fill in the key information
        IF( PRESENT(inKey) )THEN
           CALL myList % SetKey( inKey )
        ELSE
           CALL myList % SetKey( previous % key + keyInc )
        ENDIF
        
        ! Point the next to null and the tail to current
        myList % current % next => NULL( )
        
     ENDIF

 END SUBROUTINE AddToList_Node_2DList
!
!
!
 SUBROUTINE RemoveCurrent_Node_2DList( myList )
 ! S/R RemoveCurrent
 !
 !  This subroutine removes the current item in the linked list and patches together the previous
 !  and next items.
 !  
 ! =============================================================================================== !
  IMPLICIT NONE
  CLASS( NodeList_2D ) :: myList
  ! LOCAL
  TYPE( Node_2D ), POINTER :: previous, pNext
  INTEGER               :: currentKey, thisKey

     CALL myList % GetKey( currentKey )
     
     ! Check to see if this list is empty
     IF( myList % ListIsEmpty() )THEN
     
        PRINT*, 'Module EdgeClass.f90 : S/R RemoveCurrent : List is empty. Nothing to remove'
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

 END SUBROUTINE RemoveCurrent_Node_2DList
!
!
!
 SUBROUTINE MoveToNext_Node_2DList( myList )
 ! S/R MoveToNext
 !
 ! =============================================================================================== !
 ! DECLARATIONS
  IMPLICIT NONE
  CLASS( NodeList_2D ) :: myList

     myList % current => myList % current % next

 END SUBROUTINE MoveToNext_Node_2DList
!
!
!
 SUBROUTINE MoveToHead_Node_2DList( myList )
 ! S/R MoveToHead
 !
 !  
 ! ========================================================================== !
  IMPLICIT NONE
  CLASS( NodeList_2D ) :: myList

     myList % current => myList % head

 END SUBROUTINE MoveToHead_Node_2DList
!
!
!
 SUBROUTINE MoveToTail_Node_2DList( myList )
 ! S/R MoveToTail
 !
 !  
 ! ========================================================================== !
  IMPLICIT NONE
  CLASS( NodeList_2D ) :: myList

     myList % current => myList % tail

 END SUBROUTINE MoveToTail_Node_2DList
!
!
!
  SUBROUTINE GetCount_Node_2DList( myList, numberOfNodes )
 ! S/R GetCount
 !
 !
 ! ===========================================================================================s==== !
  IMPLICIT NONE
  CLASS( NodeList_2D )    :: myList
  INTEGER, INTENT(out) :: numberOfNodes

     numberOfNodes = 0 ! Initialize the number of list items
     ! Check to see if this list is empty
     IF( myList % ListIsEmpty() )THEN
     
        PRINT*, 'Module NodeClass_1D.f90 : S/R GetCount : List is empty.'
        RETURN
        
     ELSE ! the list is not empty
    
        CALL myList % MoveToHead( ) ! Rewind the list
        
        DO WHILE( ASSOCIATED( myList % current ) )
        
           numberOfNodes = numberOfNodes + 1
           CALL myList % moveToNext( )

        ENDDO

     ENDIF

 END SUBROUTINE GetCount_Node_2DList
! SUBROUTINE PrintList( myList )
 ! S/R PrintList
 !
 ! =============================================================================================== !
 ! DECLARATIONS
!  IMPLICIT NONE
!  CLASS( NodeList_2D ) :: myList

!     myList % current => myList % head

!    PRINT*, '          Data        Key'
!    DO WHILE( ASSOCIATED( myList % current ) )
 
!        PRINT*, myList % current % listData, myList % current % key

!        CALL myList % MoveToNext()

!     ENDDO

! END SUBROUTINE PrintList
!
!
!==================================================================================================!
!--------------------------------- Type-Specific Routines  ----------------------------------------!
!==================================================================================================!
!
!
 SUBROUTINE ScaleNodePosition_Node_2D( myNode, xScale, yScale )
 ! S/R ScaleNode
 ! 
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS( Node_2D ), INTENT(inout) :: myNode
   REAL(prec), INTENT(in)       :: xScale, yScale

      myNode % x = xScale*(myNode % x)
      myNode % y = yScale*(myNode % y)

 END SUBROUTINE ScaleNodePosition_Node_2D
!
!
!
END MODULE NodeClass_2D
