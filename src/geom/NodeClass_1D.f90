MODULE NodeClass_1D
! NodeClass_1D.f90
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
! nodes is unknown a'priori (e.g. Mesh Generation). One can use the Node data-structure on its own
! without the linked-list structure in cases where the number of nodes is known.
!
! =======================================================================================

USE ModelPrecision
USE ModelFlags
USE ConstantsDictionary
USE LinkedListClass


IMPLICIT NONE

   TYPE Node
      INTEGER, PRIVATE           :: key
      INTEGER, PRIVATE           :: nodeType      ! An INTEGER flag for INTERIOR or BOUNDARY
      REAL(prec), PRIVATE        :: x
      TYPE( LinkedList ), PUBLIC :: nodeToElement ! A linked list of element IDs which share this node.
      TYPE( Node ), POINTER      :: next

      CONTAINS

      PROCEDURE :: Build => Build_Node
      PROCEDURE :: Trash => Trash_Node

      PROCEDURE :: GetData => GetData_Node
      PROCEDURE :: SetData => SetData_Node
      PROCEDURE :: GetKey => GetKey_Node
      PROCEDURE :: SetKey => SetKey_Node
      PROCEDURE :: GetType => GetType_Node
      PROCEDURE :: SetType => SetType_Node
      PROCEDURE :: GetPosition => GetPosition_Node
      PROCEDURE :: SetPosition => SetPosition_Node
      
      
   END TYPE Node

   TYPE NodeList      
      TYPE( Node ), POINTER :: head, tail, current

      CONTAINS

      PROCEDURE :: Build => Build_NodeList
      PROCEDURE :: Trash => Trash_NodeList
      
      PROCEDURE :: GetData => GetCurrentData_NodeList
      PROCEDURE :: SetData => SetCurrentData_NodeList
      PROCEDURE :: GetPosition => GetCurrentPosition_NodeList
      PROCEDURE :: SetPosition => SetCurrentPosition_NodeList
      PROCEDURE :: GetKey => GetCurrentKey_NodeList
      PROCEDURE :: SetKey => SetCurrentKey_NodeList
      
      PROCEDURE :: ListIsEmpty => ListIsEmpty_NodeList
      PROCEDURE :: AddToList => AddToList_NodeList
      PROCEDURE :: RemoveCurrent => RemoveCurrent_NodeList
      PROCEDURE :: MoveToHead => MoveToHead_NodeList
      PROCEDURE :: MoveToNext => MoveToNext_NodeList
      PROCEDURE :: MoveToTail => MoveToTail_NodeList
      PROCEDURE :: GetCount => GetCount_NodeList

   END TYPE NodeList

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
 SUBROUTINE Build_Node( thisNode, xPosition )
 ! S/R Build_Node
 !  
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS( Node ), INTENT(inout) :: thisNode
   REAL(prec), INTENT(in)       :: xPosition

 
      CALL thisNode % SetPosition( xPosition )
      
      CALL thisNode % SetType( INTERIOR )
      
      CALL thisNode % nodeToElement % Build( )
      
      thisNode % next => NULL(  )
      
 END SUBROUTINE Build_Node
!
!
!
 SUBROUTINE Trash_Node( thisNode )
 ! S/R Trash_Node
 !  
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS( Node ), INTENT(inout) :: thisNode

      CALL thisNode % nodeToElement % Trash( )

 END SUBROUTINE Trash_Node
!
!
!
 SUBROUTINE Build_NodeList( myList )
 ! S/R BUILD_NodeList
 !
 !  This subroutine creates a linked-list by allocating space for the list
 !  head. The "next" attribute of the list head is nullified. 
 !  The tail is pointed to NULL and the current position in the list
 !  is pointed to the head. 
 ! 
 !   INPUT/OUTPUT : 
 !      CLASS( NodeList) :: myList
 !
 ! =============================================================================================== !
 ! DECLARATIONS
  IMPLICIT NONE
  CLASS( NodeList ) :: myList

    
     myList % head => NULL( )

     ! Point the tail to null
     myList % tail => NULL()

     ! Set the current position to Null
     myList % current => NULL( )
  
 END SUBROUTINE Build_NodeList
!
!
!  
 SUBROUTINE Trash_NodeList( myList )
 ! S/R Trash_NodeList
 ! 
 !  This subroutine traverses through the linked list and frees the memory
 !  that was allocated for the entries in the linked list.  
 !
 ! =============================================================================================== !
 ! DECLARATIONS
  IMPLICIT NONE
  CLASS( NodeList ) :: myList
  ! LOCAL
  TYPE( Node ), POINTER :: pNext

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
  
 END SUBROUTINE Trash_NodeList
!
!
!==================================================================================================!
!--------------------------------------- Accessors ------------------------------------------------!
!==================================================================================================!
!
!
 SUBROUTINE SetCurrentData_NodeList( myList, x, nodeType )
 ! S/R SetCurrentData
 !
 !
 ! ========================================================================== !
 IMPLICIT NONE
  CLASS( NodeList )      :: myList
  REAL(prec), INTENT(in) :: x
  INTEGER, INTENT(in)    :: nodeType

     CALL myList % current % SetData( x, nodeType )

 END SUBROUTINE SetCurrentData_NodeList
!
!
!
 SUBROUTINE GetCurrentData_NodeList( myList, x, nodeType )
 ! S/R GetCurrentData
 !
 !
 ! ========================================================================== !
 IMPLICIT NONE
  CLASS( NodeList )       :: myList
  REAL(prec), INTENT(out) :: x
  INTEGER, INTENT(out)    :: nodeType

     CALL myList % current % GetData( x, nodeType )

 END SUBROUTINE GetCurrentData_NodeList
!
!
!
 SUBROUTINE SetCurrentPosition_NodeList( myList, x )
 ! S/R SetCurrentPosition
 !
 !
 ! ========================================================================== !
 IMPLICIT NONE
  CLASS( NodeList )      :: myList
  REAL(prec), INTENT(in) :: x

     CALL myList % current % SetPosition( x )

 END SUBROUTINE SetCurrentPosition_NodeList
!
!
!
 SUBROUTINE GetCurrentPosition_NodeList( myList, x )
 ! S/R GetCurrentPosition
 !
 !
 ! ========================================================================== !
 IMPLICIT NONE
  CLASS( NodeList )       :: myList
  REAL(prec), INTENT(out) :: x

     CALL myList % current % GetPosition( x )

 END SUBROUTINE GetCurrentPosition_NodeList
!
!
!
 SUBROUTINE SetCurrentKey_NodeList( myList, key )
 ! S/R SetCurrentKey
 !
 !
 ! ========================================================================== !
 IMPLICIT NONE
  CLASS( NodeList )   :: myList
  INTEGER, INTENT(in) :: key

     CALL myList % current % SetKey( key )

 END SUBROUTINE SetCurrentKey_NodeList
!
!
!
 SUBROUTINE GetCurrentKey_NodeList( myList, key )
 ! S/R GetCurrentKey
 !
 !
 ! ========================================================================== !
  IMPLICIT NONE
  CLASS( NodeList )    :: myList
  INTEGER, INTENT(out) :: key

     CALL myList % current % GetKey( key )

 END SUBROUTINE GetCurrentKey_NodeList
!
!
!
 SUBROUTINE SetData_Node( thisNode, x, nodeType )
 ! S/R SetData
 !  
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS( Node ), INTENT(inout) :: thisNode
   REAL(prec), INTENT(in)       :: x
   INTEGER, INTENT(in)          :: nodeType

      CALL thisNode % SetPosition( x )
      CALL thisNode % SetType( nodeType )

 END SUBROUTINE SetData_Node
!
!
!
 SUBROUTINE GetData_Node( thisNode, x, nodeType )
 ! S/R GetData
 !  
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS( Node ), INTENT(in) :: thisNode
   REAL(prec), INTENT(out)   :: x
   INTEGER, INTENT(out)      :: nodeType

      CALL thisNode % GetPosition( x )
      CALL thisNode % GetType( nodeType )

 END SUBROUTINE GetData_Node
!
! --------------------------------------- Node Key ----------------------------------------------- !
!
 SUBROUTINE GetKey_Node( thisNode, key )
 ! S/R GetKey
 !  
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS( Node ), INTENT(in) :: thisNode
   INTEGER, INTENT(out)      :: key
   
      key = thisNode % key
   
 END SUBROUTINE GetKey_Node
!
!
!
 SUBROUTINE SetKey_Node( thisNode, key )
 ! S/R SetKey
 !  
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS( Node ), INTENT(inout) :: thisNode
   INTEGER, INTENT(in)          :: key
   
      thisNode % key = key
   
 END SUBROUTINE SetKey_Node
!
!
! --------------------------------------- Node Type ---------------------------------------------- !
!
 SUBROUTINE GetType_Node( thisNode, nodeType )
 ! S/R GetType_Node
 !  
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS( Node ), INTENT(in) :: thisNode
   INTEGER, INTENT(out)      :: nodeType
   
      nodeType = thisNode % nodeType
   
 END SUBROUTINE GetType_Node
!
!
!
 SUBROUTINE SetType_Node( thisNode, nodeType )
 ! S/R SetType_Node
 !  
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS( Node ), INTENT(inout) :: thisNode
   INTEGER, INTENT(in)          :: nodeType
   
      thisNode % nodeType = nodeType
   
 END SUBROUTINE SetType_Node
!
! --------------------------------------- Position ----------------------------------------------- !
!
 SUBROUTINE GetPosition_Node( thisNode, x )
 ! S/R GetPosition_Node
 !  
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS( Node ), INTENT(in) :: thisNode
   REAL(prec), INTENT(out)   :: x
   
      x = thisNode % x
   
 END SUBROUTINE GetPosition_Node
!
!
!
 SUBROUTINE SetPosition_Node( thisNode, x )
 ! S/R SetPosition_Node
 !  
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS( Node ), INTENT(inout) :: thisNode
   REAL(prec), INTENT(in)       :: x
   
      thisNode % x = x
   
 END SUBROUTINE SetPosition_Node
!
!
!==================================================================================================!
!----------------------------------- Linked List Routines  ----------------------------------------!
!==================================================================================================!
!
!
 FUNCTION ListIsEmpty_NodeList( myList ) RESULT( TorF )
 ! Function ListIsEmpty
 !
 !  This function checks if "myList" is an empty list. A logical is returned
 !  
 ! =============================================================================================== !
 ! DECLARATIONS
  IMPLICIT NONE
  CLASS( NodeList ) :: myList
  LOGICAL           :: TorF

     TorF = .NOT.( ASSOCIATED( myList % head  ) )
     
 END FUNCTION ListIsEmpty_NodeList
!
!
!
 SUBROUTINE AddToList_NodeList( myList, x, nodeType, inKey )
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
  CLASS( NodeList ) :: myList
  REAL(prec)        :: x
  INTEGER           :: nodeType
  INTEGER, OPTIONAL :: inKey
  ! LOCAL
  LOGICAL               :: isAssigned
  TYPE( Node ), POINTER :: previous
  INTEGER               :: allocationStatus

     ! Check to see if this list is empty
     IF( myList % ListIsEmpty() )THEN
     
        ALLOCATE( myList % head, STAT = allocationStatus )
        IF( allocationStatus /=0 )THEN
           PRINT*, 'MODULE NodeListClass.f90 : S/R AddToList : Memory not allocated for next entry in list.'
           ! An exception handler should be built to handle these problems
        ENDIF      
      
        ! Point the current position to the head
        myList % current => myList % head
        ! Set the data
        CALL myList % SetData( x, nodeType )
        
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
           PRINT*, 'MODULE NodeListClass.f90 : S/R AddToList : Memory not allocated for next entry in list.'
           ! An exception handler should be built to handle these problems
        ENDIF      
        
        !  Temporarily hold onto the tail in case we need the key (if inKey is not present)  
        previous => myList % tail
        ! Reassign the tail
        myList % tail => myList % tail % next
        
        ! Set the current to the tail
        myList % current => myList % tail
  
        ! Fill in the data
        CALL myList % SetData( x, nodeType )
        
        ! Fill in the key information
        IF( PRESENT(inKey) )THEN
           CALL myList % SetKey( inKey )
        ELSE
           CALL myList % SetKey( previous % key + keyInc )
        ENDIF
        
        ! Point the next to null and the tail to current
        myList % current % next => NULL( )
        
     ENDIF

 END SUBROUTINE AddToList_NodeList
!
!
!
 SUBROUTINE RemoveCurrent_NodeList( myList )
 ! S/R RemoveCurrent
 !
 !  This subroutine removes the current item in the linked list and patches together the previous
 !  and next items.
 !  
 ! =============================================================================================== !
  IMPLICIT NONE
  CLASS( NodeList ) :: myList
  ! LOCAL
  TYPE( Node ), POINTER :: previous, pNext
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

 END SUBROUTINE RemoveCurrent_NodeList
!
!
!
 SUBROUTINE MoveToNext_NodeList( myList )
 ! S/R MoveToNext
 !
 ! =============================================================================================== !
 ! DECLARATIONS
  IMPLICIT NONE
  CLASS( NodeList ) :: myList

     myList % current => myList % current % next

 END SUBROUTINE MoveToNext_NodeList
!
!
!
 SUBROUTINE MoveToHead_NodeList( myList )
 ! S/R MoveToHead
 !
 !  
 ! ========================================================================== !
  IMPLICIT NONE
  CLASS( NodeList ) :: myList

     myList % current => myList % head

 END SUBROUTINE MoveToHead_NodeList
!
!
!
 SUBROUTINE MoveToTail_NodeList( myList )
 ! S/R MoveToTail
 !
 !  
 ! ========================================================================== !
  IMPLICIT NONE
  CLASS( NodeList ) :: myList

     myList % current => myList % tail

 END SUBROUTINE MoveToTail_NodeList
!
!
!
  SUBROUTINE GetCount_NodeList( myList, numberOfNodes )
 ! S/R GetCount
 !
 !
 ! ===========================================================================================s==== !
  IMPLICIT NONE
  CLASS( NodeList )    :: myList
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

 END SUBROUTINE GetCount_NodeList
! SUBROUTINE PrintList( myList )
 ! S/R PrintList
 !
 ! =============================================================================================== !
 ! DECLARATIONS
!  IMPLICIT NONE
!  CLASS( NodeList ) :: myList

!     myList % current => myList % head

!    PRINT*, '          Data        Key'
!    DO WHILE( ASSOCIATED( myList % current ) )
 
!        PRINT*, myList % current % listData, myList % current % key

!        CALL myList % MoveToNext()

!     ENDDO

! END SUBROUTINE PrintList
END MODULE NodeClass_1D
