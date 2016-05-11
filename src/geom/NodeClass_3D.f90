! NodeClass_3D.f90
! 
! Copyright 2015 Joseph Schoonover <schoonover.numerics@gmail.com>
! 
! NodeClass_3D.f90 is part of the Spectral Element Libraries in Fortran (SELF).
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
 
 
MODULE NodeClass_3D
! ========================================= Logs ================================================= !
!2016-05-11  Joseph Schoonover  schoonover.numerics@gmail.com 
!
! //////////////////////////////////////////////////////////////////////////////////////////////// ! 
!
! This module defines the mesh primitive "node". The node contains a node type information, the 
! position, and a node-to-element connectivity list. Linked-List type structure and routines are
! provided for those applications which require dynamic memory allocation in which the number of 
! nodes is unknown a'priori (e.g. Mesh Generation). One can use the Node_3D data-structure on its own
! without the linked-list structure in cases where the number of nodes is known.
!
! =======================================================================================

USE ModelPrecision
USE ModelFlags
USE ConstantsDictionary
USE LinkedListClass


IMPLICIT NONE

   TYPE Node_3D
      INTEGER, PRIVATE           :: key
      INTEGER, PRIVATE           :: nodeType      ! An INTEGER flag for INTERIOR or BOUNDARY
      REAL(prec), PRIVATE        :: x, y, z
    !  TYPE( LinkedList ), PUBLIC :: nodeToElement ! A linked list of element IDs which share this node.
      TYPE( Node_3D ), POINTER      :: next

      CONTAINS

      PROCEDURE :: Build => Build_Node_3D
      PROCEDURE :: Trash => Trash_Node_3D

      PROCEDURE :: GetData => GetData_Node_3D
      PROCEDURE :: SetData => SetData_Node_3D
      PROCEDURE :: GetKey => GetKey_Node_3D
      PROCEDURE :: SetKey => SetKey_Node_3D
      PROCEDURE :: GetType => GetType_Node_3D
      PROCEDURE :: SetType => SetType_Node_3D
      PROCEDURE :: GetPosition => GetPosition_Node_3D
      PROCEDURE :: SetPosition => SetPosition_Node_3D
      
      PROCEDURE :: ScaleNodePosition => ScaleNodePosition_Node_3D
   END TYPE Node_3D

   TYPE NodeList_3D      
      TYPE( Node_3D ), POINTER :: head, tail, current

      CONTAINS

      PROCEDURE :: Build => Build_Node_3DList
      PROCEDURE :: Trash => Trash_Node_3DList
      
      PROCEDURE :: GetData => GetCurrentData_Node_3DList
      PROCEDURE :: SetData => SetCurrentData_Node_3DList
      PROCEDURE :: GetPosition => GetCurrentPosition_Node_3DList
      PROCEDURE :: SetPosition => SetCurrentPosition_Node_3DList
      PROCEDURE :: GetKey => GetCurrentKey_Node_3DList
      PROCEDURE :: SetKey => SetCurrentKey_Node_3DList
      
      PROCEDURE :: ListIsEmpty => ListIsEmpty_Node_3DList
      PROCEDURE :: AddToList => AddToList_Node_3DList
      PROCEDURE :: RemoveCurrent => RemoveCurrent_Node_3DList
      PROCEDURE :: MoveToHead => MoveToHead_Node_3DList
      PROCEDURE :: MoveToNext => MoveToNext_Node_3DList
      PROCEDURE :: MoveToTail => MoveToTail_Node_3DList
      PROCEDURE :: GetCount => GetCount_Node_3DList

      

   END TYPE NodeList_3D

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
 SUBROUTINE Build_Node_3D( thisNode, x, y, z )
 ! S/R Build_Node_3D
 !  
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS( Node_3D ), INTENT(inout) :: thisNode
   REAL(prec), INTENT(in)       :: x, y, z

 
      CALL thisNode % SetPosition( x, y, z )
      
      CALL thisNode % SetType( INTERIOR )
      
     ! CALL thisNode % nodeToElement % Build( )
      
      thisNode % next => NULL(  )
      
 END SUBROUTINE Build_Node_3D
!
!
!
 SUBROUTINE Trash_Node_3D( thisNode )
 ! S/R Trash_Node_3D
 !  
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS( Node_3D ), INTENT(inout) :: thisNode

     ! CALL thisNode % nodeToElement % Trash( )

 END SUBROUTINE Trash_Node_3D
!
!
!
 SUBROUTINE Build_Node_3DList( myList )
 ! S/R BUILD_Node_3DList
 !
 !  This subroutine creates a linked-list by allocating space for the list
 !  head. The "next" attribute of the list head is nullified. 
 !  The tail is pointed to NULL and the current position in the list
 !  is pointed to the head. 
 ! 
 !   INPUT/OUTPUT : 
 !      CLASS( NodeList_3D) :: myList
 !
 ! =============================================================================================== !
 ! DECLARATIONS
  IMPLICIT NONE
  CLASS( NodeList_3D ) :: myList

    
     myList % head => NULL( )

     ! Point the tail to null
     myList % tail => NULL()

     ! Set the current position to Null
     myList % current => NULL( )
  
 END SUBROUTINE Build_Node_3DList
!
!
!  
 SUBROUTINE Trash_Node_3DList( myList )
 ! S/R Trash_Node_3DList
 ! 
 !  This subroutine traverses through the linked list and frees the memory
 !  that was allocated for the entries in the linked list.  
 !
 ! =============================================================================================== !
 ! DECLARATIONS
  IMPLICIT NONE
  CLASS( NodeList_3D ) :: myList
  ! LOCAL
  TYPE( Node_3D ), POINTER :: pNext

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
  
 END SUBROUTINE Trash_Node_3DList
!
!
!==================================================================================================!
!--------------------------------------- Accessors ------------------------------------------------!
!==================================================================================================!
!
!
 SUBROUTINE SetCurrentData_Node_3DList( myList, x, y, z, nodeType )
 ! S/R SetCurrentData
 !
 !
 ! ========================================================================== !
 IMPLICIT NONE
  CLASS( NodeList_3D )      :: myList
  REAL(prec), INTENT(in) :: x, y, z
  INTEGER, INTENT(in)    :: nodeType

     CALL myList % current % SetData( x,  y, z, nodeType )

 END SUBROUTINE SetCurrentData_Node_3DList
!
!
!
 SUBROUTINE GetCurrentData_Node_3DList( myList, x, y, z, nodeType )
 ! S/R GetCurrentData
 !
 !
 ! ========================================================================== !
 IMPLICIT NONE
  CLASS( NodeList_3D )       :: myList
  REAL(prec), INTENT(out) :: x, y, z
  INTEGER, INTENT(out)    :: nodeType

     CALL myList % current % GetData( x, y, z, nodeType )

 END SUBROUTINE GetCurrentData_Node_3DList
!
!
!
 SUBROUTINE SetCurrentPosition_Node_3DList( myList, x, y, z )
 ! S/R SetCurrentPosition
 !
 !
 ! ========================================================================== !
 IMPLICIT NONE
  CLASS( NodeList_3D )      :: myList
  REAL(prec), INTENT(in) :: x, y, z

     CALL myList % current % SetPosition( x, y, z )

 END SUBROUTINE SetCurrentPosition_Node_3DList
!
!
!
 SUBROUTINE GetCurrentPosition_Node_3DList( myList, x, y, z )
 ! S/R GetCurrentPosition
 !
 !
 ! ========================================================================== !
 IMPLICIT NONE
  CLASS( NodeList_3D )       :: myList
  REAL(prec), INTENT(out) :: x, y, z

     CALL myList % current % GetPosition( x, y, z )

 END SUBROUTINE GetCurrentPosition_Node_3DList
!
!
!
 SUBROUTINE SetCurrentKey_Node_3DList( myList, key )
 ! S/R SetCurrentKey
 !
 !
 ! ========================================================================== !
 IMPLICIT NONE
  CLASS( NodeList_3D )   :: myList
  INTEGER, INTENT(in) :: key

     CALL myList % current % SetKey( key )

 END SUBROUTINE SetCurrentKey_Node_3DList
!
!
!
 SUBROUTINE GetCurrentKey_Node_3DList( myList, key )
 ! S/R GetCurrentKey
 !
 !
 ! ========================================================================== !
  IMPLICIT NONE
  CLASS( NodeList_3D )    :: myList
  INTEGER, INTENT(out) :: key

     CALL myList % current % GetKey( key )

 END SUBROUTINE GetCurrentKey_Node_3DList
!
!
!
 SUBROUTINE SetData_Node_3D( thisNode, x, y, z, nodeType )
 ! S/R SetData
 !  
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS( Node_3D ), INTENT(inout) :: thisNode
   REAL(prec), INTENT(in)       :: x, y, z
   INTEGER, INTENT(in)          :: nodeType

      CALL thisNode % SetPosition( x, y, z )
      CALL thisNode % SetType( nodeType )

 END SUBROUTINE SetData_Node_3D
!
!
!
 SUBROUTINE GetData_Node_3D( thisNode, x, y, z, nodeType )
 ! S/R GetData
 !  
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS( Node_3D ), INTENT(in) :: thisNode
   REAL(prec), INTENT(out)   :: x, y, z
   INTEGER, INTENT(out)      :: nodeType

      CALL thisNode % GetPosition( x, y, z )
      CALL thisNode % GetType( nodeType )

 END SUBROUTINE GetData_Node_3D
!
! --------------------------------------- Node_3D Key ----------------------------------------------- !
!
 SUBROUTINE GetKey_Node_3D( thisNode, key )
 ! S/R GetKey
 !  
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS( Node_3D ), INTENT(in) :: thisNode
   INTEGER, INTENT(out)      :: key
   
      key = thisNode % key
   
 END SUBROUTINE GetKey_Node_3D
!
!
!
 SUBROUTINE SetKey_Node_3D( thisNode, key )
 ! S/R SetKey
 !  
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS( Node_3D ), INTENT(inout) :: thisNode
   INTEGER, INTENT(in)          :: key
   
      thisNode % key = key
   
 END SUBROUTINE SetKey_Node_3D
!
!
! --------------------------------------- Node_3D Type ---------------------------------------------- !
!
 SUBROUTINE GetType_Node_3D( thisNode, nodeType )
 ! S/R GetType_Node_3D
 !  
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS( Node_3D ), INTENT(in) :: thisNode
   INTEGER, INTENT(out)      :: nodeType
   
      nodeType = thisNode % nodeType
   
 END SUBROUTINE GetType_Node_3D
!
!
!
 SUBROUTINE SetType_Node_3D( thisNode, nodeType )
 ! S/R SetType_Node_3D
 !  
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS( Node_3D ), INTENT(inout) :: thisNode
   INTEGER, INTENT(in)          :: nodeType
   
      thisNode % nodeType = nodeType
   
 END SUBROUTINE SetType_Node_3D
!
! --------------------------------------- Position ----------------------------------------------- !
!
 SUBROUTINE GetPosition_Node_3D( thisNode, x, y, z )
 ! S/R GetPosition_Node_3D
 !  
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS( Node_3D ), INTENT(in) :: thisNode
   REAL(prec), INTENT(out)   :: x, y, z
   
      x = thisNode % x
      y = thisNode % y
      z = thisNode % z
   
 END SUBROUTINE GetPosition_Node_3D
!
!
!
 SUBROUTINE SetPosition_Node_3D( thisNode, x, y, z )
 ! S/R SetPosition_Node_3D
 !  
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS( Node_3D ), INTENT(inout) :: thisNode
   REAL(prec), INTENT(in)       :: x, y, z
   
      thisNode % x = x
      thisNode % y = y
      thisNode % z = z
      
 END SUBROUTINE SetPosition_Node_3D
!
!
!==================================================================================================!
!----------------------------------- Linked List Routines  ----------------------------------------!
!==================================================================================================!
!
!
 FUNCTION ListIsEmpty_Node_3DList( myList ) RESULT( TorF )
 ! Function ListIsEmpty
 !
 !  This function checks if "myList" is an empty list. A logical is returned
 !  
 ! =============================================================================================== !
 ! DECLARATIONS
  IMPLICIT NONE
  CLASS( NodeList_3D ) :: myList
  LOGICAL           :: TorF

     TorF = .NOT.( ASSOCIATED( myList % head  ) )
     
 END FUNCTION ListIsEmpty_Node_3DList
!
!
!
 SUBROUTINE AddToList_Node_3DList( myList, x, y, z, nodeType, inKey )
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
  CLASS( NodeList_3D ) :: myList
  REAL(prec)        :: x, y, z
  INTEGER           :: nodeType
  INTEGER, OPTIONAL :: inKey
  ! LOCAL
  TYPE( Node_3D ), POINTER :: previous
  INTEGER               :: allocationStatus

     ! Check to see if this list is empty
     IF( myList % ListIsEmpty() )THEN
     
        ALLOCATE( myList % head, STAT = allocationStatus )
        IF( allocationStatus /=0 )THEN
           PRINT*, 'MODULE NodeList_3DClass.f90 : S/R AddToList : Memory not allocated for next entry in list.'
           ! An exception handler should be built to handle these problems
        ENDIF      
      
        ! Point the current position to the head
        myList % current => myList % head
        ! Set the data
        CALL myList % SetData( x, y, z, nodeType )
        
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
           PRINT*, 'MODULE NodeList_3DClass.f90 : S/R AddToList : Memory not allocated for next entry in list.'
           ! An exception handler should be built to handle these problems
        ENDIF      
        
        !  Temporarily hold onto the tail in case we need the key (if inKey is not present)  
        previous => myList % tail
        ! Reassign the tail
        myList % tail => myList % tail % next
        
        ! Set the current to the tail
        myList % current => myList % tail
  
        ! Fill in the data
        CALL myList % SetData( x, y, z, nodeType )
        
        ! Fill in the key information
        IF( PRESENT(inKey) )THEN
           CALL myList % SetKey( inKey )
        ELSE
           CALL myList % SetKey( previous % key + keyInc )
        ENDIF
        
        ! Point the next to null and the tail to current
        myList % current % next => NULL( )
        
     ENDIF

 END SUBROUTINE AddToList_Node_3DList
!
!
!
 SUBROUTINE RemoveCurrent_Node_3DList( myList )
 ! S/R RemoveCurrent
 !
 !  This subroutine removes the current item in the linked list and patches together the previous
 !  and next items.
 !  
 ! =============================================================================================== !
  IMPLICIT NONE
  CLASS( NodeList_3D ) :: myList
  ! LOCAL
  TYPE( Node_3D ), POINTER :: previous, pNext
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

 END SUBROUTINE RemoveCurrent_Node_3DList
!
!
!
 SUBROUTINE MoveToNext_Node_3DList( myList )
 ! S/R MoveToNext
 !
 ! =============================================================================================== !
 ! DECLARATIONS
  IMPLICIT NONE
  CLASS( NodeList_3D ) :: myList

     myList % current => myList % current % next

 END SUBROUTINE MoveToNext_Node_3DList
!
!
!
 SUBROUTINE MoveToHead_Node_3DList( myList )
 ! S/R MoveToHead
 !
 !  
 ! ========================================================================== !
  IMPLICIT NONE
  CLASS( NodeList_3D ) :: myList

     myList % current => myList % head

 END SUBROUTINE MoveToHead_Node_3DList
!
!
!
 SUBROUTINE MoveToTail_Node_3DList( myList )
 ! S/R MoveToTail
 !
 !  
 ! ========================================================================== !
  IMPLICIT NONE
  CLASS( NodeList_3D ) :: myList

     myList % current => myList % tail

 END SUBROUTINE MoveToTail_Node_3DList
!
!
!
  SUBROUTINE GetCount_Node_3DList( myList, numberOfNodes )
 ! S/R GetCount
 !
 !
 ! ===========================================================================================s==== !
  IMPLICIT NONE
  CLASS( NodeList_3D )    :: myList
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

 END SUBROUTINE GetCount_Node_3DList
! SUBROUTINE PrintList( myList )
 ! S/R PrintList
 !
 ! =============================================================================================== !
 ! DECLARATIONS
!  IMPLICIT NONE
!  CLASS( NodeList_3D ) :: myList

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
 SUBROUTINE ScaleNodePosition_Node_3D( myNode, xScale, yScale, zScale )
 ! S/R ScaleNode
 ! 
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS( Node_3D ), INTENT(inout) :: myNode
   REAL(prec), INTENT(in)       :: xScale, yScale, zScale

      myNode % x = xScale*(myNode % x)
      myNode % y = yScale*(myNode % y)
      myNode % z = zScale*(myNode % z)

 END SUBROUTINE ScaleNodePosition_Node_3D
!
!
!
END MODULE NodeClass_3D
