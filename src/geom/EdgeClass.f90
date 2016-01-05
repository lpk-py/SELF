MODULE EdgeClass
! EdgeClass.f90 (v2.1 - 12 Dec. 2015)
! 
! schoonover.numerics@gmail.com
! 
! o (ver 1.0) April 2014
! o (ver 2.1) Dec 2015
!
!
! The edge class stores identification information for it's starting and terminating nodes and 
! neighboring elements. In this sense, edges are viewed as the book-keeping side of the mesh 
! description. Because it is anticipated that unstructured meshes may be used in which the number
! of edges is unknown a'priori, a linked list structure is used. Utilities are provided which also
! allow for the removal of individual edges from anywhere in the list, which is useful for 
! Adaptive Mesh Refinement (AMR) and mesh generation algorithms.
!
! 
!  
! =======================================================================================
USE ModelPrecision
USE ConstantsDictionary
USE ModelFlags


IMPLICIT NONE


   TYPE Edge
   
      INTEGER, PRIVATE      :: key               ! The edge key
      INTEGER, PRIVATE      :: nodeIDs(1:2)      ! Node IDs which start and terminate this edge
      INTEGER, PRIVATE      :: elementIDs(1:2)   ! Neighboring elements IDs across the edge
      INTEGER, PRIVATE      :: elementSides(1:2) ! Local side IDs for the neighboring elements
      INTEGER, PRIVATE      :: start, inc        ! Loop start and increment for the secondary element side
      TYPE( Edge ), POINTER :: next => NULL( )


      CONTAINS
 
      PROCEDURE :: Build => Build_Edge

      PROCEDURE :: SetKey => SetKey_Edge
      PROCEDURE :: GetKey => GetKey_Edge
      PROCEDURE :: SetData => SetData_Edge
      PROCEDURE :: GetData => GetData_Edge
      PROCEDURE :: SetNodeIDs => SetNodeIDs_Edge
      PROCEDURE :: GetNodeIDs => GetNodeIDs_Edge
      PROCEDURE :: SetElementIDs => SetElementIDs_Edge
      PROCEDURE :: GetElementIDs => GetElementIDs_Edge
      PROCEDURE :: SetPrimaryElementID => SetPrimaryElementID_Edge
      PROCEDURE :: GetPrimaryElementID => GetPrimaryElementID_Edge
      PROCEDURE :: SetSecondaryElementID => SetSecondaryElementID_Edge
      PROCEDURE :: GetSecondaryElementID => GetSecondaryElementID_Edge
      PROCEDURE :: SetElementSides => SetElementSides_Edge
      PROCEDURE :: GetElementSides => GetElementSides_Edge
      PROCEDURE :: SetPrimaryElementSide => SetPrimaryElementSide_Edge
      PROCEDURE :: GetPrimaryElementSide => GetPrimaryElementSide_Edge
      PROCEDURE :: SetSecondaryElementSide => SetSecondaryElementSide_Edge
      PROCEDURE :: GetSecondaryElementSide => GetSecondaryElementSide_Edge
      PROCEDURE :: SetStart => SetStart_Edge
      PROCEDURE :: GetStart => GetStart_Edge
      PROCEDURE :: SetIncrement => SetIncrement_Edge
      PROCEDURE :: GetIncrement => GetIncrement_Edge

   END TYPE Edge

   TYPE EdgeList      
      TYPE( Edge ), POINTER :: head, tail, current

      CONTAINS

      PROCEDURE :: Build => Build_EdgeList
      PROCEDURE :: Trash => Trash_EdgeList
      
      
      PROCEDURE :: GetData => GetCurrentData_EdgeList
      PROCEDURE :: SetData => SetCurrentData_EdgeList
      PROCEDURE :: GetKey => GetCurrentKey_EdgeList
      PROCEDURE :: SetKey => SetCurrentKey_EdgeList
      
      PROCEDURE :: ListIsEmpty => ListIsEmpty_EdgeList
      PROCEDURE :: AddToList => AddToList_EdgeList
      PROCEDURE :: RemoveCurrent => RemoveCurrent_EdgeList
      PROCEDURE :: MoveToHead => MoveToHead_EdgeList
      PROCEDURE :: MoveToNext => MoveToNext_EdgeList
      PROCEDURE :: MoveToTail => MoveToTail_EdgeList
      PROCEDURE :: GetCount => GetCount_EdgeList
      
   !   PROCEDURE :: PrintList

   END TYPE EdgeList
   
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
SUBROUTINE Build_Edge( myEdge )
 ! S/R BUILD_Edge
 !
 !
 ! =============================================================================================== !
 ! DECLARATIONS
  IMPLICIT NONE
  CLASS( Edge ) :: myEdge

     myEdge % nodeIDs      = NO_NORMAL_FLOW
     myEdge % elementIDs   = NO_NORMAL_FLOW
     myEdge % elementSides = 0
     myEdge % start        = 1 
     myEdge % inc          = 1
     myEdge % next => NULL()
     
 END SUBROUTINE Build_Edge 
!
!
!
 SUBROUTINE Build_EdgeList( myList )
 ! S/R BUILD_EdgeList
 !
 !  This subroutine creates a linked-list by allocating space for the list
 !  head. The "next" attribute of the list head is nullified. 
 !  The tail is pointed to NULL and the current position in the list
 !  is pointed to the head. 
 ! 
 !   INPUT/OUTPUT : 
 !      CLASS( EdgeList) :: myList
 !
 ! =============================================================================================== !
 ! DECLARATIONS
  IMPLICIT NONE
  CLASS( EdgeList ) :: myList

    
     myList % head => NULL( )

     ! Point the tail to null
     myList % tail => NULL()

     ! Set the current position to Null
     myList % current => NULL( )
  
 END SUBROUTINE Build_EdgeList
!
!
!  
 SUBROUTINE Trash_EdgeList( myList )
 ! S/R Trash_EdgeList
 ! 
 !  This subroutine traverses through the linked list and frees the memory
 !  that was allocated for the entries in the linked list.  
 !
 ! =============================================================================================== !
 ! DECLARATIONS
  IMPLICIT NONE
  CLASS( EdgeList ) :: myList
  ! LOCAL
  TYPE( Edge ), POINTER :: pNext

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
  
 END SUBROUTINE Trash_EdgeList
!
!
!==================================================================================================!
!--------------------------------------- Accessors ------------------------------------------------!
!==================================================================================================!
!
!
 SUBROUTINE SetCurrentData_EdgeList( myList, nodeIDs, elementIDs, elementSides, key, start, inc )
 ! S/R SetCurrentData
 !
 !
 ! ========================================================================== !
 IMPLICIT NONE
  CLASS( EdgeList )   :: myList
  INTEGER, INTENT(in) :: nodeIDs(1:2), elementIDs(1:2), elementSides(1:2)
  INTEGER, INTENT(in) :: key, start, inc

     CALL myList % current % SetData( nodeIDs, elementIDs, elementSides, key, start, inc )

 END SUBROUTINE SetCurrentData_EdgeList
!
!
!
 SUBROUTINE GetCurrentData_EdgeList( myList, nodeIDs, elementIDs, elementSides, key, start, inc )
 ! S/R GetCurrentData
 !
 !
 ! ========================================================================== !
 IMPLICIT NONE
  CLASS( EdgeList )    :: myList
  INTEGER, INTENT(out) :: nodeIDs(1:2), elementIDs(1:2), elementSides(1:2)
  INTEGER, INTENT(out) :: key, start, inc

     CALL myList % current % GetData( nodeIDs, elementIDs, elementSides, key, start, inc )

 END SUBROUTINE GetCurrentData_EdgeList
!
!
!
 SUBROUTINE SetCurrentKey_EdgeList( myList, key )
 ! S/R SetCurrentKey
 !
 !
 ! ========================================================================== !
 IMPLICIT NONE
  CLASS( EdgeList )   :: myList
  INTEGER, INTENT(in) :: key

     CALL myList % current % SetKey( key )

 END SUBROUTINE SetCurrentKey_EdgeList
!
!
!
 SUBROUTINE GetCurrentKey_EdgeList( myList, key )
 ! S/R GetCurrentKey
 !
 !
 ! ========================================================================== !
  IMPLICIT NONE
  CLASS( EdgeList )    :: myList
  INTEGER, INTENT(out) :: key

     CALL myList % current % GetKey( key )

 END SUBROUTINE GetCurrentKey_EdgeList
!
!
!
 SUBROUTINE SetKey_Edge( thisEdge, key )
 ! S/R SetKey
 !  
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS( Edge ), INTENT(inout) :: thisEdge
   INTEGER, INTENT(in)          :: key

      thisEdge % key = key      

 END SUBROUTINE SetKey_Edge
!
!
!
 SUBROUTINE GetKey_Edge( thisEdge, key )
 ! S/R GetKey
 !  
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS( Edge ), INTENT(in) :: thisEdge
   INTEGER, INTENT(out)      :: key

      key = thisEdge % key

 END SUBROUTINE GetKey_Edge
!
!
!
 SUBROUTINE SetData_Edge( thisEdge, nodeIDs, elementIDs, elementSides, key, start, inc )
 ! S/R SetData
 !  
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS( Edge ), INTENT(inout) :: thisEdge
   INTEGER, INTENT(in)          :: nodeIDs(1:2), elementIDs(1:2), elementSides(1:2)
   INTEGER, INTENT(in)          :: key, start, inc

      CALL thisEdge % SetNodeIDs( nodeIDs )
      CALL thisEdge % SetElementIDs( elementIDs )
      CALL thisEdge % SetElementSides( elementSides )
      CALL thisEdge % SetKey( key )
      CALL thisEdge % SetStart( start )
      CALL thisEdge % SetIncrement( inc )

 END SUBROUTINE SetData_Edge
!
!
!
 SUBROUTINE GetData_Edge( thisEdge, nodeIDs, elementIDs, elementSides, key, start, inc )
 ! S/R GetData
 !  
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS( Edge ), INTENT(in) :: thisEdge
   INTEGER, INTENT(out)      :: nodeIDs(1:2), elementIDs(1:2), elementSides(1:2)
   INTEGER, INTENT(out)      :: key, start, inc

      CALL thisEdge % GetNodeIDs( nodeIDs )
      CALL thisEdge % GetElementIDs( elementIDs )
      CALL thisEdge % GetElementSides( elementSides )
      CALL thisEdge % GetKey( key )
      CALL thisEdge % GetStart( start )
      CALL thisEdge % GetIncrement( inc )

 END SUBROUTINE GetData_Edge
!
!
!
 SUBROUTINE SetNodeIDs_Edge( thisEdge, nodeIDs )
 ! S/R SetNodeIDs
 !  
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS( Edge ), INTENT(inout) :: thisEdge
   INTEGER, INTENT(in)          :: nodeIDs(1:2)

      thisEdge % nodeIDs = nodeIDs

 END SUBROUTINE SetNodeIDs_Edge
!
!
!
 SUBROUTINE GetNodeIDs_Edge( thisEdge, nodeIDs )
 ! S/R GetNodeIDs
 !  
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS( Edge ), INTENT(in) :: thisEdge
   INTEGER, INTENT(out)      :: nodeIDs(1:2)

      nodeIDs = thisEdge % nodeIDs

 END SUBROUTINE GetNodeIDs_Edge
!
!
!
 SUBROUTINE SetElementIDs_Edge( thisEdge, elementIDs )
 ! S/R SetElementIDs
 !  
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS( Edge ), INTENT(inout) :: thisEdge
   INTEGER, INTENT(in)          :: elementIDs(1:2)

      thisEdge % elementIDs = elementIDs

 END SUBROUTINE SetElementIDs_Edge
!
!
!
 SUBROUTINE GetElementIDs_Edge( thisEdge, elementIDs )
 ! S/R GetElementIDs
 !  
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS( Edge ), INTENT(in) :: thisEdge
   INTEGER, INTENT(out)      :: elementIDs(1:2)

      elementIDs = thisEdge % elementIDs

 END SUBROUTINE GetElementIDs_Edge
!
!
!
 SUBROUTINE SetPrimaryElementID_Edge( thisEdge, elementID )
 ! S/R SetPrimaryElementID
 !  
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS( Edge ), INTENT(inout) :: thisEdge
   INTEGER, INTENT(in)          :: elementID

      thisEdge % elementIDs(1) = elementID
      
 END SUBROUTINE SetPrimaryElementID_Edge
!
!
!
 SUBROUTINE GetPrimaryElementID_Edge( thisEdge, elementID )
 ! S/R GetPrimaryElementID
 !  
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS( Edge ), INTENT(in) :: thisEdge
   INTEGER, INTENT(out)      :: elementID

      elementID = thisEdge % elementIDs(1)

 END SUBROUTINE GetPrimaryElementID_Edge
!
!
!
 SUBROUTINE SetSecondaryElementID_Edge( thisEdge, elementID )
 ! S/R SetSecondaryElementID
 !  
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS( Edge ), INTENT(inout) :: thisEdge
   INTEGER, INTENT(in)          :: elementID

      thisEdge % elementIDs(2) = elementID
      
 END SUBROUTINE SetSecondaryElementID_Edge
!
!
!
 SUBROUTINE GetSecondaryElementID_Edge( thisEdge, elementID )
 ! S/R GetSecondaryElementID
 !  
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS( Edge ), INTENT(in) :: thisEdge
   INTEGER, INTENT(out)      :: elementID

      elementID = thisEdge % elementIDs(2)

 END SUBROUTINE GetSecondaryElementID_Edge
!
!
!
 SUBROUTINE SetElementSides_Edge( thisEdge, elementSides )
 ! S/R SetElementSides
 !  
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS( Edge ), INTENT(inout) :: thisEdge
   INTEGER, INTENT(in)          :: elementSides(1:2)

      thisEdge % elementSides = elementSides
       
 END SUBROUTINE SetElementSides_Edge
!
!
!
 SUBROUTINE GetElementSides_Edge( thisEdge, elementSides )
 ! S/R GetElementSides
 !  
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS( Edge ), INTENT(in) :: thisEdge
   INTEGER, INTENT(out)      :: elementSides(1:2)

      elementSides = thisEdge % elementSides

 END SUBROUTINE GetElementSides_Edge
!
!
!
 SUBROUTINE SetPrimaryElementSide_Edge( thisEdge, elementSide )
 ! S/R SetPrimaryElementSide
 !  
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS( Edge ), INTENT(inout) :: thisEdge
   INTEGER, INTENT(in)          :: elementSide

      thisEdge % elementSides(1) = elementSide
      
 END SUBROUTINE SetPrimaryElementSide_Edge
!
!
!
 SUBROUTINE GetPrimaryElementSide_Edge( thisEdge, elementSide )
 ! S/R GetPrimaryElementSide
 !  
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS( Edge ), INTENT(in) :: thisEdge
   INTEGER, INTENT(out)      :: elementSide

      elementSide = thisEdge % elementSides(1)

 END SUBROUTINE GetPrimaryElementSide_Edge
!
!
!
 SUBROUTINE SetSecondaryElementSide_Edge( thisEdge, elementSide )
 ! S/R SetSecondaryElementSide
 !  
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS( Edge ), INTENT(inout) :: thisEdge
   INTEGER, INTENT(in)          :: elementSide

      thisEdge % elementSides(2) = elementSide
      
 END SUBROUTINE SetSecondaryElementSide_Edge
!
!
!
 SUBROUTINE GetSecondaryElementSide_Edge( thisEdge, elementSide )
 ! S/R GetSecondaryElementSide
 !  
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS( Edge ), INTENT(in) :: thisEdge
   INTEGER, INTENT(out)      :: elementSide

      elementSide = thisEdge % elementSides(2)

 END SUBROUTINE GetSecondaryElementSide_Edge
!
!
!
 SUBROUTINE SetStart_Edge( thisEdge, start )
 ! S/R SetStart
 !  
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS( Edge ), INTENT(inout) :: thisEdge
   INTEGER, INTENT(in)          :: start
    
      thisEdge % start = start
      
 END SUBROUTINE SetStart_Edge
!
!
!
 SUBROUTINE GetStart_Edge( thisEdge, start )
 ! S/R GetStart
 !  
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS( Edge ), INTENT(in) :: thisEdge
   INTEGER, INTENT(out)      :: start

      start = thisEdge % start

 END SUBROUTINE GetStart_Edge
!
!
!
 SUBROUTINE SetIncrement_Edge( thisEdge, inc )
 ! S/R SetIncrement
 !  
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS( Edge ), INTENT(inout) :: thisEdge
   INTEGER, INTENT(in)          :: inc

      thisEdge % inc = inc

 END SUBROUTINE SetIncrement_Edge
!
!
!
 SUBROUTINE GetIncrement_Edge( thisEdge, inc )
 ! S/R GetIncrement
 !  
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS( Edge ), INTENT(in) :: thisEdge
   INTEGER, INTENT(out)      :: inc

      inc = thisEdge % inc

 END SUBROUTINE GetIncrement_Edge
!
!
!==================================================================================================!
!----------------------------------- Linked List Routines  ----------------------------------------!
!==================================================================================================!
!
!
 FUNCTION ListIsEmpty_EdgeList( myList ) RESULT( TorF )
 ! Function ListIsEmpty
 !
 !  This function checks if "myList" is an empty list. A logical is returned
 !  
 ! =============================================================================================== !
 ! DECLARATIONS
  IMPLICIT NONE
  CLASS( EdgeList ) :: myList
  LOGICAL           :: TorF

     TorF = .NOT.( ASSOCIATED( myList % head  ) )
     
 END FUNCTION ListIsEmpty_EdgeList
!
!
!
 SUBROUTINE AddToList_EdgeList( myList,  nodeIDs, elementIDs, elementSides, key, start, inc, inKey )
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
  CLASS( EdgeList ) :: myList
  INTEGER           :: nodeIDs(1:2), elementIDs(1:2), elementSides(1:2)
  INTEGER           :: key, start, inc
  INTEGER, OPTIONAL :: inKey
  ! LOCAL
  LOGICAL               :: isAssigned
  TYPE( Edge ), POINTER :: previous
  INTEGER               :: allocationStatus

     ! Check to see if this list is empty
     IF( myList % ListIsEmpty() )THEN
     
        ALLOCATE( myList % head, STAT = allocationStatus )
        IF( allocationStatus /=0 )THEN
           PRINT*, 'MODULE EdgeListClass.f90 : S/R AddToList : Memory not allocated for next entry in list.'
           ! An exception handler should be built to handle these problems
        ENDIF      
      
        ! Point the current position to the head
        myList % current => myList % head
        ! Set the data
        CALL myList % SetData(  nodeIDs, elementIDs, elementSides, key, start, inc )
        
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
           PRINT*, 'MODULE EdgeListClass.f90 : S/R AddToList : Memory not allocated for next entry in list.'
           ! An exception handler should be built to handle these problems
        ENDIF      
        
        !  Temporarily hold onto the tail in case we need the key (if inKey is not present)  
        previous => myList % tail
        ! Reassign the tail
        myList % tail => myList % tail % next
        
        ! Set the current to the tail
        myList % current => myList % tail
  
        ! Fill in the data
        CALL myList % SetData(  nodeIDs, elementIDs, elementSides, key, start, inc )
        
        ! Fill in the key information
        IF( PRESENT(inKey) )THEN
           CALL myList % SetKey( inKey )
        ELSE
           CALL myList % SetKey( previous % key + keyInc )
        ENDIF
        
        ! Point the next to null and the tail to current
        myList % current % next => NULL( )
        
     ENDIF

 END SUBROUTINE AddToList_EdgeList
!
!
!
 SUBROUTINE RemoveCurrent_EdgeList( myList )
 ! S/R RemoveCurrent
 !
 !  This subroutine removes the current item in the linked list and patches together the previous
 !  and next items.
 !  
 ! =============================================================================================== !
  IMPLICIT NONE
  CLASS( EdgeList ) :: myList
  ! LOCAL
  TYPE( Edge ), POINTER :: previous, pNext
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

 END SUBROUTINE RemoveCurrent_EdgeList
!
!
!
 SUBROUTINE MoveToNext_EdgeList( myList )
 ! S/R MoveToNext
 !
 ! =============================================================================================== !
 ! DECLARATIONS
  IMPLICIT NONE
  CLASS( EdgeList ) :: myList

     myList % current => myList % current % next

 END SUBROUTINE MoveToNext_EdgeList
!
!
!
  SUBROUTINE MoveToHead_EdgeList( myList )
 ! S/R MoveToHead
 !
 !  
 ! ========================================================================== !
  IMPLICIT NONE
  CLASS( EdgeList ) :: myList

     myList % current => myList % head

 END SUBROUTINE MoveToHead_EdgeList
!
!
!
 SUBROUTINE MoveToTail_EdgeList( myList )
 ! S/R MoveToTail
 !
 !  
 ! ========================================================================== !
  IMPLICIT NONE
  CLASS( EdgeList ) :: myList

     myList % current => myList % tail

 END SUBROUTINE MoveToTail_EdgeList
!
!
!
 
  SUBROUTINE GetCount_EdgeList( myList, numberOfEdges )
 ! S/R GetCount
 !
 !
 ! ===========================================================================================s==== !
  IMPLICIT NONE
  CLASS( EdgeList )    :: myList
  INTEGER, INTENT(out) :: numberOfEdges

     numberOfEdges = 0 ! Initialize the number of list items
     ! Check to see if this list is empty
     IF( myList % ListIsEmpty() )THEN
     
        PRINT*, 'Module EdgeClass.f90 : S/R GetCount : List is empty.'
        RETURN
        
     ELSE ! the list is not empty
    
        CALL myList % MoveToHead( ) ! Rewind the list
        
        DO WHILE( ASSOCIATED( myList % current ) )
        
           numberOfEdges = numberOfEdges + 1
           CALL myList % moveToNext( )

        ENDDO

     ENDIF

 END SUBROUTINE GetCount_EdgeList
! SUBROUTINE PrintList( myList )
 ! S/R PrintList
 !
 ! =============================================================================================== !
 ! DECLARATIONS
!  IMPLICIT NONE
!  CLASS( EdgeList ) :: myList

!     myList % current => myList % head

!    PRINT*, '          Data        Key'
!    DO WHILE( ASSOCIATED( myList % current ) )
 
!        PRINT*, myList % current % listData, myList % current % key

!        CALL myList % MoveToNext()

!     ENDDO

! END SUBROUTINE PrintList
END MODULE EdgeClass
