MODULE FaceClass
! FaceClass.f90 (v2.1 - 21 Jan 2016) * New with version 2.1
! 
! schoonover.numerics@gmail.com
! o (ver 2.1) Jan 2016
!
!
! The face class stores identification information for it's corner nodes and 
! neighboring elements. In this sense, faces are viewed as part of the book-keeping side of the mesh 
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


   TYPE Face
   
      INTEGER, PRIVATE      :: key               ! The edge key
      INTEGER, PRIVATE      :: nodeIDs(1:4)      ! Node IDs which start and terminate this edge
      INTEGER, PRIVATE      :: elementIDs(1:2)   ! Neighboring elements IDs across the edge
      INTEGER, PRIVATE      :: elementSides(1:2) ! Local side IDs for the neighboring elements
      INTEGER, PRIVATE      :: start, inc        ! Loop start and increment for the secondary element side
      LOGICAL, PRIVATE      :: swapDimensions    ! A flag for swapping the dimensions of the secondary
                                                 ! element's boundary solution/flux/etc.
                                                
      TYPE( Face ), POINTER :: next => NULL( )


      CONTAINS
 
      PROCEDURE :: Build => Build_Face

      PROCEDURE :: SetKey => SetKey_Face
      PROCEDURE :: GetKey => GetKey_Face
      PROCEDURE :: SetData => SetData_Face
      PROCEDURE :: GetData => GetData_Face
      PROCEDURE :: SetNodeIDs => SetNodeIDs_Face
      PROCEDURE :: GetNodeIDs => GetNodeIDs_Face
      PROCEDURE :: SetElementIDs => SetElementIDs_Face
      PROCEDURE :: GetElementIDs => GetElementIDs_Face
      PROCEDURE :: SetPrimaryElementID => SetPrimaryElementID_Face
      PROCEDURE :: GetPrimaryElementID => GetPrimaryElementID_Face
      PROCEDURE :: SetSecondaryElementID => SetSecondaryElementID_Face
      PROCEDURE :: GetSecondaryElementID => GetSecondaryElementID_Face
      PROCEDURE :: SetElementSides => SetElementSides_Face
      PROCEDURE :: GetElementSides => GetElementSides_Face
      PROCEDURE :: SetPrimaryElementSide => SetPrimaryElementSide_Face
      PROCEDURE :: GetPrimaryElementSide => GetPrimaryElementSide_Face
      PROCEDURE :: SetSecondaryElementSide => SetSecondaryElementSide_Face
      PROCEDURE :: GetSecondaryElementSide => GetSecondaryElementSide_Face
      PROCEDURE :: SetStart => SetStart_Face
      PROCEDURE :: GetStart => GetStart_Face
      PROCEDURE :: SetIncrement => SetIncrement_Face
      PROCEDURE :: GetIncrement => GetIncrement_Face
      PROCEDURE :: SetSwapDimensions => SetSwapDimensions_Face
      PROCEDURE :: GetSwapDimensions => GetSwapDimensions_Face

   END TYPE Face

   TYPE FaceList      
      TYPE( Face ), POINTER :: head, tail, current

      CONTAINS

      PROCEDURE :: Build => Build_FaceList
      PROCEDURE :: Trash => Trash_FaceList
      
      
      PROCEDURE :: GetData => GetCurrentData_FaceList
      PROCEDURE :: SetData => SetCurrentData_FaceList
      PROCEDURE :: GetKey => GetCurrentKey_FaceList
      PROCEDURE :: SetKey => SetCurrentKey_FaceList
      
      PROCEDURE :: ListIsEmpty => ListIsEmpty_FaceList
      PROCEDURE :: AddToList => AddToList_FaceList
      PROCEDURE :: RemoveCurrent => RemoveCurrent_FaceList
      PROCEDURE :: MoveToHead => MoveToHead_FaceList
      PROCEDURE :: MoveToNext => MoveToNext_FaceList
      PROCEDURE :: MoveToTail => MoveToTail_FaceList
      PROCEDURE :: GetCount => GetCount_FaceList
      
   !   PROCEDURE :: PrintList

   END TYPE FaceList
   
! Setting up some parameters pertaining to this module
 INTEGER, PARAMETER, PRIVATE :: keyInc   = 1 ! The default increment in the Record Key
 INTEGER, PARAMETER, PRIVATE :: keyStart = 1 ! The default starting Record key
 LOGICAL, PARAMETER, PRIVATE :: defaultSwapDim = .FALSE.

 CONTAINS
!
!
!==================================================================================================!
!------------------------------- Manual Constructors/Destructors ----------------------------------!
!==================================================================================================!
!
!
SUBROUTINE Build_Face( myFace )
 ! S/R Build_Face
 !
 !
 ! =============================================================================================== !
 ! DECLARATIONS
  IMPLICIT NONE
  CLASS( Face ) :: myFace

     myFace % nodeIDs        = NO_NORMAL_FLOW
     myFace % elementIDs     = NO_NORMAL_FLOW
     myFace % elementSides   = 0
     myFace % start          = 1 
     myFace % inc            = 1
     myFace % swapDimensions = defaultSwapDim
     myFace % next => NULL()
     
 END SUBROUTINE Build_Face 
!
!
!
 SUBROUTINE Build_FaceList( myList )
 ! S/R Build_FaceList
 !
 !  This subroutine creates a linked-list by allocating space for the list
 !  head. The "next" attribute of the list head is nullified. 
 !  The tail is pointed to NULL and the current position in the list
 !  is pointed to the head. 
 ! 
 !   INPUT/OUTPUT : 
 !      CLASS( FaceList) :: myList
 !
 ! =============================================================================================== !
 ! DECLARATIONS
  IMPLICIT NONE
  CLASS( FaceList ) :: myList

    
     myList % head => NULL( )

     ! Point the tail to null
     myList % tail => NULL()

     ! Set the current position to Null
     myList % current => NULL( )
  
 END SUBROUTINE Build_FaceList
!
!
!  
 SUBROUTINE Trash_FaceList( myList )
 ! S/R Trash_FaceList
 ! 
 !  This subroutine traverses through the linked list and frees the memory
 !  that was allocated for the entries in the linked list.  
 !
 ! =============================================================================================== !
 ! DECLARATIONS
  IMPLICIT NONE
  CLASS( FaceList ) :: myList
  ! LOCAL
  TYPE( Face ), POINTER :: pNext

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
  
 END SUBROUTINE Trash_FaceList
!
!
!==================================================================================================!
!--------------------------------------- Accessors ------------------------------------------------!
!==================================================================================================!
!
!
 SUBROUTINE SetCurrentData_FaceList( myList, nodeIDs, elementIDs, elementSides, key, start, inc )
 ! S/R SetCurrentData
 !
 !
 ! ========================================================================== !
 IMPLICIT NONE
  CLASS( FaceList )   :: myList
  INTEGER, INTENT(in) :: nodeIDs(1:4), elementIDs(1:2), elementSides(1:2)
  INTEGER, INTENT(in) :: key, start, inc

     CALL myList % current % SetData( nodeIDs, elementIDs, elementSides, key, start, inc )

 END SUBROUTINE SetCurrentData_FaceList
!
!
!
 SUBROUTINE GetCurrentData_FaceList( myList, nodeIDs, elementIDs, elementSides, key, start, inc )
 ! S/R GetCurrentData
 !
 !
 ! ========================================================================== !
 IMPLICIT NONE
  CLASS( FaceList )    :: myList
  INTEGER, INTENT(out) :: nodeIDs(1:4), elementIDs(1:2), elementSides(1:2)
  INTEGER, INTENT(out) :: key, start, inc

     CALL myList % current % GetData( nodeIDs, elementIDs, elementSides, key, start, inc )

 END SUBROUTINE GetCurrentData_FaceList
!
!
!
 SUBROUTINE SetCurrentKey_FaceList( myList, key )
 ! S/R SetCurrentKey
 !
 !
 ! ========================================================================== !
 IMPLICIT NONE
  CLASS( FaceList )   :: myList
  INTEGER, INTENT(in) :: key

     CALL myList % current % SetKey( key )

 END SUBROUTINE SetCurrentKey_FaceList
!
!
!
 SUBROUTINE GetCurrentKey_FaceList( myList, key )
 ! S/R GetCurrentKey
 !
 !
 ! ========================================================================== !
  IMPLICIT NONE
  CLASS( FaceList )    :: myList
  INTEGER, INTENT(out) :: key

     CALL myList % current % GetKey( key )

 END SUBROUTINE GetCurrentKey_FaceList
!
!
!
 SUBROUTINE SetKey_Face( thisFace, key )
 ! S/R SetKey
 !  
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS( Face ), INTENT(inout) :: thisFace
   INTEGER, INTENT(in)          :: key

      thisFace % key = key      

 END SUBROUTINE SetKey_Face
!
!
!
 SUBROUTINE GetKey_Face( thisFace, key )
 ! S/R GetKey
 !  
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS( Face ), INTENT(in) :: thisFace
   INTEGER, INTENT(out)      :: key

      key = thisFace % key

 END SUBROUTINE GetKey_Face
!
!
!
 SUBROUTINE SetData_Face( thisFace, nodeIDs, elementIDs, elementSides, key, start, inc )
 ! S/R SetData
 !  
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS( Face ), INTENT(inout) :: thisFace
   INTEGER, INTENT(in)          :: nodeIDs(1:4), elementIDs(1:2), elementSides(1:2)
   INTEGER, INTENT(in)          :: key, start, inc

      CALL thisFace % SetNodeIDs( nodeIDs )
      CALL thisFace % SetElementIDs( elementIDs )
      CALL thisFace % SetElementSides( elementSides )
      CALL thisFace % SetKey( key )
      CALL thisFace % SetStart( start )
      CALL thisFace % SetIncrement( inc )

 END SUBROUTINE SetData_Face
!
!
!
 SUBROUTINE GetData_Face( thisFace, nodeIDs, elementIDs, elementSides, key, start, inc )
 ! S/R GetData
 !  
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS( Face ), INTENT(in) :: thisFace
   INTEGER, INTENT(out)      :: nodeIDs(1:4), elementIDs(1:2), elementSides(1:2)
   INTEGER, INTENT(out)      :: key, start, inc

      CALL thisFace % GetNodeIDs( nodeIDs )
      CALL thisFace % GetElementIDs( elementIDs )
      CALL thisFace % GetElementSides( elementSides )
      CALL thisFace % GetKey( key )
      CALL thisFace % GetStart( start )
      CALL thisFace % GetIncrement( inc )

 END SUBROUTINE GetData_Face
!
!
!
 SUBROUTINE SetNodeIDs_Face( thisFace, nodeIDs )
 ! S/R SetNodeIDs
 !  
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS( Face ), INTENT(inout) :: thisFace
   INTEGER, INTENT(in)          :: nodeIDs(1:4)

      thisFace % nodeIDs = nodeIDs

 END SUBROUTINE SetNodeIDs_Face
!
!
!
 SUBROUTINE GetNodeIDs_Face( thisFace, nodeIDs )
 ! S/R GetNodeIDs
 !  
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS( Face ), INTENT(in) :: thisFace
   INTEGER, INTENT(out)      :: nodeIDs(1:4)

      nodeIDs = thisFace % nodeIDs

 END SUBROUTINE GetNodeIDs_Face
!
!
!
 SUBROUTINE SetElementIDs_Face( thisFace, elementIDs )
 ! S/R SetElementIDs
 !  
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS( Face ), INTENT(inout) :: thisFace
   INTEGER, INTENT(in)          :: elementIDs(1:2)

      thisFace % elementIDs = elementIDs

 END SUBROUTINE SetElementIDs_Face
!
!
!
 SUBROUTINE GetElementIDs_Face( thisFace, elementIDs )
 ! S/R GetElementIDs
 !  
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS( Face ), INTENT(in) :: thisFace
   INTEGER, INTENT(out)      :: elementIDs(1:2)

      elementIDs = thisFace % elementIDs

 END SUBROUTINE GetElementIDs_Face
!
!
!
 SUBROUTINE SetPrimaryElementID_Face( thisFace, elementID )
 ! S/R SetPrimaryElementID
 !  
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS( Face ), INTENT(inout) :: thisFace
   INTEGER, INTENT(in)          :: elementID

      thisFace % elementIDs(1) = elementID
      
 END SUBROUTINE SetPrimaryElementID_Face
!
!
!
 SUBROUTINE GetPrimaryElementID_Face( thisFace, elementID )
 ! S/R GetPrimaryElementID
 !  
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS( Face ), INTENT(in) :: thisFace
   INTEGER, INTENT(out)      :: elementID

      elementID = thisFace % elementIDs(1)

 END SUBROUTINE GetPrimaryElementID_Face
!
!
!
 SUBROUTINE SetSecondaryElementID_Face( thisFace, elementID )
 ! S/R SetSecondaryElementID
 !  
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS( Face ), INTENT(inout) :: thisFace
   INTEGER, INTENT(in)          :: elementID

      thisFace % elementIDs(2) = elementID
      
 END SUBROUTINE SetSecondaryElementID_Face
!
!
!
 SUBROUTINE GetSecondaryElementID_Face( thisFace, elementID )
 ! S/R GetSecondaryElementID
 !  
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS( Face ), INTENT(in) :: thisFace
   INTEGER, INTENT(out)      :: elementID

      elementID = thisFace % elementIDs(2)

 END SUBROUTINE GetSecondaryElementID_Face
!
!
!
 SUBROUTINE SetElementSides_Face( thisFace, elementSides )
 ! S/R SetElementSides
 !  
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS( Face ), INTENT(inout) :: thisFace
   INTEGER, INTENT(in)          :: elementSides(1:2)

      thisFace % elementSides = elementSides
       
 END SUBROUTINE SetElementSides_Face
!
!
!
 SUBROUTINE GetElementSides_Face( thisFace, elementSides )
 ! S/R GetElementSides
 !  
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS( Face ), INTENT(in) :: thisFace
   INTEGER, INTENT(out)      :: elementSides(1:2)

      elementSides = thisFace % elementSides

 END SUBROUTINE GetElementSides_Face
!
!
!
 SUBROUTINE SetPrimaryElementSide_Face( thisFace, elementSide )
 ! S/R SetPrimaryElementSide
 !  
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS( Face ), INTENT(inout) :: thisFace
   INTEGER, INTENT(in)          :: elementSide

      thisFace % elementSides(1) = elementSide
      
 END SUBROUTINE SetPrimaryElementSide_Face
!
!
!
 SUBROUTINE GetPrimaryElementSide_Face( thisFace, elementSide )
 ! S/R GetPrimaryElementSide
 !  
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS( Face ), INTENT(in) :: thisFace
   INTEGER, INTENT(out)      :: elementSide

      elementSide = thisFace % elementSides(1)

 END SUBROUTINE GetPrimaryElementSide_Face
!
!
!
 SUBROUTINE SetSecondaryElementSide_Face( thisFace, elementSide )
 ! S/R SetSecondaryElementSide
 !  
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS( Face ), INTENT(inout) :: thisFace
   INTEGER, INTENT(in)          :: elementSide

      thisFace % elementSides(2) = elementSide
      
 END SUBROUTINE SetSecondaryElementSide_Face
!
!
!
 SUBROUTINE GetSecondaryElementSide_Face( thisFace, elementSide )
 ! S/R GetSecondaryElementSide
 !  
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS( Face ), INTENT(in) :: thisFace
   INTEGER, INTENT(out)      :: elementSide

      elementSide = thisFace % elementSides(2)

 END SUBROUTINE GetSecondaryElementSide_Face
!
!
!
 SUBROUTINE SetStart_Face( thisFace, start )
 ! S/R SetStart
 !  
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS( Face ), INTENT(inout) :: thisFace
   INTEGER, INTENT(in)          :: start
    
      thisFace % start = start
      
 END SUBROUTINE SetStart_Face
!
!
!
 SUBROUTINE GetStart_Face( thisFace, start )
 ! S/R GetStart
 !  
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS( Face ), INTENT(in) :: thisFace
   INTEGER, INTENT(out)      :: start

      start = thisFace % start

 END SUBROUTINE GetStart_Face
!
!
!
 SUBROUTINE SetIncrement_Face( thisFace, inc )
 ! S/R SetIncrement
 !  
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS( Face ), INTENT(inout) :: thisFace
   INTEGER, INTENT(in)          :: inc

      thisFace % inc = inc

 END SUBROUTINE SetIncrement_Face
!
!
!
 SUBROUTINE GetIncrement_Face( thisFace, inc )
 ! S/R GetIncrement
 !  
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS( Face ), INTENT(in) :: thisFace
   INTEGER, INTENT(out)      :: inc

      inc = thisFace % inc

 END SUBROUTINE GetIncrement_Face
!
!
!
 SUBROUTINE SetSwapDimensions_Face( thisFace, swap )
 ! S/R SetSwapDimensions
 !  
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS( Face ), INTENT(inout) :: thisFace
   LOGICAL, INTENT(in)          :: swap

      thisFace % swapDimensions = swap

 END SUBROUTINE SetSwapDimensions_Face
!
!
!
 SUBROUTINE GetSwapDimensions_Face( thisFace, swap )
 ! S/R GetSwapDimensions
 !  
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS( Face ), INTENT(in) :: thisFace
   LOGICAL, INTENT(out)      :: swap

      swap = thisFace % swapDimensions

 END SUBROUTINE GetSwapDimensions_Face
!
!
!==================================================================================================!
!----------------------------------- Linked List Routines  ----------------------------------------!
!==================================================================================================!
!
!
 FUNCTION ListIsEmpty_FaceList( myList ) RESULT( TorF )
 ! Function ListIsEmpty
 !
 !  This function checks if "myList" is an empty list. A logical is returned
 !  
 ! =============================================================================================== !
 ! DECLARATIONS
  IMPLICIT NONE
  CLASS( FaceList ) :: myList
  LOGICAL           :: TorF

     TorF = .NOT.( ASSOCIATED( myList % head  ) )
     
 END FUNCTION ListIsEmpty_FaceList
!
!
!
 SUBROUTINE AddToList_FaceList( myList,  nodeIDs, elementIDs, elementSides, key, start, inc, inKey )
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
  CLASS( FaceList ) :: myList
  INTEGER           :: nodeIDs(1:4), elementIDs(1:2), elementSides(1:2)
  INTEGER           :: key, start, inc
  INTEGER, OPTIONAL :: inKey
  ! LOCAL
  LOGICAL               :: isAssigned
  TYPE( Face ), POINTER :: previous
  INTEGER               :: allocationStatus

     ! Check to see if this list is empty
     IF( myList % ListIsEmpty() )THEN
     
        ALLOCATE( myList % head, STAT = allocationStatus )
        IF( allocationStatus /=0 )THEN
           PRINT*, 'MODULE FaceListClass.f90 : S/R AddToList : Memory not allocated for next entry in list.'
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
           PRINT*, 'MODULE FaceListClass.f90 : S/R AddToList : Memory not allocated for next entry in list.'
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

 END SUBROUTINE AddToList_FaceList
!
!
!
 SUBROUTINE RemoveCurrent_FaceList( myList )
 ! S/R RemoveCurrent
 !
 !  This subroutine removes the current item in the linked list and patches together the previous
 !  and next items.
 !  
 ! =============================================================================================== !
  IMPLICIT NONE
  CLASS( FaceList ) :: myList
  ! LOCAL
  TYPE( Face ), POINTER :: previous, pNext
  INTEGER               :: currentKey, thisKey

     CALL myList % GetKey( currentKey )
     
     ! Check to see if this list is empty
     IF( myList % ListIsEmpty() )THEN
     
        PRINT*, 'Module FaceClass.f90 : S/R RemoveCurrent : List is empty. Nothing to remove'
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

 END SUBROUTINE RemoveCurrent_FaceList
!
!
!
 SUBROUTINE MoveToNext_FaceList( myList )
 ! S/R MoveToNext
 !
 ! =============================================================================================== !
 ! DECLARATIONS
  IMPLICIT NONE
  CLASS( FaceList ) :: myList

     myList % current => myList % current % next

 END SUBROUTINE MoveToNext_FaceList
!
!
!
  SUBROUTINE MoveToHead_FaceList( myList )
 ! S/R MoveToHead
 !
 !  
 ! ========================================================================== !
  IMPLICIT NONE
  CLASS( FaceList ) :: myList

     myList % current => myList % head

 END SUBROUTINE MoveToHead_FaceList
!
!
!
 SUBROUTINE MoveToTail_FaceList( myList )
 ! S/R MoveToTail
 !
 !  
 ! ========================================================================== !
  IMPLICIT NONE
  CLASS( FaceList ) :: myList

     myList % current => myList % tail

 END SUBROUTINE MoveToTail_FaceList
!
!
!
 
  SUBROUTINE GetCount_FaceList( myList, numberOfFaces )
 ! S/R GetCount
 !
 !
 ! ===========================================================================================s==== !
  IMPLICIT NONE
  CLASS( FaceList )    :: myList
  INTEGER, INTENT(out) :: numberOfFaces

     numberOfFaces = 0 ! Initialize the number of list items
     ! Check to see if this list is empty
     IF( myList % ListIsEmpty() )THEN
     
        PRINT*, 'Module FaceClass.f90 : S/R GetCount : List is empty.'
        RETURN
        
     ELSE ! the list is not empty
    
        CALL myList % MoveToHead( ) ! Rewind the list
        
        DO WHILE( ASSOCIATED( myList % current ) )
        
           numberOfFaces = numberOfFaces + 1
           CALL myList % moveToNext( )

        ENDDO

     ENDIF

 END SUBROUTINE GetCount_FaceList
! SUBROUTINE PrintList( myList )
 ! S/R PrintList
 !
 ! =============================================================================================== !
 ! DECLARATIONS
!  IMPLICIT NONE
!  CLASS( FaceList ) :: myList

!     myList % current => myList % head

!    PRINT*, '          Data        Key'
!    DO WHILE( ASSOCIATED( myList % current ) )
 
!        PRINT*, myList % current % listData, myList % current % key

!        CALL myList % MoveToNext()

!     ENDDO

! END SUBROUTINE PrintList
END MODULE FaceClass
