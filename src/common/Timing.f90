MODULE Timing
!
! Timing.f90 (v2.1 - 12 Nov. 2015)
!
! Joe Schoonover
!
! schoonover.numerics@gmail.com
!
!
!
! Future plans for this MODULE include routines for calculating more sophisticated timing
! statistics that can be used for benchmarking software.
! 
! ================================================================================================ !
!

USE ModelPrecision
USE ConstantsDictionary
USE CommonRoutines


 IMPLICIT NONE

    INTEGER, PARAMETER :: timerNameLength = 30
    ! ============================================================================================ !
    !
    ! RoutineTimer defines a structure for a timer which has a start and stop time, and integer ID,
    ! a name for the timer and flags for determining whether the timer has been started or stopped.
    !
    ! This is a single record for a linked list of RoutineTimers. In a linked list, this allows for
    ! easy insertion of timers into your code without having to worry about knowing the number of 
    ! timers needed a'priori.
    !
    ! ============================================================================================ !
    TYPE RoutineTimer
       LOGICAL, PRIVATE                    :: started = .FALSE., stopped = .FALSE.
       REAL(prec), PRIVATE                 :: startTime = ZERO
       REAL(prec), PRIVATE                 :: stopTime = ZERO
       REAL(prec), PRIVATE                 :: accumulatedTime = ZERO
       REAL(prec), PRIVATE                 :: nObs
       INTEGER, PRIVATE                    :: timerID
       CHARACTER(timerNameLength), PRIVATE :: whatYourTiming
       TYPE( RoutineTimer ), POINTER       :: next

       CONTAINS

       PROCEDURE :: SetName => SetTimerName
       PROCEDURE :: GetName => GetTimerName
       PROCEDURE :: SetAccumulatedTime
       PROCEDURE :: GetAccumulatedTime
  
       PROCEDURE :: StartTimer
       PROCEDURE :: StopTimer
       PROCEDURE :: ElapsedTime 

    END TYPE RoutineTimer

    TYPE MultiTimers
        TYPE( RoutineTimer ), POINTER :: head, current, tail

        CONTAINS

        PROCEDURE :: Build => Build_MultiTimers
        PROCEDURE :: Trash => Trash_MultiTimers

        PROCEDURE :: SetName => SetThisTimerName
        PROCEDURE :: GetName => GetThisTimerName
        
        PROCEDURE :: ThereAreNoTimers
        PROCEDURE :: AddTimer 
        PROCEDURE :: PointToTimer
        PROCEDURE :: MoveToNext => MoveToNextTimer
        
        PROCEDURE :: StartThisTimer
        PROCEDURE :: StopThisTimer
        PROCEDURE :: AccumulateTimings

        PROCEDURE :: Write_MultiTimers
        
    END TYPE MultiTimers

 CONTAINS
!
!
!=========================================================================!
!---------------- CONSTRUCTOR/DESTRUCTOR ROUTINES ------------------------!
!=========================================================================!
!
!
 SUBROUTINE Build_MultiTimers( theTimers  )
 !
 !
 ! ================================================ !
 ! DECLARATION
   IMPLICIT NONE
   CLASS( MultiTimers ) :: theTimers
   
     theTimers % head => NULL( )

     ! Point the tail to null
     theTimers % tail => NULL()

     ! Set the current position to Null
     theTimers % current => NULL( )


 END SUBROUTINE Build_MultiTimers
!
!
!
 SUBROUTINE Trash_MultiTimers( theTimers )
 !
 !
 ! ================================================ !
 ! DECLARATION
   IMPLICIT NONE
   CLASS( MultiTimers ) :: theTimers
   ! LOCAL 
   TYPE( RoutineTimer ), POINTER :: pNext
   
     ! Set the current position of the list to the head
     theTimers % current => theTimers % head
     
     ! Scroll through the list until the current position is nullified
     DO WHILE ( ASSOCIATED( theTimers % current ) )

        ! temporarily point to the next in the list
        pNext => theTimers % current % next 

        ! Deallocate memory pointed to by current position
        DEALLOCATE( theTimers % current ) 

        ! Update current position
        theTimers % current => pNext 

     ENDDO
      

 END SUBROUTINE Trash_MultiTimers
!
!
!=========================================================================!
!-------------------------- ACCESSOR ROUTINES ----------------------------!
!=========================================================================!
!
! 
 SUBROUTINE SetThisTimerName( theTimers, nameInput, tID )
 !
 !
 ! ================================================= !
 ! DELCARATION
   IMPLICIT NONE
   CLASS( MultiTimers )          :: theTimers
   CHARACTER(*), INTENT(IN)      :: nameInput
   INTEGER, OPTIONAL, INTENT(IN) :: tID
   ! LOCAL
   
      IF( PRESENT(tID) )THEN
      
         CALL theTimers % PointToTimer( tID ) ! Point the "current" pointer to the RoutineTimer with 
                                              ! timerID = tID.
      ENDIF
      
      CALL theTimers % current % SetName( nameInput )

 END SUBROUTINE SetThisTimerName
!
!
!
 SUBROUTINE GetThisTimerName( theTimers, nameOutput, tID )
 !
 !
 ! ================================================= !
 ! DELCARATION
   IMPLICIT NONE
   CLASS( MultiTimers )          :: theTimers
   CHARACTER(*), INTENT(OUT)     :: nameOutput
   INTEGER, OPTIONAL, INTENT(IN) :: tID
   ! LOCAL
   
      IF( PRESENT(tID) )THEN
      
         CALL theTimers % PointToTimer( tID ) ! Point the "current" pointer to the RoutineTimer with 
                                              ! timerID = tID.
      ENDIF
      
      CALL theTimers % current % GetName( nameOutput )

 END SUBROUTINE GetThisTimerName
!
!
!
 SUBROUTINE SetTimerName( theTimer, nameInput )
 !
 !
 ! ================================================= !
 ! DELCARATION
   IMPLICIT NONE
   CLASS( RoutineTimer ), INTENT(INOUT) :: theTimer
   CHARACTER(*), INTENT(IN)             :: nameInput


      theTimer % whatYourTiming = nameInput

 END SUBROUTINE SetTimerName
!
!
!
 SUBROUTINE GetTimerName( theTimer, nameOutput )
 !
 !
 ! ================================================= !
 ! DELCARATION
   IMPLICIT NONE
   CLASS( RoutineTimer ), INTENT(IN) :: theTimer
   CHARACTER(*), INTENT(OUT)         :: nameOutput

      nameOutput = theTimer % whatYourTiming

 END SUBROUTINE GetTimerName
!
!
!
 SUBROUTINE SetAccumulatedTime( theTimer, accTime )
 !
 !
 ! ================================================= !
 ! DELCARATION
   IMPLICIT NONE
   CLASS( RoutineTimer ), INTENT(INOUT)  :: theTimer
   REAL(PREC), INTENT(IN)                :: accTime


      theTimer % accumulatedTime = accTime

 END SUBROUTINE SetAccumulatedTime
!
!
!
 FUNCTION GetAccumulatedTime( theTimer ) RESULT( accTime )
 !
 !
 ! ================================================= !
 ! DELCARATION
   IMPLICIT NONE
   CLASS( RoutineTimer ) :: theTimer
   REAL(prec)            :: accTime

      accTime = theTimer % accumulatedTime

 END FUNCTION GetAccumulatedTime
!
!
!=========================================================================!
!------------------ Linked-List Type Operations ------------------- ------!
!=========================================================================!
!
! 
 FUNCTION ThereAreNoTimers( theTimers ) RESULT( TorF)
 ! Function ThereAreNoTimers
 !
 !  
 ! ========================================================================== !
  IMPLICIT NONE
  CLASS( MultiTimers ) :: theTimers
  LOGICAL              :: TorF

     TorF = .NOT.( ASSOCIATED( theTimers % head  ) )
     
 END FUNCTION ThereAreNoTimers
!
!
!
 SUBROUTINE AddTimer( theTimers, timername, inKey )
 ! S/R AddTimer
 !
 !  
 ! ========================================================================== !
  IMPLICIT NONE
  CLASS( MultiTimers ), INTENT(inout) :: theTimers
  CHARACTER(*)                        :: timername
  INTEGER                             :: inKey
  ! LOCAL
  INTEGER :: allocationStatus

     ! Check to see if this list is empty
     IF( theTimers % ThereAreNoTimers() )THEN
     
        ALLOCATE( theTimers % head, STAT = allocationStatus )
        IF( allocationStatus /=0 )THEN
           PRINT*, 'MODULE Timing.f90 : S/R AddTimer : Memory not allocated for next entry in list.'
           ! An exception handler should be built to handle these problems
        ENDIF      
      
        ! Point the current position to the head
        theTimers % current => theTimers % head
        ! Set the data
        CALL theTimers % SetName( timername )
        
        theTimers % current % timerID = inkey
        
        ! Point the next to null and the tail to current
        theTimers % current % next => NULL( )
        theTimers % tail => theTimers % current
        
     ELSE ! the list is not empty
    
        ! Then we allocate space for the next item in the list    
        ALLOCATE( theTimers % tail % next, STAT = allocationStatus )
        IF( allocationStatus /=0 )THEN
           PRINT*, 'MODULE Timing.f90 : S/R AddTimer : Memory not allocated for next entry in list.'
           ! An exception handler should be built to handle these problems
        ENDIF      
        
        ! Reassign the tail
        theTimers % tail => theTimers % tail % next
        
        ! Set the current to the tail
        theTimers % current => theTimers % tail
  
        ! Fill in the data
        CALL theTimers % SetName( timername )
        
        ! Fill in the key information
        theTimers % current % timerID = inkey

        ! Point the next to null and the tail to current
        theTimers % current % next => NULL( )
        
     ENDIF

 END SUBROUTINE AddTimer
!
!
!
 SUBROUTINE PointToTimer( theTimers, tID )
 ! S/R PointToTimer
 !
 ! Points the current timer to the timer with timerID = tID
 !
 ! ====================================================================== !
   IMPLICIT NONE
   CLASS( MultiTimers ) :: theTimers
   INTEGER, INTENT(IN)  :: tID
   
   
      theTimers % current => theTimers % head ! Point to the head of the list
      
      DO WHILE(ASSOCIATED(theTimers % current))
      
         IF( theTimers % current % timerID == tID )EXIT
         
         CALL theTimers % MoveToNext( )
      
      ENDDO
      
      
 END SUBROUTINE PointToTimer
!
!
!
 SUBROUTINE MoveToNextTimer( theTimers )
 ! S/R MoveToNextTimer
 !
 !  
 ! ========================================================================== !
  IMPLICIT NONE
  CLASS( MultiTimers ) :: theTimers

     theTimers % current => theTimers % current % next

 END SUBROUTINE MoveToNextTimer
!
!
!=========================================================================!
!------------------- Data Structure Operations ---------------------------!
!=========================================================================!
!
!
 SUBROUTINE StartThisTimer( theTimers, tID )
 !
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS( MultiTimers ), INTENT(inout) :: theTimers
   INTEGER, INTENT(in)                 :: tID

      CALL theTimers % PointToTimer( tID )
      CALL theTimers % current % StartTimer( )

 END SUBROUTINE StartThisTimer
!
!
!
 SUBROUTINE StopThisTimer( theTimers, tID )
 !
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS( MultiTimers ), INTENT(inout) :: theTimers
   INTEGER, INTENT(in)                 :: tID

      CALL theTimers % PointToTimer( tID )
      CALL theTimers % current % StopTimer( )
      CALL theTimers % AccumulateTimings( )

 END SUBROUTINE StopThisTimer
!
!
!
 SUBROUTINE StartTimer( thetimer )
 !
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS( RoutineTimer ), INTENT(inout) :: theTimer

      theTimer % started = .TRUE.
      theTimer % stopped = .FALSE.
      
      CALL CPU_TIME( theTimer % startTime )

 END SUBROUTINE StartTimer
!
!
!
 SUBROUTINE StopTimer( thetimer )
 !
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS( RoutineTimer ), INTENT(inout) :: theTimer

      CALL CPU_TIME( theTimer % stopTime )
      theTimer % stopped = .TRUE.

 END SUBROUTINE StopTimer
!
!
!
 FUNCTION ElapsedTime( theTimer ) RESULT( eTime )
 !
 !
 ! =============================================================================================== !
   IMPLICIT NONE
   CLASS( RoutineTimer ) :: theTimer
   REAL(prec)            :: eTime

      IF( theTimer % stopped )THEN
         eTime = theTimer % stopTime - theTimer % startTime 
      ELSE
         PRINT*, 'Module Timing.f90 : S/R ElapsedTime : Warning! Timer "', TRIM(theTimer % whatYourTiming),'" is not stopped'
         eTime = ZERO 
      ENDIF
 END FUNCTION ElapsedTime
!
!
!
 SUBROUTINE AccumulateTimings( theTimers, tID )
 ! S/R AccumulateTimings
 !
 ! Takes IN theTimers and adds the current measurement for
 ! theTimers % theFTTimers(k) % elapsedTime to  "totTimes(k)"
 !
 ! The attribute "nObs" is also incremented once
 !
 ! =============================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS( MultiTimers )          :: theTimers
   INTEGER, INTENT(IN), OPTIONAL :: tID
   ! LOCAL
   
    
      IF( PRESENT(tID) )THEN
      
         CALL theTimers % PointToTimer( tID ) ! Point the "current" pointer to the RoutineTimer with 
                                              ! timerID = tID.
      ENDIF
  
  
  
      theTimers % current % accumulatedTime = theTimers % current % accumulatedTime + &
                                              theTimers % current % ElapsedTime( )

      theTimers % current % nObs = theTimers % current % nObs + ONE 


 END SUBROUTINE AccumulateTimings
!
!
!
 SUBROUTINE Write_MultiTimers( theTimers )
 !
 !
 !
 ! ============================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS( MultiTimers ), INTENT(INOUT) :: theTimers 
   ! LOCAL
   INTEGER :: k, fUnit
   CHARACTER(defaultNameLength) :: tName

      OPEN( UNIT = NewUnit(fUnit), FILE = 'Timing.stats' ) 

      WRITE(fUnit,*), '====================== Timing Results ======================'
      WRITE(fUnit,*), ' '
      
      theTimers % current => theTimers % head
      k = 0
      DO WHILE( ASSOCIATED(theTimers % current) )
         k = k+1

         CALL theTimers % GetName( tName )

         WRITE(fUnit,*), tName
         WRITE(fUnit,*), 'Number of Measurements : ', theTimers % current % nObs
         WRITE(fUnit,*), 'Accumulated Time       : ', theTimers % current % accumulatedTime
         WRITE(fUnit,*), 'Average Time           : ', theTimers % current % accumulatedTime/theTimers % current % nObs
         WRITE(fUnit,*), '------------------------------------------------------------'

         CALL theTimers % MoveToNext( )

      ENDDO
      CLOSE(fUnit)


      theTimers % current => theTimers % head
      
      PRINT*, '====================== Timing Results ======================'
      PRINT*, ' '
      
      DO WHILE( ASSOCIATED(theTimers % current) )

         CALL theTimers % GetName( tName ) 

         PRINT*, tName
         PRINT*, 'Number of Measurements : ', theTimers % current % nObs
         PRINT*, 'Accumulated Time       : ', theTimers % current % accumulatedTime
         PRINT*, 'Average Time           : ', theTimers % current % accumulatedTime/theTimers % current % nObs
         PRINT*, '------------------------------------------------------------'

         CALL theTimers % MoveToNext( )

      ENDDO

      PRINT*, '------------------------------------------------------------'
      PRINT*,' Timing Results saved to  "Timing.Stats" '

 END SUBROUTINE Write_MultiTimers

END MODULE Timing
