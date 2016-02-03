! TestKeyRing.f90 ( *new with v2.1 - 2 February, 2016)
! 
! ====================================== LICENSE ================================================= !
!
!  This program is free software; you can redistribute it and/or modify
!  it under the terms of the GNU General Public License as published by
!  the Free Software Foundation; either version 2 of the License, or
!  (at your option) any later version.
!  
!  This program is distributed in the hope that it will be useful,
!  but WITHOUT ANY WARRANTY; without even the implied warranty of
!  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!  GNU General Public License for more details.
!  
!  You should have received a copy of the GNU General Public License
!  along with this program; if not, write to the Free Software
!  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
!  MA 02110-1301, USA.
! 
! ==================================== Module History ============================================ ! 
!
! o (v2.1) 2 February, 2016 
! o (newest version) module completion date (for this version)
!
! ========================================= Logs ================================================= !
! 2016-02-02  Joseph Schoonover  <jschoonover@-----> : schoonover.numerics@gmail.com
!
! On Feb. 2, 2016, this program was tested with valgrind, from which the following output is given:
!
!==13079== 
!==13079== HEAP SUMMARY:
!==13079==     in use at exit: 0 bytes in 0 blocks
!==13079==   total heap usage: 170 allocs, 170 frees, 16,781 bytes allocated
!==13079== 
!==13079== All heap blocks were freed -- no leaks are possible
!==13079== 
!==13079== For counts of detected and suppressed errors, rerun with: -v
!==13079== ERROR SUMMARY: 0 errors from 0 contexts (suppressed: 2 from 2)
! 
!   
! 
! //////////////////////////////////////////////////////////////////////////////////////////////// !


PROGRAM TestKeyRing

  USE KeyRingClass

  IMPLICIT NONE

  INTEGER, PARAMETER :: nNotches = 4
  INTEGER, PARAMETER :: keyTest = 25

  TYPE( KeyRing ) :: thisList
  
  INTEGER :: i, j, notches(1:nNotches)

  LOGICAL :: dataRetrieved

     CALL thisList % Build( )


     do i = 1, 50
 
        notches(1) = 1 + (i-1)*nNotches
        do j = 2, nNotches
           notches(j) = notches(j-1) + 1
        enddo

        CALL thisList % AddToList( i, Notches, nNotches )
   
     enddo

     CALL thisList % PrintToScreen( )
     
     notches(1) = 1 + (keyTest-1)*nNotches
     do j = 2, nNotches
        notches(j) = notches(j-1) + 1
     enddo

     CALL thisList % FindDataForNotches( notches, nNotches, i, dataRetrieved )
 
     PRINT*, ' keyTest : ', keyTest  
     PRINT*, ' nNotches :', nNotches
     PRINT*, ' notches tested : ', notches 
     PRINT*, ' Data Was Retrieved : ', dataRetrieved 
     PRINT*, ' Key Found : ', i
     
     CALL thisList % Trash( )

 
END PROGRAM TestKeyRing
