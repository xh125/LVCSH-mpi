!
! Copyright (C) 2002-2020 Quantum ESPRESSO Foundation
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!-----------------------------------------------------------------------
  FUNCTION trimcheck ( directory ,llinux)
    !-----------------------------------------------------------------------
    !! Verify if directory ends with / or \ for Windows, add one if needed.  
    !! Trim white spaces and put the result in trimcheck.
    !
    IMPLICIT NONE
    !
    CHARACTER (LEN=*), INTENT(IN) :: directory
    logical, intent(in) :: llinux
    CHARACTER (LEN=256) :: trimcheck
    character(len=1) :: separator
    INTEGER  :: l, i
    !
    
    if(llinux) then
      separator = '/'
    else
      separator = '\'
    endif
    
    l = LEN_TRIM( ADJUSTL(directory) )
    IF ( l == 0 )  CALL errore( 'trimcheck', ' input name empty', 1)
    IF ( l > LEN( trimcheck ) ) &
         CALL errore( 'trimcheck', ' input name too long', 1)
    !
    trimcheck = TRIM ( ADJUSTL(directory) )
    !
    ! for Windows: convert / to \
    !
    IF ( separator /= '/') THEN
       DO i = 1, l
          IF ( trimcheck(i:i) == '/' )  trimcheck(i:i) = separator
       END DO
    END IF
    !
    ! add final / or \ if not present - makes easier to add a file name
    !
    IF ( directory(l:l) /= separator ) THEN
       IF ( l < LEN( trimcheck ) ) THEN
          trimcheck(l+1:l+1) = separator
       ELSE
          CALL errore(  'trimcheck', ' input name too long', 2 )
       END IF
    END IF
    !
    RETURN
    !
  END FUNCTION trimcheck
  !
