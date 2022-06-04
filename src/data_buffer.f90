!
! Copyright (C) Quantum ESPRESSO group
!
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
MODULE data_buffer
    USE kinds,  ONLY : DP

    !
    IMPLICIT NONE
    !
    REAL(kind=DP), ALLOCATABLE  :: mp_buff_r(:)
    INTEGER, ALLOCATABLE        :: mp_buff_i(:)
    PUBLIC :: mp_buff_r, mp_buff_i
    !

END MODULE data_buffer