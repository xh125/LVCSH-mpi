! Copyright (C) 2002-2008 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!  Wrapper for MPI implementations that have problems with large messages
!

!  In some MPI implementation the communication subsystem
!  crashes when message exceeds a given size, so we need
!  to break down MPI communications in smaller pieces  
#define __MSGSIZ_MAX 100000
#define __BCAST_MSGSIZ_MAX 100000
!  Some implementation of MPI (OpenMPI) if it is not well tuned for the given
!  network hardware (InfiniBand) tend to lose performance or get stuck inside
!  collective routines if processors are not well synchronized
!  A barrier fixes the problem
!
!#define __USE_BARRIER  
#define __MPI  

!=----------------------------------------------------------------------------=!
!
! These routines allocate buffer spaces used in reduce_base_real_gpu.
! These should be in data_buffer.f90 but need to be here becouse size is
! depends on the __MSGSIZ_MAX definition

!module mp_base

!contains
SUBROUTINE allocate_buffers
    USE data_buffer
    IMPLICIT NONE
    INTEGER, PARAMETER :: maxb = __MSGSIZ_MAX
    !
    IF (.NOT. ALLOCATED(mp_buff_r)) ALLOCATE(mp_buff_r(maxb))
    IF (.NOT. ALLOCATED(mp_buff_i)) ALLOCATE(mp_buff_i(maxb))
    !
END SUBROUTINE allocate_buffers

SUBROUTINE deallocate_buffers
    USE data_buffer
    IMPLICIT NONE
    !
    DEALLOCATE(mp_buff_r, mp_buff_i)
    !
END SUBROUTINE deallocate_buffers

!=----------------------------------------------------------------------------=!
!

SUBROUTINE mp_synchronize( gid )
   USE mpi  
   IMPLICIT NONE
   INTEGER, INTENT(IN) :: gid
#if defined __MPI && defined __USE_BARRIER
   INTEGER :: ierr
   CALL mpi_barrier( gid, ierr )
   IF( ierr /= 0 ) CALL errore( 'mp_synchronize ', ' error in mpi_barrier ', ierr )
#endif
   RETURN
END SUBROUTINE mp_synchronize


  
  SUBROUTINE bcast_real( array, n, root, gid ) 
    use kinds,only : dp,dpc
    use mpi    
    !
    ! Broadcast Real
    !
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: n, root, gid
    REAL(DP) :: array( n )
#if defined __MPI
    INTEGER :: msgsiz_max = __BCAST_MSGSIZ_MAX    
    INTEGER :: nblk, blksiz, iblk, istart, ierr


#if defined __USE_BARRIER
    CALL mp_synchronize( gid )
#endif

    IF( n <= msgsiz_max ) THEN
       CALL MPI_BCAST( array, n, MPI_DOUBLE_PRECISION, root, gid, ierr )
       IF( ierr /= 0 ) CALL errore( ' bcast_real ', ' error in mpi_bcast 1 ', ierr )
    ELSE
       nblk   = n / msgsiz_max
       blksiz = msgsiz_max
       DO iblk = 1, nblk
          istart = (iblk-1)*msgsiz_max + 1
          CALL MPI_BCAST( array( istart ), blksiz, MPI_DOUBLE_PRECISION, root, gid, ierr )
          IF( ierr /= 0 ) CALL errore( ' bcast_real ', ' error in mpi_bcast 2 ', ierr )
       END DO
       blksiz = MOD( n, msgsiz_max )
       IF( blksiz > 0 ) THEN
          istart = nblk * msgsiz_max + 1
          CALL MPI_BCAST( array( istart ), blksiz, MPI_DOUBLE_PRECISION, root, gid, ierr )
          IF( ierr /= 0 ) CALL errore( ' bcast_real ', ' error in mpi_bcast 3 ', ierr )
       END IF
    END IF
#endif
    RETURN
  END SUBROUTINE bcast_real 

  subroutine bcast_cmpl( array, n ,root, gid)
    !
    ! Broadcast Complex
    !
    use kinds,only : dp,dpc
    use mpi    
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: n, root, gid
    Complex(kind=dpc) :: array( n )
#if defined __MPI
    INTEGER :: msgsiz_max = __BCAST_MSGSIZ_MAX    
    INTEGER :: nblk, blksiz, iblk, istart, ierr
   ! 
#if defined __USE_BARRIER
   CALL mp_synchronize(gid)
#endif
   ! 
    IF( n <= msgsiz_max ) THEN
       CALL MPI_BCAST( array, n, MPI_COMPLEX, root, gid, ierr )
       IF( ierr /= 0 ) CALL errore( ' bcast_real ', ' error in mpi_bcast 1 ', ierr )
    ELSE
       nblk   = n / msgsiz_max
       blksiz = msgsiz_max
       DO iblk = 1, nblk
          istart = (iblk-1)*msgsiz_max + 1
          CALL MPI_BCAST( array( istart ), blksiz, MPI_COMPLEX, root, gid, ierr )
          IF( ierr /= 0 ) CALL errore( ' bcast_real ', ' error in mpi_bcast 2 ', ierr )
       END DO
       blksiz = MOD( n, msgsiz_max )
       IF( blksiz > 0 ) THEN
          istart = nblk * msgsiz_max + 1
          CALL MPI_BCAST( array( istart ), blksiz, MPI_COMPLEX, root, gid, ierr )
          IF( ierr /= 0 ) CALL errore( ' bcast_real ', ' error in mpi_bcast 3 ', ierr )
       END IF
    END IF
#endif
    RETURN  
  end subroutine bcast_cmpl

  !------------------------------------------------------------------------------!
  SUBROUTINE bcast_integer(array, n, root, gid)
  !------------------------------------------------------------------------------!
  !! 
  !! Broadcast integers
  !!   
    use kinds,only : dp,dpc
    use mpi
    IMPLICIT NONE
    INTEGER, INTENT(in) :: n, root, gid
    INTEGER :: array(n)
#if defined __MPI
    INTEGER :: msgsiz_max = __BCAST_MSGSIZ_MAX    
    INTEGER :: nblk, blksiz, iblk, istart, ierr , i   
#if defined __USE_BARRIER
   CALL mp_synchronize(gid)
#endif
   !     
    IF(n <= msgsiz_max) THEN
      CALL MPI_BCAST(array, n, MPI_INTEGER, root, gid, ierr)
      IF(ierr /= 0) CALL errore(' bcast_integer ', ' error in mpi_bcast 1 ', ierr)
    ELSE
      nblk   = n / msgsiz_max
      blksiz = msgsiz_max
      DO iblk = 1, nblk
        istart = (iblk - 1) * msgsiz_max + 1
        CALL MPI_BCAST(array(istart), blksiz, MPI_INTEGER, root, gid, ierr )
        IF(ierr /= 0) CALL errore(' bcast_integer ', ' error in mpi_bcast 2 ', ierr)
      ENDDO
      blksiz = MOD( n, msgsiz_max )
      IF(blksiz > 0) THEN
        istart = nblk * msgsiz_max + 1
        CALL MPI_BCAST(array( istart ), blksiz, MPI_INTEGER, root, gid, ierr)
        IF(ierr /= 0) CALL errore(' bcast_integer ', ' error in mpi_bcast 3 ', ierr)
      ENDIF
    ENDIF
#endif    
    RETURN
  !------------------------------------------------------------------------------!
  END SUBROUTINE bcast_integer
  !------------------------------------------------------------------------------!
  ! 
  !------------------------------------------------------------------------------!
  SUBROUTINE bcast_integer8(array, n, root, gid)
  use kinds,only : dp,dpc,i8b
  use mpi
  !------------------------------------------------------------------------------!
  !! 
  !! Broadcast integers
  !!   
    ! 
    IMPLICIT NONE
    INTEGER, INTENT(in) :: n, root, gid
    INTEGER(KIND = i8b) :: array(n)
#if defined __MPI
    INTEGER :: msgsiz_max = __BCAST_MSGSIZ_MAX    
    INTEGER :: nblk, blksiz, iblk, istart, ierr    
    !
   ! 
#if defined __USE_BARRIER
   CALL mp_synchronize(gid)
#endif
   !     
    IF(n <= msgsiz_max) THEN
      CALL MPI_BCAST(array, n, MPI_INTEGER8, root, gid, ierr)
      IF(ierr /= 0) CALL errore(' bcast_integer8 ', ' error in mpi_bcast 1 ', ierr)
    ELSE
      nblk   = n / msgsiz_max
      blksiz = msgsiz_max
      DO iblk = 1, nblk
        istart = (iblk - 1) * msgsiz_max + 1
        CALL MPI_BCAST(array(istart), blksiz, MPI_INTEGER8, root, gid, ierr )
        IF(ierr /= 0) CALL errore(' bcast_integer8 ', ' error in mpi_bcast 2 ', ierr)
      ENDDO
      blksiz = MOD( n, msgsiz_max )
      IF(blksiz > 0) THEN
        istart = nblk * msgsiz_max + 1
        CALL MPI_BCAST(array( istart ), blksiz, MPI_INTEGER8, root, gid, ierr)
        IF(ierr /= 0) CALL errore(' bcast_integer8 ', ' error in mpi_bcast 3 ', ierr)
      ENDIF
    ENDIF
#endif    
    RETURN
  !------------------------------------------------------------------------------!
  END SUBROUTINE bcast_integer8
  !------------------------------------------------------------------------------!
  ! 
  !------------------------------------------------------------------------------!
  SUBROUTINE bcast_logical( array, n, root, gid )
    use kinds,only : dp,dpc
    use mpi
    use global_mpi,only : iproc
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: n, root, gid
    LOGICAL :: array( n )
#if defined __MPI
    INTEGER :: msgsiz_max = __BCAST_MSGSIZ_MAX    
    INTEGER :: nblk, blksiz, iblk, istart, ierr
#if defined __USE_BARRIER
   CALL mp_synchronize(gid)
#endif     
    IF( n <= msgsiz_max ) THEN
      CALL MPI_BCAST( array, n, MPI_LOGICAL, root, gid, ierr )
      IF( ierr /= 0 ) CALL errore( ' bcast_logical ', ' error in mpi_bcast 1 ', ierr )
    ELSE
       nblk   = n / msgsiz_max
       blksiz = msgsiz_max
       DO iblk = 1, nblk
          istart = (iblk-1)*msgsiz_max + 1
          CALL MPI_BCAST( array( istart ), blksiz, MPI_LOGICAL, root, gid, ierr )
          IF( ierr /= 0 ) CALL errore( ' bcast_logical ', ' error in mpi_bcast 2 ', ierr )
       END DO
       blksiz = MOD( n, msgsiz_max )
       IF( blksiz > 0 ) THEN
          istart = nblk * msgsiz_max + 1
          CALL MPI_BCAST( array( istart ), blksiz, MPI_LOGICAL, root, gid, ierr )
          IF( ierr /= 0 ) CALL errore( ' bcast_logical ', ' error in mpi_bcast 3 ', ierr )
       END IF
    END IF
#endif
    continue
    RETURN
  END SUBROUTINE bcast_logical

!end module mp_base
