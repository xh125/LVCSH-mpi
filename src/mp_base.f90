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
#define __USE_BARRIER
#define __USE_INPLACE_MPI  
#define __MPI  
!#define __TRACE
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


!=----------------------------------------------------------------------------=!
!  
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
       CALL MPI_BCAST( array, n, MPI_DOUBLE_COMPLEX, root, gid, ierr )
       IF( ierr /= 0 ) CALL errore( ' bcast_real ', ' error in mpi_bcast 1 ', ierr )
    ELSE
       nblk   = n / msgsiz_max
       blksiz = msgsiz_max
       DO iblk = 1, nblk
          istart = (iblk-1)*msgsiz_max + 1
          CALL MPI_BCAST( array( istart ), blksiz, MPI_DOUBLE_COMPLEX, root, gid, ierr )
          IF( ierr /= 0 ) CALL errore( ' bcast_real ', ' error in mpi_bcast 2 ', ierr )
       END DO
       blksiz = MOD( n, msgsiz_max )
       IF( blksiz > 0 ) THEN
          istart = nblk * msgsiz_max + 1
          CALL MPI_BCAST( array( istart ), blksiz, MPI_DOUBLE_COMPLEX, root, gid, ierr )
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
    use mpi
    use global_mpi,only : iproc
    IMPLICIT NONE
    INTEGER, INTENT(in) :: n, root, gid
    INTEGER :: array(n)
#if defined __MPI
    INTEGER :: msgsiz_max = __BCAST_MSGSIZ_MAX    
    INTEGER :: nblk, blksiz, iblk, istart, ierr , i   
#if defined __TRACE
   WRITE(*, *) 'BCAST_INTEGER IN'
#endif

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
#if defined __TRACE
   WRITE(*, *) 'BCAST_INTEGER OUT'
#endif
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


!
! ... "reduce"-like subroutines
!
#if defined (__USE_INPLACE_MPI)
!
!----------------------------------------------------------------------------
SUBROUTINE reduce_base_real( ldim, ps, comm, root )
  !----------------------------------------------------------------------------
  !
  ! ... sums a distributed variable ps(ldim) over the processors.
  ! ... This version uses a fixed-length buffer of appropriate (?) ldim
  !
  use kinds,only : dp
  use mpi  
  !
  IMPLICIT NONE
  !
  INTEGER,  INTENT(IN)    :: ldim     ! size of the array
  REAL(DP)                :: ps(ldim) ! array whose elements have to be reduced
  INTEGER,  INTENT(IN)    :: comm    ! communicator
  INTEGER,  INTENT(IN)    :: root    ! if root <  0 perform a reduction to all procs
                                     ! if root >= 0 perform a reduce only to root proc.
  !
#if defined (__MPI)  
  !
  INTEGER            :: info, myid
  !
#if defined __TRACE
  write(*,*) 'reduce_base_real IN'
#endif
  !
  IF ( ldim <= 0 ) GO TO 1  ! go to the end of the subroutine
  !
  ! ... synchronize processes
  !
#if defined __USE_BARRIER
  CALL mp_synchronize( comm )
#endif
  !
  IF( root >= 0 ) THEN
     CALL mpi_comm_rank( comm, myid, info )
     IF( info /= 0 ) CALL errore( 'reduce_base_real', 'error in mpi_comm_rank', info )
     IF (myid == root) THEN
        CALL MPI_REDUCE( MPI_IN_PLACE, ps, ldim, MPI_DOUBLE_PRECISION, MPI_SUM, root, comm, info )
        IF( info /= 0 ) CALL errore( 'reduce_base_real', 'error in mpi_reduce 1', info )
     ELSE
        CALL MPI_REDUCE( ps, ps, ldim, MPI_DOUBLE_PRECISION, MPI_SUM, root, comm, info )
        IF( info /= 0 ) CALL errore( 'reduce_base_real', 'error in mpi_reduce 1', info )
     ENDIF
  ELSE
     CALL MPI_ALLREDUCE( MPI_IN_PLACE, ps, ldim, MPI_DOUBLE_PRECISION, MPI_SUM, comm, info )
     IF( info /= 0 ) CALL errore( 'reduce_base_real', 'error in mpi_allreduce 1', info )
  END IF
  !
1 CONTINUE
  !
#if defined __TRACE
  write(*,*) 'reduce_base_real OUT'
#endif
  !
#endif
  !
  RETURN
  !
END SUBROUTINE reduce_base_real
!
#else
!
!----------------------------------------------------------------------------
SUBROUTINE reduce_base_real( ldim, ps, comm, root )
  !----------------------------------------------------------------------------
  !
  ! ... sums a distributed variable ps(ldim) over the processors.
  ! ... This version uses a fixed-length buffer of appropriate (?) ldim
  !
  use kinds,only : dp
  USE data_buffer, ONLY: buff => mp_buff_r
  use mpi  
  !
  IMPLICIT NONE
  !
  INTEGER,  INTENT(IN)    :: ldim     ! size of the array
  REAL(DP)                :: ps(ldim) ! array whose elements have to be reduced
  INTEGER,  INTENT(IN)    :: comm    ! communicator
  INTEGER,  INTENT(IN)    :: root    ! if root <  0 perform a reduction to all procs
                                     ! if root >= 0 perform a reduce only to root proc.
  !
#if defined (__MPI)  
  !
  INTEGER            :: info, n, nbuf, nproc, myid
  INTEGER, PARAMETER :: maxb = __MSGSIZ_MAX
  !
  !
#if defined __TRACE
  write(*,*) 'reduce_base_real IN'
#endif

  CALL mpi_comm_size( comm, nproc, info )
  IF( info /= 0 ) CALL errore( 'reduce_base_real', 'error in mpi_comm_size', info )

  CALL mpi_comm_rank( comm, myid, info )
  IF( info /= 0 ) CALL errore( 'reduce_base_real', 'error in mpi_comm_rank', info )
  !
  IF ( ldim <= 0 .OR. nproc <= 1 ) GO TO 1  ! go to the end of the subroutine
  !
  ! ... synchronize processes
  !
#if defined __USE_BARRIER
  CALL mp_synchronize( comm )
#endif
  !
  nbuf = ldim / maxb
  !
  DO n = 1, nbuf
     !
     IF( root >= 0 ) THEN
        CALL MPI_REDUCE( ps(1+(n-1)*maxb), buff, maxb, MPI_DOUBLE_PRECISION, MPI_SUM, root, comm, info )
        IF( info /= 0 ) CALL errore( 'reduce_base_real', 'error in mpi_reduce 1', info )
     ELSE
        CALL MPI_ALLREDUCE( ps(1+(n-1)*maxb), buff, maxb, MPI_DOUBLE_PRECISION, MPI_SUM, comm, info )
        IF( info /= 0 ) CALL errore( 'reduce_base_real', 'error in mpi_allreduce 1', info )
     END IF
     !                    
     IF( root < 0 ) THEN
        ps((1+(n-1)*maxb):(n*maxb)) = buff(1:maxb)
     ELSE IF( root == myid ) THEN
        ps((1+(n-1)*maxb):(n*maxb)) = buff(1:maxb)
     END IF
     !
  END DO
  !
  ! ... possible remaining elements < maxb
  !
  IF ( ( ldim - nbuf * maxb ) > 0 ) THEN
     !
     IF( root >= 0 ) THEN
        CALL MPI_REDUCE( ps(1+nbuf*maxb), buff, (ldim-nbuf*maxb), MPI_DOUBLE_PRECISION, MPI_SUM, root, comm, info )
        IF( info /= 0 ) CALL errore( 'reduce_base_real', 'error in mpi_reduce 2', info )
     ELSE
        CALL MPI_ALLREDUCE( ps(1+nbuf*maxb), buff, (ldim-nbuf*maxb), MPI_DOUBLE_PRECISION, MPI_SUM, comm, info )
        IF( info /= 0 ) CALL errore( 'reduce_base_real', 'error in mpi_allreduce 2', info )
     END IF
     !
     IF( root < 0 ) THEN
        ps((1+nbuf*maxb):ldim) = buff(1:(ldim-nbuf*maxb))
     ELSE IF( root == myid ) THEN
        ps((1+nbuf*maxb):ldim) = buff(1:(ldim-nbuf*maxb))
     END IF
     !
  END IF
  !
1 CONTINUE
  !
#if defined __TRACE
  write(*,*) 'reduce_base_real OUT'
#endif
  !
#endif
  !
  RETURN
  !
END SUBROUTINE reduce_base_real
!
#endif
!

!
#if defined (__USE_INPLACE_MPI)
!
!----------------------------------------------------------------------------
SUBROUTINE reduce_base_cmpl( ldim, ps, comm, root )
  !----------------------------------------------------------------------------
  !
  ! ... sums a distributed variable ps(ldim) over the processors.
  ! ... This version uses a fixed-length buffer of appropriate (?) ldim
  !
  use kinds,only : dp,dpc
  use mpi  
  !
  IMPLICIT NONE
  !
  INTEGER,  INTENT(IN)    :: ldim     ! size of the array
  complex(kind=dpc)       :: ps(ldim) ! array whose elements have to be reduced
  INTEGER,  INTENT(IN)    :: comm    ! communicator
  INTEGER,  INTENT(IN)    :: root    ! if root <  0 perform a reduction to all procs
                                     ! if root >= 0 perform a reduce only to root proc.
  !
#if defined (__MPI)  
  !
  INTEGER            :: info, myid
  !
#if defined __TRACE
  write(*,*) 'reduce_base_real IN'
#endif
  !
  IF ( ldim <= 0 ) GO TO 1  ! go to the end of the subroutine
  !
  ! ... synchronize processes
  !
#if defined __USE_BARRIER
  CALL mp_synchronize( comm )
#endif
  !
  IF( root >= 0 ) THEN
     CALL mpi_comm_rank( comm, myid, info )
     IF( info /= 0 ) CALL errore( 'reduce_base_real', 'error in mpi_comm_rank', info )
     IF (myid == root) THEN
        CALL MPI_REDUCE( MPI_IN_PLACE, ps, ldim, MPI_DOUBLE_complex, MPI_SUM, root, comm, info )
        IF( info /= 0 ) CALL errore( 'reduce_base_real', 'error in mpi_reduce 1', info )
     ELSE
        CALL MPI_REDUCE( ps, ps, ldim, MPI_DOUBLE_complex, MPI_SUM, root, comm, info )
        IF( info /= 0 ) CALL errore( 'reduce_base_real', 'error in mpi_reduce 1', info )
     ENDIF
  ELSE
     CALL MPI_ALLREDUCE( MPI_IN_PLACE, ps, ldim, MPI_DOUBLE_complex, MPI_SUM, comm, info )
     IF( info /= 0 ) CALL errore( 'reduce_base_real', 'error in mpi_allreduce 1', info )
  END IF
  !
1 CONTINUE
  !
#if defined __TRACE
  write(*,*) 'reduce_base_real OUT'
#endif
  !
#endif
  !
  RETURN
  !
END SUBROUTINE reduce_base_cmpl
!
#else
!
!----------------------------------------------------------------------------
SUBROUTINE reduce_base_cmpl( ldim, ps, comm, root )
  !----------------------------------------------------------------------------
  !
  ! ... sums a distributed variable ps(ldim) over the processors.
  ! ... This version uses a fixed-length buffer of appropriate (?) ldim
  !
  use kinds,only : dp,dpc
  USE data_buffer, ONLY: buff => mp_buff_c
  use mpi  
  !
  IMPLICIT NONE
  !
  INTEGER,  INTENT(IN)    :: ldim     ! size of the array
  complex(kind=dpc)       :: ps(ldim) ! array whose elements have to be reduced
  INTEGER,  INTENT(IN)    :: comm    ! communicator
  INTEGER,  INTENT(IN)    :: root    ! if root <  0 perform a reduction to all procs
                                     ! if root >= 0 perform a reduce only to root proc.
  !
#if defined (__MPI)  
  !
  INTEGER            :: info, n, nbuf, nproc, myid
  INTEGER, PARAMETER :: maxb = __MSGSIZ_MAX
  !
  !
#if defined __TRACE
  write(*,*) 'reduce_base_real IN'
#endif

  CALL mpi_comm_size( comm, nproc, info )
  IF( info /= 0 ) CALL errore( 'reduce_base_real', 'error in mpi_comm_size', info )

  CALL mpi_comm_rank( comm, myid, info )
  IF( info /= 0 ) CALL errore( 'reduce_base_real', 'error in mpi_comm_rank', info )
  !
  IF ( ldim <= 0 .OR. nproc <= 1 ) GO TO 1  ! go to the end of the subroutine
  !
  ! ... synchronize processes
  !
#if defined __USE_BARRIER
  CALL mp_synchronize( comm )
#endif
  !
  nbuf = ldim / maxb
  !
  DO n = 1, nbuf
     !
     IF( root >= 0 ) THEN
        CALL MPI_REDUCE( ps(1+(n-1)*maxb), buff, maxb, MPI_DOUBLE_complex, MPI_SUM, root, comm, info )
        IF( info /= 0 ) CALL errore( 'reduce_base_real', 'error in mpi_reduce 1', info )
     ELSE
        CALL MPI_ALLREDUCE( ps(1+(n-1)*maxb), buff, maxb, MPI_DOUBLE_complex, MPI_SUM, comm, info )
        IF( info /= 0 ) CALL errore( 'reduce_base_real', 'error in mpi_allreduce 1', info )
     END IF
     !                    
     IF( root < 0 ) THEN
        ps((1+(n-1)*maxb):(n*maxb)) = buff(1:maxb)
     ELSE IF( root == myid ) THEN
        ps((1+(n-1)*maxb):(n*maxb)) = buff(1:maxb)
     END IF
     !
  END DO
  !
  ! ... possible remaining elements < maxb
  !
  IF ( ( ldim - nbuf * maxb ) > 0 ) THEN
     !
     IF( root >= 0 ) THEN
        CALL MPI_REDUCE( ps(1+nbuf*maxb), buff, (ldim-nbuf*maxb), MPI_DOUBLE_complex, MPI_SUM, root, comm, info )
        IF( info /= 0 ) CALL errore( 'reduce_base_real', 'error in mpi_reduce 2', info )
     ELSE
        CALL MPI_ALLREDUCE( ps(1+nbuf*maxb), buff, (ldim-nbuf*maxb), MPI_DOUBLE_complex, MPI_SUM, comm, info )
        IF( info /= 0 ) CALL errore( 'reduce_base_real', 'error in mpi_allreduce 2', info )
     END IF
     !
     IF( root < 0 ) THEN
        ps((1+nbuf*maxb):ldim) = buff(1:(ldim-nbuf*maxb))
     ELSE IF( root == myid ) THEN
        ps((1+nbuf*maxb):ldim) = buff(1:(ldim-nbuf*maxb))
     END IF
     !
  END IF
  !
1 CONTINUE
  !
#if defined __TRACE
  write(*,*) 'reduce_base_real OUT'
#endif
  !
#endif
  !
  RETURN
  !
END SUBROUTINE reduce_base_cmpl
!
#endif
!


!
#if defined (__USE_INPLACE_MPI)
  !
  !----------------------------------------------------------------------------
  SUBROUTINE reduce_base_integer(ldim, ps, comm, root)
  !----------------------------------------------------------------------------
  !!
  !! Sums a distributed variable ps(ldim) over the processors.
  !! This version uses a fixed-length buffer of appropriate (?) ldim
  !!
  USE kinds,       ONLY : DP
  use mpi  
  !
  IMPLICIT NONE
  !
  INTEGER, INTENT(in)    :: ldim
  INTEGER, INTENT(inout) :: ps(ldim)
  INTEGER, INTENT(in)    :: comm    ! communicator
  INTEGER, INTENT(in)    :: root    ! if root <  0 perform a reduction to all procs
                                    ! if root >= 0 perform a reduce only to root proc.
  !
#if defined (__MPI)  
  !
  INTEGER            :: info, myid
  !
#if defined __TRACE
  WRITE(*, *) 'reduce_base_integer IN'
#endif
  !
  ! ... synchronize processes
  !
#if defined __USE_BARRIER
  CALL mp_synchronize(comm)
#endif
  !
  IF(root >= 0) THEN
    CALL mpi_comm_rank( comm, myid, info )
    IF( info /= 0 ) CALL errore( 'reduce_base_integer', 'error in mpi_comm_rank', info )
    IF (myid == root) THEN
      CALL MPI_REDUCE(MPI_IN_PLACE, ps, ldim, MPI_INTEGER, MPI_SUM, root, comm, info)
      IF(info /= 0) CALL errore('reduce_base_integer', 'error in mpi_reduce 1', info)
    ELSE
      CALL MPI_REDUCE(ps, ps, ldim, MPI_INTEGER, MPI_SUM, root, comm, info)
      IF(info /= 0) CALL errore('reduce_base_integer', 'error in mpi_reduce 1', info)
    ENDIF
  ELSE
    CALL MPI_ALLREDUCE(MPI_IN_PLACE, ps, ldim, MPI_INTEGER, MPI_SUM, comm, info)
    IF(info /= 0 ) CALL errore('reduce_base_integer', 'error in mpi_allreduce 1', info)
  ENDIF
  !
#if defined __TRACE
  WRITE(*, *) 'reduce_base_integer OUT'
#endif
  !
#endif
  !
  RETURN
  !----------------------------------------------------------------------------
  END SUBROUTINE reduce_base_integer
  !----------------------------------------------------------------------------
  !
  !----------------------------------------------------------------------------
  SUBROUTINE reduce_base_integer8(ldim, ps, comm, root)
  !----------------------------------------------------------------------------
  !!
  !! Sums a distributed variable ps(ldim) over the processors.
  !! This version uses a fixed-length buffer of appropriate (?) ldim
  !!
  USE kinds,       ONLY : DP, i8b
  use mpi
  !
  IMPLICIT NONE
  !
  INTEGER, INTENT(in)  :: ldim
  INTEGER(KIND = i8b), INTENT(inout) :: ps(ldim)
  INTEGER, INTENT(in)  :: comm    ! communicator
  INTEGER, INTENT(in)  :: root    ! if root <  0 perform a reduction to all procs
                                  ! if root >= 0 perform a reduce only to root proc.
  !   
#if defined (__MPI)  
  !
  INTEGER            :: info, myid
  !
#if defined __TRACE
  WRITE(*, *) 'reduce_base_integer IN'
#endif
  !
  ! ... synchronize processes
  !
#if defined __USE_BARRIER
  CALL mp_synchronize(comm)
#endif
  !
  IF(root >= 0) THEN
    CALL mpi_comm_rank( comm, myid, info )
    IF( info /= 0 ) CALL errore( 'reduce_base_integer8', 'error in mpi_comm_rank', info )
    IF (myid == root) THEN
      CALL MPI_REDUCE(MPI_IN_PLACE, ps, ldim, MPI_INTEGER8, MPI_SUM, root, comm, info)
      IF(info /= 0) CALL errore('reduce_base_integer8', 'error in mpi_reduce 1', info)
    ELSE
      CALL MPI_REDUCE(ps, ps, ldim, MPI_INTEGER8, MPI_SUM, root, comm, info)
      IF(info /= 0) CALL errore('reduce_base_integer8', 'error in mpi_reduce 1', info)
    ENDIF
  ELSE
    CALL MPI_ALLREDUCE(MPI_IN_PLACE, ps, ldim, MPI_INTEGER8, MPI_SUM, comm, info)
    IF(info /= 0 ) CALL errore('reduce_base_integer8', 'error in mpi_allreduce 1', info)
  ENDIF
  !
#if defined __TRACE
  WRITE(*, *) 'reduce_base_integer OUT'
#endif
  !
#endif
  !
  RETURN
  !----------------------------------------------------------------------------
  END SUBROUTINE reduce_base_integer8
  !----------------------------------------------------------------------------
#else
  !
  !----------------------------------------------------------------------------
  SUBROUTINE reduce_base_integer(ldim, ps, comm, root)
  !----------------------------------------------------------------------------
  !!
  !! Sums a distributed variable ps(ldim) over the processors.
  !! This version uses a fixed-length buffer of appropriate (?) ldim
  !!
  USE kinds,  ONLY : DP
  USE data_buffer, ONLY : buff => mp_buff_i
  use mpi  
  !
  IMPLICIT NONE
  !
  INTEGER, INTENT(in)    :: ldim
  INTEGER, INTENT(inout) :: ps(ldim)
  INTEGER, INTENT(in)    :: comm    ! communicator
  INTEGER, INTENT(in)    :: root    ! if root <  0 perform a reduction to all procs
                                    ! if root >= 0 perform a reduce only to root proc.
  !
#if defined (__MPI)  
  !
  INTEGER            :: info, n, nbuf, nproc, myid
  INTEGER, PARAMETER :: maxb = __MSGSIZ_MAX
  !
#if defined __TRACE
  WRITE(*, *) 'reduce_base_integer IN'
#endif
  !
  CALL mpi_comm_size(comm, nproc, info)
  IF(info /= 0) CALL errore('reduce_base_integer', 'error in mpi_comm_size', info)
  ! 
  CALL mpi_comm_rank(comm, myid, info)
  IF(info /= 0) CALL errore('reduce_base_integer', 'error in mpi_comm_rank', info)
  !
  IF (ldim <= 0 .OR. nproc <= 1) GOTO 1  ! go to the end of the subroutine
  !
  ! ... synchronize processes
  !
#if defined __USE_BARRIER
  CALL mp_synchronize(comm)
#endif
  !
  nbuf = ldim / maxb
  !
  DO n = 1, nbuf
    !
    IF (root >= 0) THEN
      CALL MPI_REDUCE(ps(1 + (n - 1) * maxb), buff, maxb, MPI_INTEGER, MPI_SUM, root, comm, info)
      IF (info /= 0 ) CALL errore('reduce_base_integer', 'error in mpi_reduce 1', info)
    ELSE
      CALL MPI_ALLREDUCE(ps(1 + (n - 1) * maxb), buff, maxb, MPI_INTEGER, MPI_SUM, comm, info)
      IF (info /= 0) CALL errore('reduce_base_integer', 'error in mpi_allreduce 1', info)
    ENDIF
    !
    IF (root < 0 ) THEN
      ps((1 + (n - 1) * maxb):(n * maxb)) = buff(1:maxb)
    ELSEIF( root == myid) THEN
      ps((1 + (n - 1) * maxb):(n * maxb)) = buff(1:maxb)
    ENDIF
    !
  ENDDO
  !
  ! ... possible remaining elements < maxb
  !
  IF ((ldim - nbuf * maxb) > 0) THEN
    !
    IF (root >= 0) THEN
      CALL MPI_REDUCE(ps(1 + nbuf * maxb), buff, (ldim - nbuf * maxb), MPI_INTEGER, MPI_SUM, root, comm, info)
      IF(info /= 0) CALL errore('reduce_base_integer', 'error in mpi_reduce 2', info)
    ELSE
      CALL MPI_ALLREDUCE(ps(1 + nbuf * maxb), buff, (ldim - nbuf * maxb), MPI_INTEGER, MPI_SUM, comm, info)
      IF(info /= 0) CALL errore('reduce_base_integer', 'error in mpi_allreduce 2', info)
    ENDIF
    !
    IF(root < 0) THEN
      ps((1 + nbuf * maxb):ldim) = buff(1:(ldim - nbuf * maxb))
    ELSEIF(root == myid) THEN
      ps((1 + nbuf * maxb):ldim) = buff(1:(ldim - nbuf * maxb))
    ENDIF
    !
  ENDIF
  !
1 CONTINUE
  !
#if defined __TRACE
  WRITE(*, *) 'reduce_base_integer OUT'
#endif
  !
#endif
  !
  RETURN
  !----------------------------------------------------------------------------
  END SUBROUTINE reduce_base_integer
  !----------------------------------------------------------------------------
  !
  !----------------------------------------------------------------------------
  SUBROUTINE reduce_base_integer8(ldim, ps, comm, root)
  !----------------------------------------------------------------------------
  !!
  !! Sums a distributed variable ps(ldim) over the processors.
  !! This version uses a fixed-length buffer of appropriate (?) ldim
  !!
  USE kinds,  ONLY : DP, i8b
  USE data_buffer, ONLY : buff => mp_buff_i
  use mpi  
  !
  IMPLICIT NONE
  !
  INTEGER, INTENT(in)    :: ldim
  INTEGER(KIND = i8b), INTENT(inout) :: ps(ldim)
  INTEGER, INTENT(in)    :: comm    ! communicator
  INTEGER, INTENT(in)    :: root    ! if root <  0 perform a reduction to all procs
                                    ! if root >= 0 perform a reduce only to root proc.
  !
#if defined (__MPI)  
  !
  INTEGER            :: info, n, nbuf, nproc, myid
  INTEGER, PARAMETER :: maxb = __MSGSIZ_MAX
  !
#if defined __TRACE
  WRITE(*, *) 'reduce_base_integer8 IN'
#endif
  !
  CALL mpi_comm_size(comm, nproc, info)
  IF(info /= 0) CALL errore('reduce_base_integer8', 'error in mpi_comm_size', info)
  ! 
  CALL mpi_comm_rank(comm, myid, info)
  IF(info /= 0) CALL errore('reduce_base_integer8', 'error in mpi_comm_rank', info)
  !
  IF (ldim <= 0 .OR. nproc <= 1) GOTO 1  ! go to the end of the subroutine
  !
  ! ... synchronize processes
  !
#if defined __USE_BARRIER
  CALL mp_synchronize(comm)
#endif
  !
  nbuf = ldim / maxb
  !
  DO n = 1, nbuf
    !
    IF (root >= 0) THEN
      CALL MPI_REDUCE(ps(1 + (n - 1) * maxb), buff, maxb, MPI_INTEGER8, MPI_SUM, root, comm, info)
      IF (info /= 0 ) CALL errore('reduce_base_integer8', 'error in mpi_reduce 1', info)
    ELSE
      CALL MPI_ALLREDUCE(ps(1 + (n - 1) * maxb), buff, maxb, MPI_INTEGER8, MPI_SUM, comm, info)
      IF (info /= 0) CALL errore('reduce_base_integer8', 'error in mpi_allreduce 1', info)
    ENDIF
    !
    IF (root < 0 ) THEN
      ps((1 + (n - 1) * maxb):(n * maxb)) = buff(1:maxb)
    ELSEIF( root == myid) THEN
      ps((1 + (n - 1) * maxb):(n * maxb)) = buff(1:maxb)
    ENDIF
    !
  ENDDO
  !
  ! ... possible remaining elements < maxb
  !
  IF ((ldim - nbuf * maxb) > 0) THEN
    !
    IF (root >= 0) THEN
      CALL MPI_REDUCE(ps(1 + nbuf * maxb), buff, (ldim - nbuf * maxb), MPI_INTEGER8, MPI_SUM, root, comm, info)
      IF(info /= 0) CALL errore('reduce_base_integer', 'error in mpi_reduce 2', info)
    ELSE
      CALL MPI_ALLREDUCE(ps(1 + nbuf * maxb), buff, (ldim - nbuf * maxb), MPI_INTEGER8, MPI_SUM, comm, info)
      IF(info /= 0) CALL errore('reduce_base_integer', 'error in mpi_allreduce 2', info)
    ENDIF
    !
    IF(root < 0) THEN
      ps((1 + nbuf * maxb):ldim) = buff(1:(ldim - nbuf * maxb))
    ELSEIF(root == myid) THEN
      ps((1 + nbuf * maxb):ldim) = buff(1:(ldim - nbuf * maxb))
    ENDIF
    !
  ENDIF
  !
1 CONTINUE
  !
#if defined __TRACE
  WRITE(*, *) 'reduce_base_integer8 OUT'
#endif
  !
#endif
  !
  RETURN
  !----------------------------------------------------------------------------
  END SUBROUTINE reduce_base_integer8
  !----------------------------------------------------------------------------
  !
#endif
  
  !----------------------------------------------------------------------------
  !
  ! ... "reduce"-like subroutines
  !
  !----------------------------------------------------------------------------
  SUBROUTINE reduce_base_real_to( ldim, ps, psout, comm, root )
  !----------------------------------------------------------------------------
  !
  ! ... sums a distributed variable ps(ldim) over the processors,
  ! ... and store the results in variable psout.
  ! ... This version uses a fixed-length buffer of appropriate (?) length
  !
  USE kinds, ONLY : DP
  USE mpi  
  !
  IMPLICIT NONE
  !
  INTEGER,  INTENT(IN)  :: ldim
  REAL(DP), INTENT(IN)  :: ps(ldim)
  REAL(DP)              :: psout(ldim)
  INTEGER,  INTENT(IN)  :: comm    ! communecator
  INTEGER,  INTENT(IN)  :: root    ! if root <  0 perform a reduction to all procs
                                   ! if root >= 0 perform a reduce only to root proc.
  !
#if defined (__MPI)  
  !
  INTEGER            :: info, n, nbuf, nproc, myid
  INTEGER, PARAMETER :: maxb = __MSGSIZ_MAX
  !
#if defined __TRACE
  write(*,*) 'reduce_base_real_to IN'
#endif

  CALL mpi_comm_size( comm, nproc, info )
  IF( info /= 0 ) CALL errore( 'reduce_base_real_to', 'error in mpi_comm_size', info )

  CALL mpi_comm_rank( comm, myid, info )
  IF( info /= 0 ) CALL errore( 'reduce_base_real_to', 'error in mpi_comm_rank', info )
  !
  IF ( ldim > 0 .AND. nproc <= 1 ) THEN
     psout = ps
  END IF
  IF( ldim <= 0 .OR. nproc <= 1 ) GO TO 1 ! go to the end of the subroutine
  !
  ! ... synchronize processes
  !
#if defined __USE_BARRIER
  CALL mp_synchronize( comm )
#endif
  !
  nbuf = ldim / maxb
  !
  DO n = 1, nbuf
     !
     IF( root >= 0 ) THEN
        CALL MPI_REDUCE( ps(1+(n-1)*maxb), psout(1+(n-1)*maxb), maxb, MPI_DOUBLE_PRECISION, MPI_SUM, root, comm, info )
        IF( info /= 0 ) CALL errore( 'reduce_base_real_to', 'error in mpi_reduce 1', info )
     ELSE
        CALL MPI_ALLREDUCE( ps(1+(n-1)*maxb), psout(1+(n-1)*maxb), maxb, MPI_DOUBLE_PRECISION, MPI_SUM, comm, info )
        IF( info /= 0 ) CALL errore( 'reduce_base_real_to', 'error in mpi_allreduce 1', info )
     END IF
     !                    
  END DO
  !
  ! ... possible remaining elements < maxb
  !
  IF ( ( ldim - nbuf * maxb ) > 0 ) THEN
     !
     IF( root >= 0 ) THEN
        CALL MPI_REDUCE( ps(1+nbuf*maxb), psout(1+nbuf*maxb), (ldim-nbuf*maxb), MPI_DOUBLE_PRECISION, MPI_SUM, root, comm, info )
        IF( info /= 0 ) CALL errore( 'reduce_base_real_to', 'error in mpi_reduce 2', info )
     ELSE
        CALL MPI_ALLREDUCE( ps(1+nbuf*maxb), psout(1+nbuf*maxb), (ldim-nbuf*maxb), MPI_DOUBLE_PRECISION, MPI_SUM, comm, info )
        IF( info /= 0 ) CALL errore( 'reduce_base_real_to', 'error in mpi_allreduce 2', info )
     END IF
     !
  END IF
  !
1 CONTINUE
  !
#if defined __TRACE
  write(*,*) 'reduce_base_real_to OUT'
#endif
  !
#endif
  !
  RETURN
  !
END SUBROUTINE reduce_base_real_to

  !----------------------------------------------------------------------------
  SUBROUTINE reduce_base_cmpl_to( ldim, ps, psout, comm, root )
  !----------------------------------------------------------------------------
  !
  ! ... sums a distributed variable ps(ldim) over the processors,
  ! ... and store the results in variable psout.
  ! ... This version uses a fixed-length buffer of appropriate (?) length
  !
  USE kinds, ONLY : DP,dpc
  USE mpi  
  !
  IMPLICIT NONE
  !
  INTEGER,  INTENT(IN)  :: ldim
  complex(kind=dpc), INTENT(IN)  :: ps(ldim)
  complex(kind=dpc)              :: psout(ldim)
  INTEGER,  INTENT(IN)  :: comm    ! communecator
  INTEGER,  INTENT(IN)  :: root    ! if root <  0 perform a reduction to all procs
                                   ! if root >= 0 perform a reduce only to root proc.
  !
#if defined (__MPI)  
  !
  INTEGER            :: info, n, nbuf, nproc, myid
  INTEGER, PARAMETER :: maxb = __MSGSIZ_MAX
  !
#if defined __TRACE
  write(*,*) 'reduce_base_real_to IN'
#endif

  CALL mpi_comm_size( comm, nproc, info )
  IF( info /= 0 ) CALL errore( 'reduce_base_real_to', 'error in mpi_comm_size', info )

  CALL mpi_comm_rank( comm, myid, info )
  IF( info /= 0 ) CALL errore( 'reduce_base_real_to', 'error in mpi_comm_rank', info )
  !
  IF ( ldim > 0 .AND. nproc <= 1 ) THEN
     psout = ps
  END IF
  IF( ldim <= 0 .OR. nproc <= 1 ) GO TO 1 ! go to the end of the subroutine
  !
  ! ... synchronize processes
  !
#if defined __USE_BARRIER
  CALL mp_synchronize( comm )
#endif
  !
  nbuf = ldim / maxb
  !
  DO n = 1, nbuf
     !
     IF( root >= 0 ) THEN
        CALL MPI_REDUCE( ps(1+(n-1)*maxb), psout(1+(n-1)*maxb), maxb, MPI_DOUBLE_complex, MPI_SUM, root, comm, info )
        IF( info /= 0 ) CALL errore( 'reduce_base_real_to', 'error in mpi_reduce 1', info )
     ELSE
        CALL MPI_ALLREDUCE( ps(1+(n-1)*maxb), psout(1+(n-1)*maxb), maxb, MPI_DOUBLE_complex, MPI_SUM, comm, info )
        IF( info /= 0 ) CALL errore( 'reduce_base_real_to', 'error in mpi_allreduce 1', info )
     END IF
     !                    
  END DO
  !
  ! ... possible remaining elements < maxb
  !
  IF ( ( ldim - nbuf * maxb ) > 0 ) THEN
     !
     IF( root >= 0 ) THEN
        CALL MPI_REDUCE( ps(1+nbuf*maxb), psout(1+nbuf*maxb), (ldim-nbuf*maxb), MPI_DOUBLE_complex, MPI_SUM, root, comm, info )
        IF( info /= 0 ) CALL errore( 'reduce_base_real_to', 'error in mpi_reduce 2', info )
     ELSE
        CALL MPI_ALLREDUCE( ps(1+nbuf*maxb), psout(1+nbuf*maxb), (ldim-nbuf*maxb), MPI_DOUBLE_complex, MPI_SUM, comm, info )
        IF( info /= 0 ) CALL errore( 'reduce_base_real_to', 'error in mpi_allreduce 2', info )
     END IF
     !
  END IF
  !
1 CONTINUE
  !
#if defined __TRACE
  write(*,*) 'reduce_base_real_to OUT'
#endif
  !
#endif
  !
  RETURN
  !
END SUBROUTINE reduce_base_cmpl_to


 !----------------------------------------------------------------------------
SUBROUTINE reduce_base_integer_to( ldim, ps, psout, comm, root )
  !----------------------------------------------------------------------------
  !
  ! ... sums a distributed integer variable ps(ldim) over the processors, and
  ! ... saves the result on the output variable psout.
  ! ... This version uses a fixed-length buffer of appropriate (?) length
  !
  USE kinds, ONLY : DP
  USE mpi  
  !
  IMPLICIT NONE
  !
  INTEGER,  INTENT(IN)  :: ldim
  INTEGER,  INTENT(IN)  :: ps(ldim)
  INTEGER               :: psout(ldim)
  INTEGER,  INTENT(IN)  :: comm    ! communecator
  INTEGER,  INTENT(IN)  :: root    ! if root <  0 perform a reduction to all procs
                                     ! if root >= 0 perform a reduce only to root proc.
  !
#if defined (__MPI)  
  !
  INTEGER            :: info, n, nbuf, nproc, myid
  INTEGER, PARAMETER :: maxb = __MSGSIZ_MAX
  !
#if defined __TRACE
  write(*,*) 'reduce_base_integer_to IN'
#endif

  CALL mpi_comm_size( comm, nproc, info )
  IF( info /= 0 ) CALL errore( 'reduce_base_integer_to', 'error in mpi_comm_size', info )

  CALL mpi_comm_rank( comm, myid, info )
  IF( info /= 0 ) CALL errore( 'reduce_base_integer_to', 'error in mpi_comm_rank', info )
  !
  IF ( ldim > 0 .AND. nproc <= 1 ) THEN
     psout = ps
  END IF
  IF( ldim <= 0 .OR. nproc <= 1 ) GO TO 1 ! go to the end of the subroutine
  !
  ! ... synchronize processes
  !
#if defined __USE_BARRIER
  CALL mp_synchronize( comm )
#endif
  !
  nbuf = ldim / maxb
  !
  DO n = 1, nbuf
     !
     IF( root >= 0 ) THEN
        CALL MPI_REDUCE( ps(1+(n-1)*maxb), psout( 1+(n-1)*maxb ), maxb, MPI_INTEGER, MPI_SUM, root, comm, info )
        IF( info /= 0 ) CALL errore( 'reduce_base_integer_to', 'error in mpi_reduce 1', info )
     ELSE
        CALL MPI_ALLREDUCE( ps(1+(n-1)*maxb), psout( 1+(n-1)*maxb ), maxb, MPI_INTEGER, MPI_SUM, comm, info )
        IF( info /= 0 ) CALL errore( 'reduce_base_integer_to', 'error in mpi_allreduce 1', info )
     END IF
     !                    
  END DO
  !
  ! ... possible remaining elements < maxb
  !
  IF ( ( ldim - nbuf * maxb ) > 0 ) THEN
     !
     IF( root >= 0 ) THEN
        CALL MPI_REDUCE( ps(1+nbuf*maxb), psout(1+nbuf*maxb), (ldim-nbuf*maxb), MPI_INTEGER, MPI_SUM, root, comm, info )
        IF( info /= 0 ) CALL errore( 'reduce_base_integer_to', 'error in mpi_reduce 2', info )
     ELSE
        CALL MPI_ALLREDUCE( ps(1+nbuf*maxb), psout(1+nbuf*maxb), (ldim-nbuf*maxb), MPI_INTEGER, MPI_SUM, comm, info )
        IF( info /= 0 ) CALL errore( 'reduce_base_integer_to', 'error in mpi_allreduce 2', info )
     END IF
     !
  END IF
  !
1 CONTINUE
  !
#if defined __TRACE
  write(*,*) 'reduce_base_integer_to OUT'
#endif
  !
#endif
  !
  RETURN
  !
END SUBROUTINE reduce_base_integer_to 



