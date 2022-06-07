#define __MPI
module mp
  use kinds,only : dp,dpc,i8b
  use io,only : stdout
  use mpi
  use global_mpi,only : iproc,ionode
  implicit none
  
  INTERFACE mp_bcast
    MODULE PROCEDURE mp_bcast_i1, mp_bcast_r1, mp_bcast_c1, &
      mp_bcast_z, mp_bcast_zv, &
      mp_bcast_iv, mp_bcast_rv, mp_bcast_cv, mp_bcast_l, mp_bcast_rm, &
      mp_bcast_cm, mp_bcast_im, mp_bcast_it, mp_bcast_rt, mp_bcast_lv, &
      mp_bcast_lm, mp_bcast_r4d, mp_bcast_r5d, mp_bcast_ct,  mp_bcast_c4d,&
      mp_bcast_c5d
  END INTERFACE

  INTERFACE mp_sum
    MODULE PROCEDURE mp_sum_i1, mp_sum_iv, mp_sum_im, mp_sum_it, &
      mp_sum_r1, mp_sum_rv, mp_sum_rm, mp_sum_rt, mp_sum_r4d, &
      mp_sum_c1, mp_sum_cv, mp_sum_cm, mp_sum_ct, mp_sum_c4d, &
      mp_sum_c5d, mp_sum_c6d, mp_sum_rmm, mp_sum_cmm, mp_sum_r5d
  END INTERFACE

  INTERFACE mp_root_sum
    MODULE PROCEDURE mp_root_sum_rm, mp_root_sum_cm
  END INTERFACE

  public :: mp_bcast,mp_sum,mp_root_sum
 

  character(len=80), private :: err_msg = ' '

  contains
!------------------------------------------------------------------------------!
!  mp_bcast

  SUBROUTINE mp_bcast_i1(msg,source,gid)
    IMPLICIT NONE
    INTEGER :: msg
    INTEGER :: source
    INTEGER, OPTIONAL, INTENT(IN) :: gid
    INTEGER :: group
    INTEGER :: msglen
    
#if defined(__MPI)
    msglen = 1
    group = mpi_comm_world
    IF( PRESENT( gid ) ) group = gid
    CALL bcast_integer( msg, msglen, source, group )
#endif  
  END SUBROUTINE mp_bcast_i1
!
!--------------------------------------------------------------------------!
  SUBROUTINE mp_bcast_iv(msg,source,gid)
    IMPLICIT NONE
    INTEGER :: msg(:)
    INTEGER :: source
    INTEGER, OPTIONAL, INTENT(IN) :: gid
    INTEGER :: group
    INTEGER :: msglen

#if defined(__MPI)  
    msglen = size(msg)
    group = mpi_comm_world
    IF( PRESENT( gid ) ) group = gid
    CALL bcast_integer( msg, msglen, source, group )
#endif   
  END SUBROUTINE mp_bcast_iv
!
!--------------------------------------------------------------------------!
  SUBROUTINE mp_bcast_im( msg, source, gid )
    IMPLICIT NONE
    INTEGER :: msg(:,:)
    INTEGER :: source
    INTEGER, OPTIONAL, INTENT(IN) :: gid
    INTEGER :: group
    INTEGER :: msglen

#if defined(__MPI)   
    msglen = size(msg)
    group = mpi_comm_world
    IF( PRESENT( gid ) ) group = gid
    CALL bcast_integer( msg, msglen, source, group )
#endif  
  END SUBROUTINE mp_bcast_im
!
!--------------------------------------------------------------------------!
!
  SUBROUTINE mp_bcast_it(msg,source,gid)
    IMPLICIT NONE
    INTEGER :: msg(:,:,:)
    INTEGER :: source
    INTEGER, OPTIONAL, INTENT(IN) :: gid
    INTEGER :: group
    INTEGER :: msglen

#if defined(__MPI)   
    msglen = size(msg)
    group = mpi_comm_world
    IF( PRESENT( gid ) ) group = gid
    CALL bcast_integer( msg, msglen, source, group )
#endif  
  END SUBROUTINE mp_bcast_it
!
!--------------------------------------------------------------------------!
!
  SUBROUTINE mp_bcast_r1(msg,source,gid)
    IMPLICIT NONE
    REAL (DP) :: msg
    INTEGER :: msglen, source
    INTEGER, OPTIONAL, INTENT(IN) :: gid
    INTEGER :: group

#if defined(__MPI)   
    msglen = 1
    group = mpi_comm_world
    IF( PRESENT( gid ) ) group = gid
    CALL bcast_real( msg, msglen, source, group )
#endif  
  END SUBROUTINE mp_bcast_r1
!
!--------------------------------------------------------------------------!
!
  SUBROUTINE mp_bcast_rv(msg,source,gid)
    IMPLICIT NONE
    REAL (DP) :: msg(:)
    INTEGER :: source
    INTEGER, OPTIONAL, INTENT(IN) :: gid
    INTEGER :: group
    INTEGER :: msglen

#if defined(__MPI)   
    msglen = size(msg)
    group = mpi_comm_world
    IF( PRESENT( gid ) ) group = gid
    CALL bcast_real( msg, msglen, source, group )
#endif  
  END SUBROUTINE mp_bcast_rv
!
!--------------------------------------------------------------------------!
!
  SUBROUTINE mp_bcast_rm(msg,source,gid)
    IMPLICIT NONE
    REAL (DP) :: msg(:,:)
    INTEGER :: source
    INTEGER, OPTIONAL, INTENT(IN) :: gid
    INTEGER :: group
    INTEGER :: msglen

#if defined(__MPI)   
    msglen = size(msg)
    group = mpi_comm_world
    IF( PRESENT( gid ) ) group = gid
    CALL bcast_real( msg, msglen, source, group )
#endif  
  END SUBROUTINE mp_bcast_rm
!
!--------------------------------------------------------------------------!
!
!   
!
  SUBROUTINE mp_bcast_rt(msg,source,gid)
    IMPLICIT NONE
    REAL (DP) :: msg(:,:,:)
    INTEGER :: source
    INTEGER, OPTIONAL, INTENT(IN) :: gid
    INTEGER :: group
    INTEGER :: msglen

#if defined(__MPI)   
    msglen = size(msg)
    group = mpi_comm_world
    IF( PRESENT( gid ) ) group = gid
    CALL bcast_real( msg, msglen, source, group )
#endif  
  END SUBROUTINE mp_bcast_rt
!
!--------------------------------------------------------------------------!
!
!   
!
  SUBROUTINE mp_bcast_r4d(msg, source, gid)
    IMPLICIT NONE
    REAL (DP) :: msg(:,:,:,:)
    INTEGER :: source
    INTEGER, OPTIONAL, INTENT(IN) :: gid
    INTEGER :: group
    INTEGER :: msglen

#if defined(__MPI)   
    msglen = size(msg)
    group = mpi_comm_world
    IF( PRESENT( gid ) ) group = gid
    CALL bcast_real( msg, msglen, source, group )
#endif  
  END SUBROUTINE mp_bcast_r4d

!
!--------------------------------------------------------------------------!  
!
  SUBROUTINE mp_bcast_r5d(msg, source, gid)
    IMPLICIT NONE
    REAL (DP) :: msg(:,:,:,:,:)
    INTEGER :: source
    INTEGER, OPTIONAL, INTENT(IN) :: gid
    INTEGER :: group
    INTEGER :: msglen

#if defined(__MPI)   
    msglen = size(msg)
    group = mpi_comm_world
    IF( PRESENT( gid ) ) group = gid
    CALL bcast_real( msg, msglen, source, group )
#endif  
  END SUBROUTINE mp_bcast_r5d

!--------------------------------------------------------------------------!
!
  SUBROUTINE mp_bcast_c1(msg,source,gid)
    IMPLICIT NONE
    COMPLEX (kind=dpc) :: msg
    INTEGER :: source
    INTEGER, OPTIONAL, INTENT(IN) :: gid
    INTEGER :: group
    INTEGER :: msglen

#if defined(__MPI)   
    msglen = 1
    group = mpi_comm_world
    IF( PRESENT( gid ) ) group = gid
    CALL bcast_cmpl( msg, msglen, source, group )
#endif  
  END SUBROUTINE mp_bcast_c1
!
!--------------------------------------------------------------------------!
  SUBROUTINE mp_bcast_cv(msg,source,gid)
    IMPLICIT NONE
    COMPLEX (kind=dpc) :: msg(:)
    INTEGER :: source
    INTEGER, OPTIONAL, INTENT(IN) :: gid
    INTEGER :: group
    INTEGER :: msglen

#if defined(__MPI)   
    msglen = size(msg)
    group = mpi_comm_world
    IF( PRESENT( gid ) ) group = gid
    CALL bcast_cmpl( msg, msglen, source, group )
#endif  
  END SUBROUTINE mp_bcast_cv
!
!--------------------------------------------------------------------------!
  SUBROUTINE mp_bcast_cm(msg,source,gid)
    IMPLICIT NONE
    COMPLEX (kind=dpc) :: msg(:,:)
    INTEGER :: source
    INTEGER, OPTIONAL, INTENT(IN) :: gid
    INTEGER :: group
    INTEGER :: msglen

#if defined(__MPI)   
    msglen = size(msg)
    group = mpi_comm_world
    IF( PRESENT( gid ) ) group = gid
    CALL bcast_cmpl( msg, msglen, source, group )
#endif  
  END SUBROUTINE mp_bcast_cm
!
!--------------------------------------------------------------------------!
  SUBROUTINE mp_bcast_ct(msg,source,gid)
    IMPLICIT NONE
    COMPLEX (kind=dpc) :: msg(:,:,:)
    INTEGER :: source
    INTEGER, OPTIONAL, INTENT(IN) :: gid
    INTEGER :: group
    INTEGER :: msglen

#if defined(__MPI)   
    msglen = size(msg)
    group = mpi_comm_world
    IF( PRESENT( gid ) ) group = gid
    CALL bcast_cmpl( msg, msglen, source, group )
#endif  
  END SUBROUTINE mp_bcast_ct

!
!--------------------------------------------------------------------------!
  SUBROUTINE mp_bcast_c4d(msg,source,gid)
    IMPLICIT NONE
    COMPLEX (kind=dpc) :: msg(:,:,:,:)
    INTEGER :: source
    INTEGER, OPTIONAL, INTENT(IN) :: gid
    INTEGER :: group
    INTEGER :: msglen

#if defined(__MPI)   
    msglen = size(msg)
    group = mpi_comm_world
    IF( PRESENT( gid ) ) group = gid
    CALL bcast_cmpl( msg, msglen, source, group )
#endif  
  END SUBROUTINE mp_bcast_c4d

  SUBROUTINE mp_bcast_c5d(msg,source,gid)
    IMPLICIT NONE
    COMPLEX (kind=dpc) :: msg(:,:,:,:,:)
    INTEGER :: source
    INTEGER, OPTIONAL, INTENT(IN) :: gid
    INTEGER :: group
    INTEGER :: msglen

#if defined(__MPI)   
    msglen = size(msg)
    group = mpi_comm_world
    IF( PRESENT( gid ) ) group = gid
    CALL bcast_cmpl( msg, msglen, source, group )
#endif  
  END SUBROUTINE mp_bcast_c5d

!
!--------------------------------------------------------------------------!

  SUBROUTINE mp_bcast_l(msg,source,gid)
    IMPLICIT NONE
    LOGICAL :: msg
    INTEGER :: source
    INTEGER, OPTIONAL, INTENT(IN) :: gid
    INTEGER :: group
    INTEGER :: msglen

#if defined(__MPI)   
    msglen = 1
    group = mpi_comm_world
    IF( PRESENT( gid ) ) group = gid
    CALL bcast_logical( msg, msglen, source, group )
#endif    
  END SUBROUTINE mp_bcast_l
!
!--------------------------------------------------------------------------!
!
!
  SUBROUTINE mp_bcast_lv(msg,source,gid)
    IMPLICIT NONE
    LOGICAL :: msg(:)
    INTEGER :: source
    INTEGER, OPTIONAL, INTENT(IN) :: gid
    INTEGER :: group
    INTEGER :: msglen

#if defined(__MPI)   
    msglen = size(msg)
    group = mpi_comm_world
    IF( PRESENT( gid ) ) group = gid
    CALL bcast_logical( msg, msglen, source, group )
#endif  
  END SUBROUTINE mp_bcast_lv

!--------------------------------------------------------------------------!
!
!   
!
  SUBROUTINE mp_bcast_lm(msg,source,gid)
    IMPLICIT NONE
    LOGICAL :: msg(:,:)
    INTEGER :: source
    INTEGER, OPTIONAL, INTENT(IN) :: gid
    INTEGER :: group
    INTEGER :: msglen

#if defined(__MPI)   
    msglen = size(msg)
    group = mpi_comm_world
    IF( PRESENT( gid ) ) group = gid
    CALL bcast_logical( msg, msglen, source, group )
#endif  
  END SUBROUTINE mp_bcast_lm


!
!--------------------------------------------------------------------------!
!
  SUBROUTINE mp_bcast_z(msg,source,gid)
    IMPLICIT NONE
    CHARACTER (len=*) :: msg
    INTEGER,intent(in) :: source
    INTEGER, OPTIONAL, INTENT(IN) :: gid
    INTEGER :: group
    INTEGER :: msglen, ierr, i
    INTEGER, ALLOCATABLE :: imsg(:)

#if defined(__MPI)   
    ierr = 0
    msglen = len(msg)
    group = mpi_comm_world
    IF( PRESENT( gid ) ) group = gid
    IF (ierr/=0) CALL mp_stop( 8014 )
    ALLOCATE (imsg(1:msglen), STAT=ierr)
    IF (ierr/=0) CALL mp_stop( 8015 )
    DO i = 1, msglen
      imsg(i) = ichar(msg(i:i))
    END DO
    
    CALL bcast_integer( imsg, msglen, source, group )
    DO i = 1, msglen
      msg(i:i) = char(imsg(i))
    END DO
    
    DEALLOCATE (imsg, STAT=ierr)
    IF (ierr/=0) CALL mp_stop( 8016 )
#endif
  END SUBROUTINE mp_bcast_z
!
!--------------------------------------------------------------------------!
!
!--------------------------------------------------------------------------!
!
  SUBROUTINE mp_bcast_zv(msg,source,gid)
    IMPLICIT NONE
    CHARACTER (len=*) :: msg(:)
    INTEGER :: source
    INTEGER, OPTIONAL, INTENT(IN) :: gid
    INTEGER :: group
    INTEGER :: msglen, m1, m2, ierr, i, j
    INTEGER, ALLOCATABLE :: imsg(:,:)

#if defined(__MPI)   
    ierr = 0
    m1 = LEN(msg)
    m2 = SIZE(msg)
    msglen = LEN(msg)*SIZE(msg)
    group = mpi_comm_world
    IF( PRESENT( gid ) ) group = gid
    ALLOCATE (imsg(1:m1,1:m2), STAT=ierr)
    IF (ierr/=0) CALL mp_stop( 8017 )
    DO j = 1, m2
      DO i = 1, m1
        imsg(i,j) = ichar(msg(j)(i:i))
      END DO
    END DO
    CALL bcast_integer( imsg, msglen, source, group )
    
    DO j = 1, m2
      DO i = 1, m1
        msg(j)(i:i) = char(imsg(i,j))
      END DO
    END DO
    DEALLOCATE (imsg, STAT=ierr)
    IF (ierr/=0) CALL mp_stop( 8018 )
#endif  
  END SUBROUTINE mp_bcast_zv
!
! end mp_bcast  
 
!------------------------------------------------------------------------------!
!
!..mp_sum
      SUBROUTINE mp_sum_i1(msg,gid)
        IMPLICIT NONE
        INTEGER, INTENT (INOUT) :: msg
        INTEGER, INTENT(IN) :: gid
#if defined(__MPI)
        INTEGER :: msglen
        msglen = 1
        CALL reduce_base_integer( msglen, msg, gid, -1 )
#endif
      END SUBROUTINE mp_sum_i1
      !
      !------------------------------------------------------------------------------!
      SUBROUTINE mp_sum_iv(msg, gid)
      !------------------------------------------------------------------------------!
      !! 
      !! MPI sum an integer vector from all cores and bcast the result to all.  
      !! 
      IMPLICIT NONE
      ! 
      INTEGER, INTENT(inout) :: msg(:)
      INTEGER, INTENT(in) :: gid
#if defined(__MPI)
      INTEGER :: msglen
      msglen = SIZE(msg)
      CALL reduce_base_integer(msglen, msg, gid, -1)
#endif
      !------------------------------------------------------------------------------!
      END SUBROUTINE mp_sum_iv
      !------------------------------------------------------------------------------!
      ! 
      !------------------------------------------------------------------------------!
      SUBROUTINE mp_sum_i8v(msg, gid)
      !------------------------------------------------------------------------------!
      !! 
      !! MPI sum an integer vector from all cores and bcast the result to all.  
      !! 
      IMPLICIT NONE
      ! 
      INTEGER(KIND = i8b), INTENT(inout) :: msg(:)
      INTEGER, INTENT(in) :: gid
#if defined(__MPI)
      INTEGER :: msglen
      msglen = SIZE(msg)
      CALL reduce_base_integer8(msglen, msg, gid, -1)
#endif
      !------------------------------------------------------------------------------!
      END SUBROUTINE mp_sum_i8v
      !------------------------------------------------------------------------------!
      !
      !------------------------------------------------------------------------------!
      SUBROUTINE mp_sum_im(msg,gid)
      !------------------------------------------------------------------------------!
        IMPLICIT NONE
        INTEGER, INTENT (INOUT) :: msg(:,:)
        INTEGER, INTENT(IN) :: gid
#if defined(__MPI)
        INTEGER :: msglen
        msglen = size(msg)
        CALL reduce_base_integer( msglen, msg, gid, -1 )
#endif
      END SUBROUTINE mp_sum_im
!
!------------------------------------------------------------------------------!

      SUBROUTINE mp_sum_it(msg,gid)
        IMPLICIT NONE
        INTEGER, INTENT (INOUT) :: msg(:,:,:)
        INTEGER, INTENT (IN) :: gid
#if defined(__MPI)
        INTEGER :: msglen
        msglen = size(msg)
        CALL reduce_base_integer( msglen, msg, gid, -1 )
#endif
      END SUBROUTINE mp_sum_it

!------------------------------------------------------------------------------!

      SUBROUTINE mp_sum_i4(msg,gid)
        IMPLICIT NONE
        INTEGER, INTENT (INOUT) :: msg(:,:,:,:)
        INTEGER, INTENT (IN) :: gid
#if defined(__MPI)
        INTEGER :: msglen
        msglen = size(msg)
        CALL reduce_base_integer( msglen, msg, gid, -1 )
#endif
      END SUBROUTINE mp_sum_i4

!------------------------------------------------------------------------------!

      SUBROUTINE mp_sum_i5(msg,gid)
        IMPLICIT NONE
        INTEGER, INTENT (INOUT) :: msg(:,:,:,:,:)
        INTEGER, INTENT (IN) :: gid
#if defined(__MPI)
        INTEGER :: msglen
        msglen = size(msg)
        CALL reduce_base_integer( msglen, msg, gid, -1 )
#endif
      END SUBROUTINE mp_sum_i5


!------------------------------------------------------------------------------!

      SUBROUTINE mp_sum_r1(msg,gid)
        IMPLICIT NONE
        REAL (DP), INTENT (INOUT) :: msg
        INTEGER, INTENT (IN) :: gid
#if defined(__MPI)
        INTEGER :: msglen
        msglen = 1
        CALL reduce_base_real( msglen, msg, gid, -1 )
#endif
      END SUBROUTINE mp_sum_r1

!
!------------------------------------------------------------------------------!

      SUBROUTINE mp_sum_rv(msg,gid)
        IMPLICIT NONE
        REAL (DP), INTENT (INOUT) :: msg(:)
        INTEGER, INTENT (IN) :: gid
#if defined(__MPI)
        INTEGER :: msglen
        msglen = size(msg)
        CALL reduce_base_real( msglen, msg, gid, -1 )
#endif
      END SUBROUTINE mp_sum_rv
!
!------------------------------------------------------------------------------!


      SUBROUTINE mp_sum_rm(msg, gid)
        IMPLICIT NONE
        REAL (DP), INTENT (INOUT) :: msg(:,:)
        INTEGER, INTENT (IN) :: gid
#if defined(__MPI)
        INTEGER :: msglen
        msglen = size(msg)
        CALL reduce_base_real( msglen, msg, gid, -1 )
#endif
      END SUBROUTINE mp_sum_rm


      SUBROUTINE mp_root_sum_rm( msg, res, root, gid )
        IMPLICIT NONE
        REAL (DP), INTENT (IN)  :: msg(:,:)
        REAL (DP), INTENT (OUT) :: res(:,:)
        INTEGER,   INTENT (IN)  :: root
        INTEGER,   INTENT (IN) :: gid
#if defined(__MPI)
        INTEGER :: msglen, ierr, taskid

        msglen = size(msg)

        CALL mpi_comm_rank( gid, taskid, ierr)
        IF( ierr /= 0 ) CALL mp_stop( 8059 )
        !
        IF( taskid == root ) THEN
           IF( msglen > size(res) ) CALL mp_stop( 8060 )
        END IF

        CALL reduce_base_real_to( msglen, msg, res, gid, root )

#else

        res = msg

#endif

      END SUBROUTINE mp_root_sum_rm


      SUBROUTINE mp_root_sum_cm( msg, res, root, gid )
        IMPLICIT NONE
        COMPLEX (DP), INTENT (IN)  :: msg(:,:)
        COMPLEX (DP), INTENT (OUT) :: res(:,:)
        INTEGER,   INTENT (IN)  :: root
        INTEGER,  INTENT (IN) :: gid
#if defined(__MPI)
        INTEGER :: msglen, ierr, taskid

        msglen = size(msg)

        CALL mpi_comm_rank( gid, taskid, ierr)
        IF( ierr /= 0 ) CALL mp_stop( 8061 )

        IF( taskid == root ) THEN
           IF( msglen > size(res) ) CALL mp_stop( 8062 )
        END IF

        CALL reduce_base_real_to( 2 * msglen, msg, res, gid, root )

#else

        res = msg

#endif

      END SUBROUTINE mp_root_sum_cm

!
!------------------------------------------------------------------------------!


!------------------------------------------------------------------------------!
!

      SUBROUTINE mp_sum_rmm( msg, res, root, gid )
        IMPLICIT NONE
        REAL (DP), INTENT (IN) :: msg(:,:)
        REAL (DP), INTENT (OUT) :: res(:,:)
        INTEGER, INTENT (IN) :: root
        INTEGER, INTENT (IN) :: gid
        INTEGER :: group
        INTEGER :: msglen
        INTEGER :: taskid, ierr

        msglen = size(msg)

#if defined(__MPI)

        group = gid
        !
        CALL mpi_comm_rank( group, taskid, ierr)
        IF( ierr /= 0 ) CALL mp_stop( 8063 )

        IF( taskid == root ) THEN
           IF( msglen > size(res) ) CALL mp_stop( 8064 )
        END IF
        !
        CALL reduce_base_real_to( msglen, msg, res, group, root )
        !

#else
        res = msg
#endif

      END SUBROUTINE mp_sum_rmm


!
!------------------------------------------------------------------------------!


      SUBROUTINE mp_sum_rt( msg, gid )
        IMPLICIT NONE
        REAL (DP), INTENT (INOUT) :: msg(:,:,:)
        INTEGER, INTENT(IN) :: gid
#if defined(__MPI)
        INTEGER :: msglen
        msglen = size(msg)
        CALL reduce_base_real( msglen, msg, gid, -1 )
#endif
      END SUBROUTINE mp_sum_rt

!
!------------------------------------------------------------------------------!
!
! Carlo Cavazzoni
!
      SUBROUTINE mp_sum_r4d(msg,gid)
        IMPLICIT NONE
        REAL (DP), INTENT (INOUT) :: msg(:,:,:,:)
        INTEGER, INTENT(IN) :: gid
#if defined(__MPI)
        INTEGER :: msglen
        msglen = size(msg)
        CALL reduce_base_real( msglen, msg, gid, -1 )
#endif
      END SUBROUTINE mp_sum_r4d



!------------------------------------------------------------------------------!

      SUBROUTINE mp_sum_c1(msg,gid)
        IMPLICIT NONE
        COMPLEX (DP), INTENT (INOUT) :: msg
        INTEGER, INTENT(IN) :: gid
#if defined(__MPI)
        INTEGER :: msglen
        msglen = 1
        CALL reduce_base_real( 2 * msglen, msg, gid, -1 )
#endif
      END SUBROUTINE mp_sum_c1
!
!------------------------------------------------------------------------------!

      SUBROUTINE mp_sum_cv(msg,gid)
        IMPLICIT NONE
        COMPLEX (DP), INTENT (INOUT) :: msg(:)
        INTEGER, INTENT(IN) :: gid
#if defined(__MPI)
        INTEGER :: msglen
        msglen = size(msg)
        CALL reduce_base_real( 2 * msglen, msg, gid, -1 )
#endif
      END SUBROUTINE mp_sum_cv
!
!------------------------------------------------------------------------------!

      SUBROUTINE mp_sum_cm(msg, gid)
        IMPLICIT NONE
        COMPLEX (DP), INTENT (INOUT) :: msg(:,:)
        INTEGER, INTENT (IN) :: gid
#if defined(__MPI)
        INTEGER :: msglen
        msglen = size(msg)
        CALL reduce_base_real( 2 * msglen, msg, gid, -1 )
#endif
      END SUBROUTINE mp_sum_cm
!
!------------------------------------------------------------------------------!


      SUBROUTINE mp_sum_cmm(msg, res, gid)
        IMPLICIT NONE
        COMPLEX (DP), INTENT (IN) :: msg(:,:)
        COMPLEX (DP), INTENT (OUT) :: res(:,:)
        INTEGER, INTENT (IN) :: gid
#if defined(__MPI)
        INTEGER :: msglen
        msglen = size(msg)
        CALL reduce_base_real_to( 2 * msglen, msg, res, gid, -1 )
#else
        res = msg
#endif
      END SUBROUTINE mp_sum_cmm


!
!------------------------------------------------------------------------------!
!
! Carlo Cavazzoni
!
      SUBROUTINE mp_sum_ct(msg,gid)
        IMPLICIT NONE
        COMPLEX (DP), INTENT (INOUT) :: msg(:,:,:)
        INTEGER, INTENT(IN) :: gid
#if defined(__MPI)
        INTEGER :: msglen
        msglen = SIZE(msg)
        CALL reduce_base_real( 2 * msglen, msg, gid, -1 )
#endif
      END SUBROUTINE mp_sum_ct

!
!------------------------------------------------------------------------------!
!
! Carlo Cavazzoni
!
      SUBROUTINE mp_sum_c4d(msg,gid)
        IMPLICIT NONE
        COMPLEX (DP), INTENT (INOUT) :: msg(:,:,:,:)
        INTEGER, INTENT(IN) :: gid
#if defined(__MPI)
        INTEGER :: msglen
        msglen = size(msg)
        CALL reduce_base_real( 2 * msglen, msg, gid, -1 )
#endif
      END SUBROUTINE mp_sum_c4d
!
!------------------------------------------------------------------------------!
!
! Carlo Cavazzoni
!
      SUBROUTINE mp_sum_c5d(msg,gid)
        IMPLICIT NONE
        COMPLEX (DP), INTENT (INOUT) :: msg(:,:,:,:,:)
        INTEGER, INTENT(IN) :: gid
#if defined(__MPI)
        INTEGER :: msglen
        msglen = size(msg)
        CALL reduce_base_real( 2 * msglen, msg, gid, -1 )
#endif
      END SUBROUTINE mp_sum_c5d

!------------------------------------------------------------------------------!
!
! Carlo Cavazzoni
!
      SUBROUTINE mp_sum_r5d(msg,gid)
        IMPLICIT NONE
        REAL (DP), INTENT (INOUT) :: msg(:,:,:,:,:)
        INTEGER, INTENT(IN) :: gid
#if defined(__MPI)
        INTEGER :: msglen
        msglen = size(msg)
        CALL reduce_base_real( msglen, msg, gid, -1 )
#endif
      END SUBROUTINE mp_sum_r5d


      SUBROUTINE mp_sum_r6d(msg,gid)
        IMPLICIT NONE
        REAL (DP), INTENT (INOUT) :: msg(:,:,:,:,:,:)
        INTEGER, INTENT(IN) :: gid
#if defined(__MPI)
        INTEGER :: msglen
        msglen = size(msg)
        CALL reduce_base_real( msglen, msg, gid, -1 )
#endif
      END SUBROUTINE mp_sum_r6d

!
!------------------------------------------------------------------------------!
!
! Carlo Cavazzoni
!
      SUBROUTINE mp_sum_c6d(msg,gid)
        IMPLICIT NONE
        COMPLEX (DP), INTENT (INOUT) :: msg(:,:,:,:,:,:)
        INTEGER, INTENT(IN) :: gid
#if defined(__MPI)
        INTEGER :: msglen
        msglen = size(msg)
        CALL reduce_base_real( 2 * msglen, msg, gid, -1 )
#endif
      END SUBROUTINE mp_sum_c6d





!
!------------------------------------------------------------------------------!




!------------------------------------------------------------------------------!
!
!..mp_stop
!
  SUBROUTINE mp_stop(code)
    IMPLICIT NONE
    INTEGER, INTENT (IN) :: code
    INTEGER :: ierr
    WRITE( stdout, fmt='( "*** error in Message Passing (mp) module ***")' )
    WRITE( stdout, fmt='( "*** error msg:  ",A60)' ) TRIM( err_msg )
    WRITE( stdout, fmt='( "*** error code: ",I5)' ) code

#if defined(__MPI) 
    CALL mpi_abort(mpi_comm_world,code,ierr)
#endif
    STOP
  END SUBROUTINE mp_stop 
  
  
end module mp