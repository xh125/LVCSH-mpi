module mp
  use kinds,only : dp,dpc
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

  public :: mp_bcast
 

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

    msglen = 1
    group = mpi_comm_world
    IF( PRESENT( gid ) ) group = gid
    CALL bcast_integer( msg, msglen, source, group )
  
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
  
    msglen = size(msg)
    group = mpi_comm_world
    IF( PRESENT( gid ) ) group = gid
    CALL bcast_integer( msg, msglen, source, group )
  
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
  
    msglen = size(msg)
    group = mpi_comm_world
    IF( PRESENT( gid ) ) group = gid
    CALL bcast_integer( msg, msglen, source, group )
  
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
  
    msglen = size(msg)
    group = mpi_comm_world
    IF( PRESENT( gid ) ) group = gid
    CALL bcast_integer( msg, msglen, source, group )
  
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
  
    msglen = 1
    group = mpi_comm_world
    IF( PRESENT( gid ) ) group = gid
    CALL bcast_real( msg, msglen, source, group )
  
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

  
    msglen = size(msg)
    group = mpi_comm_world
    IF( PRESENT( gid ) ) group = gid
    CALL bcast_real( msg, msglen, source, group )
  
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
  
    msglen = size(msg)
    group = mpi_comm_world
    IF( PRESENT( gid ) ) group = gid
    CALL bcast_real( msg, msglen, source, group )
  
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
  
    msglen = size(msg)
    group = mpi_comm_world
    IF( PRESENT( gid ) ) group = gid
    CALL bcast_real( msg, msglen, source, group )
  
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
  
    msglen = size(msg)
    group = mpi_comm_world
    IF( PRESENT( gid ) ) group = gid
    CALL bcast_real( msg, msglen, source, group )
  
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
  
    msglen = size(msg)
    group = mpi_comm_world
    IF( PRESENT( gid ) ) group = gid
    CALL bcast_real( msg, msglen, source, group )
  
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
  
    msglen = 1
    group = mpi_comm_world
    IF( PRESENT( gid ) ) group = gid
    CALL bcast_cmpl( msg, msglen, source, group )
  
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
  
    msglen = size(msg)
    group = mpi_comm_world
    IF( PRESENT( gid ) ) group = gid
    CALL bcast_cmpl( msg, msglen, source, group )
  
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
  
    msglen = size(msg)
    group = mpi_comm_world
    IF( PRESENT( gid ) ) group = gid
    CALL bcast_cmpl( msg, msglen, source, group )
  
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
  
    msglen = size(msg)
    group = mpi_comm_world
    IF( PRESENT( gid ) ) group = gid
    CALL bcast_cmpl( msg, msglen, source, group )
  
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
  
    msglen = size(msg)
    group = mpi_comm_world
    IF( PRESENT( gid ) ) group = gid
    CALL bcast_cmpl( msg, msglen, source, group )
  
  END SUBROUTINE mp_bcast_c4d

  SUBROUTINE mp_bcast_c5d(msg,source,gid)
    IMPLICIT NONE
    COMPLEX (kind=dpc) :: msg(:,:,:,:,:)
    INTEGER :: source
    INTEGER, OPTIONAL, INTENT(IN) :: gid
    INTEGER :: group
    INTEGER :: msglen
  
    msglen = size(msg)
    group = mpi_comm_world
    IF( PRESENT( gid ) ) group = gid
    CALL bcast_cmpl( msg, msglen, source, group )
  
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
  
    msglen = 1
    group = mpi_comm_world
    IF( PRESENT( gid ) ) group = gid
    CALL bcast_logical( msg, msglen, source, group )
    
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
  
    msglen = size(msg)
    group = mpi_comm_world
    IF( PRESENT( gid ) ) group = gid
    CALL bcast_logical( msg, msglen, source, group )
  
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
  
    msglen = size(msg)
    group = mpi_comm_world
    IF( PRESENT( gid ) ) group = gid
    CALL bcast_logical( msg, msglen, source, group )
  
  END SUBROUTINE mp_bcast_lm


!
!--------------------------------------------------------------------------!
!
  SUBROUTINE mp_bcast_z(msg,source,gid)
    IMPLICIT NONE
    CHARACTER (len=*) :: msg
    INTEGER :: source
    INTEGER, OPTIONAL, INTENT(IN) :: gid
    INTEGER :: group
    INTEGER :: msglen, ierr, i
    INTEGER, ALLOCATABLE :: imsg(:)
  
    ierr = 0
    msglen = len(msg)
    !write(*,*) "iproc",iproc," length:",msglen,trim(adjustl(msg))
    group = mpi_comm_world
    IF( PRESENT( gid ) ) group = gid
    IF (ierr/=0) CALL mp_stop( 8014 )
    ALLOCATE (imsg(1:msglen), STAT=ierr)
    IF (ierr/=0) CALL mp_stop( 8015 )
    DO i = 1, msglen
      imsg(i) = ichar(msg(i:i))
    END DO
    

    CALL bcast_integer( imsg, msglen, source, group )
    !write(*,*) "It's OK!"
    DO i = 1, msglen
      msg(i:i) = char(imsg(i))
    END DO

    !write(*,*) "iproc",iproc," length:",msglen,trim(adjustl(msg))
    
    DEALLOCATE (imsg, STAT=ierr)
    IF (ierr/=0) CALL mp_stop( 8016 )

    !call MPI_Barrier(group,ierr)   
    !write(*,*) "iproc",iproc," length:",msglen,trim(adjustl(msg))
    !call MPI_Barrier(group,ierr)

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
  
  END SUBROUTINE mp_bcast_zv
!
! end mp_bcast  
 

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

    CALL mpi_abort(mpi_comm_world,code,ierr)
    STOP
  END SUBROUTINE mp_stop 
  
  
end module mp