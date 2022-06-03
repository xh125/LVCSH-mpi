module global_mpi
  use mpi
  implicit none
  integer :: ierr
  integer :: iproc,nproc
  character(len=10) :: iproc_str
  character(len=256) :: job_path,work_path
  ! job_path   The absolute path of the folder you submit the job. This wll be got from pwd()
  ! work_path  The current working path for each processor. If DINT_TMP_DIR is defined, 
  logical :: ionode
  logical :: llinux
  integer :: filelive
  
  contains
  
  subroutine initmpi
    implicit none
    
    llinux = .true.
    ionode = .false.
    !initialize mpi
    call MPI_Init(ierr)
    call MPI_COMM_RANK(MPI_COMM_WORLD,iproc,ierr)
    call MPI_COMM_SIZE(MPI_COMM_WORLD,nproc,ierr)
    write(iproc_str,"(i0.3)") iproc
    if(iproc == 0) ionode = .true.
    
    !get the job_path
    call getcwd(job_path)
    !call MPI_Barrier(MPI_COMM_WORLD,ierr)
    call findsystem(job_path,llinux)
    !call MPI_Bcast(llinux,1,MPI_logical,ionode,MPI_COMM_WORLD,ierr)
    
    write(*,*) "Job direct is:",job_path
    
    !creat work_path
    if(llinux) then
      work_path = trim(job_path)//"/"//iproc_str
    else
      work_path = trim(job_path)//"\"//iproc_str
    endif
    
    inquire(file=trim(adjustl(work_path)),exist=filelive)
    if(filelive == 0) call system("mkdir "//trim(work_path))
    write(*,*) "Work direct for ",iproc," processor is :",work_path
    
  end subroutine initmpi
  
  subroutine findsystem(job_path,llinux)
    implicit none
    character(len=*),intent(inout) :: job_path
    logical,intent(inout) :: llinux
    
    integer :: l,i
    l = len_trim(adjustl(job_path))
    do i=1,l
      if(job_path(i:i) == "\") llinux = .false.
    enddo
    
  end subroutine findsystem
  
end module global_mpi