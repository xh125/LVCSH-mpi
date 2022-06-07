module global_mpi
  use mpi
  use io,only : io_file_unit
  use constants ,only : maxlen  
  implicit none
  integer :: ierr
  integer :: iproc,nproc,procout
  character(len=10) :: iproc_str
  character(len=256) :: job_path,work_path,procout_name
  ! job_path   The absolute path of the folder you submit the job. This wll be got from pwd()
  ! work_path  The current working path for each processor. If DINT_TMP_DIR is defined, 
  logical :: ionode
  integer :: ionode_id
  logical :: llinux
  logical :: filelive,dirlive
  integer :: ver_mpi,subver_mpi
  
  contains
  
  subroutine initmpi
    implicit none
    
    llinux = .true.
    ionode = .false.
    !initialize mpi
    call MPI_Init(ierr)
    call MPI_COMM_RANK(MPI_COMM_WORLD,iproc,ierr)
    call MPI_COMM_SIZE(MPI_COMM_WORLD,nproc,ierr)
    call allocate_buffers()
    write(iproc_str,"(i0.3)") iproc
    ionode_id = 0
    if(iproc == 0) ionode = .true.
    call MPI_get_version(ver_mpi,subver_mpi,ierr)
    !if(ionode) write(*,*) "MPI Version:",ver_mpi,".",subver_mpi
    
    !get the job_path
    call getcwd(job_path)
    !call MPI_Barrier(MPI_COMM_WORLD,ierr)
    if(ionode) call findsystem(job_path,llinux)
    call MPI_Bcast(llinux,1,MPI_logical,ionode_id,MPI_COMM_WORLD,ierr)
    
    !if(ionode) write(*,*) "Job direct is:",job_path
    
    !creat work_path
    if(llinux) then
      work_path = trim(job_path)//"/"//iproc_str
    else
      work_path = trim(job_path)//"\"//iproc_str
    endif
    
    !inquire(file=trim(adjustl(work_path)),exist=filelive)
    inquire(directory=trim(adjustl(work_path)),exist=dirlive)
    !write(*,*) "dirlive:",dirlive
    if(.not. dirlive) call system("mkdir "//trim(work_path))
    !write(*,*) "Work direct for ",iproc," processor is :",work_path
    procout = io_file_unit()
    if(llinux) then
      procout_name = trim(work_path)//"/LVCSH_"//trim(adjustl(iproc_str))//".out"
    else
      procout_name = trim(work_path)//"\LVCSH_"//trim(adjustl(iproc_str))//".out"
    endif
    !write(*,*) trim(procout_name)
    open(unit=procout,file=procout_name,status='REPLACE')
    write(procout,*) "Output of the ",iproc,"processor."
    write(procout,*) "Job_path is :",trim(job_path)
    write(procout,*) "Work_path is :",trim(work_path)
    
  end subroutine initmpi
  
  subroutine endmpi

    implicit none
    
    character(len=maxlen) :: msg
    
    close(unit=procout,iostat=ierr,iomsg=msg)
    call mpi_finalize(ierr)
    
  end subroutine endmpi
  
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