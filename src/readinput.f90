#define __MPI
module readinput
#if defined __MPI 		
  use global_mpi
  use mp
#endif
  use kinds,     only : dp,dpc
  use constants, only : maxlen
  use parameters,only : calculation,verbosity,outdir,ldecoherence,Cdecoherence,lit_gmnvkq,lit_ephonon,&
                        lreadscfout,scfoutname,lreadphout,phoutname,lreadfildyn,fildyn,epwoutname,&
                        methodsh,lfeedback,naver,nstep,nsnap,pre_nstep,pre_dt,mix_thr,&
                        gamma,ld_fric,dt,temp,l_ph_quantum,prtgmnvkq,&
                        init_kx,init_ky,init_kz,init_hband,init_eband,&
                        llaser,efield_cart,w_laser,fwhm,eps_acustic,&
                        lsetthreads,mkl_threads,lelecsh,lholesh,lehpairsh,&
                        ieband_min,ieband_max,ihband_min,ihband_max,nefre_sh,nhfre_sh,&
												nnode,ncore,savedsnap,lsortpes,l_dEa_dQ,l_dEa2_dQ2,nsample,dirsample
  use io,        only : io_file_unit,io_error,msg
  use utility,   only : utility_lowercase

  implicit none
    integer               :: in_unit,tot_num_lines,loop,in1,in2
    integer               :: num_lines,line_counter
    character(len=maxlen),allocatable:: in_data(:) 
    character(len=maxlen) :: dummy,ctmp
    integer               :: ipos
    integer               :: ios , ios2
    character(len=512)    :: line
    character, parameter  :: TABCHAR = char(9) !char(9)为制表符TAB
  
  contains  
  
  subroutine get_inputfile(filename)
    implicit none
    character(*),intent(in) :: filename
#if defined __MPI 		
    if(ionode) then
#endif
      call treat_inputfile(filename)
      call read_namelist()
#if defined __MPI 		
    endif    
    call bcast_namelist()
#endif 

    call param_in_atomicunits()
    
  end subroutine

  !=======================================!
  subroutine treat_inputfile(filename)  
  !用于将输入文件中的每一条写入字符串文件，并且(将非文件路径)改为小写，去除注释          
  !=======================================!
  !! Load the shin file into a character  
  !! array in_file, ignoring comments and  
  !! blank lines and converting everything 
  !! to lowercase characters               
  !=======================================!

    implicit none
    character(*),intent(in):: filename
    logical :: alive
		integer :: istat
    !ia = ichar('a')
    !iz = ichar('z')
    
    in_unit=io_file_unit( )
    inquire(file=trim(adjustl(filename)),exist=alive)
    if(.NOT. alive) then
      call io_error("Error:Input file "//trim(adjustl(filename))//" doesn't exist.")
    else
      open (unit=in_unit, file='LVCSH.in',form='formatted',status='old',iostat=ierr)
      if(ierr /= 0) then
        call io_error('Error: Problem opening input file LVCSH.in')
      endif
    endif
    
    num_lines=0;tot_num_lines=0;ierr=0
    do while( ierr == 0 )
      read(in_unit, '(a)', iostat = ierr ,iomsg=msg) dummy   !参考p.177 P.540
      if(ierr > 0 ) then
        call io_error('Error: Problem reading input file SHIN')
        call io_error(msg)
      elseif(ierr == 0 )then
        
        ! convert all tabulation characters to spaces
        ipos = index(dummy,TABCHAR) !查询字符串在字符串中出现的位置,并将制表符改为空格
        do while (ipos /= 0)
          dummy(ipos:ipos) = ' '
          ipos = index(dummy,TABCHAR)
        end do
        ! 
        dummy=adjustl(dummy)
        
        tot_num_lines=tot_num_lines+1
        if( dummy(1:1)/='!'  .and. dummy(1:1)/='#' ) then
          if(len_trim(adjustl(dummy)) > 0 ) num_lines=num_lines+1
        endif
      endif
    end do
    !得到SHIN文件中总的行数tot_num_lines以及非注释和空行 num_lines

    rewind(in_unit)

    allocate(in_data(num_lines),stat=ierr)  !字符串数组，内部文件 line=449
    if (ierr/=0) call io_error('Error allocating in_data in param_in_file')

    line_counter=0
    do loop=1,tot_num_lines
      read(in_unit, '(a)', iostat = ierr ,iomsg=msg) dummy
        if(ierr /= 0) then
          write(stdout,*)'Error: Problem opening input file SHIN'
          call io_error(msg)
        endif
      !I convert all tabulation characters to spaces
      ipos = index(dummy,TABCHAR)
      do while (ipos /= 0)
        dummy(ipos:ipos) = ' '
        ipos = index(dummy,TABCHAR)
      end do
      !
      !if(index(dummy,"dir")<=0) then 
      !  dummy=utility_lowercase(dummy) !将(不表示文件路径的)字符串中大写字母全部改为小写
      !endif
      dummy=trim(adjustl(dummy))     !
      if( dummy(1:1)=='!' .or.  dummy(1:1)=='#' ) cycle
      if(len(trim(dummy)) == 0 ) cycle
      if(index(dummy,'=') <=1 )  cycle  !当该行中没有‘=’ 或‘=’前没有内容则跳过该行
      line_counter=line_counter+1
      
      !去除有效行信息中的注释部分，注释可以采用 ！或者 #
      in1=index(dummy,'!')
      in2=index(dummy,'#')
      if(in1==0 .and. in2==0)  in_data(line_counter)=dummy
      !不存在'!'与'#'
      if(in1==0 .and. in2>0 )  in_data(line_counter)=dummy(:in2-1)
      if(in2==0 .and. in1>0 )  in_data(line_counter)=dummy(:in1-1)
      if(in2> 0 .and. in1>0 )  in_data(line_counter)=dummy(:min(in1,in2)-1)
      
      !如果输入参数为字符串，则给字符串加上"*"
      !itmp = index(dummy,'=')
      !ctmp = dummy(:itmp-1)
      !ctmp1= dummy(itmp+1:)
      !ctmp1= trim(adjustl(ctmp1))
      !itmp = ichar(ctmp1(1:1))
      !if(itmp>=ia .and. itmp<=iz) then
      !  dummy = ctmp//"='"//ctmp1//"'"
      !endif
      
    end do
    !得到包含有效信息的行数line_counter,和相应的数据in_data(line_counter)

    close(in_unit)

  end subroutine treat_inputfile 
  
  subroutine read_namelist()
    use parameters,only : shinput
    use io
    implicit none
    integer::incar_unit,i,long,istat
    
    !!write input file to namelist input file
    incar_unit = io_file_unit()
    open(unit=incar_unit,status='SCRATCH',iostat=ierr,iomsg=msg)
    if(ierr > 0 ) then
      call io_error('Error: Problem reading SCRATCH input namelist file')
      call io_error(msg)
    elseif(ierr == 0 )then    
      write(incar_unit,*)"&shinput" 
      do i=1,line_counter
        write(incar_unit,"(A)") trim(adjustl(in_data(i)))
      enddo
      !write(incar_unit,"(A1)") "/"
      write(incar_unit,*) "/"
    endif
    rewind(incar_unit)
    
    !   set default values for variables in namelist
		calculation   = "lvcsh"
		verbosity     = "low"
		outdir        = "./"
    methodsh      = "FSSH"
    lfeedback     = .true.
		ldecoherence  = .true.
		cdecoherence  = 0.1
		lit_gmnvkq    = 1.0d-8 !meV
		lit_ephonon   = 1.0    !meV
    eps_acustic   = 5.d0 ! cm-1
    l_ph_quantum  = .true.
    lelecsh       = .false.
    lholesh       = .false.
    lehpairsh     = .false.
    lreadscfout   = .false.
    lsortpes      = .false.
    scfoutname    = "scf.out"
    lreadphout    = .false.
    phoutname     = "ph.out"
    lreadfildyn   = .false.
    fildyn        = "prefix.dyn"
    epwoutname    = "epw.out"
    prtgmnvkq     = .false.
    ieband_min    = 0
    ieband_max    = 0
    ihband_min    = 0
    ihband_max    = 0
		nefre_sh      = 0
		nhfre_sh      = 0
    naver         = 10
    nstep         = 10
    nsnap         = 100
    pre_nstep     = 0
    pre_dt        = 1.0
		mix_thr       = 0.8
    gamma         = 0.0    ! 0.1   ! the friction coefficient 1/ps
    ld_fric       = 0.01
    dt            = 0.5    ! fs
    temp          = 300.0  ! K
    init_kx       = 0.0    ! in unit of b_x
    init_ky       = 0.0
    init_kz       = 0.0
    init_eband    = 1      ! the initial electron band
    init_hband    = 1      ! the initial hole band
    llaser        = .true.
    efield_cart   = (/ 0.0,0.0,0.0 /)  ! V/m
    w_laser       = 1.0    ! eV
    fwhm          = 10     ! fs
		nnode         = 1
		ncore         = 1
		savedsnap     = 100
    nsample       = 1
    dirsample     = "sample"
    
    l_dEa_dQ      = .true.
    l_dEa2_dQ2    = .true.
    
    lsetthreads   = .FALSE.
    mkl_threads   = 4    
    
    write(stdout,"(/,1X,A77)")   repeat("=",77)
    write(stdout,"(1X,10X,A)") "The namelist file as follows"
    write(stdout,"(1X,A77)")   repeat("=",77)
    do i=1,line_counter+2
      read(incar_unit,"(A)") ctmp
      write(stdout,"(A)") ctmp
    enddo
    rewind(incar_unit)
    read(UNIT=incar_unit,nml=shinput,iostat=ios,iomsg=msg)
    ios2 = 0
    if(ios /= 0) then
      backspace(incar_unit)
      read(incar_unit,'(A512)',iostat=ios2) line
    endif
    if(ios2 /= 0) call errore("readinput","Could not fine namelist &shinput",2)
    if(ios /= 0) then
      call errore("readinput",'Bad line in namelist &shinput:'//&
                  trim(line)//" (error could be in the previous line)",1)
      write(stdout,*)'Error: Problem reading namelist file SHIN'
      call io_error(msg)
    endif  
    close(incar_unit)
		
		long = len_trim(adjustl(outdir))
		if(long == 0) write(stdout,"(5X,A)") "Warring! the input of outdir name empty "
		outdir = trim(adjustl(outdir))
		if(outdir(long:long) /= "/") then
			outdir(long+1:long+1) = "/"
		endif
		!call system("mkdir "//trim(outdir))
		
    if(lehpairsh) then
      lelecsh = .true.
      lholesh = .true.
    endif
    write(stdout,"(/,1X,A)") "Read parameter Successful!" 
    write(stdout,"(1X,A77)")   repeat("=",77)
    
  end subroutine read_namelist
  
  subroutine param_in_atomicunits()
    ! change the input parameters into Rydberg atomic units
    use constants,only : Ryd2V_m,mev2cmm1,ryd2eV,Ry_TO_THZ,Ry_TO_fs,ryd2mev
    use lasercom,only  : fwhm_2T2
    implicit none
    
		cdecoherence = cdecoherence/2.0 ! Hartree to Rydberg
    lit_gmnvkq   = lit_gmnvkq/ryd2meV
    !
    ! from cm-1 to meV
    eps_acustic = eps_acustic / mev2cmm1
		
    gamma = gamma / Ry_TO_THZ ! change gamma in unit of (THZ) to (Ryd)
    dt    = dt / Ry_TO_fs     ! in Rydberg atomic units(1 a.u.=4.8378 * 10^-17 s)
    pre_dt= pre_dt/ Ry_TO_fs
    
    efield_cart = efield_cart / Ryd2V_m ! (in Ry a.u.;1 a.u. = 36.3609*10^10 V/m)
    w_laser = w_laser /ryd2eV
    fwhm = fwhm/ ry_to_fs  
    fwhm_2T2 = fwhm**2.0/4.0*log(2.0)
    
  end subroutine param_in_atomicunits


#if defined __MPI   
  subroutine bcast_namelist()
    implicit none
    
    call mp_bcast(calculation   ,ionode_id)  
    call mp_bcast(lfeedback     ,ionode_id)
    call mp_bcast(verbosity     ,ionode_id)
    call mp_bcast(outdir        ,ionode_id)
    call mp_bcast(methodsh      ,ionode_id)
    call mp_bcast(ldecoherence  ,ionode_id)
    call mp_bcast(cdecoherence  ,ionode_id)
    call mp_bcast(lit_gmnvkq    ,ionode_id)
    call mp_bcast(lit_ephonon   ,ionode_id)
    call mp_bcast(eps_acustic   ,ionode_id)
    call mp_bcast(l_ph_quantum  ,ionode_id)
    call mp_bcast(lelecsh       ,ionode_id)
    call mp_bcast(lholesh       ,ionode_id)
    call mp_bcast(lehpairsh     ,ionode_id)
    call mp_bcast(lreadscfout   ,ionode_id)
    call mp_bcast(lsortpes      ,ionode_id)
    call mp_bcast(scfoutname    ,ionode_id)
    call mp_bcast(lreadphout    ,ionode_id)
    call mp_bcast(phoutname     ,ionode_id)
    call mp_bcast(lreadfildyn   ,ionode_id)
    call mp_bcast(fildyn        ,ionode_id)
    call mp_bcast(epwoutname    ,ionode_id)
    call mp_bcast(prtgmnvkq     ,ionode_id)
    call mp_bcast(ieband_min    ,ionode_id)
    call mp_bcast(ieband_max    ,ionode_id)
    call mp_bcast(ihband_min    ,ionode_id)
    call mp_bcast(ihband_max    ,ionode_id)
    call mp_bcast(nefre_sh      ,ionode_id)
    call mp_bcast(nhfre_sh      ,ionode_id)
    call mp_bcast(naver         ,ionode_id)
    call mp_bcast(nstep         ,ionode_id)
    call mp_bcast(nsnap         ,ionode_id)
    call mp_bcast(pre_nstep     ,ionode_id)
    call mp_bcast(pre_dt        ,ionode_id) 
    call mp_bcast(mix_thr       ,ionode_id)
    call mp_bcast(gamma         ,ionode_id) 
    call mp_bcast(ld_fric       ,ionode_id)
    call mp_bcast(dt            ,ionode_id) 
    call mp_bcast(temp          ,ionode_id)
    call mp_bcast(init_kx       ,ionode_id) 
    call mp_bcast(init_ky       ,ionode_id)
    call mp_bcast(init_kz       ,ionode_id) 
    call mp_bcast(init_eband    ,ionode_id)
    call mp_bcast(init_hband    ,ionode_id)    
    call mp_bcast(llaser        ,ionode_id)
    call mp_bcast(efield_cart   ,ionode_id) 
    call mp_bcast(w_laser       ,ionode_id)
    call mp_bcast(fwhm          ,ionode_id) 
    call mp_bcast(nnode         ,ionode_id)
    call mp_bcast(ncore         ,ionode_id) 
    call mp_bcast(savedsnap     ,ionode_id)
    call mp_bcast(l_dEa_dQ      ,ionode_id) 
    call mp_bcast(l_dEa2_dQ2    ,ionode_id)
    call mp_bcast(lsetthreads   ,ionode_id) 
    call mp_bcast(mkl_threads   ,ionode_id)
    call mp_bcast(nsample       ,ionode_id)
    call mp_bcast(dirsample     ,ionode_id)


  if(trim(calculation) == "lvcsh" .and. verbosity  == "high" ) then
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


    write(procout,"(/,1X,A15,I5,1X,A10)") "Input file for ",iproc,"processor."   
    write(procout,*) "calculation  =",trim(calculation) 
    write(procout,*) "lfeedback    =",lfeedback     
    write(procout,*) "verbosity    =",trim(verbosity)  
    write(procout,*) "outdir       =",trim(outdir)       
    write(procout,*) "methodsh     =",methodsh   
    write(procout,*) "ldecoherence =",ldecoherence    
    write(procout,*) "cdecoherence =",cdecoherence 
    write(procout,*) "lit_gmnvkq   =",lit_gmnvkq  
    write(procout,*) "lit_ephonon  =",lit_ephonon 
    write(procout,*) "eps_acustic  =",eps_acustic    
    write(procout,*) "l_ph_quantum =",l_ph_quantum
    write(procout,*) "lelecsh      =",lelecsh     
    write(procout,*) "lholesh      =",lholesh    
    write(procout,*) "lehpairsh    =",lehpairsh     
    write(procout,*) "lreadscfout  =",lreadscfout
    write(procout,*) "lsortpes     =",lsortpes    
    write(procout,*) "scfoutname   =",trim(scfoutname) 
    write(procout,*) "lreadphout   =",lreadphout    
    write(procout,*) "phoutname    =",trim(phoutname)  
    write(procout,*) "lreadfildyn  =",lreadfildyn 
    write(procout,*) "fildyn       =",trim(fildyn)     
    write(procout,*) "epwoutname   =",trim(epwoutname)  
    write(procout,*) "prtgmnvkq    =",prtgmnvkq
    write(procout,*) "ieband_min   =",ieband_min 
    write(procout,*) "ieband_max   =",ieband_max  
    write(procout,*) "ihband_min   =",ihband_min 
    write(procout,*) "ihband_max   =",ihband_max    
    write(procout,*) "nefre_sh     =",nefre_sh   
    write(procout,*) "nhfre_sh     =",nhfre_sh    
    write(procout,*) "naver        =",naver      
    write(procout,*) "nstep        =",nstep         
    write(procout,*) "nsnap        =",nsnap      
    write(procout,*) "pre_nstep    =",pre_nstep  
    write(procout,*) "pre_dt       =",pre_dt     
    write(procout,*) "mix_thr      =",mix_thr       
    write(procout,*) "gamma        =",gamma      
    write(procout,*) "ld_fric      =",ld_fric     
    write(procout,*) "dt           =",dt         
    write(procout,*) "temp         =",temp         
    write(procout,*) "init_kx      =",init_kx      
    write(procout,*) "init_ky      =",init_ky      
    write(procout,*) "init_kz      =",init_kz      
    write(procout,*) "init_eband   =",init_eband   
    write(procout,*) "init_hband   =",init_hband   
    write(procout,*) "llaser       =",llaser       
    !write(procout,*) "efield_cart  =",efield_cart 
    write(procout,"(1X,A14,3(1X,F12.5))") "efield_cart  =",efield_cart  
    write(procout,*) "w_laser      =",w_laser      
    write(procout,*) "fwhm         =",fwhm         
    write(procout,*) "nnode        =",nnode        
    write(procout,*) "ncore        =",ncore            
    write(procout,*) "savedsnap    =",savedsnap    
    write(procout,*) "l_dEa_dQ     =",l_dEa_dQ     
    write(procout,*) "l_dEa2_dQ2   =",l_dEa2_dQ2   
    write(procout,*) "lsetthreads  =",lsetthreads  
    write(procout,*) "mkl_threads  =",mkl_threads  
    
    if(naver < nproc) then
      if(ionode) write(stdout,*) "The parameter naver need to set larger than num of processor in MPI."
      write(procout,*) "The parameter naver need to set larger than num of processor in MPI."
      stop
    endif
  
  else
    if(ionode) write(stdout,*) "Together the results in ",nsample,"samples."
    if(ionode) write(stdout,*) "The prefix directory name of each ensemble is :",trim(dirsample)
  endif  
    
  end subroutine bcast_namelist
#endif 
  
  function get_ik(kx,nkx)
    use kinds,only : dp
    implicit none
    real(kind=dp) :: kx
    integer       :: nkx
    integer       :: get_ik
    kx = MOD(kx,1.0)
    if(kx<=0.0) kx = kx + 1.0
    get_ik = Anint(kx*nkx)+1
    if(get_ik<1) then
      get_ik = get_ik +nkx
    elseif(get_ik>nkx) then
      get_ik = get_ik - nkx
    endif
    
  end function
  
end module readinput