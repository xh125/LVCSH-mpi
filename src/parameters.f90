module parameters
  !! This module contains parameters to control the actions of SCSH.
  !! Also routines to read the parameters and write them out again.
  use kinds,only : dp
  use constants,only    : maxlen
  use environments,only : mkl_threads,lsetthreads
  use lasercom,only     : llaser,lfcw,efield_cart,w_laser,fwhm
  use surfacecom,only   : methodsh,lfeedback,naver,nstep,nsnap,dt,pre_nstep,pre_dt,&
                          temp,gamma,ld_fric,l_ph_quantum,lit_gmnvkq,lit_ephonon,&
                          lelecsh,lholesh,lehpairsh,ldecoherence,Cdecoherence,eps_acustic,&
                          ieband_min,ieband_max,ihband_min,ihband_max,nefre_sh,nhfre_sh,&
                          l_dEa_dQ,l_dEa2_dQ2
  implicit none
  
  !integer,parameter :: npk = 40000 ! max number of k-points in pw.x calculation
  !integer,parameter :: ntypx = 10  ! max number of different types of atom
  character(len=maxlen) :: scfoutname,phoutname,fildyn,epwoutname

  real(kind=dp)   :: init_kx,init_ky,init_kz  !激发后初始的电子(空穴)的k坐标
  integer         :: init_hband,init_eband    !激发后初始的电子(空穴)所处的能带
  real(kind=dp)   :: init_e_en,init_h_en      !激发后初始电子(空穴)的能量
  integer         :: init_ikx,init_iky,init_ikz,init_ik 
  real(kind=dp)   :: mix_thr
	! threshold to find the mixxing states
	character(len=maxlen) :: inputfilename = "LVCSH.in"
	character(len=maxlen) :: calculation
	character(len=maxlen) :: verbosity
	character(len=maxlen) :: outdir
	integer :: nnode,ncore,naver_sum,savedsnap
  integer :: nsample,isample
  character(len=maxlen) :: dirsample
  logical :: lreadscfout,lreadphout,lreadfildyn,lsortpes
  logical :: prtgmnvkq
  
  namelist / shinput / &
           calculation,verbosity,outdir,ldecoherence,Cdecoherence,lit_gmnvkq,lit_ephonon,&
					 lreadscfout,lreadphout,scfoutname,phoutname,lreadfildyn,fildyn,&
					 epwoutname,methodsh,lfeedback,naver,nstep,nsnap,mix_thr,&
           pre_nstep,pre_dt,gamma,ld_fric,l_ph_quantum,dt,temp,&
           init_kx,init_ky,init_kz,init_hband,init_eband,&
           llaser,lfcw,efield_cart,w_laser,fwhm,eps_acustic,prtgmnvkq,&
           lsetthreads,mkl_threads,lelecsh,lholesh,lehpairsh,&
           ieband_min,ieband_max,ihband_min,ihband_max,nefre_sh,nhfre_sh,&
					 nnode,ncore,savedsnap,lsortpes,l_dEa_dQ,l_dEa2_dQ2,&
           nsample,dirsample

  contains
  

  
end module parameters   