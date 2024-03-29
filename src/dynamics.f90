#define __MPI
module dynamics
#if defined __MPI 		
  use global_mpi
#endif
  use kinds,only     : dp,dpc
  use constants,only : cone,czero
  use parameters,only: temp,verbosity
  use epwcom,only    : kqmap
  
  implicit none
  
  contains

  subroutine set_gamma(nmodes,nq,gamma,ld_fric,wf,ld_gamma)
    use constants,only : eps10
    implicit none
    integer , intent(in) :: nmodes,nq
    real(kind=dp),intent(in) :: gamma
    real(kind=dp),intent(in) :: ld_fric
    real(kind=dp),intent(in) :: wf(nmodes,nq)
    real(kind=dp),intent(out):: ld_gamma(nmodes,nq) 
    
    if(abs(gamma) > eps10) then 
      !gamma /= 0.0
      ld_gamma = gamma
    else
      ld_gamma = ld_fric * wf
    endif
    
    
  end subroutine
  
  
  !=======================================================================!
  != rk4 method to obtain coordinate and velocitie after a time interval =!
  !=======================================================================!
  != ref: http://en.wikipedia.org/wiki/runge_kutta_methods               =!
  !=======================================================================!
  subroutine rk4_nuclei(nmodes,nq,dEa_dQ,ld_gamma,wf,xx,vv,tt)
    use parameters,only : lit_ephonon
    use constants,only : ryd2mev,czero
    implicit none
    integer , intent(in) :: nmodes,nq
    real(kind=dp),intent(in)   :: ld_gamma(nmodes,nq)
    complex(kind=dpc),intent(in)   :: dEa_dQ(nmodes,nq)
    real(kind=dp),intent(in)   :: wf(nmodes,nq)
    complex(kind=dpc),intent(inout):: xx(nmodes,nq),vv(nmodes,nq)
    real(kind=dp),intent(in)   :: tt
    real(kind=dp):: tt2,tt6
    complex(kind=dpc):: xx0(nmodes,nq),dx1(nmodes,nq),dx2(nmodes,nq),dx3(nmodes,nq),dx4(nmodes,nq)
    complex(kind=dpc):: vv0(nmodes,nq),dv1(nmodes,nq),dv2(nmodes,nq),dv3(nmodes,nq),dv4(nmodes,nq)
    integer :: iq,imode
    real(kind=dp) :: womiga

    tt2=tt/2.0d0; tt6=tt/6.0d0
    call derivs_nuclei(nmodes,nq,dEa_dQ,wf,ld_gamma,xx,vv,dx1,dv1)
    xx0=xx+tt2*dx1; vv0=vv+tt2*dv1
    call derivs_nuclei(nmodes,nq,dEa_dQ,wf,ld_gamma,xx0,vv0,dx2,dv2)
    xx0=xx+tt2*dx2; vv0=vv+tt2*dv2
    call derivs_nuclei(nmodes,nq,dEa_dQ,wf,ld_gamma,xx0,vv0,dx3,dv3)
    xx0=xx+tt*dx3; vv0=vv+tt*dv3
    call derivs_nuclei(nmodes,nq,dEa_dQ,wf,ld_gamma,xx0,vv0,dx4,dv4)
    xx=xx+tt6*(dx1+2.0d0*dx2+2.0d0*dx3+dx4)
    vv=vv+tt6*(dv1+2.0d0*dv2+2.0d0*dv3+dv4)
    
    do iq=1,nq
      do imode=1,nmodes
        womiga = wf(imode,iq)
        if(womiga == 0.0 ) then
          xx(imode,iq) = czero
          vv(imode,iq) = czero
        endif
      enddo
    enddo
    
    ! ph_Q(v,-q)=ph_Q(v,q)*    (2.48)
    ! ph_P(v,-q)=ph_P(v,q)*    (2.57)
    
    call test_xv(nmodes,nq,xx,vv) 
   
  end subroutine rk4_nuclei
  
  !====================================================!
  != calculate derivative of coordinate and velocitie =!
  !====================================================!
  != ref: notebook page 630                           =!
  !====================================================!
  ! ref : Huangkun <固体物理> (3-10)
  ! ref : 1 D. M. F. M. Germana Paterlini, Chemical Physics 236 (1998) 243.
  ! ref : PPT-92
  ! ref : 1 J. Qiu, X. Bai, and L. Wang, The Journal of Physical Chemistry Letters 9 (2018) 4319.
  ! ref : "<Phonons Theory and Experiments I Lattice Dynamics and Models of Interatomic Forces by Dr. Peter Brüesch (auth.) (z-lib.org).pdf>."
  subroutine derivs_nuclei(nmodes,nq,dEa_dQ,wf,ld_gamma,xx,vv,dx,dv)
    use elph2,only : iminusq
    implicit none
    integer,intent(in) :: nmodes,nq
    complex(kind=dpc),intent(in) :: dEa_dQ(nmodes,nq)
    real(kind=dp),intent(in)  ::  wf(nmodes,nq)
    real(kind=dp),intent(in)  ::  ld_gamma(nmodes,nq)
    complex(kind=dpc),intent(in)  ::  xx(nmodes,nq),vv(nmodes,nq)
    complex(kind=dpc),intent(out) ::  dx(nmodes,nq),dv(nmodes,nq)
    
    integer :: iq,imode,ik,iband1,iband2,ikq
    
    dx = vv    ! (2.55) (2.59)
    dv = -wf**2*xx - dEa_dQ - ld_gamma*vv
    !(2.60) PPT-94
    
  endsubroutine derivs_nuclei
  

  
  subroutine derivs_electron_diabatic(nfre,HH,cc,dc)
    use f95_precision
    use blas95
    use constants,only : cmplx_i,cmplx_0
    implicit none
    integer,intent(in)           :: nfre
    complex(kind=dp),intent(in)  :: HH(nfre,nfre)
    complex(kind=dpc),intent(in) :: cc(nfre)
    complex(kind=dpc),intent(out):: dc(nfre)
    
    dc= cmplx_0
    
    !dc = MATMUL(HH,c)
    call gemv(HH,cc,dc)
  
    dc = dc *(-cmplx_i)

  endsubroutine derivs_electron_diabatic

  !===========================================================!
  != rk4 method to obtain wavefunction after a time interval =!
  !===========================================================!
  != ref: http://en.wikipedia.org/wiki/runge_kutta_methods   =!
  !===========================================================!

  subroutine rk4_electron_diabatic(nfre,HH,cc,cc0,dc1,dc2,dc3,dc4,tt)
    implicit none
    integer,intent(in)              :: nfre
    complex(kind=dpc),intent(inout) :: cc(nfre)
    complex(kind=dpc),intent(in)    :: HH(nfre,nfre)
    complex(kind=dpc),intent(out)   :: cc0(nfre),dc1(nfre),&
                        dc2(nfre),dc3(nfre),dc4(nfre)
    real(kind=dp),intent(in)        :: tt
    real(kind=dp):: tt2,tt6
    real(kind=dp):: sum_cc2
    
    tt2=tt/2.0d0; tt6=tt/6.0d0
    
    call derivs_electron_diabatic(nfre,HH,cc,dc1)
    cc0=cc+tt2*dc1
    !sum_cc2 = REAL(SUM(cc0*CONJG(cc0)))
    !cc0= cc0/sqrt(sum_cc2)
    call derivs_electron_diabatic(nfre,HH,cc0,dc2)
    cc0=cc+tt2*dc2
    !sum_cc2 = REAL(SUM(cc0*CONJG(cc0)))
    !cc0= cc0/sqrt(sum_cc2)   
    call derivs_electron_diabatic(nfre,HH,cc0,dc3)
    cc0=cc+tt*dc3     
    !sum_cc2 = REAL(SUM(cc0*CONJG(cc0)))
    !cc0= cc0/sqrt(sum_cc2)
    call derivs_electron_diabatic(nfre,HH,cc0,dc4)
    
    cc=cc+tt6*(dc1+2.0d0*dc2+2.0d0*dc3+dc4)
    
    sum_cc2 = REAL(SUM(cc*CONJG(cc)))
    cc = cc/sqrt(sum_cc2)
    
  endsubroutine rk4_electron_diabatic
      
  !===========================================================!
  != rk4 method to obtain wavefunction after a time interval =!
  !===========================================================!
  != ref: http://en.wikipedia.org/wiki/runge_kutta_methods   =!
  !===========================================================!


  subroutine get_dEa2_dQ2(nmodes,nq,nfre,nfre_sh,isurface,EE,dd,dEa2_dQ2)
    use elph2 , only : iminusq
    implicit none
    integer,intent(in) :: nmodes,nq,nfre,nfre_sh
    integer,intent(in) :: isurface
    real(kind=dp),intent(in)  :: EE(nfre)
    complex(kind=dpc),intent(in)  :: dd(2,nfre_sh,nmodes,nq)
    real(kind=dp),intent(out) :: dEa2_dQ2(nmodes,nq)
    
    integer :: iq,imode,ifre,iq_
    
    ! ref : (2.60) PPT-95
    dEa2_dQ2 = 0.0
    do iq=1,nq
      iq_ = iminusq(iq)
      do imode=1,nmodes
        do ifre=1,nfre_sh
          if(ifre /= isurface) then         
            dEa2_dQ2(imode,iq) = dEa2_dQ2(imode,iq) + &
            & 2.0*(EE(isurface)-EE(ifre))*REAL((DD(1,ifre,imode,iq)*CONJG(DD(2,ifre,imode,iq))))
          endif
        enddo
      enddo
    enddo 
    
  end subroutine get_dEa2_dQ2
  
  
  !===============================================!
  != add bath effect to coordinate and velocitie =!
  !===============================================!
  != ref: notebook page 462 and 638              =!
  !===============================================!
  ! ref: 1 D. M. F. M. Germana Paterlini, Chemical Physics 236 (1998) 243.
  SUBROUTINE ADD_BATH_EFFECT(nmodes,nq,wf,ld_gamma,temp,dEa2_dQ2,TT,l_ph_quantum,XX,VV)
    use kinds,only : dp,dpc
    use parameters,only : lit_ephonon,eps_acustic
    use elph2,only :  iminusq
    use randoms,only : GAUSSIAN_RANDOM_NUMBER_FAST
    use constants,only : KB=>K_B_Ryd,sqrt3,sqrt5,sqrt7,ryd2mev,cone,ci,tpi
    implicit none
    
    integer , intent(in)      :: nq,nmodes
    real(kind=dp), intent(in) :: wf(nmodes,nq)
    real(kind=dp), intent(in) :: ld_gamma(nmodes,nq)
    real(kind=dp), intent(in) :: temp
    real(kind=dp), intent(in) :: dEa2_dQ2(nmodes,nq)
    real(kind=dp), intent(in) :: tt
    logical,intent(in)        :: l_ph_quantum
    complex(kind=dpc), intent(inout) :: XX(nmodes,nq),VV(nmodes,nq)
    
    integer :: imode,iq,i,iq_
    real(kind=dp) :: SIGMAR
    complex(kind=dpc) :: R1,R2,R3,R4,Z1,Z2,Z3,Z4
    complex(kind=dpc) :: cplx_tmp,cplx_tmp_
    real(kind=dp) :: wwf2,wwf2_
    real(kind=dp) :: gamma,womiga,aver_E_T
    real(kind=dp) :: gamma_,womiga_,aver_E_T_
    
    
    !call test_xv(nmodes,nq,xx,vv)
    !  ph_Q(v,q)=ph_Q(v,-q)*
    !  ph_P(v,q)=ph_P(v,-q)*
    DO iq=1,nq
      iq_ = iminusq(iq)
      if( iq <= iq_ ) then
      ! ph_Q(v,q)=ph_Q(v,-q)*
        do imode=1,nmodes
          womiga = wf(imode,iq)
          womiga_= wf(imode,iq_)
          gamma = ld_gamma(imode,iq)
          gamma_= ld_gamma(imode,iq_)
          if(womiga /= womiga_) then
            write(*,*) "wf(",imode,",",iq,") /= wf(",imode,",",iq_,")"
          endif
          if(gamma /= gamma_) then
            write(*,*) "ld_gamma(",imode,",",iq,") /= ld_gamma(",imode,",",iq_,")"
          endif
          if(womiga > eps_acustic ) then
            if(l_ph_quantum) then
              aver_E_T = (bolziman(womiga,temp)+0.5)*womiga
            else
              aver_E_T = KB*TEMP
            endif

            SIGMAR=DSQRT(2.0*gamma*aver_E_T*TT)
            wwf2 = wf(imode,iq)**2+dEa2_dQ2(imode,iq)
            wwf2_= wf(imode,iq_)**2+dEa2_dQ2(imode,iq_)
            
            cplx_tmp = VV(imode,iq)/ABS(VV(imode,iq))            
            cplx_tmp_= VV(imode,iq_)/ABS(VV(imode,iq_))
            !if(cplx_tmp_ /= CONJG(cplx_tmp)) then
            !  write(*,*) "phP(",imode,",",iq,") /= CONJG:phP(",imode,",",iq_,")"
            !  write(*,*) "phP(",imode,",",iq,") =",VV(imode,iq)
            !  write(*,*) "phP(",imode,",",iq_,") =",VV(imode,iq_)
            !endif
          
            R1=GAUSSIAN_RANDOM_NUMBER_FAST(0.0D0,SIGMAR)*cplx_tmp
            R2=GAUSSIAN_RANDOM_NUMBER_FAST(0.0D0,SIGMAR)*cplx_tmp
            R3=GAUSSIAN_RANDOM_NUMBER_FAST(0.0D0,SIGMAR)*cplx_tmp
            R4=GAUSSIAN_RANDOM_NUMBER_FAST(0.0D0,SIGMAR)*cplx_tmp
            Z1=R1    ! V=phP
            Z2=TT*(R1/2.0D0+R2/SQRT3/2.0D0)  !V*T
            Z3=TT**2*(R1/6.0D0+R2*SQRT3/12.0D0+R3/SQRT5/12.0D0) ! V*T**2
            Z4=TT**3*(R1/24.0D0+R2*SQRT3/40.0D0+R3/SQRT5/24.0D0+R4/SQRT7/120.0D0) ! V*T**3
            XX(imode,iq)=XX(imode,iq)+(Z2-GAMMA*Z3+(-wwf2+GAMMA**2)*Z4)
            VV(imode,iq)=VV(imode,iq)+(Z1-GAMMA*Z2+(-wwf2+GAMMA**2)*Z3+(2.0*gamma*wwf2-GAMMA**3)*Z4)
            
            if(iq< iq_) then
              R1= CONJG(R1)
              R2= CONJG(R2)
              R3= CONJG(R3)
              R4= CONJG(R4)
              Z1=R1    ! V=phP
              Z2=TT*(R1/2.0D0+R2/SQRT3/2.0D0)  !V*T
              Z3=TT**2*(R1/6.0D0+R2*SQRT3/12.0D0+R3/SQRT5/12.0D0) ! V*T**2
              Z4=TT**3*(R1/24.0D0+R2*SQRT3/40.0D0+R3/SQRT5/24.0D0+R4/SQRT7/120.0D0) ! V*T**3
              XX(imode,iq_)=XX(imode,iq_)+(Z2-GAMMA_*Z3+(-wwf2_+GAMMA_**2)*Z4)
              VV(imode,iq_)=VV(imode,iq_)+(Z1-GAMMA_*Z2+(-wwf2_+GAMMA_**2)*Z3+(2.0*gamma_*wwf2_-GAMMA_**3)*Z4)            
            
            endif
            !XX(imode,iminusq(iq)) = CONJG(XX(imode,iq)) 
            !VV(imode,iminusq(iq)) = CONJG(VV(imode,iq))
          
          endif
        enddo
      endif
    enddo
    
    !call test_xv(nmodes,nq,xx,vv)
    
  end subroutine add_bath_effect

  
  subroutine pre_md(nmodes,nqtotf,wf,ld_gamma,temp,phQ,phP,l_ph_quantum,pre_dt)
    use surfacecom,only : dEa_dQ,dEa2_dQ2,pre_nstep
    use constants,only  : ry_to_fs,ryd2eV
    use io,only : stdout
    implicit none
    integer ,intent(in) :: nmodes,nqtotf
    real(kind=dp),intent(in) :: pre_dt,temp
    logical,intent(in) :: l_ph_quantum
    real(kind=dp),intent(in) :: wf(nmodes,nqtotf),ld_gamma(nmodes,nqtotf)
    complex(kind=dpc),intent(inout) :: phQ(nmodes,nqtotf),phP(nmodes,nqtotf)


    integer :: istep
    real(kind=dp) :: time
    character(len=2) :: ctimeunit
    
    !call test_wqv(nmodes,nqtotf,wf)
    
    !call test_xv(nmodes,nqtotf,phQ,phP)
    
    dEa_dQ = 0.0
    dEa2_dQ2 = 0.0
    do istep=1,pre_nstep
      call rk4_nuclei(nmodes,nqtotf,dEa_dQ,ld_gamma,wf,phQ,phP,pre_dt)
      call add_bath_effect(nmodes,nqtotf,wf,ld_gamma,temp,dEa2_dQ2,pre_dt,l_ph_quantum,phQ,phP)
    enddo
    !call test_xv(nmodes,nqtotf,phQ,phP)
  
    time = pre_nstep*pre_dt*ry_to_fs
    if(time<1.0E3) then
      ctimeunit = 'fs'
    elseif(time<1.0E6) then
      time = time/1.0E3
      ctimeunit = 'ps'
    elseif(time<1.0E9) then
      time = time/1.0E6
      ctimeunit = 'ns'
    endif

#if defined __MPI 		
    if(ionode) then
#endif     
    write(stdout,"(5X,A23,F6.2,A2,A19,F11.5,A4,A9,F11.5,A4)") &
    "Energy of phonon after ", time,ctimeunit," dynamica: SUM_phK=",0.5*SUM(ABS(phP)**2)*ryd2eV," eV",&
    " SUM_phU=",0.5*SUM(wf**2*ABS(phQ)**2)*ryd2eV," eV"
    write(stdout,"(5X,A23,F6.2,A2,A19,F11.5,A4)") &
    "Energy of phonon after ", time,ctimeunit," dynamica: SUM_phE="&
    ,0.5*SUM(ABS(phP)**2+wf**2*ABS(phQ)**2)*ryd2eV," eV."    
#if defined __MPI 		
    endif
		if(verbosity  == "high") then
    write(procout,"(5X,A23,F6.2,A2,A19,F11.5,A4,A9,F11.5,A4)") &
    "Energy of phonon after ", time,ctimeunit," dynamica: SUM_phK=",0.5*SUM(ABS(phP)**2)*ryd2eV," eV",&
    " SUM_phU=",0.5*SUM(wf**2*ABS(phQ)**2)*ryd2eV," eV"
    write(procout,"(5X,A23,F6.2,A2,A19,F11.5,A4)") &
    "Energy of phonon after ", time,ctimeunit," dynamica: SUM_phE="&
    ,0.5*SUM(ABS(phP)**2+wf**2*ABS(phQ)**2)*ryd2eV," eV."    
		endif
#endif 
  
  end subroutine pre_md

  
  !ref : 1 G. GRIMvall, <The electron-phonon interaction in metals by Goran Grimvall (z-lib.org).pdf> 1981),  
  !    : (3.24)  
  function bolziman(womiga,temp)
    use io,only :stdout
    use constants,only : K_B_Ryd
    implicit none
    real(kind=dp),intent(in)::womiga,temp
    real(kind=dp) :: bolziman
    if(womiga == 0.0) then
      write(stdout,*) "womiga == 0.0,phonon error"
      stop
    endif
    bolziman=1.0/(exp(womiga/(K_B_Ryd*temp))-1.0)
    !<nb>=1/(exp{hbar*w/kbT}-1)
  end function     
  
  subroutine test_wqv(nmodes,nq,wf)
    use elph2,only : iminusq
    implicit none
    integer ,intent(in) :: nmodes,nq
    real(kind=dp),intent(in) :: wf(nmodes,nq)
    
    integer :: iq,iq_
    integer :: imode
    
    do iq=1,nq
      iq_= iminusq(iq)
      do imode=1,nmodes
        if(wf(imode,iq_) /= wf(imode,iq)) then
          write(*,*) "Omiga(n,-q) /= Omiga(n,k)"
        endif
      enddo
    enddo
      
  end subroutine test_wqv
  
  subroutine test_xv(nmodes,nq,phQ,phP)
    use elph2,only :  iminusq
    implicit none
    integer, intent(in) :: nmodes,nq
    complex(kind=dpc), intent(in) :: phQ(nmodes,nq),phP(nmodes,nq)
    
    logical :: lQP_CONJG
    integer :: imode,iq,iq_
    
    lQP_CONJG = .true.
    do iq=1,nq
      iq_ = iminusq(iq)
      do imode=1,nmodes
        if(ABS(REAL(phQ(imode,iq_))/REAL(phQ(imode,iq)) - 1.0) > 1.0E-7 .or. &
           ABS(IMAG(phQ(imode,iq_))/IMAG(phQ(imode,iq)) + 1.0) > 1.0E-7) then
        !if(phQ(imode,iq_) /= CONJG(phQ(imode,iq))) then
          write(*,*) "phQ(",imode,",",iq,") /= CONJG(phQ(",imode,",",iq_,"))"
          write(*,*) "phQ(",imode,",",iq,") =",phQ(imode,iq)
          write(*,*) "phQ(",imode,",",iq_,") =",phQ(imode,iq_)
        endif
        if(ABS(REAL(phP(imode,iq_))/REAL(phP(imode,iq)) - 1.0) > 1.0E-7 .or. &
           ABS(IMAG(phP(imode,iq_))/IMAG(phP(imode,iq)) + 1.0) > 1.0E-7) then      
        !if(phP(imode,iq_) /= CONJG(phP(imode,iq))) then
          write(*,*) "phP(",imode,",",iq,") /= CONJG(phP(",imode,",",iq_,"))"
          write(*,*) "phP(",imode,",",iq,") =",phP(imode,iq)
          write(*,*) "phP(",imode,",",iq_,") =",phP(imode,iq_)
        endif        
      enddo
    enddo
    
  end subroutine test_xv

  subroutine test_v(nmodes,nq,phP)
    use elph2,only :  iminusq
    implicit none
    integer, intent(in) :: nmodes,nq
    complex(kind=dpc), intent(in) ::phP(nmodes,nq)
    
    logical :: lQP_CONJG
    integer :: imode,iq,iq_
    
    lQP_CONJG = .true.
    do iq=1,nq
      iq_ = iminusq(iq)
      do imode=1,nmodes
        if(ABS(REAL(phP(imode,iq_))/REAL(phP(imode,iq)) - 1.0) > 1.0E-7 .or. &
           ABS(IMAG(phP(imode,iq_))/IMAG(phP(imode,iq)) + 1.0) > 1.0E-7) then      
        !if(phP(imode,iq_) /= CONJG(phP(imode,iq))) then
          write(*,*) "phP(",imode,",",iq,") /= CONJG(phP(",imode,",",iq_,"))"
          write(*,*) "phP(",imode,",",iq,") =",phP(imode,iq)
          write(*,*) "phP(",imode,",",iq_,") =",phP(imode,iq_)
        endif        
      enddo
    enddo
    
  end subroutine test_v

 
end module dynamics