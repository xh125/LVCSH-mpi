#define __MPI
module surfacehopping
  use omp_lib
  use kinds, only : dp,dpc
  use constants,only : cone,czero,ci
  use epwcom,only : nkf1,nkf2,nkf3,nqf1,nqf2,nqf3,kqmap
  use elph2,only  : wf,nktotf,nbndfst,ibndmin,ibndmax
  use hamiltonian,only : nphfre,neband,nhband,nefre,nhfre,&
                         E_e,P_e,P_e_nk,E0_e,P0_e,P0_e_nk,&
                         E_h,P_h,P_h_nk,E0_h,P0_h,P0_h_nk
  use parameters, only : nsnap,naver,ncore,nnode
  use surfacecom, only : iesurface,ihsurface,esurface_type,hsurface_type,&
                         phQ,phP,phQ0,phP0,phK,phU,SUM_phU,SUM_phK,SUM_phK0,SUM_phE,&
                         phQsit,phPsit,phKsit,phUsit,ld_gamma,nqv,&
                         dEa_dQ,dEa_dQ_e,dEa_dQ_h,dEa2_dQ2,dEa2_dQ2_e,dEa2_dQ2_h,&
                         d_e,g_e,g1_e,c_e,c_e_nk,w_e,w0_e,d0_e,nefre_sh,&
                         d_h,g_h,g1_h,c_h,c_h_nk,w_h,w0_h,d0_h,nhfre_sh,&
                         pes_one_e,pes_e,csit_e,wsit_e,psit_e,&
                         pes_one_h,pes_h,csit_h,wsit_h,psit_h 
  use cc_fssh,only : S_ai_e,S_ai_h,S_bi_e,S_bi_h
	use memory_report,only : MB,GB,complex_size, real_size,int_size,ram,print_memory   
	
  implicit none
  
  complex(kind=dpc),allocatable :: cc0_e(:),dc1_e(:),dc2_e(:),dc3_e(:),dc4_e(:)
  complex(kind=dpc),allocatable :: cc0_h(:),dc1_h(:),dc2_h(:),dc3_h(:),dc4_h(:)
  real(kind=dp)    ,allocatable :: n_e(:),n_h(:)
  
  contains
  
  subroutine allocatesh(methodsh,lelecsh,lholesh,nmodes,nq)
		use io,only : msg
    implicit none
    character(len=*),intent(in) :: methodsh
    logical,intent(in):: lelecsh,lholesh
    integer,intent(in):: nmodes,nq
    integer :: ierr
    
    
    allocate(ld_gamma(nmodes,nq))
    ld_gamma = 0.0
    allocate(nqv(nmodes,nq))
    nqv = 0.0
    allocate(phQ(nmodes,nq),stat=ierr,errmsg=msg)  ! x(1:nphfre)
    if(ierr /=0) then
			call errore('surfacehopping','Error allocating phQ',1)
		endif
		phQ = 0.0
    allocate(phKsit(nmodes,nq,0:nsnap),stat=ierr,errmsg=msg)
		if(ierr /= 0 ) call errore("surfacehopping",trim(msg),1)
    allocate(phUsit(nmodes,nq,0:nsnap),stat=ierr,errmsg=msg)
		if(ierr /= 0 ) call errore("surfacehopping",trim(msg),1)
    phKsit = 0.0d0
    phUsit = 0.0d0
    ram = real_size*nmodes*nq*(1+nsnap)
		call print_memory("phKsit",ram)
		call print_memory("phUsit",ram)
		
		
    allocate(phP(nmodes,nq),stat=ierr,errmsg=msg)  ! v(1:nphfre)
    if(ierr /=0) call errore('surfacehopping','Error allocating phP',1)    
    phP = czero
    !phonons normal mode verlosity
    allocate(phQ0(nmodes,nq),stat=ierr,errmsg=msg)
    if(ierr /=0) call errore('surfacehopping','Error allocating phQ0',1)
		phQ0 = czero
    allocate(phP0(nmodes,nq),stat=ierr,errmsg=msg)
    if(ierr /=0) call errore('surfacehopping','Error allocating phP0',1)
		phP0 = czero
    allocate(phU(nmodes,nq),stat=ierr,errmsg=msg)
    if(ierr /=0) call errore('surfacehopping','Error allocating phU',1)
		phU = 0.0
    allocate(phK(nmodes,nq),stat=ierr,errmsg=msg)
    if(ierr /=0) call errore('surfacehopping','Error allocating phK',1)
		phK = 0.0
    allocate(dEa_dQ(nmodes,nq),stat=ierr,errmsg=msg)
		if(ierr /=0) call errore("surfacehopping",trim(msg),1)
		dEa_dQ = czero
    allocate(dEa2_dQ2(nmodes,nq),stat=ierr,errmsg=msg)
    if(ierr /=0) call errore("surfacehopping",trim(msg),1)
		dEa2_dQ2 = 0.0
		
    if(lelecsh) then
      allocate(E_e(1:nefre),stat=ierr,errmsg=msg)
      if(ierr /=0) call errore('surfacehopping','Error allocating E_e',1)
      allocate(P_e(nefre,nefre),stat=ierr,errmsg=msg)
      if(ierr /=0) call errore('surfacehopping','Error allocating P_e',1)
      allocate(P_e_nk(neband,nktotf,nefre),stat=ierr,errmsg=msg)
      if(ierr /=0) call errore('surfacehopping','Error allocating P_e_nk',1)
      allocate(P0_e_nk(neband,nktotf,nefre),stat=ierr,errmsg=msg)
      if(ierr /=0) call errore('surfacehopping','Error allocating P0_e_nk',1)			
      allocate(d_e(2,nefre_sh,nmodes,nq),stat=ierr,errmsg=msg) !d_ajk,d_jak
      if(ierr /=0) call errore('surfacehopping','Error allocating d_e',1)
			ram=real_size*nmodes*nq*nefre_sh*2
			call print_memory("d_e",ram)
      allocate(dEa_dQ_e(nmodes,nq))
      allocate(dEa2_dQ2_e(nmodes,nq))
      
      allocate(c_e_nk(neband,nktotf),stat=ierr,errmsg=msg)
      if(ierr /=0) call errore('surfacehopping','Error allocating c_e_nk',1)
      allocate(c_e(neband*nktotf),stat=ierr,errmsg=msg)
      if(ierr /=0) call errore('surfacehopping','Error allocating c_e',1)      
      
      allocate(pes_one_e(0:nefre,0:nsnap),pes_e(0:nefre,0:nsnap))
			ram = real_size*(1+nefre)*(1+nsnap)
			call print_memory("pes_one_e",ram)
			call print_memory("pes_e",ram)
      allocate(csit_e(nefre,0:nsnap))
      allocate(wsit_e(nefre,0:nsnap))
      allocate(psit_e(nefre,0:nsnap))
			ram = real_size*nefre*(1+nsnap)
			call print_memory("csit_e",ram)
			call print_memory("wsit_e",ram)
			call print_memory("psit_e",ram)
      pes_one_e = 0.0d0
			pes_e  = 0.0d0
      csit_e = 0.0d0
      wsit_e = 0.0d0
      psit_e = 0.0d0
			

      allocate(cc0_e(nefre),stat=ierr,errmsg=msg)
      if(ierr /=0) call errore('surfacehopping','Error allocating cc0_e',1)
      allocate(dc1_e(nefre),stat=ierr,errmsg=msg)
      if(ierr /=0) call errore('surfacehopping','Error allocating dc1_e',1)
      allocate(dc2_e(nefre),stat=ierr,errmsg=msg)
      if(ierr /=0) call errore('surfacehopping','Error allocating dc2_e',1)
      allocate(dc3_e(nefre),stat=ierr,errmsg=msg)
      if(ierr /=0) call errore('surfacehopping','Error allocating dc3_e',1)
      allocate(dc4_e(nefre),stat=ierr,errmsg=msg)
      if(ierr /=0) call errore('surfacehopping','Error allocating dc4_e',1)
      allocate(n_e(nefre),stat=ierr,errmsg=msg)
      if(ierr /=0) call errore('surfacehopping','Error allocating n_e',1)            
        
      
      allocate(w_e(nefre),stat=ierr,errmsg=msg)
      if(ierr /=0) call errore('surfacehopping','Error allocating w_e',1)  
      allocate(w0_e(nefre),stat=ierr,errmsg=msg)
      if(ierr /=0) call errore('surfacehopping','Error allocating w0_e',1)
      w_e = czero
			w0_e= czero
			
      allocate(g_e(1:nefre_sh),stat=ierr,errmsg=msg)  !g_ij
      if(ierr /=0) call errore('surfacehopping','Error allocating g_e',1)
      allocate(g1_e(1:nefre_sh),stat=ierr,errmsg=msg)
      if(ierr /=0) call errore('surfacehopping','Error allocating g1_e',1) 
      
      allocate(E0_e(1:nefre),stat=ierr,errmsg=msg)
      if(ierr /=0) call errore('surfacehopping','Error allocating E0_e',1)
      allocate(P0_e(nefre,nefre),stat=ierr,errmsg=msg)
      if(ierr /=0) call errore('surfacehopping','Error allocating P0_e',1)
      allocate(d0_e(2,nefre_sh,nmodes,nq),stat=ierr,errmsg=msg)
      if(ierr /=0) call errore('surfacehopping','Error allocating d0_e',1)
      
      if(methodsh == "CC-FSSH") then
        allocate(S_ai_e(nefre_sh))
        allocate(S_bi_e(nefre_sh))
      endif
      
    endif
    
    if(lholesh) then
      allocate(E_h(1:nhfre),stat=ierr,errmsg=msg)
      if(ierr /=0) call errore('surfacehopping','Error allocating E_h',1)
      allocate(P_h(nhfre,nhfre),stat=ierr,errmsg=msg)
      if(ierr /=0) call errore('surfacehopping','Error allocating P_h',1)
      allocate(P_h_nk(nhband,nktotf,nhfre),stat=ierr,errmsg=msg)
      if(ierr /=0) call errore('surfacehopping','Error allocating P_h_nk',1)
      allocate(P0_h_nk(nhband,nktotf,nhfre),stat=ierr,errmsg=msg)
      if(ierr /=0) call errore('surfacehopping','Error allocating P0_h_nk',1)			
      allocate(d_h(2,nhfre_sh,nmodes,nq),stat=ierr,errmsg=msg) !d_ajk
      if(ierr /=0) call errore('surfacehopping','Error allocating d_h',1)
			ram=real_size*nmodes*nq*nhfre_sh
			call print_memory("d_h",ram)
      allocate(dEa_dQ_h(nmodes,nq))
      allocate(dEa2_dQ2_h(nmodes,nq))

      allocate(c_h_nk(nhband,nktotf),stat=ierr,errmsg=msg)
      if(ierr /=0) call errore('surfacehopping','Error allocating c_h_nk',1)  
      allocate(c_h(nhband*nktotf),stat=ierr,errmsg=msg)
      if(ierr /=0) call errore('surfacehopping','Error allocating c_h',1)  
      
      allocate(pes_one_h(0:nhfre,0:nsnap),pes_h(0:nhfre,0:nsnap))
			ram = real_size*(1+nhfre)*(1+nsnap)
			call print_memory("pes_one_h",ram)
			call print_memory("pes_h",ram)			
      allocate(csit_h(nhfre,0:nsnap))
      allocate(wsit_h(nhfre,0:nsnap))
      allocate(psit_h(nhfre,0:nsnap))
			ram = real_size*nhfre*(1+nsnap)
			call print_memory("csit_h",ram)
			call print_memory("wsit_h",ram)
			call print_memory("psit_h",ram)			
			pes_one_h=0.0d0
      pes_h  = 0.0
      csit_h = 0.0
      wsit_h = 0.0
      psit_h = 0.0
			
      allocate(cc0_h(nhfre),stat=ierr,errmsg=msg)
      if(ierr /=0) call errore('surfacehopping','Error allocating cc0_h',1)
      allocate(dc1_h(nhfre),stat=ierr,errmsg=msg)
      if(ierr /=0) call errore('surfacehopping','Error allocating dc1_h',1)
      allocate(dc2_h(nhfre),stat=ierr,errmsg=msg)
      if(ierr /=0) call errore('surfacehopping','Error allocating dc2_h',1)
      allocate(dc3_h(nhfre),stat=ierr,errmsg=msg)
      if(ierr /=0) call errore('surfacehopping','Error allocating dc3_h',1)
      allocate(dc4_h(nhfre),stat=ierr,errmsg=msg)
      if(ierr /=0) call errore('surfacehopping','Error allocating dc4_h',1)
      allocate(n_h(nhfre),stat=ierr,errmsg=msg)
      if(ierr /=0) call errore('surfacehopping','Error allocating n_h',1)                  
      
      allocate(w_h(nhfre),stat=ierr,errmsg=msg)
      if(ierr /=0) call errore('surfacehopping','Error allocating w_h',1)  
      allocate(w0_h(nhfre),stat=ierr,errmsg=msg)
      if(ierr /=0) call errore('surfacehopping','Error allocating w0_h',1)
      w_h = czero
			w0_h= czero
			
      allocate(g_h(1:nhfre_sh),stat=ierr,errmsg=msg)  !g_ij
      if(ierr /=0) call errore('surfacehopping','Error allocating g_h',1)
      allocate(g1_h(1:nhfre_sh),stat=ierr,errmsg=msg)
      if(ierr /=0) call errore('surfacehopping','Error allocating g1_h',1) 
      
      allocate(E0_h(1:nhfre),stat=ierr,errmsg=msg)
      if(ierr /=0) call errore('surfacehopping','Error allocating E0_h',1)    
      
      allocate(P0_h(nhfre,nhfre),stat=ierr,errmsg=msg)
      if(ierr /=0) call errore('surfacehopping','Error allocating P0_h',1)
      allocate(d0_h(2,nhfre_sh,nmodes,nq),stat=ierr,errmsg=msg)
      if(ierr /=0) call errore('surfacehopping','Error allocating d0_h',1)
      
      if(methodsh == "CC-FSSH") then
        allocate(S_ai_h(nhfre_sh))
        allocate(S_bi_h(nhfre_sh))
      endif
    endif
    
    
  end subroutine allocatesh
  
  

  
  !=========================================================!
  != convert wavefunction from diabatix to adiabatic basis =!
  !=========================================================!  
  subroutine convert_diabatic_adiabatic(nfre,pp,cc,ww)                                        
    use f95_precision
    use blas95
    implicit none
    integer,intent(in) :: nfre 
    complex(kind=dpc),intent(in) :: pp(nfre,nfre)
    complex(kind=dpc),intent(in) :: cc(nfre)
    complex(kind=dpc),intent(out):: ww(nfre)
    real(kind=dp) :: sum_cc2,sum_ww2
    
    sum_cc2 = REAL(SUM(cc*CONJG(cc)))
    
    ww= czero
    !call gemv(pp,cc,ww)
    call gemv(pp,cc,ww,trans='T')
    
    sum_ww2 = REAL(SUM(ww*CONJG(ww)))
    
    !ww=0.0d0
    !do ik=1,nktotf
    !  do iband=1,nbndfst
    !    iefre = (ik-1)*nbndfst + iband
    !    do jk=1,nktotf
    !      do jband=1,nbndfst
    !        ww(iefre) = ww(iefre)+pp_nk(jband,jk,iefre)*c_nk(jband,jk)
    !      enddo
    !    enddo
    !  enddo
    !enddo
    !
    !do iefre=1,nefre
    !  if(ww(iefre)/=ww_(iefre)) then
    !    lgemv=.false.
    !    exit
    !  endif
    !enddo
      
    
  
  end subroutine convert_diabatic_adiabatic

  !=========================================================!
  != convert wavefunction from adiabatix to diabatic basis =!
  !=========================================================!  
  subroutine convert_adiabatic_diabatic(nfre,pp,ww,cc)                                        
    use f95_precision
    use blas95
    implicit none
    integer,intent(in) :: nfre 
    complex(kind=dpc),intent(in) :: pp(nfre,nfre)
    complex(kind=dpc),intent(in) :: ww(nfre)
    complex(kind=dpc),intent(out):: cc(nfre)
    real(kind=dp) :: sum_cc2,sum_ww2
    
    sum_ww2 = REAL(SUM(ww*CONJG(ww)))
    
    cc= czero
    !call gemv(pp,cc,ww)
    call gemv(pp,ww,cc)
    
    sum_cc2 = REAL(SUM(cc*CONJG(cc)))
  
  end subroutine convert_adiabatic_diabatic


	!ref:1	Bedard-Hearn, M. J., Larsen, R. E. & Schwartz, B. J. Mean-field dynamics with stochastic decoherence (MF-SD):
	!				a new algorithm for nonadiabatic mixed quantum/classical molecular-dynamics simulations with nuclear-induced decoherence.
	!				J Chem Phys 123, 234106, doi:10.1063/1.2131056 (2005).
	!ref:2 1	Granucci, G. & Persico, M. Critical appraisal of the fewest switches algorithm for surface hopping. 
	!				J Chem Phys 126, 134114, doi:10.1063/1.2715585 (2007).
	!ref:3 1	Zhu, C., Nangia, S., Jasper, A. W. & Truhlar, D. G.
	!				Coherent switching with decay of mixing: an improved treatment of electronic coherence for non-Born-Oppenheimer trajectories. 
	!				J Chem Phys 121, 7658-7670, doi:10.1063/1.1793991 (2004).
	!ref:4 1	Qiu, J., Bai, X. & Wang, L. Crossing Classified and Corrected Fewest Switches Surface Hopping.
	!					The Journal of Physical Chemistry Letters 9, 4319-4325, doi:10.1021/acs.jpclett.8b01902 (2018 
	subroutine add_decoherence(C,Ekin,dt,nfre,isurface,E,ww)
		implicit none
		integer ,intent(in) :: nfre,isurface
		real(kind=dp), intent(in) :: C,Ekin,dt
		real(kind=dp), intent(in) :: E(nfre)
		complex(kind=dpc),intent(inout) :: ww(nfre)
		
		integer :: ifre
		real(kind=dp) :: tau
		real(kind=dp) :: flad,factor
		
		tau = 0.0
		flad = 0.0
		do ifre=1,nfre
			if(ifre /= isurface) then
				tau = (1.0 + C/Ekin)/abs(E(ifre)-E(isurface))
				factor = exp(-1.0*dt/tau)
				ww(ifre) = ww(ifre)*exp(-1.0*dt/tau)
				flad = flad + ww(ifre)*CONJG(ww(ifre))
			endif
		enddo
		
		ww(isurface) = ww(isurface)*sqrt(1.0 - flad)/sqrt(ww(isurface)*CONJG(ww(isurface)))
		
	end subroutine add_decoherence
  
  !==================================================================!
  != calculate nonadiabatic coupling                                =!
  !==================================================================!
  ! ref : PPT-91
  ! The most time-consuming part of the program
  subroutine calculate_nonadiabatic_coupling(nmodes,nq,nband,nk,ee,p_nk,gmnvkq,nfre_sh,isurface,dd)
#if defined __MPI 		
  use global_mpi
#endif
    use kinds,only : dp
    use elph2,only : iminusq
    implicit none
    integer, intent(in)          :: nmodes,nq,nband,nk,nfre_sh,isurface
    real(kind=dp),intent(in)     :: ee(nband*nk)
    complex(kind=dpc),intent(in) :: p_nk(nband,nk,nband*nk)
    real(kind=dp),intent(in)     :: gmnvkq(nband,nband,nmodes,nk,nq)    
    complex(kind=dpc),intent(out):: dd(2,nfre_sh,nmodes,nq)
    
		integer :: nfre,ifre,iq,imode
    integer :: ik,ikq,m,n,iq_
		real(kind=dp) :: epc
    complex(kind=dpc) :: A,B
    real(kind=dp) :: C,D
    
    !nfre = nband*nk
    nfre = nfre_sh
    
    dd=czero
!$omp parallel do
		do iq=1,nq
      do ik =1 ,nk
        ikq = kqmap(ik,iq)
        do imode=1,nmodes
          do m=1,nband ! m
            do n=1,nband ! n
              epc = gmnvkq(m,n,imode,ik,iq)
              if(epc /= 0.0 ) then
                do ifre=1,nfre
                  ! dajqv
                  dd(1,ifre,imode,iq) = dd(1,ifre,imode,iq)+epc*CONJG(p_nk(m,ikq,isurface))*p_nk(n,ik,ifre)
                  ! djaqv
                  !dd(2,ifre,imode,iq) = dd(2,ifre,imode,iq)+epc*CONJG(p_nk(m,ikq,ifre))*p_nk(n,ik,isurface)
                enddo
              endif
            enddo
          enddo
        enddo
      enddo
    enddo
!$omp end parallel do


!$omp parallel do    
    do iq=1,nq
      iq_ = iminusq(iq)
      do imode=1,nmodes
        do ifre=1,nfre
          if(ifre/= isurface) then
            dd(1,ifre,imode,iq) = dd(1,ifre,imode,iq)/(ee(ifre)-ee(isurface))
            !  Let dijqv(i,j,-q,v) = -1.0*CONJG(dijqv(j,i,q,v)) (! i/=j)
            dd(2,ifre,imode,iq_) = -1.0*CONJG(dd(1,ifre,imode,iq))
          else
            dd(2,isurface,imode,iq) = dd(1,isurface,imode,iq)
          endif
        enddo
      enddo
		enddo
!$omp end parallel do
    
    ! Let dEa_dQ(imode,iq_) = dEa_dQ(imode,iq)* 
    ! dijqv(i,i,q,v) = CONJG(dijqv(i,i,-q,v))
    do iq=1,nq
      iq_=iminusq(iq)
      if(iq <= iq_) then
      do imode=1,nmodes
        if(iq_ == iq) then
          A = dd(1,isurface,imode,iq)
          B = dd(2,isurface,imode,iq)
          if(ABS(Imag(A)/Real(A)) < 1.0E-7) then
            dd(1,isurface,imode,iq) = Real(A)
          endif
          if(ABS(Imag(B)/Real(B)) < 1.0E-7) then
            dd(2,isurface,imode,iq) = Real(B)
          endif          
        else
          
          A = dd(1,isurface,imode,iq)
          B = dd(1,isurface,imode,iq_)
          if(ABS(Real(A)/REAL(B) - 1.0) < 1.0E-7 .and. ABS(Imag(A)/Imag(B) + 1.0) < 1.0E-7) then
            C = (Real(A) + Real(B))/2.0
            D = (Imag(A) - Imag(B))/2.0
            dd(1,isurface,imode,iq) = C*cone +D*ci
            dd(1,isurface,imode,iq_)= C*cone -D*ci
          endif
          
          A = dd(2,isurface,imode,iq)
          B = dd(2,isurface,imode,iq_)
          if(ABS(Real(A)/REAL(B) - 1.0) < 1.0E-7 .and. ABS(Imag(A)/Imag(B) + 1.0) < 1.0E-7) then
            C = (Real(A) + Real(B))/2.0
            D = (Imag(A) - Imag(B))/2.0
            dd(2,isurface,imode,iq) = C*cone +D*ci
            dd(2,isurface,imode,iq_)= C*cone -D*ci          
          endif
        
        endif
      enddo
      endif
    enddo

    ! Let dijqv(i,j,-q,v) = dijqv(j,i,q,v)*
    call test_dijqv(isurface,nfre,nmodes,nq,dd)
    
    !write(procout,*) "dijqv =",dd
    
  end subroutine calculate_nonadiabatic_coupling
  
  
  function bolziman(womiga,temp)
    use kinds ,only : dp
    implicit none
    real(kind=dp)::womiga,temp,bolziman

    bolziman=1.0/(exp(womiga/(temp))-1.0)
    !<nb>=1/(exp{hw/kbT}-1)
  end function bolziman
 
  !================================================================!
  != CALCULATE HOPPING PROBABILITY                                =! 
  !================================================================!
  !call calculate_hopping_probability(w0_e,phP0,d0,dt,g,g1) using FSSH in adiabatic representation
  !ref : 1 J. C. Tully, J. Chem. Phys. 93 (1990) 1061.
  !ref : 2 J. Qiu, X. Bai, and L. Wang, The Journal of Physical Chemistry Letters 9 (2018) 4319.
  !ref : eq(2)
  ! 其中gg1为原始跃迁几率，gg为采用FSSH方法，令跃迁几率小于0的部分等于0
  ! if(gg(iefre) < 0.0) gg(iefre) = 0.0
  subroutine calculate_hopping_probability(isurface,nfre,nfre_sh,nmodes,nq,WW,VV,dd,tt,gg,gg1)
    implicit none
    integer,intent(in)           :: isurface,nfre,nfre_sh,nmodes,nq
    complex(kind=dpc),intent(in) :: WW(nfre)
    complex(kind=dpc),intent(in) :: VV(nmodes,nq)
    complex(kind=dpc),intent(in) :: dd(2,nfre_sh,nmodes,nq)
    real(kind=dp),intent(in)     :: tt
    real(kind=dp),intent(out)    :: gg(nfre_sh)
    real(kind=dp),intent(out)    :: gg1(nfre_sh)
    
    complex(kind=dpc) :: sumvd
    integer :: ifre,iq,imode
    
    gg = 0.0d0
    gg1= 0.0d0
    ! FSSH
    ! ref: 1 J. Qiu, X. Bai, and L. Wang, The Journal of Physical Chemistry Letters 9 (2018) 4319.
    do ifre=1,nfre_sh
      if(ifre /= isurface) then
        sumvd = czero
        do iq=1,nq
          do imode=1,nmodes
            sumvd = sumvd+VV(imode,iq)*dd(1,ifre,imode,iq)
          enddo
        enddo
        ! in adiabatic representation：the switching probabilities from the active surface isurface to another surface iefre 
        gg(ifre)=2.0*tt*Real((CONJG(WW(isurface))*WW(ifre))*sumvd)/REAL(CONJG(WW(isurface))*WW(isurface))
        gg1(ifre) = gg(ifre)  ! 绝热表象原始的跃迁几率
        !FSSH if g_ij<0,reset to g_ij=0
        if(gg(ifre) < 0.0d0) gg(ifre) = 0.0d0
      endif
    enddo
      
  end subroutine calculate_hopping_probability

  subroutine test_deadqv(nmodes,nq,dEadQ)
    use elph2,only :  iminusq
    implicit none
    integer , intent(in) :: nmodes,nq
    complex(kind=dpc),intent(in) :: dEadQ(nmodes,nq)
    
    complex(kind=dpc) :: A,B
    integer :: imode,iq,iq_
    do iq=1,nq
      iq_=iminusq(iq)
      do imode=1,nmodes
        A = dEadQ(imode,iq)
        B = dEadQ(imode,iq_)
        if(dEadQ(imode,iq) /= CONJG(dEadQ(imode,iq_))) then
          !if(iq_ == iq) then
            !if(ABS(Imag(A)/Real(A)) > 1.0E-7) then
              !write(*,*) "dEadQ(",imode,",",iq,") is not Real."
            !endif
          !else
            !if(ABS(Real(A)/REAL(B) - 1.0) > 1.0E-7 .and. ABS(Imag(A)/Imag(B) + 1.0) > 1.0E-7) then
              write(*,*) "dEadQ(",imode,",",iq,")/=CONJG(dEadQ(",imode,",",iq_,"))"
              write(*,*) "dEadQ(",imode,",",iq,")=",dEadQ(imode,iq)
              write(*,*) "dEadQ(",imode,",",iq_,")=",dEadQ(imode,iq_)
            !endif
          !endif
        endif
      enddo
    enddo
     
  end subroutine test_deadqv

  subroutine test_dijqv(isurface,nfre,nmodes,nq,dd)
    use elph2,only : iminusq
    implicit none
    integer , intent(in) :: isurface,nfre,nmodes,nq
    complex(kind=dpc), intent(in) :: dd(2,nfre,nmodes,nq)
    
    integer :: ifre,imode,iq,iq_
    complex(kind=dpc) :: A,B
    
    ! dijqv(a,i,q,v) = dijqv(i,a,q,v)
    !do iq=1,nq
    !  do imode=1,nmodes
    !    do ifre=1,nfre
    !      if(ifre /= isurface) then
    !        if(dd(1,ifre,imode,iq) /= dd(2,ifre,imode,iq)) then
    !          write(*,*) "dajqv(",ifre,",",imode,",",iq,") /= dajqv(",ifre,",",imode,",",iq,")"
    !          write(*,*) "dajqv(",ifre,",",imode,",",iq,") =",dd(1,ifre,imode,iq)
    !          write(*,*) "djaqv(",ifre,",",imode,",",iq,") =",dd(2,ifre,imode,iq)
    !        endif
    !      endif
    !    enddo
    !  enddo
    !enddo
    
    ! dijqv(a,j,-q,v)=dijqv(a,j,q,v)*
    !do iq=1,nq
    !  iq_ = iminusq(iq)
    !  do imode=1,nmodes
    !    do ifre=1,nfre
    !      if(ifre /= isurface) then
    !        if(dd(1,ifre,imode,iq_) /= CONJG(dd(1,ifre,imode,iq))) then
    !          write(*,*) "dajqv(",ifre,",",imode,",",iq,") /= CONJG(dajqv(",ifre,",",imode,",",iq_,"))"
    !          write(*,*) "dajqv(",ifre,",",imode,",",iq,") =",dd(1,ifre,imode,iq),ABS(dd(1,ifre,imode,iq))
    !          write(*,*) "dajqv(",ifre,",",imode,",",iq_,") =",dd(1,ifre,imode,iq_),ABS(dd(1,ifre,imode,iq_))
    !        endif
    !      endif
    !    enddo
    !  enddo
    !enddo
    
    ! dijqv(i,j,q,v) = -1.0*CONJG(dijqv(j,i,-q,v)) (if i/=j)
    ! dijqv(i,i,q,v) = Real (i==j .and. -q=q+G)
    ! dijqv(i,i,q,v) = CONJG(dijqv(i,i,-q,v)) (i==j .and. -q/=q+G)
    ! 验证成功
     do iq=1,nq
      iq_ = iminusq(iq)
      do imode=1,nmodes
        do ifre=1,nfre
          if(ifre /= isurface) then
            A =  dd(2,ifre,imode,iq_)
            B =  -1.0*CONJG(dd(1,ifre,imode,iq))
            if(ABS(REAL(A)/REAL(B) - 1.0) > 1.0E-7 .or. ABS(IMAG(A)/IMAG(B) - 1.0) > 1.0E-7) then            
            !if(dd(2,ifre,imode,iq_) /= -1.0*CONJG(dd(1,ifre,imode,iq))) then
              write(*,*) "dajqv(",ifre,",",imode,",",iq,") /= -1.0*CONJG(djaqv(",ifre,",",imode,",",iq_,"))"
              write(*,*) "dajqv(",ifre,",",imode,",",iq,") =",dd(1,ifre,imode,iq)
              write(*,*) "djaqv(",ifre,",",imode,",",iq_,") =",dd(2,ifre,imode,iq_)
            endif
          else
            A =  dd(2,ifre,imode,iq_)
            B =  CONJG(dd(1,ifre,imode,iq))
            if(ABS(REAL(A)/REAL(B) - 1.0) > 1.0E-7 .or. ABS(IMAG(A)/IMAG(B) - 1.0) > 1.0E-7) then          
            !if(dd(2,isurface,imode,iq_) /= CONJG(dd(1,isurface,imode,iq))) then
              write(*,*) "daaqv(",ifre,",",imode,",",iq,") /= CONJG(daaqv(",ifre,",",imode,",",iq_,"))"
              write(*,*) "daaqv(",ifre,",",imode,",",iq,") =",dd(1,ifre,imode,iq)
              write(*,*) "daaqv(",ifre,",",imode,",",iq_,") =",dd(2,ifre,imode,iq_)
            endif            
          endif
        enddo
      enddo
    enddo   
    
    
    ! Theta(dijqv) 
    !do iq=1,nq
    !  iq_ = iminusq(iq)
    !  if(iq == iq_) then
    !    do ifre=1,nfre
    !      if(ifre /= isurface) then
    !        do imode=1,nmodes
    !          write(*,*) "Theta(dajqv(",ifre,",",imode,",",iq,"))",dd(1,ifre,imode,iq)/ABS(dd(1,ifre,imode,iq))
    !        enddo
    !      endif
    !    enddo
    !  end if
    !enddo
    
  end subroutine test_dijqv

end module surfacehopping