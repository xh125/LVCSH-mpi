#define __MPI
module getwcvk
#if defined __MPI 		
  use global_mpi
  use mp
#endif
  use kinds,only    : dp,dpc
  use parameters,only : verbosity
	use lasercom,only : W_cvk,efield_cart
  implicit none
  
  contains
  !ref : 1 S. Butscher et al., Physical Review B 72 (2005) 
  !ref : 2 <固体物理> (9-29)(9-31)
  subroutine get_Wcvk(ihband_min,ieband_max,fwhm,w_center)
    !得到光激发下垂直跃迁的跃迁几率
    use elph2,only  : vmef,nktotf  !vmef(3,nbndsub,nbndsub,nktotf)
    use readepw,only : etf,icbm
    use io,only : stdout
    use constants,only : ryd2eV,ry_to_fs
    implicit none
    integer , intent(in) :: ihband_min,ieband_max
    real(kind=dp),intent(in) :: fwhm
    real(kind=dp),intent(in) :: w_center
    
    real(kind=dp) :: fwhm_2T2
    
    ! W_cvk(cband,vband,ik)=|<E dot vmef(3,cband,vband,ik)>|^2 *f_w(w)
    !!光激发下的跃迁几率大小
    real(kind=dp) :: fcw
    complex(kind=dpc) :: Evmef
    real(kind=dp) :: E_mnk   !垂直激发能量
    integer :: ik,ikk,ibnd,jbnd,ipol
    integer :: ierr
    integer :: ivbm
    
    ivbm = icbm-1
    
    allocate(W_cvk(icbm:ieband_max,ihband_min:ivbm,nktotf),stat=ierr)
    if(ierr /=0) call errore('getmcvk','Error allocating W_cvk',1)
    
    !ref : 1 S. Butscher et al., Physical Review B 72 (2005) 
    !ref : 1 S. Fernandez-Alberti et al., The Journal of Chemical Physics 137 (2012) 
    fwhm_2T2 = fwhm**2.0/4.0*log(2.0)
    W_cvk = 0.0
    do ik=1,nktotf
      do ibnd=ihband_min,ivbm
        do jbnd=icbm,ieband_max
          Evmef = 0.0
          E_mnk = etf(jbnd,ik)-etf(ibnd,ik)
          do ipol=1,3
            Evmef =Evmef+ (efield_cart(ipol)*(vmef(ipol,jbnd,ibnd,ik)))
          enddo
          fcw = f_w(E_mnk,w_center,fwhm_2T2)
          W_cvk(jbnd,ibnd,ik) = fcw*(ABS(Evmef))**2
        enddo
      enddo
    enddo  
    
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
    !% Write laser information            %!
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
#if defined __MPI 		
    if(ionode) then
#endif    
    write(stdout,"(/,5X,A)") "In the laser obsorbtion,the Pump laser as follow:"
    write(stdout,"(5X,A22,F12.7,A4)")  "Laser centred energy :",w_center*ryd2eV," eV."
    write(stdout,"(5X,A38,F12.7,A4)")  "The full width at half-maximum:fwhm = ",fwhm*ry_to_fs," fs."    
#if defined __MPI 		
    endif
		if(verbosity  == "high") then		
    write(procout,"(/,5X,A)") "In the laser obsorbtion,the Pump laser as follow:"
    write(procout,"(5X,A22,F12.7,A4)")  "Laser centred energy :",w_center*ryd2eV," eV."
    write(procout,"(5X,A38,F12.7,A4)")  "The full width at half-maximum:fwhm = ",fwhm*ry_to_fs," fs."     
		endif
#endif    
    
  end subroutine get_Wcvk

	subroutine get_mij(nefre_sh,nhfre_sh,vij,fwhm,w_center,E_e,E_h,Mij)
		implicit none
		integer, intent(in) :: nefre_sh,nhfre_sh
		complex(kind=dpc),intent(in) :: vij(3,nefre_sh,nhfre_sh)
    real(kind=dp),intent(in) :: fwhm
    real(kind=dp),intent(in) :: w_center
		real(kind=dp),intent(in) :: E_e(nefre_sh),E_h(nhfre_sh)
		real(kind=dp),intent(out):: Mij(nefre_sh,nhfre_sh)

		real(kind=dp) :: fwhm_2T2, dE, fcw
		complex(kind=dpc) :: Evij
		integer :: i,j,ipol
		! Mij(ielec,jhole)=|<E dot vij(3,ielec,jhole)>|^2 *f_w(w)
    !!光激发下在Frocon-Condon windows下的跃迁几率大小
		
    fwhm_2T2 = fwhm**2.0/4.0*log(2.0)
    Mij = 0.0
    
    do j=1,nhfre_sh
      do i=1,nefre_sh
				dE = E_e(i)+E_h(j)
        Evij = 0.0
        do ipol=1,3
          Evij =Evij+ (efield_cart(ipol)*(vij(ipol,i,j)))
        enddo
        fcw = f_w(dE,w_center,fwhm_2T2)
        Mij(i,j) = fcw*(ABS(Evij))**2
      enddo
    enddo 			
		
	end subroutine
  
  !ref : 1 S. Fernandez-Alberti et al., The Journal of Chemical Physics 137 (2012) 
  real function f_w(w,w_laser,fwhm_2T2)
    implicit none
    real(kind=dp),intent(in) :: w,w_laser
    real(kind=dp),intent(in) :: fwhm_2T2
    real(kind=dp) :: wfwhm
    wfwhm = (w-w_laser)**2 *fwhm_2T2
    wfwhm = -0.25*wfwhm
    f_w = EXP(wfwhm)
    !f_w = exp(-1.0*(w-w_laser)**2 * fwhm_2T2/4.0)
    return
  end function

  ! <i|v|j>  where  |i> is the hole adiabatic state and |j> is the electron adiabatic state.
  ! reference : 1 S. Fernandez-Alberti et al., J. Chem. Phys. 137 (2012) 		
	subroutine get_vij(ne,mh,vij,neband,nhband,P_e,P_h)
		use elph2,only  : vmef,nktotf  !vmef(3,nbndsub,nbndsub,nktotf)
		use readepw,only : icbm
		use constants,only : czero
		implicit none
		integer,intent(in) :: ne,mh,neband,nhband
		complex(kind=dpc),intent(out) :: vij(3,ne,mh)
		complex(kind=dpc),intent(in) :: P_e(neband,nktotf,ne),P_h(nhband,nktotf,mh)
		
		integer :: i,j
		integer :: n,m,ik,n_,m_
		integer :: ivbm
		
		ivbm = icbm-1
		
		do i=1,mh
			do j=1,ne
				vij(:,j,i) = czero
				do ik=1,nktotf
					do m=1,nhband
						m_ = m+ivbm-nhband
						do n=1,neband
							n_ = n+icbm-1
							vij(:,j,i) = vij(:,j,i)+P_h(m,ik,i)*CONJG(P_e(n,ik,j))*vmef(:,n_,m_,ik)
						enddo
					enddo
				enddo
			enddo
		enddo
		
	end subroutine get_vij
	
  
end module getwcvk