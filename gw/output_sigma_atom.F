      subroutine output_sigma_atom
	use atom_mod
	use manager_mod
	use models_mod
	use parallel_mod
	use solid_mod
	use units_mod
      implicit none
      character*2 :: ch
	integer :: i_len,ind_omega,i_omega,l,i,ii,nli,li,ntl
      real*8 :: fc0
	complex*16, allocatable :: sig_out(:,:,:)
	nli=nrel*maxb+1
	fc0=evolt/2.d0
	allocate(sig_out(0:n_omega,maxntle,nli))
	sig_out=(0.d0,0.d0)
      do li=1,nrel*lmb(1)+1
	  if(irel.le.1) then
	    l=li-1
	  else
	    l=li/2
	    i=li-2*l
	    if(i.eq.0) i=-1
	  endif
        ntl=ntle(l,1)
	  do i=1,ntl
	    do ind_omega=1,ndim3_omega
	      i_omega=me3_tau*ndim3_omega+ind_omega-1
	      sig_out(i_omega,i,li)=sigc_omega_atom(i,i,li,ind_omega,1)
	    enddo
	  enddo
	enddo
c --------- Collect from all processes  --------------------------------
	if(nproc_tau/=1) then
	  call dgop(sig_out,2*maxntle*nli*(n_omega+1),'  +',comm_pnk)
	endif
	if(maswrk) then
	  i_len=len_trim(allfile)
        do li=1,nrel*lmb(1)+1
	    if(irel.le.1) then
	      l=li-1
	    else
	      l=li/2
	      i=li-2*l
	      if(i.eq.0) i=-1
	    endif
          ntl=ntle(l,1)
          if(l==0) ch='_S'
          if(l==1) ch='_P'
          if(l==2) ch='_D'
          if(l==3) ch='_F'
          if(l==4) ch='_G'
          if(l>4) cycle
	    if(ubi==' gw') then
	      open(3,file=allfile(1:i_len)//'_Im_Sigma_C'//ch//'.gw')
	      open(4,file=allfile(1:i_len)//'_Re_Sigma_C'//ch//'.gw')
	    else if(ubi==' qp') then
	      open(3,file=allfile(1:i_len)//'_Im_Sigma_C'//ch//'.qp')
	      open(4,file=allfile(1:i_len)//'_Re_Sigma_C'//ch//'.qp')
	    else if(ubi=='gwg') then
	      open(3,file=allfile(1:i_len)//'_Im_Sigma_C'//ch//'.gwg')
	      open(4,file=allfile(1:i_len)//'_Re_Sigma_C'//ch//'.gwg')
	    endif
	    do i_omega=0,n_omega
	      write(3,'(f11.4,5(1x,f12.7))')w_omega(i_omega)*fc0,
     &		               (imag(sig_out(i_omega,ii,li))*fc0,ii=1,ntl)
	      write(4,'(f11.4,5(1x,f12.7))')w_omega(i_omega)*fc0,
     &		               (real(sig_out(i_omega,ii,li))*fc0,ii=1,ntl)
	    enddo
	    close(3)
	    close(4)
	  enddo
	endif
	deallocate(sig_out)
      end