      subroutine output_dos
	use atom_mod
	use manager_mod
	use parallel_mod
	use solid_mod
	use units_mod
      implicit none
      character*2 :: s1
      character*3 :: s2
      character*5 :: ss1
      character*6 :: sss
	integer :: iedos,ispin,i_len,isort,li,k,i,ind0,lf,l,ndim,ind_dos,
     &	         ks,ie
	real*8 :: deltae,e0,consp(2),om0,fct,ds1,ds2
	integer, allocatable :: nd_dos(:),nmpi_dos(:)
	real*8, allocatable :: weigd(:,:,:),dossp(:,:),tm(:),
     &	                     weigh(:,:,:),pdos(:,:,:,:),
     &                       porb(:,:,:,:),dos_k(:,:,:)
      complex*16, allocatable :: achain(:,:,:,:,:)
	allocate(weigd(nbndf,npnt,nspin),weigh(nbndf,npnt,nspin),
     &	       dossp(nspin,0:ndos))
	dossp=0.d0
	allocate(dos_k(0:ndos,n_sym,nspin))
	allocate(pdos(nrel*maxb+1,nsort,nspin,0:ndos))
	if(irel/=2.and.iprn/=0) allocate(porb(maxel_loc,nsort,nspin,
     &	                                    0:ndos))
	allocate(tm(nbndf))
	consp(1)=-1.d0
	consp(nspin)=1.d0
      deltae=(emaxdos-emindos)/ndos
	fct=evolt/2
	if(ubi==' qp'.or.ubi==' gw') then
	  call output_z_factor
	  call output_sym_k_orb_char(e_qp,z_qp,nbndf)
	else
	  call output_sym_k_orb_char(e_bnd,z_bnd,nbndf)
	endif
c ---------------------------------------------------------------------
      if(ubi==' gw') then
	  allocate(achain(0:n_cont_frac,nbndf,nbndf,maxval(ndim3_k),
     &           nspin))
	  call a_cont_coeff(achain,nbndf)
	endif
	allocate(nd_dos(nproc_tau))
	allocate(nmpi_dos(nproc_tau))
	call size_shift_par(ndos+1,nproc_tau,nd_dos,nmpi_dos)
	dossp=0.d0
	pdos=0.d0
	dos_k=0.d0
	if(irel/=2.and.iprn/=0) porb=0.d0
	do ind_dos=1,nd_dos(me3_tau+1)
	  iedos=nmpi_dos(me3_tau+1)+ind_dos-1
        om0=emindos+iedos*deltae
        if(ubi=='dft'.or.ubi==' hf') then
          e0=om0+chem_pot
        else
          e0=om0+chem_pot_qp
        endif
c -----------------------------------------------------------------
	  if(ubi=='dft'.or.ubi==' hf'.or.ubi==' qp') then
	    call sum_up_weights(weigd,om0,nbndf)
	  else if(ubi==' gw') then
	    call a_cont_weights(achain,weigd,om0,nbndf)
	  endif
c ------ Full and Partial DOS for all variants -------------------	    
	  call dos(dossp(1,iedos),pdos(1,1,1,iedos),weigd,nbndf)
c ------ K-resolved spectral function for SYM-K -------------------
        do ispin=1,nspin
	    do k=1,npnt
	      ks=k_npnt_in_sym(k)
	      if(ks==0) cycle
	      do i=1,n_bnd(k,ispin)
              dos_k(iedos,ks,ispin)=dos_k(iedos,ks,ispin)
     &                             +weigd(i,k,ispin)/wgt(k)
          	enddo
          enddo
        enddo
      enddo
      if(nproc_tau/=1) then
	  call dgop(dossp,nspin*(1+ndos),'  +',comm_pnk)
	  call dgop(pdos,(nrel*maxb+1)*nsort*nspin*(1+ndos),'  +',
     &	          comm_pnk)
	  if(irel/=2.and.iprn/=0) 
     &	  call dgop(porb,maxel_loc*nsort*nspin*(1+ndos),'  +',
     &              comm_pnk)
	  call dgop(dos_k,nspin*n_sym*(1+ndos),'  +',comm_pnk)
	endif
      deallocate(nd_dos,nmpi_dos)
	if(ubi==' gw') deallocate(achain)
c --------- DOS, DOS_K, and PDOS  -------------
	if(maswrk) then
	  i_len=len_trim(allfile)
c ---- Full DOS (spin UP + spin DOWN ) ---------------------------	  
	  if(ubi=='dft') open(3,file=allfile(1:i_len)//'_dft.dos')
	  if(ubi==' hf') open(3,file=allfile(1:i_len)//'_hf.dos')
	  if(ubi=='  x') open(3,file=allfile(1:i_len)//'_x.dos')
	  if(ubi==' qp') open(3,file=allfile(1:i_len)//'_qp.dos')
	  if(ubi==' gw') open(3,file=allfile(1:i_len)//'_gw.dos')
	  if(irel<2) then
	    if(nspin==1) then
            do iedos=0,ndos
              om0=emindos+iedos*deltae
              ds1=dossp(1,iedos)/fct/2
              write(3,'(4(f11.5,2x))')om0*fct,ds1,ds1,ds1+ds1
            enddo
	    else if(nspin==2) then
            do iedos=0,ndos
              om0=emindos+iedos*deltae
              ds1=dossp(1,iedos)/fct/2
              ds2=dossp(2,iedos)/fct/2
              write(3,'(4(f11.5,2x))')om0*fct,ds1,ds2,ds1+ds2
            enddo
          endif
        else if(irel==2) then
          do iedos=0,ndos
            om0=emindos+iedos*deltae
            ds1=dossp(1,iedos)/fct
            write(3,'(2(f11.5,2x))')om0*fct,ds1
          enddo
        endif
	  close(3)
c ---- Full DOS_K (spin UP + spin DOWN ) ---------------------------
        do ispin=1,nspin
          if(ispin==1) s2='_dn'
          if(nspin==2.and.ispin==2) s2='_up'	 	  
	    if(ubi=='dft') open(3,file=allfile(1:i_len)//s2//'_dft.dosk')	
	    if(ubi==' hf') open(3,file=allfile(1:i_len)//s2//'_hf.dosk')
	    if(ubi=='  x') open(3,file=allfile(1:i_len)//s2//'_x.dosk')
	    if(ubi==' qp') open(3,file=allfile(1:i_len)//s2//'_qp.dosk')
	    if(ubi==' gw') open(3,file=allfile(1:i_len)//s2//'_gw.dosk')
		  write(3,'(11x,7(6x,a1,5x))')(sym0_label(i),i=1,min(n_sym,7))
          do iedos=0,ndos
            om0=emindos+iedos*deltae
            write(3,'(f11.5,7(1x,f11.5))')om0*fct,
     &       (dos_k(iedos,i,ispin)/fct,i=1,min(n_sym,7))
          enddo
	    close(3)
        enddo
c ---- Partial DOSs --------------------------------------------	  
          if(ubi=='dft'.or.ubi==' hf'.or.ubi==' qp') then
	      sss='___sum'
          else if(ubi==' gw') then
	      sss='_Chain'
	    endif
	    do isort=1,nsort
c -------- S P D F .... ------------------------------------	      
	      if(ubi=='dft') open(3,file=allfile(1:i_len)
     &		                 //txtel(isort)(1:2)//sss//'_dft.pdos')
	      if(ubi==' hf') open(3,file=allfile(1:i_len)
     &		                 //txtel(isort)(1:2)//sss//'_hf.pdos')
	      if(ubi=='  x') open(3,file=allfile(1:i_len)
     &		                 //txtel(isort)(1:2)//sss//'_x.pdos')
	      if(ubi==' gw') open(3,file=allfile(1:i_len)
     &		                 //txtel(isort)(1:2)//sss//'_gw.pdos')
	      if(ubi==' qp') open(3,file=allfile(1:i_len)
     &		                 //txtel(isort)(1:2)//sss//'_qp.pdos')
            do iedos=0,ndos
              om0=emindos+iedos*deltae
              write(3,'(9f8.3)')om0*fct,
     &			  ((pdos(li,isort,ispin,iedos)*consp(ispin)
     &                                        *nhsort(isort)/fct,
     &                  li=1,nrel*min(3,lmb(isort))+1),ispin=1,nspin)
	      enddo
	      close(3)
	    enddo   !! over isort
	  endif
c     &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
	deallocate(dossp,weigd,weigh,pdos,tm,dos_k)
	if(irel/=2.and.iprn/=0) deallocate(porb)
      call timel('***** OUTPUT_DOS finished **********')
      end
