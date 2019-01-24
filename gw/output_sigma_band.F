      subroutine output_sigma_band(f,ch1,ch2,fc)
	use atom_mod
	use manager_mod
	use parallel_mod
	use solid_mod
	use units_mod
      implicit none
      character*6, intent(in) :: ch1
      character*4, intent(in) :: ch2
      real*8, intent(in) :: f(nbndf,nbndf,2,ndim3_omega,
     &                      maxval(ndim3_k),nspin),fc
	integer :: k,i_len,m,n,ispin,l1,l2,m0,ind_k,ind_omega,i_omega,l,
     &	         ll,i,ii,n0,ks,i7,nnn
      real*8 :: fc0
	real*8, allocatable :: sig_out(:,:,:)
	nnn=8
	fc0=evolt/2.d0
	n0=nbndf
	m=5
	n=n0/m
	if(n*m<n0) n=n+1
	allocate(sig_out(0:n_omega,nbndf,n_sym))
	if(maswrk) then
	  i_len=len_trim(allfile)
	  if(ubi==' gw') then
	    open(3,file=allfile(1:i_len)//'_Im_'//ch1//'_band_k.gw')
	  else if(ubi==' qp') then
	    open(3,file=allfile(1:i_len)//'_Im_'//ch1//'_band_k.qp')
	  else if(ubi=='gwg') then
	    open(3,file=allfile(1:i_len)//'_Im_'//ch1//'_band_k.gwg')
	  endif
	endif
	do ispin=1,nspin
	  sig_out=0.d0
	  do i=1,n
          l1=m*(i-1)+1
          l2=min(m*i,n0)
	    m0=l2-l1+1
	    do ind_k=1,ndim3_k(me3_k+1)
	      k=n3_mpi_k(me3_k+1)+ind_k
	      ks=k_npnt_in_sym(k)
	      if(ks==0) cycle
	      do ind_omega=1,ndim3_omega
	        i_omega=me3_tau*ndim3_omega+ind_omega-1
	        do l=1,m0
	          ll=l1+l-1
	          if(ll>n_bnd(k,ispin)) cycle
	          sig_out(i_omega,ll,ks)=f(ll,ll,2,ind_omega,ind_k,ispin)
	        enddo
	      enddo
	    enddo
	  enddo
c --------- Collect from all processes  --------------------------------
	  if(nproc_k/=1) then
	    call dgop(sig_out,nbndf*n_sym*(n_omega+1),'  +',comm_pnt)
	  endif
	  if(nproc_tau/=1) then
	    call dgop(sig_out,nbndf*n_sym*(n_omega+1),'  +',comm_pnk)
	  endif
	  if(maswrk) write(3,*)' Ispin = ',ispin
	  do i=1,n
          l1=m*(i-1)+1
          l2=min(m*i,n0)
	    if(maswrk) write(3,*)' Bands from ',l1,' to ',l2
	    do ks=1,n_sym
	      k=ksym_in_npnt(ks)
	      if(maswrk) then
		      write(3,'(a5,3(1x,e12.5),1x,a1)')' K = ',
     &		          (pnt(i7,k),i7=1,3),sym_label(k)
	        do i_omega=0,n_omega
	          write(3,'(f11.4,5(1x,f12.7))')w_omega(i_omega)*fc0,
     &		               (sig_out(i_omega,ii,ks)*fc,ii=l1,l2)
	        enddo
	      endif
	    enddo
	  enddo
	enddo
	if(maswrk) close(3)
	if(maswrk) then
	  i_len=len_trim(allfile)
	  if(ubi==' gw') then
	    open(3,file=allfile(1:i_len)//'_Re_'//ch1//'_band_k.gw')
	    open(4,file=allfile(1:i_len)//'_Re_'//ch2//'_band_k.gw')
	  else if(ubi==' qp') then
	    open(3,file=allfile(1:i_len)//'_Re_'//ch1//'_band_k.qp')
	    open(4,file=allfile(1:i_len)//'_Re_'//ch2//'_band_k.qp')
	  else if(ubi=='gwg') then
	    open(3,file=allfile(1:i_len)//'_Re_'//ch1//'_band_k.gwg')
	    open(4,file=allfile(1:i_len)//'_Re_'//ch2//'_band_k.gwg')
	  endif
	endif
	do ispin=1,nspin
	  sig_out=0.d0
	  do i=1,n
          l1=m*(i-1)+1
          l2=min(m*i,n0)
	    m0=l2-l1+1
	    do ind_k=1,ndim3_k(me3_k+1)
	      k=n3_mpi_k(me3_k+1)+ind_k
	      ks=k_npnt_in_sym(k)
	      if(ks==0) cycle
	      do ind_omega=1,ndim3_omega
	        i_omega=me3_tau*ndim3_omega+ind_omega-1
	        do l=1,m0
	          ll=l1+l-1
	          if(ll>n_bnd(k,ispin)) cycle
	          sig_out(i_omega,ll,ks)=f(ll,ll,1,ind_omega,ind_k,ispin)
	        enddo
	      enddo
	    enddo
	  enddo
c --------- Collect from all processes  -------------------------------
	  if(nproc_k/=1) then
	    call dgop(sig_out,nbndf*n_sym*(n_omega+1),'  +',comm_pnt)
	  endif
	  if(nproc_tau/=1) then
	    call dgop(sig_out,nbndf*n_sym*(n_omega+1),'  +',comm_pnk)
	  endif
	  if(maswrk) write(3,*)' Ispin = ',ispin
	  if(maswrk) write(4,*)' Ispin = ',ispin
	  do i=1,n
          l1=m*(i-1)+1
          l2=min(m*i,n0)
	    if(maswrk) write(3,*)' Bands from ',l1,' to ',l2
	    do ks=1,n_sym
	      k=ksym_in_npnt(ks)
	      if(maswrk) then
		      write(3,'(a5,3(1x,e12.5),1x,a1)')' K = ',
     &		          (pnt(i7,k),i7=1,3),sym_label(k)
	        do i_omega=0,n_omega
	          write(3,'(f11.4,5(1x,f12.7))')w_omega(i_omega)*fc0,
     &		               (sig_out(i_omega,ii,ks)*fc,ii=l1,l2)
	        enddo
	      endif
	    enddo
	  enddo
	  if(maswrk) then
		  write(4,'(4x,8(4x,a1,5x))')(sym0_label(ks),
     &		                          ks=1,min(n_sym,nnn))
	    do i=1,nbndf
	      write(4,'(i4,8(1x,f9.4))')i,
     &		               (sig_out(0,i,ks)*fc,ks=1,min(n_sym,nnn))
	    enddo
	  endif
	enddo
	if(maswrk) close(3)
	if(maswrk) close(4)
	deallocate(sig_out)
      end
