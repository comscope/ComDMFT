      subroutine output_sigma_vrt(f,ch1,fc)
	use atom_mod
	use manager_mod
	use models_mod
	use parallel_mod
	use solid_mod
	use units_mod
      implicit none
      character*4, intent(in) :: ch1
      real*8, intent(in) :: fc
      complex*16, intent(in) :: f(nbndf_bnd,nbndf_bnd,2,ndim3_tau,
     &                            maxval(ndim3_k),nspin)
	integer :: k,i_len,m,n,ispin,l1,l2,m0,i_omega,l,ll,i,ii,n0,ks,i7,
     &	         nnn,ind_k
      real*8 :: fc0
	complex*16, allocatable :: sig_out(:,:,:),s_om(:,:,:)
	nnn=8
	fc0=evolt/2.d0
	n0=nbndf_bnd
	m=5
	n=n0/m
	if(n*m<n0) n=n+1
	allocate(sig_out(0:n_omega,nbndf_bnd,n_sym))
	allocate(s_om(nbndf_bnd,nbndf_bnd,0:n_omega))
	if(maswrk) then
	  i_len=len_trim(allfile)
	  open(3,file=allfile(1:i_len)//'_Im_'//ch1//'_band_k_VRT.gwg')
	  open(4,file=allfile(1:i_len)//'_Re_'//ch1//'_band_k_VRT.gwg')
	endif
	do ispin=1,nspin
	  sig_out=(0.d0,0.d0)
	  do ind_k=1,ndim3_k(me3_k+1)
	    k=n3_mpi_k(me3_k+1)+ind_k
	    ks=k_npnt_in_sym(k)
	    if(ks==0) cycle
	    call from_tau_to_omega_vrt(s_om,f(1,1,1,1,ind_k,ispin),n0*n0)
	    do i=1,n
            l1=m*(i-1)+1
            l2=min(m*i,n0)
	      m0=l2-l1+1
	      do i_omega=0,n_omega
	        do l=1,m0
	          ll=l1+l-1
	          if(ll>n_low_bnd(k,ispin)) cycle
	          sig_out(i_omega,ll,ks)=s_om(ll,ll,i_omega)
	        enddo
	      enddo
	    enddo
	  enddo
	  if(nproc_k/=1) call dgop(sig_out,2*(n_omega+1)*nbndf_bnd*n_sym,
     &	                         '  +',comm_pnt)
	  if(maswrk) write(3,*)' Ispin = ',ispin
	  if(maswrk) write(4,*)' Ispin = ',ispin
	  do i=1,n
          l1=m*(i-1)+1
          l2=min(m*i,n0)
	    if(maswrk) write(3,*)' Bands from ',l1,' to ',l2
	    if(maswrk) write(4,*)' Bands from ',l1,' to ',l2
	    do ks=1,n_sym
	      k=ksym_in_npnt(ks)
	      if(maswrk) then
		      write(3,'(a5,3(1x,e12.5),1x,a1)')' K = ',
     &		          (pnt(i7,k),i7=1,3),sym_label(k)
		      write(4,'(a5,3(1x,e12.5),1x,a1)')' K = ',
     &		          (pnt(i7,k),i7=1,3),sym_label(k)
	        do i_omega=0,n_omega
	          write(3,'(f11.4,5(1x,f12.7))')w_omega(i_omega)*fc0,
     &		               (imag(sig_out(i_omega,ii,ks))*fc,ii=l1,l2)
	          write(4,'(f11.4,5(1x,f12.7))')w_omega(i_omega)*fc0,
     &		               (real(sig_out(i_omega,ii,ks))*fc,ii=l1,l2)
	        enddo
	      endif
	    enddo
	  enddo
	enddo
	if(maswrk) close(3)
	if(maswrk) close(4)
	deallocate(sig_out,s_om)
      end
