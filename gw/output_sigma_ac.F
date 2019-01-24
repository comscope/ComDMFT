      subroutine output_sigma_ac(f,ch1,key)
c     key = 0 we do not transform Sigma from X to QP representation
c     key = 1 we do the transform of Sigma from X to QP representation       
	use atom_mod
	use manager_mod
	use parallel_mod
	use solid_mod
	use units_mod
      implicit none
      integer, intent(in) :: key
      character*6, intent(in) :: ch1
      real*8, intent(in) :: f(nbndf,nbndf,2,ndim3_omega,
     &                      maxval(ndim3_k),nspin)
	integer :: k,i_len,m,n,ispin,l1,l2,ind_k,ind_omega,i_omega,i,ii,
     &	         n0,ks,i7,nnn
      real*8 :: fc0,deltae
      real*8, allocatable :: erax(:)
	complex*16, allocatable :: sig(:,:,:,:),sig_out(:,:,:),tmp(:,:),
     &	                         a(:,:)
	allocate(sig(0:n_omega,nbndf,n_sym,nspin))
	allocate(sig_out(0:ndos,nbndf,n_sym))
	allocate(erax(0:ndos))
      deltae=(emaxdos-emindos)/ndos
	do i=0,ndos
	  erax(i)=emindos+deltae*i
	enddo
      sig=(0.d0,0.d0)
	do ispin=1,nspin
	  do ind_k=1,ndim3_k(me3_k+1)
	    k=n3_mpi_k(me3_k+1)+ind_k
	    ks=k_npnt_in_sym(k)
	    if(ks==0) cycle
	    n=n_bnd(k,ispin)
	    allocate(tmp(n,n))
	    allocate(a(n,n))
	    do ind_omega=1,ndim3_omega
	      i_omega=me3_tau*ndim3_omega+ind_omega-1
		    call ferm_unpack_omega(tmp,f(1,1,1,ind_omega,ind_k,ispin),
     &                             n,nbndf,n)
            if(key==1) then
 	        call zgemm('c','n',n,n,n,(1.d0,0.d0),
     &	                 q_qp(1,1,ind_k,ispin),nbndf,tmp,n,(0.d0,0.d0),
     &                   a,n)
	        call zgemm('n','n',n,n,n,(1.d0,0.d0),a,n,
     &	                 q_qp(1,1,ind_k,ispin),nbndf,(0.d0,0.d0),tmp,n)
            endif
	      do i=1,n_bnd(k,ispin)
	        sig(i_omega,i,ks,ispin)=tmp(i,i)
	      enddo
	    enddo
	    deallocate(tmp,a)
	  enddo
	enddo
c --------- Collect from all processes  --------------------------------
	if(nproc_k/=1) then
	  call dgop(sig,2*nbndf*n_sym*(n_omega+1)*nspin,'  +',comm_pnt)
	endif
	if(nproc_tau/=1) then
	  call dgop(sig,2*nbndf*n_sym*(n_omega+1)*nspin,'  +',comm_pnk)
	endif
	if(maswrk) then
	  i_len=len_trim(allfile)
	  if(ubi==' gw') then
	    open(3,file=allfile(1:i_len)//'_Re_'//ch1//'_band_ac.gw')
	    open(4,file=allfile(1:i_len)//'_Im_'//ch1//'_band_ac.gw')
	  else if(ubi==' qp') then
	    open(3,file=allfile(1:i_len)//'_Re_'//ch1//'_band_ac.qp')
	    open(4,file=allfile(1:i_len)//'_Im_'//ch1//'_band_ac.qp')
	  else if(ubi=='gwg') then
	    open(3,file=allfile(1:i_len)//'_Re_'//ch1//'_band_ac.gwg')
	    open(4,file=allfile(1:i_len)//'_Im_'//ch1//'_band_ac.gwg')
	  else if(ubi=='pt1') then
	    open(3,file=allfile(1:i_len)//'_Re_'//ch1//'_band_ac.pt1')
	    open(4,file=allfile(1:i_len)//'_Im_'//ch1//'_band_ac.pt1')
	  else if(ubi=='pt2') then
	    open(3,file=allfile(1:i_len)//'_Re_'//ch1//'_band_ac.pt2')
	    open(4,file=allfile(1:i_len)//'_Im_'//ch1//'_band_ac.pt2')
	  endif
	endif
	nnn=8
	fc0=evolt/2.d0
	n0=nbndf
	m=5
	n=n0/m
	if(n*m<n0) n=n+1
	do ispin=1,nspin
	  if(maswrk) write(3,*)' Ispin = ',ispin
	  if(maswrk) write(4,*)' Ispin = ',ispin
	  do ks=1,n_sym
          k=ksym_in_npnt(ks)
	    do i=1,n_bnd(k,ispin)
	      call a_cont_ferm(ndos+1,erax,sig(0,i,ks,ispin),
     &	                     sig_out(0,i,ks),0.d0)
	    enddo
	  enddo
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
	        do i_omega=0,ndos
	          write(3,'(f11.4,5(1x,f12.7))')erax(i_omega)*fc0,
     &		               (real(sig_out(i_omega,ii,ks))*fc0,ii=l1,l2)
	          write(4,'(f11.4,5(1x,f12.7))')erax(i_omega)*fc0,
     &		               (imag(sig_out(i_omega,ii,ks))*fc0,ii=l1,l2)
	        enddo
	      endif
	    enddo
	  enddo
	enddo
	if(maswrk) close(3)
	if(maswrk) close(4)
	deallocate(sig_out,sig,erax)
      end
