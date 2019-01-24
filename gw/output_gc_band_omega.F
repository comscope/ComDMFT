      subroutine output_gc_band_omega(ispin,k,n0)
	use atom_mod
	use etot_mod
	use manager_mod
	use parallel_mod
	use solid_mod
	use units_mod
      implicit none
	integer, intent(in) :: ispin,k,n0
	integer :: i_len,m,n,l1,l2,m0,i_omega,l,ll,i,ii,ind_omega
	real*8 :: fc
	real*8, allocatable :: g_out(:,:)
	m=5
	n=n0/m
	if(n*m<n0) n=n+1
	fc=evolt/2.d0
	allocate(g_out(0:n_omega,m))
	if(maswrk) then
	  i_len=len_trim(allfile)
	  if(ubi==' gw') then
	    open(3,file=allfile(1:i_len)//'_Im_gc_w_band_k.gw')
	  else if(ubi==' qp') then
	    open(3,file=allfile(1:i_len)//'_Im_gc_w_band_k.qp')
	  else if(ubi=='gwg') then
	    open(3,file=allfile(1:i_len)//'_Im_gc_w_band_k.gwg')
	  else if(ubi=='pt1') then
	    open(3,file=allfile(1:i_len)//'_Im_gc_w_band_k.pt1')
	  else if(ubi=='pt2') then
	    open(3,file=allfile(1:i_len)//'_Im_gc_w_band_k.pt2')
	  endif
	  if(ispin/=1.or.k/=1) call goto_end(3)
	  write(3,*)' Ispin = ',ispin,' K = ',k
	endif
	do i=1,n
        l1=m*(i-1)+1
        l2=min(m*i,n_bnd(k,ispin))
	  if(maswrk) write(3,*)' Bands from ',l1,' to ',l2
	  m0=l2-l1+1
	  if(m0<=0) cycle
	  g_out=0.d0
	  do ind_omega=1,ndim3_omega
	    i_omega=me3_tau*ndim3_omega+ind_omega-1
	    do l=1,m0
	      ll=l1+l-1
	      g_out(i_omega,l)=gc_omega(ll,ll,2,ind_omega)
	    enddo
	  enddo
	  if(nproc_tau/=1) call dgop(g_out,(n_omega+1)*m0,'  +',comm_pnk)
	  if(maswrk) then
	    do i_omega=0,n_omega
	      write(3,'(f9.3,5(1x,f12.7))')w_omega(i_omega)*fc,
     &		               (g_out(i_omega,ii),ii=1,m0)
	    enddo
	  endif
	enddo
	if(maswrk) close(3)
	if(maswrk) then
	  i_len=len_trim(allfile)
	  if(ubi==' gw') then
	    open(3,file=allfile(1:i_len)//'_Re_gc_w_band_k.gw')
	  else if(ubi==' qp') then
	    open(3,file=allfile(1:i_len)//'_Re_gc_w_band_k.qp')
	  else if(ubi=='gwg') then
	    open(3,file=allfile(1:i_len)//'_Re_gc_w_band_k.gwg')
	  else if(ubi=='pt1') then
	    open(3,file=allfile(1:i_len)//'_Re_gc_w_band_k.pt1')
	  else if(ubi=='pt2') then
	    open(3,file=allfile(1:i_len)//'_Re_gc_w_band_k.pt2')
	  endif
	  if(ispin/=1.or.k/=1) call goto_end(3)
	  write(3,*)' Ispin = ',ispin,' K = ',k
	endif
	do i=1,n
        l1=m*(i-1)+1
        l2=min(m*i,n_bnd(k,ispin))
	  if(maswrk) write(3,*)' Bands from ',l1,' to ',l2
	  m0=l2-l1+1
	  if(m0<=0) cycle
	  g_out=0.d0
	  do ind_omega=1,ndim3_omega
	    i_omega=me3_tau*ndim3_omega+ind_omega-1
	    do l=1,m0
	      ll=l1+l-1
	      g_out(i_omega,l)=gc_omega(ll,ll,1,ind_omega)
	    enddo
	  enddo
	  if(nproc_tau/=1) call dgop(g_out,(n_omega+1)*m0,'  +',comm_pnk)
	  if(maswrk) then
	    do i_omega=0,n_omega
	      write(3,'(f9.3,5(1x,f12.7))')w_omega(i_omega)*fc,
     &		               (g_out(i_omega,ii),ii=1,m0)
	    enddo
	  endif
	enddo
	if(maswrk) close(3)
	deallocate(g_out)
      end