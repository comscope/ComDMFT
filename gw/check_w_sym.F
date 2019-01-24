      subroutine check_w_sym(iq,w,n0,key)
c     key = 2 The matrix is unpacked
c     key = 1 The matrix is packed          
	use atom_mod
	use manager_mod
	use parallel_mod
	use units_mod
	use solid_mod
      implicit none
	integer, intent(in) :: key,n0,iq
	real*8, intent(in)  :: w(key*n0,n0)
	integer :: ig,j,i,n,ii,j0,j1,i0,i1,gbs_number,igtild
	real*8 :: d_mm,d_mi,d_ii,pi2,v(3),phase
	complex*16 :: cc
	complex*16, allocatable :: tmp(:,:),tmp0(:,:),tmi(:,:)
	if(itest==0) return
	n=n_pbmt+nplwgw(iq)
	pi2=pi+pi
	allocate(tmp(n,n))
	allocate(tmp0(n,n))
	allocate(tmi(n_pbmt,nplw_gw))
	d_mm=0.d0
	d_mi=0.d0
	d_ii=0.d0
	do ii=1,num_sym_k(iq)
	  ig=list_sym_k(ii,iq)
	  igtild=g_sym_0(iq,ig)
	  if(key==2) then
	    call dlacpy('a',2*n,n,w,2*n0,tmp0,2*n)
	  else
		  call boson_unpack_tau(tmp0,n,w,n0,n)
		endif
		tmp=tmp0
c ----------------- MM part ---------------------------------		
	  call sym_w_mt(ig,pnt(1,iq),tmp,n)
        do j=1,n_pbmt
          do i=1,n_pbmt
            d_mm=max(d_mm,abs(tmp(i,j)-tmp0(i,j)))
          enddo
        enddo
c ----------------- MI part ---------------------------------
        tmi(:,1:nplwgw(iq))=tmp0(1:n_pbmt,n_pbmt+1:n)
        call sym_w_mi(iq,ig,tmi)
        do j=n_pbmt+1,n
          do i=1,n_pbmt
            d_mi=max(d_mi,abs(tmi(i,j-n_pbmt)-tmp0(i,j)))
          enddo
        enddo
c ----------------- II part ---------------------------------
        do j=1,nplwgw(iq)
	    j0=indpw_gw(j,iq)  !! G' in GBS-list
	    v=gbs(:,j0)+gbs(:,igtild)  !! G'+G_A
          j1=gbs_number(v)   !! G'+G_A in GBS-list
          j1=gbs_sym(j1,ig)   !! A^(-1)*(G'+G_A) in GBS-list
          j1=iplf_gk(j1,iq)
          do i=1,nplwgw(iq)
	      i0=indpw_gw(i,iq)
	      v=gbs(:,i0)+gbs(:,igtild)  !! G+G_A
            i1=gbs_number(v)   !! G+G_A in GBS-list
            i1=gbs_sym(i1,ig)   !! A^(-1)*(G+G_A) in GBS-list
            i1=iplf_gk(i1,iq)
	      v=gbs(:,j0)-gbs(:,i0)
	      phase=pi2*dot_product(v,shift(:,ig))
	      cc=dcmplx(cos(phase),sin(phase))
	      cc=cc*tmp0(n_pbmt+i1,n_pbmt+j1)
            d_ii=max(d_ii,abs(cc-tmp0(n_pbmt+i,n_pbmt+j)))
	    enddo
	  enddo   !! over j
	enddo
	if(maswrk) then
	  write(iun,'(a5,i3,3(1x,a3,1x,e12.5))')' Q = ',iq,
     &	  'MM:',d_mm,'MI:',d_mi,'II:',d_ii 
      endif
	deallocate(tmp,tmp0,tmi)
      end