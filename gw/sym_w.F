      subroutine sym_w(iq,iq0,w_0)
	use atom_mod
	use solid_mod
      implicit none
      integer, intent(in) :: iq,iq0
	complex*16, intent(inout) :: w_0(n_pbtot,n_pbtot)
	integer :: j,i,ig
	complex*16, allocatable :: tmp(:,:)
c -------- Copy+Symmetrization of W into another k-point ---------------
	ig=k_group(iq)
      if(iq/=iq0) then
	  call sym_w_mt(ig,pnt(1,iq),w_0,n_pbtot)
	  allocate(tmp(n_pbmt,nplw_gw))
	  tmp=w_0(1:n_pbmt,n_pbmt+1:n_pbtot)
	  call sym_w_mi(iq,ig,tmp)
	  w_0(1:n_pbmt,n_pbmt+1:n_pbtot)=tmp
	  deallocate(tmp)
	  allocate(tmp(nplw_gw,nplw_gw))
	  tmp=w_0(n_pbmt+1:n_pbtot,n_pbmt+1:n_pbtot)
	  call sym_w_ii(iq,ig,tmp,nplwgw(iq0),nplw_gw)
	  w_0(n_pbmt+1:n_pbtot,n_pbmt+1:n_pbtot)=tmp
	  deallocate(tmp)
	endif
c ------- Obtaining  W_im from W_mi -------------------------------
	do j=1,n_pbmt
	  do i=n_pbmt+1,n_pbtot
	    w_0(i,j)=conjg(w_0(j,i))
	  enddo
	enddo
      end
