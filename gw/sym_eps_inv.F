      subroutine sym_eps_inv(iq,iq0,w_0)
	use atom_mod
	use solid_mod
      implicit none
      integer, intent(in) :: iq,iq0
	complex*16, intent(inout) :: w_0(n_pbtot,n_pbtot)
	integer :: ig
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
	  allocate(tmp(nplw_gw,n_pbmt))
	  tmp=w_0(n_pbmt+1:n_pbtot,1:n_pbmt)
	  call sym_w_im(iq,ig,tmp)
	  w_0(n_pbmt+1:n_pbtot,1:n_pbmt)=tmp
	  deallocate(tmp)
	  allocate(tmp(nplw_gw,nplw_gw))
	  tmp=w_0(n_pbmt+1:n_pbtot,n_pbmt+1:n_pbtot)
	  call sym_w_ii(iq,ig,tmp,nplwgw(iq0),nplw_gw)
	  w_0(n_pbmt+1:n_pbtot,n_pbmt+1:n_pbtot)=tmp
	  deallocate(tmp)
	endif
      end