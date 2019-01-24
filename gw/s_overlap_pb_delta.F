      subroutine s_overlap_pb_delta(s2,iq)
c     Calculates S2=S1*S^(-1) where S1 is not square
c     and S - interstitial overlap matrix for product basis ------
	use atom_mod
	use etot_mod
	use manager_mod
	use solid_mod
	use units_mod
      implicit none
      integer, intent(in) :: iq
      complex*16, intent(out) :: s2(nplwdif,nplw_gw)
	integer :: j,j0,i,ia,ib,ic,ind,n,iq0
	complex*16, allocatable :: s(:,:),s1(:,:)
	iq0=i_kref(iq)
	n=nplwgw(iq0)
	allocate(s(nplw_gw,nplw_gw))
	allocate(s1(nplwdif,n))
c ------ We form the inverse overlap matrix -----------------------------
      call s_overlap(1,s,nplw_gw,n,indpw_gw(1,iq))
c ------ We form extended (non-square)  overlap matrix ------------------
	do j=1,n
	  j0=indpw_gw(j,iq)
	  do i=1,nplwdif
          ia=igbs(1,j0)-igbs(1,i)
          ib=igbs(2,j0)-igbs(2,i)
          ic=igbs(3,j0)-igbs(3,i)
          ind=indplw(ia,ib,ic)
	    if(complex_ro) then
		    s1(i,j)=dcmplx(sovr(ind),sovi(ind))
	    else
		    s1(i,j)=dcmplx(sovr(ind),0.d0)
	    endif
	  enddo
	enddo
	s1=s1*amega
      call zgemm('n','n',nplwdif,n,n,(1.d0,0.d0),s1,nplwdif,s,nplw_gw,
     &		     (0.d0,0.d0),s2,nplwdif)
	deallocate(s,s1)
      end