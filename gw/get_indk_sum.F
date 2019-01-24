	subroutine get_indk_sum(ind_k)
	use solid_mod
	implicit none
	integer :: ind_k(nqdiv,npnt),iq,k,i,ind_kq
	real*8 :: qq0(3)
	do iq=1,nqdiv
	  do k=1,npnt
	    do i=1,3
	      qq0(i)=pnt(i,k)+pnt(i,iq)
	    enddo
	    call zone1_number(qq0,rb0,ndiv,ind_kq)
	    ind_k(iq,k)=index_k1(ind_kq)
        enddo
	enddo
	end