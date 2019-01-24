	subroutine get_indk1(ind_k1)
	use solid_mod
	integer :: ind_k1(nqdiv,npnt),iq,k,i,ind_kq
	real*8 :: qq0(3)
	do iq=1,npnt
	  do k=1,nqdiv
	    do i=1,3
	      qq0(i)=pnt(i,k)-pnt(i,iq)
	    enddo
	    call zone1_number(qq0,rb0,ndiv,ind_kq)
	    ind_k1(k,iq)=index_k1(ind_kq)
        enddo
	enddo
	end