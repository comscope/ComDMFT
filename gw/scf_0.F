      subroutine scf_0
      use parallel_mod
	use solid_mod
	use units_mod
      implicit none
      integer :: ind_k,ia,k
      real*8 :: fas,v(3),pi2
c     &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
	call freeatom
      call timel('**** Freeatom finished *************')
      allocate(phase_k_rr_all(nqdiv,maxval(ndim3_k),2))
      pi2=pi+pi
      do ia=1,nqdiv
        v(:)=rr_pnt(1,ia)*rbas(:,1)+rr_pnt(2,ia)*rbas(:,2)
     &                             +rr_pnt(3,ia)*rbas(:,3)
	  do ind_k=1,ndim3_k(me3_k+1)
	    k=n3_mpi_k(me3_k+1)+ind_k
	    fas=pi2*dot_product(pnt(:,k),v)
	    phase_k_rr_all(ia,ind_k,1)=cos(fas)
	    phase_k_rr_all(ia,ind_k,2)=sin(fas)
	  enddo
	enddo
      end
