      subroutine setpw_gw
	use manager_mod
	use solid_mod
      implicit none
      integer :: k,num,idum
      nplw_gw=0
	do k=1,npnt
	  call g_inside_cutoff(numplw,pnt(1,k),cutoffgw,.false.,.false.,
     &	                     idum,nplwgw(k),nplwdif)
        nplw_gw=max(nplw_gw,nplwgw(k))
      enddo
	do k=1,num_k_plot
	  call g_inside_cutoff(numplw,k_plot(1,k),cutoffgw,.false.,
     &	                     .false.,idum,nplwgw_plot(k),nplwdif)
        nplw_gw=max(nplw_gw,nplwgw_plot(k))
      enddo
	allocate(indpw_gw(nplw_gw,nqdiv))
	allocate(indpw_gw_plot(nplw_gw,num_k_plot))
	do k=1,npnt
	  call g_inside_cutoff(numplw,pnt(1,k),cutoffgw,.false.,.true.,
     &	                     indpw_gw(1,k),num,nplw_gw)
      enddo
	do k=1,num_k_plot
	  call g_inside_cutoff(numplw,k_plot(1,k),cutoffgw,.false.,.true.,
     &	                     indpw_gw_plot(1,k),num,nplw_gw)
      enddo
      end