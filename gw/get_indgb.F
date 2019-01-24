      subroutine get_indgb
	use manager_mod
	use solid_mod
	use units_mod
      implicit none
	integer :: nbndt,k
	real*8 :: cut
	do k=1,npnt
	  call g_inside_cutoff(nplwbas,pnt(1,k),cutoff,.true.,.true.,
     &	                     indgb(1,k),nbndt,nbndtm)
	enddo
	do k=1,num_k_plot
	  call g_inside_cutoff(nplwbas,k_plot(1,k),cutoff,.true.,.true.,
     &	                     indgb_plot(1,k),nbndt,nbndtm)
	enddo
      end