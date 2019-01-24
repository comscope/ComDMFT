      subroutine get_k(ierror)
	use manager_mod
	use parallel_mod
	use solid_mod
      implicit none
	integer :: ierror
c     &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      call sbzmesh(pnt,wgt,idtet,iwtet,npnt,ntet,ierror)
      if(maswrk) write(iun,1030)npnt
      if(maswrk) write(iun,1033)ntet
1030  format(1x,i4,' k-points generated')
1033  format(1x,i4,' K-Tetrahedrons generated')
      end