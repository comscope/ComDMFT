      subroutine warp
      use manager_mod
	use solid_mod
      implicit none
	integer :: iplw,ia,ib,ic
	real*8 :: gx,gy,gz,pvr,pvi
c     &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      do iplw=1,numplw
        ia=igbs(1,iplw)
        ib=igbs(2,iplw)
        ic=igbs(3,iplw)
        gx=ia*gbas(1,1)+ib*gbas(1,2)+ic*gbas(1,3)
        gy=ia*gbas(2,1)+ib*gbas(2,2)+ic*gbas(2,3)
        gz=ia*gbas(3,1)+ib*gbas(3,2)+ic*gbas(3,3)
        call plwovi(gx,gy,gz,pvr,pvi)
        sovr(iplw)=pvr
        if(complex_ro) sovi(iplw)=pvi
      enddo
      call timel('*******  warp finished *************')
      end