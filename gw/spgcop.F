      subroutine spgcop(g,g_j,ag0,h,h_j,ah,ifi,ifih)
	use manager_mod
	use solid_mod
	integer :: ifi,ifih,i
      real*8 h(maxwig),g(maxwig),ag0(3),ah(3)
	complex*16 h_j(maxwig_j),g_j(maxwig_j)
      ifih=ifi
      do i=1,maxwig
         h(i)=g(i)
         if(dabs(h(i)).lt.1.d-10) h(i)=0.d0
      enddo
      do i=1,3
         ah(i)=ag0(i)
         if(dabs(ah(i)).lt.1.d-10) ah(i)=0.d0
      enddo
	if(irel.eq.2) then
        do i=1,maxwig_j
          h_j(i)=g_j(i)
          if(abs(h_j(i)).lt.1.d-10) h_j(i)=(0.d0,0.d0)
        enddo
	endif
      end
