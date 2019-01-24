      subroutine spgprd(g1,g1_j,a1,g2,g2_j,a2,ifi1,ifi2,g,g_j,a,ifi,
     &	              lmax,jmax,su1,su1_j)
	use manager_mod
	use solid_mod
	implicit none
	integer :: lmax,jmax,ifi1,ifi2,l,iwig,i,ist,ii,jj,ifi
      real*8 g1(maxwig),g2(maxwig),g(maxwig),a1(3),a2(3),a(3),
     &       su1(2*lmax+1,2*lmax+1)
	complex*16 g1_j(maxwig_j),g2_j(maxwig_j),g_j(maxwig_j),
     &	       su1_j(jmax+1,jmax+1)
      ifi=ifi1*ifi2
      do l=0,lmax
        iwig=l*(2*l-1)*(2*l+1)/3+1
        call dcopy((2*l+1)**2,g1(iwig),1,su1,1)
        call prd(su1,g2(iwig),g(iwig),l)
      enddo
      call rotate(a2(1),a2(2),a2(3),a(1),a(2),a(3),g1(2),2)
      do i=1,3
        a(i)=a1(i)+a(i)
      enddo
	if(irel.eq.2) then
        do l=0,lmax
	    ist=-1
	    if(l.eq.0) ist=1
	    do ii=ist,1,2
	      jj=l+l+ii
            iwig=4*l*(2*l*l+1)/3+2*(ii+1)*l*l+1
            call zcopy((jj+1)**2,g1_j(iwig),1,su1_j,1)
            call prd_j(su1_j,g2_j(iwig),g_j(iwig),jj)
          enddo
	  enddo
	endif
      end