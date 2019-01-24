      subroutine reduc(a,b,c,rin,rout,key)
*************************************************************
*       subroutine for reducing vector r in first (key=1) , *
*       central (key=0) lattice cell.if (key=-1) return the *
*       r in basis coordinate.                              *
*************************************************************
	use manager_mod
	use parallel_mod
      implicit none
	integer :: i,j,key,i0,j0,k0,l,k
      real*8 :: rin(3),rout(3),rw(3),a(3),b(3),c(3),abc(3,3),abci(3,3),
     &	      det,s,rc,rmod
c start up:
      do i=1,3
        abc(i,1)=a(i)
        abc(i,2)=b(i)
        abc(i,3)=c(i)
	enddo
      abci(1,1)=b(2)*c(3)-b(3)*c(2)
      abci(1,2)=b(3)*c(1)-b(1)*c(3)
      abci(1,3)=b(1)*c(2)-b(2)*c(1)
      det=a(1)*abci(1,1)+a(2)*abci(1,2)+a(3)*abci(1,3)
      if (abs(det).le.1.d-08) then
        if(maswrk) write(iun,*)'invalid basis in reduc'
        call ending
      endif
      abci(2,1)=a(3)*c(2)-a(2)*c(3)
      abci(2,2)=a(1)*c(3)-a(3)*c(1)
      abci(2,3)=a(2)*c(1)-a(1)*c(2)
      abci(3,1)=a(2)*b(3)-a(3)*b(2)
      abci(3,2)=a(3)*b(1)-a(1)*b(3)
      abci(3,3)=a(1)*b(2)-a(2)*b(1)
      do i=1,3
        do j=1,3
          abci(i,j)=abci(i,j)/det
	  enddo
	enddo
c
c rin in a,b,c basis:
      do i=1,3
        s=0
        do j=1,3
          s=s+abci(i,j)*rin(j)
	  enddo
        rout(i)=s
	enddo
      if(key.lt.0) return
c
      rw(1)=rout(1)-aint(rout(1))
      rw(2)=rout(2)-aint(rout(2))
      rw(3)=rout(3)-aint(rout(3))
c first cell:
      if(rw(1).lt.0.) rw(1)=rw(1)+1
      if(rw(2).lt.0.) rw(2)=rw(2)+1
      if(rw(3).lt.0.) rw(3)=rw(3)+1
c return to old basis:
      do i=1,3
        s=0
        do j=1,3
          s=s+abc(i,j)*rw(j)
	  enddo
        rout(i)=s
	enddo
      if(key.eq.1) return
c central cell:
      i0=0
      j0=0
      k0=0
      rc=rout(1)**2+rout(2)**2+rout(3)**2
      do i=0,1
        do j=0,1
          do k=0,1
            rmod=0
            do l=1,3
              rmod=rmod+(rout(l)-i*a(l)-j*b(l)-k*c(l))**2
	      enddo
            if(rmod.lt.rc) then
              i0=i
              j0=j
              k0=k
              rc=rmod
            endif
	    enddo
        enddo
      enddo
      do l=1,3
         rout(l)=rout(l)-i0*a(l)-j0*b(l)-k0*c(l)
	enddo
      end