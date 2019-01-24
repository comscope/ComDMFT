      logical function spgeql(g1,a1,g2,a2,ifi1,ifi2,qb,maxwig)
c  determines whether g1 is equal to g2
      implicit none
	integer :: ifi1,ifi2,maxwig,m,iq,iac
      real*8 :: g1(maxwig),g2(maxwig),a1(3),a2(3),qb(3,3),c,ca,dc
      spgeql=.true.
      if(ifi1.ne.ifi2) then
        spgeql=.false.
        return
      endif
      do m=1,maxwig
        if(dabs(g1(m)-g2(m)).gt.1.d-5) then
          spgeql=.false.
          return
        endif
      enddo
      do iq=1,3
        c=(a1(1)-a2(1))*qb(1,iq)+(a1(2)-a2(2))*qb(2,iq)
     &                          +(a1(3)-a2(3))*qb(3,iq)
        ca=dabs(c)
        iac=ca+0.5d0
        dc=ca-iac
        if(dabs(dc).gt.1.d-5) then
          spgeql=.false.
          return
        endif
      enddo
      end