      subroutine angle(ll,mm,l,i,mj,l1,i1,mj1,key,angr,angi)
	use manager_mod
      implicit none
	integer, intent(in) :: ll,mm,l,i,mj,l1,i1,mj1,key
	real*8, intent(out) :: angr,angi
	integer :: lll,lll1,ii,ii1,lmj,lmj1,mu,mjm,mjm1,lm,l1m1,mmu,lm0,
     &	         l1m10,mjm2,mjm3
	real*8 :: sqr2,sckg,g,sckg1,sckg2,g1,g2
      data sqr2/1.41421356237310d0/
      angr=0.d0
      angi=0.d0
      if(ll.gt.l+l1.or.ll.lt.iabs(l-l1).or.mod(ll+l+l1,2).eq.1) return
      lll=l*(l+1)+1
      lll1=l1*(l1+1)+1
      ii=(3+i)/2
      ii1=(3+i1)/2
      lmj=l*(l+2)+(mj+3)/2         !!!  IV-37-7
      lmj1=l1*(l1+2)+(mj1+3)/2
      if(key.eq.0) then             !!! IV-37-6
        do mu=-1,1,2
          mjm=(mj-mu)/2
          if(iabs(mjm).gt.l) goto 1
          mjm1=(mj1-mu)/2
          if(iabs(mjm1).gt.l1) goto 1
          if(mm+mjm1.ne.mjm.and.mjm1-mm.ne.mjm) goto 1
          lm=lll+mjm
          l1m1=lll1+mjm1
          mmu=(3+mu)/2
          sckg=ckg2(ii,mmu,lmj)*ckg2(ii1,mmu,lmj1)
          g=gaunt(ll/2,lm,l1m1)
          if(mm.lt.0) then
            if(mj-mj1.ge.0) g=-g*(-1)**mm
            angi=angi+sckg*g/sqr2
          else if(mm.eq.0) then
            angr=angr+sckg*g
          else if(mm.gt.0) then
            if(mj-mj1.gt.0) g=g*(-1)**mm
            angr=angr+sckg*g/sqr2
          endif
1         continue
        enddo !!! over mu
      else if(key.eq.1) then          !!! $ IV-21
        mjm=(mj+1)/2
        mjm1=(mj1-1)/2
        sckg1=0.d0
        if(mjm1+mm.eq.mjm) sckg1=ckg2(ii,1,lmj)*ckg2(ii1,2,lmj1)
        if(mjm1-mm.eq.mjm) sckg1=ckg2(ii,1,lmj)*ckg2(ii1,2,lmj1)
        lm0=lll+mjm
        l1m10=lll1+mjm1
        mjm2=(mj-1)/2
        mjm3=(mj1+1)/2
        sckg2=0.d0
        if(mjm3+mm.eq.mjm2) sckg2=ckg2(ii,2,lmj)*ckg2(ii1,1,lmj1)
        if(mjm3-mm.eq.mjm2) sckg2=ckg2(ii,2,lmj)*ckg2(ii1,1,lmj1)
        lm=lll+mjm2
        l1m1=lll1+mjm3
        if(mm.lt.0) then
          if(iabs(mjm).gt.l.or.iabs(mjm1).gt.l1) then
            g1=0.d0
          else
            g1=gaunt(ll/2,lm0,l1m10)
            if((mj-mj1)/2+1.gt.0) g1=-g1*(-1)**mm
          endif
          if(iabs(mjm2).gt.l.or.iabs(mjm3).gt.l1) then
            g2=0.d0
          else
            g2=gaunt(ll/2,lm,l1m1)
            if((mj-mj1)/2-1.gt.0) g2=-g2*(-1)**(-mm)
          endif
          angi=angi+(sckg1*g1+sckg2*g2)/sqr2
        else if(mm.eq.0) then
          if(iabs(mjm).gt.l.or.iabs(mjm1).gt.l1) then
            g1=0.d0
          else
            g1=gaunt(ll/2,lm0,l1m10)
          endif
          if(iabs(mjm2).gt.l.or.iabs(mjm3).gt.l1) then
            g2=0.d0
          else
            g2=gaunt(ll/2,lm,l1m1)
          endif
          angr=angr+sckg1*g1+sckg2*g2
        else if(mm.gt.0) then
          if(iabs(mjm).gt.l.or.iabs(mjm1).gt.l1) then
            g1=0.d0
          else
            g1=gaunt(ll/2,lm0,l1m10)
            if((mj-mj1)/2+1.gt.0) g1=g1*(-1)**mm
          endif
          if(iabs(mjm2).gt.l.or.iabs(mjm3).gt.l1) then
            g2=0.d0
          else
            g2=gaunt(ll/2,lm,l1m1)
            if((mj-mj1)/2-1.gt.0) g2=g2*(-1)**mm
          endif
          angr=angr+(sckg1*g1+sckg2*g2)/sqr2
        endif
      else if(key.eq.2) then
        mjm=(mj+1)/2
        lm0=lll+mjm
        mjm1=(mj1-1)/2
        l1m10=lll1+mjm1
        sckg1=0.d0
        if(mjm1+mm.eq.mjm) sckg1=ckg2(ii,1,lmj)*ckg2(ii1,2,lmj1)
        if(mjm1-mm.eq.mjm) sckg1=ckg2(ii,1,lmj)*ckg2(ii1,2,lmj1)
        mjm2=(mj-1)/2
        lm=lll+mjm2
        mjm3=(mj1+1)/2
        l1m1=lll1+mjm3
        sckg2=0.d0
        if(mjm3+mm.eq.mjm2) sckg2=ckg2(ii,2,lmj)*ckg2(ii1,1,lmj1)
        if(mjm3-mm.eq.mjm2) sckg2=ckg2(ii,2,lmj)*ckg2(ii1,1,lmj1)
        if(mm.lt.0) then
          if(iabs(mjm).gt.l.or.iabs(mjm1).gt.l1) then
            g1=0.d0
          else
            g1=gaunt(ll/2,lm0,l1m10)
            if((mj-mj1)/2+1.gt.0) g1=-g1*(-1)**(-mm)
          endif
          if(iabs(mjm2).gt.l.or.iabs(mjm3).gt.l1) then
            g2=0.d0
          else
            g2=gaunt(ll/2,lm,l1m1)
            if((mj-mj1)/2-1.gt.0) g2=-g2*(-1)**(-mm)
          endif
          angr=angr-(sckg1*g1-sckg2*g2)/sqr2
        else if(mm.eq.0) then
          if(iabs(mjm).gt.l.or.iabs(mjm1).gt.l1) then
            g1=0.d0
          else
            g1=gaunt(ll/2,lm0,l1m10)
          endif
          if(iabs(mjm2).gt.l.or.iabs(mjm3).gt.l1) then
            g2=0.d0
          else
            g2=gaunt(ll/2,lm,l1m1)
          endif
          angi=angi+sckg1*g1-sckg2*g2
        else if(mm.gt.0) then
          if(iabs(mjm).gt.l.or.iabs(mjm1).gt.l1) then
            g1=0.d0
          else
            g1=gaunt(ll/2,lm0,l1m10)
            if((mj-mj1)/2+1.gt.0) g1=g1*(-1)**(mm)
          endif
          if(iabs(mjm2).gt.l.or.iabs(mjm3).gt.l1) then
            g2=0.d0
          else
            g2=gaunt(ll/2,lm,l1m1)
            if((mj-mj1)/2-1.gt.0) g2=g2*(-1)**(mm)
          endif
          angi=angi+(sckg1*g1-sckg2*g2)/sqr2
        endif
      else if(key.eq.3) then
        mjm=(mj-1)/2
        lm0=lll+mjm
        mjm1=(mj1-1)/2
        l1m10=lll1+mjm1
        sckg1=0.d0
        if(mjm1+mm.eq.mjm) sckg1=ckg2(ii,2,lmj)*ckg2(ii1,2,lmj1)
        if(mjm1-mm.eq.mjm) sckg1=ckg2(ii,2,lmj)*ckg2(ii1,2,lmj1)
        mjm2=(mj+1)/2
        lm=lll+mjm2
        mjm3=(mj1+1)/2
        l1m1=lll1+mjm3
        sckg2=0.d0
        if(mjm3+mm.eq.mjm2) sckg2=ckg2(ii,1,lmj)*ckg2(ii1,1,lmj1)
        if(mjm3-mm.eq.mjm2) sckg2=ckg2(ii,1,lmj)*ckg2(ii1,1,lmj1)
        if(mm.lt.0) then
          if(iabs(mjm).gt.l.or.iabs(mjm1).gt.l1) then
            g1=0.d0
          else
            g1=gaunt(ll/2,lm0,l1m10)
            if((mj-mj1)/2.gt.0) g1=-g1*(-1)**(-mm)
          endif
          if(iabs(mjm2).gt.l.or.iabs(mjm3).gt.l1) then
            g2=0.d0
          else
            g2=gaunt(ll/2,lm,l1m1)
            if((mj-mj1)/2.gt.0) g2=-g2*(-1)**(-mm)
          endif
          angi=angi+(sckg1*g1-sckg2*g2)/sqr2
        else if(mm.eq.0) then
          if(iabs(mjm).gt.l.or.iabs(mjm1).gt.l1) then
            g1=0.d0
          else
            g1=gaunt(ll/2,lm0,l1m10)
          endif
          if(iabs(mjm2).gt.l.or.iabs(mjm3).gt.l1) then
            g2=0.d0
          else
            g2=gaunt(ll/2,lm,l1m1)
          endif
          angr=angr+sckg1*g1-sckg2*g2
        else if(mm.gt.0) then
          if(iabs(mjm).gt.l.or.iabs(mjm1).gt.l1) then
            g1=0.d0
          else
            g1=gaunt(ll/2,lm0,l1m10)
            if((mj-mj1)/2.gt.0) g1=g1*(-1)**(mm)
          endif
          if(iabs(mjm2).gt.l.or.iabs(mjm3).gt.l1) then
            g2=0.d0
          else
            g2=gaunt(ll/2,lm,l1m1)
            if((mj-mj1)/2.gt.0) g2=g2*(-1)**mm
          endif
          angr=angr+(sckg1*g1-sckg2*g2)/sqr2
        endif
      endif
      end