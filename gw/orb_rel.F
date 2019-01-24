      subroutine orb_rel(l,i,mj,l1,i1,mj1,aorbx,aorby,aorbz,aorbxi,
     &	                 aorbyi,aorbzi)
	use manager_mod
      implicit none
	integer :: l,i,mj,l1,mj1,i1,ii,ii1,lmj,lmj1,mu,mjm,mjm1,mmu
	real*8 :: aorbx,aorby,aorbz,aorbxi,aorbyi,aorbzi,sckg,con
      aorbx=0.d0
      aorby=0.d0
      aorbz=0.d0
      aorbxi=0.d0
      aorbyi=0.d0
      aorbzi=0.d0
      if(l.ne.l1) return
      if(iabs(mj-mj1).gt.2) return
      ii=(3+i)/2
      ii1=(3+i1)/2
      lmj=l*(l+2)+(mj+3)/2         !!!  IV-37-7
      lmj1=l1*(l1+2)+(mj1+3)/2
      do mu=-1,1,2
         mjm=(mj-mu)/2
         if(iabs(mjm).gt.l) goto 1
         mjm1=(mj1-mu)/2
         if(iabs(mjm1).gt.l1) goto 1
         mmu=(3+mu)/2
         sckg=ckg2(ii,mmu,lmj)*ckg2(ii1,mmu,lmj1)
         if(mjm.eq.mjm1-1) then
            con=0.5d0*sckg*sqrt(dfloat(l1+mjm1)*dfloat(l1-mjm1+1))
            aorbx=aorbx+con
            aorbyi=aorbyi+con
         else if(mjm.eq.mjm1+1) then
            con=0.5d0*sckg*sqrt(dfloat(l1-mjm1)*dfloat(l1+mjm1+1))
            aorbx=aorbx+con
            aorbyi=aorbyi-con
         else if(mjm.eq.mjm1) then
            aorbz=aorbz+sckg*dfloat(mjm1)
         endif
1        continue
      enddo   !!! over mu
      end