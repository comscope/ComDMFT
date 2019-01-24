      subroutine sp_rel(l,i,mj,l1,i1,mj1,spx,spy,spz,spxi,spyi,spzi)
	use manager_mod
      implicit none
	integer :: l,i,mj,l1,i1,mj1,ii1,lmj,lmj1,ii
	real*8 :: spx,spy,spz,spxi,spyi,spzi,sckg11,sckg21,sckg12,sckg22
c -----------  See IV-21 ---------------------------------------------
      spx=0.d0
      spy=0.d0
      spz=0.d0
      spxi=0.d0
      spyi=0.d0
      spzi=0.d0
      if(l.ne.l1) return
      if(iabs(mj-mj1).gt.2) return
      ii=(3+i)/2
      ii1=(3+i1)/2
      lmj=l*(l+2)+(mj+3)/2         !!!  IV-37-7
      lmj1=l1*(l1+2)+(mj1+3)/2
      sckg11=ckg2(ii,1,lmj)
      sckg21=ckg2(ii,2,lmj)
      sckg12=ckg2(ii1,1,lmj1)
      sckg22=ckg2(ii1,2,lmj1)
	if(mj.eq.mj1-2) then
	  spx=spx+0.5d0*sckg11*sckg22
	  spyi=spyi+0.5d0*sckg11*sckg22
	else if(mj.eq.mj1+2) then
	  spx=spx+0.5d0*sckg21*sckg12
	  spyi=spyi-0.5d0*sckg21*sckg12
	else if(mj.eq.mj1) then
	  spz=spz+0.5d0*(sckg21*sckg22-sckg11*sckg12)
	endif
      end