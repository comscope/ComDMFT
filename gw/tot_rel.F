      subroutine tot_rel(jj,mj,jj1,mj1,totx,toty,totz,totxi,totyi,totzi)
      implicit none
	integer :: jj,mj,jj1,mj1
	real*8 :: totx,toty,totz,totxi,totyi,totzi
c -----------  See V-79 ---------------------------------------------
      totx=0.d0
      toty=0.d0
      totz=0.d0
      totxi=0.d0
      totyi=0.d0
      totzi=0.d0
      if(jj.eq.jj1) then
        if(mj1.eq.mj) totz=dfloat(mj)/2.d0
        if(mj.eq.mj1+2) then
          totx=0.25d0*sqrt(dfloat(jj1-mj1)*dfloat(jj1+mj1+2))
          totyi=-totx
        else if(mj.eq.mj1-2) then
          totx=0.25d0*sqrt(dfloat(jj1+mj1)*dfloat(jj1-mj1+2))
          totyi=totx
        endif
      endif
      end