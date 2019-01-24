      subroutine getlimj(limj,l,i,mj,li,key)
c$$$      l: l quantum number
c$$$      i: (j-l)*2 quantum mumber
c$$$      mj: mj*2
c$$$      li: hilbert space index by (l, j) quantum number 
c     key=0 l,i,mj,li from limj
c     key=1 limj from l,i,mj
      
      implicit none
      integer, intent(in) :: key
c     integer, intent(inout) :: limj,l,i,mj,li
      integer :: limj,l,i,mj,li
      if(key.eq.0) then
        l=sqrt((dfloat(limj)-0.1d0)/2.d0)
        if(limj-2*l*l.le.l+l) then
          i=-1
        else
          i=1
        endif
        li=(5+i)/2+2*(l-1)
        mj=2*(limj-l*(2*l+2+i)-1)-i
      else
        limj=l*(2*l+2+i)+(mj+i)/2+1
      endif
      end

      subroutine getlm(lm,l,m,key)
c$$$      l: l quantum number
c$$$      m: m starting from 1
c     key=0 l,i,mj,li from limj
c     key=1 limj from l,i,mj
      
      implicit none
      integer, intent(in) :: key
      integer :: lm,l,m
      integer, external :: lget
      if(key.eq.0) then
        l=lget(lm)
        m=lm-l**2
      else
        lm=l**2+m
      endif
      end  
      

      
