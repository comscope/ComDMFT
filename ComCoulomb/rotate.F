
      subroutine rotate1(x,y,l,u0,key)

      use solid_mod
      implicit none
c     key=1 Y(m)=Sum_m1 D_m_m1 X_m1
c     key=2 Y(m)=Sum_m1 D_m1_m X_m1
c$$$  u0 point group matrix
      integer :: l,key,nls,lm,m,nls1,m1,l1m1,iwig,l1,l2m,ll2m
      double precision x(2*l+1),y(2*l+1),u0(maxwig),sum

      l1=l+1
      l2m=l1+l
      ll2m=l2m*l
      nls=(l2m-2)*ll2m/3+l1+ll2m
      lm=0
      if(key.eq.1) then
        do m=-l,l
          lm=lm+1
          nls1=nls+m
          sum=0.d0
          do m1=-l,l
            l1m1=l1+m1
            iwig=nls1+l2m*m1
            sum=sum+u0(iwig)*x(l1m1)
          enddo
          y(lm)=sum
        enddo
      else if(key.eq.2) then
        do m=-l,l
          lm=lm+1
          nls1=nls+l2m*m
          sum=0.d0
          do m1=-l,l
            l1m1=l1+m1
            iwig=nls1+m1
            sum=sum+u0(iwig)*x(l1m1)
          enddo
          y(lm)=sum
        enddo
      endif
      end


      subroutine rotate2(x,y,l,ii,u0,key,key1)
c$$$  complex version of rotate1
      use solid_mod
      implicit none

c     key=1 Y(m)=Sum_m1 D_m_m1 X_m1
c     key=2 Y(m)=Sum_m1 D_m1_m X_m1
c     key1=1 conjugated rotation  (D*)
c     key1=2 usual                (D )
      integer :: key,key1,l,ii,jj,iwig0,lm,m,nls1,l1m1,m1,iwig
      complex*16 x(*),y(*),u0(maxwig_j),sum
      jj=l+l+ii
      if(key1.eq.1) then
        iwig0=(8*l**3+12*l**2+10*l+3)/3+l*(2*l+1)*ii
        if(key.eq.1) then
          lm=0
          do m=-jj,jj,2
            lm=lm+1
            nls1=iwig0+(ii+m)/2
            sum=(0.d0,0.d0)
            l1m1=0
            do m1=-jj,jj,2
              l1m1=l1m1+1
              iwig=nls1+(2*l+ii+1)*(ii+m1)/2
              sum=sum+conjg(u0(iwig))*x(l1m1)
            enddo
            y(lm)=sum
          enddo
        else if(key.eq.2) then
          lm=0
          do m=-jj,jj,2
            lm=lm+1
            nls1=iwig0+(2*l+ii+1)*(ii+m)/2
            sum=(0.d0,0.d0)
            l1m1=0
            do m1=-jj,jj,2
              l1m1=l1m1+1
              iwig=nls1+(ii+m1)/2
              sum=sum+conjg(u0(iwig))*x(l1m1)
            enddo
            y(lm)=sum
          enddo
        endif
      else if(key1.eq.2) then
        iwig0=(8*l**3+12*l**2+10*l+3)/3+l*(2*l+1)*ii
        if(key.eq.1) then
          lm=0
          do m=-jj,jj,2
            lm=lm+1
            nls1=iwig0+(ii+m)/2
            sum=(0.d0,0.d0)
            l1m1=0
            do m1=-jj,jj,2
              l1m1=l1m1+1
              iwig=nls1+(2*l+ii+1)*(ii+m1)/2
              sum=sum+u0(iwig)*x(l1m1)
            enddo
            y(lm)=sum
          enddo
        else if(key.eq.2) then
          lm=0
          do m=-jj,jj,2
            lm=lm+1
            nls1=iwig0+(2*l+ii+1)*(ii+m)/2
            sum=(0.d0,0.d0)
            l1m1=0
            do m1=-jj,jj,2
              l1m1=l1m1+1
              iwig=nls1+(ii+m1)/2
              sum=sum+u0(iwig)*x(l1m1)
            enddo
            y(lm)=sum
          enddo
        endif
      endif
      end
