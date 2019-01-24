
      subroutine rotate(x2,y2,z2,x1,y1,z1,ur,key)

! given 3 by 3 rotation matrix u, totate {x2y2z2} vector to {x1,y1,z1}

      use comwann_mod
      implicit none
      integer :: key
      double precision :: ur(3,3),yg(3),yg1(3),rr,x1,y1,z1,x2,y2,z2,s3
      rr=x2*x2+y2*y2+z2*z2

      if(rr.lt.1.0d-8)then
        x1=0.d0
        y1=0.d0
        z1=0.d0
        return
      endif
      s3=sqrt(0.75d0/pi/rr)
      yg(1)=s3*y2
      yg(2)=s3*z2
      yg(3)=s3*x2
      if(key.eq.1)call prd1(yg,ur,yg1) !!! inversional ig-operation
      if(key.eq.2)call prd2(yg,ur,yg1) !!!      direct ig-operation
      x1=yg1(3)/s3
      y1=yg1(1)/s3
      z1=yg1(2)/s3
      end


      subroutine rotate1(x,y,l,u0,key)

      use comwann_mod
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
      use comwann_mod
      implicit none

c     key=1 Y(m)=Sum_m1 D_m_m1 X_m1
c     key=2 Y(m)=Sum_m1 D_m1_m X_m1
c     key1=1 conjugated rotation  (D*)
c     key1=2 usual                (D )
      integer :: key,key1,l,ii,jj,iwig0,lm,m,nls1,l1m1,m1,iwig
      complex*16 x(*),y(*),u0(maxwig_j),sumval
      jj=l+l+ii
      if(key1.eq.1) then
        iwig0=(8*l**3+12*l**2+10*l+3)/3+l*(2*l+1)*ii
        if(key.eq.1) then
          lm=0
          do m=-jj,jj,2
            lm=lm+1
            nls1=iwig0+(ii+m)/2
            sumval=(0.d0,0.d0)
            l1m1=0
            do m1=-jj,jj,2
              l1m1=l1m1+1
              iwig=nls1+(2*l+ii+1)*(ii+m1)/2
              sumval=sumval+conjg(u0(iwig))*x(l1m1)
            enddo
            y(lm)=sumval
          enddo
        else if(key.eq.2) then
          lm=0
          do m=-jj,jj,2
            lm=lm+1
            nls1=iwig0+(2*l+ii+1)*(ii+m)/2
            sumval=(0.d0,0.d0)
            l1m1=0
            do m1=-jj,jj,2
              l1m1=l1m1+1
              iwig=nls1+(ii+m1)/2
              sumval=sumval+conjg(u0(iwig))*x(l1m1)
            enddo
            y(lm)=sumval
          enddo
        endif
      else if(key1.eq.2) then
        iwig0=(8*l**3+12*l**2+10*l+3)/3+l*(2*l+1)*ii
        if(key.eq.1) then
          lm=0
          do m=-jj,jj,2
            lm=lm+1
            nls1=iwig0+(ii+m)/2
            sumval=(0.d0,0.d0)
            l1m1=0
            do m1=-jj,jj,2
              l1m1=l1m1+1
              iwig=nls1+(2*l+ii+1)*(ii+m1)/2
              sumval=sumval+u0(iwig)*x(l1m1)
            enddo
            y(lm)=sumval
          enddo
        else if(key.eq.2) then
          lm=0
          do m=-jj,jj,2
            lm=lm+1
            nls1=iwig0+(2*l+ii+1)*(ii+m)/2
            sumval=(0.d0,0.d0)
            l1m1=0
            do m1=-jj,jj,2
              l1m1=l1m1+1
              iwig=nls1+(ii+m1)/2
              sumval=sumval+u0(iwig)*x(l1m1)
            enddo
            y(lm)=sumval
          enddo
        endif
      endif
      end
