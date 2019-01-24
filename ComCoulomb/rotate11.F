      subroutine rotate(x2,y2,z2,x1,y1,z1,u,key)

! given 3 by 3 rotation matrix u, totate {x2y2z2} vector to {x1,y1,z1}

      use units_mod
      implicit none
      integer :: key
      double precision :: u(3,3),yg(3),yg1(3),r,x1,y1,z1,x2,y2,z2,s3
      r=x2*x2+y2*y2+z2*z2

      if(r.lt.1.0d-8)then
        x1=0.d0
        y1=0.d0
        z1=0.d0
        return
      endif
      s3=sqrt(0.75d0/pi/r)
      yg(1)=s3*y2
      yg(2)=s3*z2
      yg(3)=s3*x2
      if(key.eq.1)call prd1(yg,u,yg1) !!! inversional ig-operation
      if(key.eq.2)call prd2(yg,u,yg1) !!!      direct ig-operation
      x1=yg1(3)/s3
      y1=yg1(1)/s3
      z1=yg1(2)/s3
      end



      subroutine rotate11(x,y,nlm,u0,key)
      use solid_mod
      implicit none
c     key=1 Y(m)=Sum_m1 D_m_m1 X_m1
c     key=2 Y(m)=Sum_m1 D_m1_m X_m1


      integer :: nlm,key,lmsf,l,nls,m,lm,nls1,m1,l1m1,iwig,lget
      double precision :: x(nlm),y(nlm),u0(maxwig),sum
      lmsf=lget(nlm)
      do l=0,lmsf
        nls=l*(2*l-1)*(2*l+1)/3+l+1
        if(key.eq.1) then
          do m=-l,l
            lm=l*(l+1)+m+1
            nls1=nls+m
            sum=0.d0
            do m1=-l,l
              l1m1=l*(l+1)+m1+1
              iwig=nls1+(2*l+1)*(l+m1)
              sum=sum+u0(iwig)*x(l1m1)
            enddo
            y(lm)=sum
          enddo
        else if(key.eq.2) then
          do m=-l,l
            lm=l*(l+1)+m+1
            nls1=nls+(2*l+1)*(l+m)
            sum=0.d0
            do m1=-l,l
              l1m1=l*(l+1)+m1+1
              iwig=nls1+m1
              sum=sum+u0(iwig)*x(l1m1)
            enddo
            y(lm)=sum
          enddo
        endif
      enddo
      end
