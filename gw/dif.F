      subroutine dif(f,n,fd)
      implicit none
	integer :: n,irad
      real*8 :: f(0:n),fd(0:n),con1,con2,con3,c2,c3
      data con3,con2,con1/0.08730158730159d0,-0.26587301587302d0,
     &                    -0.23015873015873d0/
      fd(3)=con3*(f(0)-f(6))+con2*(f(1)-f(5))+con1*(f(2)-f(4))
      c3=-49.d0*fd(3)/397.d0-(27.d0*(f(0)-f(6))+8.d0*(f(1)-f(5))+f(2)-
     &   f(4))/1588.d0
      c2=(5.d0*(f(0)+f(6))-3.d0*(f(2)+f(4))-4.d0*f(3))/84.d0
      fd(0)=fd(3)-6.d0*c2+27.d0*c3
      fd(1)=fd(3)-4.d0*c2+12.d0*c3
      fd(2)=fd(3)-2.d0*c2+3.d0*c3
      do irad=4,n-4
         fd(irad)=con3*(f(irad-3)-f(irad+3))+
     &            con2*(f(irad-2)-f(irad+2))+
     &            con1*(f(irad-1)-f(irad+1))
      enddo
      fd(n-3)=con3*(f(n-6)-f(n))+con2*(f(n-5)-f(n-1))+
     &        con1*(f(n-4)-f(n-2))
      c3=-49.d0*fd(n-3)/397.d0-(27.d0*(f(n-6)-f(n))+8.d0*(f(n-5)-f(n-1))
     &                          +f(n-4)-f(n-2))/1588.d0
      c2=(5.d0*(f(n-6)+f(n))-3.d0*(f(n-4)+f(n-2))-4.d0*f(n-3))/84.d0
      fd(n-2)=fd(n-3)+2.d0*c2+3.d0*c3
      fd(n-1)=fd(n-3)+4.d0*c2+12.d0*c3
      fd(n)=fd(n-3)+6.d0*c2+27.d0*c3
      end