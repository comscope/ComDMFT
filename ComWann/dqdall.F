      double precision function dqdall(h,y,n)
      
c$$$  y =f(x) dx
c$$$  returns int_0^1 f(x) dx      

      
      implicit none
      integer, intent(in) :: n
      double precision, intent(in) :: y(0:*),h
      integer :: n0,i
      double precision :: fac,sum06,sum15,sum24,sum3

      if(n.eq.1) then
        dqdall=h/2.d0*(y(0)+y(1))
        return
      else if(n.eq.2) then
        dqdall=h/3.d0*(y(0)+y(2)+4.d0*y(1))
        return
      else if(n.eq.3) then
        dqdall=h/8.d0*(3.d0*(y(0)+y(3))+9.d0*(y(1)+y(2)))
        return
      else if(n.eq.4) then
        dqdall=h/45.d0*(14.d0*(y(0)+y(4))+64.d0*(y(1)+y(3))+24.d0*y(2))
        return
      else if(n.eq.5) then
        dqdall=h/288.d0*(95.d0*(y(0)+y(5))+375.d0*(y(1)+y(4))+
     &    250.d0*(y(2)+y(3)))
        return
      endif

      n0=mod(n,6)
      if(n0.eq.0) then
        dqdall=0.d0
      else if(n0.eq.1) then
        fac=h/60480.d0
        dqdall=fac*(19087.d0*y(0)+65112.d0*y(1)-46461.d0*y(2)+
     &    37504.d0*y(3)-20211.d0*y(4)+6312.d0*y(5)-
     &    863.d0*y(6))
      else if(n0.eq.2) then
        fac=h/3780.d0
        dqdall=fac*(1139.d0*y(0)+5640.d0*y(1)+33.d0*y(2)+
     &    1328.d0*y(3)-807.d0*y(4)+264.d0*y(5)-
     &    37.d0*y(6))
      else if(n0.eq.3) then
        fac=h/2240.d0
        dqdall=fac*(685.d0*y(0)+3240.d0*y(1)+1161.d0*y(2)+
     &    2176.d0*y(3)-729.d0*y(4)+216.d0*y(5)-
     &    29.d0*y(6))
      else if(n0.eq.4) then
        fac=h/945.d0
        dqdall=fac*(286.d0*y(0)+1392.d0*y(1)+384.d0*y(2)+
     &    1504.d0*y(3)+174.d0*y(4)+48.d0*y(5)-
     &    8.d0*y(6))
      else if(n0.eq.5) then
        fac=h/12096.d0
        dqdall=fac*(3715.d0*y(0)+17400.d0*y(1)+6375.d0*y(2)+
     &    16000.d0*y(3)+11625.d0*y(4)+5640.d0*y(5)-
     &    275.d0*y(6))
      endif

      fac=h/140.d0
      sum06=0.d0
      sum15=0.d0
      sum24=0.d0
      sum3=0.d0
      do i=n0,n-6,6
        sum06=sum06+y(i)+y(i+6)
        sum15=sum15+y(i+1)+y(i+5)
        sum24=sum24+y(i+2)+y(i+4)
        sum3=sum3+y(i+3)
      enddo
      dqdall=dqdall+fac*(41.d0*sum06+216.d0*sum15+27.d0*sum24+
     &  272.d0*sum3)
      end
