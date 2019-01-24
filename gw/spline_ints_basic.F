      subroutine spline_ints_basic(t0,t1,h0,ti1,t,sign,aa)
	use manager_mod
	use parallel_mod
	use units_mod
      implicit none
	real*8, intent(in) :: t0,t1,h0,ti1,t,sign
	real*8, intent(out) :: aa(0:3,0:3)
	real*8 :: b,a,dx1,dx2,dx3,dx4,dx5,dx6,dx7
	b=t-h0
	a=ti1
	dx1=t1-t0
	dx2=dx1*(t1+t0)
	dx3=dx1*(t1**2+t1*t0+t0**2)
	dx4=dx2*(t1**2+t0**2)
	dx5=dx1*(t1**4+t1**3*t0+t1**2*t0**2+t1*t0**3+t0**4)
	dx6=dx1*(t1+t0)*(t1*t1-t1*t0+t0*t0)*(t1**2+t1*t0+t0**2)
	dx7=dx1*(t1**6+t1**5*t0+t1**4*t0**2+t1**3*t0**3+t1**2*t0**4
     &	    +t1*t0**5+t0**6)
	aa(0,0)=dx1
	aa(1,0)=b*dx1-0.5d0*dx2
	aa(2,0)=dx3/3-b*dx2+b**2*dx1
	aa(3,0)=-0.25d0*(-4*b**3*dx1+6*b**2*dx2-4*b*dx3+dx4)
	aa(0,1)=0.5d0*(dx2-2*a*dx1)
	aa(1,1)=-dx3/3+0.5d0*(a+b)*dx2-a*b*dx1
	aa(2,1)=0.25d0*dx4-(a+2*b)*dx3/3+0.5d0*b*(2*a+b)*dx2-a*b**2*dx1
	aa(3,1)=-0.2d0*dx5+0.25d0*(a+3*b)*dx4-b*(a+b)*dx3
     &	    +0.5d0*b*b*(3*a+b)*dx2-a*b**3*dx1
	aa(0,2)=dx3/3-a*dx2+a**2*dx1
	aa(1,2)=-0.25d0*dx4+(2*a+b)*dx3/3-0.5d0*a*(a+2*b)*dx2+a*a*b*dx1
	aa(2,2)=0.2d0*dx5-0.5d0*(a+b)*dx4+(a*a+4*b*a+b*b)*dx3/3
     &	   -a*b*(a+b)*dx2+a*a*b*b*dx1
	aa(3,2)=-dx6/6+0.2d0*(2*a+3*b)*dx5-0.25d0*(a*a+6*a*b+3*b*b)*dx4
     &	    +b*(3*a*a+6*b*a+b*b)*dx3/3-0.5d0*a*b*b*(3*a+2*b)*dx2
     &        +a*a*b**3*dx1
	aa(0,3)=0.25d0*(-4*a**3*dx1+6*a**2*dx2-4*a*dx3+dx4)
	aa(1,3)=-0.2d0*dx5+0.25d0*(3*a+b)*dx4-a*(a+b)*dx3
     &	    +0.5d0*a*a*(a+3*b)*dx2-a**3*b*dx1
	aa(2,3)=dx6/6-0.2d0*(3*a+2*b)*dx5+0.25d0*(3*a*a+6*b*a+b*b)*dx4
     &	   -a*(a*a+6*b*a+3*b*b)*dx3/3+0.5d0*a*a*b*(2*a+3*b)*dx2
     &       -a**3*b*b*dx1
	aa(3,3)=-dx7/7+0.5d0*(a+b)*dx6-0.6d0*(a*a+3*b*a+b*b)*dx5
     &	    +0.25d0*(a**3+9*b*a*(a+b)+b**3)*dx4
     &        -a*b*(a*a+3*a*b+b*b)*dx3+1.5d0*a*a*b*b*(a+b)*dx2
     &        -(a*b)**3*dx1
	aa=sign*aa
      end