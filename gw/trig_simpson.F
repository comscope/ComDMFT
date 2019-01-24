      subroutine trig_simpson(x,y,n,k,sn,cs,info)
      implicit none
	integer, intent(in) :: n
	real*8, intent(in) :: x(0:n),y(0:n),k
	integer, intent(out) :: info
	real*8, intent(out) :: sn,cs
	integer :: i
	real*8 :: k2,k3,c,s,a,b,fs1,gs1,hs1,fs3,gs3,hs3,fc1,gc1,hc1,fc3,
     &	      gc3,hc3,x1,x2,x3,a1,a2,a3,aa,bb,cc
	k2=2.d0/(k*k)
	k3=k2/k
	info=0
	sn=0.d0
	cs=0.d0
	c=k*x(0)
	s=sin(c)
	c=cos(c)
	a=k3-x(0)**2/k
	b=k2*x(0)
	fs3=a*c+b*s
	gs3=(s/k-x(0)*c)/k
	hs3=-c/k
	fc3=-a*s+b*c
	gc3=(c/k+x(0)*s)/k
	hc3=s/k
	do i=2,n,2
	  fs1=fs3
	  gs1=gs3
	  hs1=hs3
	  fc1=fc3
	  gc1=gc3
	  hc1=hc3
	  c=k*x(i)
	  s=sin(c)
	  c=cos(c)
	  a=k3-x(i)**2/k
	  b=k2*x(i)
	  fs3=a*c+b*s
	  gs3=(s/k-x(i)*c)/k
	  hs3=-c/k
	  fc3=-a*s+b*c
	  gc3=(c/k+x(i)*s)/k
	  hc3=s/k
	  x1=x(i-2)
	  x2=x(i-1)
	  x3=x(i)
	  a1=y(i-2)/(x1-x2)/(x1-x3)
	  a2=y(i-1)/(x2-x1)/(x2-x3)
	  a3=y(i)/(x3-x1)/(x3-x2)
	  aa=a1+a2+a3
	  bb=-(a2+a3)*x1-(a1+a3)*x2-(a1+a2)*x3
	  cc=a1*x2*x3+a2*x1*x3+a3*x1*x2
	  sn=sn+aa*(fs3-fs1)+bb*(gs3-gs1)+cc*(hs3-hs1)
	  cs=cs+aa*(fc3-fc1)+bb*(gc3-gc1)+cc*(hc3-hc1)
	enddo
	if(mod(n,2)/=0) info=1
      end