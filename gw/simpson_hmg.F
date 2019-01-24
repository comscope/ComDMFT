      subroutine simpson_hmg(a,b,n,f,res)
      implicit none
	integer, intent(in) :: n
      real*8, intent(in) :: a,b,f(0:n)
      real*8, intent(out) :: res
	integer :: i
	res=f(0)+f(n)
	do i=1,n-1,2
	  res=res+4.d0*f(i)
	enddo
	do i=2,n-2,2
	  res=res+2.d0*f(i)
	enddo
	res=res*(b-a)/3.d0/n
      end
      
      subroutine simpson_ihmg_weights(x,k,w)
      implicit none
	integer, intent(in) :: k
      real*8, intent(in) :: x(0:k)
      real*8, intent(out) :: w(0:k)
	integer :: i,i0,i1,i2
	real*8 :: h1,h2,h12
	if(k==1) then
	  h1=x(1)-x(0)
	  w(0)=h1/2
	  w(1)=w(0)
	  return
	else if(mod(k,2)==1) then
	  h1=x(1)-x(0)
	  h2=x(2)-x(1)
	  h12=h1+h2
	  w(0)=h1*(2*h1+3*h2)/6.d0/h12
	  w(1)=h1*(h1+3*h2)/6.d0/h2
	  w(2)=-h1**3/6.d0/h12/h2
	  w(3:k)=0.d0
	  i0=1
	else
	  w=0.d0
	  i0=0
	endif
	do i=i0,k-2,2
	  i1=i+1
	  i2=i+2
	  h1=x(i1)-x(i)
	  h2=x(i2)-x(i1)
	  h12=h1+h2
	  w(i)=w(i)+h12*(2*h1-h2)/6.d0/h1
	  w(i1)=w(i1)+h12**3/6.d0/h1/h2
	  w(i2)=w(i2)+h12*(2*h2-h1)/6.d0/h2
	enddo
      end