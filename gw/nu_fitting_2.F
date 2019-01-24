      subroutine nu_fitting_2(f,fp,nb_fit,err,w0,w1,n)
	use manager_mod
      implicit none
      integer, intent(in) :: nb_fit,n
	real*8, intent(in) :: f(0:n_nu),w0(2),w1(2)
	real*8, intent(out) :: fp(2,nb_fit),err
	integer :: i_nu,i1,i2
	real*8 :: d,s,hh1,hh2,g,a,c1,c2,e1,e2,b,st,tt,t,green0_boson_nu
	real*8, allocatable :: g1(:),g2(:)
	allocate(g1(0:n_nu))
	allocate(g2(0:n_nu))
	hh1=(w1(1)-w0(1))/n
	hh2=(w1(2)-w0(2))/n
	d=1.d5
	do i1=0,n
	  e1=w0(1)+hh1*i1
	  do i_nu=0,n_nu
	    g1(i_nu)=green0_boson_nu(e1,w_nu(i_nu))
	  enddo
	  do i2=1,n
	    e2=w0(2)+hh2*i2
	    if(e2-e1<2*hh2) cycle
	    do i_nu=0,n_nu
	      g2(i_nu)=green0_boson_nu(e2,w_nu(i_nu))
	    enddo
	    a=e1**2*f(0)
	    b=(e1/e2)**2
	    st=0.d0
	    tt=0.d0
	    do i_nu=1,n_nu-1
	      s=f(i_nu)-a*g1(i_nu)
	      t=b*g1(i_nu)-g2(i_nu)
	      st=st+s*t
	      tt=tt+t*t
	    enddo
	    c2=-st/tt
	    c1=a-b*c2
	    s=0.d0
	    do i_nu=0,n_nu
	      g=c1*g1(i_nu)+c2*g2(i_nu)
	      s=s+(f(i_nu)-g)**2
	    enddo
	    if(s<d) then
	      d=s
	      fp(2,1)=e1
	      fp(2,2)=e2
	      fp(1,1)=c1
	      fp(1,2)=c2
	    endif
	  enddo
	enddo
	deallocate(g1,g2)
	err=sqrt(d/(n_nu+1))
      end