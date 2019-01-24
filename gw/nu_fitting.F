      subroutine nu_fitting(f,fp,err)
c     Fitting F(v)=A/(v^2+B^2)    
	use manager_mod
      implicit none
	real*8, intent(in) :: f(0:n_nu)
	real*8, intent(out) :: fp(2),err
	integer :: i_nu,i1_nu,nu_fit
	real*8 :: d,s,kap2,w2,w22,g,green0_boson_nu
	d=1.d5
	do i_nu=1,n_nu
	  w2=w_nu(i_nu)**2
	  kap2=f(i_nu)*w2/(f(0)-f(i_nu))
	  s=0.d0
	  do i1_nu=1,n_nu
	    if(i1_nu==i_nu) cycle
	    w22=w_nu(i1_nu)**2
	    s=s+abs(f(i1_nu)-f(0)*kap2/(w22+kap2))
	  enddo
	  if(s<d) then
	    d=s
	    nu_fit=i_nu
	    fp(1)=f(0)*kap2
	    fp(2)=sqrt(kap2)
	  endif
	enddo
	s=0.d0
	do i_nu=0,n_nu
	  g=fp(1)*green0_boson_nu(fp(2),w_nu(i_nu))
	  s=s+(f(i_nu)-g)**2
	enddo
	err=sqrt(s/(n_nu+1))
      end