      subroutine interp_ferm_tau(n,coef,t,f_tau)
c ------- Chebyshev's interpolation, given coef(L) ----------
	use manager_mod
	use parallel_mod
	use solid_mod
      implicit none
      integer, intent(in) :: n
	real*8, intent(in) :: t,coef(n,0:n_tau/2-1,2)
	real*8, intent(out) :: f_tau(n)
	integer :: l
	real*8 :: x,t1,t2
	real*8, allocatable :: cheb(:),a(:),b(:)
	allocate(a(n))
	allocate(b(n))
	allocate(cheb(0:n_tau-2))
      t1=t
      if(t<0.d0) t1=t+betta_t
      if(t>betta_t) t1=t-betta_t
      t2=t1
      if(t2>betta_t/2.d0) t2=betta_t-t2
	x=((t2+t2)/betta_t)**(1.d0/3.d0)-1.d0
	call dcheb_recurr(x,cheb,n_tau-2)
c ---------- Odd part -----------------------------------------
      a=0.d0
	do l=0,n_tau/2-2
	  a=a+coef(:,l,1)*cheb(l+l+1)
	enddo
c ---------- Even part -----------------------------------------
      b=0.d0
	do l=0,n_tau/2-1
	  b=b+coef(:,l,2)*cheb(l+l)
	enddo
	if(t1<=betta_t/2) f_tau=(a+b)/2
	if(t1>betta_t/2) f_tau=(b-a)/2
      if(t<0.d0) f_tau=-f_tau
      if(t>betta_t) f_tau=-f_tau
	deallocate(cheb,a,b)
      end