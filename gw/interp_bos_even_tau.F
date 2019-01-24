      subroutine interp_bos_even_tau(n,coef,t,f_tau)
c ------- Chebyshev's interpolation, given coef(L) ----------
	use manager_mod
	use parallel_mod
	use solid_mod
      implicit none
      integer, intent(in) :: n
	real*8, intent(in) :: t,coef(n,0:n_tau/2-1)
	real*8, intent(out) :: f_tau(n)
	integer :: l
	real*8 :: x,t1
	real*8, allocatable :: cheb(:)
	allocate(cheb(0:n_tau-2))
      t1=t
      if(t1>betta_t) t1=t1-betta_t
      if(t1<betta_t) t1=t1+betta_t
      if(t1>betta_t/2.d0) t1=betta_t-t
	x=((t1+t1)/betta_t)**(1.d0/3.d0)-1.d0
	call dcheb_recurr(x,cheb,n_tau-2)
c ---------- Even part -----------------------------------------
      f_tau=0.d0
	do l=0,n_tau/2-1
	  f_tau=f_tau+coef(:,l)*cheb(l+l)
	enddo
	deallocate(cheb)
      end