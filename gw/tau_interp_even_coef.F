      subroutine tau_interp_even_coef(fcoef,t0,fout,n,n0,n1)
c ------- Chebyshev's interpolation, given f_coef(L) ----------
	use manager_mod
	use parallel_mod
	use solid_mod
      implicit none
      integer, intent(in) :: n,n0,n1
	real*8, intent(in) :: fcoef(n0,n0,0:n_tau/2-1),t0
	real*8, intent(out) :: fout(n1,n1)
	integer :: l
	real*8 :: x,t1
	real*8, allocatable :: cheb(:)
	allocate(cheb(0:n_tau-2))
      t1=t0
      if(t1>betta_t/2.d0) t1=betta_t-t0
	x=((t1+t1)/betta_t)**(1.d0/3.d0)-1.d0
	call dcheb_recurr(x,cheb,n_tau-2)
	fout=0.d0
	do l=0,n_tau/2-1
	  fout(1:n,1:n)=fout(1:n,1:n)+fcoef(1:n,1:n,l)*cheb(l+l)
	enddo
	deallocate(cheb)
      end