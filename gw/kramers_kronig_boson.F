      subroutine kramers_kronig_boson(h_xi_rax,nrax,n,wnu,s,s1,chk)
	use atom_mod
	use manager_mod
	use models_mod
	use parallel_mod
	use solid_mod
	use units_mod
      implicit none
      integer, intent(in) :: nrax,n
      real*8, intent(in) :: wnu(0:n),h_xi_rax
	complex*16, intent(in) :: s(0:nrax),s1(0:nrax)
	complex*16, intent(out) :: chk(0:n)
	integer :: i,i_nu
	real*8 :: vi2,x,v,dqdall
	real*8, allocatable :: f(:),g(:)
	allocate(f(0:nrax))
	allocate(g(0:nrax))
	do i_nu=0,n
	  vi2=wnu(i_nu)**2
	  if(nspin==1) then
	    f=0.d0
	    do i=0,nrax
	      if(i_nu==0.and.i==0) cycle
	      x=i*h_xi_rax
	      v=x*x
	      f(i)=x*v*imag(s(i))/(v*v+vi2)
	    enddo
	    chk(i_nu)=4.d0*dqdall(h_xi_rax,f,nrax)/pi
	  else if(nspin==2) then
	    f=0.d0
	    g=0.d0
	    do i=0,nrax
	      if(i_nu==0.and.i==0) cycle
	      x=i*h_xi_rax
	      v=x*x
	      f(i)=x*v*(imag(s(i))-imag(s1(i)))/(v*v+vi2)
	      g(i)=-x*v*(real(s(i))-real(s1(i)))/(v*v+vi2)
	    enddo
	    x=2.d0*dqdall(h_xi_rax,f,nrax)/pi
	    v=2.d0*dqdall(h_xi_rax,g,nrax)/pi
	    chk(i_nu)=dcmplx(x,v)
	  endif
	enddo
	deallocate(f,g)
      end