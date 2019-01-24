      subroutine a_cont_ferm(nrax,erax,f_in,f_out,esm)
	use atom_mod
	use etot_mod
	use manager_mod
	use parallel_mod
	use solid_mod
	use units_mod
      implicit none
      integer, intent(in) :: nrax
      real*8, intent(in) :: erax(nrax),esm
      complex*16, intent(in) :: f_in(0:n_omega)
      complex*16, intent(out) :: f_out(nrax)
	integer :: i
	real*8 :: chk
	complex*16 :: om0
      real*8, allocatable :: s0(:,:),s0_cheb(:,:)
      complex*16, allocatable :: a(:)
      chk=maxval(abs(f_in))
      if(chk<1.d-12) then
        f_out=(0.d0,0.d0)
      else
        allocate(a(0:n_cont_frac))
	  if(n_cont_frac<=n_omega_exa) then
          call pade_fermion_chain(f_in,1,1,a,n_cont_frac)
        else
	    allocate(s0(2,0:n_omega))
	    allocate(s0_cheb(2,0:n_omega_asy-2))
	    s0(1,:)=real(f_in)
	    s0(2,:)=imag(f_in)
	    call dgemm('n','t',2,n_omega_asy-1,n_omega_asy-1,1.d0,
     &               s0(1,n_omega_exa+1),2,cheb_omega(0,1),
     &               n_omega_asy-1,0.d0,s0_cheb(1,0),2)
          call ac_fermion_chain(s0,s0_cheb,1,1,a,n_cont_frac)
          deallocate(s0,s0_cheb)
        endif
        do i=1,nrax
          om0=dcmplx(erax(i),esm)
		  call pade_fermion_chain_out(1,1,a,om0,f_out(i),n_cont_frac)
        enddo
	  deallocate(a)
	endif
      end