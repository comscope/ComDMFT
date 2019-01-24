      subroutine psi_function_to_real_mesh(k,nbas,psi_g,psi_r)
	use atom_mod
	use manager_mod
	use parallel_mod
	use solid_mod
	use units_mod
      implicit none
	integer, intent(in) :: nbas,k
	complex*16, intent(in) :: psi_g(nbas)
	complex*16, intent(out) :: psi_r(nr_full)
	integer :: i,i0,ind
      psi_r=(0.d0,0.d0)
	do i=1,nbas
	  i0=indgb(i,k)
	  call zone1_number(gbs(:,i0),rbas,nrdiv,ind)
	  psi_r(ind)=psi_g(i)
	enddo
	call fft3(nrdiv(1),nrdiv(2),nrdiv(3),1,psi_r,1)
      end