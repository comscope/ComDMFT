      subroutine sig_omega0_g1(n,s_co,z_0,z_fac)
	use atom_mod
	use manager_mod
	use parallel_mod
	use solid_mod
	use units_mod
      implicit none
	integer, intent(in) :: n
	complex*16, intent(in) :: s_co(maxel_loc,maxel_loc,ndim3_omega)
	complex*16, intent(out) :: z_0(n,n),z_fac(n,n)
	integer :: i,j,ind_omega,i_omega
	complex*16, allocatable :: tmp(:,:,:),tmp0(:,:)
	allocate(tmp(n,n,2))
	allocate(tmp0(n,0:n_omega))
	tmp=(0.d0,0.d0)
	tmp0=(0.d0,0.d0)
	do ind_omega=1,ndim3_omega
	  i_omega=me3_tau*ndim3_omega+ind_omega-1
	  do j=1,n
	    tmp0(j,i_omega)=s_co(j,j,ind_omega)
        enddo
	  if(i_omega<=1) then
	    i=i_omega+1
	    tmp(:,:,i)=s_co(1:n,1:n,ind_omega)
        endif
      enddo   !! over ind_omega
	if(nproc_tau/=1) then
	  call dgop(tmp0,2*n*(n_omega+1),'  +',comm_pnk)
	  call dgop(tmp,4*n*n,'  +',comm_pnk)
	endif
	call sigma_0_z_ren(n,n,tmp,tmp0)
	z_0=tmp(:,:,1)
	z_fac=tmp(:,:,2)
	deallocate(tmp,tmp0)
      end