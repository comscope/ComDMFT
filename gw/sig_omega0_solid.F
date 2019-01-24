      subroutine sig_omega0_solid(n,ind_k,ispin)
	use atom_mod
	use manager_mod
	use parallel_mod
	use solid_mod
	use units_mod
      implicit none
	integer, intent(in) :: n,ind_k,ispin
	integer :: i,j,ind_omega,i_omega,n0
	complex*16, allocatable :: tmp(:,:,:),tmp0(:,:)
	n0=nbndf
	allocate(tmp(n0,n0,2))
	allocate(tmp0(n0,0:n_omega))
	tmp=(0.d0,0.d0)
	tmp0=(0.d0,0.d0)
	do ind_omega=1,ndim3_omega
	  i_omega=me3_tau*ndim3_omega+ind_omega-1
	  do j=1,n
	    tmp0(j,i_omega)=dcmplx(sig_c_omega(j,j,1,ind_omega,ind_k,
     &                                       ispin),
     &                           sig_c_omega(j,j,2,ind_omega,ind_k,
     &                                       ispin))
        enddo
	  if(i_omega<=1) then
	    i=i_omega+1
	    call ferm_unpack_omega(tmp(1,1,i),
     &                           sig_c_omega(1,1,1,ind_omega,ind_k,
     &                                       ispin),n,n0,n0)
        endif
      enddo   !! over ind_omega
	if(nproc_tau/=1) then
	  call dgop(tmp0,2*n0*(n_omega+1),'  +',comm_pnk)
	  call dgop(tmp,2*n0**2*2,'  +',comm_pnk)
	endif
	call sigma_0_z_ren(n,n0,tmp,tmp0)
	call pack_hermit(tmp,sig_c_0(1,1,ind_k,ispin),n,n0,n0,
     &                 (0.d0,0.d0),(1.d0,0.d0))
	call pack_hermit(tmp(1,1,2),z_ren(1,1,ind_k,ispin),n,n0,
     &                 n0,(0.d0,0.d0),(1.d0,0.d0))
	deallocate(tmp,tmp0)
      end