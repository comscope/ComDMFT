      subroutine sigx_k_solid
	use atom_mod
	use manager_mod
	use parallel_mod
	use solid_mod
      implicit none
	integer :: ispin,k,ind_k
	real*8, allocatable :: v_rs_mm(:,:,:)
	sigx_solid_k=0.d0
	allocate(v_rs_mm(n_pbmtm,n_pbmtm,ndim_nrr(me3_k+1)))
	call v_mm_from_q(v_rs_mm)
	do ispin=1,nspin
	  call sigma_prelim(ispin,v_rs_mm)
	  do ind_k=1,ndim3_k(me3_k+1)
	    k=n3_mpi_k(me3_k+1)+ind_k
	    call sigma_loc_fin(ind_k,k,ispin)
	    call sigma_mt_mt_fin(ind_k,k,ispin)
	    call sigma_mt_int_fin(ind_k,k,ispin)
	    call sigma_int_int_fin(ind_k,k,ispin)
	  enddo
	  deallocate(pvv,pvg,pgg)
	enddo
	deallocate(v_rs_mm)
	if(ubi==' gw'.or.ubi=='gwg') call etot_x_from_pw(1)
      end
